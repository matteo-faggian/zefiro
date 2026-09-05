"""Il livello di servizio: tutto cio' che l'API sa fare.

Sta separato dalle rotte HTTP di proposito. Le funzioni qui dentro non sanno
nulla di richieste, risposte o codici di stato: prendono dati e restituiscono
dati. Il risultato pratico e' che sono chiamabili dai test senza un client HTTP,
e che una futura interfaccia diversa (una CLI arricchita, uno script) userebbe
esattamente lo stesso codice.

Nessuna fisica vive qui. Ogni funzione chiama `derive`, `evaluate_l0`,
`build_and_export` — le stesse della riga di comando.
"""
from __future__ import annotations

import logging
import math
from datetime import datetime, timezone
from hashlib import blake2b
from json import dumps
from typing import Any, Callable

from zefiro.cli import runs_root
from zefiro.config import load_design_vector, load_material, load_operating_point
from zefiro.doe import latin_hypercube
from zefiro.feed import HP_TO_W, burst_duration, compressor_bounds, tank_blowdown
from zefiro.geometry import BuildOptions, build_and_export
from zefiro.geometry.build import wall_profile_mm
from zefiro.geometry.parameters import (
    DESIGN_BOUNDS,
    INTEGER_PARAMETERS,
    check_manufacturability,
    derive,
)
from zefiro.l0.equilibrium import chamber_state
from zefiro.l0.mixture import MixtureModel
from zefiro.opt.objectives import objectives_l0
from zefiro.schemas import DesignVector, OperatingPoint, ZefiroError
from zefiro.store import RunStore, env_fingerprint, make_run_id
from zefiro.store.runid import git_revision
from zefiro.thermal import adiabatic_wall_temperature, bartz_h_gas, transient_slab
from zefiro.units import BAR, DEG, MM

log = logging.getLogger("zefiro.web")

RHO_316L_NOMINAL = 7980.0   # kg/m3, PLACEHOLDER: solo per stimare una massa

#: I 12 parametri liberi, con unita' di VISUALIZZAZIONE e fattore verso SI.
#: Insieme a config.py e' l'unico posto del progetto dove esistono unita' non SI.
DESIGN_FIELDS: list[dict[str, Any]] = [
    {"name": "p_c", "label": "Pressione di camera", "unit": "bar", "si": BAR, "step": 0.01,
     "group": "ciclo",
     "help": "Limitata dall'alimentazione: p_c piu' il Dp dell'iniettore deve stare sotto la pressione a monte a fine raffica."},
    {"name": "phi_core", "label": "Rapporto di equivalenza", "unit": "–", "si": 1.0, "step": 0.005,
     "group": "ciclo",
     "help": "Del core, che vede tutta l'aria. L'ottimo non e' 1: il massimo di c* sta poco sopra."},
    {"name": "f_film", "label": "Frazione al film cooling", "unit": "–", "si": 1.0, "step": 0.005,
     "group": "ciclo",
     "help": "A 5 s la camera non ne ha bisogno, e il film iniettato in testa non arriva alla gola."},
    {"name": "Dc_over_Dt", "label": "Rapporto di contrazione", "unit": "Dc/Dt", "si": 1.0, "step": 0.05,
     "group": "camera",
     "help": "Alto significa bassa velocita' in camera e piu' tempo di residenza, al prezzo di massa e area da raffreddare."},
    {"name": "Lc_over_Dc", "label": "Lunghezza di camera", "unit": "Lc/Dc", "si": 1.0, "step": 0.05,
     "group": "camera", "help": "Governa il tempo di residenza."},
    {"name": "conv_half_angle", "label": "Semiangolo del convergente", "unit": "°", "si": DEG, "step": 0.5,
     "group": "camera", "help": "Perdite di pressione totale contro lunghezza."},
    {"name": "N_inj", "label": "Elementi di iniezione", "unit": "n", "si": 1.0, "step": 1.0,
     "group": "iniezione",
     "help": "Intero. Piu' elementi danno mixing piu' fine, ma peggiorano Dp e fabbricabilita'."},
    {"name": "d_ox_ratio", "label": "Diametro foro aria", "unit": "d/Dc", "si": 1.0, "step": 0.002,
     "group": "iniezione", "help": "Determina la velocita' di iniezione dell'aria."},
    {"name": "fuel_vel_ratio", "label": "Velocita' GPL / aria", "unit": "–", "si": 1.0, "step": 0.05,
     "group": "iniezione",
     "help": "E' IL parametro che governa il mixing, attraverso il rapporto delle quantita' di moto."},
    {"name": "theta_swirl", "label": "Angolo di swirl", "unit": "°", "si": DEG, "step": 1.0,
     "group": "iniezione", "help": "Accorcia la fiamma, costa pressione totale."},
    {"name": "plug_trunc", "label": "Troncamento del plug", "unit": "–", "si": 1.0, "step": 0.01,
     "group": "ugello",
     "help": "1.0 e' il plug completo. Troncare accorcia e alleggerisce pagando prestazione."},
    {"name": "t_wall", "label": "Spessore di parete", "unit": "mm", "si": MM, "step": 0.1,
     "group": "ugello",
     "help": "A pozzo termico ispessire aiuta finche' satura; con raffreddamento sul dorso e' solo resistenza in piu'."},
]

OPERATING_FIELDS: list[dict[str, Any]] = [
    {"name": "mdot_air_max", "label": "Portata d'aria", "unit": "g/s", "si": 1e-3,
     "help": "Massa estraibile dal serbatoio divisa per la durata."},
    {"name": "burn_time", "label": "Durata", "unit": "s", "si": 1.0,
     "help": "Decide se la parete puo' fare da pozzo termico o va raffreddata."},
    {"name": "T_air_in", "label": "Temperatura aria", "unit": "K", "si": 1.0,
     "help": "A valle del riduttore. Non e' la temperatura ambiente: l'espansione raffredda."},
    {"name": "T_fuel_in", "label": "Temperatura GPL", "unit": "K", "si": 1.0,
     "help": "A valle del riduttore."},
    {"name": "cd_injector_ox", "label": "Cd iniettore aria", "unit": "–", "si": 1.0,
     "help": "Senza, i vincoli sul Dp non sono valutabili e restano fuori dal calcolo."},
    {"name": "cd_injector_fuel", "label": "Cd iniettore GPL", "unit": "–", "si": 1.0,
     "help": "Idem."},
]

_DESIGN_BY_NAME = {f["name"]: f for f in DESIGN_FIELDS}
_OPERATING_BY_NAME = {f["name"]: f for f in OPERATING_FIELDS}


# --------------------------------------------------------------------------- #
# Conversioni al confine
# --------------------------------------------------------------------------- #
def design_to_si(values: dict[str, float]) -> dict[str, float]:
    out: dict[str, float] = {}
    mancanti = [n for n in _DESIGN_BY_NAME if n not in values]
    if mancanti:
        raise ValueError(f"parametri di progetto mancanti: {', '.join(sorted(mancanti))}")
    for name, f in _DESIGN_BY_NAME.items():
        out[name] = float(values[name]) * f["si"]
    for name in INTEGER_PARAMETERS:
        out[name] = float(round(out[name]))
    return out


def operating_from(overrides: dict[str, float | None]) -> OperatingPoint:
    """Punto operativo di `config/`, con le sovrascritture dell'interfaccia.

    Le sovrascritture NON scrivono sui file. Cambiare il progetto si fa in
    `config/`, sotto controllo di versione, dove il diff si vede in revisione.
    Quello che l'interfaccia offre e' esplorazione, non modifica del progetto.
    """
    base = load_operating_point()
    kw: dict[str, Any] = {}
    for name, f in _OPERATING_BY_NAME.items():
        v = overrides.get(name)
        if v not in (None, ""):
            kw[name] = float(v) * f["si"]
    return OperatingPoint(
        p_amb=base.p_amb, p_air_supply=base.p_air_supply,
        p_fuel_supply=base.p_fuel_supply, fuel=base.fuel,
        T_air_in=kw.get("T_air_in", base.T_air_in),
        T_fuel_in=kw.get("T_fuel_in", base.T_fuel_in),
        mdot_air_max=kw.get("mdot_air_max", base.mdot_air_max),
        cd_injector_ox=kw.get("cd_injector_ox", base.cd_injector_ox),
        cd_injector_fuel=kw.get("cd_injector_fuel", base.cd_injector_fuel),
        burn_time=kw.get("burn_time", base.burn_time),
    )


def schema() -> dict[str, Any]:
    x = load_design_vector()
    op = load_operating_point()
    mat = load_material()
    design = []
    for f in DESIGN_FIELDS:
        lo, hi = DESIGN_BOUNDS[f["name"]]
        design.append({**{k: v for k, v in f.items() if k != "si"},
                       "min": round(lo / f["si"], 6), "max": round(hi / f["si"], 6),
                       "value": round(x.values[f["name"]] / f["si"], 6),
                       "integer": f["name"] in INTEGER_PARAMETERS})
    operating = []
    for f in OPERATING_FIELDS:
        v = getattr(op, f["name"])
        operating.append({**{k: v2 for k, v2 in f.items() if k != "si"},
                          "value": None if v is None else round(v / f["si"], 6),
                          "missing": v is None})
    return {
        "design": design,
        "operating": operating,
        "constants": {
            "p_amb_bar": op.p_amb / BAR,
            "p_air_supply_bar": op.p_air_supply / BAR,
            "p_fuel_supply_bar": op.p_fuel_supply / BAR,
            "fuel": dict(op.fuel.composition),
            "thermo_source": op.fuel.thermo_source,
        },
        "material_is_placeholder": mat["density_kg_m3"] is None,
        "min_feature_size_m": mat["process"]["min_feature_size_m"],
        "code_rev": git_revision(),
        "runs_root": str(runs_root()),
    }


# --------------------------------------------------------------------------- #
# Valutazione
# --------------------------------------------------------------------------- #
def _section(params) -> dict[str, Any]:
    """Sezione meridiana per il disegno.

    NON e' una ricostruzione approssimata a scopo grafico: e' esattamente la
    poligonale che viene rivoluzionata per generare il solido. Se lo spaccato a
    schermo e' sbagliato, lo e' anche il pezzo che esce dalla stampante.
    """
    d = params.derived
    return {
        "profile_mm": [[round(p[0], 4), round(p[2], 4)] for p in wall_profile_mm(params)],
        "stations_mm": {
            "L_c": d["L_c"] * 1e3, "L_conv": d["L_conv"] * 1e3,
            "x_lip": (d["L_c"] + d["L_conv"]) * 1e3,
            "R_c": d["R_c"] * 1e3, "R_lip": d["R_lip"] * 1e3,
            "r_centerbody": d["r_centerbody"] * 1e3,
            "t_wall": d["t_wall"] * 1e3, "t_face": d["t_face"] * 1e3,
            "x_throat": (d["L_c"] + d["L_conv"] + d["x_throat_plug"]) * 1e3,
            "x_tip": (d["L_c"] + d["L_conv"] + d["x_tip_full"]) * 1e3,
            "D_t_eq": d["D_t_eq"] * 1e3,
        },
        "injectors": {
            "N": int(d["N_inj"]), "R_inj_mm": d["R_inj"] * 1e3,
            "d_ox_mm": d["d_ox"] * 1e3, "d_fuel_mm": d["d_fuel"] * 1e3,
            "d_film_mm": d["d_film"] * 1e3, "R_film_mm": d["R_film"] * 1e3,
            "swirl_deg": math.degrees(params.free["theta_swirl"]),
        },
    }


def _persist(run_id: str, x: DesignVector, l0, obj) -> bool:
    try:
        env = env_fingerprint()
        with RunStore(runs_root() / "runs.db") as store:
            store.insert(run_id, x, l0, datetime.now(timezone.utc).isoformat(),
                         git_revision(),
                         blake2b(dumps(env, sort_keys=True).encode(), digest_size=8).hexdigest(),
                         objectives=obj, replace=True)
        return True
    except Exception:  # noqa: BLE001
        log.exception("run %s non salvata nel database", run_id)
        return False


def evaluate(design: dict[str, float], operating: dict[str, float | None]) -> dict[str, Any]:
    """Valutazione L0 completa. Costa ~50 ms: si serve nella richiesta."""
    op = operating_from(operating)
    x = DesignVector(values=design_to_si(design), bounds=DESIGN_BOUNDS)
    params, l0 = derive(x, op)
    mat = load_material()
    min_feature = mat["process"]["min_feature_size_m"]
    run_id = make_run_id(x, op)
    obj = objectives_l0(run_id, l0, op, x.values["p_c"],
                        min_feature_size=min_feature, derived=params.derived)
    return {
        "run_id": run_id,
        "stored": _persist(run_id, x, l0, obj),
        "l0": l0.to_dict(),
        "derived": dict(params.derived),
        "section": _section(params),
        "objectives": {"f": dict(obj.f), "g": dict(obj.g), "feasible": obj.feasible},
        "manufacturability": list(check_manufacturability(params, min_feature)),
        "warnings": list(l0.warnings),
        "assumptions": list(l0.assumptions),
    }


def build_geometry(design: dict[str, float], operating: dict[str, float | None],
                   sector: bool = False,
                   progress: Callable[[float, str], None] | None = None) -> dict[str, Any]:
    """Costruisce il solido ed esporta STEP e STL. Costa ~2 s: va in coda."""
    def step(f: float, m: str = "") -> None:
        if progress:
            progress(f, m)

    step(0.05, "valutazione termochimica")
    op = operating_from(operating)
    x = DesignVector(values=design_to_si(design), bounds=DESIGN_BOUNDS)
    params, l0 = derive(x, op)
    run_id = make_run_id(x, op)

    step(0.25, "costruzione del solido")
    art = build_and_export(params, runs_root() / run_id, run_id,
                           BuildOptions(sector=sector))
    step(0.9, "verifica di tenuta")
    (runs_root() / run_id / "geometry_params.json").write_text(
        params.to_json(), encoding="utf-8")
    step(1.0, "fatto")
    return {
        "run_id": run_id,
        "stl_url": f"/api/v1/geometry/{run_id}/stl",
        "step_url": f"/api/v1/geometry/{run_id}/step",
        "volume_mm3": art.volume * 1e9,
        "wetted_area_mm2": art.wetted_area * 1e6,
        "mass_g": art.volume * RHO_316L_NOMINAL * 1e3,
        "mass_is_placeholder": True,
        "valid_brep": art.is_valid_brep,
        "watertight": art.is_watertight_mesh,
        "n_triangles": art.n_triangles,
        "sha256_stl": art.sha256_stl,
        "sector": sector,
    }


def thermal(design: dict[str, float], operating: dict[str, float | None],
            k_wall: float, rho_wall: float, cp_wall: float,
            water_cooled: bool, h_coolant: float, T_coolant: float) -> dict[str, Any]:
    """Carico termico e risposta della parete, con e senza raffreddamento."""
    op = operating_from(operating)
    x = DesignVector(values=design_to_si(design), bounds=DESIGN_BOUNDS)
    params, l0 = derive(x, op)
    op.require("burn_time")

    model = MixtureModel.from_fuel(op.fuel)
    st = chamber_state(model, x.values["p_c"], op.T_air_in, op.T_fuel_in,
                       l0.mdot_air, l0.mdot_fuel_core)
    g = model.gas
    g.TPX = st.T, st.p, st.X
    g.equilibrate("TP")
    mu, cpg = g.viscosity, g.cp_mass
    Pr = mu * cpg / g.thermal_conductivity
    D_t = 2.0 * math.sqrt(l0.A_t / math.pi)
    A_c = math.pi * (params.derived["R_c"] ** 2 - params.derived["r_centerbody"] ** 2)
    t_end = float(op.burn_time)
    t_axis: list[float] = []

    stations = []
    for label, ar, M in (("camera", A_c / l0.A_t, 0.1), ("gola", 1.0, 1.0)):
        T_aw = adiabatic_wall_temperature(st.T, M, st.gamma, Pr)
        h = bartz_h_gas(D_t, x.values["p_c"], l0.c_star, mu, cpg, Pr, st.gamma,
                        ar, M, 800.0, st.T)
        sink = transient_slab(x.values["t_wall"], rho_wall, cp_wall, k_wall,
                              h, T_aw, t_end, n_steps=400)
        if not t_axis:
            t_axis = [round(v, 4) for v in sink.t.tolist()]
        entry: dict[str, Any] = {
            "label": label, "h_gas": h, "T_aw": T_aw,
            "q0_MW_m2": h * (T_aw - 293.15) / 1e6,
            "biot": sink.biot, "penetration_mm": sink.penetration_depth * 1e3,
            "sink": [round(v, 1) for v in sink.T_hot.tolist()],
        }
        if water_cooled:
            cooled = transient_slab(x.values["t_wall"], rho_wall, cp_wall, k_wall,
                                    h, T_aw, t_end, n_steps=400,
                                    h_back=h_coolant, T_back=T_coolant)
            entry["cooled"] = [round(v, 1) for v in cooled.T_hot.tolist()]
            entry["cooled_back"] = [round(v, 1) for v in cooled.T_cold.tolist()]
        stations.append(entry)

    return {
        "t": t_axis, "stations": stations, "T_ad": st.T,
        "burn_time": t_end,
        "material_is_placeholder": load_material()["density_kg_m3"] is None,
        "correlation": "Bartz 1957 — empirica, ±30 %, tarata su gole molto piu' grandi",
    }


def sweep(n: int, seed: int, operating: dict[str, float | None],
          progress: Callable[[float, str], None] | None = None) -> dict[str, Any]:
    """Piano sperimentale Latin Hypercube su L0, salvato nel database."""
    op = operating_from(operating)
    mat = load_material()
    min_feature = mat["process"]["min_feature_size_m"]
    campione = latin_hypercube(n, DESIGN_BOUNDS, seed, INTEGER_PARAMETERS)
    rev = git_revision()
    env = env_fingerprint()
    env_hash = blake2b(dumps(env, sort_keys=True).encode(), digest_size=8).hexdigest()
    now = datetime.now(timezone.utc).isoformat()

    ok = scartati = 0
    with RunStore(runs_root() / "runs.db") as store:
        for i, valori in enumerate(campione):
            try:
                x = DesignVector(values=valori, bounds=DESIGN_BOUNDS)
                params, l0 = derive(x, op)
                run_id = make_run_id(x, op, rev)
                obj = objectives_l0(run_id, l0, op, valori["p_c"],
                                    min_feature_size=min_feature, derived=params.derived)
                store.insert(run_id, x, l0, now, rev, env_hash, objectives=obj, replace=True)
                ok += 1
            except (ZefiroError, ValueError, NotImplementedError):
                # un candidato impossibile non e' un errore del piano: e' un dato
                scartati += 1
            if progress and (i % 10 == 0 or i == len(campione) - 1):
                progress((i + 1) / len(campione), f"{i+1}/{len(campione)} · {ok} validi")
        fattibili = store.query("SELECT COUNT(*) c FROM runs WHERE feasible = 1")[0]["c"]
        totale = store.count("runs")
    return {"valutati": ok, "scartati": scartati, "seed": seed,
            "fattibili_totali": fattibili, "run_totali": totale}


def plant(power_hp: float, tank_l: float, p_max_bar: float, p_min_bar: float,
          burn_s: float, T_amb_C: float) -> dict[str, Any]:
    T = T_amb_C + 273.15
    b = compressor_bounds(power_hp * HP_TO_W, 101325.0, p_max_bar * BAR, T)
    bd = tank_blowdown(tank_l * 1e-3, p_max_bar * BAR, p_min_bar * BAR, T)
    curva = []
    for mdot in (0.010, 0.020, 0.040, 0.060, 0.080, 0.100, 0.150, 0.200):
        ta, ti = burst_duration(bd, mdot)
        curva.append({"mdot_g_s": mdot * 1e3,
                      "adiabatico_s": None if math.isinf(ta) else round(ta, 2),
                      "isotermo_s": None if math.isinf(ti) else round(ti, 2)})
    return {
        "compressore": {
            "potenza_W": b.shaft_power,
            "mdot_isotermo_g_s": b.mdot_isothermal * 1e3,
            "mdot_adiabatico_g_s": b.mdot_adiabatic * 1e3,
        },
        "serbatoio": {
            "massa_iniziale_g": bd.mass_initial * 1e3,
            "utilizzabile_adiabatico_g": bd.mass_usable_adiabatic * 1e3,
            "utilizzabile_isotermo_g": bd.mass_usable_isothermal * 1e3,
            "T_finale_C": bd.T_final_adiabatic - 273.15,
        },
        "per_durata": {"secondi": burn_s,
                       "mdot_g_s": bd.mass_usable_adiabatic / burn_s * 1e3},
        "curva": curva,
    }


def runs(limit: int = 500) -> dict[str, Any]:
    db = runs_root() / "runs.db"
    if not db.exists():
        return {"rows": [], "total": 0, "feasible": 0}
    with RunStore(db) as store:
        rows = store.query(
            "SELECT run_id, thrust, Isp_s, T_ad, p_c, phi_core, f_film, epsilon, "
            "c_star, feasible, created_utc FROM runs ORDER BY created_utc DESC LIMIT ?",
            (limit,))
        fattibili = store.query("SELECT COUNT(*) c FROM runs WHERE feasible = 1")[0]["c"]
        return {"rows": [dict(r) for r in rows], "total": store.count("runs"),
                "feasible": fattibili}
