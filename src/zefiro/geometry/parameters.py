"""Separazione fra parametri LIBERI e DERIVATI, e la derivazione stessa.

Regola: un derivato non e' mai un ingresso. Se ti serve un derivato come
ingresso, la parametrizzazione e' sbagliata e va cambiata qui, non aggirata.

Ordine di calcolo (docs/architettura.md sezione 4.5): la termochimica L0 non
dipende dalla geometria, quindi il ciclo apparente si spezza e derive() e' una
funzione PURA, senza iterazioni ne' tolleranze.
"""
from __future__ import annotations

import dataclasses
import functools
import math

from zefiro.geometry.aerospike import area_ratio, plug_contour
from zefiro.geometry.profile import (
    meridian_polygon,
    revolved_volume,
    wetted_area_from_contour,
)
from zefiro.l0.cycle import evaluate_l0
from zefiro.schemas import DesignVector, GeometryParams, L0Result, OperatingPoint
from zefiro.units import DEG, R_UNIVERSAL

#: Bounds dei 12 parametri liberi. Unita' SI. Vedi docs/architettura.md 4.1
#: per la motivazione di ciascuno.
DESIGN_BOUNDS: dict[str, tuple[float, float]] = {
    "p_c":             (3.0e5, 6.0e5),      # Pa
    "phi_core":        (0.70, 1.15),        # -
    "f_film":          (0.00, 0.35),        # -
    "Dc_over_Dt":      (2.00, 5.00),        # -
    "Lc_over_Dc":      (0.80, 4.00),        # -
    "N_inj":           (6.0, 24.0),         # - (INTERO: vedi nota sotto)
    "theta_swirl":     (0.0, 60.0 * DEG),   # rad
    "d_ox_ratio":      (0.02, 0.12),        # - (d_ox / D_c)
    "fuel_vel_ratio":  (0.30, 3.00),        # - (u_fuel / u_ox)
    "conv_half_angle": (20.0 * DEG, 45.0 * DEG),
    "plug_trunc":      (0.20, 1.00),        # -
    "t_wall":          (0.8e-3, 4.0e-3),    # m
}

#: `N_inj` e' discreto. L'ottimizzatore deve trattarlo come intero (pymoo lo
#: supporta nativamente). Arrotondare un float a valle produrrebbe due vettori
#: x diversi con lo stesso run_id, che rompe il versionamento (sezione 5.1).
INTEGER_PARAMETERS: frozenset[str] = frozenset({"N_inj"})


def default_design_vector(**overrides: float) -> DesignVector:
    """Punto centrale dei bounds, con eventuali sovrascritture."""
    values = {k: 0.5 * (lo + hi) for k, (lo, hi) in DESIGN_BOUNDS.items()}
    values["N_inj"] = float(round(values["N_inj"]))
    values.update(overrides)
    return DesignVector(values=values, bounds=DESIGN_BOUNDS)


def _ideal_gas_density(p: float, T: float, MW_kg_per_kmol: float) -> float:
    return p * MW_kg_per_kmol / (R_UNIVERSAL * 1000.0 * T)


def derive(
    x: DesignVector,
    op: OperatingPoint,
    mdot_air: float | None = None,
) -> tuple[GeometryParams, L0Result]:
    """Da (x, punto operativo) a geometria completa. Deterministico e puro."""
    v = dict(x.values)
    for name in INTEGER_PARAMETERS:
        if v[name] != round(v[name]):
            raise ValueError(f"{name} deve essere intero, vale {v[name]!r}")

    # ---- 1. termochimica: non dipende dalla geometria -------------------- #
    l0 = evaluate_l0(
        op, p_c=v["p_c"], phi_core=v["phi_core"], f_film=v["f_film"], mdot_air=mdot_air
    )

    A_t = l0.A_t
    D_t_eq = math.sqrt(4.0 * A_t / math.pi)      # diametro di gola EQUIVALENTE
    gamma = l0.gamma_c

    # ---- 2. aerospike ----------------------------------------------------- #
    # eps coerente con gamma (non quello dell'equilibrio spostato di L0): e' la
    # condizione perche' il contorno chiuda esattamente sull'asse.
    eps_gamma = area_ratio(l0.M_e, gamma)
    R_lip = math.sqrt(A_t * eps_gamma / math.pi)
    contour = plug_contour(gamma, l0.M_e, R_lip, truncation=v["plug_trunc"])

    # ---- 3. camera anulare ------------------------------------------------ #
    R_c = 0.5 * v["Dc_over_Dt"] * D_t_eq
    if R_c <= R_lip:
        raise ValueError(
            f"Dc_over_Dt = {v['Dc_over_Dt']:.2f} da' R_c = {R_c*1e3:.2f} mm <= R_lip = "
            f"{R_lip*1e3:.2f} mm: la camera non contiene l'ugello. Alza Dc_over_Dt."
        )
    L_c = v["Lc_over_Dc"] * (2.0 * R_c)
    L_conv = (R_c - R_lip) / math.tan(v["conv_half_angle"])
    r_cb = contour.r_throat                       # raggio del corpo centrale in camera

    V_chamber = math.pi * (R_c**2 - r_cb**2) * L_c
    V_conv = math.pi / 3.0 * L_conv * (R_c**2 + R_c * R_lip + R_lip**2) \
        - math.pi * r_cb**2 * L_conv
    V_c = V_chamber + V_conv
    L_star = V_c / A_t

    # ---- 4. iniezione ----------------------------------------------------- #
    N = int(v["N_inj"])
    d_ox = v["d_ox_ratio"] * (2.0 * R_c)
    A_ox_tot = N * math.pi / 4.0 * d_ox**2

    MW_air = 28.9649           # kg/kmol, coerente con units.AIR_MOLE_FRACTIONS
    rho_air = _ideal_gas_density(v["p_c"], float(op.T_air_in), MW_air)
    u_ox = l0.mdot_air / (rho_air * A_ox_tot)

    if op.fuel.phase_at_injection != "gas":
        raise NotImplementedError(
            "phase_at_injection = 'liquid': serve un modello di flash/vaporizzazione "
            "a monte di L0 (TODO n.3 in docs/architettura.md sezione 8). "
            "Non esiste una densita' di iniezione da assumere."
        )
    MW_fuel = _fuel_molar_mass(op)
    rho_fuel = _ideal_gas_density(v["p_c"], float(op.T_fuel_in), MW_fuel)
    u_fuel = v["fuel_vel_ratio"] * u_ox
    A_f_tot = l0.mdot_fuel_core / (rho_fuel * u_fuel)
    d_f = math.sqrt(4.0 * A_f_tot / (N * math.pi))

    # rapporto delle quantita' di moto: e' IL numero che governa il mixing
    J = (rho_fuel * u_fuel**2) / (rho_air * u_ox**2)

    dp_ox = _injector_dp(l0.mdot_air, op.cd_injector_ox, A_ox_tot, rho_air)
    dp_f = _injector_dp(l0.mdot_fuel_core, op.cd_injector_fuel, A_f_tot, rho_fuel)

    # ---- film cooling: FORI DISCRETI, non una fessura anulare ------------ #
    # Criterio di velocita': il film entra alla stessa velocita' del gas di
    # camera (co-flusso). E' il criterio che minimizza lo shear all'interfaccia,
    # quindi massimizza la lunghezza su cui il film sopravvive. Non c'e' nulla
    # da tabellare: la velocita' di camera si calcola.
    rho_c = v["p_c"] * l0.MW_c / (R_UNIVERSAL * 1000.0 * l0.T_ad)
    A_annulus = math.pi * (R_c**2 - r_cb**2)
    mdot_tot = l0.mdot_air + l0.mdot_fuel_core + l0.mdot_fuel_film
    u_chamber = mdot_tot / (rho_c * A_annulus)

    if l0.mdot_fuel_film > 0.0:
        A_film_tot = l0.mdot_fuel_film / (rho_fuel * u_chamber)
        d_film = math.sqrt(4.0 * A_film_tot / (N * math.pi))
    else:
        A_film_tot, d_film = 0.0, 0.0
    # I fori stanno appena dentro la parete di camera, con un setto (land) fra
    # foro e parete pari a mezzo diametro: e' quel setto il punto debole per SLM.
    R_film = R_c - d_film
    film_land = 0.5 * d_film

    # Setto fra fori d'iniezione ADIACENTI, sull'arco del settore a R_inj.
    # Non e' la stessa cosa del diametro minimo: due fori possono essere
    # entrambi fabbricabili e non starci comunque uno accanto all'altro.
    # Trovato dal generatore di mesh su un punto del fronte di Pareto: 18 fori
    # d'aria da 2.62 mm su un arco di 4.11 mm lasciavano 0.23 mm di materiale.
    R_inj = 0.5 * (r_cb + R_c)
    arco_settore = R_inj * 2.0 * math.pi / N
    diametri = [d_ox, d_f] + ([d_film] if d_film > 0.0 else [])
    injector_land = (arco_settore - sum(diametri)) / (len(diametri) + 1)

    # Volume del pezzo, ESATTO: e' la rivoluzione del poligono meridiano, e
    # quel volume ha forma chiusa (geometry/profile.py). Non e' una stima di
    # parete sottile - quella sbagliava del 7 % perche' contava due volte il
    # materiale agli spigoli concavi - e non costa OCCT.
    area_bagnata = wetted_area_from_contour(R_c, R_lip, L_c, L_conv, r_cb,
                                            contour.x, contour.r)
    t_face = max(2.0 * v["t_wall"], 2.0e-3)
    derived = {
        "A_t": A_t,
        "wetted_area": area_bagnata,
        "D_t_eq": D_t_eq, "gamma_c": gamma,
        "M_e": l0.M_e, "epsilon_gamma": eps_gamma, "epsilon_l0": l0.epsilon,
        "epsilon_mismatch": abs(eps_gamma - l0.epsilon) / l0.epsilon,
        "R_lip": R_lip, "r_throat_plug": contour.r_throat,
        "x_throat_plug": contour.x_throat, "x_tip_full": contour.x_tip_full,
        "nu_e": contour.nu_e,
        "R_c": R_c, "L_c": L_c, "L_conv": L_conv, "r_centerbody": r_cb,
        "V_c": V_c, "L_star": L_star,
        "N_inj": float(N), "sector_angle": 2.0 * math.pi / N,
        "d_ox": d_ox, "d_fuel": d_f, "A_ox_tot": A_ox_tot, "A_fuel_tot": A_f_tot,
        "R_inj": R_inj, "injector_land": injector_land,
        "injector_pitch_arc": arco_settore,
        "u_ox": u_ox, "u_fuel": u_fuel, "rho_air_inj": rho_air, "rho_fuel_inj": rho_fuel,
        "momentum_flux_ratio_J": J,
        "rho_chamber": rho_c, "u_chamber": u_chamber,
        "A_film_tot": A_film_tot, "d_film": d_film,
        "R_film": R_film, "film_land": film_land,
        "t_wall": v["t_wall"], "t_face": t_face,
        "mdot_air": l0.mdot_air,
        "mdot_fuel_core": l0.mdot_fuel_core,
        "mdot_fuel_film": l0.mdot_fuel_film,
    }
    if dp_ox is not None:
        derived["dp_inj_ox"] = dp_ox
    if dp_f is not None:
        derived["dp_inj_fuel"] = dp_f

    wall_volume = revolved_volume(
        meridian_polygon(derived, contour.x, contour.r)
    )
    derived["wall_volume"] = wall_volume
    l0 = dataclasses.replace(l0, wall_volume=wall_volume,
                             tau_res=(rho_c * V_c / mdot_tot))

    params = GeometryParams(
        free=v,
        derived=derived,
        plug_contour_x=contour.x,
        plug_contour_r=contour.r,
    )
    return params, l0


def _injector_dp(mdot: float, cd: float | None, area: float, rho: float) -> float | None:
    """Dp = (mdot/(Cd A))^2 / (2 rho). None se Cd non e' noto: non si inventa."""
    if cd is None or area <= 0.0:
        return None
    return (mdot / (cd * area)) ** 2 / (2.0 * rho)


def _fuel_molar_mass(op: OperatingPoint) -> float:
    """Massa molare del combustibile [kg/kmol], dal meccanismo dichiarato."""
    return _molar_mass_cached(
        tuple(sorted(op.fuel.composition.items())), op.fuel.thermo_source
    )


@functools.lru_cache(maxsize=64)
def _molar_mass_cached(composition: tuple, mechanism: str) -> float:
    """Il risultato dipende solo da (composizione, meccanismo), che in una run
    non cambiano mai: si memorizza il NUMERO, non l'oggetto Cantera, cosi' non
    c'e' nessuno stato condiviso di cui preoccuparsi."""
    from zefiro.l0.mixture import solution

    gas = solution(mechanism)
    return sum(frac * gas.molecular_weights[gas.species_index(name)]
               for name, frac in composition)


def check_manufacturability(
    params: GeometryParams, min_feature_size: float | None
) -> tuple[str, ...]:
    """Vincoli di processo SLM sulle dimensioni caratteristiche.

    `min_feature_size` e' il TODO n.7 (docs/architettura.md sezione 8): finche'
    non lo fornisci, questa funzione NON restituisce "tutto ok", restituisce
    esplicitamente che il controllo non e' stato eseguito. La differenza fra
    "verificato" e "non verificabile" non va nascosta.
    """
    if min_feature_size is None:
        return (
            "min_feature_size non fornito: i vincoli di fabbricabilita' SLM NON sono "
            "stati verificati (TODO n.7).",
        )
    d = params.derived
    checks = {
        "diametro foro ossidante": d["d_ox"],
        "diametro foro combustibile": d["d_fuel"],
        "diametro foro film": d["d_film"],
        "setto fra foro film e parete": d["film_land"],
        "spessore di parete": d["t_wall"],
        "spessore del labbro": d["t_wall"],
    }
    return tuple(
        f"{name} = {value*1e3:.3f} mm < min_feature_size = {min_feature_size*1e3:.3f} mm"
        for name, value in checks.items()
        if 0.0 < value < min_feature_size
    )
