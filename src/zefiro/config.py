"""Caricamento e validazione della configurazione YAML.

Unico punto del codice in cui esistono unita' non SI: qui si legge `p_c_bar` e
si scrive `p_c`. Da qui in poi tutto e' SI.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from zefiro.schemas import DesignVector, FuelSpec, MissingDatum, OperatingPoint
from zefiro.units import BAR, DEG, MM

CONFIG_DIR = Path(__file__).resolve().parents[2] / "config"


def _load_yaml(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)
    if not isinstance(data, dict):
        raise MissingDatum(f"{path} non contiene una mappa YAML.")
    return data


def load_operating_point(path: Path | None = None) -> OperatingPoint:
    """Legge il punto operativo. I `null` restano `None` e falliranno all'uso.

    Deliberatamente NON si sostituiscono i null con default: un default fisico
    inventato e' peggio di un errore, perche' produce numeri credibili e falsi.
    """
    path = path or CONFIG_DIR / "operating_point.yaml"
    d = _load_yaml(path)
    f = d["fuel"]
    fuel = FuelSpec(
        composition=dict(f["composition_mole"]),
        phase_at_injection=f["phase_at_injection"],
        thermo_source=f["thermo_source"],
    )
    p_fuel, p_air = _supply(d, fuel)
    return OperatingPoint(
        p_amb=float(d["p_amb_bar"]) * BAR,
        p_air_supply=p_air,
        p_fuel_supply=p_fuel,
        fuel=fuel,
        T_air_in=_opt(d.get("T_air_in_K")),
        T_fuel_in=_opt(d.get("T_fuel_in_K")),
        mdot_air_max=_opt(d.get("mdot_air_max_kg_s")),
        cd_injector_ox=_opt(d.get("cd_injector_ox")),
        cd_injector_fuel=_opt(d.get("cd_injector_fuel")),
        burn_time=_opt(d.get("burn_time_s")),
    )


def load_design_vector(path: Path | None = None) -> DesignVector:
    from zefiro.geometry.parameters import DESIGN_BOUNDS

    path = path or CONFIG_DIR / "design_default.yaml"
    d = _load_yaml(path)["free"]
    values = {
        "p_c": float(d["p_c_bar"]) * BAR,
        "phi_core": float(d["phi_core"]),
        "f_film": float(d["f_film"]),
        "Dc_over_Dt": float(d["Dc_over_Dt"]),
        "Lc_over_Dc": float(d["Lc_over_Dc"]),
        "N_inj": float(d["N_inj"]),
        "theta_swirl": float(d["theta_swirl_deg"]) * DEG,
        "d_ox_ratio": float(d["d_ox_ratio"]),
        "fuel_vel_ratio": float(d["fuel_vel_ratio"]),
        "conv_half_angle": float(d["conv_half_angle_deg"]) * DEG,
        "plug_trunc": float(d["plug_trunc"]),
        "t_wall": float(d["t_wall_mm"]) * MM,
    }
    return DesignVector(values=values, bounds=DESIGN_BOUNDS)


def load_material(name: str = "aisi316l") -> dict[str, Any]:
    """Proprieta' del materiale. I `null` sono TODO n.6 e restano `null`."""
    return _load_yaml(CONFIG_DIR / "materials" / f"{name}.yaml")


def _opt(v: Any) -> float | None:
    return None if v is None else float(v)


DERIVED = "derived"


def _supply(d: dict[str, Any], fuel: FuelSpec) -> tuple[float, float]:
    """Risolve `p_fuel_supply_bar` e `p_air_supply_bar`.

    Il valore ammesso e' `derived`: le due pressioni si calcolano dalla
    temperatura di progetto della bombola (vedi `zefiro.feed.supply_pressures`
    per il perche' fisico di ciascuna).

    Un numero esplicito resta ammesso perche' una misura di targa e' un dato
    legittimo, ma viene **rifiutato se supera p_sat alla temperatura di
    progetto**: quella e' una pressione che la termodinamica non produce, e un
    file di configurazione non deve poterla esprimere. E' esattamente l'errore
    che questa funzione esiste per rendere impossibile.
    """
    from zefiro.feed import saturation_pressure, supply_pressures

    bottle = d["fuel"].get("bottle") or {}
    T_design = bottle.get("T_design_K")
    if T_design is None:
        raise MissingDatum(
            "fuel.bottle.T_design_K assente. E' la temperatura MINIMA a cui il "
            "motore verra' acceso: fissa p_sat, quindi il tetto di p_c, quindi "
            "l'intera geometria. Senza, non c'e' punto operativo."
        )
    p_tank_max = float(d["plant"]["air_tank"]["p_max_bar"]) * BAR
    sup = supply_pressures(dict(fuel.composition), float(T_design), p_tank_max,
                           _DP_FRACTION)

    p_sat = saturation_pressure(dict(fuel.composition), float(T_design))
    p_fuel = _pressione(d, "p_fuel_supply_bar", sup.p_fuel_supply)
    p_air = _pressione(d, "p_air_supply_bar", sup.p_air_supply)
    if p_fuel > p_sat * (1.0 + 1.0e-9):
        raise MissingDatum(
            f"p_fuel_supply_bar = {p_fuel/BAR:.2f} bar, ma a "
            f"{float(T_design)-273.15:.1f} C la bombola satura a {p_sat/BAR:.2f} bar. "
            "Nessuna bombola eroga piu' della propria tensione di vapore: o la "
            "temperatura di progetto e' piu' alta, o la miscela contiene piu' "
            "propano di quanto dichiarato (misura pressione E temperatura "
            "insieme e usa feed.composition_from_pressure), o il numero e' sbagliato."
        )
    if p_air > p_tank_max * (1.0 + 1.0e-9):
        raise MissingDatum(
            f"p_air_supply_bar = {p_air/BAR:.2f} bar supera il massimo del "
            f"serbatoio ({p_tank_max/BAR:.2f} bar)."
        )
    return p_fuel, p_air


#: Frazione minima di p_c richiesta come salto d'iniezione. Duplicata qui e non
#: importata da opt.objectives per non far dipendere il caricamento della
#: configurazione dall'ottimizzatore; il test test_config.py verifica che i due
#: valori coincidano.
_DP_FRACTION = 0.15


def _pressione(d: dict[str, Any], chiave: str, derivata: float) -> float:
    v = d[chiave]
    if isinstance(v, str):
        if v.strip().lower() != DERIVED:
            raise MissingDatum(f"{chiave}: atteso un numero o '{DERIVED}', letto {v!r}.")
        return derivata
    return float(v) * BAR
