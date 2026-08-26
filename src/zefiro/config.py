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
    return OperatingPoint(
        p_amb=float(d["p_amb_bar"]) * BAR,
        p_air_supply=float(d["p_air_supply_bar"]) * BAR,
        p_fuel_supply=float(d["p_fuel_supply_bar"]) * BAR,
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
