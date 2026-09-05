"""Tempi caratteristici chimici: ritardo di ignizione e velocita' di fiamma laminare.

NON usa un valore tabellato di L*: il criterio di dimensionamento della camera
si costruisce confrontando tau_res (geometrico) con tau_chem (calcolato).
Questo richiede un meccanismo CINETICO valido per il GPL: gri30 non lo e'
(tarato su metano, non validato oltre C3). Senza il meccanismo esterno queste
funzioni sollevano MissingThermoData e L0Result.tau_chem resta None.
"""
from __future__ import annotations

import cantera as ct

from zefiro.l0.mixture import MixtureModel
from zefiro.schemas import MissingThermoData

# Meccanismi la cui CINETICA non e' valida per il GPL, anche se la loro
# termodinamica lo e'. Elencati per nome perche' l'errore va detto, non subito.
KINETICS_NOT_VALID_FOR_LPG = ("gri30.yaml", "gri30_highT.yaml", "gri30_ion.yaml")


def _check_kinetics(model: MixtureModel) -> None:
    src = model.fuel.thermo_source
    if src in KINETICS_NOT_VALID_FOR_LPG:
        raise MissingThermoData(
            f"{src!r} e' GRI-Mech 3.0: cinetica ottimizzata su metano e NON validata "
            "oltre C3. Usarla per il GPL darebbe tempi di ignizione privi di significato. "
            "Esegui scripts/fetch_mechanism.py per ottenere il San Diego mech."
        )


def ignition_delay(
    model: MixtureModel, T0: float, p: float, phi: float, t_max: float = 0.1
) -> float:
    """Ritardo di ignizione a volume costante, definito come istante di max dT/dt.

    Massimo gradiente e non soglia fissa: la soglia dipende dalle condizioni
    iniziali, il massimo no.
    """
    _check_kinetics(model)
    gas = model.gas
    gas.TP = T0, p
    gas.set_equivalence_ratio(phi, model.fuel_string, model.air_string, basis="mole")
    reactor = ct.IdealGasReactor(gas)
    net = ct.ReactorNet([reactor])
    t, t_prev, T_prev, best_slope, t_ign = 0.0, 0.0, reactor.T, 0.0, float("nan")
    while t < t_max:
        t = net.step()
        slope = (reactor.T - T_prev) / max(t - t_prev, 1.0e-12)
        if slope > best_slope:
            best_slope, t_ign = slope, t
        t_prev, T_prev = t, reactor.T
    return t_ign


def laminar_flame_speed(model: MixtureModel, T0: float, p: float, phi: float) -> float:
    """Velocita' di fiamma laminare S_L [m/s], fiamma 1D premiscelata libera."""
    _check_kinetics(model)
    gas = model.gas
    gas.TP = T0, p
    gas.set_equivalence_ratio(phi, model.fuel_string, model.air_string, basis="mole")
    flame = ct.FreeFlame(gas, width=0.03)
    flame.set_refine_criteria(ratio=3.0, slope=0.07, curve=0.14)
    flame.solve(loglevel=0, auto=True)
    return float(flame.velocity[0])
