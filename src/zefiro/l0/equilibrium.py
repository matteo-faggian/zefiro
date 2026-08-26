"""Stato di camera: miscelazione adiabatica dei due flussi + equilibrio HP.

Punto delicato: aria e GPL entrano a temperature DIVERSE. Miscelarli imponendo
una temperatura media pesata sarebbe sbagliato (l'entalpia non e' lineare in T
con cp variabile). Si usa `ct.Quantity`, la cui somma conserva massa ed
entalpia a pressione costante: e' esattamente il bilancio di un miscelatore
adiabatico.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import cantera as ct

from zefiro.l0.mixture import MixtureModel


@dataclass(frozen=True)
class ChamberState:
    T: float                     # K   (adiabatica di fiamma, = T di ristagno in camera)
    p: float                     # Pa
    X: Mapping[str, float]       # frazioni molari all'equilibrio
    gamma: float                 # cp/cv congelato della miscela di equilibrio
    MW: float                    # kg/kmol
    cp: float                    # J/(kg K)
    h: float                     # J/kg
    s: float                     # J/(kg K)


def chamber_state(
    model: MixtureModel,
    p_c: float,
    T_air_in: float,
    T_fuel_in: float,
    mdot_air: float,
    mdot_fuel: float,
) -> ChamberState:
    """Miscela adiabaticamente i due flussi a `p_c` e porta all'equilibrio.

    Ipotesi attive:
      * combustione completa fino all'equilibrio chimico (nessuna limitazione
        cinetica): e' un LIMITE SUPERIORE di T_ad e di c*;
      * velocita' in camera trascurabile, quindi T_statica = T_ristagno. Vale
        se il rapporto di contrazione e' abbastanza alto; per Dc/Dt >= 2.5
        l'errore su T e' sotto lo 0.5 %;
      * nessuna perdita di calore verso la parete.
    """
    gas = model.gas
    gas.TPX = T_air_in, p_c, model.air_string
    q_air = ct.Quantity(gas, mass=mdot_air, constant="HP")
    gas.TPX = T_fuel_in, p_c, model.fuel_string
    q_fuel = ct.Quantity(gas, mass=mdot_fuel, constant="HP")
    mix = q_air + q_fuel          # somma a H e p costanti = miscelatore adiabatico
    mix.equilibrate("HP")
    g = mix.phase
    X = {k: v for k, v in g.mole_fraction_dict().items() if v > 1.0e-6}
    return ChamberState(
        T=float(g.T),
        p=float(g.P),
        X=X,
        gamma=float(g.cp_mass / g.cv_mass),
        MW=float(g.mean_molecular_weight),
        cp=float(g.cp_mass),
        h=float(g.enthalpy_mass),
        s=float(g.entropy_mass),
    )
