"""Espansione 1D isentropica e prestazione dell'ugello.

Metodo: nessuna formula chiusa con gamma costante. Si integra l'espansione
isentropica reale in Cantera, cercando la gola come **massimo del flusso di
massa specifico** rho*u lungo l'isentropica. Questo:
  * non richiede di scegliere un gamma "rappresentativo";
  * funziona identico per espansione congelata e in equilibrio spostato;
  * la condizione sonica esce come risultato, non come ipotesi.

Il metodo si rompe se l'espansione non e' isentropica (urti, attrito, scambio
termico), che e' proprio cio' che la CFD di L1 dovra' quantificare.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import cantera as ct

from zefiro.l0.equilibrium import ChamberState

FROZEN = "frozen"
SHIFTING = "shifting"

# Sotto questo rapporto di pressione l'ugello e' appena supersonico e il
# vantaggio di compensazione di quota dell'aerospike e' irrilevante.
AEROSPIKE_MEANINGFUL_PR = 8.0


@dataclass(frozen=True)
class NozzleSolution:
    c_star: float        # m/s
    A_t_over_mdot: float # m^2 / (kg/s)   -> A_t = mdot * questo
    epsilon: float       # A_e / A_t
    u_e: float           # m/s
    M_e: float
    p_e: float           # Pa
    T_e: float           # K
    C_F: float
    mode: str


def _state_at_pressure(gas: ct.Solution, s0: float, h0: float, p: float, mode: str):
    """Stato isentropico a pressione p; ritorna (rho*u, u, rho, T, a)."""
    gas.SP = s0, p
    if mode == SHIFTING:
        gas.equilibrate("SP")
    dh = h0 - gas.enthalpy_mass
    u = math.sqrt(2.0 * dh) if dh > 0.0 else 0.0
    return gas.density * u, u, gas.density, gas.T, gas.sound_speed


def _find_throat(gas: ct.Solution, s0: float, h0: float, p_c: float, mode: str) -> float:
    """Pressione di gola = argmax di rho*u. Sezione aurea, bracket ampio.

    Per gamma in [1.1, 1.4] il rapporto critico p*/p_c sta fra 0.53 e 0.58;
    il bracket [0.20, 0.95] lo contiene con ampio margine e la funzione e'
    unimodale su di esso.
    """
    lo, hi = 0.20 * p_c, 0.95 * p_c
    inv_phi = (math.sqrt(5.0) - 1.0) / 2.0
    a, b = lo, hi
    c = b - inv_phi * (b - a)
    d = a + inv_phi * (b - a)
    fc = _state_at_pressure(gas, s0, h0, c, mode)[0]
    fd = _state_at_pressure(gas, s0, h0, d, mode)[0]
    for _ in range(80):
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - inv_phi * (b - a)
            fc = _state_at_pressure(gas, s0, h0, c, mode)[0]
        else:
            a, c, fc = c, d, fd
            d = a + inv_phi * (b - a)
            fd = _state_at_pressure(gas, s0, h0, d, mode)[0]
        if abs(b - a) < 1.0e-8 * p_c:
            break
    return 0.5 * (a + b)


def expand_to_pressure(
    gas: ct.Solution,
    chamber: ChamberState,
    p_e: float,
    p_amb: float,
    mode: str = SHIFTING,
) -> NozzleSolution:
    """Espande dallo stato di camera fino a `p_e`.

    Per un aerospike al punto di progetto si passa p_e = p_amb (espansione
    adattata): il termine di pressione si annulla e C_F = u_e / c*.
    """
    if mode not in (FROZEN, SHIFTING):
        raise ValueError(f"mode non valido: {mode!r}")
    if p_e >= chamber.p:
        raise ValueError("p_e deve essere minore della pressione di camera.")

    gas.TPX = chamber.T, chamber.p, chamber.X
    if mode == SHIFTING:
        gas.equilibrate("TP")
    s0, h0 = gas.entropy_mass, gas.enthalpy_mass

    p_t = _find_throat(gas, s0, h0, chamber.p, mode)
    G_t = _state_at_pressure(gas, s0, h0, p_t, mode)[0]

    G_e, u_e, _rho_e, T_e, a_e = _state_at_pressure(gas, s0, h0, p_e, mode)

    c_star = chamber.p / G_t                     # = p_c A_t / mdot
    A_t_over_mdot = 1.0 / G_t
    epsilon = G_t / G_e
    A_e_over_mdot = 1.0 / G_e
    thrust_over_mdot = u_e + (p_e - p_amb) * A_e_over_mdot
    C_F = thrust_over_mdot / (chamber.p * A_t_over_mdot)
    return NozzleSolution(
        c_star=c_star,
        A_t_over_mdot=A_t_over_mdot,
        epsilon=epsilon,
        u_e=u_e,
        M_e=u_e / a_e,
        p_e=p_e,
        T_e=T_e,
        C_F=C_F,
        mode=mode,
    )
