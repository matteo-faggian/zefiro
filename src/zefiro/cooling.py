"""Raffreddamento a liquido di una parete: dimensionamento del canale.

Contesto: su un banco a terra l'acqua di rete e' il refrigerante piu' semplice
disponibile, e il raffreddamento ad acqua e' la soluzione classica per i
combustori da prova. La domanda non e' pero' "quanta potenza devo togliere"
(quella e' facile), e' **se il flusso termico LOCALE riesce ad attraversare
l'interfaccia parete-acqua senza far collassare l'ebollizione**.

Il criterio di progetto adottato qui e' il piu' conservativo possibile:
mantenere la parete lato acqua **sotto la temperatura di saturazione**, cioe'
evitare del tutto l'ebollizione. Cosi' non serve alcuna correlazione di crisi
termica per il dimensionamento, e la CHF resta solo come verifica di margine.

Proprieta' dell'acqua da CoolProp (EOS di riferimento IAPWS-95).
Le correlazioni usate sono dichiarate una per una nel docstring: Dittus-Boelter
e Blasius sono EMPIRICHE, Zuber e' derivata da instabilita' idrodinamica.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

from CoolProp.CoolProp import PropsSI

G0 = 9.80665

DITTUS_BOELTER_REFERENCE = (
    "Dittus & Boelter, Univ. California Publ. Eng. 2:443, 1930. Valido per "
    "Re > 10^4, 0.6 < Pr < 160, L/D > 10, differenze di temperatura moderate."
)
ZUBER_REFERENCE = (
    "Zuber, N., 'On the Stability of Boiling Heat Transfer', Trans. ASME 80, "
    "1958. CHF in pool boiling saturo, derivata da instabilita' di "
    "Helmholtz del getto di vapore."
)


def saturation_temperature(p: float, fluid: str = "Water") -> float:
    return PropsSI("T", "P", p, "Q", 0, fluid)


def dittus_boelter_h(
    hydraulic_diameter: float,
    velocity: float,
    T_bulk: float,
    p: float,
    fluid: str = "Water",
    heating: bool = True,
) -> float:
    """Coefficiente di scambio lato liquido [W/(m^2 K)].

        Nu = 0.023 Re^0.8 Pr^n,   n = 0.4 se il fluido si scalda

    CORRELAZIONE EMPIRICA. Fuori dal suo campo di validita' (in particolare per
    Re < 10^4 o per differenze di temperatura grandi) l'errore puo' superare il
    25 %. Qui la si usa per dimensionare, non per verificare.
    """
    rho = PropsSI("D", "T", T_bulk, "P", p, fluid)
    mu = PropsSI("V", "T", T_bulk, "P", p, fluid)
    k = PropsSI("L", "T", T_bulk, "P", p, fluid)
    cp = PropsSI("C", "T", T_bulk, "P", p, fluid)
    re = rho * velocity * hydraulic_diameter / mu
    pr = mu * cp / k
    nu = 0.023 * re**0.8 * pr ** (0.4 if heating else 0.3)
    return nu * k / hydraulic_diameter


def reynolds(hydraulic_diameter: float, velocity: float, T: float, p: float,
             fluid: str = "Water") -> float:
    rho = PropsSI("D", "T", T, "P", p, fluid)
    mu = PropsSI("V", "T", T, "P", p, fluid)
    return rho * velocity * hydraulic_diameter / mu


def channel_pressure_drop(
    hydraulic_diameter: float, length: float, velocity: float, T: float, p: float,
    fluid: str = "Water",
) -> float:
    """Perdita di carico distribuita [Pa], attrito di Blasius per tubo liscio.

        f = 0.316 Re^-0.25   (CORRELAZIONE EMPIRICA, Re < 10^5)
        dp = f (L/D) rho V^2 / 2

    Non include le perdite concentrate (imbocco, curve, collettore), che su un
    canale corto e tortuoso possono valere quanto il distribuito.
    """
    rho = PropsSI("D", "T", T, "P", p, fluid)
    re = reynolds(hydraulic_diameter, velocity, T, p, fluid)
    f = 0.316 * re**-0.25
    return f * (length / hydraulic_diameter) * rho * velocity**2 / 2.0


def zuber_chf(p: float, fluid: str = "Water") -> float:
    """Flusso termico critico in pool boiling saturo [W/m^2] (Zuber).

        q_CHF = 0.131 h_fg rho_v^0.5 [sigma g (rho_l - rho_v)]^0.25

    E' un **pavimento**, non il limite reale del caso in esame: con convezione
    forzata e liquido sottoraffreddato la CHF sale di quasi un ordine di
    grandezza. Serve a capire quanto si e' lontani dal caso peggiore.
    """
    T_sat = saturation_temperature(p, fluid)
    rho_l = PropsSI("D", "T", T_sat, "Q", 0, fluid)
    rho_v = PropsSI("D", "T", T_sat, "Q", 1, fluid)
    h_fg = PropsSI("H", "T", T_sat, "Q", 1, fluid) - PropsSI("H", "T", T_sat, "Q", 0, fluid)
    sigma = PropsSI("I", "T", T_sat, "Q", 0, fluid)
    return 0.131 * h_fg * math.sqrt(rho_v) * (sigma * G0 * (rho_l - rho_v)) ** 0.25


def coolant_mass_flow(heat_load: float, delta_T: float, T_in: float, p: float,
                      fluid: str = "Water") -> float:
    """Portata [kg/s] per assorbire `heat_load` con un salto `delta_T`."""
    cp = PropsSI("C", "T", T_in + 0.5 * delta_T, "P", p, fluid)
    return heat_load / (cp * delta_T)


@dataclass(frozen=True)
class WallTemperatures:
    T_gas_side: float
    T_coolant_side: float
    T_sat: float
    margin_to_boiling: float     # K; positivo = niente ebollizione
    conduction_drop: float       # K attraverso la parete
    film_drop: float             # K attraverso il film liquido


def steady_wall_temperatures(
    q: float, thickness: float, k_wall: float, h_coolant: float,
    T_coolant: float, p_coolant: float, fluid: str = "Water",
) -> WallTemperatures:
    """Parete piana in regime STAZIONARIO, raffreddata sul dorso.

        T_lato_acqua = T_acqua + q/h
        T_lato_gas   = T_lato_acqua + q t / k

    Regime stazionario e non transitorio: con l'acqua la parete raggiunge
    l'equilibrio in frazioni di secondo, quindi il transitorio non e' il caso
    dimensionante. E' l'opposto della parete a pozzo termico, dove il
    transitorio E' il problema.
    """
    film = q / h_coolant
    cond = q * thickness / k_wall
    T_w = T_coolant + film
    T_sat = saturation_temperature(p_coolant, fluid)
    return WallTemperatures(
        T_gas_side=T_w + cond,
        T_coolant_side=T_w,
        T_sat=T_sat,
        margin_to_boiling=T_sat - T_w,
        conduction_drop=cond,
        film_drop=film,
    )


def required_h_for_no_boiling(
    q: float, T_coolant: float, p_coolant: float, safety_margin_K: float = 20.0,
    fluid: str = "Water",
) -> float:
    """h minimo perche' la parete lato acqua resti sotto saturazione."""
    dT = saturation_temperature(p_coolant, fluid) - T_coolant - safety_margin_K
    if dT <= 0.0:
        raise ValueError(
            "L'acqua entra troppo vicina alla saturazione: nessun margine."
        )
    return q / dT
