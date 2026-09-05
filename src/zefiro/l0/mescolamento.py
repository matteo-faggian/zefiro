"""Che cosa costa, in spinta, una miscela non uniforme.

PERCHE' QUESTO MODULO ESISTE. Il modello L0 calcola la prestazione mettendo le
due portate in un miscelatore adiabatico e portando il risultato all'equilibrio
(`zefiro.l0.equilibrium.chamber_state`). Dentro quella riga c'e' l'ipotesi piu'
pesante di tutta la catena delle prestazioni - **il combustibile e l'aria sono
perfettamente mescolati** - e fino a questo modulo non era ne' dichiarata ne'
quantificata. Isp, c*, T_ad, e quindi le portate, l'area di gola e tutto il
dimensionamento, erano il limite di mescolamento perfetto.

La CFD di mescolamento misura una disuniformita'

    U = sqrt( < (Y - <Y>)^2 > ) / <Y>

sui piani a valle dei getti. Questo modulo traduce quella U in una perdita di
spinta, in modo che le due meta' del progetto - la fluidodinamica e la
prestazione - parlino dello stesso motore.

IL MODELLO, E PERCHE' QUESTO. Si divide il flusso in TUBI DI FLUSSO che non si
scambiano massa: ognuno brucia al proprio rapporto e si espande per conto suo,
e le spinte si sommano. E' il trattamento classico della maldistribuzione del
rapporto di miscela, e qui e' anche l'unico onesto: se il fluido si mescolasse
lungo l'ugello, la disuniformita' misurata all'ingresso della camera non
significherebbe niente.

La CHIUSURA A DUE PARCELLE. Di una distribuzione si conoscono solo i primi due
momenti, media e varianza: <Y> dalla portata, U dalla CFD. Serve una
distribuzione, e se ne sceglie una a DUE valori perche' e' l'unica che si scrive
senza inventare nient'altro:

    una frazione di portata f porta tutto il combustibile a Y_r,
    il resto e' aria pura.

Imponendo media e varianza si ricava, in forma chiusa e senza parametri liberi,

    Y_r = <Y> (1 + U^2)      f = <Y> / Y_r = 1 / (1 + U^2)

Nota che per U > 1 una distribuzione SIMMETRICA attorno alla media sarebbe
impossibile (darebbe frazioni di massa negative): U > 1 significa gia' di per
se' che gran parte della portata e' aria quasi pura e il combustibile sta tutto
in una minoranza di flusso molto ricca. Non e' un'ipotesi del modello, e' una
proprieta' del numero.

E' un LIMITE INFERIORE della prestazione a parita' di U: due parcelle sono la
segregazione massima compatibile con quella varianza. La realta' sta fra questa
curva e il mescolamento perfetto. Detto altrimenti: se questa curva dice che si
puo' vivere con la U misurata, ci si puo' vivere davvero.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import cantera as ct

from zefiro.l0.equilibrium import chamber_state
from zefiro.l0.mixture import MixtureModel
from zefiro.l0.nozzle import SHIFTING, expand_to_pressure


@dataclass(frozen=True)
class PerditaDiMescolamento:
    """Il conto di una disuniformita' `U` sulla spinta."""
    U: float
    Y_ricco: float          # frazione di massa di combustibile nella parcella ricca
    phi_ricco: float        # rapporto di equivalenza della parcella ricca
    frazione_di_portata: float   # f: quanta portata sta nella parcella ricca
    T_ricco: float          # K, temperatura di quella parcella
    spinta: float           # N
    spinta_mescolata: float  # N, riferimento a mescolamento perfetto
    Isp: float              # s

    @property
    def rendimento(self) -> float:
        return self.spinta / self.spinta_mescolata


def _tubo(model: MixtureModel, p_c: float, p_amb: float, T_air: float,
          T_fuel: float, m_air: float, m_fuel: float) -> tuple[float, float]:
    """(velocita' d'uscita, temperatura) di un tubo di flusso che brucia da solo."""
    if m_fuel <= 0.0:
        #: ARIA PURA: non si chiama l'equilibrio. Non e' un'ottimizzazione, e'
        #: correttezza - a 293 K si finirebbe sotto il limite inferiore delle
        #: tabelle NASA-7 (300 K) e Cantera avvisa. Aria fredda che si espande
        #: e' un'espansione isentropica e basta, senza chimica da risolvere.
        gas = model.gas
        gas.TPX = T_air, p_c, model.air_string
        s0, h0 = gas.entropy_mass, gas.enthalpy_mass
        gas.SP = s0, p_amb
        dh = h0 - gas.enthalpy_mass
        return (math.sqrt(2.0 * dh) if dh > 0.0 else 0.0), T_air
    ch = chamber_state(model, p_c, T_air, T_fuel, m_air, m_fuel)
    noz = expand_to_pressure(model.gas, ch, p_amb, p_amb, mode=SHIFTING)
    return noz.u_e, ch.T


def perdita_di_mescolamento(
    model: MixtureModel, U: float, p_c: float, p_amb: float,
    T_air: float, T_fuel: float, mdot_air: float, mdot_fuel: float,
) -> PerditaDiMescolamento:
    """Spinta con disuniformita' `U`, contro quella a mescolamento perfetto."""
    if U < 0.0:
        raise ValueError("U e' uno scarto relativo: non puo' essere negativa.")
    m_tot = mdot_air + mdot_fuel
    Y_med = mdot_fuel / m_tot
    u_id, T_id = _tubo(model, p_c, p_amb, T_air, T_fuel, mdot_air, mdot_fuel)
    F_id = m_tot * u_id
    if U == 0.0:
        return PerditaDiMescolamento(
            U=0.0, Y_ricco=Y_med, phi_ricco=model.phi_from_flows(mdot_air, mdot_fuel),
            frazione_di_portata=1.0, T_ricco=T_id, spinta=F_id,
            spinta_mescolata=F_id, Isp=u_id / 9.80665)

    Y_r = Y_med * (1.0 + U * U)
    if Y_r >= 1.0:
        raise ValueError(
            f"U = {U:.2f} porta la parcella ricca a Y = {Y_r:.3f} >= 1, cioe' a "
            "combustibile puro: la chiusura a due parcelle non si applica piu'. "
            "Serve una distribuzione con piu' di due valori, e per sceglierla "
            "servono piu' di due momenti - cioe' l'istogramma dalla CFD.")
    f = Y_med / Y_r
    m_r, m_l = f * m_tot, (1.0 - f) * m_tot
    mf_r, ma_r = m_r * Y_r, m_r * (1.0 - Y_r)
    u_r, T_r = _tubo(model, p_c, p_amb, T_air, T_fuel, ma_r, mf_r)
    u_l, _ = _tubo(model, p_c, p_amb, T_air, T_fuel, m_l, 0.0)
    F = m_r * u_r + m_l * u_l
    return PerditaDiMescolamento(
        U=U, Y_ricco=Y_r, phi_ricco=model.phi_from_flows(ma_r, mf_r),
        frazione_di_portata=f, T_ricco=T_r, spinta=F, spinta_mescolata=F_id,
        Isp=F / (m_tot * 9.80665))
