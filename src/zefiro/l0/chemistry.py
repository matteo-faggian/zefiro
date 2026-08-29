"""Tempo chimico di camera: blowout di un reattore perfettamente miscelato (PSR).

PERCHE' NON IL RITARDO DI AUTOACCENSIONE
----------------------------------------
La grandezza che si usa di solito per "quanto ci mette a bruciare" e' il
ritardo di autoaccensione a volume/pressione costante. Qui sarebbe la
grandezza SBAGLIATA, e non di poco: i reagenti entrano in camera a ~300 K, e
a 300 K una miscela propano/aria non si autoaccende in alcun tempo
significativo. Il ritardo calcolato sarebbe enorme e la conclusione ("la
camera non brucia") falsa.

In una camera di combustione l'accensione non e' un'autoaccensione: e'
sostenuta dal ricircolo dei prodotti caldi. La domanda corretta e' quindi
"dato che l'accensione c'e', quanto deve restare il gas per bruciare?", che e'
esattamente la definizione del **tempo di blowout di un PSR**: il minimo tempo
di residenza per cui esiste ancora una soluzione stazionaria accesa.

Sotto tau_blowout la reazione si spegne per diluizione con reagenti freddi;
sopra, la conversione e' completa. E' un tempo caratteristico della sola
CINETICA, calcolato dal meccanismo dichiarato: non c'e' nessun numero
tabellato in questo modulo.

COSTO
-----
Una bisezione costa ~2.5 s: troppo per il ciclo interno di un ottimizzatore
(migliaia di valutazioni). Per questo il modulo espone anche una
tabulazione su griglia (p_c, phi) con interpolazione logaritmica, il cui
errore viene MISURATO dal test, non assunto.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import cantera as ct
import numpy as np

from zefiro.l0.mixture import MixtureModel

#: Frazione dell'innalzamento di temperatura di equilibrio sotto la quale il
#: reattore e' considerato SPENTO. Il salto fra ramo acceso e ramo spento a
#: cavallo del blowout e' netto (centinaia di K), quindi il risultato e'
#: insensibile a questa soglia in un intervallo largo: il test lo verifica.
EXTINCTION_FRACTION = 0.5

#: Tolleranza relativa della bisezione sul tempo di blowout.
BISECTION_RTOL = 0.02

#: Tempo di residenza di partenza [s]: se il PSR non e' acceso qui, non lo e'
#: da nessuna parte (e' ~30 volte il piu' lungo blowout trovato nel box).
TAU_UPPER = 1.0e-2


@dataclass(frozen=True)
class PSRPoint:
    tau_blowout: float | None   # s ; None = non si accende nemmeno a TAU_UPPER
    T_inlet: float              # K ; miscela fredda, entalpia conservata
    T_equilibrium: float        # K


def _inlet_state(
    model: MixtureModel, p: float, phi: float, T_air: float, T_fuel: float
) -> tuple[float, np.ndarray]:
    """Reagenti MISCELATI ma non bruciati, a entalpia ed elementi conservati."""
    g = model.gas
    afr = model.afr_stoichiometric()
    g.TPX = T_air, p, model.air_string
    q_air = ct.Quantity(g, mass=1.0, constant="HP")
    g.TPX = T_fuel, p, model.fuel_string
    q_fuel = ct.Quantity(g, mass=phi / afr, constant="HP")
    mix = q_air + q_fuel
    return float(mix.phase.T), np.array(mix.phase.Y, dtype=float)


def _psr_temperature(mech: str, p: float, T_in: float, Y_in: np.ndarray, tau: float) -> float:
    """Temperatura stazionaria del PSR, innescato dal ramo acceso.

    L'inizializzazione all'equilibrio non e' un trucco: il PSR ha due rami
    stazionari (acceso e spento) e noi cerchiamo il limite di ESISTENZA di
    quello acceso. Partire dall'equilibrio e' il modo di trovarlo se c'e'.
    """
    g = ct.Solution(mech)
    g.TPY = T_in, p, Y_in
    upstream = ct.Reservoir(g, clone=True)
    g.equilibrate("HP")
    reactor = ct.IdealGasReactor(g, energy="on", clone=True)
    downstream = ct.Reservoir(g, clone=True)
    mfc = ct.MassFlowController(upstream, reactor, mdot=reactor.mass / tau)
    ct.PressureController(reactor, downstream, primary=mfc, K=1.0e-5)
    net = ct.ReactorNet([reactor])
    net.advance_to_steady_state()
    return float(reactor.T)


def psr_blowout(
    model: MixtureModel, p: float, phi: float, T_air: float, T_fuel: float
) -> PSRPoint:
    """Tempo di blowout per bisezione geometrica sul tempo di residenza."""
    mech = model.fuel.thermo_source
    T_in, Y_in = _inlet_state(model, p, phi, T_air, T_fuel)
    g = ct.Solution(mech)
    g.TPY = T_in, p, Y_in
    g.equilibrate("HP")
    T_eq = float(g.T)
    threshold = T_in + EXTINCTION_FRACTION * (T_eq - T_in)

    def lit(tau: float) -> bool:
        return _psr_temperature(mech, p, T_in, Y_in, tau) > threshold

    hi = TAU_UPPER
    if not lit(hi):
        return PSRPoint(None, T_in, T_eq)
    lo = hi
    for _ in range(40):
        lo *= 0.5
        if not lit(lo):
            break
    else:                                   # pragma: no cover - fuori scala fisica
        return PSRPoint(lo, T_in, T_eq)
    while hi / lo > 1.0 + BISECTION_RTOL:
        mid = math.sqrt(lo * hi)
        if lit(mid):
            hi = mid
        else:
            lo = mid
    return PSRPoint(hi, T_in, T_eq)


@dataclass(frozen=True)
class BlowoutTable:
    """Tabulazione di tau_blowout(p, phi), con interpolazione bilineare in log.

    In log perche' tau varia di un fattore ~3 sul box e la dipendenza da p e'
    prossima a una potenza: in log la funzione e' quasi un piano, e l'errore
    di interpolazione crolla. Il test lo misura contro il calcolo esatto.
    """

    p_grid: tuple[float, ...]
    phi_grid: tuple[float, ...]
    log_tau: tuple[tuple[float, ...], ...]
    T_air: float
    T_fuel: float
    mech: str

    @classmethod
    def build(
        cls,
        model: MixtureModel,
        p_grid,
        phi_grid,
        T_air: float,
        T_fuel: float,
    ) -> "BlowoutTable":
        rows = []
        for p in p_grid:
            row = []
            for phi in phi_grid:
                pt = psr_blowout(model, p, phi, T_air, T_fuel)
                if pt.tau_blowout is None:
                    raise ValueError(
                        f"nessuna soluzione accesa a p = {p/1e5:.2f} bar, phi = {phi:.3f}: "
                        "la tabella non copre un punto spento e non deve inventarne uno"
                    )
                row.append(math.log(pt.tau_blowout))
            rows.append(tuple(row))
        return cls(tuple(p_grid), tuple(phi_grid), tuple(rows), T_air, T_fuel,
                   model.fuel.thermo_source)

    def __call__(self, p: float, phi: float) -> float:
        """tau_blowout interpolato [s]. Fuori griglia: errore, non estrapolazione."""
        p_lo = self.p_grid[0]
        p_hi = self.p_grid[-1]
        f_lo = self.phi_grid[0]
        f_hi = self.phi_grid[-1]
        if not (p_lo <= p <= p_hi and f_lo <= phi <= f_hi):
            raise ValueError(
                f"(p = {p/1e5:.3f} bar, phi = {phi:.3f}) fuori dalla griglia tabulata "
                f"[{p_lo/1e5:.2f}, {p_hi/1e5:.2f}] bar x [{f_lo:.2f}, {f_hi:.2f}]. "
                "Estrapolare un tempo chimico non e' ammesso: ricostruisci la tabella."
            )
        i = min(np.searchsorted(self.p_grid, p) - 1, len(self.p_grid) - 2)
        j = min(np.searchsorted(self.phi_grid, phi) - 1, len(self.phi_grid) - 2)
        i = max(i, 0)
        j = max(j, 0)
        p0, p1 = self.p_grid[i], self.p_grid[i + 1]
        f0, f1 = self.phi_grid[j], self.phi_grid[j + 1]
        s = 0.0 if p1 == p0 else (p - p0) / (p1 - p0)
        t = 0.0 if f1 == f0 else (phi - f0) / (f1 - f0)
        L = self.log_tau
        val = ((1 - s) * (1 - t) * L[i][j] + s * (1 - t) * L[i + 1][j]
               + (1 - s) * t * L[i][j + 1] + s * t * L[i + 1][j + 1])
        return math.exp(val)
