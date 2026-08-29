"""Ottimizzazione multi-obiettivo su L0: NSGA-II sui 12 parametri liberi.

PERCHE' QUESTI OBIETTIVI E NON ALTRI
------------------------------------
La scelta non e' di gusto: e' misurata da `scripts/objective_screening.py` su
un campione LHS di 400 punti del box di progetto, con correlazione di RANGO
(Spearman), che e' la statistica giusta perche' la dominanza di Pareto e' una
relazione d'ordine.

Il risultato (runs/objective_screening.json, seed 20260829):

    neg_thrust    <-> neg_isp_total     rho = +0.979
    neg_isp_total <-> q_throat          rho = -0.994
    neg_thrust    <-> q_throat          rho = -0.956
    neg_tau_res   <-> neg_damkohler     rho = +0.954

Entrambi gli estremi rovinano un fronte, in modi opposti e ugualmente inutili:

  * rho -> +1  (spinta e Isp): un obiettivo e' una copia dell'altro. Il fronte
    collassa su un punto e si e' sprecata una dimensione.
  * rho -> -1  (Isp e q di gola): OGNI punto e' non dominato. Il fronte e'
    tutto lo spazio di progetto, NSGA-II perde ogni pressione selettiva e
    restituisce rumore. Fisicamente la ragione e' che entrambe sono quasi
    funzioni monotone della sola p_c: quel "compromesso" e' una scansione a
    un parametro, non un problema multi-obiettivo, e va risolto come tale.

Restano tre metriche mutuamente quasi indipendenti (|rho| < 0.21 su tutte e
tre le coppie), che e' esattamente cio' che serve a un fronte tridimensionale
con informazione vera:

    neg_thrust    la spinta: e' cio' per cui il motore esiste
    wall_volume   la massa: tempo macchina e polvere in SLM, e il suo budget
    neg_isp_fuel  spinta per kg di GPL: l'aria e' gratis (compressore), il GPL
                  no (bombola finita). E' L'efficienza economica di Zefiro,
                  ed e' quasi ortogonale alla spinta (rho = 0.081) perche' la
                  spinta vuole lo stechiometrico e l'economia vuole il magro.

Il flusso termico di gola NON e' un obiettivo: e' un VINCOLO, al livello che
il raffreddamento ad acqua dichiarato riesce a estrarre. Cosi' la coppia
degenere (spinta, q) diventa "massimizza la spinta restando raffreddabile",
che e' il problema di progetto vero.
"""
from __future__ import annotations

import json
import math
import os
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np

from zefiro.geometry.parameters import DESIGN_BOUNDS, INTEGER_PARAMETERS, derive
from zefiro.opt.objectives import CONSTRAINT_NAMES, objectives_l0
from zefiro.schemas import DesignVector, OperatingPoint, ZefiroError

#: Insiemi di obiettivi ammessi. Un nome, non una lista sparsa nel codice: cosi'
#: una run dice in una parola su cosa ha ottimizzato, e il confronto fra run e'
#: possibile. Aggiungerne uno e' una scelta di progetto, non un dettaglio.
OBJECTIVE_SETS: dict[str, tuple[str, ...]] = {
    #: Il default, giustificato sopra dallo screening.
    "zefiro": ("neg_thrust", "wall_volume", "neg_isp_fuel"),
    #: Solo il compromesso prestazione/massa, per un fronte 2D leggibile.
    "thrust_mass": ("neg_thrust", "wall_volume"),
    #: Il compromesso termico esplicito, se si vuole vederlo come fronte
    #: invece che come vincolo. Attenzione: e' quasi degenere (rho = -0.956).
    "thrust_heat": ("neg_thrust", "q_throat"),
}


@dataclass
class OptimizationSetup:
    """Tutto cio' che serve a rendere una run ripetibile, in un oggetto solo."""

    op: OperatingPoint
    objective_set: str = "zefiro"
    tau_chem: Callable[[float, float], float] | None = None
    q_removable: float | None = None
    min_feature_size: float | None = None
    bounds: Mapping[str, tuple[float, float]] = field(default_factory=lambda: DESIGN_BOUNDS)

    @property
    def objective_names(self) -> tuple[str, ...]:
        if self.objective_set not in OBJECTIVE_SETS:
            raise ValueError(
                f"insieme di obiettivi {self.objective_set!r} sconosciuto; "
                f"quelli definiti sono {sorted(OBJECTIVE_SETS)}"
            )
        return OBJECTIVE_SETS[self.objective_set]

    @property
    def variable_names(self) -> tuple[str, ...]:
        return tuple(sorted(self.bounds))


#: Valore assegnato a un obiettivo quando la valutazione fallisce.
#:
#: NON e' "un numero grande a caso": e' il modo in cui NSGA-II rappresenta un
#: punto inammissibile senza fermarsi. Un `derive` che solleva (camera piu'
#: piccola dell'ugello, per esempio) e' un punto che non esiste, e deve essere
#: dominato da qualunque punto che esiste. Le run fallite vengono CONTATE e
#: riportate: se sono tante, il box e' mal posto e va corretto, non nascosto.
FAILED = 1.0e12


def _evaluate_one(
    x_vec: Sequence[float], setup: OptimizationSetup
) -> tuple[list[float], list[float], str | None]:
    names = setup.variable_names
    values = {n: float(v) for n, v in zip(names, x_vec)}
    for n in INTEGER_PARAMETERS:
        if n in values:
            values[n] = float(round(values[n]))
    try:
        x = DesignVector(values=values, bounds=dict(setup.bounds))
        params, l0 = derive(x, setup.op)
        tau = None
        if setup.tau_chem is not None:
            tau = setup.tau_chem(values["p_c"], values["phi_core"])
        obj = objectives_l0(
            run_id="", l0=l0, op=setup.op, p_c=values["p_c"],
            geometry=None, min_feature_size=setup.min_feature_size,
            derived=params.derived, tau_chem=tau, q_removable=setup.q_removable,
        )
    except (ZefiroError, ValueError, ArithmeticError, NotImplementedError) as e:
        n_f = len(setup.objective_names)
        return [FAILED] * n_f, [FAILED] * len(CONSTRAINT_NAMES), f"{type(e).__name__}: {e}"

    f = []
    for n in setup.objective_names:
        v = obj.f.get(n)
        if v is None or not math.isfinite(v):
            return ([FAILED] * len(setup.objective_names),
                    [FAILED] * len(CONSTRAINT_NAMES),
                    f"obiettivo {n!r} non calcolabile")
        f.append(float(v))
    # Un vincolo NON valutabile resta fuori da obj.g. Metterlo a 0 direbbe
    # "soddisfatto"; qui si mette a 0 solo perche' pymoo vuole un vettore
    # pieno, ed e' per questo che `unevaluated_constraints` lo dichiara.
    g = [float(obj.g.get(n, 0.0)) for n in CONSTRAINT_NAMES]
    return f, g, None


def _worker(args):
    return _evaluate_one(args[0], args[1])


class L0Problem:
    """Adattatore verso pymoo. Sta qui e non in `objectives` di proposito:
    la definizione fisica del problema non deve dipendere dall'ottimizzatore."""

    def __new__(cls, setup: OptimizationSetup, n_jobs: int = 1):
        from pymoo.core.problem import Problem

        names = setup.variable_names
        xl = np.array([setup.bounds[n][0] for n in names])
        xu = np.array([setup.bounds[n][1] for n in names])

        class _P(Problem):
            def __init__(self):
                super().__init__(n_var=len(names), n_obj=len(setup.objective_names),
                                 n_ieq_constr=len(CONSTRAINT_NAMES), xl=xl, xu=xu)
                self.setup = setup
                self.n_failed = 0
                self.failures: dict[str, int] = {}

            def _evaluate(self, X, out, *_a, **_k):
                if n_jobs > 1:
                    from concurrent.futures import ProcessPoolExecutor
                    with ProcessPoolExecutor(max_workers=n_jobs) as ex:
                        res = list(ex.map(_worker, [(row, setup) for row in X],
                                          chunksize=max(1, len(X) // (4 * n_jobs))))
                else:
                    res = [_evaluate_one(row, setup) for row in X]
                F, G = [], []
                for f, g, err in res:
                    F.append(f)
                    G.append(g)
                    if err is not None:
                        self.n_failed += 1
                        key = err.split(":")[0]
                        self.failures[key] = self.failures.get(key, 0) + 1
                out["F"] = np.array(F, dtype=float)
                out["G"] = np.array(G, dtype=float)

        return _P()


@dataclass
class ParetoRun:
    objective_names: tuple[str, ...]
    variable_names: tuple[str, ...]
    X: np.ndarray
    F: np.ndarray
    G: np.ndarray
    n_eval: int
    n_failed: int
    failures: dict[str, int]
    seed: int
    pop_size: int
    n_gen: int
    wall_time_s: float
    history_igd: list[float] = field(default_factory=list)

    def to_json(self) -> dict:
        return {
            "objective_names": list(self.objective_names),
            "variable_names": list(self.variable_names),
            "constraint_names": list(CONSTRAINT_NAMES),
            "X": self.X.tolist(), "F": self.F.tolist(), "G": self.G.tolist(),
            "n_eval": self.n_eval, "n_failed": self.n_failed,
            "failures": self.failures, "seed": self.seed,
            "pop_size": self.pop_size, "n_gen": self.n_gen,
            "wall_time_s": self.wall_time_s,
        }


def run_nsga2(
    setup: OptimizationSetup,
    pop_size: int,
    n_gen: int,
    seed: int,
    n_jobs: int = 1,
    verbose: bool = False,
    checkpoint: Path | None = None,
) -> ParetoRun:
    """NSGA-II sul modello L0.

    `seed` e' obbligatorio: senza, la run non e' riproducibile e non deve
    entrare nel database (docs/architettura.md sezione 5.1).

    `N_inj` e' intero. Si usa `RoundingRepair`, cioe' si arrotonda DENTRO
    l'algoritmo, non a valle: due individui che l'ottimizzatore considera
    diversi ma che arrotondano allo stesso intero avrebbero lo stesso run_id
    con due x diverse, e il versionamento delle run si romperebbe.
    """
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.operators.crossover.sbx import SBX
    from pymoo.operators.mutation.pm import PM
    from pymoo.operators.sampling.lhs import LHS
    from pymoo.optimize import minimize

    names = setup.variable_names
    interi = [i for i, n in enumerate(names) if n in INTEGER_PARAMETERS]

    problem = L0Problem(setup, n_jobs=n_jobs)
    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=LHS(),                     # copertura iniziale, non casuale pura
        crossover=SBX(prob=0.9, eta=15, vtype=float,
                      repair=_int_repair(interi) if interi else None),
        mutation=PM(eta=20, vtype=float,
                    repair=_int_repair(interi) if interi else None),
        eliminate_duplicates=True,
    )
    t0 = time.perf_counter()
    res = minimize(problem, algorithm, ("n_gen", n_gen), seed=seed,
                   verbose=verbose, save_history=False)
    dt = time.perf_counter() - t0

    X = np.atleast_2d(res.X)
    F = np.atleast_2d(res.F)
    G = np.atleast_2d(res.G) if res.G is not None else np.zeros((len(X), 0))
    for j in interi:
        X[:, j] = np.round(X[:, j])

    out = ParetoRun(
        objective_names=setup.objective_names, variable_names=names,
        X=X, F=F, G=G, n_eval=int(res.algorithm.evaluator.n_eval),
        n_failed=problem.n_failed, failures=dict(problem.failures),
        seed=seed, pop_size=pop_size, n_gen=n_gen, wall_time_s=dt,
    )
    if checkpoint is not None:
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
        checkpoint.write_text(json.dumps(out.to_json(), indent=1))
    return out


def _int_repair(columns: Sequence[int]):
    """Riparatore che arrotonda SOLO le colonne intere.

    `RoundingRepair` di pymoo fa `np.around(X).astype(int)` su TUTTO il
    vettore. Qui distruggerebbe undici variabili continue: p_c in pascal
    diventerebbe intero (innocuo), ma t_wall in metri vale ~2e-3 e diventerebbe
    ZERO. Un motore a spessore di parete nullo passerebbe silenziosamente in
    ottimizzazione, con volume minimo: sarebbe l'ottimo globale del problema
    sbagliato.
    """
    from pymoo.core.repair import Repair

    cols = list(columns)

    class _IntRepair(Repair):
        def _do(self, problem, X, **kwargs):
            X = np.asarray(X, dtype=float).copy()
            if X.ndim == 1:
                X[cols] = np.round(X[cols])
            else:
                X[:, cols] = np.round(X[:, cols])
            return X

    return _IntRepair()


def default_n_jobs() -> int:
    """Un processo per thread logico, meno uno per non affamare la macchina."""
    return max(1, (os.cpu_count() or 1) - 1)


# --- qualita' del fronte ----------------------------------------------------- #

def front_quality(runs: Sequence[ParetoRun]) -> dict:
    """Un fronte da un solo seme non e' un risultato: e' un campione.

    NSGA-II e' stocastico. Due semi diversi danno due fronti diversi, e la
    domanda che conta non e' "quanto e' bello questo fronte" ma "quanto
    dipende dal seme". Se due semi danno fronti molto diversi, il budget e'
    insufficiente e le conclusioni di progetto non reggono.

    Si misura con l'IPERVOLUME, normalizzato: gli obiettivi qui hanno unita' e
    scale incomparabili (spinta ~1e2 N, volume ~1e-5 m^3, Isp ~1e3 s) e un
    ipervolume su numeri grezzi sarebbe dominato dalla sola spinta. La
    normalizzazione usa ideale e nadir dell'UNIONE dei fronti, cosi' i semi si
    misurano sulla stessa scala.
    """
    if not runs:
        raise ValueError("nessuna run da confrontare")
    nomi = runs[0].objective_names
    for r in runs:
        if r.objective_names != nomi:
            raise ValueError("fronti con obiettivi diversi non sono confrontabili")

    tutti = np.vstack([r.F for r in runs])
    ideale = tutti.min(axis=0)
    nadir = tutti.max(axis=0)
    span = np.where(nadir > ideale, nadir - ideale, 1.0)

    from pymoo.indicators.hv import HV

    # riferimento appena oltre il nadir normalizzato: include tutto il fronte
    hv = HV(ref_point=np.full(len(nomi), 1.1))
    valori = [float(hv(np.atleast_2d((r.F - ideale) / span))) for r in runs]
    v = np.array(valori)
    return {
        "objective_names": list(nomi),
        "hypervolume": valori,
        "hv_mean": float(v.mean()),
        "hv_spread": float((v.max() - v.min()) / v.mean()) if v.mean() else float("nan"),
        "n_points": [int(len(r.F)) for r in runs],
        "seeds": [r.seed for r in runs],
        "ideal": ideale.tolist(),
        "nadir": nadir.tolist(),
    }
