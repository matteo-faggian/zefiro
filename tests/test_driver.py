"""Verifica del DRIVER, non di pymoo.

Distinzione che conta: pymoo e' testato a monte. Cio' che qui puo' rompersi e'
il modo in cui lo si usa - operatori, riparazione degli interi, segno dei
vincoli, ordine delle colonne. Per questo i test girano su problemi di cui si
conosce la soluzione ESATTA in forma chiusa: se il fronte calcolato non ci
finisce sopra, il difetto e' nel cablaggio.
"""
from __future__ import annotations

import numpy as np
import pytest

from zefiro.opt.driver import _int_repair


def _nsga2(problem, pop=80, gen=120, seed=1):
    from pymoo.algorithms.moo.nsga2 import NSGA2
    from pymoo.operators.sampling.lhs import LHS
    from pymoo.optimize import minimize

    return minimize(problem, NSGA2(pop_size=pop, sampling=LHS(),
                                   eliminate_duplicates=True),
                    ("n_gen", gen), seed=seed, verbose=False)


@pytest.mark.slow
def test_ritrova_il_fronte_analitico_di_zdt1():
    """ZDT1 ha fronte noto: f2 = 1 - sqrt(f1), con f1 in [0, 1].

    Tre criteri, tutti necessari, perche' ciascuno da solo si supera barando:
      * ACCURATEZZA: i punti stanno sulla curva.
      * COPERTURA: la coprono da un capo all'altro. Cento copie dello stesso
        punto passerebbero il primo criterio e sarebbero inutili.
      * UNIFORMITA': nessun buco largo. Due grumi agli estremi passerebbero i
        primi due.

    Le soglie non sono scelte a occhio: vengono dalla misura al variare del
    budget, con margine. Errore mediano dal fronte esatto, seed 1:

        pop  80 x 120 gen ( 9 600 valutazioni) -> 1.6e-2
        pop 100 x 250 gen (25 000 valutazioni) -> 8.4e-4
        pop 100 x 400 gen (40 000 valutazioni) -> 1.1e-4

    cioe' converge di due ordini di grandezza col budget, che e' il
    comportamento atteso di NSGA-II e la vera prova che il cablaggio e' giusto.
    Qui si usa il budget piu' alto e si chiede meno di quanto misurato.
    """
    from pymoo.problems import get_problem

    res = _nsga2(get_problem("zdt1"), pop=100, gen=400)
    f1, f2 = res.F[:, 0], res.F[:, 1]
    errore = np.abs(f2 - (1.0 - np.sqrt(f1)))
    assert np.median(errore) < 2.0e-3, f"errore mediano {np.median(errore):.2e}"
    assert errore.max() < 2.0e-2, f"errore massimo {errore.max():.4f}"
    assert f1.min() < 0.02 and f1.max() > 0.98, f"copertura [{f1.min():.3f}, {f1.max():.3f}]"
    buchi = np.diff(np.sort(f1))
    assert buchi.max() < 0.05, f"buco massimo {buchi.max():.3f}"


@pytest.mark.slow
def test_rispetta_un_vincolo_con_ottimo_noto():
    """Problema costruito perche' il vincolo MORDA e la risposta sia esatta.

        min (f1, f2) = (x1, x2)   su   x in [0, 1]^2
        vincolo: x1 + x2 >= 1     ->   g = 1 - x1 - x2 <= 0

    Senza vincolo l'ottimo sarebbe il solo punto (0, 0); col vincolo e'
    esattamente il segmento x1 + x2 = 1. La verifica DECISIVA sul segno e' che
    il vincolo risulti ATTIVO: g deve arrivare a zero da sotto. Col segno
    invertito la popolazione finirebbe nell'angolo opposto e g resterebbe
    lontano da zero, cosa che si vedrebbe subito.

    Sull'ultimo decimale NSGA-II lascia qualche individuo indietro (max |s-1|
    si assesta a ~1.2e-2 anche raddoppiando le generazioni, mentre la mediana
    sta a 5.5e-4): si verifica quindi la distribuzione, non il caso peggiore.
    """
    from pymoo.core.problem import Problem

    class P(Problem):
        def __init__(self):
            super().__init__(n_var=2, n_obj=2, n_ieq_constr=1,
                             xl=np.zeros(2), xu=np.ones(2))

        def _evaluate(self, X, out, *a, **k):
            out["F"] = X.copy()
            out["G"] = (1.0 - X.sum(axis=1)).reshape(-1, 1)

    res = _nsga2(P(), pop=100, gen=300)
    scarto = np.abs(res.F.sum(axis=1) - 1.0)
    assert (res.G <= 0.0).all(), "punti non ammissibili nel risultato"
    assert res.G.max() > -1.0e-6, f"vincolo non attivo: g_max = {res.G.max():.2e}"
    assert np.median(scarto) < 5.0e-3, f"mediana |s-1| = {np.median(scarto):.2e}"
    assert np.percentile(scarto, 95) < 2.0e-2, f"p95 |s-1| = {np.percentile(scarto,95):.2e}"


def test_repair_tocca_solo_le_colonne_intere():
    """Il difetto che questo test esclude e' reale: `RoundingRepair` di pymoo
    arrotonda TUTTO il vettore, e t_wall (in metri, ~2e-3) diventerebbe 0."""
    rep = _int_repair([1])
    X = np.array([[4.2e5, 11.7, 2.4e-3], [3.1e5, 6.2, 8.0e-4]])
    Y = rep._do(None, X)
    assert (Y[:, 1] == np.array([12.0, 6.0])).all()
    assert (Y[:, 0] == X[:, 0]).all() and (Y[:, 2] == X[:, 2]).all()


def test_repair_su_vettore_singolo():
    rep = _int_repair([0])
    assert rep._do(None, np.array([7.6, 1.5e-3]))[0] == 8.0


def test_front_quality_vede_la_differenza_fra_fronti():
    """`front_quality` deve dare uno spread PICCOLO su fronti simili e GRANDE
    su fronti diversi: se non distinguesse i due casi non servirebbe a nulla."""
    from zefiro.opt.driver import ParetoRun, front_quality

    def finto(F, seed):
        return ParetoRun(("a", "b"), ("x",), np.zeros((len(F), 1)), np.array(F),
                         np.zeros((len(F), 0)), len(F), 0, {}, seed, 10, 10, 0.0)

    t = np.linspace(0, 1, 25)
    buono = np.c_[t, 1 - np.sqrt(t)]                 # fronte convergito
    quasi = np.c_[t, 1 - np.sqrt(t) + 0.002]         # praticamente uguale
    scarso = np.c_[t, 1 - np.sqrt(t) + 0.30]         # nettamente peggiore

    simili = front_quality([finto(buono, 1), finto(quasi, 2)])
    diversi = front_quality([finto(buono, 1), finto(scarso, 2)])
    assert simili["hv_spread"] < 0.02, simili["hv_spread"]
    assert diversi["hv_spread"] > 0.10, diversi["hv_spread"]
