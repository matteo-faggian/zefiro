"""Il DOE deve essere deterministico e davvero un Latin Hypercube."""
from __future__ import annotations

import pytest

from zefiro.doe import latin_hypercube, projection_is_stratified
from zefiro.geometry.parameters import DESIGN_BOUNDS, INTEGER_PARAMETERS

N = 32


@pytest.fixture(scope="module")
def sample():
    return latin_hypercube(N, DESIGN_BOUNDS, seed=0, integer_names=INTEGER_PARAMETERS)


def test_stesso_seme_stesso_campione(sample):
    assert latin_hypercube(N, DESIGN_BOUNDS, 0, INTEGER_PARAMETERS) == sample


def test_semi_diversi_campioni_diversi(sample):
    assert latin_hypercube(N, DESIGN_BOUNDS, 1, INTEGER_PARAMETERS) != sample


def test_tutti_i_punti_sono_dentro_i_bounds(sample):
    for point in sample:
        for name, v in point.items():
            lo, hi = DESIGN_BOUNDS[name]
            assert lo <= v <= hi, f"{name} = {v} fuori da [{lo}, {hi}]"


@pytest.mark.parametrize("name", [n for n in DESIGN_BOUNDS if n not in INTEGER_PARAMETERS])
def test_ogni_dimensione_continua_e_stratificata(sample, name):
    """La proprieta' che DEFINISCE un Latin Hypercube: un punto per strato su
    ogni proiezione monodimensionale. Se salta, non e' un LHS."""
    assert projection_is_stratified(sample, DESIGN_BOUNDS, name)


def test_le_dimensioni_intere_sono_intere(sample):
    for point in sample:
        assert point["N_inj"] == int(point["N_inj"])
        lo, hi = DESIGN_BOUNDS["N_inj"]
        assert lo <= point["N_inj"] <= hi


def test_un_solo_punto_e_ammesso():
    assert len(latin_hypercube(1, DESIGN_BOUNDS, 0)) == 1


def test_zero_punti_e_un_errore():
    with pytest.raises(ValueError):
        latin_hypercube(0, DESIGN_BOUNDS, 0)


def test_lhs_copre_meglio_del_casuale_puro():
    """Verifica che valga la pena usare LHS invece di campionare a caso.

    Su una dimensione, la discrepanza rispetto alla distribuzione uniforme
    (statistica di Kolmogorov-Smirnov a un campione) deve essere piu' piccola
    per LHS che per un campione casuale, in modo sistematico.
    """
    import numpy as np

    b = {"u": (0.0, 1.0)}
    def ks(vals):
        v = np.sort(np.array(vals))
        n = len(v)
        return float(np.max(np.abs(v - (np.arange(1, n + 1) - 0.5) / n)))

    lhs_worse = 0
    for seed in range(20):
        lhs = ks([p["u"] for p in latin_hypercube(50, b, seed)])
        rnd = ks(np.random.Generator(np.random.PCG64(seed)).random(50))
        if lhs >= rnd:
            lhs_worse += 1
    assert lhs_worse <= 2, f"LHS peggiore del casuale in {lhs_worse}/20 semi"
