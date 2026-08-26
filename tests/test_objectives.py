"""Obiettivi e vincoli: convenzione dei segni e adimensionalizzazione.

La convenzione (minimizzare `f`, `g <= 0` fattibile) e' facile da rompere e
impossibile da accorgersene guardando i numeri: un segno invertito produce un
ottimizzatore che converge tranquillamente sulla soluzione peggiore.
"""
from __future__ import annotations

import dataclasses

import pytest

from zefiro.geometry.parameters import derive
from zefiro.opt.objectives import (
    CONSTRAINT_NAMES,
    MIN_INJECTOR_DP_FRACTION,
    OBJECTIVE_NAMES,
    objectives_l0,
    objectives_l1,
)
from zefiro.schemas import FEMResult


@pytest.fixture
def evaluated(operating_point, design_vector):
    params, l0 = derive(design_vector, operating_point)
    return params, l0, design_vector.values["p_c"]


def test_gli_obiettivi_da_massimizzare_entrano_col_segno_invertito(
    operating_point, evaluated
):
    params, l0, p_c = evaluated
    o = objectives_l0("r", l0, operating_point, p_c, derived=params.derived)
    assert o.f["neg_thrust"] == -l0.thrust
    assert o.f["neg_isp_total"] == -l0.Isp_s
    # spinta maggiore -> obiettivo minore: e' cio' che l'ottimizzatore cerca
    migliore = dataclasses.replace(l0, thrust=l0.thrust * 2.0)
    o2 = objectives_l0("r", migliore, operating_point, p_c, derived=params.derived)
    assert o2.f["neg_thrust"] < o.f["neg_thrust"]


def test_i_nomi_prodotti_stanno_nel_registro(operating_point, evaluated):
    """Il registro sono le colonne del database: un nome fuori registro
    verrebbe silenziosamente perso all'inserimento."""
    params, l0, p_c = evaluated
    o = objectives_l0("r", l0, operating_point, p_c, min_feature_size=0.4e-3,
                      derived=params.derived)
    assert set(o.f) <= set(OBJECTIVE_NAMES)
    assert set(o.g) <= set(CONSTRAINT_NAMES)


def test_vincolo_di_stabilita_delliniezione(operating_point, evaluated):
    """g = (0.15 p_c - Dp) / p_c: negativo quando il Dp e' sufficiente."""
    params, l0, p_c = evaluated
    dp = params.derived["dp_inj_fuel"]
    o = objectives_l0("r", l0, operating_point, p_c, derived=params.derived)
    atteso = (MIN_INJECTOR_DP_FRACTION * p_c - dp) / p_c
    assert o.g["fuel_dp_stability"] == pytest.approx(atteso, rel=1e-12)
    assert (o.g["fuel_dp_stability"] <= 0.0) == (dp >= MIN_INJECTOR_DP_FRACTION * p_c)


def test_vincolo_di_pressione_di_alimentazione(operating_point, evaluated):
    """p_c + Dp deve stare sotto la pressione di bombola."""
    params, l0, p_c = evaluated
    o = objectives_l0("r", l0, operating_point, p_c, derived=params.derived)
    fattibile = p_c + params.derived["dp_inj_fuel"] <= operating_point.p_fuel_supply
    assert (o.g["fuel_supply_pressure"] <= 0.0) == fattibile


def test_i_vincoli_sono_adimensionali(operating_point, evaluated):
    """Sommare un Dp in pascal con un margine strutturale puro darebbe alla
    pressione un peso arbitrario di 10^5 nella violazione aggregata."""
    params, l0, p_c = evaluated
    o = objectives_l0("r", l0, operating_point, p_c, min_feature_size=0.4e-3,
                      derived=params.derived)
    for name, v in o.g.items():
        assert abs(v) < 1.0e3, f"{name} = {v}: sembra dimensionale"


def test_un_vincolo_non_valutabile_resta_fuori(operating_point, evaluated):
    """Metterlo a zero direbbe 'soddisfatto', che e' falso. La differenza fra
    'verificato' e 'non verificabile' non va nascosta."""
    params, l0, p_c = evaluated
    senza_cd = dataclasses.replace(operating_point, cd_injector_fuel=None,
                                   cd_injector_ox=None)
    o = objectives_l0("r", l0, senza_cd, p_c, derived={})
    assert "fuel_dp_stability" not in o.g
    assert "min_feature" not in o.g


def test_fabbricabilita_come_vincolo(operating_point, evaluated):
    """Il foro piu' piccolo deve superare min_feature_size."""
    params, l0, p_c = evaluated
    d = params.derived
    quota_min = min(d[k] for k in ("d_ox", "d_fuel", "d_film", "film_land", "t_wall")
                    if d[k] > 0.0)
    o = objectives_l0("r", l0, operating_point, p_c, min_feature_size=0.4e-3,
                      derived=d)
    assert o.g["min_feature"] == pytest.approx((0.4e-3 - quota_min) / 0.4e-3, rel=1e-12)
    # a questa scala il vincolo e' ATTIVO: e' il risultato di progetto, non un bug
    assert o.g["min_feature"] > 0.0


def test_il_margine_del_fem_entra_col_segno_invertito(operating_point, evaluated):
    """FEMResult.margin_yield e' un margine (positivo = buono), `g` vuole il
    contrario. L'inversione avviene una volta sola, in objectives_l1."""
    params, l0, p_c = evaluated
    fem = FEMResult(run_id="r", cfd_id="c", T_wall_max=900.0, von_mises_max=1.0e8,
                    margin_yield=1.5, margin_temp=0.3, solver_version="test")
    o = objectives_l1("r", l0, operating_point, p_c, fem, derived=params.derived)
    assert o.g["margin_yield"] == -1.5
    assert o.g["margin_temp"] == -0.3
    assert o.fidelity == "L1"


def test_feasible_richiede_tutti_i_vincoli_non_positivi(operating_point, evaluated):
    params, l0, p_c = evaluated
    o = objectives_l0("r", l0, operating_point, p_c, derived=params.derived)
    assert o.feasible == all(v <= 0.0 for v in o.g.values())
