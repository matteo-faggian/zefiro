"""Stime strutturali analitiche: servono anche a verificare il FEM in fase 4."""
from __future__ import annotations

import pytest

from zefiro.structural import (
    biot_number,
    fully_constrained_thermal_stress,
    lame_hoop_stress_inner,
    thin_wall_hoop_stress,
)


def test_parete_sottile_e_lame_coincidono_quando_la_parete_e_sottile():
    p, r = 5.0e5, 0.0277
    for t, tol in ((0.001, 0.02), (0.0024, 0.05)):
        sottile = thin_wall_hoop_stress(p, r, t)
        esatta = lame_hoop_stress_inner(p, r, r + t)
        assert sottile == pytest.approx(esatta, rel=tol)


def test_parete_sottile_sbaglia_quando_la_parete_e_spessa():
    """Il criterio t/r < 0.1 non e' decorativo."""
    p, r, t = 5.0e5, 0.0277, 0.020        # t/r = 0.72
    sottile = thin_wall_hoop_stress(p, r, t)
    esatta = lame_hoop_stress_inner(p, r, r + t)
    assert abs(sottile - esatta) / esatta > 0.30


def test_lame_diverge_per_parete_evanescente():
    p, r = 5.0e5, 0.03
    assert lame_hoop_stress_inner(p, r, r * 1.0001) > 100.0 * p


def test_spessori_non_validi_sono_errori():
    with pytest.raises(ValueError):
        thin_wall_hoop_stress(1e5, 0.03, 0.0)
    with pytest.raises(ValueError):
        lame_hoop_stress_inner(1e5, 0.03, 0.03)


def test_tensione_termica_lineare_nel_gradiente():
    E, a, nu = 190e9, 17e-6, 0.30
    s1 = fully_constrained_thermal_stress(E, a, 100.0, nu)
    s2 = fully_constrained_thermal_stress(E, a, 200.0, nu)
    assert s2 == pytest.approx(2.0 * s1, rel=1e-12)


def test_a_bassa_pressione_il_termico_domina_di_ordini_di_grandezza():
    """RISULTATO DI PROGETTO: a 5 bar di camera la pressione produce ~6 MPa,
    mentre un gradiente di 100 K in parete impedita ne produce ~460.

    Conseguenza per la fase 4: il FEM di Zefiro e' un problema TERMO-elastico
    in cui il carico di pressione e' quasi irrilevante. Dimensionare sullo
    spessore a pressione sarebbe rispondere alla domanda sbagliata.
    """
    p_stress = thin_wall_hoop_stress(5.22e5, 0.0277, 0.0024)
    t_stress = fully_constrained_thermal_stress(190e9, 17e-6, 100.0, 0.30)
    assert p_stress / 1e6 == pytest.approx(6.0, abs=0.5)
    assert t_stress / p_stress > 50.0


def test_biot():
    assert biot_number(2400.0, 0.0024, 16.0) == pytest.approx(0.36, abs=0.01)
    assert biot_number(273.0, 0.0024, 16.0) < 0.1        # camera: quasi isoterma
