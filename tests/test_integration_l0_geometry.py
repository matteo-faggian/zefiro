"""Percorso completo: vettore di progetto -> L0 -> geometria -> STEP/STL."""
from __future__ import annotations

import math

import pytest

from zefiro.geometry import build_and_export
from zefiro.geometry.parameters import check_manufacturability, derive
from zefiro.store import make_run_id


def test_percorso_completo(operating_point, design_vector, tmp_path):
    params, l0 = derive(design_vector, operating_point)
    run_id = make_run_id(design_vector, operating_point, code_rev="test")
    art = build_and_export(params, tmp_path, run_id)

    assert art.is_valid_brep and art.is_watertight_mesh
    assert art.run_id == run_id
    assert art.sha256_step and art.sha256_stl

    # coerenza fra i due blocchi: l'area di gola usata dalla geometria e'
    # quella calcolata da L0, non una sua copia rimaneggiata
    assert params.derived["A_t"] == pytest.approx(l0.A_t, rel=1e-12)

    # il disaccordo fra eps a gamma costante ed eps in equilibrio spostato e'
    # un indicatore diagnostico, e a queste condizioni deve restare piccolo
    assert params.derived["epsilon_mismatch"] < 0.05


def test_conservazione_della_portata_fra_l0_e_geometria(operating_point, design_vector):
    """mdot = p_c * A_t / c*: identita' che lega i due blocchi."""
    params, l0 = derive(design_vector, operating_point)
    mdot_tot = l0.mdot_air + l0.mdot_fuel_core + l0.mdot_fuel_film
    assert design_vector.values["p_c"] * l0.A_t / l0.c_star == pytest.approx(mdot_tot, rel=1e-9)


def test_la_spinta_e_coerente_con_c_star_e_cf(operating_point, design_vector):
    """F = C_F * p_c * A_t, e per un ugello adattato C_F = u_e / c*."""
    _, l0 = derive(design_vector, operating_point)
    assert l0.thrust == pytest.approx(l0.C_F * design_vector.values["p_c"] * l0.A_t, rel=1e-9)


def test_gli_avvisi_sono_presenti_quando_devono(operating_point, design_vector):
    """A p_c/p_amb ~ 4 il codice deve dire che l'aerospike non e' nel suo regime."""
    _, l0 = derive(design_vector, operating_point)
    assert any("compensazione di quota" in w for w in l0.warnings)
    assert any("equilibrio chimico completo" in a for a in l0.assumptions)


def test_fabbricabilita_non_verificata_e_dichiarata(operating_point, design_vector):
    """Senza min_feature_size il risultato NON e' 'tutto ok', e' 'non verificato'."""
    params, _ = derive(design_vector, operating_point)
    msgs = check_manufacturability(params, None)
    assert len(msgs) == 1 and "NON sono stati verificati" in msgs[0]
