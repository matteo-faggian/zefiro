"""I contratti dati devono rifiutare gli input incoerenti, non aggiustarli."""
from __future__ import annotations

import pytest

from zefiro.geometry.parameters import DESIGN_BOUNDS, default_design_vector
from zefiro.schemas import (
    SCHEMA_VERSION,
    ContractViolation,
    DesignVector,
    FuelSpec,
    MissingDatum,
    Objectives,
    OperatingPoint,
)


def test_design_vector_fuori_dai_bounds_e_un_errore():
    """Volutamente un errore e non un clamp: un ottimizzatore che esce dai
    bounds ha un bug, e il clamp lo renderebbe invisibile."""
    vals = {k: 0.5 * (lo + hi) for k, (lo, hi) in DESIGN_BOUNDS.items()}
    vals["p_c"] = 99.0e5
    with pytest.raises(ContractViolation, match="p_c"):
        DesignVector(values=vals, bounds=DESIGN_BOUNDS)


def test_design_vector_con_chiavi_mancanti_e_un_errore():
    vals = {k: 0.5 * (lo + hi) for k, (lo, hi) in DESIGN_BOUNDS.items()}
    del vals["t_wall"]
    with pytest.raises(ContractViolation, match="t_wall"):
        DesignVector(values=vals, bounds=DESIGN_BOUNDS)


def test_composizione_combustibile_deve_sommare_a_uno():
    with pytest.raises(ContractViolation, match="1.0"):
        FuelSpec(composition={"C3H8": 0.6, "C4H10": 0.2},
                 phase_at_injection="gas", thermo_source="gri30.yaml")


def test_fase_di_iniezione_non_valida():
    with pytest.raises(ContractViolation):
        FuelSpec(composition={"C3H8": 1.0},
                 phase_at_injection="supercritical", thermo_source="gri30.yaml")


def test_dati_mancanti_elencati_tutti_insieme(propane_fuel):
    """Scoprire i TODO uno per volta a ogni esecuzione fallita e' il modo piu'
    lento possibile di raccogliere dati: `require` li elenca tutti."""
    op = OperatingPoint(p_amb=1e5, p_air_supply=10e5, p_fuel_supply=8e5, fuel=propane_fuel)
    with pytest.raises(MissingDatum) as exc:
        op.require("T_air_in", "T_fuel_in", "mdot_air_max")
    msg = str(exc.value)
    assert "T_air_in" in msg and "T_fuel_in" in msg and "mdot_air_max" in msg


def test_objectives_convenzione_di_fattibilita():
    o = Objectives(run_id="x", fidelity="L0", f={"neg_thrust": -10.0}, g={"a": -1.0, "b": 0.0})
    assert o.feasible
    assert not Objectives(run_id="x", fidelity="L0", f={}, g={"a": 1e-9}).feasible


def test_serializzazione_round_trip_esatta_sui_float():
    x = default_design_vector()
    d = x.to_dict()
    assert d["schema_version"] == SCHEMA_VERSION
    assert d["values"]["p_c"] == x.values["p_c"]
