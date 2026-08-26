"""Fixture condivise."""
from __future__ import annotations

import pytest

from zefiro.geometry.parameters import default_design_vector
from zefiro.schemas import FuelSpec, OperatingPoint

#: Punto operativo di PROVA, non quello reale di Zefiro.
#: Serve a rendere i test eseguibili mentre i TODO di config/ sono ancora
#: aperti. I valori di portata e temperatura sono dichiaratamente arbitrari e
#: NON vanno usati per concludere nulla sul motore.
TEST_MDOT_AIR = 0.010          # kg/s, valore di comodo
TEST_T_AIR = 300.0             # K
TEST_T_FUEL = 290.0            # K


@pytest.fixture
def propane_fuel() -> FuelSpec:
    return FuelSpec(
        composition={"C3H8": 1.0},
        phase_at_injection="gas",
        thermo_source="gri30.yaml",
    )


@pytest.fixture
def operating_point(propane_fuel: FuelSpec) -> OperatingPoint:
    return OperatingPoint(
        p_amb=101325.0,
        p_air_supply=10.0e5,
        p_fuel_supply=8.0e5,
        fuel=propane_fuel,
        T_air_in=TEST_T_AIR,
        T_fuel_in=TEST_T_FUEL,
        mdot_air_max=TEST_MDOT_AIR,
        cd_injector_ox=0.8,
        cd_injector_fuel=0.8,
    )


@pytest.fixture
def design_vector():
    return default_design_vector(
        p_c=4.0e5, phi_core=1.0, f_film=0.15, N_inj=12.0, plug_trunc=0.8
    )
