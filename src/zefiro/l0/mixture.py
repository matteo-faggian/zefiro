"""Composizione della miscela: AFR stechiometrico, portate, rapporto di equivalenza.

Nessun valore tabellato: AFR viene dal bilancio degli elementi calcolato da
Cantera sui polinomi NASA del meccanismo dichiarato.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping

import cantera as ct

from zefiro.schemas import FuelSpec, MissingThermoData
from zefiro.units import air_composition_string

# Sopra questo phi, l'insieme di specie di gri30 (nessuna C2+ pesante, nessun
# precursore di fuliggine) non descrive piu' l'equilibrio ricco.
RICH_EQUILIBRIUM_LIMIT = 1.5


@dataclass
class MixtureModel:
    """Modello di miscela legato a un preciso file di meccanismo Cantera."""

    fuel: FuelSpec
    gas: ct.Solution

    @classmethod
    def from_fuel(cls, fuel: FuelSpec) -> "MixtureModel":
        gas = ct.Solution(fuel.thermo_source)
        missing = [s for s in fuel.composition if s not in gas.species_names]
        if missing:
            raise MissingThermoData(
                f"Il meccanismo {fuel.thermo_source!r} non contiene {missing}. "
                "Non esiste alcuna sostituzione automatica: procurati un meccanismo "
                "che le contenga (scripts/fetch_mechanism.py scarica il San Diego mech, "
                "che include n-C4H10) e aggiorna FuelSpec.thermo_source."
            )
        return cls(fuel=fuel, gas=gas)

    # -- composizioni come stringhe Cantera ------------------------------- #
    @property
    def fuel_string(self) -> str:
        return ", ".join(f"{k}:{v!r}" for k, v in self.fuel.composition.items())

    @property
    def air_string(self) -> str:
        return air_composition_string(list(self.gas.species_names))

    # -- grandezze derivate ------------------------------------------------ #
    def afr_stoichiometric(self) -> float:
        """Rapporto aria/combustibile stechiometrico in MASSA.

        Verifica a mano per propano puro e aria (O2 .20946 / N2 .78084 /
        AR .00934 / CO2 .00036):
            C3H8 + 5 O2 -> 3 CO2 + 4 H2O
            moli aria / mole fuel = 5 / 0.20946 = 23.870
            MW_aria = 28.949 kg/kmol ; MW_C3H8 = 44.096 kg/kmol
            AFR = 23.870 * 28.949 / 44.096 = 15.67
        """
        g = self.gas
        g.TP = 300.0, ct.one_atm
        g.set_equivalence_ratio(1.0, self.fuel_string, self.air_string, basis="mole")
        y_fuel = sum(g.mass_fraction_dict().get(s, 0.0) for s in self.fuel.composition)
        if y_fuel <= 0.0:
            raise MissingThermoData("Frazione di combustibile nulla: composizione incoerente.")
        return (1.0 - y_fuel) / y_fuel

    def mass_flows(
        self, mdot_air: float, phi_core: float, f_film: float
    ) -> tuple[float, float, float]:
        """(mdot_fuel_core, mdot_fuel_film, phi_global).

        Definizione adottata (docs/architettura.md sezione 4.1): la variabile
        libera e' `phi_core`, cioe' il rapporto di equivalenza visto dal core,
        che riceve TUTTA l'aria e la frazione (1 - f_film) del combustibile.
        Ne segue phi_global = phi_core / (1 - f_film).
        """
        if not (0.0 <= f_film < 1.0):
            raise ValueError(f"f_film deve stare in [0, 1), vale {f_film!r}")
        afr = self.afr_stoichiometric()
        mdot_fuel_core = mdot_air / afr * phi_core
        mdot_fuel_total = mdot_fuel_core / (1.0 - f_film)
        return mdot_fuel_core, mdot_fuel_total - mdot_fuel_core, phi_core / (1.0 - f_film)

    def phi_from_flows(self, mdot_air: float, mdot_fuel: float) -> float:
        return (mdot_air / mdot_fuel) ** -1 * self.afr_stoichiometric()

    def equilibrium_validity_warning(self, phi: float) -> str | None:
        if phi > RICH_EQUILIBRIUM_LIMIT:
            return (
                f"phi = {phi:.2f} > {RICH_EQUILIBRIUM_LIMIT}: l'insieme di specie del "
                f"meccanismo {self.fuel.thermo_source!r} probabilmente non contiene gli "
                "idrocarburi pesanti e i precursori di fuliggine che dominano l'equilibrio "
                "ricco. Il T_ad calcolato e' sovrastimato."
            )
        return None
