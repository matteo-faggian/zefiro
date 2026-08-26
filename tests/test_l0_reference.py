"""Verifica di L0 su un caso di riferimento CALCOLABILE A MANO.

Il caso: propano puro + aria secca, stechiometrico, 300 K, 1 atm.
Scelto perche' entrambe le grandezze verificate hanno un riscontro indipendente:
il rapporto aria/combustibile si ricava con carta e penna dalla stechiometria,
e la temperatura adiabatica di fiamma di questa miscela e' un valore
ampiamente riportato in letteratura (circa 2267 K).
"""
from __future__ import annotations

import math

import cantera as ct
import pytest

from zefiro.l0.equilibrium import chamber_state
from zefiro.l0.mixture import MixtureModel
from zefiro.l0.nozzle import FROZEN, SHIFTING, _find_throat, expand_to_pressure
from zefiro.schemas import FuelSpec, MissingThermoData
from zefiro.units import AIR_MOLE_FRACTIONS


def test_afr_stechiometrico_coincide_col_conto_a_mano(propane_fuel):
    """CONTO A MANO.

        C3H8 + 5 O2 -> 3 CO2 + 4 H2O

    L'aria secca usata dal progetto (units.AIR_MOLE_FRACTIONS, frazioni molari
    rinormalizzate) ha x_O2 = 0.20946/(somma) e massa molare

        MW_aria = SOMMA_i x_i MW_i = 28.9657 kg/kmol

    Servono 5 kmol di O2 per kmol di C3H8, quindi

        n_aria/n_fuel = 5 / x_O2         = 23.8709
        AFR           = 23.8709 * 28.9657 / 44.0956 = 15.6799

    Il test rifa' esattamente questo conto con i pesi molecolari presi dal
    meccanismo, e lo confronta col valore prodotto dal modulo.
    """
    model = MixtureModel.from_fuel(propane_fuel)
    afr_code = model.afr_stoichiometric()

    gas = ct.Solution(propane_fuel.thermo_source)
    tot = sum(AIR_MOLE_FRACTIONS.values())
    x = {k: v / tot for k, v in AIR_MOLE_FRACTIONS.items()}
    mw = {k: gas.molecular_weights[gas.species_index(k)] for k in x}
    mw_air = sum(x[k] * mw[k] for k in x)
    mw_fuel = gas.molecular_weights[gas.species_index("C3H8")]
    afr_hand = (5.0 / x["O2"]) * mw_air / mw_fuel

    assert afr_hand == pytest.approx(15.6799, abs=1e-3)      # il conto a mano
    assert afr_code == pytest.approx(afr_hand, rel=1e-9)     # il codice ci coincide


def test_temperatura_adiabatica_di_fiamma_valore_di_letteratura(propane_fuel):
    """T_ad di propano/aria stechiometrico a 300 K e 1 atm ~ 2267 K.

    Tolleranza +-25 K: il valore di letteratura dipende leggermente dalla
    composizione dell'aria assunta e dall'insieme di specie del meccanismo, e
    una tolleranza piu' stretta verificherebbe il meccanismo, non il codice.
    """
    model = MixtureModel.from_fuel(propane_fuel)
    afr = model.afr_stoichiometric()
    state = chamber_state(model, ct.one_atm, 300.0, 300.0, mdot_air=afr, mdot_fuel=1.0)

    assert state.T == pytest.approx(2267.0, abs=25.0)
    assert 1.20 < state.gamma < 1.30
    # prodotti dominanti attesi da C3H8 + aria
    assert state.X["N2"] > 0.65
    assert state.X["H2O"] == pytest.approx(0.148, abs=0.02)
    assert state.X["CO2"] == pytest.approx(0.103, abs=0.02)


def test_la_miscelazione_conserva_lentalpia(propane_fuel):
    """Bilancio del miscelatore adiabatico: H_out = H_air + H_fuel.

    Verifica che aria a 300 K e GPL a 290 K vengano miscelati per ENTALPIA e
    non per temperatura media, che sarebbe sbagliato con cp variabile.
    """
    model = MixtureModel.from_fuel(propane_fuel)
    gas = model.gas
    p, m_a, m_f = 4.0e5, 1.0, 0.05

    gas.TPX = 300.0, p, model.air_string
    h_air = gas.enthalpy_mass
    gas.TPX = 290.0, p, model.fuel_string
    h_fuel = gas.enthalpy_mass

    state = chamber_state(model, p, 300.0, 290.0, m_a, m_f)
    h_expected = (m_a * h_air + m_f * h_fuel) / (m_a + m_f)
    assert state.h == pytest.approx(h_expected, rel=1e-9)


def test_gola_sonica_e_condizione_di_massimo_flusso(propane_fuel):
    """La gola trovata come massimo di rho*u deve risultare sonica: M = 1.

    E' il controllo che il metodo usato nell'ugello (nessun gamma imposto,
    nessuna formula chiusa, la gola come argmax del flusso specifico) produca
    davvero la condizione critica. La verifica si fa in modo CONGELATO perche'
    e' li' che la velocita' del suono di Cantera (congelata) e' quella giusta:
    in equilibrio spostato il Mach critico va riferito alla velocita' del suono
    di equilibrio, che e' una grandezza diversa.
    """
    model = MixtureModel.from_fuel(propane_fuel)
    afr = model.afr_stoichiometric()
    state = chamber_state(model, 4.0e5, 300.0, 290.0, afr, 1.0)

    gas = model.gas
    gas.TPX = state.T, state.p, state.X
    s0, h0 = gas.entropy_mass, gas.enthalpy_mass
    p_t = _find_throat(gas, s0, h0, state.p, FROZEN)

    gas.SP = s0, p_t
    u = math.sqrt(max(2.0 * (h0 - gas.enthalpy_mass), 0.0))
    assert u / gas.sound_speed == pytest.approx(1.0, abs=2e-3)

    # e il rapporto critico deve essere quello isentropico con QUEL gamma
    g = gas.cp_mass / gas.cv_mass
    assert p_t / state.p == pytest.approx((2.0 / (g + 1.0)) ** (g / (g - 1.0)), rel=5e-3)


def test_equilibrio_spostato_da_prestazione_maggiore_del_congelato(propane_fuel):
    """L'espansione in equilibrio spostato deve dare u_e >= congelato.

    Motivo fisico: ricombinando, la miscela restituisce al flusso parte
    dell'energia di dissociazione. Se il codice desse il contrario, ci sarebbe
    un errore di segno da qualche parte.
    """
    model = MixtureModel.from_fuel(propane_fuel)
    afr = model.afr_stoichiometric()
    state = chamber_state(model, 4.0e5, 300.0, 290.0, afr, 1.0)
    frozen = expand_to_pressure(model.gas, state, 101325.0, 101325.0, mode=FROZEN)
    shifting = expand_to_pressure(model.gas, state, 101325.0, 101325.0, mode=SHIFTING)
    assert shifting.u_e >= frozen.u_e
    assert shifting.c_star >= frozen.c_star * 0.98


def test_butano_fallisce_invece_di_essere_sostituito():
    """gri30 non contiene C4H10: il codice deve DIRLO, non ripiegare sul propano."""
    fuel = FuelSpec(
        composition={"C3H8": 0.7, "C4H10": 0.3},
        phase_at_injection="gas",
        thermo_source="gri30.yaml",
    )
    with pytest.raises(MissingThermoData, match="C4H10"):
        MixtureModel.from_fuel(fuel)
