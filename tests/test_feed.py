"""Impianto di alimentazione: limiti del compressore, serbatoio, bombola GPL.

Le proprieta' del GPL vengono da CoolProp (equazioni di stato di riferimento).
I test verificano che il MODULO usi correttamente quei dati, non che i dati
siano giusti: quelli sono pubblicati e citati in feed.py.
"""
from __future__ import annotations

import math

import pytest
from CoolProp.CoolProp import PropsSI

from zefiro.feed import (
    CV_TO_W,
    HP_TO_W,
    R_AIR,
    autorefrigeration,
    burst_duration,
    composition_from_pressure,
    compressor_bounds,
    mdot_from_fad,
    saturation_pressure,
    tank_blowdown,
)


# --------------------------------------------------------------------------- #
# Compressore
# --------------------------------------------------------------------------- #
def test_lavoro_isotermo_coincide_col_conto_a_mano():
    """w_iso = R T ln(p2/p1). Per aria a 293.15 K e rapporto 10:

        287.05 * 293.15 * ln(10) = 193.8 kJ/kg
    """
    b = compressor_bounds(2237.0, 101325.0, 1.01325e6, 293.15)
    atteso = R_AIR * 293.15 * math.log(10.0)
    assert b.w_isothermal == pytest.approx(atteso, rel=1e-12)
    assert atteso / 1000.0 == pytest.approx(193.8, abs=0.5)


def test_isotermo_e_il_limite_inferiore_di_lavoro():
    """La compressione isoterma e' il lavoro MINIMO possibile: se il codice
    desse il contrario, il secondo principio sarebbe violato."""
    b = compressor_bounds(2237.0, 101325.0, 10e5, 293.15)
    assert b.w_isothermal < b.w_adiabatic
    assert b.mdot_isothermal > b.mdot_adiabatic


def test_portata_massima_di_un_3hp_e_dellordine_di_dieci_grammi_al_secondo():
    """Limite superiore ASSOLUTO, non prestazione di una macchina reale."""
    b = compressor_bounds(3 * HP_TO_W, 101325.0, 10e5, 293.15)
    assert 0.010 < b.mdot_isothermal < 0.013
    assert b.mdot_adiabatic < b.mdot_isothermal


def test_hp_e_cv_non_sono_la_stessa_cosa():
    """In Italia '3 cavalli' puo' voler dire 3 CV (735.5 W) o 3 HP (745.7 W).
    La differenza e' l'1.4 %: irrilevante qui, ma il codice non deve
    confonderli in silenzio."""
    assert HP_TO_W / CV_TO_W == pytest.approx(1.0139, abs=1e-3)


def test_fad_valutata_alle_condizioni_di_aspirazione():
    """La FAD e' un volume ASPIRATO: va convertita con la densita' ambiente,
    non con quella di mandata. E' l'errore piu' comune su questo numero."""
    rho_amb = 101325.0 / (R_AIR * 293.15)
    assert mdot_from_fad(300.0) == pytest.approx(300e-3 / 60.0 * rho_amb, rel=1e-12)
    assert mdot_from_fad(300.0) == pytest.approx(0.0060, abs=2e-4)


# --------------------------------------------------------------------------- #
# Serbatoio
# --------------------------------------------------------------------------- #
def test_massa_iniziale_nel_serbatoio_conto_a_mano():
    """m = p V / (R T) = 10e5 * 0.1 / (287.05 * 293.15) = 1.188 kg."""
    bd = tank_blowdown(0.100, 10e5, 6e5, 293.15)
    assert bd.mass_initial == pytest.approx(10e5 * 0.1 / (R_AIR * 293.15), rel=1e-12)
    assert bd.mass_initial == pytest.approx(1.188, abs=0.005)


def test_lo_svuotamento_adiabatico_rende_meno_massa():
    """Espandendosi il gas residuo si raffredda, quindi ne resta di piu' dentro
    a parita' di pressione finale. Usare il limite isotermo SOVRASTIMA la
    durata della raffica."""
    bd = tank_blowdown(0.100, 10e5, 6e5, 293.15)
    assert bd.mass_usable_adiabatic < bd.mass_usable_isothermal
    assert bd.T_final_adiabatic < 293.15


def test_raffica_infinita_se_il_compressore_tiene_il_passo():
    bd = tank_blowdown(0.100, 10e5, 6e5, 293.15)
    assert burst_duration(bd, 0.005, mdot_recharge=0.006) == (math.inf, math.inf)


def test_durata_raffica_scala_come_inverso_della_portata():
    bd = tank_blowdown(0.100, 10e5, 6e5, 293.15)
    t1, _ = burst_duration(bd, 0.100)
    t2, _ = burst_duration(bd, 0.200)
    assert t1 / t2 == pytest.approx(2.0, rel=1e-9)


def test_pressione_finale_maggiore_di_quella_iniziale_e_un_errore():
    with pytest.raises(ValueError):
        tank_blowdown(0.100, 6e5, 10e5, 293.15)


# --------------------------------------------------------------------------- #
# Bombola GPL
# --------------------------------------------------------------------------- #
def test_propano_puro_riproduce_la_tensione_di_vapore_di_riferimento():
    """A 20 C il propano ha p_sat = 8.37 bar (EOS di Lemmon 2009).

    E' il numero che rende leggibile il dato 'la bombola arriva a 8 bar'.
    """
    p = saturation_pressure({"C3H8": 1.0}, 293.15)
    assert p == pytest.approx(PropsSI("P", "T", 293.15, "Q", 0, "Propane"), rel=1e-12)
    assert p / 1e5 == pytest.approx(8.37, abs=0.05)


def test_il_butano_abbassa_molto_la_pressione():
    """n-butano a 20 C sta a ~2.08 bar: una miscela 70/30 non arriva a 8 bar.
    E' il motivo per cui la pressione misurata identifica la composizione."""
    p_mix = saturation_pressure({"C3H8": 0.7, "C4H10": 0.3}, 293.15)
    assert p_mix < saturation_pressure({"C3H8": 1.0}, 293.15)
    assert p_mix / 1e5 == pytest.approx(6.48, abs=0.15)


@pytest.mark.parametrize("x_prop", [0.5, 0.7, 0.9, 1.0])
@pytest.mark.parametrize("T", [283.15, 293.15, 303.15])
def test_inversione_di_raoult_e_consistente(x_prop, T):
    """composition_from_pressure deve invertire esattamente saturation_pressure."""
    p = saturation_pressure({"C3H8": x_prop, "C4H10": 1.0 - x_prop}, T)
    assert composition_from_pressure(p, T) == pytest.approx(x_prop, rel=1e-9)


def test_una_misura_incompatibile_lo_dice():
    """8 bar a 10 C darebbe x_propano > 1: impossibile per una binaria
    propano/butano. Il codice restituisce il valore fuori range invece di
    troncarlo a 1: e' un'informazione sulla misura, non un errore da nascondere."""
    assert composition_from_pressure(8.0e5, 283.15) > 1.0


def test_autorefrigerazione_chiude_il_bilancio_energetico():
    """Q * t = C * dT: e' il bilancio da cui nasce la formula."""
    a = autorefrigeration(0.0075, 293.15, {"C3H8": 1.0},
                          liquid_mass=6.0, shell_mass=8.5,
                          shell_specific_heat=470.0, burn_time=10.0)
    assert a.heat_rate * 10.0 == pytest.approx(a.thermal_capacity * a.dT_over_burn, rel=1e-9)
    assert a.p_end < a.p_start


def test_autorefrigerazione_trascurabile_su_una_raffica_breve():
    """Con 7.5 g/s di propano per 10 s la bombola perde circa 1 K.

    Conclusione di progetto: l'autorefrigerazione NON e' il vincolo che limita
    il funzionamento a raffica. Lo sono la capacita' del serbatoio d'aria e la
    portata ammessa dal riduttore.
    """
    a = autorefrigeration(0.0075, 293.15, {"C3H8": 1.0}, 6.0, 8.5, 470.0, 10.0)
    assert a.dT_over_burn < 2.0
    assert a.p_end / 1e5 > 8.0
