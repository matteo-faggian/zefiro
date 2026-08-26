"""Raffreddamento a liquido: dimensionamento del canale.

Come per il termico, gli statuti sono diversi e i test lo rispettano:
  * bilanci e conduzione sono ESATTI -> tolleranze strette;
  * Dittus-Boelter e Blasius sono EMPIRICHE -> si verificano solo le
    proprieta' di scala che le correlazioni dichiarano;
  * Zuber e' derivata da instabilita' idrodinamica -> si verifica il
    comportamento qualitativo con la pressione e il confronto con CoolProp.
"""
from __future__ import annotations

import math

import pytest
from CoolProp.CoolProp import PropsSI

from zefiro.cooling import (
    channel_pressure_drop,
    coolant_mass_flow,
    dittus_boelter_h,
    required_h_for_no_boiling,
    reynolds,
    saturation_temperature,
    steady_wall_temperatures,
    zuber_chf,
)

T_IN = 288.15
P_W = 4.0e5
K_WALL = 16.0          # PLACEHOLDER 316L (TODO n.J)
Q_THROAT = 3.81e6      # W/m^2, dal calcolo con Bartz in gola


def test_saturazione_coincide_con_coolprop():
    """A 1 atm l'acqua bolle a 100 C: e' il controllo che il modulo stia
    davvero interrogando l'EOS e non una tabella a memoria."""
    assert saturation_temperature(101325.0) == pytest.approx(373.12, abs=0.1)
    assert saturation_temperature(P_W) == pytest.approx(
        PropsSI("T", "P", P_W, "Q", 0, "Water"), rel=1e-12
    )


def test_la_pressione_alza_la_saturazione_e_il_margine():
    """E' il motivo per cui avere buona pressione in rete non serve solo a
    spingere l'acqua: sposta la soglia di ebollizione."""
    assert saturation_temperature(2e5) < saturation_temperature(4e5) < saturation_temperature(6e5)


# --------------------------------------------------------------------------- #
# Bilanci e conduzione: esatti
# --------------------------------------------------------------------------- #
def test_le_due_cadute_di_temperatura_sono_esatte():
    """T_gas = T_acqua + q/h + q t/k. Sono due resistenze in serie, punto."""
    h = 45000.0
    t = 1.5e-3
    w = steady_wall_temperatures(Q_THROAT, t, K_WALL, h, T_IN, P_W)
    assert w.film_drop == pytest.approx(Q_THROAT / h, rel=1e-12)
    assert w.conduction_drop == pytest.approx(Q_THROAT * t / K_WALL, rel=1e-12)
    assert w.T_gas_side == pytest.approx(T_IN + w.film_drop + w.conduction_drop, rel=1e-12)


def test_la_parete_piu_sottile_e_piu_fredda_lato_gas():
    """Al contrario del pozzo termico, dove ispessire aiutava (fino a un punto):
    con raffreddamento sul dorso lo spessore e' solo resistenza in piu'."""
    h = 45000.0
    T = [steady_wall_temperatures(Q_THROAT, t, K_WALL, h, T_IN, P_W).T_gas_side
         for t in (0.8e-3, 1.5e-3, 2.4e-3)]
    assert T[0] < T[1] < T[2]


def test_portata_dal_bilancio_energetico():
    Q, dT = 15.0e3, 20.0
    m = coolant_mass_flow(Q, dT, T_IN, P_W)
    cp = PropsSI("C", "T", T_IN + 0.5 * dT, "P", P_W, "Water")
    assert m * cp * dT == pytest.approx(Q, rel=1e-12)
    assert m * 60.0 == pytest.approx(10.8, abs=0.5)      # L/min, acqua ~1 kg/L


def test_h_richiesto_e_coerente_con_le_temperature_di_parete():
    """required_h_for_no_boiling deve essere l'inverso esatto di
    steady_wall_temperatures: se uso quell'h, il margine e' il margine chiesto."""
    margine = 20.0
    h = required_h_for_no_boiling(Q_THROAT, T_IN, P_W, safety_margin_K=margine)
    w = steady_wall_temperatures(Q_THROAT, 1.5e-3, K_WALL, h, T_IN, P_W)
    assert w.margin_to_boiling == pytest.approx(margine, abs=1e-6)


def test_acqua_gia_calda_non_ha_margine():
    with pytest.raises(ValueError, match="margine"):
        required_h_for_no_boiling(Q_THROAT, saturation_temperature(P_W) - 5.0, P_W)


# --------------------------------------------------------------------------- #
# Correlazioni empiriche: si verificano le SCALE dichiarate
# --------------------------------------------------------------------------- #
def test_dittus_boelter_scala_come_velocita_alla_08():
    a = dittus_boelter_h(1.0e-3, 5.0, T_IN, P_W)
    b = dittus_boelter_h(1.0e-3, 10.0, T_IN, P_W)
    assert b / a == pytest.approx(2.0**0.8, rel=1e-9)


def test_dittus_boelter_premia_i_canali_stretti():
    """h ~ D^-0.2 a velocita' fissata: stringere il canale alza lo scambio.

    E' il motivo per cui il canale di gola va stretto e non largo, e per cui a
    parita' di portata conviene un gap piccolo."""
    a = dittus_boelter_h(2.0e-3, 10.0, T_IN, P_W)
    b = dittus_boelter_h(1.0e-3, 10.0, T_IN, P_W)
    assert b / a == pytest.approx(2.0**0.2, rel=1e-9)


def test_perdita_di_carico_scala_come_v_alla_175():
    """Blasius: f ~ Re^-0.25, quindi dp ~ V^2 V^-0.25 = V^1.75.

    Conseguenza pratica: raddoppiare la velocita' costa 3.4 volte la pressione,
    non 4. E' il motivo per cui alzare la velocita' resta conveniente finche'
    la pressione di rete lo permette."""
    a = channel_pressure_drop(1.0e-3, 0.04, 5.0, T_IN, P_W)
    b = channel_pressure_drop(1.0e-3, 0.04, 10.0, T_IN, P_W)
    assert b / a == pytest.approx(2.0**1.75, rel=1e-3)


def test_reynolds_turbolento_nel_campo_di_progetto():
    """Dittus-Boelter vuole Re > 10^4. A gap 0.5 mm servono almeno ~10 m/s."""
    assert reynolds(1.0e-3, 10.0, T_IN, P_W) > 8000.0
    assert reynolds(1.0e-3, 3.0, T_IN, P_W) < 4000.0     # qui la correlazione non vale


# --------------------------------------------------------------------------- #
# CHF: pavimento, non limite reale
# --------------------------------------------------------------------------- #
def test_chf_di_zuber_cresce_con_la_pressione():
    assert zuber_chf(1e5) < zuber_chf(4e5) < zuber_chf(8e5)


def test_il_flusso_di_gola_supera_la_chf_in_pool_boiling():
    """RISULTATO DI PROGETTO, e va tenuto ben presente: il flusso di gola e'
    circa il doppio della CHF in pool boiling saturo.

    Il progetto funziona SOLO perche' si tiene il liquido sottoraffreddato e in
    convezione forzata, cioe' si evita del tutto l'ebollizione. Non c'e'
    margine per un calo di portata: se l'acqua rallenta, si innesca
    l'ebollizione a film, il flusso crolla e la parete brucia in pochi secondi.
    Da qui la necessita' di un interblocco sulla portata."""
    assert Q_THROAT / zuber_chf(P_W) > 1.5


# --------------------------------------------------------------------------- #
# Il punto di progetto proposto
# --------------------------------------------------------------------------- #
def test_il_canale_di_progetto_non_bolle_e_sta_nella_rete_domestica():
    """gap 0.5 mm, 10 m/s, acqua a 4 bar e 15 C, parete di gola 1.5 mm."""
    gap, V = 0.5e-3, 10.0
    D_h = 2.0 * gap
    h = dittus_boelter_h(D_h, V, T_IN, P_W)
    w = steady_wall_temperatures(Q_THROAT, 1.5e-3, K_WALL, h, T_IN, P_W)
    dp = channel_pressure_drop(D_h, 0.040, V, T_IN, P_W)
    portata_Lmin = 2.0 * math.pi * 0.008 * gap * V * 1000.0 * 60.0

    assert w.margin_to_boiling > 30.0            # niente ebollizione, con margine
    assert w.T_gas_side < 273.15 + 550.0         # il 316L a 460 C ha tutta la sua resistenza
    assert dp < 1.0e5                            # meno di 1 bar nel canale
    assert 10.0 < portata_Lmin < 20.0            # compatibile con un rubinetto da giardino
