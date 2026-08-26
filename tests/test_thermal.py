"""Carico termico e risposta della parete.

Attenzione allo statuto epistemico diverso dei due blocchi:
  * `transient_slab` e' fisica esatta -> si verifica contro soluzioni
    ANALITICHE, con tolleranze strette;
  * `bartz_h_gas` e' una correlazione empirica -> si verificano solo le sue
    proprieta' di SCALA, che sono quelle che l'equazione dichiara. Verificarne
    il valore assoluto significherebbe verificare i dati del 1957, non il
    codice.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from zefiro.thermal import (
    STEFAN_BOLTZMANN,
    adiabatic_wall_temperature,
    bartz_h_gas,
    gas_radiation_flux,
    mean_beam_length,
    recovery_factor,
    semi_infinite_surface_temperature,
    transient_slab,
)

# Proprieta' di PROVA, non del 316L reale (TODO n.J). Servono solo a far girare
# i test: nessuna conclusione fisica su Zefiro dipende da questi numeri.
RHO, CP, K = 8000.0, 500.0, 16.0


# --------------------------------------------------------------------------- #
# Parete transitoria: fisica esatta
# --------------------------------------------------------------------------- #
def test_il_solutore_conserva_lenergia():
    """Il calore entrato dalla faccia calda deve trovarsi tutto immagazzinato
    nella parete: il dorso e' adiabatico, non c'e' altra via d'uscita."""
    r = transient_slab(0.004, RHO, CP, K, h_gas=2000.0, T_aw=2200.0, t_end=5.0)
    assert r.energy_in == pytest.approx(r.energy_stored, rel=1e-6)


def test_coincide_col_semispazio_ai_tempi_brevi():
    """Finche' la penetrazione e' molto minore dello spessore, la parete non
    "sa" di essere finita e deve comportarsi da semispazio.

    Confronto a flusso quasi costante: si usa h piccolo, cosi' la superficie si
    scalda poco e il flusso convettivo resta vicino al suo valore iniziale.
    """
    h, T_aw, t = 150.0, 2000.0, 1.0
    r = transient_slab(0.100, RHO, CP, K, h, T_aw, t_end=t, n_nodes=401, n_steps=4000)
    assert r.penetration_depth < 0.1 * 0.100      # davvero semispazio
    q0 = h * (T_aw - 293.15)
    analitico = semi_infinite_surface_temperature(q0, t, RHO, CP, K)
    # il numerico sta SOTTO l'analitico a flusso costante, perche' il flusso
    # vero cala man mano che la superficie si scalda
    assert r.T_hot[-1] < analitico
    assert r.T_hot[-1] == pytest.approx(analitico, rel=0.02)


def test_coincide_col_modello_concentrato_ai_tempi_lunghi():
    """Con Biot piccolo e tempo lungo la parete e' isoterma e tende
    all'equilibrio T_parete -> T_aw."""
    r = transient_slab(0.001, RHO, CP, K, h_gas=50.0, T_aw=1500.0,
                       t_end=600.0, n_steps=6000)
    assert r.biot < 0.01
    assert r.T_hot[-1] == pytest.approx(r.T_cold[-1], rel=1e-4)   # isoterma
    assert r.T_hot[-1] == pytest.approx(1500.0, rel=0.02)         # all'equilibrio


def test_il_gradiente_cresce_col_numero_di_biot():
    caldo = []
    for h in (200.0, 2000.0, 20000.0):
        r = transient_slab(0.005, RHO, CP, K, h, 2200.0, t_end=3.0)
        caldo.append(r.T_hot[-1] - r.T_cold[-1])
    assert caldo[0] < caldo[1] < caldo[2]


def test_lo_spessore_satura(): 
    """RISULTATO DI PROGETTO, non solo di codice: oltre una certa profondita'
    aumentare lo spessore non abbassa piu' la temperatura della faccia calda.

    Motivo fisico: superata la profondita' di penetrazione, il materiale in
    piu' non viene raggiunto dal calore nel tempo disponibile, quindi non
    partecipa. Il limite diventa la CONDUCIBILITA', non la capacita' termica.
    """
    T = [transient_slab(s, RHO, CP, K, 2400.0, 2230.0, t_end=5.0).T_hot[-1]
         for s in (0.0024, 0.008, 0.012, 0.020, 0.040)]
    assert T[0] > T[1] > T[2]                       # all'inizio serve
    assert T[-1] == pytest.approx(T[-2], rel=0.01)  # poi non serve piu'


def test_spessore_o_proprieta_non_positive_sono_errori():
    for kwargs in ({"thickness": 0.0}, {"rho": -1.0}, {"k": 0.0}):
        base = dict(thickness=0.004, rho=RHO, cp=CP, k=K)
        base.update(kwargs)
        with pytest.raises(ValueError):
            transient_slab(h_gas=100.0, T_aw=1000.0, t_end=1.0, **base)


# --------------------------------------------------------------------------- #
# Temperatura di parete adiabatica: fisica esatta
# --------------------------------------------------------------------------- #
def test_a_mach_zero_la_parete_adiabatica_e_al_ristagno():
    assert adiabatic_wall_temperature(2200.0, 0.0, 1.25, 0.7) == pytest.approx(2200.0)


def test_a_mach_uno_la_parete_adiabatica_sta_fra_statica_e_ristagno():
    T0, g, Pr = 2200.0, 1.25, 0.7
    T_aw = adiabatic_wall_temperature(T0, 1.0, g, Pr)
    T_static = T0 / (1.0 + 0.5 * (g - 1.0))
    assert T_static < T_aw < T0
    # con r = Pr^(1/3) < 1 il recupero e' parziale, ed e' proprio questo che
    # rende T_aw minore della temperatura di ristagno
    assert recovery_factor(Pr) == pytest.approx(Pr ** (1 / 3))
    assert recovery_factor(Pr, turbulent=False) == pytest.approx(math.sqrt(Pr))


# --------------------------------------------------------------------------- #
# Bartz: si verificano solo le proprieta' di SCALA dichiarate dall'equazione
# --------------------------------------------------------------------------- #
def _h(**kw):
    args = dict(D_t=0.016, p_c=5.0e5, c_star=1250.0, mu=7.0e-5, cp=1450.0,
                prandtl=0.7, gamma=1.25, area_ratio_local=1.0, mach_local=1.0,
                T_wall=800.0, T_c=2250.0)
    args.update(kw)
    return bartz_h_gas(**args)


def test_scala_come_pressione_alla_08():
    assert _h(p_c=1.0e6) / _h(p_c=5.0e5) == pytest.approx(2.0**0.8, rel=1e-9)


def test_scala_come_diametro_di_gola_alla_meno_02():
    assert _h(D_t=0.032) / _h(D_t=0.016) == pytest.approx(2.0**-0.2, rel=1e-9)


def test_scala_come_rapporto_daree_alla_meno_09():
    """(A_t/A)^0.9: e' il termine per cui la camera vede molto meno flusso
    della gola. Con contrazione 11.7 il rapporto e' circa 9."""
    assert _h(area_ratio_local=10.0) / _h(area_ratio_local=1.0) == \
        pytest.approx(10.0**-0.9, rel=1e-9)


def test_la_gola_e_il_punto_piu_caldo():
    assert _h(area_ratio_local=1.0) > _h(area_ratio_local=5.0) > _h(area_ratio_local=12.0)


def test_una_parete_piu_fredda_riceve_piu_flusso():
    """Il fattore sigma di Bartz: la parete fredda addensa il gas nello strato
    limite e alza lo scambio."""
    assert _h(T_wall=400.0) > _h(T_wall=1200.0)


# --------------------------------------------------------------------------- #
# Irraggiamento: fisica esatta
# --------------------------------------------------------------------------- #
def test_stefan_boltzmann_esatta():
    q = gas_radiation_flux(1.0, 1000.0, 0.0)
    assert q == pytest.approx(STEFAN_BOLTZMANN * 1000.0**4, rel=1e-12)


def test_flusso_radiativo_nullo_allequilibrio():
    assert gas_radiation_flux(1.0, 1200.0, 1200.0) == pytest.approx(0.0, abs=1e-9)


def test_il_radiativo_conta_in_camera_e_non_in_gola():
    """A 2258 K il radiativo e' una frazione importante del convettivo di
    camera (basso) e trascurabile rispetto a quello di gola (alto)."""
    q_rad = gas_radiation_flux(0.15, 2258.0, 700.0)
    q_conv_camera = 273.0 * (2258.0 - 700.0)
    q_conv_gola = 2412.0 * (2230.0 - 700.0)
    assert 0.10 < q_rad / q_conv_camera < 1.0
    assert q_rad / q_conv_gola < 0.10


def test_lunghezza_media_del_raggio():
    """Sfera di raggio R: V/A = R/3, quindi L = 3.6 R/3 = 1.2 R."""
    R = 0.05
    V = 4.0 / 3.0 * math.pi * R**3
    A = 4.0 * math.pi * R**2
    assert mean_beam_length(V, A) == pytest.approx(1.2 * R, rel=1e-12)
