"""Verifiche del contorno dell'aerospike (metodo di Angelino).

Ogni test controlla una proprieta' che DEVE valere per costruzione: se una
salta, l'errore e' nella derivazione, non nella tolleranza.
"""
from __future__ import annotations

import math

import pytest

from zefiro.geometry.aerospike import (
    area_ratio,
    mach_from_area_ratio,
    mach_from_pressure_ratio,
    plug_contour,
    prandtl_meyer,
)

GAMMAS = [1.15, 1.25, 1.30, 1.40]


@pytest.mark.parametrize("gamma", GAMMAS)
def test_prandtl_meyer_si_annulla_a_mach_uno(gamma):
    assert prandtl_meyer(1.0, gamma) == 0.0
    assert prandtl_meyer(1.0000001, gamma) < 1e-4


@pytest.mark.parametrize("gamma", GAMMAS)
def test_prandtl_meyer_tende_al_limite_teorico(gamma):
    """nu(M -> inf) = (pi/2)(sqrt((g+1)/(g-1)) - 1). E' un limite esatto."""
    limite = 0.5 * math.pi * (math.sqrt((gamma + 1) / (gamma - 1)) - 1.0)
    assert prandtl_meyer(1.0e4, gamma) == pytest.approx(limite, rel=1e-3)


@pytest.mark.parametrize("gamma", GAMMAS)
def test_area_ratio_minimo_a_mach_uno(gamma):
    assert area_ratio(1.0, gamma) == pytest.approx(1.0, rel=1e-12)
    assert area_ratio(0.5, gamma) > 1.0
    assert area_ratio(2.0, gamma) > 1.0


@pytest.mark.parametrize("gamma", GAMMAS)
@pytest.mark.parametrize("mach", [1.2, 1.6, 2.5, 4.0])
def test_inversione_area_ratio(gamma, mach):
    assert mach_from_area_ratio(area_ratio(mach, gamma), gamma) == pytest.approx(mach, rel=1e-9)


def test_mach_da_rapporto_di_pressione():
    """A p0/p = 1.8929 con gamma = 1.4 corrisponde M = 1 (rapporto critico)."""
    g = 1.4
    pr_crit = ((g + 1.0) / 2.0) ** (g / (g - 1.0))
    assert mach_from_pressure_ratio(pr_crit, g) == pytest.approx(1.0, rel=1e-9)


@pytest.mark.parametrize("gamma", [1.20, 1.25, 1.33])
@pytest.mark.parametrize("M_e", [1.4, 2.0, 3.0])
def test_il_plug_chiude_esattamente_sullasse(gamma, M_e):
    """Al Mach di progetto il raggio del plug deve annullarsi.

    Non e' una coincidenza numerica: sostituendo M = M_e nella formula si ha
    A/A* = eps e sin(alpha) = sin(mu_e) = 1/M_e, quindi l'argomento della
    radice diventa 1 - M_e * eps * (1/M_e) / eps = 0. E' la verifica che
    l'implementazione riproduce l'algebra.
    """
    c = plug_contour(gamma, M_e, R_lip=5.0e-3)
    assert c.r[-1] == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("gamma", [1.20, 1.25, 1.33])
@pytest.mark.parametrize("M_e", [1.4, 2.0, 3.0])
def test_contorno_monotono(gamma, M_e):
    """x strettamente crescente, r non crescente: e' un contorno di ugello."""
    c = plug_contour(gamma, M_e, R_lip=5.0e-3)
    assert all(c.x[i] < c.x[i + 1] for i in range(len(c.x) - 1))
    assert all(c.r[i] >= c.r[i + 1] - 1e-15 for i in range(len(c.r) - 1))


@pytest.mark.parametrize("gamma", [1.20, 1.25, 1.33])
@pytest.mark.parametrize("M_e", [1.4, 2.0, 3.0])
def test_area_di_gola_coerente_con_epsilon(gamma, M_e):
    """L'area di gola ricavata dalla GEOMETRIA deve valere pi*R_lip^2/eps.

    La gola e' inclinata di nu_e rispetto al piano radiale, quindi l'area
    normale al flusso e' l'anello diviso cos(nu_e), non l'anello. Questo test
    e' il controllo incrociato fra formula del contorno e definizione di eps.
    """
    R = 5.0e-3
    c = plug_contour(gamma, M_e, R_lip=R)
    A_geom = math.pi * (R**2 - c.r_throat**2) / math.cos(c.nu_e)
    A_def = math.pi * R**2 / c.epsilon
    assert A_geom == pytest.approx(A_def, rel=1e-9)


def test_il_punto_di_gola_e_a_monte_del_labbro():
    """x_throat < 0: la gola inclinata cade a MONTE del labbro. E' corretto."""
    c = plug_contour(1.25, 2.0, R_lip=5.0e-3)
    assert c.x_throat < 0.0


def test_troncamento_accorcia_senza_cambiare_la_gola():
    """Troncare il plug non deve toccare la gola: e' un taglio a valle."""
    full = plug_contour(1.25, 2.0, R_lip=5.0e-3, truncation=1.0)
    trunc = plug_contour(1.25, 2.0, R_lip=5.0e-3, truncation=0.5)
    assert trunc.r_throat == pytest.approx(full.r_throat, rel=1e-12)
    assert trunc.x[-1] < full.x[-1]
    assert trunc.r[-1] > 0.0            # base tronca, non punta


def test_mach_di_uscita_subsonico_rifiutato():
    with pytest.raises(ValueError):
        plug_contour(1.25, 0.9, R_lip=5.0e-3)
