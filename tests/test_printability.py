"""Stampabilita': sbalzi, supporti impossibili, evacuazione della polvere.

Il test piu' importante di questo file non riguarda il motore: riguarda il
VERSO delle normali. Se marching cubes le restituisse entranti invece che
uscenti, tutta l'analisi degli sbalzi sarebbe capovolta - direbbe "tutto a
posto" esattamente dove il pezzo non si stampa. E' un'assunzione che non si
puo' lasciare implicita, quindi si verifica su una sfera, dove la risposta si
scrive a mano.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from zefiro.sdf.core import Grid
from zefiro.sdf.meshing import isosurface
from zefiro.sdf.printability import (
    SELF_SUPPORT_ANGLE_DEG,
    analizza_sbalzi,
    passo_elica_minimo,
)
from zefiro.sdf.shapes import cylinder, revolve_polygon, sphere


def test_le_normali_puntano_fuori_dal_materiale():
    """Su una sfera centrata nell'origine, una normale uscente soddisfa
    n . p > 0 in ogni punto. Se fallisse, `analizza_sbalzi` scambierebbe
    soffitti e pavimenti."""
    g = Grid.bounding([-0.03] * 3, [0.03] * 3, 0.0008)
    v, f = isosurface(sphere(g, (0.0, 0.0, 0.0), 0.02))
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    n = np.cross(b - a, c - a)
    n /= np.linalg.norm(n, axis=1, keepdims=True)
    centri = v[f].mean(axis=1)
    assert ((n * centri).sum(axis=1) > 0).mean() > 0.999


def test_un_cilindro_verticale_non_ha_sbalzi():
    """Costruito lungo il proprio asse, un cilindro ha solo pareti verticali e
    due facce piane: quella inferiore poggia sulla piastra, quella superiore
    guarda in alto. Niente da supportare."""
    g = Grid.bounding([-0.004, -0.02, -0.02], [0.044, 0.02, 0.02], 0.0005)
    v, f = isosurface(cylinder(g, 0.015, 0.0, 0.04))
    r = analizza_sbalzi(v, f, (1.0, 0.0, 0.0))
    assert r.area_da_supportare < 0.02 * r.area_totale, r.area_da_supportare
    assert r.area_sulla_piastra > 0.0, "la faccia appoggiata non e' stata riconosciuta"


def test_una_sporgenza_orizzontale_viene_vista():
    """Il controllo deve trovare ROBA VERA: un disco che sporge da un albero
    ha una faccia inferiore piatta, ed e' esattamente cio' che non si stampa.
    Un test che passa solo su geometrie sane non dimostra niente."""
    g = Grid.bounding([-0.004, -0.025, -0.025], [0.044, 0.025, 0.025], 0.0005)
    profilo = [(0.0, 0.0), (0.0, 0.006), (0.020, 0.006), (0.020, 0.020),
               (0.024, 0.020), (0.024, 0.006), (0.040, 0.006), (0.040, 0.0)]
    v, f = isosurface(revolve_polygon(g, profilo))
    r = analizza_sbalzi(v, f, (1.0, 0.0, 0.0))
    atteso = math.pi * (0.020**2 - 0.006**2)      # corona inferiore del disco
    assert r.area_da_supportare > 0.7 * atteso, (r.area_da_supportare, atteso)
    assert r.angolo_minimo < 5.0


def test_passo_minimo_e_reversibile():
    """`passo_elica_minimo` e la sua inversa devono chiudere il cerchio.

    Non e' pedanteria: la prima stesura scriveva 2 pi R / tan(alpha) invece di
    2 pi R tan(alpha). A 45 gradi tan vale 1 e le due COINCIDONO, quindi il
    caso di prova piu' ovvio non le distingue. Qui si verifica a 30 e 60."""
    R = 0.0132
    for angolo in (30.0, 45.0, 60.0):
        p = passo_elica_minimo(R, angolo)
        indietro = math.degrees(math.atan(p / (2.0 * math.pi * R)))
        assert indietro == pytest.approx(angolo, abs=1e-9), (angolo, indietro)
    assert passo_elica_minimo(R, 45.0) == pytest.approx(2.0 * math.pi * R)
    # piu' margine si vuole, piu' passo serve: se fosse invertita, calerebbe
    assert passo_elica_minimo(R, 60.0) > passo_elica_minimo(R, 45.0)


def test_la_soglia_e_dichiarata():
    assert 30.0 <= SELF_SUPPORT_ANGLE_DEG <= 50.0
