"""La verifica topologica della pelle.

Il punto di questi test e' che ciascuno usa una forma di cui il GENERE si
conosce a priori: sfera 0, toro 1, due tori 2. Se il calcolo non ritrova quei
numeri interi, non serve a niente su una geometria vera.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from zefiro.sdf.core import Field, Grid
from zefiro.sdf.meshing import isosurface
from zefiro.sdf.shapes import cylinder, revolve_polygon, sphere
from zefiro.sdf.topology import analizza, salda_e_pulisci


def _toro(grid, R, r, xp=np):
    """Toro attorno all'asse x: la rivoluzione di un cerchio."""
    n = 96
    t = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
    poligono = [(r * math.cos(a), R + r * math.sin(a)) for a in t]
    return revolve_polygon(grid, poligono)


def test_la_sfera_ha_genere_zero():
    g = Grid.bounding([-0.03] * 3, [0.03] * 3, 0.0008)
    v, f = salda_e_pulisci(*isosurface(sphere(g, (0, 0, 0), 0.02)))
    r = analizza(v, f, genere_atteso=0)
    assert r.euler == 2, r.euler
    assert r.genere == pytest.approx(0.0)
    assert r.sana, r.note


def test_il_toro_ha_genere_uno():
    """Il controllo decisivo: un toro E' un buco passante. Se il calcolo non
    lo vede, non vedra' nemmeno un canale che ha sfondato."""
    g = Grid.bounding([-0.012, -0.035, -0.035], [0.012, 0.035, 0.035], 0.0006)
    v, f = salda_e_pulisci(*isosurface(_toro(g, 0.022, 0.007)))
    r = analizza(v, f, genere_atteso=1)
    assert r.euler == 0, r.euler
    assert r.genere == pytest.approx(1.0)
    assert r.sana, r.note


def test_due_buchi_danno_genere_due():
    """Due tori separati: due componenti, genere totale 2. Verifica che la
    formula regga anche quando le componenti sono piu' di una."""
    g = Grid.bounding([-0.045, -0.035, -0.035], [0.045, 0.035, 0.035], 0.0008)
    # due tori a quote x diverse: si costruiscono traslando il POLIGONO, non
    # arrotolando l'array, che avvolgerebbe il dominio su se stesso
    def toro_a(x0):
        n = 96
        t = np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)
        return revolve_polygon(g, [(x0 + 0.007 * math.cos(u),
                                    0.022 + 0.007 * math.sin(u)) for u in t])
    v, f = isosurface(toro_a(-0.030).union(toro_a(0.030)))
    v, f = salda_e_pulisci(v, f)
    r = analizza(v, f)
    assert r.n_componenti == 2
    assert r.genere == pytest.approx(2.0)


def test_un_buco_non_voluto_viene_denunciato():
    """La prova che conta per il motore: si prende un cilindro pieno (genere 0)
    e ci si fa passare un foro trasversale, come farebbe un canale che sfonda.
    Il controllo deve dire che il genere e' salito di uno E dire di quanto."""
    g = Grid.bounding([-0.004, -0.02, -0.02], [0.034, 0.02, 0.02], 0.0005)
    pieno = cylinder(g, 0.014, 0.0, 0.030)
    v, f = salda_e_pulisci(*isosurface(pieno))
    assert analizza(v, f, genere_atteso=0).sana

    # foro passante lungo y: il campo di un cilindro ruotato
    X, Y, Z = g.coords()
    d = np.sqrt((X - 0.015) ** 2 + Z**2) - 0.004
    foro = Field(g, np.broadcast_to(d, g.shape).astype(np.float32).copy())
    bucato = pieno.difference(foro)
    v2, f2 = salda_e_pulisci(*isosurface(bucato))
    r = analizza(v2, f2, genere_atteso=0)
    assert r.genere == pytest.approx(1.0), r.genere
    assert not r.sana
    assert any("GENERE" in n for n in r.note), r.note


def test_una_superficie_aperta_viene_vista():
    """Se la pelle e' aperta, il genere non ha senso e va detto invece di
    stampare un numero."""
    g = Grid.bounding([-0.01] * 3, [0.01] * 3, 0.0005)
    v, f = isosurface(sphere(g, (0, 0, 0), 0.0105))   # tagliata dal bordo
    r = analizza(v, f)
    assert r.spigoli_di_bordo > 0
    assert not r.sana
    assert any("APERTA" in n for n in r.note)


def test_le_facce_degeneri_vengono_contate():
    v = np.array([[0.0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    f = np.array([[0, 1, 2], [0, 1, 1]])     # la seconda ha area nulla
    assert analizza(v, f).facce_degeneri == 1


def test_saldare_e_pulire_non_apre_buchi():
    """La riparazione deve TOGLIERE i degeneri senza cambiare la topologia.

    Il caso di prova e' costruito a mano: un tetraedro chiuso a cui si aggiunge
    un vertice DUPLICATO e due facce che lo usano. E' la forma esatta del
    difetto che marching cubes produceva prima di `allow_degenerate=False`.
    Una pulizia che apra un buco e' peggio del problema che risolve, quindi si
    verifica che il genere non cambi e che non compaiano bordi."""
    v = np.array([[0.0, 0, 0], [1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0],
                  [1.0, 0, 0]])                      # il 4 duplica il vertice 1
    f = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3],
                  [1, 4, 2]])                        # faccia di area nulla
    prima = analizza(v, f)
    assert prima.facce_degeneri >= 1
    v2, f2 = salda_e_pulisci(v, f)
    dopo = analizza(v2, f2)
    assert dopo.facce_degeneri == 0
    assert dopo.spigoli_di_bordo == 0
    assert dopo.euler == 2 and dopo.genere == pytest.approx(0.0)
    assert dopo.n_facce == 4
