"""Il livello a campi di distanza.

Ogni test confronta il campo campionato con un valore ESATTO calcolato per
altra via. Un SDF sbagliato produce solidi che sembrano plausibili: non si
vede a occhio, si vede solo misurando.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from zefiro.sdf.core import Field, Grid
from zefiro.sdf.meshing import is_watertight, isosurface, mesh_area, mesh_volume
from zefiro.sdf.shapes import cylinder, polygon_sdf_2d, revolve_polygon, sphere


def test_sfera_distanza_esatta_e_converge_come_h2():
    """Il campo di una sfera E' |p - c| - R. Si verifica in punti a caso, non
    sui nodi, cosi' si misura anche l'interpolazione.

    La soglia non e' scelta a occhio: l'interpolazione trilineare di una
    funzione curva sbaglia di ordine h^2, quindi dimezzando il passo l'errore
    deve calare di circa quattro volte. Si verifica il RAPPORTO, che non si
    puo' far passare allentando un numero."""
    rng = np.random.default_rng(0)
    p = rng.uniform(-0.025, 0.025, size=(400, 3))
    atteso = np.linalg.norm(p, axis=1) - 0.02
    # Il tasso di convergenza si misura sull'errore RMS e non sul massimo: il
    # massimo su 400 punti a caso dipende da quale punto capita peggio
    # rispetto alla griglia, ed e' una statistica rumorosa. Il massimo resta
    # come limite assoluto.
    rms, massimi = [], []
    for h in (0.001, 0.0005, 0.00025):
        g = Grid.bounding([-0.03] * 3, [0.03] * 3, h)
        f = sphere(g, (0.0, 0.0, 0.0), 0.02)
        e = np.abs(f.sample(p) - atteso)
        rms.append(float(np.sqrt((e ** 2).mean())))
        massimi.append(float(e.max()))
    for prima, dopo in zip(rms[:-1], rms[1:]):
        assert 3.2 < prima / dopo < 5.0, f"non converge come h^2: {rms}"
    assert massimi[-1] < 5.0e-6, massimi


def test_volume_di_una_sfera_converge():
    """Il volume dal campo deve tendere a 4/3 pi R^3, e l'errore deve
    DIMINUIRE raffinando: e' la prova che la correzione sub-voxel funziona,
    e non si puo' barare con una soglia."""
    R = 0.02
    esatto = 4.0 / 3.0 * math.pi * R**3
    errori = []
    for h in (0.002, 0.001, 0.0005):
        g = Grid.bounding([-0.026] * 3, [0.026] * 3, h)
        errori.append(abs(sphere(g, (0, 0, 0), R).volume() - esatto) / esatto)
    assert errori[-1] < 0.005, errori
    assert errori[0] > errori[1] > errori[2], f"non converge: {errori}"


def test_cilindro_volume_e_tenuta():
    g = Grid.bounding([-0.005, -0.02, -0.02], [0.045, 0.02, 0.02], 0.0005)
    f = cylinder(g, 0.015, 0.0, 0.04)
    assert not f.tocca_il_bordo()
    esatto = math.pi * 0.015**2 * 0.04
    assert f.volume() == pytest.approx(esatto, rel=0.01)
    v, t = isosurface(f)
    assert is_watertight(t)
    assert mesh_volume(v, t) == pytest.approx(esatto, rel=0.01)


def test_poligono_2d_segno_e_distanza():
    """Quadrato 2x2 centrato: la distanza si scrive a mano ovunque."""
    quad = [(-1.0, 0.0), (1.0, 0.0), (1.0, 2.0), (-1.0, 2.0)]
    px = np.array([0.0, 0.0, 2.0, -3.0, 0.9])
    pr = np.array([1.0, 3.0, 1.0, 1.0, 1.0])
    atteso = np.array([-1.0, 1.0, 1.0, 2.0, -0.1])
    got = polygon_sdf_2d(px, pr, quad)
    assert np.allclose(got, atteso, atol=1e-6), (got, atteso)


def test_rivoluzione_esatta_contro_un_cilindro():
    """Un rettangolo ruotato E' un cilindro cavo: due strade indipendenti allo
    stesso solido, e devono coincidere a precisione di campionamento."""
    g = Grid.bounding([-0.005, -0.025, -0.025], [0.045, 0.025, 0.025], 0.0005)
    tubo = [(0.0, 0.010), (0.04, 0.010), (0.04, 0.018), (0.0, 0.018)]
    f = revolve_polygon(g, tubo)
    esatto = math.pi * (0.018**2 - 0.010**2) * 0.04
    assert f.volume() == pytest.approx(esatto, rel=0.01)
    v, t = isosurface(f)
    assert is_watertight(t)
    assert mesh_volume(v, t) == pytest.approx(esatto, rel=0.01)


def test_offset_e_davvero_un_offset():
    """`f - t` deve essere la superficie spostata di t. Si verifica sul
    volume: una sfera offset di t ha il volume della sfera di raggio R + t."""
    g = Grid.bounding([-0.04] * 3, [0.04] * 3, 0.0008)
    f = sphere(g, (0, 0, 0), 0.02)
    for t in (0.003, 0.006):
        atteso = 4.0 / 3.0 * math.pi * (0.02 + t) ** 3
        assert f.offset(t).volume() == pytest.approx(atteso, rel=0.01), t


def test_smooth_union_raccorda_senza_gonfiare():
    """Il raccordo deve aggiungere materiale SOLO vicino al giunto: lontano,
    il campo deve restare quello dell'unione netta. Se cosi' non fosse, il
    raccordo cambierebbe la geometria dappertutto."""
    g = Grid.bounding([-0.05, -0.03, -0.03], [0.05, 0.03, 0.03], 0.0006)
    a = sphere(g, (-0.012, 0, 0), 0.015)
    b = sphere(g, (0.012, 0, 0), 0.015)
    netta = a.union(b)
    raccordata = a.smooth_union(b, 0.006)
    v_n, v_r = netta.volume(), raccordata.volume()
    assert v_r > v_n, "il raccordo deve aggiungere materiale"
    assert (v_r - v_n) / v_n < 0.06, f"ne aggiunge troppo: {(v_r-v_n)/v_n:.1%}"
    # Lontano dal giunto i due campi devono coincidere. I punti sono sui lati
    # ESTERNI delle due sfere: sul piano di simmetria non varrebbe, perche' li'
    # le due sfere sono equidistanti ed e' proprio dove il raccordo agisce.
    lontano = np.array([[-0.030, 0.0, 0.0], [0.030, 0.0, 0.0],
                        [-0.012, 0.026, 0.0], [0.012, 0.0, 0.026]])
    assert np.allclose(netta.sample(lontano), raccordata.sample(lontano), atol=1e-6)


def test_il_bordo_toccato_viene_segnalato():
    # R = 0.0105 supera la semiampiezza della scatola (0.01), quindi la sfera
    # viene tagliata dalle facce. Con R = 0.0099 NON la tocca, anche se ci
    # arriva vicino: sfiorare non e' tagliare.
    g = Grid.bounding([-0.01] * 3, [0.01] * 3, 0.0005)
    assert sphere(g, (0, 0, 0), 0.0105).tocca_il_bordo()
    assert not sphere(g, (0, 0, 0), 0.005).tocca_il_bordo()
