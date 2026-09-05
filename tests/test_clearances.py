"""Spessori di parete fra vuoti, isole e sacche.

In SLM la domanda che decide se un pezzo si stampa o perde non e' "e' chiuso?"
ma "quanto materiale resta fra due cavita' che non devono comunicare?". Un
pezzo puo' essere chiuso, in un blocco solo, col circuito sigillato, e avere un
setto da 0.15 mm che in stampa esce poroso.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from zefiro.sdf.clearances import distanza_minima
from zefiro.sdf.core import Field, Grid
from zefiro.sdf.printability import isole_di_materiale, riempi_sacche_chiuse
from zefiro.sdf.shapes import sphere


def test_distanza_minima_e_esatta_su_due_sfere():
    """Due sfere di raggio R a distanza D fra i centri distano D - 2R.
    Se la misura fosse sbagliata, ogni verdetto sugli spessori lo sarebbe."""
    g = Grid.bounding([-0.03, -0.015, -0.015], [0.03, 0.015, 0.015], 0.0004)
    a = sphere(g, (-0.012, 0, 0), 0.006)
    b = sphere(g, (+0.012, 0, 0), 0.006)
    atteso = 0.024 - 2 * 0.006
    d = distanza_minima(np.asarray(a.a) < 0, np.asarray(b.a) < 0, g.spacing)
    assert d == pytest.approx(atteso, abs=2.5 * g.spacing), (d, atteso)


def test_due_corpi_che_si_toccano_danno_zero():
    g = Grid.bounding([-0.02, -0.012, -0.012], [0.02, 0.012, 0.012], 0.0004)
    a = sphere(g, (-0.004, 0, 0), 0.006)
    b = sphere(g, (+0.004, 0, 0), 0.006)
    assert distanza_minima(np.asarray(a.a) < 0, np.asarray(b.a) < 0, g.spacing) == 0.0


def test_le_isole_di_materiale_vengono_trovate_e_tolte():
    """Un frammento staccato in SLM non e' un dettaglio: e' materiale
    sinterizzato che si stacca e, dentro un circuito, lo ostruisce."""
    g = Grid.bounding([-0.02] * 3, [0.02] * 3, 0.0005)
    grosso = sphere(g, (-0.008, 0, 0), 0.006)
    scheggia = sphere(g, (0.010, 0, 0), 0.002)
    insieme = grosso.union(scheggia)
    v_prima = insieme.volume()
    pulito, n, vol = isole_di_materiale(insieme, volume_minimo=1.0e-6)
    assert n == 1
    assert vol == pytest.approx(4 / 3 * math.pi * 0.002**3, rel=0.25)
    assert pulito.volume() < v_prima * 0.98


def test_le_sacche_chiuse_vengono_riempite_solo_se_piccole():
    """Una bolla dentro il materiale e' polvere che non esce. Le minuscole
    sono artefatti degli spigoli e si riempiono; una grande e' un errore di
    progetto e deve RESTARE, per essere vista."""
    g = Grid.bounding([-0.02] * 3, [0.02] * 3, 0.0004)
    blocco = sphere(g, (0, 0, 0), 0.015)

    piccola = blocco.difference(sphere(g, (0.005, 0, 0), 0.0012))
    _, n, vol, rimaste = riempi_sacche_chiuse(piccola, volume_massimo=(3.0e-3) ** 3)
    assert n == 1 and rimaste == 0

    grande = blocco.difference(sphere(g, (0.005, 0, 0), 0.005))
    _, n2, _, rimaste2 = riempi_sacche_chiuse(grande, volume_massimo=(1.0e-3) ** 3)
    assert n2 == 0 and rimaste2 == 1, "una sacca grande non va nascosta"
