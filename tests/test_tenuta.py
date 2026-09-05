"""Prova di tenuta numerica: su geometrie in cui la risposta si sa a mano.

Un controllo che non e' mai stato visto fallire su un difetto costruito
apposta non e' un controllo: e' una funzione che ritorna sempre True. Qui i
difetti si costruiscono di proposito, uno alla volta.
"""
from __future__ import annotations

import numpy as np
import pytest

from zefiro.sdf.core import Field, Grid
from zefiro.sdf.shapes import cylinder, sphere
from zefiro.sdf.tenuta import Sonda, confronta, verifica_tenuta

PASSO = 5.0e-4


def blocco(grid, lo, hi):
    """Parallelepipedo come campo di distanza (esatto fuori, sottostimato
    dentro: basta, perche' qui serve solo il segno)."""
    x, y, z = grid.coords()
    d = None
    for c, a, b in zip((x, y, z), lo, hi):
        e = np.maximum(a - c, c - b)
        d = e if d is None else np.maximum(d, e)
    return Field(grid, d.astype(np.float32) + 0.0 * (x + y + z).astype(np.float32))


#: Tutte le forme di `zefiro.sdf.shapes` sono di rivoluzione attorno all'ASSE X,
#: come il motore. Le cavita' di prova stanno quindi su y = z = 0.
@pytest.fixture(scope="module")
def griglia():
    return Grid.bounding((-0.01, -0.015, -0.015), (0.06, 0.015, 0.015), PASSO,
                         margin=0.002)


def _due_cavita(grid, foro_tra=False, foro_fuori=False, terza_cavita=False):
    """Un blocco con due sfere vuote dentro. Le opzioni aggiungono i difetti."""
    solido = blocco(grid, (0.0, -0.010, -0.010), (0.050, 0.010, 0.010))
    a = sphere(grid, (0.012, 0.0, 0.0), 0.006)
    b = sphere(grid, (0.038, 0.0, 0.0), 0.006)
    vuoti = a.union(b)
    if terza_cavita:
        vuoti = vuoti.union(sphere(grid, (0.025, 0.0, 0.0065), 0.002))
    if foro_tra:
        vuoti = vuoti.union(cylinder(grid, 0.0008, 0.010, 0.040))
    if foro_fuori:
        # un canale radiale che sbuca dalla faccia y = -10 mm: il buco nella
        # pelle vero e proprio.
        x, y, z = grid.coords()
        r = np.sqrt((x - 0.012) ** 2 + z ** 2)
        d = np.maximum(r - 0.0008, np.maximum(-0.020 - y, y - 0.0))
        vuoti = vuoti.union(Field(grid, np.broadcast_to(
            d, grid.shape).astype(np.float32).copy()))
    return solido.difference(vuoti)


SONDE = [
    Sonda("cavita_A", (0.012, 0.0, 0.0)),
    Sonda("cavita_B", (0.038, 0.0, 0.0)),
    Sonda("esterno", (-0.008, 0.014, 0.014)),
]


def test_due_cavita_separate_sono_stagne(griglia):
    r = verifica_tenuta(_due_cavita(griglia), griglia, SONDE)
    assert r.stagno, r.riassunto()
    assert len({r.etichette[s.nome] for s in SONDE}) == 3
    # e i volumi tornano: due sfere da 6 mm di raggio
    atteso = 4.0 / 3.0 * np.pi * 0.006 ** 3
    for nome in ("cavita_A", "cavita_B"):
        assert r.volumi[r.etichette[nome]] == pytest.approx(atteso, rel=0.05)


def test_un_foro_fra_le_due_cavita_e_una_perdita(griglia):
    """Il difetto che la caratteristica di Eulero NON distingue: un manico in
    piu' e un manico spostato danno lo stesso genere."""
    r = verifica_tenuta(_due_cavita(griglia, foro_tra=True), griglia, SONDE)
    assert not r.stagno
    assert ("cavita_A", "cavita_B") in r.inattese
    assert ("cavita_A", "esterno") not in r.comunicazioni


def test_un_foro_nella_pelle_mette_in_comunicazione_con_l_esterno(griglia):
    r = verifica_tenuta(_due_cavita(griglia, foro_fuori=True), griglia, SONDE)
    assert not r.stagno
    assert ("cavita_A", "esterno") in r.inattese
    assert ("cavita_B", "esterno") not in r.comunicazioni


def test_una_comunicazione_dichiarata_non_e_una_perdita(griglia):
    """I fori d'iniezione mettono in comunicazione GPL e aria di proposito.
    Dichiararli e' obbligatorio, altrimenti il motore funzionante risulta
    difettoso e si smette di guardare i risultati."""
    r = verifica_tenuta(_due_cavita(griglia, foro_tra=True), griglia, SONDE,
                        attese={("cavita_A", "cavita_B")})
    assert r.stagno, r.riassunto()
    assert not r.inattese


def test_una_comunicazione_voluta_che_manca_e_un_difetto(griglia):
    """Un motore in cui il GPL non arriva in camera e' rotto quanto uno che
    perde, e nessun controllo di tenuta classico se ne accorge."""
    r = verifica_tenuta(_due_cavita(griglia), griglia, SONDE,
                        attese={("cavita_A", "cavita_B")})
    assert not r.stagno
    assert ("cavita_A", "cavita_B") in r.mancate


def test_un_vuoto_che_nessuno_rivendica_viene_segnalato(griglia):
    """Una sacca chiusa non dichiarata: polvere che non esce piu'."""
    r = verifica_tenuta(_due_cavita(griglia, terza_cavita=True), griglia, SONDE)
    assert not r.stagno
    assert len(r.non_rivendicate) == 1
    assert r.volume_non_rivendicato == pytest.approx(
        4.0 / 3.0 * np.pi * 0.002 ** 3, rel=0.2)


def test_una_sonda_nel_pieno_e_un_errore_non_un_successo(griglia):
    """Se la sonda finisce nel materiale la prova non ha misurato niente.
    Ritornare 'stagno' sarebbe il modo piu' facile di superare il controllo
    per sbaglio."""
    r = verifica_tenuta(_due_cavita(griglia), griglia,
                        [Sonda("dentro_il_metallo", (0.025, 0.0, 0.0))])
    assert not r.stagno
    assert r.sonde_fuori == ["dentro_il_metallo"]


def test_una_sonda_fuori_dalla_griglia_e_un_errore_esplicito(griglia):
    with pytest.raises(ValueError, match="fuori dalla griglia"):
        verifica_tenuta(_due_cavita(griglia), griglia, [Sonda("x", (9.0, 0.0, 0.0))])


def test_i_voxel_che_si_toccano_per_uno_spigolo_non_comunicano():
    """Perche' la connettivita' e' a 6 facce e non a 26.

    Due cavita' separate da una parete spessa un voxel, sfalsate di un voxel in
    diagonale: con la connettivita' a 26 risulterebbero comunicanti attraverso
    uno spigolo, cioe' attraverso una parete che c'e'. E' un falso allarme, e i
    falsi allarmi su un controllo di tenuta sono costosi quanto i mancati.
    """
    g = Grid(origin=(0.0, 0.0, 0.0), spacing=1.0e-3, shape=(9, 9, 9))
    pieno = np.full(g.shape, -1.0e-3, dtype=np.float32)
    pieno[2, 2, 2] = 1.0e-3
    pieno[3, 3, 2] = 1.0e-3          # tocca il primo solo per uno spigolo
    campo = Field(g, pieno)
    r = verifica_tenuta(campo, g, [Sonda("a", (0.002, 0.002, 0.002)),
                                   Sonda("b", (0.003, 0.003, 0.002))])
    assert r.comunicazioni == []
    assert r.n_componenti == 2


def test_il_confronto_fra_risoluzioni_vede_un_risultato_instabile(griglia):
    """Un verdetto che cambia col passo non e' un verdetto. Questo e'
    esattamente il caso in cui ci si e' gia' trovati con il genere del motore:
    setti da 0.6 mm su una griglia da 0.15 mm.
    """
    # Griglie costruite a mano e non con `bounding`: il punto del test e' che
    # il setto cada FRA due nodi della griglia grossa e SOPRA un nodo di quella
    # fine, e questo dipende dall'allineamento, non solo dal passo.
    fine = Grid(origin=(-0.010, -0.015, -0.015), spacing=2.5e-4,
                shape=(281, 121, 121))
    grosso = Grid(origin=(-0.010, -0.015, -0.015), spacing=1.0e-3,
                  shape=(71, 31, 31))

    def con_setto_sottile(g):
        """Due cavita' separate da 0.4 mm di pieno, fra x = 24.5 e 24.9 mm.
        A 0.25 mm il setto contiene nodi ed esiste; a 1.0 mm i nodi cadono a
        24.0 e 25.0 mm, entrambi nel vuoto, e il setto SPARISCE."""
        solido = blocco(g, (0.0, -0.010, -0.010), (0.050, 0.010, 0.010))
        a = cylinder(g, 0.004, 0.010, 0.0245)
        b = cylinder(g, 0.004, 0.0249, 0.040)
        return solido.difference(a.union(b))

    sonde = [Sonda("A", (0.015, 0.0, 0.0)), Sonda("B", (0.035, 0.0, 0.0)),
             Sonda("esterno", (-0.008, 0.014, 0.014))]
    rf = verifica_tenuta(con_setto_sottile(fine), fine, sonde)
    rg = verifica_tenuta(con_setto_sottile(grosso), grosso, sonde)
    assert ("A", "B") not in rf.comunicazioni
    assert ("A", "B") in rg.comunicazioni, (
        "il setto sotto-risolto doveva sparire: se non sparisce il test non "
        "sta piu' provando quello che dice di provare")
    assert confronta(rf, rg), "il confronto deve accorgersi della differenza"
