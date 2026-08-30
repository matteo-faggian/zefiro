"""Il motore da 50 N in geometria implicita.

I difetti che questi test escludono sono TUTTI difetti realmente incontrati
costruendo questo pezzo, e nessuno era visibile nei volumi: il pezzo restava
plausibile mentre era sbagliato. Si vedevano solo in sezione o per via
topologica, ed e' per questo che i controlli sono topologici.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

from zefiro.sdf.engine import costruisci                          # noqa: E402
from zefiro.sdf.meshing import (                                  # noqa: E402
    cavities_are_closed, connected_components, is_watertight, isosurface,
    mesh_volume,
)

#: 0.13 mm, non 0.24. Il setto fra andata e ritorno e' 0.6 mm: a 0.24 mm sono
#: due voxel e mezzo, e la mesh si frantuma - non perche' il pezzo sia
#: sbagliato, ma perche' la griglia non lo risolve. Il fatto che serva 0.13 mm
#: per modellare queste quote e' esso stesso un'informazione di progetto: sono
#: quote al limite anche per il processo.
#: 0.13 mm: sotto questa finezza il pezzo si frantuma nel MODELLO, perche' i
#: setti da 0.6 mm sono meno di quattro voxel. Non e' un difetto del controllo:
#: e' la prova che quelle quote sono al limite anche per il processo.
#: Serve ~4 GB: in un container piccolo la suite intera va in OOM.
PASSO = 1.3e-4

#: Memoria necessaria, misurata: la griglia a 0.13 mm e' ~22 Mvoxel e il
#: costruttore tiene in vita una decina di campi float32, cioe' ~1 GB, piu' la
#: mesh e le trasformate di distanza. Sotto ~6 GB liberi il processo viene
#: ucciso a meta'. Meglio saltare dichiarandolo che finire in OOM e far
#: sembrare che la suite sia passata.
def _memoria_sufficiente(minimo_gb: float = 12.0) -> bool:
    try:
        with open("/proc/meminfo") as fh:
            for riga in fh:
                if riga.startswith("MemAvailable:"):
                    return int(riga.split()[1]) / 1024**2 >= minimo_gb
    except OSError:
        pass
    return True


pesante = pytest.mark.skipif(
    not _memoria_sufficiente(),
    reason="serve almeno 12 GB liberi: la griglia a 0.13 mm che risolve setti "
           "da 0.6 mm non ci sta in meno",
)


@pytest.fixture(scope="module")
def motore():
    from build_50N import CAMERA, GOLA, punto_operativo
    op, params, l0 = punto_operativo()
    m = costruisci(params.derived, params.plug_contour_x, params.plug_contour_r,
                   [CAMERA, GOLA], PASSO, raccordo=5.0e-4)
    return params, l0, m


@pesante
@pytest.mark.slow
def test_e_un_pezzo_solo_e_chiuso(motore):
    """Due difetti in un test, perche' si presentano insieme.

    Un pezzo in DUE componenti non si stampa in un colpo: e' successo davvero,
    perche' il collettore piazzato a x = 0 sconfinava nella piastra e vi
    scavava un anello, staccando un disco da 5432 triangoli. Una superficie
    APERTA non e' stampabile affatto: e' successo dove i due circuiti si
    accavallavano e lasciavano pareti piu' sottili della griglia.

    Nessuno dei due si vedeva nel volume, che restava plausibile.
    """
    _, _, m = motore
    v, f = isosurface(m.solido)
    componenti = connected_components(v, f)
    assert len(componenti) == 1, f"il pezzo e' in {len(componenti)} parti: {componenti[:5]}"
    assert is_watertight(f), "superficie aperta: non stampabile"


@pesante
@pytest.mark.slow
def test_i_canali_non_sbucano_nel_gas(motore):
    """Un canale che sbuca in camera scarica acqua nel gas. E' il modo piu'
    rapido di spegnere il motore, e non si vede da fuori."""
    _, _, m = motore
    chiusa, frazione = cavities_are_closed(m.cavita_gas, m.canali)
    assert chiusa, f"il {frazione:.2%} del volume dei canali sbuca nella cavita' gassosa"


@pesante
@pytest.mark.slow
def test_i_canali_non_svuotano_il_corpo_centrale(motore):
    """La condizione |G - profondita| <= h/2 e' soddisfatta su ENTRAMBI i lati
    della parete: un campo di distanza non ha un verso. Senza il vincolo
    esplicito, i canali si scavavano anche dentro il plug.

    Si verifica che lungo l'asse, dove non ci deve essere nessun canale, il
    materiale sia pieno da x = 0 fino alla punta."""
    params, _, m = motore
    d = params.derived
    x_base = d["L_c"] + d["L_conv"] + params.plug_contour_x[-1]
    xs = np.linspace(0.002, x_base - 0.002, 40)
    punti = np.stack([xs, np.zeros_like(xs), np.zeros_like(xs)], axis=1)
    valori = m.solido.sample(punti)
    assert (valori < 0).all(), (
        f"{(valori >= 0).sum()} punti sull'asse sono vuoti: il corpo centrale "
        "e' stato svuotato dai canali"
    )


@pesante
@pytest.mark.slow
def test_il_volume_dal_campo_e_quello_dalla_mesh(motore):
    """Due misure indipendenti: conteggio sub-voxel del campo contro il
    teorema della divergenza sui triangoli."""
    _, _, m = motore
    v, f = isosurface(m.solido)
    assert mesh_volume(v, f) == pytest.approx(m.solido.volume(), rel=0.02)


@pesante
@pytest.mark.slow
def test_la_sezione_dei_canali_e_quella_progettata(motore):
    """Il volume dei canali deve corrispondere a sezione x lunghezza: se il
    campo li tagliasse, la portata reale non sarebbe quella dimensionata e
    tutto il calcolo termico cadrebbe."""
    _, _, m = motore
    atteso = 0.0
    for c in m.circuiti:
        raggio_medio = 0.012                     # ordine di grandezza, vedi sotto
        lunghezza = (c.x_fine - c.x_inizio) * math.sqrt(
            1.0 + (2.0 * math.pi * raggio_medio / c.passo_elica) ** 2)
        atteso += c.n_canali * c.lato**2 * lunghezza
    v = m.canali.volume()
    # tolleranza larga: i collettori aggiungono volume e il raggio medio e'
    # una stima. Serve a prendere errori di ORDINE DI GRANDEZZA - un canale
    # non scavato, o scavato al doppio - non a validare il terzo decimale.
    assert 0.5 * atteso < v < 3.0 * atteso, (v, atteso)


def test_il_collettore_non_puo_strozzare():
    """Se l'attacco e' piu' stretto della somma dei canali, la portata la
    decide l'attacco e non il progetto. Il costruttore deve dirlo."""
    from build_50N import CAMERA, GOLA
    from zefiro.sdf.engine import PORTA_DIAMETRO, PORTE_PER_COLLETTORE

    # Sono TRE attacchi per collettore, non uno: conta la loro somma.
    sezione = PORTE_PER_COLLETTORE * math.pi / 4.0 * PORTA_DIAMETRO**2
    for c in (CAMERA, GOLA):
        assert sezione >= c.n_canali * c.lato**2, (
            f"attacchi {sezione*1e6:.2f} mm2 contro canali "
            f"{c.n_canali*c.lato**2*1e6:.2f} mm2 sul ramo {c.nome}")


@pesante
@pytest.mark.slow
def test_la_gola_non_e_strozzata_dal_raccordo(motore):
    """Il difetto piu' grave incontrato, e l'unico che nessun altro controllo
    vedeva: con un raccordo da 1.5 mm il materiale entrava nell'anello di gola
    e l'area libera scendeva del 25 %, cioe' 25 % di spinta in meno. Il pezzo
    restava chiuso, in un blocco solo, col circuito sigillato e un volume del
    tutto plausibile.

    Un raccordo e' un'operazione LOCALE solo se i due corpi distano piu' del
    suo raggio. Qui plug e labbro distano 1.78 mm.
    """
    from zefiro.sdf.engine import area_di_gola

    params, _, m = motore
    d = params.derived
    x_lip = d["L_c"] + d["L_conv"]
    misurata, teorica = area_di_gola(m, x_lip, params.plug_contour_r[0], d["R_lip"])
    # la tolleranza e' larga perche' su una griglia grossa il difetto di
    # discretizzazione e' esso stesso qualche percento: serve a prendere una
    # strozzatura VERA (che vale decine di percento), non il terzo decimale.
    assert misurata > 0.90 * teorica, (
        f"area di gola {misurata*1e6:.2f} mm2 contro {teorica*1e6:.2f} attesa "
        f"({misurata/teorica-1:+.1%}): qualcosa ostruisce il passaggio"
    )


@pesante
@pytest.mark.slow
def test_nessun_vuoto_comunica_con_un_altro_che_non_deve(motore):
    """Il controllo che chiude il cerchio in ottica SLM: fra due cavita' che
    non devono comunicare deve restare materiale, e sopra il minimo di processo.

    Due difetti veri trovati cosi', entrambi invisibili in ogni altra verifica:
      * i due attacchi del circuito a U aprivano nello STESSO collettore, quindi
        l'acqua entrava e usciva senza passare per i canali - un cortocircuito;
      * l'attacco d'uscita della camera e quello d'ingresso della gola, a 4 mm
        di distanza con fori da 5.6 mm, si compenetravano: i due rami paralleli
        del raffreddamento diventavano un ramo solo.
    """
    import itertools

    from zefiro.sdf.clearances import MIN_WALL_SLM, rapporto_spessori
    from zefiro.sdf.engine import PORTE_PER_COLLETTORE

    _, _, m = motore
    # connessioni VOLUTE: i fori d'iniezione si aprono in camera (e' la loro
    # funzione), e ogni attacco si apre nel proprio strato
    voluti = {("fori_iniezione", "gas")}
    for c in m.circuiti:
        for k in (0, 1):
            strato = f"{c.nome}_ritorno" if (c.ritorno and k == 1) else f"{c.nome}_andata"
            for j in range(PORTE_PER_COLLETTORE):
                voluti.add(tuple(sorted((f"{c.nome}_attacco_{k}_{j}", strato))))

    nomi = sorted(m.parti)
    coppie = [(a, b, MIN_WALL_SLM) for a, b in itertools.combinations(nomi, 2)
              if tuple(sorted((a, b))) not in voluti]
    assert len(coppie) > 50, "il motore non ha esposto le sue parti"

    esiti = [x for x in rapporto_spessori(m, coppie) if x.minima != float("inf")]
    compenetrano = [x for x in esiti if x.esito == "COMPENETRANO"]
    assert not compenetrano, [
        f"{x.a}<->{x.b}" for x in compenetrano]
    # nessuna parete sotto il minimo, al netto della tolleranza di griglia
    sotto = [x for x in esiti if x.esito == "SOTTO IL MINIMO"]
    assert not sotto, [f"{x.a}<->{x.b} = {x.minima*1e3:.3f} mm" for x in sotto]


@pesante
@pytest.mark.slow
def test_niente_polvere_intrappolata_ne_frammenti(motore):
    """Due difetti SLM speculari: una sacca di vuoto chiusa e' polvere che non
    esce; un'isola di materiale e' un frammento che si stacca. Nel circuito di
    raffreddamento entrambi finiscono per ostruire un canale da 0.6 mm."""
    from zefiro.sdf.meshing import connected_components, isosurface
    from zefiro.sdf.printability import polvere_evacuabile

    _, _, m = motore
    ok, _, sacche = polvere_evacuabile(m.canali, m.solido, m.grid)
    assert ok, f"{sacche} sacche di vuoto chiuse: polvere che non esce"
    v, f = isosurface(m.solido)
    assert len(connected_components(v, f)) == 1
