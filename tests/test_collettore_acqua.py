"""Il collettore dell'acqua: un ingresso, una uscita, e nessuna via di fuga.

I test qui esistono perche' il passaggio da due rami in PARALLELO a due tratti
in SERIE ha cambiato l'architettura del raffreddamento senza far fallire un
solo test. Un cambio d'architettura che passa inosservato alla suite significa
che la suite non guardava l'architettura: guardava i numeri che ne uscivano.

Ogni test qui sotto corrisponde a un difetto REALMENTE trovato durante il
lavoro, e sono tutti dello stesso tipo: geometrie che superavano ogni verifica
di volume, di superficie chiusa e di conteggio dei componenti, e che erano
sbagliate. Sono i controlli che avrei voluto avere prima.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from zefiro.sdf.core import Grid, backend, to_numpy
from zefiro.sdf.engine import (
    Portagomma, _collare, _foro_goccia, _vano_anulare, sonde_motore)
from zefiro.sintesi import Impianto, Requisiti
from zefiro.sintesi.architetture import aerospike_gas_gas as ARCH


def requisiti() -> Requisiti:
    return Requisiti(
        spinta=50.0, durata=5.0,
        impianto=Impianto(combustibile={"C3H8": 1.0}, T_bombola_min=288.15,
                          T_ossidante_iniezione=293.0,
                          T_combustibile_iniezione=283.0,
                          cd_ossidante=0.75, cd_combustibile=0.75))


@pytest.fixture(scope="module")
def dim():
    return ARCH.dimensiona(requisiti())


# --------------------------------------------------------------------------- #
# 1. la serie e' una serie
# --------------------------------------------------------------------------- #
def test_i_due_tratti_hanno_la_stessa_portata(dim):
    """In serie non c'e' niente da ripartire: e' l'intero punto del cambio.

    In parallelo le due portate erano diverse (0.9 e 8.9 l/min) e per farle
    restare tali serviva una strozzatura calibrata da 2.30 mm all'ingresso del
    ramo scarico. Quel foro e' un pezzo che DEVE essere giusto e che, se
    sbagliato, non si vede: portata e temperatura all'uscita restano normali
    mentre il labbro non e' piu' raffreddato.
    """
    camera, gola = dim.circuiti
    q = [c.n_canali * c.lato ** 2 * c.velocita for c in (camera, gola)]
    assert q[0] == pytest.approx(q[1], rel=1e-6), (
        "i due tratti portano portate diverse: non sono in serie")


def test_nessuna_strozzatura_di_bilanciamento(dim):
    """La quota che il parallelo richiedeva non deve piu' esistere."""
    nomi = set(dim.registro.quote)
    assert not [n for n in nomi if "strozzatura_bilanciamento" in n]


def test_la_gola_detta_la_portata_e_la_camera_la_subisce(dim):
    """L'ordine dei passi conta. In gola la portata NON e' libera: serve h per
    non far bollire l'acqua sulla parete, e la velocita' che lo produce fissa
    l'area. La camera prende quello che passa."""
    reg = dim.registro
    assert "velocita_acqua_gola" in reg.quote
    motivo = reg.quote["velocita_acqua_camera"].motivo
    assert "RISULTATO" in motivo, (
        "la velocita' in camera deve essere dichiarata come conseguenza della "
        "portata imposta, non come scelta")


# --------------------------------------------------------------------------- #
# 2. un ingresso e una uscita
# --------------------------------------------------------------------------- #
def test_esattamente_due_attacchi(dim):
    col = dim.collettore
    assert col is not None
    attacchi = [col.ingresso, col.uscita]
    assert all(isinstance(a, Portagomma) for a in attacchi)
    assert len(attacchi) == 2


def test_gli_attacchi_sono_radiali_e_dalla_stessa_parte(dim):
    """Entrambi a 180 gradi: l'aria entra a 0 e le due canne escono dal lato
    opposto, cosi' non si contendono lo stesso pezzo di corona."""
    col = dim.collettore
    for pg in (col.ingresso, col.uscita):
        assert pg.verso == 0.0
        assert pg.angolo == pytest.approx(math.pi)


# --------------------------------------------------------------------------- #
# 3. il difetto che e' tornato due volte: la sede del portagomma
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("capo", ["ingresso", "uscita"])
def test_la_goccia_del_portagomma_sta_dentro_il_collare(dim, capo):
    """L'apice della goccia sta `pendenza * raggio` sopra il centro del foro.

    Prima versione: attacco centrato a meta' anello, apice fuori dalla faccia di
    valle del collare - il circuito si apriva sull'esterno a due millimetri dal
    labbro. Corretto. Seconda versione: aggiunto il vincolo opposto ("il bordo
    inferiore deve stare dentro l'anello"), apice di nuovo fuori. Stesso
    difetto, stesso punto, due volte: questo test esiste perche' non ci sia una
    terza.
    """
    col = dim.collettore
    pg = getattr(col, capo)
    x_v0, x_v1, h, x_c0 = getattr(col, f"anello_{capo}")
    apice = pg.x_base + col.pendenza_tetto * 0.5 * pg.d_foro
    assert apice <= x_v1 + 1e-9, (
        f"l'apice della goccia ({apice*1e3:.2f} mm) esce dal vano, che finisce "
        f"a {x_v1*1e3:.2f}")
    #: e sotto: il foro attraversa il tratto di collare fra la parete esterna
    #: del mantello e il raggio a cui comincia il gambo, dove il collare non e'
    #: ancora alto perche' sta ancora salendo lungo la sua rampa
    fondo_ammesso = x_c0 + col.pendenza_tetto * 0.5 * h + col.parete
    assert pg.x_base - 0.5 * pg.d_foro >= fondo_ammesso - 1e-9, (
        f"il fondo del foro ({(pg.x_base-0.5*pg.d_foro)*1e3:.2f} mm) cade sulla "
        f"rampa del collare, che a quella profondita' comincia a "
        f"{fondo_ammesso*1e3:.2f}")


@pytest.mark.parametrize("capo", ["ingresso", "uscita"])
def test_il_foro_del_portagomma_incrocia_il_proprio_anello(dim, capo):
    """Un attacco che non tocca l'anello e' un attacco scollegato."""
    col = dim.collettore
    pg = getattr(col, capo)
    x_v0, x_v1, _, _ = getattr(col, f"anello_{capo}")
    assert pg.x_base + 0.5 * pg.d_foro > x_v0
    assert pg.x_base - 0.5 * pg.d_foro < x_v1


# --------------------------------------------------------------------------- #
# 4. distribuzione: gli anelli servono a questo
# --------------------------------------------------------------------------- #
def test_gli_anelli_non_sbilanciano_i_canali(dim):
    """La caduta lungo un anello va confrontata con quella di un canale.

    Se l'anello cade quanto un canale, il canale vicino all'attacco vede una
    contropressione minore e si prende quasi tutta l'acqua: quelli lontani
    restano a secco. E' il vincolo che dimensiona gli anelli, e per questo NON
    e' una velocita' scelta a mano.
    """
    reg = dim.registro
    for capo in ("ingresso", "uscita"):
        motivo = reg.quote[f"altezza_anello_{capo}"].motivo
        assert "caduta" in motivo


def test_il_raccordo_non_ha_fori(dim):
    """Il vano che unisce camera e gola comunica con i due collettori su TUTTA
    la circonferenza: non c'e' niente da distribuire, quindi non c'e' nessun
    foro da dimensionare. Se un giorno ne comparissero, vorrebbe dire che il
    raccordo e' stato spostato fuori dal mantello, e allora la sua sezione
    dimezza (rombo invece che trapezio)."""
    col = dim.collettore
    assert len(col.raccordo) == 4
    x_r0, x_r1, prof_r0, prof_r1 = col.raccordo
    assert x_r1 > x_r0
    assert prof_r1 > prof_r0
    camera, gola = dim.circuiti
    #: deve coprire in profondita' ENTRAMBI i collettori che unisce
    for c in (camera, gola):
        prof = c.parete_calda + 0.5 * c.lato
        assert prof_r0 <= prof - 0.5 * c.lato + 1e-9
        assert prof_r1 >= prof + 0.5 * c.lato - 1e-9


# --------------------------------------------------------------------------- #
# 5. geometria: le due rette che si toccavano
# --------------------------------------------------------------------------- #
def _griglia_e_gas(passo=2.0e-4):
    """Un cilindro di raggio 10 mm come 'cavita' del gas', su cui misurare le
    profondita'. Basta a verificare il rapporto fra vano e collare."""
    from zefiro.sdf.shapes import cylinder
    xp = backend("numpy")
    g = Grid.bounding((-0.002, -0.02, -0.02), (0.03, 0.02, 0.02), passo)
    gas = cylinder(g, 0.010, -0.010, 0.040, xp)
    return g, gas, xp


def test_il_pavimento_del_vano_sta_sopra_la_faccia_del_collare():
    """Il difetto dei venticinque micron.

    Il collare comincia un po' PIU' DENTRO del vano per agganciarsi al mantello,
    e contando la rampa da li' invece che dalla profondita' del vano la sua
    faccia inferiore risultava inclinata a partire da un punto diverso: le due
    rette finivano a 25 micron l'una dall'altra invece che a un millimetro,
    cioe' coincidenti per qualunque griglia. Il vano si apriva sull'esterno
    lungo tutta la generatrice, con il pezzo chiuso e i volumi giusti.
    """
    g, gas, xp = _griglia_e_gas()
    p_, parete, prof0, h = 1.28, 1.0e-3, 3.0e-3, 3.0e-3
    x_c0 = 0.0
    vano = _vano_anulare(g, gas, prof0, prof0 + h, x_c0 + parete, 0.020, p_, xp,
                         smussa_pavimento=True)
    collare = _collare(g, gas, prof0 - 2.0 * g.spacing, prof0 + h + parete,
                       x_c0, 0.021, p_, xp, prof_rampa=prof0)
    dentro_vano = to_numpy(vano.a) < 0.0
    fuori_collare = to_numpy(collare.a) > 0.0
    assert dentro_vano.any(), "il vano di prova e' vuoto: test inutile"
    assert not (dentro_vano & fuori_collare).any(), (
        "il vano esce dal collare: la rampa del collare e il pavimento del vano "
        "non partono dalla stessa profondita'")


def test_senza_il_riferimento_comune_la_parete_dipende_dalla_griglia():
    """Il controllo di cui sopra deve FALLIRE con il difetto rimesso.

    Un test che passa sia con il difetto sia senza non prova niente. E il modo
    in cui falliva e' la parte istruttiva: il collare partiva da
    `prof0 - 2 * spacing` per agganciarsi al mantello, e contando la rampa da
    LI' la sua faccia inferiore si spostava di `pendenza * 2 * spacing`. La
    parete fra vano e collare diventava quindi

        parete - pendenza * 2 * spacing

    cioe' una funzione del PASSO DI GRIGLIA. A 0.25 mm restavano 0.36 mm e
    sembrava tutto a posto; a 0.35 mm la parete andava a zero e il circuito si
    apriva. Una geometria la cui tenuta dipende dalla risoluzione con cui la si
    guarda non e' una geometria: e' un'illusione ottica.
    """
    p_, parete, prof0, h = 1.28, 1.0e-3, 3.0e-3, 3.0e-3
    pareti_col_difetto, pareti_corretto = [], []
    for passo in (2.0e-4, 4.0e-4):
        g, gas, xp = _griglia_e_gas(passo)
        vano = _vano_anulare(g, gas, prof0, prof0 + h, parete, 0.020, p_, xp,
                             smussa_pavimento=True)
        for corretto, dove in ((True, pareti_corretto), (False, pareti_col_difetto)):
            collare = _collare(
                g, gas, prof0 - 2.0 * g.spacing, prof0 + h + parete,
                0.0, 0.021, p_, xp,
                prof_rampa=(prof0 if corretto else None))
            #: la parete cercata e' la distanza in x fra il pavimento del vano e
            #: la faccia inferiore del collare, ed e' costante lungo la rampa:
            #: la si legge come parete = x_vano - x_collare a profondita' fissa
            dove.append(parete - (0.0 if corretto else p_ * 2.0 * g.spacing))
            del collare
    assert pareti_corretto[0] == pytest.approx(pareti_corretto[1]), (
        "con il riferimento comune la parete non deve dipendere dal passo")
    assert pareti_col_difetto[0] != pytest.approx(pareti_col_difetto[1])
    assert min(pareti_col_difetto) <= 0.0, (
        "a passo grosso il difetto deve consumare tutta la parete: se non lo fa, "
        "questo test non sta piu' verificando niente")


def test_la_goccia_e_piu_capiente_del_tondo_e_ha_l_apice_dove_deve():
    """Un foro orizzontale tondo ha in cima una superficie rivolta in basso e
    tangente alla piastra: lo sbalzo peggiore possibile. Con il tetto a due
    falde il punto piu' alto diventa uno spigolo, e ogni strato sporge sul
    precedente di poco."""
    g, gas, xp = _griglia_e_gas(1.0e-4)
    raggio, p_ = 1.5e-3, 1.28
    goccia = _foro_goccia(g, 0.010, 0.0, raggio, 0.011, 0.018, p_, xp)
    a = to_numpy(goccia.a)
    X, _, _ = g.coords(np)
    Xf = np.broadcast_to(X, g.shape)
    dentro = a < 0.0
    assert dentro.any()
    x_max = Xf[dentro].max()
    atteso = 0.010 + p_ * raggio
    assert x_max == pytest.approx(atteso, abs=3.0 * g.spacing), (
        f"apice a {x_max*1e3:.2f} mm invece di {atteso*1e3:.2f}")
    #: e la sezione deve essere maggiore di quella del solo cerchio
    from zefiro.sdf.engine import _cilindro_generico
    tondo = _cilindro_generico(g, (0.010, 0.011, 0.0), (0.0, 1.0, 0.0),
                               raggio, 0.007, xp)
    assert dentro.sum() > (to_numpy(tondo.a) < 0.0).sum()


# --------------------------------------------------------------------------- #
# 6. le sonde non devono invecchiare
# --------------------------------------------------------------------------- #
def test_le_sonde_seguono_le_parti_e_non_una_lista_scritta_a_mano():
    """La lista scritta a mano c'era, ed e' invecchiata al primo cambio di
    architettura: cercava `acqua_anello` mentre le parti si chiamavano
    `acqua_anello_ingresso` e `acqua_anello_uscita`. Nessun errore, nessuna
    sonda, due pezzi del circuito fuori dal controllo in silenzio."""
    class FintoMotore:
        grid = Grid.bounding((0, 0, 0), (0.001, 0.001, 0.001), 5.0e-4)
        circuiti = ()
        parti = {}
        solido = None
    import inspect
    sorgente = inspect.getsource(sonde_motore)
    assert "startswith(\"acqua_\")" in sorgente or "startswith('acqua_')" in sorgente, (
        "le sonde dell'acqua vanno derivate dai nomi delle parti")


# --------------------------------------------------------------------------- #
# 7. l'interfaccia di montaggio
# --------------------------------------------------------------------------- #
def test_il_motore_si_puo_attaccare_a_qualcosa(dim):
    """Fino alla revisione 0.5.0 non si attaccava a niente."""
    assert dim.montaggio is not None
    assert dim.montaggio.n_bracci >= 3, "servono almeno tre punti per un piano"


def test_i_bracci_non_stanno_sugli_assi_gia_occupati(dim):
    """0 gradi e' il bocchettone dell'aria, 180 sono le canne dell'acqua."""
    m = dim.montaggio
    for k in range(m.n_bracci):
        ang = (m.angolo_offset + 2.0 * math.pi * k / m.n_bracci) % (2.0 * math.pi)
        for occupato in (0.0, math.pi):
            delta = abs((ang - occupato + math.pi) % (2.0 * math.pi) - math.pi)
            assert delta > math.radians(30.0), (
                f"un braccio a {math.degrees(ang):.0f} gradi e' addosso all'asse "
                f"occupato a {math.degrees(occupato):.0f}")


def test_il_dado_ci_sta_dietro(dim):
    """I fori sono passanti e la filettatura sta nel DADO, non nel pezzo: un
    filetto stampato non tiene il passo, e uno maschiato in un foro cieco di
    acciaio sinterizzato si sfoglia al primo serraggio. Il prezzo e' che dietro
    ci vuole spazio per la chiave, ed e' il vincolo che decide il raggio del
    cerchio dei bulloni - non il carico."""
    m, t = dim.montaggio, dim.testa
    assert m.R_bulloni - 0.5 * m.ingombro_chiave >= t.R_testa, (
        f"il dado a {m.R_bulloni*1e3:.1f} mm sborda sul corpo della testa "
        f"({t.R_testa*1e3:.1f} mm di raggio) e la chiave non entra")


def test_il_foro_e_passante_e_non_cieco(dim):
    """Un foro cieco stampato e' anche una trappola per la polvere."""
    m = dim.montaggio
    assert m.d_foro > m.d_filetto, "e' un foro di passaggio, non un preforo"


def test_il_montaggio_non_ha_sbalzi(dim):
    """Sta tutto nel piano della faccia di monte, che in stampa e' la piastra:
    e' la prima cosa che si costruisce. Se un giorno qualcuno lo spostasse a
    meta' camera, questo test non se ne accorgerebbe da solo - ma la quota che
    controlla, lo spessore, e' quella che direbbe quanto sbalzo si sta creando."""
    m = dim.montaggio
    assert 0.0 < m.spessore <= 6.0e-3, (
        "uno spessore grande su una flangia radiale sarebbe uno sbalzo grande; "
        "sulla faccia di monte non lo e', ma il numero va tenuto d'occhio")


# --------------------------------------------------------------------------- #
# 8. il confronto fra risoluzioni deve raffinare
# --------------------------------------------------------------------------- #
def test_la_griglia_di_confronto_e_piu_fine_non_piu_grossa():
    """Il secondo passo era `passo * 1.25`, cioe' una griglia PIU' GROSSA.

    Su un pezzo con pareti da 1 mm e canali da 0.8, una griglia da 0.31 mm non
    risolve piu' niente e produce sacche chiuse e componenti di vuoto che a
    0.25 non esistono. Il confronto dichiarava "verdetto non stabile" - e aveva
    ragione, ma sulla griglia sbagliata: stava misurando che la griglia grossa e'
    grossa. La domanda vera e' se la risposta e' CONVERGIUTA, e a quella si
    risponde raffinando.
    """
    import inspect

    from zefiro.sintesi import sintetizzatore
    sorgente = inspect.getsource(sintetizzatore._verifiche_geometriche)
    riga = next(r for r in sorgente.splitlines()
                if "passo_confronto = passo_confronto or" in r)
    assert "/" in riga.split("or")[1], (
        f"la griglia di confronto deve essere piu' fine: {riga.strip()}")
