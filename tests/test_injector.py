"""Iniettore a getti trasversali: efflusso comprimibile e regola di Holdeman."""
from __future__ import annotations

import math

import pytest

from zefiro.injector import (
    C_OTTIMO, DISPOSIZIONI_A_DUE_FILE, R_UNIVERSALE, disposizione_realizzabile,
    efflusso_gas, getti_massimi, progetta, spazio_di_scelta,
)

ARIA = dict(massa_molare=28.9649, gamma=1.400)
GPL = dict(massa_molare=44.0956, gamma=1.130)


def test_l_orifizio_si_strozza_al_rapporto_critico():
    """Sotto il rapporto critico la portata non dipende piu' dalla pressione a
    valle. Il numero di controllo e' analitico: il flusso di massa strozzato
    vale p0 sqrt(gamma/(R T0)) (2/(g+1))^((g+1)/(2(g-1)))."""
    p0, T0 = 10.0e5, 300.0
    e = efflusso_gas(1.0e-3, p0, 1.0e5, T0, **ARIA, cd=1.0)
    assert e.strozzato
    assert e.mach == pytest.approx(1.0, rel=1e-6)
    g = ARIA["gamma"]
    R = R_UNIVERSALE / ARIA["massa_molare"]
    G = p0 * math.sqrt(g / (R * T0)) * (2.0 / (g + 1.0)) ** ((g + 1.0) / (2.0 * (g - 1.0)))
    assert 1.0e-3 / e.area == pytest.approx(G, rel=1e-6)


def test_sotto_il_rapporto_critico_l_efflusso_e_subsonico():
    e = efflusso_gas(1.0e-3, 7.315e5, 6.361e5, 293.0, **ARIA, cd=0.75)
    assert not e.strozzato
    assert 0.0 < e.mach < 1.0


def test_il_modello_incomprimibile_sbaglia_l_area_di_alcuni_punti():
    """La ragione per cui questo modulo non usa mdot = Cd A sqrt(2 rho Dp).

    Non e' un dettaglio accademico: la differenza sull'area si trasferisce
    tutta sul diametro dei fori, cioe' sul numero che decide se il pezzo e'
    stampabile."""
    p0, p1, T0, cd = 7.315e5, 6.361e5, 293.0, 0.75
    e = efflusso_gas(2.0e-3, p0, p1, T0, **ARIA, cd=cd)
    R = R_UNIVERSALE / ARIA["massa_molare"]
    rho0 = p0 / (R * T0)
    a_incomp = 2.0e-3 / (cd * math.sqrt(2.0 * rho0 * (p0 - p1)))
    scarto = abs(a_incomp - e.area) / e.area
    assert 0.02 < scarto < 0.20, scarto


def test_la_portata_torna_dalla_geometria_ricavata():
    """Controllo di chiusura: rifacendo il conto al contrario dai fori si deve
    ritrovare la portata di partenza."""
    i = _iniettore()
    a_tot = i.n_getti * math.pi / 4.0 * i.d_getto ** 2
    mdot = 0.75 * i.rho_getto * i.V_getto * a_tot
    assert mdot == pytest.approx(1.926e-3, rel=1e-9)


def _iniettore(R=3.10e-3, n=None, disp="un_lato"):
    return progetta(
        mdot_aria=33.56e-3, mdot_gpl=1.926e-3, p_c=6.361e5,
        p_aria_monte=7.315e5, p_gpl_monte=7.315e5, T_aria=293.0, T_gpl=283.0,
        MW_aria=ARIA["massa_molare"], MW_gpl=GPL["massa_molare"],
        gamma_aria=ARIA["gamma"], gamma_gpl=GPL["gamma"], cd=0.75,
        R_iniezione=R, n_getti=n, disposizione=disp)


def test_C_e_coerente_con_S_H_e_J():
    i = _iniettore()
    assert i.C == pytest.approx((i.passo / i.altezza_anello) * math.sqrt(i.J))


def test_salti_di_pressione_uguali_danno_J_vicino_a_uno():
    """Non e' una coincidenza ed e' bene che sia scritto: con la stessa
    frazione di p_c come salto e densita' confrontabili, il rapporto dei flussi
    di quantita' di moto viene vicino a 1 da solo. J NON e' un parametro che si
    sceglie in questo motore: e' una conseguenza della bombola."""
    assert _iniettore().J == pytest.approx(1.0, abs=0.05)


def test_l_area_dell_anello_e_quella_esatta_non_quella_di_corona_sottile():
    """A = pi (Re^2 - Ri^2), non 2 pi R H.

    Con H/R = 0.5 l'approssimazione di corona sottile sbaglia l'area del 25 %,
    e quell'errore si propaga al quadrato: velocita' dell'aria +25 %, J -56 %,
    passo dei getti -25 %. E' un errore che c'era, ed e' il tipo di errore che
    non si vede perche' tutti i numeri restano plausibili.
    """
    i = _iniettore()
    Ri, Re = i.R_iniezione, i.R_medio
    assert math.pi * (Re ** 2 - Ri ** 2) == pytest.approx(i.area_aria, rel=1e-12)
    sottile = 2.0 * math.pi * Ri * i.altezza_anello
    assert abs(sottile / i.area_aria - 1.0) > 0.15


def test_il_diametro_dei_getti_e_proporzionale_all_altezza_del_condotto():
    """L'invariante che governa tutto il dimensionamento.

    A portate e salti fissati vale  d = k H  con k che non dipende dal raggio:
    combinando  n = 2 pi R sqrt(J) / (C H),  A_getti = n pi d^2/4  e
    H = A_aria / (2 pi R), il raggio si semplifica. La conseguenza pratica e'
    controintuitiva e va verificata, non creduta: **strizzare il corpo
    centrale fa fori piu' grandi**.
    """
    coppie = [(r, _iniettore(R=r)) for r in (2.8e-3, 3.2e-3, 3.6e-3, 4.0e-3)]
    k = [i.d_getto / i.altezza_anello for _, i in coppie]
    # l'arrotondamento di n a intero introduce qualche punto percentuale
    assert max(k) / min(k) < 1.20
    # e il verso e' quello: raggio piu' piccolo, foro piu' grande
    d = [i.d_getto for _, i in coppie]
    assert d == sorted(d, reverse=True)


def test_l_area_dell_aria_non_dipende_dal_raggio_di_iniezione():
    """E' fissata dalla portata e dal salto: il raggio ne cambia solo la forma
    (piu' stretta e alta, o piu' larga e bassa), non l'area."""
    a = [_iniettore(R=r).area_aria for r in (2.8e-3, 4.0e-3)]
    assert a[0] == pytest.approx(a[1], rel=1e-12)


def test_la_fila_singola_punta_al_proprio_C():
    """E l'errore residuo e' SOLO quello dell'arrotondamento di n a intero.

    Il limite si scrive: se n_ideale = 2 pi R sqrt(J)/(C* H) e n e' il suo
    intero piu' vicino, allora |C/C* - 1| <= 1/(2n). Verificarlo cosi' invece
    che con una tolleranza a occhio serve a distinguere un errore di
    arrotondamento (inevitabile) da un errore di formula (grave)."""
    i = _iniettore(disp="un_lato")
    assert i.C_obiettivo == C_OTTIMO["un_lato"]
    assert abs(i.C / i.C_obiettivo - 1.0) <= 1.0 / (2.0 * i.n_getti) + 1e-9


@pytest.mark.parametrize("disp", sorted(DISPOSIZIONI_A_DUE_FILE))
def test_una_disposizione_a_due_file_si_rifiuta_invece_di_fingere(disp):
    """IL DIFETTO CHE QUESTO TEST BLOCCA, ed e' il piu' insidioso di tutto il
    modulo perche' non falliva: produceva numeri.

    Chiedendo `opposti_allineati` la funzione restituiva UNA fila sola di 8
    getti con C = 1.25, etichettata con l'ottimo della disposizione opposta.
    Ma 1.25 per una fila SOLA e' esattamente il caso che Holdeman giudica
    pessimo - getti troppo fitti che restano attaccati alla parete. Il modello
    diceva "ottimo" di una geometria che, cosi' com'era costruita, e' la
    peggiore. Nessuna avvertenza, nessun errore: solo un motore sbagliato.
    """
    with pytest.raises(ValueError, match="due file|opposta"):
        _iniettore(disp=disp)


def test_su_un_anello_le_file_allineate_non_stanno_entrambe_sull_ottimo():
    """Un fatto GEOMETRICO, non una limitazione del codice.

    C = (S/H) sqrt(J): J e H sono comuni alle due file, quindi C sta al passo.
    "Allineate" vuol dire stessi azimut, quindi stesso numero di getti su due
    circonferenze diverse, quindi S_e/S_i = R_e/R_i e C_e/C_i pure. Su questo
    anello il rapporto vale 1.6: una delle due file e' fuori del 60 %.
    Holdeman studia un condotto rettangolare, dove le due pareti sono lunghe
    uguale e il problema non esiste - per questo la conclusione non e' nella
    correlazione e va ricavata qui."""
    R_i, H = 2.646e-3, 1.654e-3
    ok, motivo = disposizione_realizzabile("opposti_allineati", R_i, R_i + H)
    assert not ok and "SFALSATA" in motivo
    #: sfalsate invece non chiedono gli stessi azimut: ogni fila sceglie il
    #: proprio numero di getti e sta sul proprio ottimo.
    assert disposizione_realizzabile("opposti_sfalsati", R_i, R_i + H)[0]
    #: e su un condotto quasi piano (R_e ~ R_i) il problema sparisce, che e'
    #: il caso di Holdeman.
    assert disposizione_realizzabile("opposti_allineati", 1.0, 1.05)[0]


def test_il_numero_di_getti_ha_un_tetto_che_lo_fissa_il_foro_minimo():
    """Il vincolo che decide tutta l'architettura d'iniezione, e che non si
    aggira scegliendo una disposizione.

    L'area totale dei getti e' portata / (densita' * velocita'), e la velocita'
    la fissa il salto di pressione disponibile: e' un dato, non una scelta.
    Quell'area diviso l'area del foro minimo realizzabile da' quanti getti si
    possono avere. Punto. Se il tetto e' 5, nessuna disposizione che ne chiede
    21 e' raggiungibile, per quanto sia migliore."""
    i = _iniettore()
    #: si verifica la FORMULA, non un numero: un numero dipenderebbe dal punto
    #: di prova e cambierebbe a ogni ritocco della fixture, che e' il modo piu'
    #: sicuro di avere un test che non dice piu' niente.
    A_tot = math.pi / 4.0 * i.d_getto ** 2 * i.n_getti
    for d_min in (0.60e-3, 0.45e-3, 0.30e-3):
        assert getti_massimi(i, d_min) == int(A_tot / (math.pi / 4.0 * d_min ** 2))
    #: il tetto va come 1/d_min^2: dimezzare la soglia quadruplica i getti
    #: possibili. E' il motivo per cui `spazio_di_scelta` prende d_min come
    #: argomento invece di leggerlo - il numero non e' mai stato misurato e il
    #: risultato cambia QUALITATIVAMENTE nell'intervallo plausibile.
    assert getti_massimi(i, 0.30e-3) >= 3 * getti_massimi(i, 0.60e-3)


def test_lo_spazio_di_scelta_boccia_cio_che_non_e_raggiungibile():
    i = _iniettore()
    alt = {a.disposizione: a for a in spazio_di_scelta(i)}
    assert set(alt) == set(C_OTTIMO)
    #: la disposizione migliore in un condotto rettangolare (C* = 1.25, ogni
    #: getto attraversa mezza altezza) qui e' fuori due volte: geometria e foro.
    a = alt["opposti_allineati"]
    assert not a.realizzabile_geometricamente
    assert not a.stampabile(0.60e-3)
    assert a.n_totale > getti_massimi(i, 0.60e-3)
    #: quella scelta e' l'unica che sta dentro tutti e due i vincoli con
    #: l'ipotesi di fila periodica ancora in piedi.
    assert alt["un_lato"].realizzabile_geometricamente
    assert alt["un_lato"].stampabile(0.60e-3)


def test_imporre_pochi_getti_fa_sovra_penetrare():
    """Meno getti dello stretto necessario = passo troppo largo = C troppo
    grande = getti che attraversano il condotto e sbattono dall'altra parte."""
    pochi = _iniettore(n=2)
    assert pochi.C > 2.0 * pochi.C_obiettivo
    assert pochi.d_getto > _iniettore().d_getto


def test_il_punto_scelto_per_il_motore_da_50N():
    """Congelamento del punto di progetto. Se qualcuno cambia una pressione o
    una portata a monte, questo test si accorge che l'iniettore non e' piu'
    quello disegnato."""
    i = _iniettore()
    assert i.n_getti == 5
    assert i.d_getto == pytest.approx(0.666e-3, abs=5e-6)
    assert i.altezza_anello == pytest.approx(1.540e-3, abs=5e-6)
    assert i.R_medio == pytest.approx(4.640e-3, abs=5e-6)
    assert abs(i.scarto_C) < 0.01
    assert i.stampabile


def test_p_valle_maggiore_di_p_monte_e_un_errore():
    with pytest.raises(ValueError, match="nessun efflusso"):
        efflusso_gas(1e-3, 1.0e5, 2.0e5, 300.0, **ARIA, cd=1.0)


def test_disposizione_sconosciuta_e_un_errore():
    with pytest.raises(ValueError):
        _iniettore(disp="a_caso")
