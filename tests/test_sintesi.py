"""Il sintetizzatore: progetta(requisiti) -> Progetto | Rifiuto.

I test qui rispondono a tre domande diverse, e la seconda e' la piu' importante:

  1. il modello ritrova cio' che era stato ricavato a mano? (regressione)
  2. il modello RIFIUTA quando deve? (un modello che consegna sempre qualcosa
     inventa i dati che non ha)
  3. e' un modello o e' un progetto travestito? (cambiando i requisiti le quote
     devono muoversi, e muoversi nel verso giusto)
"""
from __future__ import annotations

import math

import pytest

from zefiro.sintesi import Impianto, Processo, Requisiti, Rifiuto, progetta
from zefiro.sintesi.architetture import aerospike_gas_gas as ARCH
from zefiro.sintesi.esito import Esito
from zefiro.sintesi.provenienza import Origine, Registro, Scelta


def requisiti_zefiro(**k) -> Requisiti:
    """I requisiti di Zefiro cosi' come sono stati decisi a mano."""
    campi = dict(spinta=50.0, durata=5.0)
    campi.update({a: b for a, b in k.items() if a not in ("impianto",)})
    imp = dict(combustibile={"C3H8": 1.0}, T_bombola_min=288.15,
               T_ossidante_iniezione=293.0, T_combustibile_iniezione=283.0,
               cd_ossidante=0.75, cd_combustibile=0.75)
    imp.update(k.get("impianto", {}))
    return Requisiti(impianto=Impianto(**imp), **campi)


@pytest.fixture(scope="module")
def zefiro():
    return progetta(requisiti_zefiro(), geometria=False)


# --------------------------------------------------------------------------- #
# 1. regressione: il modello deve ritrovare il progetto fatto a mano
# --------------------------------------------------------------------------- #
def test_ritrova_la_pressione_di_camera_ricavata_a_mano(zefiro):
    """p_c viene solo dalla bombola: p_sat del propano a 15 C diviso 1.15.
    Non dipende da nessuna scelta di progetto a valle, quindi e' il pezzo del
    calcolo a mano che DEVE restare identico qualunque cosa cambi dopo."""
    assert zefiro.inviluppo.p_camera == pytest.approx(6.361e5, rel=1e-3)
    assert zefiro.inviluppo.spinta == pytest.approx(50.0, rel=1e-9)


def test_il_punto_di_progetto_corrente(zefiro):
    """Congelamento del punto attuale.

    NON coincide piu' con quello ricavato a mano, ed e' un progresso, non una
    deriva: il progetto a mano dava Isp 143.7 s ed eps 1.622 con un corpo
    centrale che FONDEVA a 2.8 s. Deviando il 9 % del combustibile al film che
    lo protegge, quel combustibile brucia comunque - piu' a valle e peggio
    miscelato, ma brucia - e l'Isp SALE a 148.2 s.

    Il vecchio numero era migliore solo perche' descriveva un motore che non
    arrivava a fine raffica.
    """
    i = zefiro.inviluppo
    assert i.isp == pytest.approx(148.2, abs=0.3)
    assert i.epsilon == pytest.approx(1.636, abs=0.003)
    assert i.durata_raffica == pytest.approx(7.27, abs=0.06)
    assert zefiro.registro.valore("frazione_film_plug") == pytest.approx(0.09, abs=0.006)


def test_l_iniettore_segue_il_film(zefiro):
    """L'iniettore e' cambiato da solo quando il film gli ha tolto il 9 % del
    combustibile: meno portata, meno area di efflusso, e la regola di scelta
    del raggio si sposta su 4 getti piu' grossi invece di 5.

    E' il comportamento voluto - il raggio d'iniezione non e' una costante, e'
    la soluzione di un problema che dipende dalle portate - ma va verificato,
    perche' un iniettore che NON si muove quando cambiano le portate sarebbe un
    iniettore scritto a mano.
    """
    inj = zefiro.geometria.iniettore
    assert inj.n_getti == 4
    assert inj.d_getto == pytest.approx(0.731e-3, abs=6e-6)
    assert inj.J == pytest.approx(0.99, abs=0.02)
    assert abs(inj.scarto_C) < 0.02
    # e resta sopra il minimo stampabile, con piu' margine di prima
    assert inj.d_getto > zefiro.requisiti.processo.margine_foro()


def test_sceglie_le_stesse_filettature(zefiro):
    """G3/8 per l'aria e G1/8 per il GPL non erano una scelta di catalogo: era
    un bilancio di pressione. Il modello deve rifarlo e arrivarci."""
    from zefiro.injector import BSPP

    reg = zefiro.registro
    assert reg.valore("filettatura_aria") == pytest.approx(BSPP["G3/8"]["d_esterno"])
    assert reg.valore("filettatura_GPL") == pytest.approx(BSPP["G1/8"]["d_esterno"])


#: Difetti NOTI e non ancora risolti. Un test che li elenca serve a due cose:
#: che nessuno se ne dimentichi, e che nessun ALTRO difetto passi inosservato
#: nascosto in mezzo a quelli gia' noti.
DIFETTI_APERTI: set[str] = set()


def test_le_uniche_verifiche_che_falliscono_sono_quelle_note(zefiro):
    """L'elenco e' vuoto, e ci e' voluto un film per svuotarlo.

    Il corpo centrale senza film fonde a 2.8 s contro i 5 richiesti. Non e'
    raffreddabile in nessun altro modo: un gas dall'interno da' h di 3 kW/m2K
    contro i 29 dell'acqua, l'acqua non ci sta nel corpo centrale strizzato,
    troncare cambia 0.2 s e ingrandire il motore non basta. Resta il film.
    """
    fallite = {v.nome for v in zefiro.verifiche if v.esito is Esito.FALLITA}
    assert fallite == DIFETTI_APERTI, (
        f"inattese: {fallite - DIFETTI_APERTI}; "
        f"risolte senza aggiornare l'elenco: {DIFETTI_APERTI - fallite}")
    non_concl = [v.nome for v in zefiro.verifiche if v.esito is Esito.NON_CONCLUSIVA]
    assert not non_concl, f"verifiche non conclusive: {non_concl}"


def test_il_film_del_plug_e_derivato_e_non_scelto(zefiro):
    """La frazione di film e' la MINIMA che tiene il plug sotto soglia anche se
    il film rende meta' di quanto la correlazione promette. E' un punto fisso:
    il film sposta phi del nucleo, quindi T_ad, quindi il carico sul plug,
    quindi il film che serve."""
    import zefiro.sintesi.architetture.aerospike_gas_gas as A

    f = zefiro.registro.valore("frazione_film_plug")
    assert 0.0 < f < 0.35
    T = zefiro.registro.valore("T_plug_fine_raffica")
    assert T <= A.T_LIMITE_PLUG
    # e senza film il plug NON ce la farebbe: il vincolo e' attivo
    st, _ = A.stato_plug(zefiro.requisiti, zefiro.geometria.op,
                         zefiro.geometria.params, zefiro.geometria.l0,
                         zefiro.inviluppo.p_camera, 0.0)
    assert st.T_finale > A.T_LIMITE_PLUG


def test_una_raffica_troppo_lunga_non_e_salvabile_col_film():
    """C'e' un limite: il film costa combustibile, e il combustibile e' la
    risorsa piu' scarsa del motore. Oltre una certa durata il modello deve
    rifiutare invece di deviare meta' del GPL."""
    r = progetta(requisiti_zefiro(durata=60.0), geometria=False)
    assert isinstance(r, Rifiuto)


# --------------------------------------------------------------------------- #
# 2. il rifiuto: la parte del contratto che rende il modello onesto
# --------------------------------------------------------------------------- #
def test_rifiuta_se_mancano_le_temperature_a_valle_dei_riduttori():
    """L'espansione attraverso un riduttore e' quasi isentalpica e raffredda il
    gas: la temperatura d'iniezione NON e' quella ambiente e va misurata.
    Inventarla sposta la densita' e quindi tutte le aree di efflusso."""
    r = progetta(requisiti_zefiro(impianto={"T_ossidante_iniezione": None}),
                 geometria=False)
    assert isinstance(r, Rifiuto)
    assert any("T_ossidante" in d for d in r.dati_mancanti)


def test_rifiuta_se_mancano_i_coefficienti_di_efflusso():
    r = progetta(requisiti_zefiro(impianto={"cd_ossidante": None}), geometria=False)
    assert isinstance(r, Rifiuto)
    assert any("cd_" in d for d in r.dati_mancanti)


def test_elenca_TUTTI_i_dati_mancanti_insieme():
    """Uno alla volta significa far ripartire chi progetta in laboratorio
    cinque volte."""
    r = progetta(requisiti_zefiro(impianto={"cd_ossidante": None,
                                            "T_ossidante_iniezione": None,
                                            "T_combustibile_iniezione": None}),
                 geometria=False)
    assert isinstance(r, Rifiuto)
    assert len(r.dati_mancanti) >= 2


def test_rifiuta_una_spinta_fuori_dall_inviluppo_verificato():
    r = progetta(requisiti_zefiro(spinta=50000.0), geometria=False)
    assert isinstance(r, Rifiuto)
    assert r.architetture_provate
    assert "spinta" in r.architetture_provate[0][1]


def test_rifiuta_un_ossidante_che_non_sa_trattare():
    with pytest.raises(NotImplementedError, match="gas-gas"):
        requisiti_zefiro(impianto={"ossidante": "N2O"})


def test_rifiuta_quando_il_raffreddamento_non_ce_la_fa():
    """Con una parete minima di processo troppo grossa il salto termico
    attraverso la parete calda supera il limite e non esiste parete possibile.
    Il modello deve dirlo, non stampare una parete che fondera'."""
    r = progetta(requisiti_zefiro(), geometria=False)
    assert not isinstance(r, Rifiuto)
    grosso = Processo(parete_minima=3.0e-3, foro_minimo=3.0e-3)
    r2 = progetta(Requisiti(spinta=50.0, durata=5.0,
                            impianto=requisiti_zefiro().impianto, processo=grosso),
                  geometria=False)
    assert isinstance(r2, Rifiuto), "una parete da 3 mm in gola non puo' funzionare"


def test_un_rifiuto_non_e_consegnabile():
    r = progetta(requisiti_zefiro(impianto={"cd_ossidante": None}), geometria=False)
    assert r.consegnabile is False
    assert "RIFIUTO" in r.referto()


# --------------------------------------------------------------------------- #
# 3. e' un modello, non un progetto travestito
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("spinta", [50.0, 100.0, 200.0])
def test_la_spinta_richiesta_esce_esatta(spinta):
    p = progetta(requisiti_zefiro(spinta=spinta), geometria=False)
    assert not isinstance(p, Rifiuto), getattr(p, "architetture_provate", "")
    assert p.inviluppo.spinta == pytest.approx(spinta, rel=1e-9)


def test_il_plug_non_raffreddato_e_un_limite_di_tutta_la_famiglia():
    """Non e' un caso particolare di Zefiro: il plug scala male. La massa va
    come il cubo della scala e il calore come il quadrato, quindi un motore
    PIU' GRANDE dovrebbe reggere meglio - e infatti il tempo alla fusione
    cresce con la spinta. Serve sapere di quanto, perche' dice se la cura e'
    ingrandire il motore o cambiare architettura."""
    from zefiro.sintesi.sintetizzatore import _durata_del_plug

    tempi = []
    for spinta in (50.0, 200.0):
        p = progetta(requisiti_zefiro(spinta=spinta), geometria=False)
        tempi.append(_durata_del_plug(p.geometria)[1])
    assert tempi[1] > tempi[0], "il plug deve durare di piu' su un motore piu' grande"


def test_sotto_una_certa_spinta_i_getti_non_sono_piu_stampabili():
    """Un limite VERO dell'architettura, e il modello lo dice con un numero.

    Il diametro dei getti va come la radice della portata di combustibile, che
    va come la spinta: sotto una certa spinta i fori scendono sotto il minimo
    stampabile e non c'e' raggio d'iniezione che li salvi. A 25 N il massimo
    ottenibile e' 0.61 mm contro 0.66 richiesti.

    Il rifiuto non si limita a dire di no: dice quale spinta funzionerebbe e
    quale sarebbe l'alternativa (fori realizzati dopo la stampa).
    """
    r = progetta(requisiti_zefiro(spinta=25.0), geometria=False)
    assert isinstance(r, Rifiuto)
    motivo = r.architetture_provate[0][1]
    assert "N" in motivo and "elettroerosione" in motivo


def test_sopra_una_certa_spinta_non_esiste_raccordo_che_porti_l_aria():
    """L'altro capo dell'inviluppo: a 500 N servono 336 g/s d'aria, cioe' 39
    litri al secondo, e nemmeno un G1 li porta senza mangiarsi il salto
    d'iniezione. Non e' il motore a non stare in piedi: e' l'impianto."""
    r = progetta(requisiti_zefiro(spinta=500.0, durata=0.5), geometria=False)
    assert isinstance(r, Rifiuto)
    assert "filettatura" in r.architetture_provate[0][1]


def test_una_bombola_piu_calda_alza_p_c_e_accorcia_la_raffica():
    """Il risultato controintuitivo che dimensiona il motore, e che il modello
    deve riprodurre da solo: la bombola va tenuta FRESCA."""
    freddo = progetta(requisiti_zefiro(impianto={"T_bombola_min": 283.15}),
                      geometria=False)
    caldo = progetta(requisiti_zefiro(impianto={"T_bombola_min": 293.15}),
                     geometria=False)
    assert caldo.inviluppo.p_camera > freddo.inviluppo.p_camera
    assert caldo.inviluppo.isp > freddo.inviluppo.isp
    assert caldo.inviluppo.durata_raffica < freddo.inviluppo.durata_raffica


def test_una_spinta_maggiore_chiede_piu_raffreddamento():
    piccolo = progetta(requisiti_zefiro(spinta=50.0), geometria=False)
    grande = progetta(requisiti_zefiro(spinta=200.0), geometria=False)
    assert (grande.registro.valore("potenza_termica")
            > piccolo.registro.valore("potenza_termica"))
    assert (grande.registro.valore("portata_acqua_totale")
            > piccolo.registro.valore("portata_acqua_totale"))


def test_una_macchina_piu_fine_fa_fori_piu_piccoli_e_piu_numerosi():
    """Il diametro dei getti e' vincolato dal minimo stampabile: se la macchina
    ne fa di piu' piccoli, il modello usa piu' getti, che e' meglio perche' un
    foro otturato fa meno danno."""
    base = progetta(requisiti_zefiro(), geometria=False)
    fine = progetta(Requisiti(
        spinta=50.0, durata=5.0, impianto=requisiti_zefiro().impianto,
        processo=Processo(foro_minimo=3.0e-4, misurati=frozenset({"foro_minimo"}))),
        geometria=False)
    assert fine.geometria.iniettore.n_getti > base.geometria.iniettore.n_getti
    assert fine.geometria.iniettore.d_getto < base.geometria.iniettore.d_getto


def test_il_raccordo_cresce_con_la_portata():
    """L'aria e' 3.9 l/s a 50 N: un G1/4 si mangerebbe un terzo del salto
    d'iniezione. A spinta maggiore serve una filettatura maggiore, e con essa
    una testa piu' lunga: la ferramenta risale la catena."""
    piccolo = progetta(requisiti_zefiro(spinta=50.0), geometria=False)
    grande = progetta(requisiti_zefiro(spinta=200.0), geometria=False)
    assert (grande.registro.valore("filettatura_aria")
            > piccolo.registro.valore("filettatura_aria"))
    assert (grande.registro.valore("L_plenum_aria")
            >= piccolo.registro.valore("filettatura_aria"))


# --------------------------------------------------------------------------- #
# La provenienza: il meccanismo che rende leggibile tutto il resto
# --------------------------------------------------------------------------- #
def test_ogni_quota_ha_un_motivo():
    with pytest.raises(ValueError, match="non ha motivo"):
        Scelta("x", 1.0, "m", Origine.DECISA, "   ")


def test_una_quota_mancante_non_puo_avere_un_valore_finito():
    """Un segnaposto numerico si propaga nei conti e produce un risultato
    credibile e falso. E' il modo piu' comune in cui muore un modello."""
    with pytest.raises(ValueError, match="MANCANTE"):
        Scelta("d_min", 6.0e-4, "m", Origine.MANCANTE, "da confermare col service")
    Scelta("d_min", float("nan"), "m", Origine.MANCANTE, "da confermare col service")


def test_due_regole_in_conflitto_sullo_stesso_nome_sono_un_errore():
    r = Registro()
    r.aggiungi(Scelta("a", 1.0, "m", Origine.DERIVATA, "prima regola"))
    with pytest.raises(ValueError, match="conflitto"):
        r.aggiungi(Scelta("a", 2.0, "m", Origine.DERIVATA, "seconda regola"))


def test_la_maggior_parte_delle_quote_e_derivata(zefiro):
    """La misura di maturita' del modello. Non c'e' un valore giusto, ma c'e'
    un verso: le quote DECISE sono quelle che qualcuno un giorno rimettera' in
    discussione, e devono essere poche e ben motivate."""
    m = zefiro.registro.maturita()
    assert m["mancante"] == 0.0
    assert m["derivata"] > 0.7, zefiro.registro.referto()


def test_le_quote_decise_sono_elencabili_una_per_una(zefiro):
    decise = zefiro.registro.per_origine()[Origine.DECISA]
    assert 1 <= len(decise) <= 12
    for s in decise:
        assert len(s.motivo) > 40, f"{s.nome}: motivo troppo scarno"
