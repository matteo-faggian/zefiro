"""Il caso OpenFOAM del mescolamento: le scelte che non si possono perdere.

Ognuno di questi test corrisponde a un modo in cui il caso e' gia' andato
storto, e in cui andrebbe storto di nuovo se qualcuno cambiasse una riga senza
sapere perche' c'era. Nessuno di questi difetti fa fallire il solutore in modo
evidente: due su tre lo fanno girare per ore producendo numeri.
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts" / "cfd"))


def quote_finte() -> dict:
    """Un dizionario di quote plausibile, per scrivere il caso senza dover
    dimensionare tutto il motore."""
    return {
        "n_getti": 4, "R_getti": 2.65e-3, "R_anello": 4.30e-3,
        "d_getto": 7.31e-4, "x_getti": -3.31e-3, "r_bore": 1.65e-3,
        "x_ingresso": -6.62e-3, "R_c": 12.08e-3, "r_centerbody": 3.91e-3,
        "x_rampa": 5.95e-3,
        "V_aria": 151.8, "rho_aria": 7.56, "T_aria": 293.0,
        "V_gpl": 121.6, "rho_gpl": 11.9, "T_gpl": 283.0,
        "p_c": 6.361e5, "mdot_aria": 0.0324, "mdot_gpl": 0.00186,
        "H": 1.654e-3, "C": 2.498, "J": 0.988, "L_mescolamento": 3.31e-3,
        "piani_x": [-3.31e-3 + f * 3.31e-3 for f in (0.25, 0.5, 1.0, 1.5)]
                   + [0.0, 3.0e-3, 8.0e-3],
    }


@pytest.fixture(scope="module")
def caso(tmp_path_factory):
    pytest.importorskip("cantera")
    from caso_mescolamento import scrivi
    d = tmp_path_factory.mktemp("caso") / "prova"
    scrivi(d, quote_finte(), 1.0e-4, 1.0e-5, 4, campionamento=1.0e-6, binario=True)
    return d


def leggi(caso, *parti) -> str:
    return (caso.joinpath(*parti)).read_text(encoding="utf-8")


def voce(testo: str, chiave: str) -> str:
    """Il valore di una voce di dizionario OpenFOAM, cercata a INIZIO RIGA.

    Cercarla come sottostringa qualunque non funziona: `stopAt endTime;` viene
    prima di `endTime 0.00013;` e la prima ricerca trova la parola dentro il
    valore di un'altra voce."""
    for riga in testo.splitlines():
        pezzi = riga.strip().rstrip(";").split(None, 1)
        if len(pezzi) == 2 and pezzi[0] == chiave:
            return pezzi[1].strip()
    raise KeyError(chiave)


# --------------------------------------------------------------------------- #
# 1. lo sbocco non puo' essere uno specchio acustico
# --------------------------------------------------------------------------- #
def test_lo_sbocco_non_e_a_pressione_imposta(caso):
    """IL DIFETTO PIU' COSTOSO DI TUTTA LA CATENA, e il piu' invisibile.

    Con `p: fixedValue` allo sbocco il contorno riflette perfettamente. Il
    dominio parte fermo e gli ingressi cominciano a soffiare a 152 m/s: l'onda
    di compressione che ne nasce percorre i 18.6 mm, arriva allo sbocco e torna
    indietro tutta intera. Il calcolo e' esploso DUE VOLTE allo stesso istante -
    56.6 e 56.2 microsecondi - e un transito acustico del dominio vale
    18.6 mm / 343 m/s = 54 microsecondi. Non era il getto: era il contorno.

    Il sintomo non e' un crollo: e' il passo temporale che va a 2.8 picosecondi
    con il Courant MEDIO a 3e-6 e il MASSIMO a 0.4. Il solutore continua a
    girare, non stampa nessun errore, e in nove ore non arriva a meta'.

    QUESTO TEST E' STATO RISCRITTO DOPO LA CORSA 5, e vale la pena dire
    perche'. Diceva `assert "waveTransmissive" in blocco`: cioe' fissava la
    SOLUZIONE di allora invece del REQUISITO. La soluzione di allora ha poi
    causato il guasto successivo - uno sbocco trasparente non vincola la
    pressione media, e la corsa 5 e' andata alla deriva fino a strozzare il
    getto. Il requisito vero non e' "usare waveTransmissive": e' **non
    riflettere un'onda di avvio**, e lo si puo' soddisfare anche togliendo
    l'onda invece del riflettore. Un test che nomina l'implementazione
    impedisce proprio la correzione che serve.
    """
    p = leggi(caso, "0", "p")
    u = leggi(caso, "0", "U")
    sbocco = blocco(p, "uscita")
    riflette = "fixedValue" in sbocco or "fixedMean" in sbocco
    if riflette:
        #: se il contorno riflette, l'onda NON deve nascere: gli ingressi
        #: devono salire in rampa. Le due cose sono ammesse solo insieme.
        for patch in ("ingresso_aria", "ingresso_gpl"):
            assert "table" in blocco(u, patch), (
                "lo sbocco riflette (fixedValue/fixedMean) ma gli ingressi "
                "partono a gradino: e' la combinazione che ha fatto esplodere "
                "le corse 1 e 2 a un transito acustico esatto.")
    else:
        assert "waveTransmissive" in sbocco, (
            "lo sbocco non riflette e non e' nemmeno non-riflettente "
            "dichiarato: che condizione e'?")


def test_lo_sbocco_ha_tutti_i_dati_che_la_condizione_chiede(caso):
    """Una condizione allo sbocco a cui manca un dato non parte, e il messaggio
    di OpenFOAM non e' fra i piu' chiari. Le chiavi richieste dipendono da
    quale condizione si e' scelta, quindi si guarda quella che c'e'."""
    sbocco = blocco(leggi(caso, "0", "p"), "uscita")
    richieste = {
        "waveTransmissive": ("field", "psi", "gamma", "fieldInf", "lInf", "value"),
        "fixedMean": ("meanValue", "value"),
        "fixedValue": ("value",),
    }
    for tipo, chiavi in richieste.items():
        if tipo in sbocco:
            for chiave in chiavi:
                assert chiave in sbocco, (
                    f"manca '{chiave}': {tipo} non parte senza.")
            return
    raise AssertionError(f"condizione allo sbocco non riconosciuta:\n{sbocco}")


def test_le_specie_rientranti_portano_dentro_valori_sensati(caso):
    """Con `zeroGradient` un ricircolo che rientra dallo sbocco porta dentro
    qualunque cosa si trovi fuori. Non dovrebbe capitare, ma se capita la
    risposta dev'essere sensata lo stesso."""
    for campo in ("T", "C3H8", "O2", "N2"):
        blocco = leggi(caso, "0", campo).split("uscita", 1)[1].split("}", 1)[0]
        assert "inletOutlet" in blocco, f"{campo} allo sbocco e' ancora zeroGradient"


# --------------------------------------------------------------------------- #
# 2. i limiti numerici esistono e si dichiarano
# --------------------------------------------------------------------------- #
def test_i_limiti_ci_sono_e_sono_larghi(caso):
    """I limiti tagliano valori che NON ESISTONO: in un caso non reattivo con
    ingressi a 283 e 293 K e 152 m/s, la temperatura fisica sta fra 270 e 310 K.
    Se fossero stretti taglierebbero fisica invece che numerica."""
    fv = leggi(caso, "constant", "fvOptions")
    assert "limitTemperature" in fv and "limitVelocity" in fv
    minimo = float(voce(fv, "min"))
    massimo = float(voce(fv, "max"))
    assert minimo <= 260.0 and massimo >= 450.0, (
        f"i limiti ({minimo}-{massimo} K) sono troppo stretti: taglierebbero la "
        "soluzione invece degli artefatti")


def test_il_polinomio_termodinamico_non_lavora_sul_proprio_bordo(caso):
    """Tlow era a 200 K e la soluzione oscillava a 203: un'estrapolazione
    valutata sul proprio bordo e' il posto peggiore dove stare, perche' cp e
    comprimibilita' diventano rumorosi proprio dove il campo e' gia' malmesso.
    Sotto ci sono i limiti dichiarati, che sono un'altra cosa: si contano."""
    termo = leggi(caso, "constant", "thermophysicalProperties")
    tlow = float(voce(termo, "Tlow"))
    assert tlow <= 180.0, f"Tlow = {tlow} K e' troppo vicino alla zona di lavoro"


# --------------------------------------------------------------------------- #
# 3. lo schema non deve fabbricare il risultato
# --------------------------------------------------------------------------- #
def test_le_frazioni_di_massa_non_sono_in_upwind(caso):
    """L'upwind diffonde numericamente, e una diffusione numerica in un caso di
    MESCOLAMENTO produce esattamente il risultato che si sta cercando di
    misurare. Sarebbe la peggiore delle conferme."""
    schemi = leggi(caso, "system", "fvSchemes")
    riga = next(r for r in schemi.splitlines() if "div(phi,Yi_h)" in r)
    assert "upwind" not in riga, riga
    assert "limitedLinear" in riga


# --------------------------------------------------------------------------- #
# 4. il campionamento deve guardare dove serve
# --------------------------------------------------------------------------- #
def test_si_campiona_anche_il_piano_di_fine_mescolamento(caso):
    """E' l'unico piano su cui la correlazione di Holdeman fa una promessa
    verificabile. Campionarne altri sei e non quello sarebbe un modo elaborato
    di non rispondere."""
    d = quote_finte()
    controllo = leggi(caso, "system", "controlDict")
    assert "surfaces" in controllo and "meridiano_getto" in controllo
    x_fine = d["x_getti"] + d["L_mescolamento"]
    assert any(abs(x - x_fine) < 1e-9 for x in d["piani_x"])
    assert controllo.count("cuttingPlane") >= len(d["piani_x"]) + 2


def test_si_scrive_abbastanza_spesso_da_poter_ripartire(caso):
    """La macchina virtuale si e' spenta tre volte sotto il calcolo. Con
    `startFrom latestTime` la frequenza di scrittura E' la granularita' con cui
    una corsa interrotta si riprende: a 40 microsecondi si perdevano quaranta
    minuti, a 10 se ne perdono quattro."""
    controllo = leggi(caso, "system", "controlDict")
    assert "latestTime" in controllo
    intervallo = float(voce(controllo, "writeInterval"))
    fine = float(voce(controllo, "endTime"))
    assert intervallo <= fine / 8.0, (
        f"scritture ogni {intervallo:.1e} s su una corsa di {fine:.1e}: una "
        "interruzione costa troppo")


# --------------------------------------------------------------------------- #
# Corsa 4: la portata imposta al posto della velocita' imposta.
#
# I tre test che seguono corrispondono a un difetto reale che ha fatto perdere
# tre corse. Il difetto era invisibile: il caso girava, la mesh era sana,
# checkMesh diceva "Mesh OK", e il solutore produceva numeri per ore prima di
# fermarsi. Erano solo i numeri di un altro motore.
# --------------------------------------------------------------------------- #
def test_gli_ingressi_impongono_la_portata_non_la_velocita(caso):
    """Con `U: fixedValue` + `p: zeroGradient` il flusso di massa entrante vale
    rho U S con rho = psi p: se p scende entra meno massa, e se entra meno massa
    p scende ancora. E' un anello di reazione positivo, e su questo caso si e'
    innescato tre volte sullo stesso labbro di foro, portando una cella a 1e-4
    bar e il passo temporale a 3 picosecondi."""
    u = (caso / "0" / "U").read_text()
    for ingresso in ("ingresso_aria", "ingresso_gpl"):
        blocco = u.split(ingresso, 1)[1].split("}", 1)[0]
        assert "flowRateInletVelocity" in blocco, (
            f"{ingresso} torna a imporre la velocita': l'anello p->rho->mdot->p "
            "si richiude e la corsa si ripianta.")
        assert "massFlowRate" in blocco
        #: `rho rho` usa il CAMPO di densita'. Con un `rhoInlet` fisso il flusso
        #: di massa tornerebbe a dipendere da p, cioe' al problema di prima con
        #: un altro nome.
        assert "rho             rho;" in blocco


def test_la_portata_imposta_e_quella_del_dimensionamento(caso):
    """Il difetto vero: le aree del sintetizzatore sono GEOMETRICHE (contengono
    gia' Cd), le velocita' sono IDEALI (attraverso l'area efficace Cd*A).
    Imporre le seconde sulle prime immette mdot/Cd, cioe' il 33 % di portata in
    piu' su entrambi i flussi."""
    d = quote_finte()
    u = (caso / "0" / "U").read_text()
    atteso = {"ingresso_aria": d["mdot_aria"] / d["n_getti"],
              "ingresso_gpl": d["mdot_gpl"] / d["n_getti"]}
    for ingresso, m in atteso.items():
        b = u.split(ingresso, 1)[1].split("massFlowRate", 1)[1]
        #: la portata ora e' una TABELLA (rampa d'avvio): il valore che conta e'
        #: quello di regime, cioe' l'ultimo. Leggere la prima cifra che capita
        #: prenderebbe lo zero iniziale della rampa e il test passerebbe
        #: sempre - contro un motore che non soffia.
        numeri = re.findall(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", b.split(");", 1)[0])
        letto = float(numeri[-1])
        assert letto == pytest.approx(m, rel=1e-6), (
            f"{ingresso}: la CFD soffia {letto/m:.3f} volte la portata di "
            "progetto. Non e' un dettaglio: sposta pressioni, velocita' in "
            "camera e tempo di residenza.")


def test_la_pressione_e_limitata_e_la_velocita_no(caso):
    """`limitVelocity` taglia U senza correggere phi ne' p: in una corsa che
    diverge non contiene il problema, lo nasconde. Il limite va messo dove la
    variabile diverge davvero, cioe' sulla pressione."""
    #: SI TOLGONO I COMMENTI PRIMA DI CERCARE. Il dizionario spiega a lungo
    #: perche' `limitVelocity` e' stato tolto, e la spiegazione contiene la
    #: parola: cercarla nel testo grezzo fa fallire il test proprio quando la
    #: cosa e' fatta bene. Un test che punisce la documentazione e' un test
    #: rotto, non un progetto disciplinato.
    def senza_commenti(testo: str) -> str:
        return "\n".join(r for r in testo.splitlines()
                          if not r.lstrip().startswith("//"))

    fv = senza_commenti((caso / "constant" / "fvOptions").read_text())
    assert "limitVelocity" not in fv, (
        "il limitatore di velocita' e' tornato: mascherava la divergenza "
        "tenendo |U| inchiodata a 400 m/s mentre p andava a zero indisturbata.")
    assert "limitTemperature" in fv, (
        "il limite di temperatura invece serve e va tenuto: e' su una "
        "grandezza fisicamente limitata da 270 a 310 K in questo caso.")
    fs = senza_commenti((caso / "system" / "fvSolution").read_text())
    assert "pMin" in fs and "pMax" in fs, (
        "senza pMin nulla impedisce a rho = psi*p di andare a zero in una "
        "cella e portarsi dietro il passo temporale di tutte le altre.")
    #: il limite deve stare BEN SOTTO qualunque pressione fisica del problema,
    #: altrimenti non e' una rete di sicurezza ma un vincolo che falsa il campo.
    riga = [r for r in fs.splitlines() if "pMin" in r][0]
    assert float(riga.split()[-1].rstrip(";")) <= 0.25 * quote_finte()["p_c"]


# --------------------------------------------------------------------------- #
# Il tempo di lavaggio: il numero che decide QUANDO un risultato e' un
# risultato. Senza di lui la corsa 3 e' stata letta per ore come se dicesse
# qualcosa sul mescolamento in camera, dove non aveva ancora rinnovato il 5 %
# del fluido.
# --------------------------------------------------------------------------- #
def test_il_tempo_di_lavaggio_torna_in_forma_chiusa_dove_l_area_e_costante():
    """Nel tratto d'anello l'area e' costante, quindi il lavaggio DEVE valere
    L/V esattamente. Se l'integrale non torna li', non torna da nessuna parte."""
    from sezioni import tempo_di_lavaggio
    d = quote_finte()
    A = math.pi * (d["R_anello"] ** 2 - d["R_getti"] ** 2)
    V = (d["mdot_gpl"] + d["mdot_aria"]) / (d["rho_aria"] * A)
    atteso = abs(d["x_ingresso"]) / V
    assert tempo_di_lavaggio(d, 0.0) == pytest.approx(atteso, rel=5e-3)


def test_il_tempo_di_lavaggio_cresce_e_esplode_dopo_il_gradino():
    """Il salto di area a x = 0 (anello -> camera, 11 volte) e' il motivo per
    cui un verdetto unico su tutti i piani non ha senso: se questo test smette
    di passare, e' cambiata la geometria e vanno rifatti i conti sui tempi."""
    from sezioni import tempo_di_lavaggio
    d = quote_finte()
    t = [tempo_di_lavaggio(d, x) for x in (-2.5e-3, 0.0, 3.0e-3, 12.0e-3)]
    assert all(b > a for a, b in zip(t, t[1:])), "il lavaggio deve crescere con x"
    #: i primi 12 mm di camera costano piu' di dieci volte tutto l'anello
    assert (t[3] - t[1]) > 10.0 * t[1]


def test_un_piano_immaturo_non_produce_un_verdetto(caso):
    """Il difetto che questo blocca: leggere l'ultimo fotogramma di un piano
    che il flusso non ha ancora lavato, e chiamarlo 'disuniformita' della
    miscela'. Il numero e' vero; non e' una risposta."""
    import sezioni
    assert sezioni.LAVAGGI_MINIMI >= 3.0, (
        "meno di tre lavaggi non e' una soglia, e' una speranza.")


def test_i_contorni_si_misurano_mentre_gira(caso):
    """"La simulazione gira con le pressioni giuste, da ogni lato?" era una
    domanda senza risposta: sugli ingressi non si impone nessuna pressione, si
    impone la portata, e quella che si sviluppa nessuno la guardava. Adesso e'
    nel log a ogni scrittura."""
    cd = (caso / "system" / "controlDict").read_text()
    for patch in ("ingresso_aria", "ingresso_gpl", "uscita"):
        assert f"p_{patch}" in cd, (
            f"nessuna misura di pressione su {patch}: il salto d'iniezione "
            "resta un numero del modello 0-D che nessuno confronta con la CFD.")
        assert f"portata_{patch}" in cd, (
            f"nessuna misura di portata su {patch}: se "
            "`flowRateInletVelocity` non facesse il suo mestiere non si "
            "vedrebbe da nessuna parte.")
    #: la portata si misura come somma di phi, che E' il flusso di massa: una
    #: media d'area di U non lo sarebbe, e sarebbe l'errore facile da fare.
    assert cd.count("operation       sum;") >= 3


def test_i_piani_di_misura_non_si_ripetono():
    """Il difetto: `x_getti` vale ESATTAMENTE meno la lunghezza di
    mescolamento, quindi il piano a f = 1.0 cade su x = 0, che era gia' in
    coda alla lista. Uscivano due `sez_XX` con le stesse coordinate, campionati
    due volte e stampati in tabella come due righe diverse.

    Non falsava nessun numero - erano identici - ma due righe uguali in una
    tabella si leggono come due misure che si confermano a vicenda, ed erano
    la stessa misura scritta due volte. La CFD costa ore: campionare due volte
    lo stesso piano e' anche disco e tempo buttati."""
    pytest.importorskip("cantera")
    from mesh_iniettore import quote_dal_progetto
    piani = quote_dal_progetto()["piani_x"]
    for i, x in enumerate(piani):
        for y in piani[i + 1:]:
            assert abs(x - y) > 1e-9, (
                f"due piani di misura coincidono a x = {x*1e3:.3f} mm")


def test_il_piano_di_fianco_non_cade_sul_bordo_di_simmetria(caso):
    """IL DIFETTO CHE HA PRODOTTO MEZZA FIGURA BIANCA PER TRE CORSE.

    Il settore ha i due fianchi dichiarati `symmetryPlane`, e uno di essi E' il
    piano z = 0. Un `cuttingPlane` con normale (0 0 1) passante per l'origine ci
    cade esattamente sopra: non ci sono celle da intersecare e OpenFOAM
    restituisce quattro triangoli degeneri. Nessun errore, nessun avviso: solo
    un piano di misura vuoto, in una tabella accanto a uno pieno.

    Il piano va spostato dentro al dominio - poco, perche' deve restare "meta'
    strada fra due getti", ma abbastanza da avere celle sotto."""
    d = quote_finte()
    cd = (caso / "system" / "controlDict").read_text()
    blocco = cd.split("meridiano_fianco", 1)[1].split("}", 1)[0]
    #: si cerca "point (" e non "point": la riga `planeType pointAndNormal;`
    #: contiene la parola e viene prima. Un test che pesca la riga sbagliata
    #: fallisce per il motivo sbagliato, che e' peggio di non averlo.
    riga = [r for r in blocco.splitlines() if "point (" in r][0]
    z = float(riga.split("(")[1].split(")")[0].split()[2])
    assert z > 0.0, (
        "il piano di fianco e' tornato sull'origine, cioe' sul bordo di "
        "simmetria: uscira' vuoto e nessuno se ne accorgera'.")
    #: e deve restare meta' strada: lo scarto angolare al raggio d'iniezione
    #: va confrontato con i gradi che separano il piano dal getto.
    ang = math.degrees(math.asin(z / d["R_getti"]))
    assert ang < 0.25 * (180.0 / d["n_getti"]), (
        f"spostato di {ang:.1f} gradi: non e' piu' meta' strada fra due getti.")


def test_il_foro_del_getto_non_scava_una_tasca_nel_metallo():
    """IL DIFETTO: un "+2 mm" giusto per un taglio, applicato a una fusione.

    Il cilindro del getto veniva costruito lungo (R_anello - r_bore) + 2 mm e
    FUSO con il dominio fluido. L'eccedenza di 2 mm serve a non lasciare facce
    coincidenti nell'operazione booleana - e per un TAGLIO e' corretta, perche'
    finisce dentro il solido. Per una FUSIONE no: l'eccedenza diventa fluido, e
    davanti a ogni getto restava una tasca cieca profonda 2 mm scavata nel
    metallo, che nel pezzo vero non c'e'.

    Misurata sulla corsa 4 conteneva Y_GPL = 0.0000 su 927 punti, quindi non
    aveva falsato il mescolamento. Si toglie perche' e' geometria che il motore
    non ha - non perche' avesse fatto danno."""
    from mesh_iniettore import tratto_del_getto
    r_bore, R_getti, R_anello = 1.646e-3, 2.646e-3, 4.300e-3
    r0, L = tratto_del_getto(r_bore, R_getti, R_anello)
    #: PARTE ESATTAMENTE dalla parete del plenum. Questa riga diceva
    #: `assert r0 < r_bore` - "il foro deve pescare nel GPL" - e la
    #: motivazione era sbagliata: il plenum NON e' nel dominio fluido (il
    #: raggio sotto R_getti viene tolto), quindi la base del cilindro non
    #: pesca in niente: E' la faccia d'ingresso del combustibile. Spostarla
    #: piu' in dentro allungava il foro e faceva sparire la patch.
    #: Il test verificava dunque il difetto, non il requisito.
    assert r0 == pytest.approx(r_bore)
    #: sconfina sopra la parete interna: il foro deve sboccare nell'anello
    assert r0 + L > R_getti
    #: ma NON arriva alla parete esterna: oltre quella c'e' metallo
    assert r0 + L < R_anello, (
        f"il cilindro arriva a {(r0+L)*1e3:.3f} mm con la parete esterna a "
        f"{R_anello*1e3:.3f}: e' tornata la tasca cieca.")


def test_una_geometria_incoerente_si_rifiuta_invece_di_costruire():
    """Se l'anello fosse cosi' sottile che lo sconfinamento lo attraversa,
    la funzione tornerebbe a scavare nel metallo. Meglio fermarsi."""
    from mesh_iniettore import tratto_del_getto
    with pytest.raises(ValueError, match="tasca cieca"):
        tratto_del_getto(1.646e-3, 2.646e-3, 2.746e-3)
    with pytest.raises(ValueError, match="incoerenti"):
        tratto_del_getto(3.0e-3, 2.646e-3, 4.300e-3)


def blocco(testo: str, nome: str) -> str:
    """Il corpo del blocco `nome` di un dizionario OpenFOAM.

    Cercare per sottostringa non basta: `vorticity` compare anche nell'elenco
    dei campi campionati, e si finirebbe a leggere le impostazioni del blocco
    sbagliato. Qui si parte dalla riga che E' il nome e si conta le graffe."""
    righe = testo.splitlines()
    for i, r in enumerate(righe):
        if r.strip() == nome and righe[i + 1].strip().startswith("{"):
            liv, pezzo = 0, []
            for r2 in righe[i + 1:]:
                liv += r2.count("{") - r2.count("}")
                pezzo.append(r2)
                if liv == 0:
                    return "\n".join(pezzo)
    raise AssertionError(f"blocco '{nome}' assente")


def test_le_sonde_di_bordo_scrivono_una_curva_non_due_punti(caso):
    """IL DIFETTO CHE QUESTO TEST BLOCCA. Nel run 4 le sonde di bordo erano
    scritte con `adjustableRunTime`: in 250 microsecondi hanno prodotto DUE
    campioni, a 20 e 50 us. Con due punti non si dimostra che una pressione si
    e' assestata, e la domanda "gira con le pressioni giuste?" e' rimasta
    senza risposta misurata per una corsa intera - mentre sembrava di averla
    strumentata.

    `timeStep` non dipende dalla scrittura dei volumi, sopravvive ai riavvii e
    non chiede al solutore di centrare istanti esatti (cosa che gli deforma il
    passo temporale proprio mentre lo si sta misurando)."""
    cd = (caso / "system" / "controlDict").read_text()
    for nome in ("p_uscita", "portata_uscita", "p_ingresso_aria",
                 "massa_nel_dominio"):
        b = blocco(cd, nome)
        assert "writeControl    timeStep;" in b, (
            f"{nome} non campiona a passi: con `adjustableRunTime` il run 4 ha "
            "prodotto due soli campioni in tutta la corsa.")
        n = int(voce(b, "writeInterval"))
        assert 1 <= n <= 100, (
            f"{nome} campiona ogni {n} passi: troppo rado per essere una curva "
            "di convergenza.")


def test_la_massa_contenuta_si_misura(caso):
    """PERCHE' NON BASTA LO SBILANCIO DELLE PORTATE. Portata entrante e uscente
    sono due numeri grandi e quasi uguali: la loro differenza, vicino alla
    convergenza, e' sotto il rumore del calcolo. L'integrale di volume di rho
    e' invece la grandezza che il transitorio di riempimento sta cambiando, e
    la sua derivata va a zero solo quando il dominio ha finito di riempirsi.
    E' questa la curva che stabilisce quando le pressioni sono giuste."""
    b = blocco((caso / "system" / "controlDict").read_text(), "massa_nel_dominio")
    assert "volFieldValue" in b and "volIntegrate" in b, (
        "la massa contenuta va integrata sul volume, non mediata su un bordo.")
    assert "(rho)" in b


def test_la_vorticita_esiste_quando_i_piani_la_campionano(caso):
    """IL MODO SILENZIOSO DI SBAGLIARE. Il blocco `superfici` campiona la
    vorticita' molto piu' spesso di quanto si scrivano i volumi. Se la
    vorticita' fosse calcolata solo a `writeTime`, quasi tutti i campionamenti
    troverebbero un campo che non esiste nel registro - e OpenFOAM non
    solleverebbe un errore: scriverebbe piani senza quel campo, e la figura
    uscirebbe vuota senza che niente lo segnali. E' esattamente il modo in cui
    il piano di fianco e' uscito vuoto per tre corse di fila."""
    cd = (caso / "system" / "controlDict").read_text()
    v = blocco(cd, "vorticity")
    assert "executeControl  timeStep;" in v and "executeInterval 1;" in v, (
        "la vorticita' non viene calcolata a ogni passo: i piani la "
        "campionerebbero quando non c'e'.")
    sup = blocco(cd, "superfici")
    assert "vorticity" in voce(sup, "fields"), (
        "si calcola la vorticita' e non la si campiona: lavoro sprecato.")


def test_le_medie_temporali_non_partono_dentro_il_transitorio(caso):
    """Mediare durante il riempimento significa mediare un motore che non
    esiste: il campo passa da fermo a regime, e la media pesa allo stesso modo
    i due. Il valore di default e' un terzo della corsa, che e' una scelta
    prudente e comunque sovrascrivibile - ma NON puo' essere zero."""
    m = blocco((caso / "system" / "controlDict").read_text(), "medie")
    t0 = float(voce(m, "timeStart"))
    fine = float(voce((caso / "system" / "controlDict").read_text(), "endTime"))
    assert t0 > 0.0, "le medie partono da zero: dentro il transitorio."
    assert t0 <= 0.5 * fine, (
        f"le medie partono a {t0*1e6:.0f} us su {fine*1e6:.0f} di corsa: resta "
        "troppo poca finestra per una media temporale sensata.")
    assert "prime2Mean on" in m, (
        "senza prime2Mean non si sa quanta unsteadiness sia RISOLTA, e non si "
        "puo' dire se il campo URANS sia di fatto stazionario o stia oscillando.")


def test_il_condotto_del_getto_ha_la_lunghezza_del_pezzo_vero():
    """DUE DIFETTI, NELLO STESSO PARAMETRO, IN DUE CORREZIONI DI FILA.

    Il primo: `lung = (R_anello - r_bore) + 2 mm`, che sfondava di 2 mm oltre
    la parete esterna e lasciava una tasca cieca nel metallo.

    Il secondo, introdotto CORREGGENDO il primo: lo sconfinamento messo anche
    all'estremita' interna, a r_bore - 0.3 mm. Sembra la scelta simmetrica e
    prudente, e invece allunga il condotto del 60 %: il foro vero va dalla
    parete del plenum a quella dell'anello, e il suo L/d - da cui dipende se
    il getto riattacca dopo la vena contracta - passava da 1.37 a 2.19.

    Questo test fissa la lunghezza al valore geometrico, non a un valore
    comodo: l'unico sconfinamento ammesso e' quello esterno, che finisce
    dentro fluido gia' esistente e non aggiunge volume."""
    from mesh_iniettore import SCONFINAMENTO_GETTO, tratto_del_getto
    r_bore, R_getti, R_anello = 1.646e-3, 2.646e-3, 4.300e-3
    r0, lung = tratto_del_getto(r_bore, R_getti, R_anello)

    assert r0 == pytest.approx(r_bore), (
        f"il condotto parte da {r0*1e3:.3f} mm invece che dalla parete del "
        f"plenum a {r_bore*1e3:.3f}: e' un foro piu' lungo di quello stampato.")
    spessore = R_getti - r_bore
    assert lung == pytest.approx(spessore + SCONFINAMENTO_GETTO), (
        "la lunghezza non e' lo spessore da attraversare piu' il solo "
        "sconfinamento esterno.")
    #: e l'eccedenza resta dentro fluido gia' esistente, cioe' non scava niente
    assert r0 + lung > R_getti, "senza sconfinare, le facce coincidono"
    assert r0 + lung < R_anello, "sconfina oltre la parete esterna: tasca cieca"


def test_la_faccia_d_ingresso_del_gpl_sta_dove_la_si_cerca():
    """IL MODO IN CUI IL DIFETTO E' RIUSCITO A PASSARE INOSSERVATO.

    `_marca_bordi` riconosce l'ingresso del combustibile per POSIZIONE: e' la
    faccia a raggio r_bore vicino a x_getti. Se la base del cilindro si sposta,
    la faccia non viene piu' riconosciuta, finisce nel gruppo 'pareti', e la
    mesh esce con cinque bordi invece di sei - **senza errori**: checkMesh dice
    "Mesh OK", perche' geometricamente lo e'. Il combustibile semplicemente non
    entra piu' nel motore.

    Questo test lega le due cose che devono restare d'accordo: dove la
    geometria mette la faccia, e dove la classificazione la cerca."""
    from mesh_iniettore import tratto_del_getto
    r_bore, R_getti, R_anello = 1.646e-3, 2.646e-3, 4.300e-3
    r0, _ = tratto_del_getto(r_bore, R_getti, R_anello)
    #: la tolleranza e' quella scritta in `_marca_bordi`
    assert abs(r0 - r_bore) < 5e-5, (
        f"la base del getto sta a {r0*1e3:.3f} mm ma la classificazione cerca "
        f"l'ingresso GPL entro 0.05 mm da {r_bore*1e3:.3f}: la patch "
        "'ingresso_gpl' non verrebbe creata e il GPL non entrerebbe.")


def test_la_pressione_media_allo_sbocco_e_ancorata(caso):
    """COME E' MORTA LA CORSA 5, e perche' non poteva non morire.

    Con `flowRateInletVelocity` su ENTRAMBI gli ingressi, nessun contorno
    impone una pressione: gli ingressi impongono massa, le pareti e i piani di
    simmetria sono a gradiente nullo. L'unico vincolo sul LIVELLO di pressione
    del dominio era il rilassamento di `waveTransmissive` verso `fieldInf` su
    lInf = 0.05 m - 2.7 volte la lunghezza del dominio, cioe' quasi nessun
    vincolo. Il commento nel generatore lo diceva gia': "piu' e' lunga, meno
    vincola la pressione media".

    Per 300 us ha retto. Poi lo squilibrio di massa ha vinto: camera da 6.53 a
    3.38 bar, alimentazione da 7.64 a 13.40, il getto ha strozzato e spingere
    una portata FISSA attraverso una gola strozzata a densita' calante ha
    portato la pressione a monte a crescere senza limite. Il 9 % delle celle
    contro il limitatore a 500 K.

    Questo test impone che il dominio abbia SEMPRE almeno un punto in cui la
    pressione e' ancorata."""
    p = (caso / "0" / "p").read_text()
    assert "fixedMean" in p, (
        "nessuna condizione ancora la pressione media: con la portata imposta "
        "su tutti gli ingressi il livello di pressione e' libero di andare "
        "alla deriva, ed e' esattamente cosi' che e' morta la corsa 5.")
    #: e deve essere ancorata al valore del ciclo, non a un numero qualsiasi
    b = blocco(p, "uscita")
    assert abs(float(voce(b, "meanValue")) - quote_finte()["p_c"]) < 1.0, (
        "lo sbocco e' ancorato a una pressione che non e' quella di camera.")


def test_gli_ingressi_partono_in_rampa_e_non_a_gradino(caso):
    """L'ALTRA META' DELLO STESSO GUASTO.

    Lo sbocco trasparente non era un capriccio: serviva a non riflettere
    l'onda di compressione dell'avvio impulsivo, che aveva fatto esplodere le
    corse 1 e 2 a 56.6 e 56.2 us - un transito acustico esatto del dominio.
    Ancorare la pressione senza togliere quell'onda significherebbe tornare a
    farla rimbalzare.

    Un gradino di portata e' un pistone che parte di colpo: genera
    dp ~ rho*c*V = 7.87*343*114 = 3.1 bar. Con una rampa lunga qualche transito
    acustico (L/c = 54 us) il campo si adatta quasi-staticamente e l'onda
    coerente non nasce. Le due correzioni stanno o cadono insieme: questo test
    impedisce di rimettere il gradino lasciando `fixedMean`."""
    from caso_mescolamento import T_RAMPA
    u = (caso / "0" / "U").read_text()
    for patch in ("ingresso_aria", "ingresso_gpl"):
        b = blocco(u, patch)
        assert "table" in b, (
            f"{patch} parte a gradino: con lo sbocco che ora riflette, l'onda "
            "d'avvio tornerebbe indietro tutta intera.")
        assert "(0            0)" in b, (
            f"{patch} non parte da portata nulla: non e' una rampa.")
    #: la rampa deve durare piu' di un transito acustico, altrimenti non serve
    L_dominio = 18.6e-3
    c_suono = 343.0
    assert T_RAMPA > L_dominio / c_suono, (
        f"la rampa dura {T_RAMPA*1e6:.0f} us contro un transito acustico di "
        f"{L_dominio/c_suono*1e6:.0f} us: troppo corta per togliere il gradino.")


def test_la_rampa_arriva_esattamente_alla_portata_di_progetto():
    """Una rampa che sbaglia il valore finale non si vedrebbe: il transitorio
    d'avvio la nasconderebbe, e il motore girerebbe per sempre alla portata
    sbagliata. Si controlla il valore, non solo la forma."""
    from caso_mescolamento import _rampa
    testo = _rampa(8.088459291e-03, 1.0e-4)
    assert "8.088459291e-03" in testo
    #: e resta a quel valore anche molto dopo la fine della rampa
    assert testo.count("8.088459291e-03") == 2, (
        "la tabella non tiene il valore dopo la rampa: `table` interpola fra i "
        "punti dati, e senza un punto finale lontano il comportamento oltre "
        "l'ultimo istante dipende dal clamping, non dalla nostra intenzione.")
    with pytest.raises(ValueError):
        _rampa(1.0e-3, 0.0)
