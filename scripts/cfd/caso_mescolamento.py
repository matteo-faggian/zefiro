#!/usr/bin/env python3
"""Genera il caso OpenFOAM del mescolamento aria/GPL nell'iniettore.

NON REATTIVO, e la scelta e' deliberata: si separa il mescolamento dalla
combustione, perche' mettendoli insieme e ottenendo un risultato brutto non si
saprebbe quale dei due modelli incolpare. La domanda a cui questo caso deve
rispondere e' una sola:

    il GPL diventa davvero uniforme entro la lunghezza di mescolamento che la
    correlazione di Holdeman promette?

Se la risposta e' no, il C = 2.5 su cui poggia tutto l'iniettore non vale
niente e la geometria va rifatta.

LE PROPRIETA' DEI GAS VENGONO DA CANTERA, non da tabelle ricopiate: i
coefficienti NASA a 7 termini si estraggono dallo stesso meccanismo (gri30) che
il resto del progetto usa per la termochimica. Cosi' non esistono due verita'
sulle stesse specie, che e' il modo piu' silenzioso di far divergere due parti
di un modello.
"""
from __future__ import annotations

import argparse
import math
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

SPECIE = ("C3H8", "O2", "N2")
#: Composizione dell'ARIA in frazioni di massa, coerente con units.AIR_MOLE_FRACTIONS.
Y_ARIA = {"O2": 0.23292, "N2": 0.76708, "C3H8": 0.0}
#: Intensita' di turbolenza agli ingressi. 5 % e' il valore convenzionale per un
#: condotto; qui a monte c'e' una contrazione 6:1, che la turbolenza la SCHIACCIA,
#: quindi 5 % e' prudente in eccesso.
INTENSITA = 0.05


def intestazione(cls: str, oggetto: str, loc: str) -> str:
    return f"""FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    location    "{loc}";
    object      {oggetto};
}}
"""


def _termo_da_cantera() -> str:
    """Coefficienti NASA-7 e Sutherland per specie, estratti da gri30.

    IL TRASPORTO E' SUTHERLAND, e non e' una preferenza: questo OpenFOAM
    accetta `janaf` solo accoppiato a `sutherland` (con `const` vuole `hConst`,
    cioe' cp costante, che a 2200 K sarebbe sbagliato del 30 %). Le due
    costanti As e Ts NON sono prese da una tabella: si ricavano imponendo che
    la formula di Sutherland passi per le viscosita' che Cantera da' a 300 e a
    800 K, cioe' agli estremi dell'intervallo di questo caso.

        mu = As sqrt(T) / (1 + Ts/T) = As T^1.5 / (T + Ts)

    Fuori da quell'intervallo l'estrapolazione peggiora, ed e' il motivo per
    cui il caso reattivo (fino a 2200 K) andra' rifittato prima di girare.

    Tlow e' abbassata a 150 K contro i 300 dichiarati da gri30. Gli ingressi
    stanno a 293 e 283 K, cioe' 7 e 17 K SOTTO il limite del polinomio, e senza
    abbassarla OpenFOAM si ferma. E' un'estrapolazione, e va detta: su specie
    stabili il cp del ramo basso e' quasi lineare e 17 K costano meno dello
    0.5 %, ma resta un'estrapolazione e non una valutazione.

    PERCHE' 150 E NON 200. Con 200 il limite del polinomio cadeva dentro
    l'intervallo in cui la soluzione oscillava (la prima corsa e' scesa a 203 K
    per sottoelongazione numerica), e un'estrapolazione valutata sul proprio
    bordo e' il posto peggiore dove stare: cp e psi diventano rumorosi proprio
    dove il campo e' gia' malmesso. Portandola a 150 il bordo si allontana dalla
    zona di lavoro, e la temperatura viene invece TAGLIATA a 250 K da un limite
    dichiarato in `fvOptions` - che e' una cosa diversa e onesta, perche' si
    conta quante celle tocca.
    """
    import cantera as ct

    gas = ct.Solution("gri30.yaml")
    blocchi = []
    for nome in SPECIE:
        sp = gas.species(nome)
        c = sp.thermo.coeffs        # [T_mid, alte(7), basse(7)]
        T_mid, alte, basse = c[0], c[1:8], c[8:15]
        #: fit di Sutherland su due punti, 300 e 800 K
        T1, T2 = 300.0, 800.0
        gas.TPX = T1, 101325.0, {nome: 1.0}
        mu1 = gas.viscosity
        gas.TPX = T2, 101325.0, {nome: 1.0}
        mu2 = gas.viscosity
        a1, a2 = T1 ** 1.5, T2 ** 1.5
        Ts = (mu2 * a1 * T2 - mu1 * a2 * T1) / (mu1 * a2 - mu2 * a1)
        As = mu1 * (T1 + Ts) / a1
        blocchi.append(f"""
{nome}
{{
    specie
    {{
        molWeight   {sp.molecular_weight:.6f};
    }}
    thermodynamics
    {{
        Tlow 150.0;
        Thigh       {sp.thermo.max_temp:.1f};
        Tcommon     {T_mid:.1f};
        highCpCoeffs ( {' '.join(f'{v:.9e}' for v in alte)} );
        lowCpCoeffs  ( {' '.join(f'{v:.9e}' for v in basse)} );
    }}
    transport
    {{
        As          {As:.6e};
        Ts          {Ts:.4f};
    }}
}}""")
    return "\n".join(blocchi)


#: Durata della rampa di portata all'avvio, e il numero si ricava, non si
#: sceglie a occhio.
#:
#: Un avvio a gradino e' un pistone che parte di colpo a V: genera un'onda di
#: compressione di ampiezza
#:      dp ~ rho c V = 7.87 * 343 * 114 = 3.1 bar
#: che e' l'ordine di grandezza dei picchi visti davvero nelle corse (9-11 bar
#: contro 6.36 di riferimento). Con uno sbocco che ancora la pressione, e
#: quindi RIFLETTE, un'onda cosi' non e' accettabile.
#:
#: La grandezza che conta e' il rapporto fra il tempo di salita e il tempo
#: acustico del dominio, L/c = 18.6 mm / 343 m/s = 54 us. Se la rampa dura
#: MOLTO PIU' di un transito, l'informazione attraversa il dominio molte volte
#: mentre la portata sale, il campo si adatta quasi-staticamente e l'onda
#: coerente non si forma: resta solo la pressurizzazione regolare fino al
#: salto d'iniezione di progetto, ~1 bar. Se dura molto MENO, si torna al
#: gradino.
#:
#: 100 us sono 1.9 transiti: abbastanza da togliere il gradino senza spendere
#: mezzo transitorio in rampa. Non e' gratis - i primi 100 us non sono motore,
#: sono avviamento - ed e' il motivo per cui `--media-da` va messo ben oltre.
T_RAMPA = 1.0e-4


def _rampa(mdot: float, t_rampa: float) -> str:
    """Portata che sale da zero al valore di progetto in `t_rampa`.

    PERCHE' ESISTE. Le corse 1 e 2 sono esplose a 56.6 e 56.2 us, cioe' a un
    transito acustico esatto del dominio: il campo partiva fermo e gli
    ingressi cominciavano a soffiare a piena portata nello stesso istante. Quel
    gradino genera un'onda di compressione che percorre i 18.6 mm, trova lo
    sbocco e torna indietro.

    La risposta di allora fu rendere lo sbocco TRASPARENTE (`waveTransmissive`
    con lInf lungo) perche' l'onda uscisse. Ha funzionato contro l'onda, e ha
    creato il guasto della corsa 5: un contorno trasparente non vincola la
    pressione media, e con la portata imposta su entrambi gli ingressi non
    restava nessun punto in cui la pressione fosse fissata.

    Questa e' la cura dell'altra meta' del problema: si toglie il GRADINO, e
    l'onda non nasce. Allora lo sbocco puo' tornare ad ancorare la pressione.
    """
    if t_rampa <= 0.0:
        raise ValueError("la rampa deve durare un tempo positivo")
    return (f"        massFlowRate    table\n"
            f"        (\n"
            f"            (0            0)\n"
            f"            ({t_rampa:.6e} {mdot:.9e})\n"
            f"            (1            {mdot:.9e})\n"
            f"        );")


def _controlli_di_contorno(ogni: int) -> str:
    """Sonde che misurano pressione, portata e massa contenuta nel dominio.

    IL CAMPIONAMENTO E' A PASSI, NON A TEMPO, e la differenza non e' di gusto.
    Nel run 4 questi stessi blocchi erano scritti con `adjustableRunTime`: in
    250 microsecondi hanno prodotto DUE campioni in tutto, a 20 e 50 us. Con
    due punti non si dimostra che una pressione si e' assestata - non si
    dimostra nemmeno in che verso stia andando - e la domanda "la simulazione
    gira con le pressioni giuste?" e' rimasta senza risposta misurata per una
    corsa intera. `timeStep` non dipende dalla scrittura dei volumi, non viene
    riazzerato dai riavvii e non chiede al solutore di centrare istanti esatti
    (che gli deforma il passo temporale): scrive e basta, ogni `ogni` passi.

    LA MASSA CONTENUTA E' LA SONDA VERA. Lo sbilancio fra portata entrante e
    uscente e' la differenza fra due numeri grandi e quasi uguali: al 2 % di
    convergenza e' gia' sotto il rumore del calcolo. L'integrale di volume di
    rho e' invece la grandezza che il transitorio di riempimento sta davvero
    cambiando, e la sua derivata va a zero SOLO quando il dominio ha finito di
    riempirsi. E' questa la curva che dice quando le pressioni sono giuste.
    """
    pezzi = []
    for patch in ("ingresso_aria", "ingresso_gpl", "uscita"):
        pezzi.append(f"""    p_{patch}
    {{
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        regionType      patch;
        name            {patch};
        operation       areaAverage;
        fields          (p T);
        writeControl    timeStep;
        writeInterval   {ogni:d};
        writeFields     no;
        log             no;
    }}
    portata_{patch}
    {{
        type            surfaceFieldValue;
        libs            (fieldFunctionObjects);
        regionType      patch;
        name            {patch};
        operation       sum;
        fields          (phi);
        writeControl    timeStep;
        writeInterval   {ogni:d};
        writeFields     no;
        log             no;
    }}""")
    pezzi.append(f"""    massa_nel_dominio
    {{
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        regionType      all;
        operation       volIntegrate;
        fields          (rho);
        writeControl    timeStep;
        writeInterval   {ogni:d};
        writeFields     no;
        log             no;
    }}
    pressione_media
    {{
        type            volFieldValue;
        libs            (fieldFunctionObjects);
        regionType      all;
        operation       volAverage;
        fields          (p T k);
        writeControl    timeStep;
        writeInterval   {ogni:d};
        writeFields     no;
        log             no;
    }}""")
    return "\n".join(pezzi)


def _controlli_di_turbolenza(media_da: float, passo_sup: float) -> str:
    """Quello che serve per GIUDICARE il campo turbolento, non solo per vederlo.

    Tre cose, e ognuna risponde a una domanda che altrimenti resterebbe
    un'opinione.

    `yPlus` - LE FUNZIONI DI PARETE SONO LECITE? k-omega SST con le wall
    function chiede y+ fuori dalla zona cuscinetto (5 < y+ < 30), o sotto 1 se
    si vuole risolvere il sottostrato. La mesh e' stata fatta senza misurare
    questo numero: se y+ cade in mezzo, l'attrito a parete - e quindi la
    perdita di pressione e la produzione di turbolenza vicino al muro - non e'
    ne' risolto ne' modellato, e va detto.

    `vorticity` - CALCOLATA DAL SOLUTORE, non ricostruita dopo. Finora la
    vorticita' delle figure veniva da un gradiente ai nodi fatto a mano sulla
    superficie di taglio, cioe' da un campo gia' interpolato una volta:
    derivare un'interpolazione amplifica il rumore dell'interpolazione. Qui il
    rotore si prende sulle celle, con lo stesso schema con cui il solutore
    calcola tutto il resto, e solo dopo si taglia il piano.

    `fieldAverage` - QUANTA UNSTEADINESS C'E' DAVVERO. In URANS convivono due
    turbolenze: quella MODELLATA (k, che il modello porta) e quella RISOLTA
    (le fluttuazioni che il campo medio compie nel tempo, che il modello NON
    sa di avere). prime2Mean misura la seconda. Se e' trascurabile rispetto a
    k il calcolo e' di fatto stazionario e ha senso leggerlo come tale; se non
    lo e', il getto sta oscillando e una media temporale e' obbligatoria prima
    di dire qualsiasi numero sul mescolamento. Le medie partono da `media_da`,
    cioe' DOPO l'assestamento: mediare dentro il transitorio di riempimento
    significa mediare un motore che non esiste.
    """
    return f"""    yPlus
    {{
        type            yPlus;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        log             yes;
    }}
    vorticity
    {{
        type            vorticity;
        libs            (fieldFunctionObjects);
        //: SI CALCOLA A OGNI PASSO ANCHE SE SI SCRIVE DI RADO. Il blocco
        //: `superfici` campiona la vorticita' ogni {passo_sup*1e6:.1f} us, e a
        //: quell'istante il campo deve ESISTERE nel registro: se lo si
        //: calcolasse solo quando si scrivono i volumi, il campionamento
        //: troverebbe il vuoto quasi sempre e fallirebbe in silenzio.
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
        log             no;
    }}
    medie
    {{
        type            fieldAverage;
        libs            (fieldFunctionObjects);
        timeStart       {media_da:.8f};
        //: la media si ACCUMULA a ogni passo e si SCRIVE quando si scrivono i
        //: volumi: accumularla di rado sarebbe una media di campioni, non una
        //: media temporale, e su un getto che oscilla a ~50 kHz darebbe un
        //: risultato che dipende da quando si e' guardato.
        executeControl  timeStep;
        executeInterval 1;
        writeControl    writeTime;
        restartOnRestart false;
        fields
        (
            U     {{ mean on; prime2Mean on; base time; }}
            p     {{ mean on; prime2Mean on; base time; }}
            T     {{ mean on; prime2Mean on; base time; }}
            C3H8  {{ mean on; prime2Mean on; base time; }}
            k     {{ mean on; prime2Mean off; base time; }}
        );
    }}"""


def scrivi(caso: Path, d: dict, tempo_finale: float, scrittura: float,
           n_proc: int, campionamento: float | None = None,
           binario: bool = False, passo_sonde: int = 20,
           media_da: float | None = None) -> None:
    for sub in ("0", "constant", "system"):
        (caso / sub).mkdir(parents=True, exist_ok=True)

    ang = math.pi / d["n_getti"]
    n = d["n_getti"]

    #: ------------------------------------------------------------------ #
    #: SI IMPONE LA PORTATA, NON LA VELOCITA'. Due ragioni, e nessuna delle
    #: due e' estetica.
    #:
    #: (1) COERENZA CON IL DIMENSIONAMENTO. Il sintetizzatore calcola le aree
    #:     con  A = mdot / (Cd rho V)  e Cd = 0.75: `d_getto` e l'altezza
    #:     dell'anello sono aree GEOMETRICHE, mentre `V_getto` e `V_aria` sono
    #:     le velocita' ideali attraverso l'area EFFICACE Cd*A. La mesh usa
    #:     l'area geometrica. Imporre la velocita' ideale su di essa immette
    #:     mdot/Cd = 1.333 volte la portata di progetto - misurato, non
    #:     stimato - su entrambi i flussi. Il motore simulato non era questo.
    #:
    #:     Imponendo invece la portata, le velocita' scendono di Cd su TUTTI E
    #:     DUE i flussi, quindi J e C - che sono rapporti - restano esattamente
    #:     quelli di progetto. Non si perde la fisica che interessa: si perde
    #:     solo un errore del 33 % sulla portata, sulla pressione e sul tempo
    #:     di residenza.
    #:
    #:     Che poi la velocita' giusta sia mdot/(rho A) e non quella ideale, per
    #:     l'anello e' certo: l'anello e' un CONDOTTO, l'aria lo riempie tutto,
    #:     non c'e' vena contratta. Per il getto il foro e' lungo 1.00 mm con
    #:     d = 0.731 mm (L/d = 1.37), cioe' al limite del riattacco: e' un TODO
    #:     da misurare, ed e' dichiarato come tale.
    #:
    #: (2) STABILITA'. Con `U: fixedValue` e `p: zeroGradient` il flusso di massa
    #:     entrante vale rho U S con rho = psi p preso dalla cella interna:
    #:     se p scende entra meno massa, e se entra meno massa p scende ancora.
    #:     E' un anello di reazione POSITIVO. Su questo caso si e' innescato sul
    #:     labbro del foro (x = -3.07 mm, r = 2.66 mm) e ha portato quella cella
    #:     da 2.86 bar a 1e-4 bar in 12 microsecondi: rho -> 0, Courant -> inf,
    #:     passo temporale a 3 picosecondi. Tre corse su tre si sono fermate li'.
    #:     `flowRateInletVelocity` fissa phi = mdot indipendentemente da p, e
    #:     l'anello si spezza alla radice.
    #: ------------------------------------------------------------------ #
    A_getto = math.pi / 4.0 * d["d_getto"] ** 2
    A_aria = math.pi * (d["R_anello"] ** 2 - d["R_getti"] ** 2) / n
    mdot_gpl_settore = d["mdot_gpl"] / n
    mdot_aria_settore = d["mdot_aria"] / n
    V_gpl_eff = mdot_gpl_settore / (d["rho_gpl"] * A_getto)
    V_aria_eff = mdot_aria_settore / (d["rho_aria"] * A_aria)

    #: la turbolenza d'ingresso si stima sulla velocita' EFFETTIVA, non su
    #: quella ideale: usare due velocita' diverse nello stesso contorno e' il
    #: genere di incoerenza che poi non si ritrova piu'.
    Dh_aria = 2.0 * d["H"]
    k_aria = 1.5 * (V_aria_eff * INTENSITA) ** 2
    om_aria = math.sqrt(k_aria) / (0.09 ** 0.25 * 0.07 * Dh_aria)
    k_gpl = 1.5 * (V_gpl_eff * INTENSITA) ** 2
    om_gpl = math.sqrt(k_gpl) / (0.09 ** 0.25 * 0.07 * d["d_getto"])
    J_eff = ((d["rho_gpl"] * V_gpl_eff ** 2)
             / (d["rho_aria"] * V_aria_eff ** 2))
    print(f"portate imposte per settore: aria {mdot_aria_settore*1e3:.4f} g/s, "
          f"GPL {mdot_gpl_settore*1e3:.4f} g/s")
    print(f"velocita' che ne seguono:    aria {V_aria_eff:.2f} m/s "
          f"(ideale {d['V_aria']:.2f}), GPL {V_gpl_eff:.2f} m/s "
          f"(ideale {d['V_gpl']:.2f})")
    print(f"J con queste velocita': {J_eff:.5f}  (di progetto {d['J']:.5f})")
    #: LA VERIFICA CHE RENDE LECITO IL CAMBIO. Passare dalla velocita' imposta
    #: alla portata imposta e' innocuo per il mescolamento SOLO SE i due flussi
    #: scalano dello stesso fattore: e' quello che tiene J, e quindi C, ai
    #: valori dimensionati. Se un giorno il sintetizzatore usasse due Cd diversi
    #: per aria e combustibile, la condizione cadrebbe in silenzio e il caso
    #: girerebbe lo stesso su un iniettore che non e' quello progettato. La
    #: tolleranza e' del 3 % perche' J entra nella penetrazione sotto radice:
    #: 3 % su J e' 1.5 % sulla penetrazione, sotto il rumore di tutto il resto.
    if abs(J_eff - d["J"]) > 0.03 * d["J"]:
        raise ValueError(
            f"J e' cambiato passando alla portata imposta: {J_eff:.5f} contro "
            f"{d['J']:.5f} ({(J_eff/d['J']-1)*100:+.1f} %). I due flussi non "
            "scalano dello stesso fattore, quindi la fisica del getto in flusso "
            "trasversale NON e' piu' quella dimensionata: il caso non va "
            "lanciato finche' non si capisce da dove viene la differenza.")

    def campo(nome, cls, dim, interno, bordi):
        righe = [intestazione(cls, nome, "0"), f"dimensions      {dim};",
                 f"internalField   uniform {interno};", "boundaryField", "{"]
        for b, spec in bordi.items():
            righe.append(f"    {b}\n    {{\n{spec}\n    }}")
        righe.append("}")
        (caso / "0" / nome).write_text("\n".join(righe) + "\n", encoding="utf-8")

    zg = "        type            zeroGradient;"
    sym = "        type            symmetryPlane;"
    campo("U", "volVectorField", "[0 1 -1 0 0 0 0]", "(0 0 0)", {
        #: `rho rho` dice al contorno di usare il CAMPO di densita' invece di
        #: un valore fisso: cosi' phi = mdot esattamente, qualunque cosa faccia
        #: la pressione. Con `rhoInlet` fisso si tornerebbe a un flusso di massa
        #: che dipende da p, cioe' al problema di prima con un altro nome.
        "ingresso_aria": (f"        type            flowRateInletVelocity;\n"
                          f"{_rampa(mdot_aria_settore, T_RAMPA)}\n"
                          f"        rho             rho;\n"
                          f"        value           uniform (0 0 0);"),
        "ingresso_gpl": (f"        type            flowRateInletVelocity;\n"
                         f"{_rampa(mdot_gpl_settore, T_RAMPA)}\n"
                         f"        rho             rho;\n"
                         f"        value           uniform (0 0 0);"),
        "uscita": "        type            pressureInletOutletVelocity;\n        value           uniform (0 0 0);",
        "simmetria_1": sym, "simmetria_2": sym, "pareti": "        type            noSlip;"})
    #: LO SBOCCO NON PUO' ESSERE A PRESSIONE IMPOSTA.
    #:
    #: Con `fixedValue` il contorno di uscita e' uno SPECCHIO ACUSTICO perfetto:
    #: l'onda di compressione che parte all'avvio - il dominio e' inizializzato
    #: fermo e gli ingressi cominciano a soffiare a 152 m/s - percorre i 18.6 mm
    #: del dominio, arriva allo sbocco e torna indietro tutta intera.
    #:
    #: E' successo due volte, allo STESSO ISTANTE. Il dominio e' lungo 18.6 mm e
    #: il suono in aria a 293 K fa 343 m/s: un transito sono 54 microsecondi. La
    #: prima corsa e' esplosa a t = 56.6 us, la seconda - con i limitatori
    #: numerici - a 56.2. Non era una divergenza casuale: era un'onda che
    #: tornava, sempre alla stessa ora.
    #:
    #: `waveTransmissive` risolveva l'equazione delle caratteristiche sul
    #: contorno e lasciava uscire l'onda invece di rimandarla dentro. `lInf`
    #: era la lunghezza su cui la pressione veniva rilassata verso `fieldInf`:
    #: piu' e' lunga, piu' il contorno e' trasparente e meno vincola la
    #: pressione media. Era 0.05 m, cioe' 2.7 volte il dominio.
    #:
    #: E QUELLA FRASE ERA GIA' LA DIAGNOSI DEL GUASTO SUCCESSIVO, scritta
    #: senza accorgersene: "meno vincola la pressione media". La corsa 5 e'
    #: morta esattamente li'. Con la PORTATA imposta su tutti e due gli
    #: ingressi, il dominio non ha nessun altro punto in cui la pressione sia
    #: fissata: l'unica ancora era il rilassamento debole dello sbocco. Fino a
    #: 300 us ha tenuto; poi lo squilibrio di massa ha vinto e la camera si e'
    #: svuotata senza ritorno -
    #:
    #:      finestra   camera   alimentazione   salto
    #:        100 us   6.53 bar     7.64 bar    1.11
    #:        300 us   6.13         6.53        0.39
    #:        500 us   4.65         9.79        5.14
    #:        700 us   3.38        13.40       10.02
    #:
    #: A camera bassa il rapporto di pressione sull'iniettore supera il
    #: critico, il getto STROZZA, e da li' in poi spingere una portata FISSA
    #: attraverso una gola strozzata a densita' calante richiede una pressione
    #: a monte che cresce senza limite. Il 9 % delle celle e' finito contro il
    #: limitatore a 500 K e k e' andata a 4.6e5 m2/s2. La corsa 4 non l'aveva
    #: mai visto solo perche' si era fermata a 250 us, cioe' PRIMA.
    #:
    #: LA CORREZIONE NON E' TARARE `lInf`. Il contorno trasparente serviva a
    #: non riflettere l'onda di compressione dell'avvio impulsivo (corse 1 e 2,
    #: esplose a 56.6 e 56.2 us, un transito acustico esatto). Ma quell'onda e'
    #: un problema di AVVIO, e si toglie alla radice: si RAMPA la portata
    #: d'ingresso da zero al valore di progetto in 20 us - vedi `_rampa` - e
    #: l'onda non nasce. Tolta l'onda, non serve piu' un contorno trasparente,
    #: e la pressione media si puo' finalmente ancorare.
    #:
    #: `fixedMean` impone che la MEDIA D'AREA sullo sbocco valga p_c, lasciando
    #: al profilo la liberta' di svilupparsi. E' anche la condizione onesta:
    #: p_c = 6.361 bar e' un dato del ciclo 0-D, non una previsione della CFD
    #: (il dominio finisce a 12 mm, molto prima della gola, e con la chimica
    #: spenta una gola nel dominio strozzerebbe a 2.3 bar). Fissarla dichiara
    #: quell'ipotesi invece di lasciarla galleggiare.
    campo("p", "volScalarField", "[1 -1 -2 0 0 0 0]", f"{d['p_c']:.1f}", {
        "ingresso_aria": zg, "ingresso_gpl": zg,
        "uscita": (f"        type            fixedMean;\n"
                   f"        meanValue       {d['p_c']:.1f};\n"
                   f"        value           uniform {d['p_c']:.1f};"),
        "simmetria_1": sym, "simmetria_2": sym, "pareti": zg})
    campo("T", "volScalarField", "[0 0 0 1 0 0 0]", f"{d['T_aria']:.1f}", {
        "ingresso_aria": f"        type            fixedValue;\n        value           uniform {d['T_aria']:.1f};",
        "ingresso_gpl": f"        type            fixedValue;\n        value           uniform {d['T_gpl']:.1f};",
        #: se un ricircolo fa RIENTRARE gas dallo sbocco, `zeroGradient` gli
        #: lascia portare dentro qualunque valore si trovi fuori. `inletOutlet`
        #: si comporta da zeroGradient in uscita e impone il valore di
        #: riferimento in entrata: e' la stessa condizione, ma con una risposta
        #: sensata anche nel caso che non dovrebbe capitare.
        "uscita": (f"        type            inletOutlet;\n"
                   f"        inletValue      uniform {d['T_aria']:.1f};\n"
                   f"        value           uniform {d['T_aria']:.1f};"),
        "simmetria_1": sym, "simmetria_2": sym, "pareti": zg})
    for sp in SPECIE:
        campo(sp, "volScalarField", "[0 0 0 0 0 0 0]", f"{Y_ARIA[sp]:.5f}", {
            "ingresso_aria": f"        type            fixedValue;\n        value           uniform {Y_ARIA[sp]:.5f};",
            "ingresso_gpl": f"        type            fixedValue;\n        value           uniform {1.0 if sp == 'C3H8' else 0.0:.5f};",
            "uscita": (f"        type            inletOutlet;\n"
                       f"        inletValue      uniform {Y_ARIA[sp]:.5f};\n"
                       f"        value           uniform {Y_ARIA[sp]:.5f};"),
            "simmetria_1": sym, "simmetria_2": sym, "pareti": zg})
    campo("k", "volScalarField", "[0 2 -2 0 0 0 0]", f"{k_aria:.5f}", {
        "ingresso_aria": f"        type            fixedValue;\n        value           uniform {k_aria:.5f};",
        "ingresso_gpl": f"        type            fixedValue;\n        value           uniform {k_gpl:.5f};",
        "uscita": zg, "simmetria_1": sym, "simmetria_2": sym,
        "pareti": "        type            kqRWallFunction;\n        value           uniform 1e-8;"})
    campo("omega", "volScalarField", "[0 0 -1 0 0 0 0]", f"{om_aria:.2f}", {
        "ingresso_aria": f"        type            fixedValue;\n        value           uniform {om_aria:.2f};",
        "ingresso_gpl": f"        type            fixedValue;\n        value           uniform {om_gpl:.2f};",
        "uscita": zg, "simmetria_1": sym, "simmetria_2": sym,
        "pareti": "        type            omegaWallFunction;\n        value           uniform 1e3;"})
    for nome, dim in (("alphat", "[1 -1 -1 0 0 0 0]"), ("nut", "[0 2 -1 0 0 0 0]")):
        wf = "compressible::alphatWallFunction" if nome == "alphat" else "nutkWallFunction"
        campo(nome, "volScalarField", dim, "0", {
            "ingresso_aria": "        type            calculated;\n        value           uniform 0;",
            "ingresso_gpl": "        type            calculated;\n        value           uniform 0;",
            "uscita": "        type            calculated;\n        value           uniform 0;",
            "simmetria_1": sym, "simmetria_2": sym,
            "pareti": f"        type            {wf};\n        value           uniform 0;"})

    (caso / "constant" / "thermophysicalProperties").write_text(
        intestazione("dictionary", "thermophysicalProperties", "constant") + f"""
thermoType
{{
    type            hePsiThermo;
    mixture         multiComponentMixture;
    transport       sutherland;
    thermo          janaf;
    energy          sensibleEnthalpy;
    equationOfState perfectGas;
    specie          specie;
}}

//: L'ordine conta: OpenFOAM ricava l'ULTIMA specie per differenza, quindi
//: l'inerte va per ultimo. Con C3H8 per ultimo, l'errore di chiusura finirebbe
//: tutto sulla specie che si sta cercando di misurare.
species ( C3H8 O2 N2 );
inertSpecie N2;
{_termo_da_cantera()}
""", encoding="utf-8")

    (caso / "constant" / "turbulenceProperties").write_text(
        intestazione("dictionary", "turbulenceProperties", "constant") + """
simulationType      RAS;
RAS
{
    RASModel        kOmegaSST;
    turbulence      on;
    printCoeffs     on;
}
""", encoding="utf-8")
    (caso / "constant" / "chemistryProperties").write_text(
        intestazione("dictionary", "chemistryProperties", "constant") + """
//: CHIMICA SPENTA. Questo caso misura il MESCOLAMENTO, non la combustione:
//: accendendola, un risultato brutto non direbbe quale dei due modelli e'
//: sbagliato.
chemistryType { solver noChemistrySolver; }
chemistry           off;
initialChemicalTimeStep 1e-07;
""", encoding="utf-8")
    (caso / "constant" / "fvOptions").write_text(
        intestazione("dictionary", "fvOptions", "constant") + """
//: LIMITI NUMERICI DICHIARATI, non nascosti.
//:
//: In questo caso la chimica e' spenta e i due ingressi stanno a 283 e 293 K:
//: con 152 m/s la differenza fra temperatura totale e statica vale una decina
//: di gradi, quindi la temperatura FISICA non puo' uscire da 270-310 K. La
//: prima corsa e' arrivata a 203 e 334 K, e poi una singola cella e' esplosa:
//: il passo temporale e' crollato da 1.5e-8 a 2.7e-12 secondi con il Courant
//: MEDIO a 3e-6 e il MASSIMO a 0.4, cioe' una cella su un milione decideva
//: tutto. Nove ore di calcolo per arrivare a meta' del transitorio.
//:
//: Questi limiti non "aggiustano" il risultato: tagliano valori che non
//: esistono. Se pero' intervengono su molte celle il risultato NON e' valido, e
//: per questo il solutore stampa quante ne tocca a ogni passo: e' un numero da
//: guardare, non da ignorare.
limitaT
{
    type            limitTemperature;
    active          yes;
    selectionMode   all;
    min             250;
    max             500;
}
//: IL LIMITATORE DI VELOCITA' E' STATO TOLTO, e vale la pena dire perche'.
//: `limitVelocity` taglia U dopo la soluzione della quantita' di moto, ma NON
//: corregge il flusso phi ne' la pressione: la massa che entra in una cella
//: smette di essere quella che ne esce. In una corsa che sta gia' divergendo
//: questo non contiene il problema, lo sposta - e infatti nella corsa 3 la
//: velocita' massima restava inchiodata a 400.000 m/s (cioe' il limitatore
//: lavorava a ogni passo) mentre la pressione di una cella scendeva a 1e-4 bar
//: indisturbata. Mascherava il sintomo e lasciava correre la causa. Il posto
//: giusto dove fermare una divergenza e' l'equazione della pressione: `pMin`
//: nel PIMPLE, che e' un limite sulla variabile che sta davvero divergendo.
""", encoding="utf-8")
    (caso / "constant" / "combustionProperties").write_text(
        intestazione("dictionary", "combustionProperties", "constant") + """
combustionModel     none;
""", encoding="utf-8")

    #: PIANI DI TAGLIO. Il settore va da 0 a 2*ang in y-z, con il getto sulla
    #: bisettrice. Servono due piani meridiani diversi e non intercambiabili:
    #:   * quello che CONTIENE l'asse del getto (normale perpendicolare al
    #:     getto), dove si vede il getto penetrare;
    #:   * il fianco a meta' strada fra due getti (z = 0), dove si vede se il
    #:     GPL arriva anche dove non e' stato sparato. E' li' che una miscela
    #:     disuniforme si tradisce.
    n_getto = (0.0, -math.sin(ang), math.cos(ang))
    piani = "\n".join(
        f"""        sez_{i:02d}
        {{
            type            cuttingPlane;
            planeType       pointAndNormal;
            pointAndNormalDict {{ point ({x:.6f} 0 0); normal (1 0 0); }}
            interpolate     true;
        }}""" for i, x in enumerate(d["piani_x"]))

    #: NOTA DI METODO. Il campionamento delle superfici e' separato dalla
    #: scrittura dei campi di volume: le superfici costano poco e vanno scritte
    #: fitte (e' il film), i volumi costano molto e bastano poche istantanee.
    #: Scriverli con lo stesso passo significa scegliere fra un film a scatti e
    #: un disco pieno.
    passo_sup = campionamento if campionamento else scrittura
    contorni = _controlli_di_contorno(passo_sonde)
    #: se non lo si dice, si mediano gli ultimi due terzi della corsa: il
    #: primo terzo e' assestamento in ogni caso realistico, e mediarlo
    #: sporcherebbe la media con un transitorio che non e' il motore.
    t_medie = media_da if media_da is not None else tempo_finale / 3.0
    turbolenza = _controlli_di_turbolenza(t_medie, passo_sup)
    #: quanto spostare il piano di fianco: una frazione dell'altezza del
    #: condotto, abbastanza da uscire dal bordo e abbastanza poco da restare
    #: a meta' strada fra due getti. Si dichiara in gradi perche' e' cosi'
    #: che va giudicato.
    scarto_fianco = 0.12 * d["H"]
    ang_fianco = math.degrees(math.asin(min(scarto_fianco / d["R_getti"], 1.0)))
    ang_getto = 180.0 / d["n_getti"]
    if ang_fianco > 0.25 * ang_getto:
        raise ValueError(
            f"il piano di fianco e' spostato di {ang_fianco:.1f} gradi contro i "
            f"{ang_getto:.0f} che lo separano dal getto: non e' piu' 'meta' strada'.")
    formato = "binary" if binario else "ascii"
    (caso / "system" / "controlDict").write_text(
        intestazione("dictionary", "controlDict", "system") + f"""
application         reactingFoam;
startFrom           latestTime;
startTime           0;
stopAt              endTime;
endTime             {tempo_finale:.8f};
deltaT              1e-08;
writeControl        adjustableRunTime;
writeInterval       {scrittura:.8f};
//: SI TENGONO GLI ULTIMI ISTANTI DI VOLUME, NON TUTTI. Non e' avarizia di
//: disco: `startFrom latestTime` riprende dall'ultimo scritto, quindi la
//: frequenza di scrittura E' la granularita' con cui una corsa interrotta si
//: puo' riprendere. La prima corsa e' stata interrotta due volte da cause
//: esterne al solutore (la macchina virtuale si e' spenta sotto di lui) e ogni
//: volta ha perso tutto il tratto dall'ultima scrittura: con un intervallo di
//: 40 microsecondi erano quaranta minuti buttati. I piani campionati, che sono
//: il risultato vero, non vengono mai cancellati.
//: DODICI istanti, non quattro. Con quattro e una scrittura ogni 20 us il
//: disco conserva 80 us di storia: la corsa 5 e' degenerata a partire da
//: 387 us ed e' stata trovata a 634, quando gli unici istanti rimasti
//: (560-620) erano gia' tutti contaminati. Non c'era piu' nessuno stato sano
//: da cui ripartire, e sei ore sono state buttate invece che recuperate.
//: Dodici istanti sono 240 us di storia: piu' del tempo che passa fra il
//: primo sintomo e il momento in cui uno se ne accorge.
purgeWrite          12;
writeFormat         {formato};
writePrecision      8;
runTimeModifiable   true;
adjustTimeStep      yes;
maxCo               0.4;
maxDeltaT           1e-06;

functions
{{
//: I CONTORNI SI MISURANO, NON SI DANNO PER BUONI.
//:
//: Domanda a cui questi tre blocchi rispondono con un numero invece che con un
//: ragionamento: "la simulazione gira con le pressioni giuste, da ogni lato?"
//: Prima non c'era risposta, perche' sugli ingressi NON si impone nessuna
//: pressione (si impone la portata) e quella che si sviluppa e' un risultato -
//: un risultato che nessuno guardava. Il dimensionamento dice: alimentazione
//: 7.315 bar, camera 6.361, salto d'iniezione 0.954 bar = 15.0 % di p_c. Se la
//: CFD sviluppa un salto diverso, o il modello 0-D dell'iniettore e' ottimista
//: o la geometria non e' quella che il modello crede. In entrambi i casi si
//: vuole saperlo mentre gira, non a fine corsa.
//:
//: `sum(phi)` sugli ingressi e' la verifica della verifica: deve valere
//: esattamente la portata imposta. Se `flowRateInletVelocity` non stesse
//: facendo il suo mestiere - per esempio se `rho` non fosse il campo giusto -
//: si vedrebbe qui e da nessun'altra parte.
{contorni}
{turbolenza}
    superfici
    {{
        type            surfaces;
        libs            (sampling);
        writeControl    adjustableRunTime;
        writeInterval   {passo_sup:.8f};
        surfaceFormat   vtk;
        writeFormat     binary;
        interpolationScheme cellPoint;
        //: `vorticity`, `k` e `nut` viaggiano con gli altri perche' le
        //: domande sul mescolamento e quelle sulla turbolenza si guardano
        //: sullo STESSO istante e sullo stesso piano. Campionarle in un
        //: secondo momento vorrebbe dire riaprire i volumi, che sono scritti
        //: rado apposta, e confrontare istanti diversi.
        fields          (C3H8 O2 T U rho vorticity k nut);
        surfaces
        {{
            meridiano_getto
            {{
                type            cuttingPlane;
                planeType       pointAndNormal;
                pointAndNormalDict
                {{ point (0 0 0); normal ({n_getto[0]:.6f} {n_getto[1]:.6f} {n_getto[2]:.6f}); }}
                interpolate     true;
            }}
            //: IL PIANO DI FIANCO E' SPOSTATO DENTRO IL DOMINIO, e non e' un
            //: dettaglio: il settore ha i due fianchi dichiarati piani di
            //: SIMMETRIA, e uno di essi E' il piano z = 0. Un `cuttingPlane`
            //: che ci cade sopra non ha celle da intersecare e restituisce una
            //: manciata di triangoli degeneri - il piano usciva VUOTO in tre
            //: corse su tre, e la figura mostrava mezza tavola bianca senza
            //: che niente segnalasse un errore.
            //: Lo scarto vale {scarto_fianco*1e3:.1f} mm, cioe' {ang_fianco:.1f} gradi al raggio
            //: d'iniezione, contro i {ang_getto:.0f} che separano questo piano dal getto:
            //: resta a tutti gli effetti "meta' strada fra due getti".
            meridiano_fianco
            {{
                type            cuttingPlane;
                planeType       pointAndNormal;
                pointAndNormalDict {{ point (0 0 {scarto_fianco:.6e}); normal (0 0 1); }}
                interpolate     true;
            }}
{piani}
        }}
    }}
}}
""", encoding="utf-8")

    (caso / "system" / "fvSchemes").write_text(
        intestazione("dictionary", "fvSchemes", "system") + """
ddtSchemes      { default Euler; }
gradSchemes     { default cellLimited Gauss linear 1; }
divSchemes
{
    default         none;
    div(phi,U)      Gauss limitedLinearV 1;
    //: limitedLinear sulle frazioni di massa e non upwind puro: l'upwind
    //: diffonde numericamente, e una diffusione numerica in un caso di
    //: MESCOLAMENTO produce esattamente il risultato che si sta cercando di
    //: misurare. Sarebbe la peggiore delle conferme.
    div(phi,Yi_h)   Gauss limitedLinear01 1;
    div(phi,K)      Gauss limitedLinear 1;
    div(phi,k)      Gauss limitedLinear 1;
    div(phi,omega)  Gauss limitedLinear 1;
    div(phid,p)     Gauss limitedLinear 1;
    div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear corrected; }
interpolationSchemes { default linear; }
snGradSchemes   { default corrected; }

//: k-omega SST ha bisogno della distanza dalla parete, e OpenFOAM non la
//: calcola se non gli si dice come. meshWave e' esatto su mesh non strutturate.
wallDist        { method meshWave; }
""", encoding="utf-8")

    (caso / "system" / "fvSolution").write_text(
        intestazione("dictionary", "fvSolution", "system") + """
solvers
{
    "rho.*"         { solver diagonal; }
    p
    {
        solver          PCG;
        preconditioner  DIC;
        tolerance       1e-07;
        relTol          0.01;
    }
    pFinal          { $p; relTol 0; }
    "(U|h|k|omega|Yi)"
    {
        solver          PBiCGStab;
        preconditioner  DILU;
        tolerance       1e-08;
        relTol          0.01;
    }
    "(U|h|k|omega|Yi)Final" { $U; relTol 0; }
}

PIMPLE
{
    momentumPredictor   yes;
    nOuterCorrectors    1;
    nCorrectors         2;
    //: DUE correttori non-ortogonali, non uno. La mesh e' tetraedrica con
    //: non-ortogonalita' massima 57.9 gradi (checkMesh): e' esattamente il
    //: caso per cui il correttore esiste, e costa solo una soluzione di p in
    //: piu' per passo.
    nNonOrthogonalCorrectors 2;
    //: LIMITI SULLA PRESSIONE, dichiarati. Il campo lavora fra 6.36 bar
    //: (sbocco) e 7.32 bar (alimentazione): 1 bar non e' raggiungibile da
    //: nessun processo fisico dentro questo dominio, che e' tutto subsonico e
    //: senza espansioni. Una cella sotto pMin non e' fisica, e' il solutore
    //: che sta divergendo - ed e' quello che e' successo tre volte. Il limite
    //: non aggiusta un risultato: impedisce a rho = psi*p di andare a zero e
    //: portarsi dietro il passo temporale di tutti.
    //: SCALARI NUDI, senza nome e senza dimensioni: `pressureControl` li
    //: legge con `readScalar`, e la forma `pMin pMin [1 -1 -2 ...] 1e5` fa
    //: abortire tutti e 24 i processi al primo passo con "Wrong token type".
    pMin                1.0e5;
    pMax                5.0e6;
}
""", encoding="utf-8")

    (caso / "system" / "decomposeParDict").write_text(
        intestazione("dictionary", "decomposeParDict", "system") + f"""
numberOfSubdomains  {n_proc};
method              scotch;
""", encoding="utf-8")

    #: NIENTE BLOCCO `functions`, ed e' un aggiramento, non una scelta.
    #: Questa build di OpenFOAM 1912 (pacchetto Debian) cade su QUALUNQUE
    #: functionObject con un errore su uno stream chiamato "sha1", prima
    #: ancora del primo passo temporale. Provate e scartate: sintassi a
    #: dizionario e a lista per `surfaces`, `libs` in forma corta e lunga,
    #: runTimeModifiable a on e off. Cade sempre.
    #:
    #: Il campionamento si fa quindi DOPO, sui campi scritti, con foamToVTK e
    #: `scripts/cfd/analizza_mescolamento.py`. Costa disco e un passaggio in
    #: piu', ma ha un vantaggio: le misure si possono rifare e correggere senza
    #: rilanciare il calcolo.


def main() -> int:
    from mesh_iniettore import costruisci, quote_dal_progetto

    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", type=Path, default=Path("runs/cfd/mescolamento"))
    ap.add_argument("--fine", type=float, default=6.0e-5)
    ap.add_argument("--grossa", type=float, default=4.0e-4)
    ap.add_argument("--tempo", type=float, default=2.0e-3, help="s simulati")
    ap.add_argument("--scrittura", type=float, default=1.0e-4)
    ap.add_argument("--proc", type=int, default=8)
    ap.add_argument("--solo-caso", action="store_true", help="non rifare la mesh")
    ap.add_argument("--campionamento", type=float, default=None,
                    help="passo di scrittura dei piani di taglio [s]")
    ap.add_argument("--binario", action="store_true",
                    help="campi di volume in binario (mesh grandi)")
    ap.add_argument("--passo-sonde", type=int, default=20,
                    help="ogni quanti passi si misurano pressioni, portate e "
                         "massa contenuta (e' la curva di convergenza)")
    ap.add_argument("--media-da", type=float, default=None,
                    help="istante [s] da cui accumulare le medie temporali; "
                         "va messo DOPO l'assestamento, non prima")
    a = ap.parse_args()

    d = quote_dal_progetto()
    print(f"getti {d['n_getti']} x {d['d_getto']*1e3:.3f} mm   condotto H "
          f"{d['H']*1e3:.3f} mm   C {d['C']:.3f}   J {d['J']:.3f}")
    print(f"aria {d['V_aria']:.1f} m/s   GPL {d['V_gpl']:.1f} m/s   "
          f"p_c {d['p_c']/1e5:.3f} bar")
    print(f"lunghezza di mescolamento promessa dalla correlazione: "
          f"{d['L_mescolamento']*1e3:.2f} mm")

    if a.caso.exists() and not a.solo_caso:
        shutil.rmtree(a.caso)
    msh = a.caso.parent / "iniettore.msh"
    if not a.solo_caso:
        n = costruisci(d, msh, a.fine, a.grossa)
        print(f"mesh: {n} tetraedri")
    scrivi(a.caso, d, a.tempo, a.scrittura, a.proc, a.campionamento, a.binario,
           a.passo_sonde, a.media_da)
    (a.caso / "quote.txt").write_text(
        "\n".join(f"{k} = {v}" for k, v in sorted(d.items())) + "\n", encoding="utf-8")
    print(f"caso scritto in {a.caso}")
    print(f"\nda lanciare (dentro l'ambiente OpenFOAM):")
    print(f"  gmshToFoam {msh} -case {a.caso}")
    print(f"  # i gruppi di gmsh arrivano come 'patch': i due fianchi vanno")
    print(f"  # dichiarati symmetryPlane, altrimenti il settore non e' un settore")
    print(f"  decomposePar -case {a.caso}")
    print(f"  mpirun -np {a.proc} reactingFoam -parallel -case {a.caso}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
