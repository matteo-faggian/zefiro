"""Il motore da 50 N come geometria implicita, con raffreddamento conforme.

L'architettura del raffreddamento non e' una scelta di stile: viene dalla
mappa del calore (`scripts/heat_map_50N.py`) e dai vincoli idraulici.

  * il calore NON e' concentrato in gola. La camera prende il 62 % dei 3.69 kW
    totali, perche' ha l'area; il plug ne prende il 10 %. Quindi i canali
    servono soprattutto in camera, non solo alla gola;
  * un circuito UNICO non funziona. Serve h alto solo al labbro, ma la
    velocita' che lo produce va pagata in perdita di carico su TUTTA la
    lunghezza: 2.5 bar su una rete che ne da' 4, con 14 K di margine
    all'ebollizione. Due circuiti in PARALLELO da un collettore comune
    scendono a 0.77 bar, perche' il ramo veloce e' anche quello corto.

E' questo che produce i collettori: non sono decorazione, sono la conseguenza
di due zone con richieste opposte alimentate dalla stessa acqua.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field as dc_field
from typing import Sequence

import numpy as np

from zefiro.sdf.core import Field, Grid
from zefiro.sdf.printability import passo_elica_minimo
from zefiro.sdf.shapes import cylinder, helical_channels, revolve_polygon


#: Larghezza assiale della gola anulare di collettore [m]. Presa pari a circa
#: due volte il lato del canale piu' grande: serve che la sezione di passaggio
#: del collettore sia molto maggiore di quella di un singolo canale, altrimenti
#: e' il collettore a strozzare e la distribuzione fra i canali non e' uniforme.
COLLETTORE_LARGHEZZA = 3.0e-3

#: Pendenza del soffitto del collettore, in avanzamento assiale per unita' di
#: profondita'. 1.0 darebbe esattamente 45 gradi, cioe' il limite: con la
#: discretizzazione a voxel una parte delle faccette finiva appena sotto.
#: 2.4 da' circa 67 gradi sul tratto CILINDRICO. Serve tutto quel margine
#: perche' lo smusso e' definito rispetto alla NORMALE alla parete, e sul
#: convergente la normale e' gia' inclinata di 37 gradi rispetto all'asse:
#: quel che sul cilindro sono 67 gradi, sul cono ne diventa una trentina.
#: Con 1.4 (54 gradi sul cilindro) i collettori sul cono finivano sotto
#: soglia, e sono tre su sei.
COLLETTORE_PENDENZA = 3.2

#: Angolo di autosostentamento richiesto in progetto, piu' severo della soglia
#: fisica di 45 gradi. Stessa logica dello smusso: si progetta con margine
#: sulla soglia di processo, non uguale ad essa.
MARGINE_ANGOLO_DEG = 52.0

#: Diametro dell'attacco radiale [m] e quanti attacchi per collettore.
#:
#: TRE attacchi da 3.2 mm invece di uno da 5.6. Stessa sezione complessiva
#: (24 contro 25 mm2), ma tre vantaggi che uno solo non da':
#:   * un foro orizzontale di 5.6 mm ha una calotta superiore che in SLM cede;
#:     3.2 mm e' una luce che il processo attraversa senza supporti;
#:   * su un motore lungo 47 mm due attacchi da 5.6 mm posti a 4 mm di distanza
#:     SI COMPENETRANO, e i due rami paralleli del raffreddamento diventano uno
#:     solo. Con attacchi piccoli e sfalsati angolarmente il problema sparisce;
#:   * tre ingressi a 120 gradi distribuiscono l'acqua nel collettore meglio di
#:     uno solo, che alimenterebbe di piu' i canali che gli stanno davanti.
PORTA_DIAMETRO = 3.2e-3
PORTE_PER_COLLETTORE = 3

#: Sfasamento angolare fra collettori consecutivi, in frazione del passo fra
#: attacchi dello stesso collettore. Serve a garantire che due collettori
#: vicini in x non abbiano MAI attacchi allineati.
SFALSAMENTO_COLLETTORI = 0.5

#: Quanto il collettore e' piu' profondo, per lato, dei canali che unisce [m].
COLLETTORE_MARGINE = 1.5e-4


@dataclass(frozen=True)
class CircuitoRaffreddamento:
    """Un ramo del raffreddamento. Tutti i numeri vengono dal dimensionamento
    idraulico e termico, nessuno e' scelto per estetica."""
    nome: str
    n_canali: int
    lato: float            # m, sezione quadrata
    parete_calda: float    # m, spessore fra gas e canale
    parete_fredda: float   # m, spessore fra canale ed esterno
    x_inizio: float        # m
    x_fine: float          # m
    passo_elica: float     # m, avanzamento assiale per giro completo
    velocita: float        # m/s (per la documentazione: non entra in geometria)
    #: Circuito a U: l'acqua scende fino a `x_fine`, gira nel collettore di
    #: fondo e torna indietro su un SECONDO strato di canali piu' esterno,
    #: cosi' entrambi gli attacchi stanno a `x_inizio`.
    #:
    #: Non e' un vezzo. Verso il labbro la parete e' conica, e nell'ingombro
    #: assiale di un attacco radiale il raggio cambia di 2 mm: non esiste
    #: nessuna profondita' che insieme raggiunga il collettore e non buchi la
    #: parete. Misurato: l'attacco arrivava a 0.18 mm dal gas e in parte lo
    #: bucava. Con la U il problema sparisce, perche' gli attacchi tornano
    #: sulla parte cilindrica.
    ritorno: bool = False
    #: Pieno fra lo strato di andata e quello di ritorno. 0.6 mm e non 0.4:
    #: 0.4 sta sotto il minimo di parete dichiarato per l'SLM, e una parete
    #: sotto il minimo non e' sottile, e' porosa - cioe' i due strati
    #: comunicano e il circuito si cortocircuita da solo.
    setto_ritorno: float = 0.6e-3
    #: Di quanto lo strato di RITORNO parte piu' a valle di quello di andata.
    #: Non e' un dettaglio: l'attacco che alimenta lo strato interno (andata)
    #: arriva da fuori e dovrebbe attraversare quello esterno (ritorno),
    #: mettendoli in comunicazione. Facendo partire il ritorno piu' a valle si
    #: lascia un tratto in cui lo strato esterno NON c'e', e li' si mette
    #: l'attacco dell'andata.
    sfalsamento_ritorno: float = 3.0e-3



#: Diametro interno della CANNA dell'acqua [m].
#:
#: DATO NON MISURATO. 1/2 pollice (12.7 mm) e' il formato da giardinaggio piu'
#: diffuso, ma la canna non e' stata misurata: questo numero e' una SCELTA
#: dichiarata, non una derivazione, e se la canna e' da 5/8 o 3/4 cambia il
#: portagomma e nient'altro. Sta qui e non sparso nel codice proprio perche'
#: cambiarlo sia un'operazione sola.
CANNA_DIAMETRO_INTERNO = 12.7e-3

#: Quanto la cresta del portagomma e' piu' larga del diametro interno della
#: canna [m]. E' l'interferenza che, con una fascetta, fa la tenuta.
PORTAGOMMA_INTERFERENZA = 0.6e-3
PORTAGOMMA_CRESTE = 3
PORTAGOMMA_PASSO = 4.0e-3
PORTAGOMMA_LUNGHEZZA = 16.0e-3


@dataclass(frozen=True)
class Portagomma:
    """Attacco per la canna dell'acqua.

    ASSIALE, non radiale, e non e' una preferenza estetica. In stampa si
    costruisce lungo +x: un cilindro con l'asse lungo x non ha nessuna
    superficie rivolta in basso, mentre un attacco radiale e' un cilindro
    ORIZZONTALE, e sopra i 3 mm circa la sua calotta superiore cede. Un
    portagomma da mezzo pollice radiale non si stampa senza supporti, e dentro
    un condotto i supporti non si tolgono.
    """
    x_base: float          # m, dove il gambo poggia sul pezzo
    #: +1 sporge verso +x, -1 verso -x, 0 = RADIALE (sporge verso l'esterno).
    #: Un portagomma radiale ha il foro ORIZZONTALE rispetto alla costruzione:
    #: si stampa con il tetto a due falde ("a goccia") e non tondo, altrimenti
    #: la calotta superiore e' uno sbalzo e cede.
    verso: float
    R_posizione: float     # m, raggio dell'asse
    angolo: float          # rad
    d_canna: float         # m, diametro INTERNO della canna
    d_foro: float          # m, foro di passaggio
    lunghezza: float       # m, sporgenza del gambo
    n_creste: int
    passo_creste: float
    interferenza: float    # m, quanto la cresta e' piu' larga della canna


@dataclass(frozen=True)
class FessuraFilm:
    """La fessura che protegge il corpo centrale, e come ci arriva il gas.

    IL PEZZO DI MOTORE CHE NESSUN CIRCUITO PUO' RAGGIUNGERE. Il plug e' un cono
    lungo sette millimetri, largo quattro alla base e appeso nel getto: non c'e'
    posto per un canale d'acqua, non c'e' parete a cui appoggiarlo, e senza film
    fonde a 2.8 s contro i 5 richiesti. La fessura non e' un miglioramento, e'
    l'unica cosa che tiene insieme il progetto.

    TRE PEZZI, E OGNUNO RISOLVE UN PROBLEMA DIVERSO:

      * il **foro di dosatura**, assiale, subito a valle del plenum del GPL. E'
        lui a decidere quanto combustibile va al film, e non la fessura: la
        fessura vede lo stesso salto di pressione dei getti principali, quindi
        se fosse lei a dosare si porterebbe via una frazione confrontabile della
        portata invece del 18 %. Il suo diametro e' il piu' piccolo che il
        processo sa fare, e la frazione di film e' la CONSEGUENZA di quel
        diametro, non un desiderio;
      * il **condotto**, che porta il gas dalla testa fino alla base del plug -
        quaranta millimetri dentro un corpo centrale largo otto. E' assiale,
        quindi in stampa non ha soffitto;
      * la **fessura** vera, un anello inclinato di `angolo` sulla direzione di
        costruzione che sbocca raso sulla superficie del plug. Inclinata e non
        radiale per due ragioni che coincidono: un film deve uscire TANGENTE
        alla parete che protegge, e un condotto inclinato entro l'angolo di
        processo e' autosostentato.

    IL LABBRO NON PUO' ESSERE UNA LAMA. All'uscita la fessura lascia sopra di
    se' un anello di materiale che si assottiglia fino a zero: si tronca a
    `spessore_labbro`, perche' un labbro piu' sottile della parete minima non e'
    sottile, e' assente.
    """
    x_uscita: float          # m, dove la fessura sbocca
    r_uscita: float          # m, raggio della superficie del plug li'
    altezza: float           # m, apertura della fessura misurata normale
    angolo: float            # rad, inclinazione sull'asse del motore
    r_condotto: float        # m, raggio del condotto assiale
    x_condotto_0: float      # m, dove il condotto comincia (a valle del plenum)
    d_dosatura: float        # m, foro che decide la portata
    x_dosatura: float        # m
    spessore_labbro: float   # m
    avvertenze: tuple[str, ...] = ()


@dataclass(frozen=True)
class InterfacciaMontaggio:
    """Come il motore si attacca a qualcosa.

    Fino a questa revisione non si attaccava a niente: era un pezzo bellissimo
    da tenere in mano. Un motore che non si puo' vincolare non si puo' provare,
    e non si puo' nemmeno accendere in sicurezza.

    PERCHE' SULLA FACCIA DI MONTE, E NON DI FIANCO. Tre ragioni, in ordine di
    peso:

      1. **Stampa.** Si costruisce lungo +x con la faccia d'iniezione sulla
         piastra: un rilievo su quella faccia e' la PRIMA cosa che si costruisce
         e non ha, per definizione, nessuno sbalzo. Una flangia radiale a meta'
         camera avrebbe invece una corona rivolta verso il basso larga quanto la
         flangia stessa, e per smussarla a 52 gradi servirebbero 12 mm di corsa
         assiale che fra i due collari dell'acqua non ci sono (ce ne sono 9).
      2. **Carico.** La spinta e' diretta verso monte: appoggiando la faccia a
         una piastra il vincolo lavora in COMPRESSIONE, che e' il modo in cui un
         'interfaccia si rompe meno.
      3. **Temperatura.** E' l'estremita' fredda. La parete della camera sta a
         378 K perche' e' raffreddata; la testa sta all'ambiente.

    IL PREZZO, ed e' un prezzo vero: il raccordo del GPL e' assiale al centro
    della stessa faccia. La piastra di banco deve avere un foro centrale da cui
    esce la linea del combustibile. E' la disposizione classica di un motore
    bullonato a una piastra di spinta, ma va detta perche' vincola il banco.

    QUATTRO PIAZZOLE E NON UNA FLANGIA ANULARE: la corona continua dovrebbe
    passare sopra il bocchettone dell'aria, che sporge a 18.4 mm di raggio e ha
    bisogno di spazio per una chiave. Quattro piazzole a 45, 135, 225 e 315
    gradi lasciano libero l'asse dell'aria (0 gradi) e quello delle canne
    dell'acqua (180 gradi).
    """
    #: FORI PASSANTI E DADO DIETRO, non filetti ricavati nel pezzo.
    #:
    #: Un filetto stampato non tiene il passo, e uno maschiato in un foro cieco
    #: di acciaio sinterizzato si sfoglia al primo serraggio - e se il maschio si
    #: spezza dentro, il pezzo e' perso. Con i fori passanti la filettatura sta
    #: nel dado, che e' un pezzo di commercio e costa dieci centesimi. Il prezzo
    #: e' che dietro la flangia deve esserci spazio per la chiave: e' il vincolo
    #: che ha spinto il cerchio dei bulloni fuori dal corpo della testa.
    n_bracci: int
    angolo_offset: float     # rad, sfasamento del primo braccio
    R_bulloni: float         # m, raggio del cerchio dei bulloni
    d_foro: float            # m, foro passante di passaggio
    d_piazzola: float        # m, diametro della piazzola attorno al foro
    larghezza_braccio: float # m
    spessore: float          # m, quanto flangia e bracci sporgono dalla faccia
    r_attacco: float         # m, da dove il braccio nasce sul corpo della testa
    d_filetto: float         # m, diametro nominale del bullone
    ingombro_chiave: float   # m, spazio libero che serve dietro per il dado
    carico_per_bullone: float   # N, quello che deve reggere
    avvertenze: tuple[str, ...] = ()


@dataclass(frozen=True)
class CollettoreAcqua:
    """Il collettore che rende il raffreddamento una SERIE con un solo ingresso
    e una sola uscita.

    PERCHE' ESISTE. L'impianto ha una canna sola. Due rami in parallelo vogliono
    due ingressi e due uscite, e in piu' una strozzatura di bilanciamento: due
    rami in parallelo hanno per definizione la stessa caduta, quindi con cadute
    di progetto diverse l'acqua se ne va quasi tutta nel ramo che oppone meno
    resistenza - via dalla gola, che e' il ramo che ne ha bisogno. Quel foro
    calibrato e' un pezzo che DEVE essere giusto e che, se sbagliato, non si
    vede: portata e temperatura all'uscita restano normali mentre il labbro non
    e' piu' raffreddato. In serie non c'e' niente da bilanciare.

    IL PERCORSO:

        canna -> anello d'ingresso -> camera -> raccordo -> gola
              -> inversione al labbro -> ritorno -> anello d'uscita -> canna

    Ogni elemento e' adiacente al successivo, e questo NON e' una comodita': le
    stazioni assiali dei collettori si intrecciano, e qualunque condotto che
    unisca due stazioni non adiacenti scavalca il collettore di un altro tratto.

    PERCHE' DUE ANELLI ESTERNI E NON DUE FORI. Un collettore anulare alimentato
    in un solo punto deve portare meta' della portata verso ciascun lato: la
    velocita' li' dentro diventa maggiore che nei canali, e la caduta lungo
    l'anello - diversa per il canale vicino all'attacco e per quello opposto -
    li alimenta in modo diverso. La cura non e' approfondire il collettore
    (costerebbe spessore di mantello su tutta la lunghezza, cioe' massa): e'
    mettere la distribuzione FUORI, in un anello che comunica con il collettore
    in piu' punti.

    PERCHE' NON NELLA TESTA. La prima stesura metteva l'ingresso dentro la
    piastra: un plenum anulare, fori assiali, portagomma sulla faccia. Due
    difetti, entrambi trovati dalla prova di tenuta e non a occhio.
    Primo: il raccordo dell'aria e' un foro RADIALE da 14 mm che attraversa la
    piastra, e un plenum anulare completo ci passa dentro. Secondo, e fatale:
    la testa NON e' un disco del diametro del mantello - e' un cilindro da 10.5
    mm di raggio che si allarga a cono solo nell'ultimo tratto. Il plenum, messo
    al raggio del collettore della camera, stava semplicemente in aria.
    """
    #: ingresso e uscita hanno la STESSA forma: anello esterno, fori radiali,
    #: portagomma radiale. Non e' pigrizia: e' il fatto che il vincolo (arrivare
    #: a un collettore anulare da fuori senza sbalzi) e' lo stesso ai due capi.
    ingresso: Portagomma
    #: (x0, x1, altezza radiale, x d'inizio del collare)
    anello_ingresso: tuple[float, float, float, float]
    x_collettore_ingresso: float
    n_fori_ingresso: int
    d_foro_ingresso: float
    #: anello di raccordo: vano anulare INTERNO al mantello che unisce il
    #: collettore d'uscita della camera a quello d'ingresso della gola su tutta
    #: la circonferenza. Non ha fori: i due collettori condividono un anello.
    #: (x0, x1, profondita' minima, profondita' massima) dalla parete del gas.
    raccordo: tuple[float, float, float, float]
    uscita: Portagomma
    anello_uscita: tuple[float, float, float, float]
    x_collettore_uscita: float
    n_fori_uscita: int
    d_foro_uscita: float
    #: quanto il tetto (e, nel collare, il pavimento) di un vano anulare arretra
    #: in x per ogni unita' di estensione radiale: e' tan(angolo di progetto)
    pendenza_tetto: float
    parete: float
    dp_totale: float
    avvertenze: tuple[str, ...] = ()


@dataclass(frozen=True)
class TestaIniezione:
    """La testa d'iniezione a getti trasversali.

    Sostituisce il coassiale a taglio, che per una coppia gas-gas e' l'elemento
    con l'efficienza di mescolamento piu' bassa fra quelli provati (NASA
    CR-121234). Il dimensionamento sta in `zefiro.injector`; qui c'e' solo la
    geometria che ne consegue, e le quote assiali che la letteratura non da'.

    Percorso dell'aria, da monte a valle:
      plenum anulare -> contrazione -> condotto a sezione costante -> camera.
    Il plenum serve a distribuire in circonferenza l'aria che arriva da UN solo
    attacco radiale; la contrazione la accelera da qualche m/s ai 152 m/s che
    il salto di pressione consente, e lo fa con una parete continua invece che
    con uno spigolo, perche' uno spigolo vivo da' una vena contratta e un Cd di
    0.6 invece di 0.95, cioe' il 35 % di area in piu' da trovare.

    Percorso del GPL:
      foro assiale sull'asse -> plenum dentro il corpo centrale -> `n_getti`
    fori radiali che attraversano il condotto d'aria.

    QUOTE ASSIALI, e perche' quelle:
      * i getti stanno a `x_getti` < 0, cioe' DENTRO il condotto e non sulla
        faccia. La correlazione di Holdeman vale per getti in un condotto
        CONFINATO: se i getti sparassero direttamente in camera non ci sarebbe
        nessun condotto e nessuna penetrazione da governare. Fra i getti e lo
        sbocco deve restare almeno la lunghezza di mescolamento (circa 2 H).
      * il corpo centrale e' STRIZZATO a `R_getti` nel condotto e torna al suo
        raggio di camera entro `x_rampa`. La strizione non e' estetica: a
        portate e salti fissati vale d_getto = 0.386 H e H = A/(2 pi R), quindi
        **abbassare il raggio di iniezione e' il solo modo di ingrandire i
        fori** senza toccare nessuna pressione. Da 3.93 a 3.10 mm i fori
        passano da 0.53 a 0.67 mm, cioe' da sotto a sopra la soglia di stampa.
    """
    R_getti: float          # m, raggio strizzato del corpo centrale
    R_anello: float         # m, raggio esterno del condotto d'aria
    n_getti: int
    d_getto: float          # m
    x_getti: float          # m, negativo
    L_condotto: float       # m, tratto a sezione costante
    L_contrazione: float    # m
    R_plenum: float         # m
    L_plenum: float         # m
    #: Raggio del corpo centrale NEL PLENUM. Non e' R_getti: li' il corpo
    #: centrale deve contenere la filettatura del GPL, che ha un nocciolo di
    #: 8.8 mm. Si assottiglia poi a R_getti attraverso la contrazione, e la
    #: contrazione avviene quindi su ENTRAMBE le pareti, come in un vero
    #: ugello anulare. Fra plenum e condotto si passa da 224 a 37 mm2.
    R_corpo_plenum: float   # m
    t_monte: float          # m, parete della faccia di monte
    t_testa: float          # m, parete fra plenum ed esterno
    r_bore_gpl: float       # m, raggio del plenum del GPL nel corpo centrale
    x_rampa: float          # m, positivo: dove il corpo centrale torna a r_cb
    d_porta_aria: float     # m, foro STAMPATO, sottomisura: si fora e si maschia
    R_boss_aria: float      # m, raggio della faccia piana del bocchettone
    d_boss_aria: float      # m, diametro della bozza
    x_porta_aria_fissa: float  # m
    d_filetto_gpl: float    # m, foro STAMPATO per il G1/8
    L_filetto_gpl: float    # m
    L_cono: float           # m, raccordo esterno fra testa e mantello

    @property
    def x_condotto(self) -> float:
        return -self.L_condotto

    @property
    def x_contrazione(self) -> float:
        return self.x_condotto - self.L_contrazione

    @property
    def x_monte(self) -> float:
        return self.x_contrazione - self.L_plenum

    @property
    def x_faccia(self) -> float:
        """Faccia di monte del pezzo: e' la faccia che poggia sulla piastra di
        costruzione, quindi e' gia' spianata dopo il taglio dal supporto. E'
        li' che va la filettatura del GPL, e non serve nessuna bozza."""
        return self.x_monte - self.t_monte

    @property
    def R_testa(self) -> float:
        return self.R_plenum + self.t_testa

    @property
    def x_porta_aria(self) -> float:
        return self.x_porta_aria_fissa

    def lunghezza_confinata(self) -> float:
        """Quanto condotto resta a valle dei getti."""
        return -self.x_getti

    def area_anello(self) -> float:
        return math.pi * (self.R_anello ** 2 - self.R_getti ** 2)

    def area_getti(self) -> float:
        return self.n_getti * math.pi / 4.0 * self.d_getto ** 2

    def area_plenum_gpl(self) -> float:
        return math.pi * self.r_bore_gpl ** 2

    def verifica(self, lunghezza_mescolamento: float) -> list[str]:
        """Le condizioni che la geometria deve rispettare, e che il
        dimensionamento da solo non garantisce."""
        n = []
        if self.lunghezza_confinata() < lunghezza_mescolamento:
            n.append(
                f"i getti stanno a {-self.x_getti*1e3:.2f} mm dallo sbocco ma la "
                f"lunghezza di mescolamento e' {lunghezza_mescolamento*1e3:.2f} mm: "
                "il mescolamento finirebbe in camera, dove non c'e' piu' condotto "
                "e la correlazione non vale piu'")
        if self.x_getti < self.x_condotto:
            n.append("i getti cadono nella contrazione, dove la sezione cambia: "
                     "H non e' piu' definita")
        #: Il plenum deve essere abbastanza grande da non essere lui a decidere
        #: la portata. Regola: sezione di passaggio del plenum almeno doppia di
        #: quella che alimenta, altrimenti la velocita' nel plenum non e'
        #: trascurabile e la distribuzione in circonferenza non e' uniforme.
        passaggio = 2.0 * self.L_plenum * (self.R_plenum - self.R_corpo_plenum)
        if passaggio < 2.0 * self.area_anello():
            n.append(
                f"il plenum dell'aria offre {passaggio*1e6:.0f} mm2 di passaggio "
                f"contro {self.area_anello()*1e6:.0f} mm2 di condotto: sarebbe il "
                "plenum a distribuire male, e i getti vedrebbero un flusso "
                "trasversale diverso da un settore all'altro")
        if self.area_plenum_gpl() < 4.0 * self.area_getti():
            n.append(
                f"il plenum del GPL ({self.area_plenum_gpl()*1e6:.1f} mm2) e' meno "
                f"di 4 volte l'area dei getti ({self.area_getti()*1e6:.2f} mm2): "
                "i getti piu' lontani dall'ingresso riceverebbero meno")
        semi = 0.5 * self.d_boss_aria
        if self.x_porta_aria - semi < self.x_faccia - 1e-9 or \
                self.x_porta_aria + semi > 1e-9:
            n.append(
                f"la bozza dell'aria (diametro {self.d_boss_aria*1e3:.1f} mm centrata "
                f"a x = {self.x_porta_aria*1e3:.2f}) sborda dalla testa, che va da "
                f"{self.x_faccia*1e3:.2f} a 0: la faccia piana per la guarnizione non "
                "sarebbe appoggiata su niente")
        if 0.5 * self.d_filetto_gpl + 5.0e-4 > self.R_corpo_plenum:
            n.append(
                f"la filettatura del GPL (foro {self.d_filetto_gpl*1e3:.1f} mm) non ci "
                f"sta nel corpo centrale del plenum (raggio "
                f"{self.R_corpo_plenum*1e3:.1f} mm)")
        if self.R_getti - self.r_bore_gpl < 5.0e-4:
            n.append(
                f"parete fra GPL e aria di {(self.R_getti-self.r_bore_gpl)*1e3:.2f} "
                "mm: sotto il minimo dichiarato per l'SLM, e li' un poro mette il "
                "GPL nell'aria a monte dei getti")
        return n


def contorno_aria(t: TestaIniezione) -> list[tuple[float, float]]:
    """Poligono meridiano della cavita' d'aria: plenum, contrazione, condotto.

    La contrazione avviene su ENTRAMBE le pareti: da (6.0, 9.5) mm nel plenum a
    (3.10, 4.64) mm nel condotto. Il corpo centrale e' grosso nel plenum perche'
    li' dentro ci passa la filettatura del GPL, e si assottiglia dove servono i
    getti, che e' il solo punto in cui il raggio deve essere minimo.

    Nel TRATTO DEI GETTI le due pareti tornano cilindriche e restano tali fino
    allo sbocco: la correlazione di Holdeman presuppone un condotto a sezione
    costante, e i fori radiali su una superficie cilindrica non hanno
    sottosquadri.
    """
    return [
        (t.x_monte, t.R_corpo_plenum),
        (t.x_monte, t.R_plenum),
        (t.x_contrazione, t.R_plenum),
        (t.x_condotto, t.R_anello),
        (0.0, t.R_anello),
        (0.0, t.R_getti),
        (t.x_condotto, t.R_getti),
        (t.x_contrazione, t.R_corpo_plenum),
    ]


def contorno_testa(t: TestaIniezione, R_out: float) -> list[tuple[float, float]]:
    """Poligono meridiano del SOLIDO della testa, prima di scavarci dentro.

    Il raccordo conico all'esterno serve a passare dal raggio della testa a
    quello del mantello raffreddato. Costruendo lungo +x con la faccia di monte
    sulla piastra, e' una superficie rivolta in basso e va supportata: sono
    supporti ESTERNI e accessibili, quindi ammessi. L'alternativa - testa dello
    stesso diametro del mantello - non ha sbalzi ma pesa il doppio, ed e'
    materiale messo dove la pressione e' 7 bar e la tensione di cerchio 3 MPa.
    """
    x_cono = -t.L_cono
    return [
        (t.x_faccia, 0.0),
        (t.x_faccia, t.R_testa),
        (x_cono, t.R_testa),
        (0.0, R_out),
        (0.0, 0.0),
    ]


def _getti_gpl(grid, t: TestaIniezione, xp):
    """I fori radiali del GPL che attraversano il condotto d'aria.

    Partono DENTRO il plenum del GPL e arrivano oltre la parete esterna del
    condotto, cosi' la lunghezza non dipende da nessun conto: il taglio con la
    cavita' d'aria e col plenum lo fanno gli altri campi.
    """
    campo = None
    lunghezza = t.R_anello - t.r_bore_gpl + 4.0 * grid.spacing
    for k in range(t.n_getti):
        ang = 2.0 * math.pi * k / t.n_getti
        base = (t.x_getti, t.r_bore_gpl * math.cos(ang),
                t.r_bore_gpl * math.sin(ang))
        c = _cilindro_generico(grid, base, (0.0, math.cos(ang), math.sin(ang)),
                               0.5 * t.d_getto, lunghezza, xp)
        campo = c if campo is None else campo.union(c)
    return campo


def cavita_gpl(grid, t: TestaIniezione, xp):
    """Sede della filettatura G1/8, plenum nel corpo centrale, getti radiali.

    Il foro e' a GRADINI e non cilindrico: i primi `L_filetto_gpl` dalla faccia
    di monte sono il preforo della filettatura (stampato sottomisura, poi
    forato a 8.8 e maschiato G1/8), il resto e' il plenum vero e proprio.

    La filettatura sta sulla FACCIA DI MONTE, che e' la faccia appoggiata alla
    piastra di costruzione: dopo il taglio dal supporto e' gia' piana e
    perpendicolare all'asse, cioe' e' gia' la sede della guarnizione a legare.
    Non serve nessuna bozza e non si spende niente per averla.
    """
    sede = cylinder(grid, 0.5 * t.d_filetto_gpl,
                    t.x_faccia - 4.0 * grid.spacing,
                    t.x_faccia + t.L_filetto_gpl, xp)
    bore = cylinder(grid, t.r_bore_gpl,
                    t.x_faccia + t.L_filetto_gpl - grid.spacing,
                    t.x_getti + 0.5 * t.d_getto + 5.0e-4, xp)
    return sede.union(bore).union(_getti_gpl(grid, t, xp))


def bozza_aria(grid, t: TestaIniezione, xp):
    """La bozza (boss) su cui si maschia il G3/8 dell'aria.

    PERCHE' E' COSI' GRANDE, che e' la prima domanda che viene guardandola.
    L'aria sono 33.6 g/s a 8.7 kg/m3, cioe' 3.9 LITRI AL SECONDO. In un
    passaggio da 7.5 mm (un G1/4) fanno 87 m/s e 0.33 bar di perdita, che e'
    un TERZO del salto d'iniezione disponibile - e il salto d'iniezione qui e'
    tutto quello che c'e', perche' p_c e' stata fissata a p_sat/1.15 senza
    margini nascosti. Con un G3/8 (10 mm) si scende a 49 m/s e 0.105 bar,
    l'11 %. Vedi injector.perdita_raccordo.

    La bozza e' grande perche' il raccordo e' grande, e il raccordo e' grande
    perche' la portata volumetrica lo e'. Non c'e' niente da limare.

    La faccia e' piana e perpendicolare al raggio: e' la sede della guarnizione
    a legare, che e' il modo piu' semplice e affidabile di tenere 7 bar su una
    filettatura cilindrica.
    """
    dirz = (0.0, -1.0, 0.0)
    base = (t.x_porta_aria, t.R_boss_aria, 0.0)
    return _cilindro_generico(grid, base, dirz, 0.5 * t.d_boss_aria,
                              t.R_boss_aria - t.R_corpo_plenum, xp)


def cavita_aria(grid, t: TestaIniezione, xp):
    """Plenum + contrazione + condotto, piu' il foro filettato dell'aria."""
    interno = revolve_polygon(grid, contorno_aria(t), xp)
    porta = _cilindro_generico(
        grid,
        (t.x_porta_aria, (t.R_boss_aria + 4.0 * grid.spacing), 0.0),
        (0.0, -1.0, 0.0), 0.5 * t.d_porta_aria,
        t.R_boss_aria - t.R_plenum + 8.0 * grid.spacing, xp)
    return interno.union(porta)



def spessore_mantello(c: CircuitoRaffreddamento) -> float:
    """Spessore di mantello che un ramo richiede, dalla parete calda in fuori.

    DIFETTO GRAVE CORRETTO QUI, e trovato dalla prova di tenuta topologica.

    Fino alla revisione 0.4.0 questa somma non contava lo strato di RITORNO:
    valeva `parete_calda + lato + parete_fredda`, cioe' 2.8 mm per il ramo di
    gola, mentre il ramo di gola e' a U e ne occupa 0.8 + 0.6 + setto + 0.6 +
    1.0. Il mantello veniva quindi costruito spesso 3.4 mm (il massimo, dettato
    dal ramo di camera) e lo strato di ritorno, che ne chiedeva 4.0, veniva
    TAGLIATO dall'intersezione con il mantello.

    Tagliato dove? Sulla sua faccia esterna, cioe' proprio sulla pelle del
    pezzo. **Il circuito di raffreddamento della gola era aperto
    sull'atmosfera per tutta la sua lunghezza.**

    Perche' nessun controllo se n'era accorto: la superficie restava chiusa
    (`is_watertight` = True), il pezzo restava in una sola componente
    connessa, i volumi restavano plausibili, e il controllo della polvere
    diceva "0 sacche chiuse su 1 componente di vuoto" - che sembra un successo
    ed era invece il sintomo: UNA sola componente di vuoto significa che il
    circuito dell'acqua e l'atmosfera sono la stessa cosa.

    La prova di tenuta di `zefiro.sdf.tenuta` lo vede al primo colpo, perche'
    e' l'unica che chiede "chi comunica con chi" invece di "e' chiusa".
    """
    strati = c.parete_calda + c.lato + c.parete_fredda
    if c.ritorno:
        strati += c.setto_ritorno + c.lato
    return strati


@dataclass
class MotoreSDF:
    """Il pezzo completo, con i campi intermedi tenuti per la verifica."""
    grid: Grid
    solido: Field
    cavita_gas: Field
    canali: Field
    vuoti: Field          # canali + fori d'iniezione + attacchi: TUTTO il vuoto
    circuiti: tuple[CircuitoRaffreddamento, ...]
    testa: "TestaIniezione | None" = None
    collettore: "CollettoreAcqua | None" = None
    montaggio: "InterfacciaMontaggio | None" = None
    #: I vuoti tenuti SEPARATI, per poter misurare gli spessori di parete
    #: a coppie. Un unico campo unito direbbe solo che il pezzo e' chiuso,
    #: non quanto materiale resta fra due cavita' che non devono toccarsi.
    parti: dict = dc_field(default_factory=dict)
    note: list[str] = dc_field(default_factory=list)


def contorno_gas_chiuso(d, contour_x, contour_r,
                        testa: "TestaIniezione | None" = None) -> list[tuple[float, float]]:
    """Poligono della CAVITA' GASSOSA interna, chiusa sul piano del labbro.

    Chiusa e non aperta perche' serve un campo di distanza da una regione, e
    una regione aperta non ne ha uno. La chiusura sul piano del labbro viene
    poi ritagliata via: il getto esterno non e' parte del pezzo.

    Con la testa a getti trasversali il corpo centrale entra in camera
    STRIZZATO a `testa.R_getti` e riprende il suo raggio entro `testa.x_rampa`.
    La rampa e' dentro la camera, dove di area ce n'e' da vendere: mettere li'
    l'allargamento invece che nel condotto d'aria evita di far cambiare sezione
    al condotto proprio dove i getti devono trovarla costante.
    """
    R_c, R_lip, L_c, L_conv = d["R_c"], d["R_lip"], d["L_c"], d["L_conv"]
    r_cb = d["r_centerbody"]
    x_lip = L_c + L_conv
    r0 = r_cb if testa is None else testa.R_getti
    punti = [
        (0.0, r0),            # faccia d'iniezione, dal corpo centrale in fuori
        (0.0, R_c),
        (L_c, R_c),           # parete di camera
        (x_lip, R_lip),       # convergente fino al labbro
        (x_lip, r_cb),        # chiusura sul piano del labbro
    ]
    if testa is not None and testa.x_rampa > 0.0:
        punti.append((testa.x_rampa, r_cb))
    return punti


def contorno_plug(d, contour_x, contour_r,
                  testa: "TestaIniezione | None" = None) -> list[tuple[float, float]]:
    """Corpo centrale cilindrico piu' il contorno del plug, come solido pieno.

    Con la testa, il corpo centrale comincia alla faccia di monte con il raggio
    strizzato e si allarga a r_cb entro x_rampa: e' lo stesso profilo della
    parete interna della cavita' gassosa, visto dall'altra parte.
    """
    L_c, L_conv, r_cb = d["L_c"], d["L_conv"], d["r_centerbody"]
    x_lip = L_c + L_conv
    if testa is None:
        pts = [(0.0, 0.0), (0.0, r_cb)]
    else:
        pts = [(testa.x_faccia, 0.0), (testa.x_faccia, testa.R_getti),
               (0.0, testa.R_getti), (testa.x_rampa, r_cb)]
    x_throat = x_lip + contour_x[0]
    pts.append((x_throat, r_cb))
    for xx, rr in zip(contour_x[1:], contour_r[1:]):
        pts.append((x_lip + xx, rr))
    if contour_r[-1] > 1.0e-9:
        pts.append((x_lip + contour_x[-1], 0.0))
    return pts


def _raggio_parete(d, x: float) -> float:
    """Raggio della parete lato gas alla quota x: costante in camera, lineare
    nel convergente. Serve agli attacchi radiali, che devono partire dalla
    superficie esterna VERA e non da un raggio unico."""
    R_c, R_lip, L_c, L_conv = d["R_c"], d["R_lip"], d["L_c"], d["L_conv"]
    if x <= L_c:
        return R_c
    t = min((x - L_c) / max(L_conv, 1e-12), 1.0)
    return R_c + t * (R_lip - R_c)


def _cilindro_generico(grid, base, direzione, raggio, lunghezza, xp):
    """Cilindro di raggio `raggio` da `base` lungo `direzione` (versore)."""
    X, Y, Z = grid.coords(xp)
    d = np.asarray(direzione, dtype=float)
    d = d / np.linalg.norm(d)
    px, py, pz = X - base[0], Y - base[1], Z - base[2]
    t = px * d[0] + py * d[1] + pz * d[2]
    qx, qy, qz = px - t * d[0], py - t * d[1], pz - t * d[2]
    dr = xp.sqrt(qx * qx + qy * qy + qz * qz) - raggio
    dt = xp.maximum(-t, t - lunghezza)
    dentro = xp.maximum(dr, dt)
    fuori = xp.sqrt(xp.maximum(dr, 0.0) ** 2 + xp.maximum(dt, 0.0) ** 2)
    a = xp.where(dentro > 0, fuori, dentro)
    return Field(grid, xp.broadcast_to(a, grid.shape).astype(xp.float32).copy(), xp)


def _fori_iniezione(grid, d, t_face, xp):
    """I fori d'aria e di GPL attraverso la piastra, alle stesse quote della
    parametrizzazione L0: non sono un dettaglio grafico, sono le sezioni di
    passaggio che fissano la velocita' d'iniezione."""
    N = int(round(d["N_inj"]))
    R_inj = d["R_inj"]
    prof = t_face + 2.0 * grid.spacing
    campo = None
    for k in range(N):
        ang = 2.0 * math.pi * k / N
        # sfasati di mezzo passo fra loro, come nella disposizione di L1
        for raggio, dang in ((0.5 * d["d_ox"], 0.0),
                             (0.5 * d["d_fuel"], 0.5 * 2.0 * math.pi / N)):
            y = R_inj * math.cos(ang + dang)
            z = R_inj * math.sin(ang + dang)
            c = _cilindro_generico(grid, (-t_face - grid.spacing, y, z),
                                   (1.0, 0.0, 0.0), raggio, prof, xp)
            campo = c if campo is None else campo.union(c)
    return campo


def _porta_acqua(grid, x_centro, r_esterno_locale, penetrazione, raggio_foro,
                 sporgenza, xp, angolo: float = 0.0):
    """Attacco radiale che sbuca nel collettore.

    La LUNGHEZZA e' quella che serve a bucare la parete esterna ed entrare nel
    collettore, non di piu'. La prima versione passava una lunghezza pari al
    raggio: il foro trapanava fino all'asse e tagliava in due il mantello.
    Non si e' visto nei volumi (il pezzo restava plausibile), si e' visto in
    sezione.

    `r_esterno_locale` e' il raggio esterno del pezzo A QUELLA quota x, che nel
    convergente e' molto minore che in camera: un valore unico farebbe partire
    il foro dal vuoto o dentro il materiale.
    """
    rr = r_esterno_locale + sporgenza
    base = (x_centro, rr * math.cos(angolo), rr * math.sin(angolo))
    direzione = (0.0, -math.cos(angolo), -math.sin(angolo))
    return _cilindro_generico(grid, base, direzione, raggio_foro,
                              sporgenza + penetrazione, xp)


def _cilindrico(grid, angolo: float, xp):
    """(q, p) = coordinate nel piano y-z ruotate di `angolo`.

    `q` e' la distanza dall'asse misurata NELLA direzione dell'angolo (cioe' il
    raggio proiettato), `p` e' quella perpendicolare. Servono a costruire fori e
    attacchi radiali senza ruotare la griglia.
    """
    _, Y, Z = grid.coords(xp)
    ca, sa = math.cos(angolo), math.sin(angolo)
    return Y * ca + Z * sa, -Y * sa + Z * ca


def _vano_anulare(grid, gas, prof_min: float, prof_max: float,
                  x0: float, x1: float, pendenza: float, xp,
                  smussa_pavimento: bool = False):
    """Vano anulare fra due profondita' misurate dalla parete del gas.

    IL TETTO E' L'UNICA FACCIA CHE VA SMUSSATA, se il vano sta DENTRO materiale
    pieno. Costruendo lungo +x, la faccia a x maggiore e' la superficie inferiore
    del materiale sovrastante: uno sbalzo. La faccia a x minore e' invece la
    superficie superiore del materiale sottostante, e non chiede niente.

    `smussa_pavimento` serve per i vani dentro un COLLARE che sporge dal corpo:
    li' anche il pavimento e' sospeso, perche' sotto c'e' la rampa con cui il
    collare nasce. Un vano con due facce inclinate ha sezione a rombo, cioe'
    meta' area a parita' di lunghezza: e' il motivo per cui il raccordo si fa
    dentro il mantello e non in un collare.
    """
    X, _, _ = grid.coords(xp)
    radiale = xp.maximum(prof_min - gas.a, gas.a - prof_max)
    affondo = xp.maximum(gas.a - prof_min, 0.0)
    monte = (x0 + pendenza * affondo if smussa_pavimento else x0) - X
    valle = X - (x1 - pendenza * affondo)
    a = xp.maximum(radiale, xp.maximum(monte, valle))
    return Field(grid, xp.broadcast_to(a, grid.shape).astype(xp.float32).copy(), xp)


def _collare(grid, gas, prof_min: float, prof_max: float,
             x0: float, x1: float, pendenza: float, xp,
             prof_rampa: float | None = None):
    """Materiale aggiunto fuori dal mantello, per ospitare un vano anulare.

    LA RAMPA. Un collare che sporge in fuori ha una faccia rivolta verso la
    piastra: quella di monte. Va inclinata di `pendenza`, cioe' il collare
    comincia piu' tardi man mano che si allontana dal corpo. La faccia di valle
    non ha lo stesso problema: guarda dalla parte opposta.
    """
    X, _, _ = grid.coords(xp)
    radiale = xp.maximum(prof_min - gas.a, gas.a - prof_max)
    #: LA RAMPA SI MISURA DALLA STESSA PROFONDITA' DEL VANO CHE DEVE CONTENERE,
    #: non dalla propria. Il collare comincia un po' PIU' DENTRO del vano per
    #: agganciarsi al mantello, e contando la rampa da li' la sua faccia
    #: inferiore risultava inclinata a partire da un punto diverso: le due rette
    #: - pavimento del vano e faccia inferiore del collare - finivano a
    #: venticinque micron l'una dall'altra invece che a un millimetro, cioe'
    #: coincidenti per qualunque griglia. Il vano si apriva sull'esterno lungo
    #: tutta la generatrice.
    riferimento = prof_min if prof_rampa is None else prof_rampa
    affondo = xp.maximum(gas.a - riferimento, 0.0)
    monte = (x0 + pendenza * affondo) - X
    valle = X - x1
    a = xp.maximum(radiale, xp.maximum(monte, valle))
    return Field(grid, xp.broadcast_to(a, grid.shape).astype(xp.float32).copy(), xp)


def _foro_goccia(grid, x_centro: float, angolo: float, raggio: float,
                 q0: float, q1: float, pendenza: float, xp):
    """Foro RADIALE a sezione di goccia: cerchio piu' tetto a due falde.

    Un foro orizzontale tondo ha, in cima, una superficie rivolta in basso e
    tangente alla piastra: e' lo sbalzo peggiore possibile, e sopra i 3 mm circa
    cede. Aggiungendogli un tetto a due falde inclinate di `pendenza` il punto
    piu' alto diventa uno spigolo, e ogni strato sporge sul precedente di poco:
    il foro si stampa senza supporti, e dentro un foro i supporti non si tolgono.

    L'apice sta a `x_centro + pendenza * raggio`; la sezione utile e' quella del
    cerchio piu' quella del triangolo.
    """
    X, _, _ = grid.coords(xp)
    q, pp = _cilindrico(grid, angolo, xp)
    dx = X - x_centro
    lungo = xp.maximum(q0 - q, q - q1)
    tondo = xp.maximum(xp.sqrt(pp * pp + dx * dx) - raggio, lungo)
    #: il triangolo: |p| <= raggio - dx / pendenza, con dx >= 0
    falde = xp.maximum(xp.abs(pp) - xp.maximum(raggio - dx / pendenza, 0.0), -dx)
    tetto = xp.maximum(falde, lungo)
    a = xp.minimum(tondo, tetto)
    return Field(grid, xp.broadcast_to(a, grid.shape).astype(xp.float32).copy(), xp)


def _gambo(grid, pg: "Portagomma", xp, raggio_extra: float = 0.0):
    """Il gambo del portagomma piu' le sue creste, come solido."""
    R, ang = pg.R_posizione, pg.angolo
    r_gambo = 0.5 * pg.d_canna + raggio_extra
    r_cresta = 0.5 * (pg.d_canna + pg.interferenza) + raggio_extra
    campo = None

    def aggiungi(base, direzione, raggio, lunghezza):
        nonlocal campo
        c = _cilindro_generico(grid, base, direzione, raggio, lunghezza, xp)
        campo = c if campo is None else campo.union(c)

    if pg.verso == 0.0:
        #: radiale: sporge verso l'esterno lungo l'angolo
        dirz = (0.0, math.cos(ang), math.sin(ang))
        base = (pg.x_base, R * math.cos(ang), R * math.sin(ang))
        aggiungi(base, dirz, r_gambo, pg.lunghezza)
        for k in range(pg.n_creste):
            q = pg.lunghezza - (k + 0.5) * pg.passo_creste
            if q <= 0.0:
                break
            b = (pg.x_base, (R + q) * math.cos(ang), (R + q) * math.sin(ang))
            aggiungi(b, dirz, r_cresta, 1.2e-3)
    else:
        dirz = (pg.verso, 0.0, 0.0)
        base = (pg.x_base, R * math.cos(ang), R * math.sin(ang))
        aggiungi(base, dirz, r_gambo, pg.lunghezza)
        for k in range(pg.n_creste):
            q = pg.lunghezza - (k + 0.5) * pg.passo_creste
            if q <= 0.0:
                break
            aggiungi((pg.x_base + pg.verso * q, base[1], base[2]), dirz,
                     r_cresta, 1.2e-3)
    return campo


def _mensola(grid, pg: "Portagomma", altezza: float, pendenza: float, xp):
    """Il rinforzo sotto un portagomma RADIALE.

    Un cilindro orizzontale lungo 16 mm e largo 13 ha, sotto, una generatrice
    tangente alla piastra: e' uno sbalzo che non si sostiene. La mensola e' una
    lama nel piano che contiene l'asse del motore e quello dell'attacco, con il
    bordo inclinato di `pendenza`: sale col pezzo e regge il gambo. Non e'
    ornamento e non e' un supporto di stampa da togliere - e' materiale del
    pezzo, e regge anche il tiro della canna.
    """
    X, _, _ = grid.coords(xp)
    q, pp = _cilindrico(grid, pg.angolo, xp)
    spessore = 2.0e-3
    dq = q - pg.R_posizione
    #: il bordo inferiore parte dal corpo e sale: x >= x_base - (altezza - dq)/1
    bordo = (pg.x_base - altezza) + dq / pendenza - X
    a = xp.maximum(xp.abs(pp) - 0.5 * spessore,
                   xp.maximum(bordo, xp.maximum(-dq, dq - pg.lunghezza)))
    a = xp.maximum(a, X - pg.x_base)
    return Field(grid, xp.broadcast_to(a, grid.shape).astype(xp.float32).copy(), xp)


def _collettore(grid, gas, mantello, c, x_centro: float, xp,
                strato: str = "andata"):
    """Gola anulare che raccoglie tutti i canali di un ramo, alla stessa
    profondita' dei canali. E' l'elemento che rende il circuito un circuito."""
    X, _, _ = grid.coords(xp)
    # Con la U il collettore deve unire ANCHE lo strato di ritorno, quindi si
    # estende in profondita' fino a coprirli entrambi.
    prof_and = c.parete_calda + 0.5 * c.lato
    prof_rit = c.parete_calda + c.lato + c.setto_ritorno + 0.5 * c.lato
    if strato == "andata":
        profondita, semi = prof_and, 0.5 * c.lato
    elif strato == "ritorno":
        profondita, semi = prof_rit, 0.5 * c.lato
    elif strato == "entrambi":
        profondita = 0.5 * (prof_and + prof_rit)
        semi = 0.5 * (prof_rit - prof_and) + 0.5 * c.lato
    else:
        raise ValueError(f"strato {strato!r} sconosciuto")
    # Il collettore e' leggermente PIU' PROFONDO dei canali che unisce, cosi'
    # la banda del canale ci sta dentro tutta e i due bordi non si incontrano
    # tangenti. All'incrocio tangente marching cubes lasciava sacche di vuoto
    # chiuse da un voxel: polvere che non esce.
    radiale = xp.abs(gas.a - profondita) - (semi + COLLETTORE_MARGINE)
    # IL SOFFITTO DI UNA CAVA STA A x MAGGIORE, non minore.
    #
    # Costruendo lungo +x, il materiale sopra il vuoto e' quello a x piu'
    # grande, e la sua faccia esposta guarda verso la piastra: normale -x,
    # sbalzo. La faccia a x minore e' invece la superficie SUPERIORE del
    # materiale sottostante, e non ha bisogno di niente.
    #
    # Alla prima correzione avevo smussato la faccia sbagliata, e i 2.4 cm2 di
    # tetto piatto erano rimasti esattamente dov'erano. Il verso delle normali
    # di marching cubes non e' un dettaglio da assumere: e' verificato su una
    # sfera in tests/test_printability.py, dove la risposta si conosce.
    #
    # Il soffitto si smussa a 45 gradi facendolo arretrare man mano che si va
    # in profondita': in sezione meridiana il collettore diventa un cuneo.
    # il riferimento dello smusso deve essere il bordo VERO del collettore,
    # margine compreso: usando quello nominale lo smusso partiva un margine
    # piu' in la' e il primo tratto di soffitto restava piatto.
    profondita_min = profondita - semi - COLLETTORE_MARGINE
    monte = (x_centro - 0.5 * COLLETTORE_LARGHEZZA) - X
    valle = X - (x_centro + 0.5 * COLLETTORE_LARGHEZZA
                 - COLLETTORE_PENDENZA * (gas.a - profondita_min))
    a = xp.maximum(radiale, xp.maximum(monte, valle))
    campo = Field(grid, xp.broadcast_to(a, grid.shape).astype(xp.float32).copy(), xp)
    return campo.intersection(mantello)


def _capo_del_collettore(grid, gas, mantello, d, col: "CollettoreAcqua", t_max: float,
                         anello, x_collettore, n_fori, d_foro, prof_collettore,
                         pg: "Portagomma", nome: str, xp,
                         pavimento_inclinato: bool, inclinati: bool = False):
    """Un capo del collettore: collare, vano anulare, fori radiali, portagomma.

    Ingresso e uscita sono lo stesso oggetto costruito due volte, con quote
    diverse. Scriverlo una volta sola non e' economia di righe: e' la garanzia
    che un difetto trovato a un capo sia corretto anche all'altro. Il primo
    difetto di questo tipo - l'apice della goccia che usciva dal collare - era
    stato corretto solo all'uscita.
    """
    p_ = col.pendenza_tetto
    x_v0, x_v1, h_an, x_c0 = anello
    prof0, prof1 = t_max, t_max + h_an
    collare = _collare(grid, gas, prof0 - 2.0 * grid.spacing, prof1 + col.parete,
                       x_c0, x_v1 + col.parete, p_, xp, prof_rampa=prof0)
    gambo = _gambo(grid, pg, xp, raggio_extra=col.parete)
    mensola = _mensola(grid, pg, 3.0 * col.parete + h_an, p_, xp)
    solido = collare.union(gambo).union(mensola)

    #: OGNI VUOTO SI RITAGLIA SUL MATERIALE CHE LO DEVE CONTENERE.
    #:
    #: Non e' prudenza: e' l'unico modo di essere sicuri. Un foro radiale ha una
    #: lunghezza calcolata, e la lunghezza giusta dipende dalla forma del
    #: contenitore alla quota in cui il foro la incontra. Il collare non ha un
    #: raggio esterno costante - nasce dal corpo con una rampa a 52 gradi, e in
    #: cima ci arriva sei millimetri piu' a valle. I fori d'ingresso, tagliati
    #: alla profondita' massima dell'anello, uscivano dal collare dalla parte di
    #: monte, dove il collare a quella profondita' ancora non c'e': il circuito
    #: dell'acqua si apriva sull'esterno a un decimo di millimetro dalla faccia
    #: d'iniezione. Superficie chiusa, un pezzo solo, volumi giusti.
    contenitore = collare.union(gambo).union(mantello)
    parti = {}
    parti[f"acqua_anello_{nome}"] = _vano_anulare(
        grid, gas, prof0, prof1, x_v0, x_v1, p_, xp,
        smussa_pavimento=pavimento_inclinato).intersection(collare)

    R_par = _raggio_parete(d, x_collettore)
    q0 = R_par + prof_collettore
    #: IL FORO FINISCE A META' DELL'ANELLO, non al suo bordo esterno. Portandolo
    #: fino a `prof1` il suo raggio lo faceva sporgere oltre la parete esterna
    #: del collare: il ritaglio sul contenitore lo teneva dentro, ma di un
    #: quarto di millimetro, cioe' meno del passo di griglia. Una tenuta che
    #: dipende dalla risoluzione non e' una tenuta.
    q1 = R_par + prof0 + 0.5 * h_an
    fori = None
    for k in range(n_fori):
        ang = 2.0 * math.pi * k / n_fori
        if inclinati:
            #: FORI INCLINATI, non radiali, e la ragione e' geometrica.
            #:
            #: A monte il collare nasce dal corpo con una rampa a 52 gradi, e
            #: alla quota del collettore d'ingresso - un millimetro e mezzo
            #: dalla faccia - la rampa e' appena cominciata: sopra al mantello
            #: non c'e' ancora niente. Un foro radiale li' esce dal fianco della
            #: rampa e apre il circuito sull'esterno. Un foro che sale con la
            #: STESSA pendenza della rampa le resta parallelo, quindi dentro, e
            #: per giunta e' autosostentato per definizione: la sua inclinazione
            #: e' esattamente quella che il processo ammette.
            dr = q1 - q0
            dx = p_ * dr
            dirz = (dx, dr * math.cos(ang), dr * math.sin(ang))
            lung = math.hypot(dx, dr)
            f = _cilindro_generico(
                grid, (x_collettore - grid.spacing,
                       q0 * math.cos(ang), q0 * math.sin(ang)),
                dirz, 0.5 * d_foro, lung + 2.0 * grid.spacing, xp)
        else:
            f = _foro_goccia(grid, x_collettore, ang, 0.5 * d_foro,
                             q0 - grid.spacing, q1 + grid.spacing, p_, xp)
        fori = f if fori is None else fori.union(f)
    parti[f"acqua_fori_{nome}"] = fori.intersection(contenitore)

    parti[f"acqua_{nome}"] = _foro_goccia(
        grid, pg.x_base, pg.angolo, 0.5 * pg.d_foro,
        _raggio_parete(d, pg.x_base) + t_max,
        pg.R_posizione + pg.lunghezza + 2.0 * grid.spacing, p_, xp
    ).intersection(contenitore)
    return solido, parti


def _braccio(grid, angolo, q0, q1, larghezza, x0, x1, xp):
    """Una lama radiale piatta, dal corpo alla piazzola.

    Sta tutta nel piano della faccia, che in stampa e' la PIASTRA: e' la prima
    cosa che si costruisce e non ha nessuna superficie rivolta in basso. Una
    flangia anulare continua peserebbe trentacinque grammi; quattro bracci ne
    pesano dieci e reggono lo stesso carico, perche' il carico e' assiale e i
    bracci lavorano a compressione nel proprio piano.
    """
    X, _, _ = grid.coords(xp)
    q, p_ = _cilindrico(grid, angolo, xp)
    a = xp.maximum(xp.abs(p_) - 0.5 * larghezza,
                   xp.maximum(q0 - q, q - q1))
    a = xp.maximum(a, xp.maximum(x0 - X, X - x1))
    return Field(grid, xp.broadcast_to(a, grid.shape).astype(xp.float32).copy(), xp)


def costruisci_film(grid, f: "FessuraFilm", xp):
    """Condotto, dosatura e fessura, come un unico vuoto.

    Tutto quello che serve e' assialsimmetrico tranne il foro di dosatura, e per
    la parte assialsimmetrica si usa un poligono rivoltato: e' ESATTO, non
    discretizzato in facce, ed e' il modo giusto di disegnare una fessura alta
    mezzo millimetro su una griglia che di millimetri ne risolve un quarto.
    """
    ca, sa = math.cos(f.angolo), math.sin(f.angolo)
    #: L'USCITA E' RASO SULLA SUPERFICIE, non sotto.
    #:
    #: La prima versione faceva finire la fessura `spessore_labbro` sotto il
    #: raggio del plug, per non lasciare un labbro a lama. Il risultato era che
    #: la fessura non sbucava affatto: restava un condotto cieco da 217 mm3
    #: dentro il corpo centrale - polvere che non esce, e film che non esce
    #: nemmeno lui. La prova di tenuta l'ha visto subito come "vuoto non
    #: rivendicato", che e' la sua diagnosi per una cavita' in cui non c'e'
    #: nessuna sonda perche' non e' collegata a niente.
    #:
    #: Il labbro sottile resta, e si accetta: essendo la fessura inclinata di
    #: `angolo`, il labbro raggiunge lo spessore di progetto gia'
    #: `spessore_labbro / tan(angolo)` prima del bordo. Gli ultimi decimi di
    #: millimetro usciranno arrotondati dalla stampa, come esce arrotondato
    #: qualunque spigolo vivo, e per un labbro di fessura va benissimo.
    r_est = f.r_uscita
    h_r = f.altezza / ca                      # apertura misurata in raggio
    r_int = r_est - h_r
    #: quanto indietro va la fessura per scendere dal labbro al condotto
    L = (r_int - f.r_condotto) / sa if sa > 1e-9 else 0.0
    x0 = f.x_uscita - L * ca
    #: la fessura e' una BANDA fra due rette parallele inclinate di `angolo`:
    #: quattro vertici, non di piu'. Il bordo esterno va da (x_uscita, r_est)
    #: indietro fino a (x0, r_condotto + h_r); quello interno gli sta sotto di
    #: `h_r` in raggio, cioe' di `altezza` misurata normale alla fessura.
    fessura = revolve_polygon(grid, [
        (f.x_uscita, r_est),
        (f.x_uscita, r_int),
        (x0, f.r_condotto),
        (x0, f.r_condotto + h_r),
    ], xp)
    condotto = revolve_polygon(grid, [
        (f.x_condotto_0, 0.0),
        (x0 + 2.0 * grid.spacing, 0.0),
        (x0 + 2.0 * grid.spacing, f.r_condotto),
        (f.x_condotto_0, f.r_condotto),
    ], xp)
    dosatura = _cilindro_generico(
        grid, (f.x_dosatura - 2.0 * grid.spacing, 0.0, 0.0), (1.0, 0.0, 0.0),
        0.5 * f.d_dosatura, (f.x_condotto_0 - f.x_dosatura) + 4.0 * grid.spacing, xp)
    return fessura.union(condotto).union(dosatura)


def costruisci_montaggio(grid, testa, m: "InterfacciaMontaggio", xp):
    """Bracci, piazzole e fori passanti. (solido, vuoto)."""
    solido = vuoto = None
    x0, x1 = testa.x_faccia, testa.x_faccia + m.spessore
    for k in range(m.n_bracci):
        ang = m.angolo_offset + 2.0 * math.pi * k / m.n_bracci
        y, z = m.R_bulloni * math.cos(ang), m.R_bulloni * math.sin(ang)
        pad = _cilindro_generico(grid, (x0, y, z), (1.0, 0.0, 0.0),
                                 0.5 * m.d_piazzola, m.spessore, xp)
        braccio = _braccio(grid, ang, m.r_attacco, m.R_bulloni,
                           m.larghezza_braccio, x0, x1, xp)
        pezzo = pad.union(braccio)
        solido = pezzo if solido is None else solido.union(pezzo)
        #: il foro e' PASSANTE: attraversa braccio e piazzola e si apre
        #: sull'esterno da entrambe le parti, quindi non intrappola polvere.
        foro = _cilindro_generico(
            grid, (x0 - 2.0 * grid.spacing, y, z), (1.0, 0.0, 0.0),
            0.5 * m.d_foro, m.spessore + 4.0 * grid.spacing, xp)
        vuoto = foro if vuoto is None else vuoto.union(foro)
    return solido, vuoto


def costruisci_collettore(grid, gas, mantello, d, circuiti, testa,
                          col: "CollettoreAcqua", t_max: float, xp):
    """Il collettore dell'acqua: che cosa si aggiunge e che cosa si toglie.

    Restituisce (solido_da_aggiungere, parti_vuote). Le parti vuote entrano nei
    CANALI e non fra i vuoti generici, perche' devono passare le stesse
    verifiche: un vano che buca la parete calda sprizza acqua nel gas
    esattamente come un canale.
    """
    p_ = col.pendenza_tetto
    cam = next(c for c in circuiti if c.nome == "camera")
    gola = next(c for c in circuiti if c.nome == "gola")

    # --- raccordo camera -> gola: vano anulare DENTRO il mantello ------------
    x_r0, x_r1, prof_r0, prof_r1 = col.raccordo
    parti = {"acqua_raccordo": _vano_anulare(
        grid, gas, prof_r0, prof_r1, x_r0, x_r1, p_, xp).intersection(mantello)}

    #: INGRESSO. Il suo collare comincia a x = 0, dove comincia il mantello:
    #: prima non c'e' corpo a cui attaccarsi, c'e' il cono della testa. La rampa
    #: con cui nasce e' quindi tutta a valle, e il pavimento del vano deve
    #: seguirla - sezione a rombo invece che a trapezio, meta' area a parita' di
    #: lunghezza. E' il prezzo di stare all'estremita' di monte del pezzo.
    prof_cam = cam.parete_calda + 0.5 * cam.lato
    s_in, p_in = _capo_del_collettore(
        grid, gas, mantello, d, col, t_max, col.anello_ingresso, col.x_collettore_ingresso,
        col.n_fori_ingresso, col.d_foro_ingresso, prof_cam, col.ingresso,
        "ingresso", xp, pavimento_inclinato=True, inclinati=True)

    #: USCITA. Qui il collare puo' cominciare molto prima del vano, perche' a
    #: monte c'e' mantello in abbondanza: la rampa finisce prima che il vano
    #: cominci, e il pavimento puo' essere piatto.
    prof_rit = (gola.parete_calda + gola.lato + gola.setto_ritorno
                + 0.5 * gola.lato)
    s_out, p_out = _capo_del_collettore(
        grid, gas, mantello, d, col, t_max, col.anello_uscita, col.x_collettore_uscita,
        col.n_fori_uscita, col.d_foro_uscita, prof_rit, col.uscita,
        "uscita", xp, pavimento_inclinato=False)

    parti.update(p_in)
    parti.update(p_out)
    return s_in.union(s_out), parti


def costruisci(
    d,
    contour_x: Sequence[float],
    contour_r: Sequence[float],
    circuiti: Sequence[CircuitoRaffreddamento],
    passo_griglia: float,
    raccordo: float = 0.0,
    xp=np,
    testa: "TestaIniezione | None" = None,
    collettore: "CollettoreAcqua | None" = None,
    montaggio: "InterfacciaMontaggio | None" = None,
    film: "FessuraFilm | None" = None,
) -> MotoreSDF:
    """Assembla il motore. Ogni passaggio lascia il proprio campo, cosi' la
    verifica puo' misurare i pezzi e non solo il risultato."""
    R_c = d["R_c"]
    L_c, L_conv = d["L_c"], d["L_conv"]
    x_lip = L_c + L_conv
    x_base = x_lip + contour_x[-1]
    t_max = max(spessore_mantello(c) for c in circuiti)
    t_face = max(2.0 * d["t_wall"], 2.0e-3)

    R_out = R_c + t_max
    x_monte = -t_face if testa is None else testa.x_faccia
    #: La griglia deve contenere anche la BOZZA dell'aria, che sporge oltre il
    #: raggio del mantello. Senza, il solido toccava il bordo, marching cubes
    #: tagliava la superficie e la mesh usciva aperta - con un errore di volume
    #: del 6 % che era esattamente il pezzo di bozza tagliato via.
    R_grid = R_out if testa is None else max(R_out, testa.R_boss_aria)
    if montaggio is not None:
        R_grid = max(R_grid, montaggio.R_bulloni + 0.5 * montaggio.d_piazzola
                     + 2.0 * passo_griglia)
    #: LA SCATOLA NON E' UN CUBO. I portagomma sporgono da UNA sola parte, e
    #: prendere la loro punta come semilato in tutte le direzioni quadruplica i
    #: voxel: a 0.25 mm sono 137 MB per campo invece di 42, e il processo muore
    #: per esaurimento di memoria senza scrivere niente. Si costruisce quindi la
    #: scatola vera, che e' asimmetrica come il pezzo.
    lo_y = hi_y = lo_z = hi_z = R_grid
    if collettore is not None:
        R_grid = max(R_grid,
                     max(collettore.anello_ingresso[2], collettore.anello_uscita[2])
                     + R_out + collettore.parete)
        lo_y = hi_y = lo_z = hi_z = R_grid
        for pg in (collettore.ingresso, collettore.uscita):
            r_cresta = 0.5 * (pg.d_canna + pg.interferenza) + collettore.parete
            if pg.verso == 0.0:
                #: RADIALE: sporge di tutta la sua lunghezza oltre l'anello. Se
                #: la griglia non lo contiene, marching cubes lo taglia e la
                #: mesh esce aperta proprio dove c'e' un foro passante.
                punta = pg.R_posizione + pg.lunghezza + r_cresta
                cy, cz = math.cos(pg.angolo), math.sin(pg.angolo)
                hi_y = max(hi_y, punta * cy + r_cresta)
                lo_y = max(lo_y, -(punta * cy) + r_cresta)
                hi_z = max(hi_z, punta * cz + r_cresta)
                lo_z = max(lo_z, -(punta * cz) + r_cresta)
            else:
                lo_y = hi_y = max(hi_y, pg.R_posizione + r_cresta)
                lo_z = hi_z = max(hi_z, pg.R_posizione + r_cresta)
                x_monte = min(x_monte, pg.x_base + pg.verso * pg.lunghezza)
    lo = (x_monte - 3 * passo_griglia, -lo_y, -lo_z)
    hi = (x_base + 3 * passo_griglia, hi_y, hi_z)
    grid = Grid.bounding(lo, hi, passo_griglia, margin=4 * passo_griglia)

    gas = revolve_polygon(grid, contorno_gas_chiuso(d, contour_x, contour_r, testa), xp)
    plug = revolve_polygon(grid, contorno_plug(d, contour_x, contour_r, testa), xp)

    # --- mantello: guscio esterno alla cavita', ritagliato al labbro -------- #
    X, _, _ = grid.coords(xp)
    piano_labbro = Field(grid, xp.broadcast_to(X - x_lip, grid.shape)
                         .astype(xp.float32).copy(), xp)
    # Il guscio di una regione avvolge TUTTO il suo contorno, e il contorno
    # della cavita' gassosa comprende anche la superficie del corpo centrale:
    # il guscio "esterno" includeva quindi uno strato DENTRO il plug, e i
    # canali ci finivano dentro svuotandolo. Si toglie esplicitamente il plug.
    # Visto in sezione, non nei numeri: il volume totale era plausibile.
    # Il mantello esiste SOLO fra la faccia d'iniezione e il labbro. Senza il
    # taglio a monte, il collettore piazzato a x = 0 sconfinava dentro la
    # piastra e vi scavava un anello attorno alla radice del corpo centrale,
    # staccando un disco di materiale: 5432 triangoli in una seconda
    # componente connessa. Il volume totale non se ne accorgeva.
    piano_iniezione = Field(grid, xp.broadcast_to(-X, grid.shape)
                            .astype(xp.float32).copy(), xp)
    mantello = (gas.shell(t_max, verso="esterno")
                .intersection(piano_labbro)
                .intersection(piano_iniezione)
                .difference(plug))

    # --- testa d'iniezione --------------------------------------------------- #
    # Senza `testa` resta la vecchia piastra piena: serve ai test che
    # verificano il mantello e non hanno bisogno dell'iniettore.
    if testa is None:
        piastra = cylinder(grid, R_out, -t_face, 0.0, xp)
    else:
        piastra = revolve_polygon(grid, contorno_testa(testa, R_out), xp)
        piastra = piastra.union(bozza_aria(grid, testa, xp))

    corpo = mantello.union(piastra)
    corpo = corpo.smooth_union(plug, raccordo) if raccordo > 0 else corpo.union(plug)

    # --- canali: definiti RISPETTO al campo del gas, quindi conformi -------- #
    #
    # DIFETTO TROVATO GUARDANDO LA SEZIONE, non i numeri. La condizione
    # `|G - profondita| <= h/2` e' soddisfatta su ENTRAMBI i lati della parete:
    # un punto a 1.2 mm dalla cavita' gassosa lo e' sia andando verso l'esterno
    # (dove il canale ci vuole) sia andando verso l'interno, dentro il corpo
    # centrale (dove svuota il plug). Il campo di distanza non sa da che lato
    # sta la parete, e non puo' saperlo: e' una distanza, non ha un verso.
    #
    # Il rimedio e' dirglielo: il canale esiste solo DENTRO il mantello. Cosi'
    # il vincolo geometrico e' esplicito invece che sperato.
    canali = None
    parti: dict = {}
    note_passo: list[str] = []
    #: (nome, x_centro) di ogni collettore, per verificare a posteriori che
    #: nessuno sporga oltre i due piani che limitano il mantello.
    collettori_x: list[tuple[str, float]] = []
    for c in circuiti:
        ch = helical_channels(
            grid, gas,
            profondita=c.parete_calda + 0.5 * c.lato,
            larghezza=c.lato, altezza=c.lato,
            n_canali=c.n_canali, passo=c.passo_elica,
            x_inizio=c.x_inizio, x_fine=c.x_fine, xp=xp,
        )
        ch = ch.intersection(mantello)
        parti[f"{c.nome}_andata"] = ch
        # Il raggio che conta e' quello del CANALE, non quello esterno del
        # pezzo: e' li' che sta la superficie da sostenere. Con il raggio
        # esterno il controllo bocciava progetti sani.
        prof_esterna = (c.parete_calda + c.lato + c.setto_ritorno + 0.5 * c.lato
                        if c.ritorno else c.parete_calda + 0.5 * c.lato)
        raggio_canale = max(_raggio_parete(d, c.x_inizio),
                            _raggio_parete(d, c.x_fine)) + prof_esterna
        passo_minimo = passo_elica_minimo(raggio_canale, MARGINE_ANGOLO_DEG)
        if abs(c.passo_elica) < passo_minimo:
            note_passo.append(
                f"il ramo '{c.nome}' ha passo {abs(c.passo_elica)*1e3:.0f} mm contro i "
                f"{passo_minimo*1e3:.0f} mm minimi a raggio {raggio_canale*1e3:.1f} mm "
                f"(soglia con margine: {MARGINE_ANGOLO_DEG:.0f} gradi): il tetto dei "
                "canali e' uno sbalzo, e dentro un canale non si possono mettere supporti"
            )
        if c.ritorno:
            prof_rit = (c.parete_calda + c.lato + c.setto_ritorno + 0.5 * c.lato)
            rit = helical_channels(
                grid, gas, profondita=prof_rit,
                larghezza=c.lato, altezza=c.lato,
                # STESSO passo, non opposto. Con eliche controrotanti andata e
                # ritorno si incrociano, e a ogni incrocio il setto fra i due
                # strati si riduce a una lamella: marching cubes vi lasciava
                # 116 isole di materiale da un voxel, cioe' polvere sinterizzata
                # che si stacca. Concordi, il ritorno corre esattamente sopra
                # l'andata e il setto resta spesso uguale ovunque.
                n_canali=c.n_canali, passo=c.passo_elica,
                x_inizio=c.x_inizio + c.sfalsamento_ritorno, x_fine=c.x_fine, xp=xp,
            ).intersection(mantello)
            parti[f"{c.nome}_ritorno"] = rit
            ch = ch.union(rit)
        # Collettori. Per un circuito a U ne servono TRE, non due:
        #   * a x_inizio, alla profondita' della sola ANDATA;
        #   * a x_inizio + sfalsamento, alla profondita' del solo RITORNO;
        #   * a x_fine, che li unisce entrambi ed e' l'inversione a U.
        # Con un collettore unico a monte, che li univa tutti e due, l'acqua
        # entrava da un attacco e usciva dall'altro senza passare per i canali:
        # un cortocircuito idraulico invisibile in ogni vista.
        if c.ritorno:
            ch = ch.union(_collettore(grid, gas, mantello, c, c.x_inizio, xp,
                                      strato="andata"))
            collettori_x.append((f"{c.nome}_andata", c.x_inizio))
            ch = ch.union(_collettore(grid, gas, mantello, c,
                                      c.x_inizio + c.sfalsamento_ritorno, xp,
                                      strato="ritorno"))
            collettori_x.append((f"{c.nome}_ritorno",
                                 c.x_inizio + c.sfalsamento_ritorno))
            # L'inversione a U sta TUTTA a monte di x_fine, non a cavallo.
            #
            # Un collettore e' largo COLLETTORE_LARGHEZZA e centrato su
            # x_centro: centrandolo su x_fine sporgeva di 1.5 mm oltre, e
            # x_fine sta gia' a soli 0.5 mm dal labbro. Il mantello finisce al
            # labbro, quindi il collettore veniva TAGLIATO dal piano del labbro
            # e si apriva sulla faccia: il circuito della gola comunicava con
            # l'atmosfera da li'. Trovato dalla prova di tenuta, invisibile a
            # tutto il resto (superficie chiusa, un solo pezzo, volumi giusti).
            #
            # L'alternativa - arretrare x_fine di 2 mm - avrebbe tolto
            # raffreddamento proprio al labbro, dove q vale 4.3 MW/m2. Cosi'
            # invece i canali arrivano dove arrivavano e a spostarsi e' solo
            # l'inversione.
            x_u = c.x_fine - 0.5 * COLLETTORE_LARGHEZZA
            ch = ch.union(_collettore(grid, gas, mantello, c, x_u, xp,
                                      strato="entrambi"))
            collettori_x.append(("gola_U", x_u))
        else:
            for x_col in (c.x_inizio, c.x_fine):
                ch = ch.union(_collettore(grid, gas, mantello, c, x_col, xp,
                                          strato="andata"))
                collettori_x.append((f"{c.nome}@{x_col*1e3:.1f}", x_col))
        canali = ch if canali is None else canali.union(ch)

    # --- attacchi dell'acqua ------------------------------------------------ #
    # Entrano in `canali` e non fra i vuoti generici, perche' devono
    # sottostare alla stessa verifica: un attacco che buca la parete sprizza
    # acqua nel gas esattamente come un canale. La prima versione li teneva
    # fuori dal controllo, e l'attacco di uscita del ramo gola bucava davvero.
    for i_c, c in enumerate(circuiti if collettore is None else ()):
        prof_and = c.parete_calda + 0.5 * c.lato
        prof_rit = c.parete_calda + c.lato + c.setto_ritorno + 0.5 * c.lato
        if c.ritorno:
            # L'attacco dell'ANDATA sta a x_inizio, dove lo strato di ritorno
            # non e' ancora cominciato: cosi' non lo attraversa e i due strati
            # restano separati. Quello del RITORNO sta piu' a valle e penetra
            # solo fino allo strato esterno.
            stazioni = [(c.x_inizio, prof_and),
                        (c.x_inizio + c.sfalsamento_ritorno, prof_rit)]
        else:
            stazioni = [(c.x_inizio, prof_and), (c.x_fine, prof_and)]

        for k, (x_col, prof) in enumerate(stazioni):
            # r_est generoso: con il taglio sul campo, partire piu' fuori del
            # necessario non fa danno, mentre partire troppo dentro si'.
            r_est = _raggio_parete(d, x_col) + 1.4 * t_max
            penetrazione = 1.4 * t_max - prof + 0.5 * c.lato + grid.spacing
            # Gli attacchi di collettori consecutivi sono SFALSATI: due
            # collettori vicini in x, con attacchi allineati, si compenetrano.
            # Su questo motore succedeva davvero fra l'uscita della camera e
            # l'ingresso della gola, a 4 mm di distanza con fori da 5.6 mm: i
            # due rami paralleli diventavano un ramo solo.
            base = (2.0 * math.pi / PORTE_PER_COLLETTORE) * SFALSAMENTO_COLLETTORI \
                * (2 * i_c + k)
            # L'attacco viene TAGLIATO CON IL CAMPO alla profondita' del suo
            # strato, invece di fidarsi della lunghezza calcolata.
            #
            # Il motivo e' che su una parete CONICA la distanza radiale non e'
            # quella normale: a 37 gradi di semiapertura, 4 mm di parete misurati
            # normalmente sono 5 mm misurati in raggio. L'attacco partiva percio'
            # un millimetro DENTRO il materiale e sfondava di altrettanto,
            # arrivando nello strato di andata e cortocircuitando la U.
            # Tagliandolo col campo, si ferma alla profondita' giusta su
            # qualunque forma di parete, cono o cilindro che sia.
            # ATTENZIONE AL VERSO: `gas.a` e' la profondita' misurata DALLA
            # parete verso l'esterno, quindi cresce allontanandosi dal gas.
            # L'attacco arriva da fuori (profondita' grande) e scende: per
            # fermarlo al proprio strato si tiene `gas.a >= limite`, non <=.
            # Con il verso sbagliato il taglio teneva proprio la parte che
            # doveva togliere, e l'attacco del ritorno continuava a sfondare
            # nello strato di andata.
            limite = prof - 0.5 * c.lato
            clip = Field(grid, xp.broadcast_to(limite - gas.a, grid.shape)
                         .astype(xp.float32).copy(), xp)
            for j in range(PORTE_PER_COLLETTORE):
                angolo = base + 2.0 * math.pi * j / PORTE_PER_COLLETTORE
                porta = _porta_acqua(grid, x_col, r_est, penetrazione,
                                     0.5 * PORTA_DIAMETRO, 3.0 * grid.spacing, xp,
                                     angolo=angolo).intersection(clip)
                parti[f"{c.nome}_attacco_{k}_{j}"] = porta
                canali = porta if canali is None else canali.union(porta)

    # --- fessura di film del corpo centrale ----------------------------------- #
    if film is not None:
        #: IL FILM E' UN VUOTO CHE SI APRE IN CAMERA, come i getti: non entra
        #: nei canali, che sono i vuoti che NON devono aprirsi da nessuna parte.
        parti["film"] = costruisci_film(grid, film, xp)

    # --- interfaccia di montaggio -------------------------------------------- #
    if montaggio is not None and testa is not None:
        pad, fori_m = costruisci_montaggio(grid, testa, montaggio, xp)
        corpo = corpo.union(pad)
        parti["montaggio"] = fori_m

    # --- collettore dell'acqua: un ingresso e una uscita ---------------------- #
    if collettore is not None:
        piu_solido, parti_col = costruisci_collettore(
            grid, gas, mantello, d, circuiti, testa, collettore, t_max, xp)
        if piu_solido is not None:
            corpo = corpo.union(piu_solido)
        for nome_p, campo_p in parti_col.items():
            parti[nome_p] = campo_p
            canali = campo_p if canali is None else canali.union(campo_p)

    # --- iniezione ----------------------------------------------------------- #
    # Questi vuoti SI' devono aprirsi in camera: e' la loro funzione. Restano
    # quindi fuori dal controllo degli spessori, in campi separati.
    if testa is None:
        fori = _fori_iniezione(grid, d, t_face, xp)
        parti["fori_iniezione"] = fori
    else:
        aria = cavita_aria(grid, testa, xp).difference(plug)
        gpl = cavita_gpl(grid, testa, xp)
        parti["aria"] = aria
        parti["gpl"] = gpl
        fori = aria.union(gpl)
        note_testa = testa.verifica(2.0 * (testa.R_anello - testa.R_getti))
    parti["gas"] = gas
    if "montaggio" in parti:
        fori = fori.union(parti["montaggio"])
    if "film" in parti:
        fori = fori.union(parti["film"])
    vuoti = fori if canali is None else canali.union(fori)

    solido = corpo.difference(vuoti)

    # Sacche di vuoto chiuse piu' piccole di un cubo da mezzo millimetro:
    # artefatti degli spigoli acuti, si riempiono. Le piu' grandi restano e
    # vengono segnalate, perche' quelle sono errori di progetto.
    from zefiro.sdf.printability import isole_di_materiale, riempi_sacche_chiuse
    solido, n_riempite, vol_riempito, n_rimaste = riempi_sacche_chiuse(
        solido, volume_massimo=(5.0e-4) ** 3)
    # e, simmetricamente, le isole di MATERIALE staccate: frammenti
    # sinterizzati che si staccano e girano nel circuito finche' non
    # ostruiscono un canale da 0.6 mm.
    solido, n_isole, vol_isole = isole_di_materiale(
        solido, volume_minimo=(5.0e-4) ** 3)

    note = list(note_passo)
    #: Nessun collettore deve sporgere oltre i piani che limitano il mantello
    #: (la faccia d'iniezione a x = 0 e il piano del labbro): oltre, viene
    #: tagliato e si apre sulla superficie del pezzo.
    for nome_col, x_col in collettori_x:
        lo_c = x_col - 0.5 * COLLETTORE_LARGHEZZA
        hi_c = x_col + 0.5 * COLLETTORE_LARGHEZZA
        if lo_c < 0.0 or hi_c > x_lip:
            note.append(
                f"il collettore '{nome_col}' occupa x da {lo_c*1e3:.2f} a "
                f"{hi_c*1e3:.2f} mm, fuori dal mantello (0 - {x_lip*1e3:.2f}): "
                "viene tagliato e apre il circuito sulla superficie")
    if testa is not None:
        note.extend(note_testa)
    if n_rimaste:
        note.append(
            f"{n_rimaste} sacche di vuoto CHIUSE sopra la soglia: polvere che non "
            "esce dal pezzo. Vanno collegate a un percorso di scarico o eliminate.")
    #: REGOLA RITIRATA, e vale la pena dire perche'. Finche' i due rami erano in
    #: PARALLELO, chiedere che la sezione degli attacchi superasse quella dei
    #: canali aveva senso: due rami che si contendono la stessa acqua vengono
    #: alimentati in proporzione alle loro resistenze, e un attacco stretto
    #: sposta portata da un ramo all'altro. In SERIE non c'e' niente da
    #: contendere: la stessa acqua passa da tutto, e un attacco stretto non
    #: cambia la ripartizione, aggiunge solo caduta. Il criterio giusto e'
    #: quindi il bilancio di pressione, che sta nel dimensionamento e non qui.
    if collettore is None:
        for c in circuiti:
            sezione_canali = c.n_canali * c.lato ** 2
            sezione_porte = PORTE_PER_COLLETTORE * math.pi / 4.0 * PORTA_DIAMETRO ** 2
            if sezione_porte < sezione_canali:
                note.append(
                    f"gli attacchi del ramo '{c.nome}' ({sezione_porte*1e6:.2f} mm2 in "
                    f"{PORTE_PER_COLLETTORE}) sono piu' stretti della somma dei suoi "
                    f"canali ({sezione_canali*1e6:.2f} mm2): sarebbero gli attacchi a "
                    "decidere la portata")
    if solido.tocca_il_bordo():
        note.append("il solido tocca il bordo della griglia: la mesh uscira' aperta")
    if n_isole:
        dg = getattr(solido, "diagnostica_isole", {})
        rimasta = dg.get("volume_massima_rimasta", 0.0)
        principale = dg.get("volume_principale", 0.0)
        note.append(
            f"{n_isole} isole di materiale staccate, {vol_isole*1e9:.3f} mm3 in "
            f"tutto; ne sono state TOLTE {dg.get('n_tolte', 0)} "
            f"({dg.get('volume_tolto', 0.0)*1e9:.3f} mm3), le altre restano. "
            f"Corpo principale {principale*1e6:.2f} cm3, "
            f"isola rimasta piu' grande {rimasta*1e9:.1f} mm3")
        if rimasta > 0.02 * principale:
            note.append(
                "un'isola vale piu' del 2 % del corpo: NON e' un frammento. O la "
                "geometria si e' spezzata davvero, o la griglia non la risolve. "
                "Rilanciare a passo dimezzato e confrontare prima di concludere.")
    if n_riempite:
        note.append(
            f"riempite {n_riempite} sacche minuscole ({vol_riempito*1e9:.3f} mm3 "
            f"in tutto, {vol_riempito/max(solido.volume(), 1e-30):.1e} del pezzo)")
    return MotoreSDF(grid=grid, solido=solido, cavita_gas=gas, canali=canali,
                     vuoti=vuoti, parti=parti, circuiti=tuple(circuiti),
                     testa=testa, collettore=collettore, montaggio=montaggio,
                     note=note)



# --------------------------------------------------------------------------- #
# Prova di tenuta: sonde e tappi
# --------------------------------------------------------------------------- #
def sonde_motore(motore: "MotoreSDF", d) -> list:
    """Dove mettere il gas traccia, un punto per ogni dominio fluido.

    Le coordinate NON sono scritte a mano: per ogni dominio si prende il punto
    piu' interno del suo campo, cioe' il centro del vuoto. Vedi
    `tenuta.punto_piu_interno` per il perche' una sonda scritta a mano e' un
    modo eccellente di far passare una prova sbagliata.
    """
    from zefiro.sdf.core import to_numpy
    from zefiro.sdf.tenuta import Sonda, punto_piu_interno

    g = motore.grid
    #: Dove una sonda ha diritto di stare: nel vuoto del pezzo FINITO e dentro
    #: il suo involucro. I campi delle singole parti sono costruiti per
    #: intersezioni e non sono campi di distanza validi lontano dalla loro
    #: superficie: senza questa maschera l'argomento del minimo puo' finire in
    #: aria libera. E' successo con il ramo di gola.
    #: L'involucro e' il solido PIU' tutte le sue cavita': la camera non sta
    #: in `vuoti` (non e' mai stata materiale da togliere, e' lo spazio fra
    #: mantello e plug), quindi usare solo solido+vuoti lasciava la sonda del
    #: gas senza nessun voxel ammesso.
    inv = motore.solido
    for campo in motore.parti.values():
        inv = inv.union(campo)
    dentro = (to_numpy(motore.solido.a) > 0.0) & (to_numpy(inv.a) < 0.0)
    del inv
    sonde = []
    nomi = [("gas", "gas"), ("aria", "aria"), ("gpl", "gpl")]
    for c in motore.circuiti:
        for strato in ("andata", "ritorno"):
            nomi.append((f"acqua_{c.nome}_{strato}", f"{c.nome}_{strato}"))
    #: Una sonda per ogni pezzo del collettore. Non e' ridondanza: se il
    #: raccordo non tocca uno dei due collettori, o se i fori d'uscita non
    #: arrivano nell'anello, il circuito si spezza in due domini e senza una
    #: sonda per parte la prova non se ne accorge - vedrebbe solo un dominio
    #: sano e uno che non ha nessuno dentro.
    #: I NOMI SI PRENDONO DALLE PARTI, non da una lista scritta a mano. La lista
    #: scritta a mano c'era, ed e' invecchiata al primo cambio di architettura:
    #: le sonde cercavano `acqua_anello` mentre le parti si chiamavano
    #: `acqua_anello_ingresso` e `acqua_anello_uscita`. Nessun errore, nessuna
    #: sonda: due pezzi del circuito uscivano dal controllo in silenzio.
    for chiave in sorted(k for k in motore.parti if k.startswith("acqua_")):
        nomi.append((chiave, chiave))
    for nome, chiave in nomi:
        campo = motore.parti.get(chiave)
        if campo is not None:
            sonde.append(Sonda(nome, punto_piu_interno(campo, g, dentro=dentro)))
    # L'esterno: uno spigolo della griglia, il piu' lontano possibile dal pezzo.
    sonde.append(Sonda("esterno", tuple(
        o + s_ * (n_ - 2) for o, s_, n_ in
        zip(g.origin, (g.spacing,) * 3, g.shape))))
    return sonde


def tappi(motore: "MotoreSDF", d, xp=np):
    """I tappi che si avvitano sugli attacchi prima di mettere in pressione.

    In collaudo non si prova un motore con i raccordi aperti: si tappa tutto e
    si guarda il manometro. Qui e' uguale, e serve a rendere la prova
    significativa: senza tappi ogni circuito comunica con l'esterno dal proprio
    attacco, quindi comunicano tutti fra loro attraverso l'esterno, e il
    controllo non distingue piu' un motore sano da uno bucato.

    I tappi sono cilindri di raggio generoso posti FUORI dal pezzo e affondati
    quel tanto che basta a chiudere la bocca di ciascun attacco.
    """
    from zefiro.sdf.clearances import _maschera_esterno  # noqa: F401  (documenta il verso)

    grid = motore.grid
    t = motore.testa
    campo = None

    def aggiungi(c):
        nonlocal campo
        campo = c if campo is None else campo.union(c)

    R_out = d["R_c"] + max(spessore_mantello(cc) for cc in motore.circuiti)
    col = motore.collettore
    if col is not None:
        #: CON IL COLLETTORE I TAPPI SONO DUE, non dodici: e' il senso di tutto
        #: il lavoro. Ciascuno chiude la bocca del proprio portagomma, e la
        #: chiude DOVE FINISCE IL GAMBO, non dove comincia: un tappo appoggiato
        #: sulla superficie del pezzo lascerebbe aperta tutta la canna del
        #: gambo, e la prova segnalerebbe una perdita che e' del tappo.
        for pg in (col.ingresso, col.uscita):
            if pg.verso == 0.0:
                q = pg.R_posizione + pg.lunghezza
                dirz = (0.0, math.cos(pg.angolo), math.sin(pg.angolo))
                base = (pg.x_base, q * math.cos(pg.angolo), q * math.sin(pg.angolo))
            else:
                dirz = (pg.verso, 0.0, 0.0)
                base = (pg.x_base + pg.verso * pg.lunghezza,
                        pg.R_posizione * math.cos(pg.angolo),
                        pg.R_posizione * math.sin(pg.angolo))
            #: il tappo affonda di un paio di voxel nel gambo e sporge di
            #: qualcuno in fuori: cosi' chiude la bocca qualunque sia il verso
            #: con cui marching cubes la trova.
            partenza = tuple(b - 2.0 * grid.spacing * v for b, v in zip(base, dirz))
            aggiungi(_cilindro_generico(
                grid, partenza, dirz,
                0.5 * (pg.d_canna + pg.interferenza) + 4.0 * grid.spacing,
                6.0 * grid.spacing, xp))
    # attacchi dell'acqua: radiali, a qualunque angolo -> si tappa tutta la
    # corona esterna con un guscio cilindrico chiuso
    for c in (motore.circuiti if col is None else ()):
        stazioni = [c.x_inizio, c.x_inizio + c.sfalsamento_ritorno] if c.ritorno \
            else [c.x_inizio, c.x_fine]
        for x in stazioni:
            #: Il tappo deve poggiare sulla superficie ESTERNA VERA a quella
            #: quota, che sul convergente e' molto piu' dentro che in camera:
            #: a x = 35.5 mm il pezzo finisce a 14.9 mm di raggio, non a 16.1.
            #: Con il raggio massimo il tappo galleggiava nel vuoto 0.7 mm
            #: fuori dal pezzo e la bocca dell'attacco restava aperta: la prova
            #: di tenuta segnalava una perdita che era del TAPPO, non del
            #: motore. Un falso allarme su un controllo di tenuta costa quanto
            #: un mancato allarme, perche' insegna a non fidarsi.
            r_loc = _raggio_parete(d, x) + spessore_mantello(c)
            #: Il tappo copre la bocca dell'attacco e NIENTE ALTRO. Con
            #: semilarghezza pari al diametro invece che al raggio, gli anelli
            #: diventavano manicotti lunghi che murava-no interi tratti di
            #: superficie esterna, e isolavano dall'"esterno" regioni che
            #: esterne lo sono: la prova segnalava 313 mm3 di "vuoto non
            #: rivendicato" che erano semplicemente aria attorno al motore.
            semi = 0.5 * PORTA_DIAMETRO + 2.0 * grid.spacing
            anello = cylinder(grid, R_out + 6.0 * grid.spacing,
                              x - semi, x + semi, xp)
            interno = cylinder(grid, r_loc - 2.0 * grid.spacing,
                               x - 2.0 * semi, x + 2.0 * semi, xp)
            aggiungi(anello.difference(interno))
    if t is not None:
        # attacco dell'aria: radiale, sulla testa
        x = t.x_porta_aria
        semi = 0.5 * t.d_porta_aria + 2.0 * grid.spacing
        anello = cylinder(grid, t.R_testa + 6.0 * grid.spacing,
                          x - semi, x + semi, xp)
        interno = cylinder(grid, t.R_plenum + 2.0 * grid.spacing,
                           x - 2.0 * semi, x + 2.0 * semi, xp)
        aggiungi(anello.difference(interno))
        # attacco del GPL: assiale, sulla faccia di monte. Il raggio da
        # tappare e' quello della SEDE FILETTATA (3.8 mm), non quello del
        # plenum a valle (2.1): con il secondo il tappo era piu' piccolo del
        # foro e lasciava una corona aperta tutt'intorno.
        aggiungi(cylinder(grid, 0.5 * t.d_filetto_gpl + 6.0 * grid.spacing,
                          t.x_faccia - 6.0 * grid.spacing,
                          t.x_faccia + 2.0 * grid.spacing, xp))
    return campo


def area_di_gola(motore: "MotoreSDF", x_gola: float, r_min: float, r_max: float,
                 n_stazioni: int = 5) -> tuple[float, float]:
    """Area di passaggio LIBERA misurata sul solido costruito, alla gola.

    PERCHE' ESISTE. Il raccordo fra mantello e corpo centrale, impostato a
    1.5 mm, gettava materiale dentro l'anello di gola: l'area libera usciva
    del **12 % piu' piccola** di quella di progetto, cioe' 12 % di spinta in
    meno. E non lo vedeva nessuno degli altri controlli - il pezzo era chiuso,
    in un solo blocco, con il circuito sigillato e un volume del tutto
    plausibile. Un raccordo e' un'operazione LOCALE solo se i due corpi sono
    lontani piu' del suo raggio; qui plug e labbro distano 1.78 mm e il
    raccordo li ha uniti attraverso il getto.

    L'area di gola e' la grandezza che fissa la portata e quindi la spinta:
    su un motore va misurata sul pezzo, non data per scontata dal disegno.

    Ritorna (area_misurata, area_teorica) in m^2.
    """
    import numpy as np

    from zefiro.sdf.core import to_numpy

    g = motore.grid
    campo = to_numpy(motore.solido.a)
    ay, az = np.array(g.axes()[1]), np.array(g.axes()[2])
    Y, Z = np.meshgrid(ay, az, indexing="ij")
    R = np.hypot(Y, Z)
    anello = (R > r_min) & (R < r_max)

    # piu' stazioni attorno alla gola: il minimo e' cio' che strozza davvero
    aree = []
    for k in range(n_stazioni):
        xq = x_gola - (n_stazioni - 1 - k) * g.spacing
        i = int(round((xq - g.origin[0]) / g.spacing))
        if not (0 <= i < g.shape[0]):
            continue
        aree.append(float((campo[i][anello] >= 0).sum()) * g.spacing**2)
    misurata = min(aree) if aree else 0.0
    teorica = math.pi * (r_max**2 - r_min**2)
    return misurata, teorica
