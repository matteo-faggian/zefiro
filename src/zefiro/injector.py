"""Iniettore a getti trasversali (jet in crossflow) per la coppia aria / GPL.

PERCHE' NON UN COASSIALE A TAGLIO
---------------------------------
Il motore ha tempi chimici di 30-100 us (zefiro.l0.chemistry, tempo di
spegnimento di un PSR) e un tempo di permanenza di alcuni millisecondi. Il
fattore mille fra i due dice che **la combustione non e' limitata dalla
chimica: e' limitata dal mescolamento**. Cio' che decide le prestazioni non e'
quanto in fretta brucia, e' quanto in fretta i due gas si trovano.

Per una coppia gas-gas la classifica sperimentale degli elementi e' nota e non
e' opinabile (NASA CR-121234, Calhoon, Ito, Kors, 1973, caratterizzazione di
iniettori gas-gas): premiscelato ed elementi a impatto danno efficienze di
mescolamento alte, il coassiale con swirl sta in mezzo, e **il coassiale a
taglio semplice e' l'ultimo**. Il motivo e' fisico: un coassiale affida il
mescolamento alla crescita di uno strato di taglio, che e' lenta e cresce
linearmente. Un getto trasversale genera invece una coppia di vortici
controrotanti che avvolge il flusso principale, e l'area d'interfaccia cresce
in modo molto piu' rapido.

Il premiscelato, che sarebbe il migliore, e' scartato: mettere propano e aria
gia' miscelati a monte di un iniettore non strozzato significa avere una
camera di combustione dentro il collettore. Vedi la nota in fondo al modulo.

LA REGOLA DI DIMENSIONAMENTO
----------------------------
Per getti multipli in un flusso trasversale confinato la penetrazione ottima e'
governata da un solo numero adimensionale (Holdeman, "Mixing of multiple jets
with a confined subsonic crossflow", Progress in Energy and Combustion Science
19 (1993) 31-70):

    C = (S / H) * sqrt(J)

dove S e' il passo fra i getti, H l'altezza del condotto attraversato e

    J = (rho_getto * V_getto^2) / (rho_flusso * V_flusso^2)

il rapporto dei flussi di quantita' di moto. Il valore ottimo, cioe' quello per
cui i getti arrivano al centro del condotto e i profili di miscela sono i piu'
uniformi, e':

    C = 2.5   per getti da UN SOLO lato
    C = 1.25  per due file opposte allineate
    C = 5.0   per due file opposte sfalsate

Sono valori pubblicati, non calibrati qui. Il senso fisico e' immediato: se
C e' piccolo i getti sono troppo fitti o troppo deboli e restano attaccati alla
parete; se e' grande sono troppo radi o troppo forti, attraversano il condotto e
vanno a sbattere dall'altra parte. In mezzo c'e' un solo valore.

COSA QUESTO MODULO CALCOLA E COSA NO
------------------------------------
Calcola: le aree di efflusso dalle portate e dai salti di pressione (efflusso
comprimibile isentropico, non incomprimibile: a rapporti di pressione di 0.87
la differenza sulla densita' e' del 10 % e va sull'area), J, il passo S, il
numero di getti e il loro diametro. Verifica che il diametro sia realizzabile.

NON calcola l'efficienza di mescolamento: quella richiede o una prova o una
CFD, ed e' il passo successivo. Qui si sceglie una geometria che la letteratura
dice essere quella giusta, e la si rende verificabile.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

R_UNIVERSALE = 8314.462618      # J/(kmol K)

#: Holdeman 1993, tabella dei valori ottimi di C = (S/H) sqrt(J).
C_OTTIMO = {"un_lato": 2.5, "opposti_allineati": 1.25, "opposti_sfalsati": 5.0}

#: Diametro minimo di un foro passante ottenibile in SLM su 316L senza
#: ripresa meccanica. E' un TODO da confermare con il service (vedi TODO 7 in
#: ROADMAP.md): finche' non e' misurato, il valore qui e' una SOGLIA DI
#: PROGETTO conservativa, non un dato. Sotto questa soglia il foro va forato o
#: elettroerroso dopo la stampa, il che e' possibile ma va deciso, non subito.
D_FORO_MINIMO_STAMPATO = 6.0e-4


@dataclass(frozen=True)
class Efflusso:
    """Stato di un gas che esce da un orifizio, efflusso isentropico."""
    area: float          # m2, area geometrica (Cd gia' incluso nel calcolo)
    velocita: float      # m/s
    densita: float       # kg/m3, alla sezione di uscita
    mach: float
    strozzato: bool
    T_uscita: float      # K
    flusso_qdm: float    # Pa, rho V^2


def efflusso_gas(mdot: float, p_monte: float, p_valle: float, T_monte: float,
                 massa_molare: float, gamma: float, cd: float) -> Efflusso:
    """Area di un orifizio per un gas, con espansione isentropica.

    Il modello incomprimibile (mdot = Cd A sqrt(2 rho Dp)) sbaglia qui: con
    p_valle/p_monte = 0.87 il gas si espande, la densita' all'uscita e' l'88 %
    di quella a monte e la velocita' e' quella corrispondente al salto
    entalpico, non a sqrt(2Dp/rho). Sull'area la differenza e' di alcuni punti
    percentuali, cioe' esattamente l'ordine di grandezza su cui si decide se un
    foro e' stampabile o no.

    Se il rapporto scende sotto quello critico l'orifizio si strozza e la
    portata non dipende piu' da p_valle: e' una condizione DESIDERABILE (isola
    la camera dall'alimentazione) ma qui non si verifica, perche' il GPL non ha
    abbastanza pressione di bombola per pagarsi un salto critico.
    """
    if p_valle >= p_monte:
        raise ValueError("p_valle >= p_monte: nessun efflusso.")
    R = R_UNIVERSALE / massa_molare
    rapporto_critico = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    r = p_valle / p_monte
    strozzato = r <= rapporto_critico
    if strozzato:
        r = rapporto_critico
    T = T_monte * r ** ((gamma - 1.0) / gamma)
    cp = gamma * R / (gamma - 1.0)
    V = math.sqrt(max(2.0 * cp * (T_monte - T), 0.0))
    rho = p_monte * r / (R * T)
    area = mdot / (cd * rho * V)
    return Efflusso(area=area, velocita=V, densita=rho,
                    mach=V / math.sqrt(gamma * R * T), strozzato=strozzato,
                    T_uscita=T, flusso_qdm=rho * V * V)


@dataclass(frozen=True)
class Iniettore:
    """Geometria completa dell'iniettore a getti trasversali."""
    # -- flusso trasversale (aria) --
    R_medio: float          # m, raggio ESTERNO dell'anello d'aria
    altezza_anello: float   # m, H: altezza radiale del condotto attraversato
    area_aria: float        # m2
    V_aria: float
    rho_aria: float
    # -- getti (GPL) --
    n_getti: int
    d_getto: float          # m
    passo: float            # m, S misurato sulla circonferenza di iniezione
    R_iniezione: float      # m, raggio da cui partono i getti
    V_getto: float
    rho_getto: float
    # -- adimensionali --
    J: float
    C: float
    C_obiettivo: float
    # -- verifiche --
    dp_aria: float
    dp_gpl: float
    lunghezza_mescolamento: float   # m, stima x/H ~ 1 (Holdeman)
    avvertenze: tuple[str, ...] = ()

    @property
    def stampabile(self) -> bool:
        return self.d_getto >= D_FORO_MINIMO_STAMPATO

    @property
    def scarto_C(self) -> float:
        return self.C / self.C_obiettivo - 1.0


#: Su un ANELLO le due file opposte non stanno sulla stessa circonferenza: la
#: esterna e' piu' lunga di R_e/R_i. Da qui segue un fatto geometrico che la
#: correlazione, ricavata su un condotto rettangolare, non poteva prevedere -
#: e che decide quali disposizioni sono realizzabili qui. Vedi
#: `disposizione_realizzabile`.
DISPOSIZIONI_A_DUE_FILE = ("opposti_allineati", "opposti_sfalsati")


def disposizione_realizzabile(disposizione: str, R_interno: float,
                              R_esterno: float, tolleranza: float = 0.10):
    """(realizzabile, motivo) per una disposizione su un condotto ANULARE.

    IL PUNTO GEOMETRICO. Holdeman studia un condotto rettangolare, dove le due
    pareti opposte sono lunghe uguale e due file opposte hanno per forza lo
    stesso passo. Su un anello no: la circonferenza esterna e' piu' lunga di
    quella interna del rapporto R_e/R_i, che qui vale 1.63. Allora:

      * ALLINEATI vuol dire stessi azimut sulle due file, quindi lo STESSO
        numero di getti, quindi passi diversi: S_e/S_i = R_e/R_i. E siccome
        C = (S/H) sqrt(J) con J e H comuni, le due file hanno C diversi dello
        stesso rapporto. **Non possono stare tutte e due sull'ottimo.** Una
        delle due sara' fuori del 63 %, cioe' sei volte la tolleranza con cui
        questo modulo accetta uno scarto.
      * SFALSATI non chiede gli stessi azimut, quindi ogni fila puo' avere il
        proprio numero di getti e stare sul proprio ottimo. **Realizzabile.**

    Non e' un dettaglio di implementazione: e' il motivo per cui su un iniettore
    anulare la disposizione opposta sensata e' quella sfalsata. Un condotto
    rettangolare non lo direbbe.
    """
    if disposizione not in DISPOSIZIONI_A_DUE_FILE:
        return True, ""
    rapporto = R_esterno / R_interno
    if disposizione == "opposti_sfalsati":
        return True, ""
    if abs(rapporto - 1.0) <= tolleranza:
        return True, ""
    return False, (
        f"su un anello le due file allineate hanno per forza lo stesso numero di "
        f"getti su circonferenze diverse ({R_interno*1e3:.2f} e {R_esterno*1e3:.2f} mm): "
        f"i loro C differiscono del {(rapporto-1.0)*100:.0f} %, contro una tolleranza "
        f"del {tolleranza*100:.0f} %. Su un anello la disposizione opposta che sta "
        "sull'ottimo su entrambe le file e' quella SFALSATA.")


def progetta(
    mdot_aria: float, mdot_gpl: float,
    p_c: float, p_aria_monte: float, p_gpl_monte: float,
    T_aria: float, T_gpl: float,
    MW_aria: float, MW_gpl: float, gamma_aria: float, gamma_gpl: float,
    cd: float,
    R_iniezione: float,
    n_getti: int | None = None,
    disposizione: str = "un_lato",
) -> Iniettore:
    """Dimensiona l'iniettore.

    L'ordine dei passi non e' arbitrario, ed e' il punto di tutto il modulo:

    1. le due velocita' NON si scelgono: le fissa il salto di pressione
       disponibile, che a sua volta lo fissa la bombola. Sono un dato.
    2. da velocita' e portate seguono le aree, quindi H (l'aria passa in una
       corona di raggio medio dato) e l'area totale dei getti.
    3. J segue dalle velocita' e dalle densita': **non e' un parametro libero**.
       Con salti di pressione uguali in frazione di p_c e densita' simili, J
       viene vicino a 1 da solo.
    4. il numero di getti e' l'unica cosa che resta libera, e la fissa la
       regola di Holdeman: n = 2 pi R / S con S = C H / sqrt(J).
    5. il diametro segue. Se e' troppo piccolo per la macchina, non si "aggiusta
       un parametro": si cambia il raggio di iniezione, perche' e' il raggio a
       decidere quanti getti stanno su una circonferenza.

    `n_getti` si puo' imporre per esplorare il compromesso fra diametro
    realizzabile e valore di C: e' il solo modo onesto di allontanarsi
    dall'ottimo, cioe' sapendo di quanto.
    """
    if disposizione not in C_OTTIMO:
        raise ValueError(f"disposizione: {sorted(C_OTTIMO)}")
    if disposizione in DISPOSIZIONI_A_DUE_FILE:
        #: SI RIFIUTA invece di costruire. Chiedendo `opposti_allineati` questa
        #: funzione restituiva, in silenzio, UNA fila sola di 8 getti con
        #: C = 1.25 - cioe' esattamente la geometria che per una fila sola
        #: Holdeman giudica pessima (getti troppo fitti, restano attaccati alla
        #: parete), etichettata con l'ottimo di una disposizione che non c'era.
        #: Non falliva: produceva numeri. Costruire una fila opposta per davvero
        #: chiede due cose che qui non ci sono - il raggio della seconda fila e
        #: la ripartizione della portata fra le due - e finche' non ci sono, il
        #: rifiuto e' l'unica risposta onesta. `spazio_di_scelta` dice quanto
        #: costerebbe averle.
        raise ValueError(
            f"disposizione '{disposizione}': questo modulo dimensiona UNA fila "
            "di getti su una circonferenza. Una disposizione opposta ha due file "
            "su due raggi diversi e una ripartizione di portata fra loro, che "
            "questa firma non prende. Usa `spazio_di_scelta` per sapere quanti "
            "getti e che diametro servirebbero.")
    aria = efflusso_gas(mdot_aria, p_aria_monte, p_c, T_aria, MW_aria, gamma_aria, cd)
    gpl = efflusso_gas(mdot_gpl, p_gpl_monte, p_c, T_gpl, MW_gpl, gamma_gpl, cd)

    # H dall'area anulare. ESATTA, non con A = 2 pi R H: quella e'
    # l'approssimazione di corona sottile, e qui H/R vale 0.4, dove sbaglia
    # l'area del 25 %. Su un'area sbagliata del 25 % si sbaglia la velocita'
    # dell'aria del 25 %, quindi J del 56 %, quindi il passo dei getti del 25 %.
    #     A = pi (R_e^2 - R_i^2)  ->  R_e = sqrt(R_i^2 + A/pi)
    R_esterno = math.sqrt(R_iniezione ** 2 + aria.area / math.pi)
    H = R_esterno - R_iniezione
    J = gpl.flusso_qdm / aria.flusso_qdm
    C_obj = C_OTTIMO[disposizione]

    circonferenza = 2.0 * math.pi * R_iniezione
    if n_getti is None:
        S = C_obj * H / math.sqrt(J)
        n = max(int(round(circonferenza / S)), 1)
    else:
        n = int(n_getti)
    S = circonferenza / n
    C = (S / H) * math.sqrt(J)
    d = math.sqrt(4.0 * gpl.area / (n * math.pi))

    avvertenze = []
    if n < 4:
        #: La correlazione di Holdeman e' ricavata su una FILA PERIODICA di
        #: getti. Con pochi getti su una circonferenza intera la periodicita'
        #: e' una finzione, e l'arrotondamento di n a intero sposta C fino a
        #: 1/(2n), cioe' oltre il 12 % gia' a n = 4. Il numero esce ancora, ma
        #: non significa piu' quello che dice di significare.
        avvertenze.append(
            f"solo {n} getti sull'intera circonferenza: la correlazione "
            "presuppone una fila periodica, qui l'ipotesi e' debole")
    if abs((S / H) * math.sqrt(J) / C_obj - 1.0) > 0.10:
        avvertenze.append(
            f"C = {C:.2f} contro {C_obj}: fuori dal 10 %, i getti "
            + ("attraversano il condotto" if C > C_obj else "restano attaccati alla parete"))
    if d < D_FORO_MINIMO_STAMPATO:
        avvertenze.append(
            f"foro da {d*1e3:.3f} mm sotto la soglia di stampabilita' "
            f"({D_FORO_MINIMO_STAMPATO*1e3:.2f} mm): va forato dopo la stampa")

    return Iniettore(
        avvertenze=tuple(avvertenze),
        R_medio=R_esterno, altezza_anello=H, area_aria=aria.area,
        V_aria=aria.velocita, rho_aria=aria.densita,
        n_getti=n, d_getto=d, passo=S, R_iniezione=R_iniezione,
        V_getto=gpl.velocita, rho_getto=gpl.densita,
        J=J, C=C, C_obiettivo=C_obj,
        dp_aria=p_aria_monte - p_c, dp_gpl=p_gpl_monte - p_c,
        #: Holdeman 1993: con spaziatura ottima i profili di miscela sono
        #: sostanzialmente uniformi entro un'altezza di condotto a valle.
        #: Si prende x/H = 2 per prudenza. E' una CORRELAZIONE, non un calcolo:
        #: la verifica vera e' una CFD o una prova.
        lunghezza_mescolamento=2.0 * H,
    )


# --------------------------------------------------------------------------- #
# Nota sul premiscelato, che sarebbe l'elemento migliore
# --------------------------------------------------------------------------- #
#
# CR-121234 mette il premiscelato in cima alla classifica: se i due gas sono
# gia' mescolati quando entrano, l'efficienza di mescolamento e' per
# definizione unitaria. E' scartato qui, e la ragione va scritta perche' e' una
# decisione, non una dimenticanza.
#
# Un collettore premiscelato e' un volume di propano-aria a phi 0.9 a 7 bar. Il
# ritorno di fiamma si ferma solo se il passaggio e' piu' stretto della distanza
# di spegnimento, che per propano-aria stechiometrico vale circa 2 mm a 1 bar e
# scala come l'inverso della pressione: a 6.4 bar sono circa 0.3 mm. Bisognerebbe
# quindi che OGNI passaggio del premiscelatore fosse sotto i 0.3 mm, cioe' sotto
# la soglia di stampabilita' dichiarata sopra. Non e' realizzabile con questa
# tecnologia.
#
# La transizione a detonazione e' invece improbabile (la dimensione di cella del
# propano-aria e' dell'ordine di decine di millimetri, molto piu' grande di
# qualunque passaggio qui dentro), quindi il rischio non e' la detonazione: e'
# una deflagrazione stabilizzata dentro il collettore, che fonde l'iniettore e
# risale verso la bombola. Basta e avanza.


# --------------------------------------------------------------------------- #
# Raccordi: la perdita di carico di un attacco NON e' trascurabile qui
# --------------------------------------------------------------------------- #
#: Filettature cilindriche BSPP (ISO 228-1). Per ciascuna:
#: (diametro esterno, passo, diametro del nocciolo, punta da maschiare,
#:  lunghezza di avvitamento tipica del maschio, diametro esterno della
#:  guarnizione a legare "bonded seal").
#: Sono quote di norma, non stime.
BSPP = {
    "G1/8": dict(d_esterno=9.728e-3, passo=0.907e-3, d_nocciolo=8.566e-3,
                 punta=8.80e-3, avvitamento=8.0e-3, d_guarnizione=13.9e-3),
    "G1/4": dict(d_esterno=13.157e-3, passo=1.337e-3, d_nocciolo=11.445e-3,
                 punta=11.80e-3, avvitamento=11.0e-3, d_guarnizione=17.9e-3),
    "G3/8": dict(d_esterno=16.662e-3, passo=1.337e-3, d_nocciolo=14.950e-3,
                 punta=15.25e-3, avvitamento=12.0e-3, d_guarnizione=21.2e-3),
    "G1/2": dict(d_esterno=20.955e-3, passo=1.814e-3, d_nocciolo=18.631e-3,
                 punta=19.00e-3, avvitamento=15.0e-3, d_guarnizione=26.7e-3),
    #: G3/4 e G1 aggiunte dopo che il sintetizzatore ha rifiutato un motore da
    #: 200 N: 134 g/s d'aria sono 15 litri al secondo, e nemmeno un G1/2 li
    #: porta senza mangiarsi il salto d'iniezione. La tabella si fermava dove
    #: si era fermato il progetto da 50 N, ed era il modello a essere corto,
    #: non la fisica.
    "G3/4": dict(d_esterno=26.441e-3, passo=1.814e-3, d_nocciolo=24.117e-3,
                 punta=24.50e-3, avvitamento=16.0e-3, d_guarnizione=32.4e-3),
    "G1":   dict(d_esterno=33.249e-3, passo=2.309e-3, d_nocciolo=30.291e-3,
                 punta=30.75e-3, avvitamento=18.0e-3, d_guarnizione=39.5e-3),
}

#: Sottomisura di stampa di un foro che verra' forato e maschiato dopo. Un foro
#: SLM esce sempre piu' piccolo del nominale e con la parete rugosa: si stampa
#: sotto e si porta a misura con la punta, che e' anche l'unico modo di avere
#: un foro tondo abbastanza da maschiarci dentro.
SOTTOMISURA_STAMPA = 1.2e-3


def perdita_raccordo(mdot: float, rho: float, d_passaggio: float,
                     k: float = 1.0) -> tuple[float, float]:
    """(velocita', perdita di carico) attraverso il passaggio di un raccordo.

    Perche' questo conto esiste. Il salto d'iniezione dell'aria vale 0.954 bar
    e **e' tutto il salto disponibile**: p_c e' stata fissata a p_sat/1.15, non
    c'e' un margine nascosto. Ogni pascal speso a monte dell'iniettore e' un
    pascal che l'iniettore non ha, e con esso se ne va la stabilita' che il
    salto d'iniezione serve a comprare.

    Un raccordo non e' un tubo: e' una strozzatura seguita da un allargamento
    brusco nel plenum, e l'allargamento brusco non recupera niente. La perdita
    e' quindi dell'ordine dell'intera pressione dinamica nel passaggio,
    Dp = k * rho V^2 / 2 con k vicino a 1.

    Il risultato e' brutale e vale la pena vederlo: 33.6 g/s d'aria a 8.7
    kg/m3 fanno 3.9 litri al secondo, e in un G1/4 (passaggio 7.5 mm) sono 87
    m/s. Un terzo del salto d'iniezione buttato nel raccordo.
    """
    area = math.pi / 4.0 * d_passaggio ** 2
    v = mdot / (rho * area)
    return v, k * 0.5 * rho * v * v


@dataclass(frozen=True)
class Alternativa:
    """Una disposizione di Holdeman valutata sulla geometria di un progetto."""
    disposizione: str
    C_obiettivo: float
    n_fila_interna: int
    n_fila_esterna: int          # 0 per le disposizioni a una fila
    d_getto: float               # m
    realizzabile_geometricamente: bool
    motivo: str                  # vuoto se realizzabile

    @property
    def n_totale(self) -> int:
        return self.n_fila_interna + self.n_fila_esterna

    def stampabile(self, d_min: float) -> bool:
        return self.d_getto >= d_min


def spazio_di_scelta(inj: "Iniettore", d_min: float = D_FORO_MINIMO_STAMPATO
                     ) -> list[Alternativa]:
    """Le tre disposizioni di Holdeman valutate sul progetto gia' fatto.

    A COSA SERVE. "Si puo' mescolare meglio?" non e' una domanda a cui si
    risponde ritoccando C: C sta gia' sul suo ottimo. Si risponde cambiando
    DISPOSIZIONE, e allora la domanda diventa quale disposizione sia
    raggiungibile - e la risposta la danno due vincoli che non si negoziano:

      1. LA GEOMETRIA ANULARE, che rende impossibile mettere due file allineate
         entrambe sull'ottimo (vedi `disposizione_realizzabile`);
      2. IL FORO MINIMO REALIZZABILE, che e' un tetto duro sul NUMERO di getti:
         l'area totale dei getti la fissa la portata diviso la velocita', e la
         velocita' la fissa il salto di pressione disponibile. Quell'area
         divisa per l'area del foro piu' piccolo che la macchina sa fare da'
         quanti getti si possono avere, punto.

    Il secondo vincolo e' quello che decide, ed e' appeso a un numero che
    NESSUNO HA MISURATO (vedi `D_FORO_MINIMO_STAMPATO`). Per questo la funzione
    prende `d_min` come argomento invece di leggerlo: il risultato cambia
    qualitativamente fra 0.60 e 0.35 mm, e chi guarda deve poterlo vedere.

    Le velocita', e quindi J, restano quelle del progetto: cambiare
    disposizione non cambia il salto di pressione disponibile.
    """
    H = inj.altezza_anello
    R_i = inj.R_iniezione
    R_e = R_i + H
    #: l'area totale dei getti non dipende dalla disposizione: e' portata
    #: diviso (densita' * velocita'), e tutte e due le fissa l'alimentazione.
    A_tot = math.pi / 4.0 * inj.d_getto ** 2 * inj.n_getti
    fuori = []
    for nome, C_obj in sorted(C_OTTIMO.items(), key=lambda kv: kv[1]):
        S = C_obj * H / math.sqrt(inj.J)
        due_file = nome in DISPOSIZIONI_A_DUE_FILE
        n_i = max(int(round(2.0 * math.pi * R_i / S)), 1)
        n_e = max(int(round(2.0 * math.pi * R_e / S)), 1) if due_file else 0
        n_tot = n_i + n_e
        ok, motivo = disposizione_realizzabile(nome, R_i, R_e)
        fuori.append(Alternativa(
            disposizione=nome, C_obiettivo=C_obj,
            n_fila_interna=n_i, n_fila_esterna=n_e,
            d_getto=math.sqrt(4.0 * A_tot / (n_tot * math.pi)),
            realizzabile_geometricamente=ok, motivo=motivo))
    return fuori


def getti_massimi(inj: "Iniettore", d_min: float = D_FORO_MINIMO_STAMPATO) -> int:
    """Quanti getti stanno nell'area disponibile, al foro minimo dato.

    E' il tetto duro di tutta l'architettura d'iniezione e vale la pena averlo
    in una riga: l'area totale dei getti e' fissata dall'alimentazione, quindi
    il numero di getti e' quell'area diviso l'area del foro piu' piccolo che si
    sa fare. Nessuna scelta di disposizione puo' superarlo.
    """
    A_tot = math.pi / 4.0 * inj.d_getto ** 2 * inj.n_getti
    return int(A_tot / (math.pi / 4.0 * d_min ** 2))
