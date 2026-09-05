"""Film cooling a fessura tangenziale, e il transitorio del corpo centrale.

PERCHE' QUESTO MODULO ESISTE. Il plug di un aerospike e' un corpo isolato in
mezzo alla camera: nessun circuito lo raggiunge, perche' portarci dentro
dell'acqua vorrebbe dire attraversare il getto. La sua sola difesa e' la
capacita' termica, e quella si esaurisce in due secondi.

Le alternative sono state provate e scartate con dei numeri, non a intuito:

  * RAFFREDDARLO CON UN GAS dall'interno. Con propano o aria in un foro da 2 mm
    si arriva a h = 3 kW/m2K, contro i 29 dell'acqua nello stesso foro. A 3
    MW/m2 il salto attraverso il film sarebbe di MILLE gradi. Un gas non
    raffredda una gola.
  * PORTARCI L'ACQUA. Alla stazione dei getti il corpo centrale ha 3.09 mm di
    raggio e dovrebbe contenere andata, ritorno, l'anello del GPL e tre pareti:
    servono 3.9 mm. Allargarlo a 4 mm ci farebbe stare tutto, ma il diametro
    dei getti va come 0.386 H e H come 1/R: i fori scenderebbero a 0.53 mm,
    sotto il minimo stampabile. L'iniettore si romperebbe per salvare il plug.
  * TRONCARE il plug. Massa e area calano insieme: da trunc 0.8 a 0.25 il tempo
    alla fusione passa da 2.03 a 2.26 s. Non e' una cura.
  * INGRANDIRE il motore. Da 50 a 200 N si va da 2.03 a 4.66 s. Nemmeno.

Resta il film: una fessura anulare alla base del plug che stende un velo di
combustibile fra il gas e la parete. E' la tecnica classica degli aerospike, e
qui il combustibile e' gia' dentro il corpo centrale perche' ci passa per
alimentare i getti.

CORRELAZIONE. Fessura tangenziale, Stollery & El-Ehwany (Int. J. Heat Mass
Transfer 8 (1965) 55-65):

    eta = 1                per xi <= 8      (film ancora integro)
    eta = 3.09 xi^-0.8     per xi > 8       (film ormai mescolato)
    xi  = (x / (M s)) Re_s^-0.25

con M = rho_c u_c / (rho_g u_g) rapporto dei flussi di massa, s altezza della
fessura, Re_s costruito sulla fessura. E' una correlazione da esperimenti su
lastra piana con gas simili: qui la parete e' conica, il gas e' molto piu' caldo
del refrigerante e il numero di Mach non e' quello degli esperimenti. Vale come
DIMENSIONAMENTO, non come verifica: la verifica e' una CFD o una prova.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

#: Soglia di xi oltre la quale il film comincia a degradare.
XI_CRITICO = 8.0
#: Costante della correlazione nel tratto degradato.
C_STOLLERY = 3.09


def efficacia_film(x: float, altezza_fessura: float, M: float, Re_s: float) -> float:
    """Efficacia adiabatica del film alla distanza `x` dalla fessura."""
    if x <= 0.0 or M <= 0.0 or altezza_fessura <= 0.0:
        return 1.0
    xi = (x / (M * altezza_fessura)) * Re_s ** -0.25
    if xi <= XI_CRITICO:
        return 1.0
    return min(1.0, C_STOLLERY * xi ** -0.8)


@dataclass(frozen=True)
class Film:
    frazione: float          # del combustibile totale
    altezza_fessura: float   # m
    mdot: float              # kg/s
    M: float                 # rapporto dei flussi di massa
    Re_s: float
    velocita: float          # m/s nella fessura
    T_refrigerante: float    # K
    #: Perche' il numero che segue NON e' una previsione. Un modello che
    #: restituisce un valore anche fuori dal proprio campo di validita' e' il
    #: modo piu' comodo di sbagliarsi: il numero c'e', ha le unita' giuste, e
    #: nessuno chiede da dove venga.
    avvertenze: tuple[str, ...] = ()


#: Reynolds di fessura sotto il quale la correlazione di Stollery & El-Ehwany
#: non e' piu' applicabile.
#:
#: La correlazione e' TURBOLENTA: e' ricavata su fessure tangenziali con
#: Re_s dell'ordine di 10^4 - 10^5, e la sua dipendenza x^-0.8 viene proprio
#: dallo strato limite turbolento. Sotto, il film resta laminare, si mescola
#: molto meno, e l'efficacia REALE puo' essere sia molto maggiore (il film
#: sopravvive piu' a lungo perche' non si mescola) sia molto minore (si stacca).
#: La correlazione non lo sa, e continua a restituire un numero.
RE_S_MINIMO_CORRELAZIONE = 1.0e4

#: Sotto questo rapporto di soffiaggio la fessura tangenziale e' fuori dal campo
#: in cui e' stata caratterizzata.
M_MINIMO_RAGIONEVOLE = 0.5


def dimensiona_film(mdot_film: float, altezza_fessura: float, raggio_fessura: float,
                    rho_c: float, mu_c: float, rho_g: float, u_g: float,
                    T_c: float, frazione: float) -> Film:
    area = 2.0 * math.pi * raggio_fessura * altezza_fessura
    u_c = mdot_film / (rho_c * area)
    Re_s = rho_c * u_c * altezza_fessura / mu_c
    avvertenze: list[str] = []
    if 0.0 < Re_s < RE_S_MINIMO_CORRELAZIONE:
        #: NOTA CHE VALE PIU' DEL NUMERO. Re_s NON dipende dall'altezza della
        #: fessura: Re_s = rho u s / mu con u = mdot / (rho 2 pi R s), quindi
        #: Re_s = mdot / (2 pi R mu). Stringere la fessura alza la velocita' e
        #: abbassa s nella stessa misura. L'unico modo di alzare Re_s e' mandare
        #: piu' portata, e la portata e' combustibile sottratto alla spinta.
        #: Un film laminare su questo motore non e' una scelta: e' una
        #: conseguenza della sua taglia.
        avvertenze.append(
            f"Re_s = {Re_s:.0f}, contro i {RE_S_MINIMO_CORRELAZIONE:.0e} sotto cui la "
            "correlazione di Stollery & El-Ehwany (turbolenta) non e' applicabile. "
            "E Re_s non dipende dall'altezza della fessura - vale "
            "mdot/(2 pi R mu) - quindi non si aggiusta stringendola: si alza solo "
            "mandando piu' combustibile al film. L'efficacia calcolata qui e' "
            "un'estrapolazione fuori dal campo di validita', non una previsione")
    M = (rho_c * u_c) / (rho_g * u_g)
    if 0.0 < M < M_MINIMO_RAGIONEVOLE:
        avvertenze.append(
            f"rapporto di soffiaggio M = {M:.2f}: le fessure "
            f"tangenziali lavorano tipicamente fra {M_MINIMO_RAGIONEVOLE:.1f} e 2. "
            "Sotto, il film si consuma in fretta e la lunghezza protetta e' corta")
    return Film(frazione=frazione, altezza_fessura=altezza_fessura, mdot=mdot_film,
                M=M, Re_s=Re_s,
                velocita=u_c, T_refrigerante=T_c, avvertenze=tuple(avvertenze))





@dataclass(frozen=True)
class TransitorioPlug:
    T_finale: float          # K raggiunta a fine raffica
    t_fusione: float         # s, inf se non fonde
    massa: float             # kg
    potenza_iniziale: float  # W
    potenza_finale: float    # W


def transitorio_plug(stazioni, massa: float, durata: float,
                     film: Film | None, cp: float, T_fusione: float,
                     emissivita: float = 0.3, T_iniziale: float = 300.0,
                     dt: float = 1.0e-3,
                     fattore_efficacia: float = 1.0) -> TransitorioPlug:
    """Integra la temperatura del plug trattato come corpo termicamente sottile.

    `stazioni` e' una lista di (x_dalla_fessura, area, h_gas, T_aw) prese dalla
    mappa di calore. Il flusso NON e' costante: dipende dalla temperatura del
    plug, che cresce, e il motore si autolimita mentre si scalda. Fissare q e
    dividere l'energia per la potenza (come faceva la prima versione di questo
    controllo) sottostima la durata.

    Con il film, T_aw diventa T_aw - eta (T_aw - T_c): la parete non vede piu'
    il gas, vede una miscela.

    `fattore_efficacia` moltiplica eta ed e' il modo onesto di usare una
    correlazione ricavata su lastra piana, con gas simili e senza il salto di
    temperatura che c'e' qui. Dimensionare a fattore 1 significa credere alla
    correlazione al punto da giocarci il pezzo; a 0.5 significa chiedere che il
    progetto regga anche se il film rende meta' di quanto promette. Costa
    combustibile, e il combustibile e' la risorsa piu' scarsa del motore: e'
    proprio per questo che il fattore va scritto e non lasciato implicito.
    """
    SIGMA = 5.670374419e-8
    area_tot = sum(a for _, a, _, _ in stazioni)
    T = T_iniziale
    t = 0.0
    q0 = None
    t_fus = float("inf")
    while t < durata:
        Q = 0.0
        for x, area, h, T_aw in stazioni:
            if film is not None:
                eta = fattore_efficacia * efficacia_film(
                    x, film.altezza_fessura, film.M, film.Re_s)
                T_aw = T_aw - eta * (T_aw - film.T_refrigerante)
            Q += h * area * (T_aw - T)
        Q -= emissivita * SIGMA * area_tot * (T ** 4 - T_iniziale ** 4)
        if q0 is None:
            q0 = Q
        T += Q / (massa * cp) * dt
        t += dt
        if T >= T_fusione and t_fus == float("inf"):
            t_fus = t
            break
    return TransitorioPlug(T_finale=T, t_fusione=t_fus, massa=massa,
                           potenza_iniziale=q0 or 0.0, potenza_finale=Q)
