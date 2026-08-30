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


@dataclass
class MotoreSDF:
    """Il pezzo completo, con i campi intermedi tenuti per la verifica."""
    grid: Grid
    solido: Field
    cavita_gas: Field
    canali: Field
    vuoti: Field          # canali + fori d'iniezione + attacchi: TUTTO il vuoto
    circuiti: tuple[CircuitoRaffreddamento, ...]
    #: I vuoti tenuti SEPARATI, per poter misurare gli spessori di parete
    #: a coppie. Un unico campo unito direbbe solo che il pezzo e' chiuso,
    #: non quanto materiale resta fra due cavita' che non devono toccarsi.
    parti: dict = dc_field(default_factory=dict)
    note: list[str] = dc_field(default_factory=list)


def contorno_gas_chiuso(d, contour_x, contour_r) -> list[tuple[float, float]]:
    """Poligono della CAVITA' GASSOSA interna, chiusa sul piano del labbro.

    Chiusa e non aperta perche' serve un campo di distanza da una regione, e
    una regione aperta non ne ha uno. La chiusura sul piano del labbro viene
    poi ritagliata via: il getto esterno non e' parte del pezzo.
    """
    R_c, R_lip, L_c, L_conv = d["R_c"], d["R_lip"], d["L_c"], d["L_conv"]
    r_cb = d["r_centerbody"]
    x_lip = L_c + L_conv
    return [
        (0.0, r_cb),          # faccia d'iniezione, dal corpo centrale in fuori
        (0.0, R_c),
        (L_c, R_c),           # parete di camera
        (x_lip, R_lip),       # convergente fino al labbro
        (x_lip, r_cb),        # chiusura sul piano del labbro
    ]


def contorno_plug(d, contour_x, contour_r) -> list[tuple[float, float]]:
    """Corpo centrale cilindrico piu' il contorno del plug, come solido pieno."""
    L_c, L_conv, r_cb = d["L_c"], d["L_conv"], d["r_centerbody"]
    x_lip = L_c + L_conv
    pts = [(0.0, 0.0), (0.0, r_cb)]
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


def costruisci(
    d,
    contour_x: Sequence[float],
    contour_r: Sequence[float],
    circuiti: Sequence[CircuitoRaffreddamento],
    passo_griglia: float,
    raccordo: float = 0.0,
    xp=np,
) -> MotoreSDF:
    """Assembla il motore. Ogni passaggio lascia il proprio campo, cosi' la
    verifica puo' misurare i pezzi e non solo il risultato."""
    R_c = d["R_c"]
    L_c, L_conv = d["L_c"], d["L_conv"]
    x_lip = L_c + L_conv
    x_base = x_lip + contour_x[-1]
    t_max = max(c.parete_calda + c.lato + c.parete_fredda for c in circuiti)
    t_face = max(2.0 * d["t_wall"], 2.0e-3)

    R_out = R_c + t_max
    lo = (-t_face - 3 * passo_griglia, -R_out, -R_out)
    hi = (x_base + 3 * passo_griglia, R_out, R_out)
    grid = Grid.bounding(lo, hi, passo_griglia, margin=4 * passo_griglia)

    gas = revolve_polygon(grid, contorno_gas_chiuso(d, contour_x, contour_r), xp)
    plug = revolve_polygon(grid, contorno_plug(d, contour_x, contour_r), xp)

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

    # --- piastra d'iniezione ------------------------------------------------ #
    piastra = cylinder(grid, R_out, -t_face, 0.0, xp)

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
            ch = ch.union(_collettore(grid, gas, mantello, c,
                                      c.x_inizio + c.sfalsamento_ritorno, xp,
                                      strato="ritorno"))
            ch = ch.union(_collettore(grid, gas, mantello, c, c.x_fine, xp,
                                      strato="entrambi"))
        else:
            for x_col in (c.x_inizio, c.x_fine):
                ch = ch.union(_collettore(grid, gas, mantello, c, x_col, xp,
                                          strato="andata"))
        canali = ch if canali is None else canali.union(ch)

    # --- attacchi dell'acqua ------------------------------------------------ #
    # Entrano in `canali` e non fra i vuoti generici, perche' devono
    # sottostare alla stessa verifica: un attacco che buca la parete sprizza
    # acqua nel gas esattamente come un canale. La prima versione li teneva
    # fuori dal controllo, e l'attacco di uscita del ramo gola bucava davvero.
    for i_c, c in enumerate(circuiti):
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

    # --- fori d'iniezione ---------------------------------------------------- #
    # Questi SI' devono aprirsi in camera: e' la loro funzione. Restano quindi
    # fuori dal controllo di tenuta, e in un campo separato per non confonderli.
    fori = _fori_iniezione(grid, d, t_face, xp)
    parti["fori_iniezione"] = fori
    parti["gas"] = gas
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
    if n_rimaste:
        note.append(
            f"{n_rimaste} sacche di vuoto CHIUSE sopra la soglia: polvere che non "
            "esce dal pezzo. Vanno collegate a un percorso di scarico o eliminate.")
    for c in circuiti:
        sezione_canali = c.n_canali * c.lato ** 2
        sezione_porte = PORTE_PER_COLLETTORE * math.pi / 4.0 * PORTA_DIAMETRO ** 2
        if sezione_porte < sezione_canali:
            note.append(
                f"gli attacchi del ramo '{c.nome}' ({sezione_porte*1e6:.2f} mm2 in "
                f"{PORTE_PER_COLLETTORE}) sono piu' stretti della somma dei suoi canali "
                f"({sezione_canali*1e6:.2f} mm2): sarebbero gli attacchi a decidere la portata")
    if solido.tocca_il_bordo():
        note.append("il solido tocca il bordo della griglia: la mesh uscira' aperta")
    if n_isole:
        note.append(
            f"tolte {n_isole} isole di materiale staccate ({vol_isole*1e9:.3f} mm3): "
            "frammenti che si staccherebbero nel circuito")
    if n_riempite:
        note.append(
            f"riempite {n_riempite} sacche minuscole ({vol_riempito*1e9:.3f} mm3 "
            f"in tutto, {vol_riempito/max(solido.volume(), 1e-30):.1e} del pezzo)")
    return MotoreSDF(grid=grid, solido=solido, cavita_gas=gas, canali=canali,
                     vuoti=vuoti, parti=parti, circuiti=tuple(circuiti), note=note)


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
