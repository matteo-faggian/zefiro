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
from zefiro.sdf.shapes import cylinder, helical_channels, revolve_polygon


#: Larghezza assiale della gola anulare di collettore [m]. Presa pari a circa
#: due volte il lato del canale piu' grande: serve che la sezione di passaggio
#: del collettore sia molto maggiore di quella di un singolo canale, altrimenti
#: e' il collettore a strozzare e la distribuzione fra i canali non e' uniforme.
COLLETTORE_LARGHEZZA = 3.0e-3

#: Diametro dell'attacco radiale che porta l'acqua nel collettore [m].
#: Dimensionato perche' la sua sezione superi quella di tutti i canali del ramo
#: messi insieme: se strozzasse li', la portata la deciderebbe l'attacco e non
#: il progetto dei canali.
PORTA_DIAMETRO = 5.6e-3


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
    setto_ritorno: float = 0.4e-3   # m, pieno fra andata e ritorno


@dataclass
class MotoreSDF:
    """Il pezzo completo, con i campi intermedi tenuti per la verifica."""
    grid: Grid
    solido: Field
    cavita_gas: Field
    canali: Field
    vuoti: Field          # canali + fori d'iniezione + attacchi: TUTTO il vuoto
    circuiti: tuple[CircuitoRaffreddamento, ...]
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


def _collettore(grid, gas, mantello, c, x_centro: float, xp):
    """Gola anulare che raccoglie tutti i canali di un ramo, alla stessa
    profondita' dei canali. E' l'elemento che rende il circuito un circuito."""
    X, _, _ = grid.coords(xp)
    # Con la U il collettore deve unire ANCHE lo strato di ritorno, quindi si
    # estende in profondita' fino a coprirli entrambi.
    if c.ritorno:
        prof_rit = c.parete_calda + c.lato + c.setto_ritorno + 0.5 * c.lato
        centro = 0.5 * (c.parete_calda + 0.5 * c.lato + prof_rit)
        mezza = 0.5 * (prof_rit - (c.parete_calda + 0.5 * c.lato)) + 0.5 * c.lato
        profondita, semi = centro, mezza
    else:
        profondita, semi = c.parete_calda + 0.5 * c.lato, 0.5 * c.lato
    radiale = xp.abs(gas.a - profondita) - semi
    assiale = xp.abs(X - x_centro) - 0.5 * COLLETTORE_LARGHEZZA
    a = xp.maximum(radiale, assiale)
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
    for c in circuiti:
        ch = helical_channels(
            grid, gas,
            profondita=c.parete_calda + 0.5 * c.lato,
            larghezza=c.lato, altezza=c.lato,
            n_canali=c.n_canali, passo=c.passo_elica,
            x_inizio=c.x_inizio, x_fine=c.x_fine, xp=xp,
        )
        ch = ch.intersection(mantello)
        if c.ritorno:
            prof_rit = (c.parete_calda + c.lato + c.setto_ritorno + 0.5 * c.lato)
            rit = helical_channels(
                grid, gas, profondita=prof_rit,
                larghezza=c.lato, altezza=c.lato,
                n_canali=c.n_canali, passo=-c.passo_elica,
                x_inizio=c.x_inizio, x_fine=c.x_fine, xp=xp,
            ).intersection(mantello)
            ch = ch.union(rit)
        # collettori: due gole anulari che uniscono tutti i canali del ramo.
        # Senza, i canali non sono collegati a nulla e il pezzo non e' un
        # circuito ma una serie di buchi ciechi.
        for x_col in (c.x_inizio, c.x_fine):
            ch = ch.union(_collettore(grid, gas, mantello, c, x_col, xp))
        canali = ch if canali is None else canali.union(ch)

    # --- attacchi dell'acqua ------------------------------------------------ #
    # Entrano in `canali` e non fra i vuoti generici, perche' devono
    # sottostare alla stessa verifica: un attacco che buca la parete sprizza
    # acqua nel gas esattamente come un canale. La prima versione li teneva
    # fuori dal controllo, e l'attacco di uscita del ramo gola bucava davvero.
    for c in circuiti:
        prof_col = c.parete_calda + 0.5 * c.lato      # profondita' del collettore
        # con la U entrambi gli attacchi stanno a monte, sulla parte cilindrica
        stazioni = (c.x_inizio, c.x_inizio) if c.ritorno else (c.x_inizio, c.x_fine)
        for k, x_col in enumerate(stazioni):
            r_est = _raggio_parete(d, x_col) + t_max
            # si ferma appena dentro il bordo INTERNO del collettore, che sta a
            # (prof_col - lato/2) dalla parete calda. Andare oltre non serve a
            # niente e mangia la parete: sul ramo gola arrivava a 0.18 mm dal gas.
            prof_max = (c.parete_calda + c.lato + c.setto_ritorno + c.lato
                        if c.ritorno else prof_col + 0.5 * c.lato)
            penetrazione = t_max - prof_max + grid.spacing
            angolo = math.pi * k                     # i due attacchi opposti
            porta = _porta_acqua(grid, x_col, r_est, penetrazione,
                                 0.5 * PORTA_DIAMETRO, 3.0 * grid.spacing, xp,
                                 angolo=angolo)
            canali = porta if canali is None else canali.union(porta)

    # --- fori d'iniezione ---------------------------------------------------- #
    # Questi SI' devono aprirsi in camera: e' la loro funzione. Restano quindi
    # fuori dal controllo di tenuta, e in un campo separato per non confonderli.
    fori = _fori_iniezione(grid, d, t_face, xp)
    vuoti = fori if canali is None else canali.union(fori)

    solido = corpo.difference(vuoti)
    note = []
    for c in circuiti:
        sezione_canali = c.n_canali * c.lato ** 2
        sezione_porta = math.pi / 4.0 * PORTA_DIAMETRO ** 2
        if sezione_porta < sezione_canali:
            note.append(
                f"l'attacco del ramo '{c.nome}' ({sezione_porta*1e6:.2f} mm2) e' piu' "
                f"stretto della somma dei suoi canali ({sezione_canali*1e6:.2f} mm2): "
                "sarebbe l'attacco a decidere la portata")
    if solido.tocca_il_bordo():
        note.append("il solido tocca il bordo della griglia: la mesh uscira' aperta")
    return MotoreSDF(grid=grid, solido=solido, cavita_gas=gas, canali=canali,
                     vuoti=vuoti, circuiti=tuple(circuiti), note=note)
