"""Il poligono meridiano e il suo volume di rivoluzione, in forma CHIUSA.

Sta in un modulo suo, senza OCCT, per una ragione precisa: il volume del pezzo
e' uno degli obiettivi dell'ottimizzatore, che lo valuta decine di migliaia di
volte. Costruire il solido con OCCT costa ~100 ms; la formula qui sotto costa
microsecondi ed e' ESATTA, non approssimata.

Esatta e non approssimata perche' il solido *e'* la rivoluzione di questo
poligono: non c'e' nessuna ipotesi di parete sottile da fare.

    V = 2 pi * (integrale di r su tutta l'area meridiana)

e con il teorema di Green, prendendo P = -r^2/2 e Q = 0,

    int int r dA = contorno di (-r^2/2) dx

che su un lato da (x0, r0) a (x1, r1) vale in forma chiusa

    int r^2 dx = (x1 - x0) (r0^2 + r0 r1 + r1^2) / 3

da cui  V = (pi/3) * | somma dei (x1 - x0)(r0^2 + r0 r1 + r1^2) |.

L'unica differenza rispetto al solido CAD sono i fori d'iniezione, che vengono
tagliati dopo: il test la misura invece di assumerla.
"""
from __future__ import annotations

import math
from typing import Mapping, Sequence


def meridian_polygon(
    derived: Mapping[str, float],
    contour_x: Sequence[float],
    contour_r: Sequence[float],
) -> list[tuple[float, float]]:
    """Poligono chiuso e semplice nel semipiano (x, r), in METRI.

    Percorso: faccia posteriore della piastra -> mantello esterno -> esterno
    del convergente -> faccia del labbro -> lato gas del convergente -> parete
    di camera -> faccia di iniezione -> corpo centrale -> contorno del plug ->
    chiusura sull'asse. Ne risulta UN SOLO solido connesso: piastra, mantello e
    plug sono lo stesso pezzo, che e' anche il modo in cui verra' stampato.
    """
    t = derived["t_wall"]
    tf = derived["t_face"]
    R_c = derived["R_c"]
    R_lip = derived["R_lip"]
    L_c = derived["L_c"]
    L_conv = derived["L_conv"]
    r_cb = derived["r_centerbody"]
    x_lip = L_c + L_conv
    x_plug_start = x_lip + contour_x[0]
    if x_plug_start <= 0.0:
        raise ValueError(
            "Il punto di gola del plug cade a monte della faccia di iniezione: "
            "camera troppo corta. Alza Lc_over_Dc o conv_half_angle."
        )

    pts: list[tuple[float, float]] = [
        (-tf, 0.0),
        (-tf, R_c + t),
        (L_c, R_c + t),
        (x_lip, R_lip + t),      # esterno del convergente (spessore RADIALE)
        (x_lip, R_lip),          # faccia del labbro, smussata di spessore t
        (L_c, R_c),              # lato gas del convergente, verso monte
        (0.0, R_c),              # parete di camera
        (0.0, r_cb),             # faccia di iniezione, verso l'asse
        (x_plug_start, r_cb),    # corpo centrale cilindrico
    ]
    # Si salta contour[0], che coincide con l'ultimo punto gia' inserito: un
    # punto duplicato genera uno spigolo di lunghezza nulla e OCCT lo rifiuta.
    pts += [(x_lip + xi, ri) for xi, ri in zip(contour_x[1:], contour_r[1:])]
    if contour_r[-1] > 1.0e-9:                  # plug troncato: faccia di base
        pts.append((x_lip + contour_x[-1], 0.0))
    return pts


def revolved_volume(polygon: Sequence[tuple[float, float]]) -> float:
    """Volume esatto del solido di rivoluzione del poligono chiuso [m^3].

    Il poligono e' chiuso implicitamente: l'ultimo punto si ricongiunge al
    primo. Il risultato e' preso in valore assoluto, cosi' non dipende dal
    verso di percorrenza.
    """
    n = len(polygon)
    if n < 3:
        raise ValueError(f"servono almeno 3 vertici, ce ne sono {n}")
    acc = 0.0
    for i in range(n):
        x0, r0 = polygon[i]
        x1, r1 = polygon[(i + 1) % n]
        acc += (x1 - x0) * (r0 * r0 + r0 * r1 + r1 * r1)
    return abs(math.pi / 3.0 * acc)


def revolved_lateral_area(polygon: Sequence[tuple[float, float]]) -> float:
    """Area totale della superficie di rivoluzione del contorno [m^2]."""
    n = len(polygon)
    total = 0.0
    for i in range(n):
        x0, r0 = polygon[i]
        x1, r1 = polygon[(i + 1) % n]
        total += math.pi * (r0 + r1) * math.hypot(x1 - x0, r1 - r0)
    return total


def wetted_area_from_contour(
    R_c: float, R_lip: float, L_c: float, L_conv: float, r_cb: float,
    contour_x: Sequence[float], contour_r: Sequence[float],
) -> float:
    """Superficie bagnata DAI GAS [m^2], per il teorema di Pappo-Guldino.

    E' un sottoinsieme del contorno meridiano, non tutto: il mantello esterno e
    la faccia posteriore della piastra non vedono i gas. Serve per il carico
    termico, dove conta solo cio' che scambia calore col gas caldo.

    Calcolata analiticamente dal profilo, non dal CAD: e' quindi una verifica
    INDIPENDENTE dalla tassellazione e dalle booleane.
    """
    x_lip = L_c + L_conv
    segs = [
        (0.0, R_c, L_c, R_c),                             # parete di camera
        (L_c, R_c, x_lip, R_lip),                         # convergente
        (0.0, R_c, 0.0, r_cb),                            # faccia di iniezione
        (0.0, r_cb, x_lip + contour_x[0], r_cb),          # corpo centrale
    ]
    segs += [(x_lip + contour_x[i], contour_r[i],
              x_lip + contour_x[i + 1], contour_r[i + 1])
             for i in range(len(contour_x) - 1)]
    return sum(math.pi * (r0 + r1) * math.hypot(x1 - x0, r1 - r0)
               for x0, r0, x1, r1 in segs)


#: Estensione del campo lontano, in multipli del RAGGIO DEL LABBRO.
#:
#: Il labbro e non la lunghezza del motore: la camera sta a monte e non ha
#: nulla a che vedere con quanto deve essere grande la regione di scarico.
#: La scala del getto e' il raggio di uscita, quindi e' quella che comanda.
#:
#: NON sono valori da manuale. Sono un punto di partenza dichiarato, e
#: l'indipendenza del risultato da questi due numeri va VERIFICATA con uno
#: studio: si allarga il dominio e si guarda se la spinta e la pressione di
#: base cambiano. Finche' quello studio non c'e', ogni risultato CFD porta
#: questa incertezza, e va detto invece che dimenticato.
FARFIELD_RADIAL_FACTOR = 8.0
FARFIELD_AXIAL_FACTOR = 10.0


def fluid_polygon(
    derived: Mapping[str, float],
    contour_x: Sequence[float],
    contour_r: Sequence[float],
    radial_factor: float = FARFIELD_RADIAL_FACTOR,
    axial_factor: float = FARFIELD_AXIAL_FACTOR,
) -> tuple[list[tuple[float, float]], list[str], list[tuple[int, int]]]:
    """Poligono meridiano del DOMINIO FLUIDO: punti, nome di ogni lato, e
    quali tratti sono CURVE LISCE invece che spezzate.

    Non e' il solido: e' il suo complemento dentro un contenitore che include
    la regione di scarico. L'aerospike espande all'ESTERNO, quindi il contorno
    del getto non e' una parete e il dominio deve arrivare fino all'ambiente.

    Ritorna `(punti, nomi, curve)`. `len(nomi) == len(punti)`: il lato `i` va
    dal vertice `i` al vertice `(i+1) % n` e si chiama `nomi[i]`. Il poligono e'
    chiuso implicitamente.

    `curve` elenca gli intervalli di vertici `[i, j]` che vanno costruiti come
    UNA spline e non come una spezzata. **Non e' cosmesi geometrica, e' la
    differenza fra una mesh utilizzabile e una da buttare.** Il contorno del
    plug arriva da `plug_contour` come ~140 punti su una decina di millimetri:
    se lo si da' a Gmsh come 140 segmenti, ogni vertice diventa un nodo
    obbligato e forza celle da 0.07 mm accanto a celle da 0.8 mm. Il risultato
    misurato erano schegge con non-ortogonalita' di 89 gradi, che OpenFOAM
    rifiuta. Come spline unica, Gmsh e' libero di scegliere i nodi secondo il
    campo di dimensione.

    Solo il plug e' liscio: gli altri tratti hanno spigoli VERI (la faccia di
    base, il labbro, gli angoli del campo lontano) e una spline li
    arrotonderebbe, cambiando la geometria invece che descriverla meglio.

    **Perche' i nomi si assegnano qui.** Dopo la rivoluzione si potrebbe
    riconoscere ogni superficie dalla posizione del suo baricentro, ed e' il
    modo in cui si fa di solito. E' anche il modo in cui si sbaglia: una
    condizione al contorno applicata alla faccia sbagliata non fa fallire
    niente e non si vede da nessuna parte, produce solo un risultato falso.
    Qui i nomi nascono in codice puro, senza gmsh, dove un test li puo'
    verificare uno per uno contro la geometria.

    **Frontiere canoniche assenti, e perche'.** `wall_throat` e `wall_cowl` non
    vengono prodotte: su un aerospike a espansione esterna con labbro di
    spessore nullo la gola e' delimitata dal plug (dentro) e dal solo SPIGOLO
    del labbro (fuori), che non e' una superficie. Inventarle come patch vuote
    darebbe a valle l'impressione di aver imposto una condizione che non
    esiste. Il flusso termico in gola si estrae dalla coordinata x su
    `wall_plug`: le patch servono alle condizioni al contorno, non al
    post-processing.
    """
    R_c = derived["R_c"]
    R_lip = derived["R_lip"]
    L_c = derived["L_c"]
    L_conv = derived["L_conv"]
    r_cb = derived["r_centerbody"]
    x_lip = L_c + L_conv

    x_throat = x_lip + contour_x[0]
    if x_throat <= 0.0:
        raise ValueError(
            "La gola cade a monte della faccia di iniezione: camera troppo corta."
        )
    x_base = x_lip + contour_x[-1]
    r_base = contour_r[-1]
    troncato = r_base > 1.0e-9

    R_far = radial_factor * R_lip
    # a valle: il piu' grande fra "tante volte il raggio di uscita" e "il
    # doppio del plug", cosi' un plug lungo non si ritrova il fondo addosso
    x_far = x_base + max(axial_factor * R_lip, 2.0 * (x_base - x_throat))
    if R_far <= R_c:
        raise ValueError(
            f"campo lontano R = {R_far*1e3:.1f} mm dentro la camera "
            f"(R_c = {R_c*1e3:.1f} mm): alza radial_factor"
        )

    pts: list[tuple[float, float]] = []
    nomi: list[str] = []
    curve: list[tuple[int, int]] = []

    def lato(punto: tuple[float, float], nome: str) -> None:
        """Aggiunge il vertice `punto` e dichiara che il lato che ne PARTE si
        chiama `nome`."""
        pts.append(punto)
        nomi.append(nome)

    lato((0.0, r_cb), "wall_faceplate")          # faccia anulare d'iniezione
    lato((0.0, R_c), "wall_chamber")
    lato((L_c, R_c), "wall_convergent")
    lato((x_lip, R_lip), "outlet_far")           # dal labbro verso l'ambiente
    lato((x_lip, R_far), "outlet_far")           # cilindro esterno
    lato((x_far, R_far), "outlet_far")           # sezione di uscita
    lato((x_far, 0.0), "axis")
    if troncato:
        lato((x_base, 0.0), "wall_plug")         # faccia di base, spigolo VERO
        inizio_curva = len(pts)
        lato((x_base, r_base), "wall_plug")
    else:
        inizio_curva = len(pts)
        lato((x_base, 0.0), "wall_plug")
    # contorno del plug a ritroso, dalla base verso la gola: e' liscio
    for i in range(len(contour_x) - 2, 0, -1):
        lato((x_lip + contour_x[i], contour_r[i]), "wall_plug")
    # il vertice di gola chiude la curva; il lato che ne parte e' il corpo
    # centrale, retto, che chiude il poligono tornando al primo vertice
    lato((x_throat, contour_r[0]), "wall_plug")
    curve.append((inizio_curva, len(pts) - 1))

    if len(pts) != len(nomi):                    # pragma: no cover - invariante
        raise AssertionError("un lato per vertice: invariante rotta")
    for i, j in curve:
        if not (0 <= i < j < len(pts)):          # pragma: no cover - invariante
            raise AssertionError(f"intervallo di curva assurdo: ({i}, {j})")
        if len(set(nomi[i:j])) != 1:             # pragma: no cover - invariante
            raise AssertionError(
                f"la curva ({i}, {j}) attraversa piu' frontiere: {set(nomi[i:j])}"
            )
    return pts, nomi, curve
