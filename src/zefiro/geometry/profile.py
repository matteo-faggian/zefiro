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
