"""Forme implicite: rivoluzione esatta di un poligono, e canali conformi.

La rivoluzione e' ESATTA, non approssimata, e vale la pena dire perche'.
Per un solido di rivoluzione il punto di superficie piu' vicino a P giace nel
semipiano meridiano che contiene P. Si vede in una riga: con P a (x_p, r_p, 0)
e Q a (x_q, r_q, theta),

    d^2 = (x_p - x_q)^2 + r_p^2 + r_q^2 - 2 r_p r_q cos(theta)

e con r_p, r_q >= 0 il minimo su theta e' in theta = 0, dove
d^2 = (x_p - x_q)^2 + (r_p - r_q)^2. Quindi **la distanza 3D dal solido di
rivoluzione e' la distanza 2D dal poligono nel piano (x, r)**, senza
approssimazioni. L'unico errore che resta e' quello di campionamento sulla
griglia, e i test lo misurano.
"""
from __future__ import annotations

import math
from typing import Sequence

import numpy as np

from zefiro.sdf.core import Field, Grid


def polygon_sdf_2d(px, pr, poligono: Sequence[tuple[float, float]], xp=np):
    """Distanza con segno da un poligono chiuso, nel piano (x, r).

    Negativa dentro. Il segno viene dal numero di attraversamenti di una
    semiretta (crossing number), non dall'orientamento del poligono: cosi' il
    risultato non dipende da come e' stato costruito il poligono, che e'
    esattamente il tipo di dipendenza nascosta da evitare.
    """
    P = xp.asarray(poligono, dtype=xp.float32)
    n = len(poligono)
    d2 = xp.full(px.shape, xp.float32(np.inf), dtype=xp.float32)
    dentro = xp.zeros(px.shape, dtype=bool)

    for i in range(n):
        ax, ar = P[i, 0], P[i, 1]
        bx, br = P[(i + 1) % n, 0], P[(i + 1) % n, 1]
        ex, er = bx - ax, br - ar
        wx, wr = px - ax, pr - ar
        ll = ex * ex + er * er
        t = xp.clip((wx * ex + wr * er) / (ll if ll > 0 else 1.0), 0.0, 1.0)
        cx, cr = wx - t * ex, wr - t * er
        d2 = xp.minimum(d2, cx * cx + cr * cr)

        # crossing number: il lato attraversa la semiretta orizzontale a destra?
        cond1 = (pr >= ar) & (pr < br)
        cond2 = (pr >= br) & (pr < ar)
        lato = ex * (pr - ar) - er * (px - ax)
        attraversa = ((cond1 & (lato > 0)) | (cond2 & (lato < 0)))
        dentro = dentro ^ attraversa

    d = xp.sqrt(d2)
    return xp.where(dentro, -d, d)


def revolve_polygon(grid: Grid, poligono: Sequence[tuple[float, float]],
                    xp=np) -> Field:
    """SDF esatto del solido ottenuto ruotando `poligono` attorno all'asse x."""
    X, Y, Z = grid.coords(xp)
    R = xp.sqrt(Y * Y + Z * Z)
    Xb, Rb = xp.broadcast_arrays(X.astype(xp.float32), R.astype(xp.float32))
    return Field(grid, polygon_sdf_2d(Xb, Rb, poligono, xp).astype(xp.float32), xp)


def cylinder(grid: Grid, raggio: float, x0: float, x1: float, xp=np) -> Field:
    """Cilindro sull'asse x. Serve soprattutto ai test: la sua distanza esatta
    si scrive a mano, quindi permette di misurare l'errore del campionamento."""
    X, Y, Z = grid.coords(xp)
    R = xp.sqrt(Y * Y + Z * Z)
    dr = R - raggio
    dx = xp.maximum(x0 - X, X - x1)
    dentro = xp.maximum(dr, dx)
    fuori = xp.sqrt(xp.maximum(dr, 0.0) ** 2 + xp.maximum(dx, 0.0) ** 2)
    a = xp.where(dentro > 0, fuori, dentro)
    a = xp.broadcast_to(a, grid.shape).astype(xp.float32).copy()
    return Field(grid, a, xp)


def sphere(grid: Grid, centro, raggio: float, xp=np) -> Field:
    X, Y, Z = grid.coords(xp)
    d = xp.sqrt((X - centro[0]) ** 2 + (Y - centro[1]) ** 2 + (Z - centro[2]) ** 2)
    a = xp.broadcast_to(d - raggio, grid.shape).astype(xp.float32).copy()
    return Field(grid, a, xp)


def helical_channels(
    grid: Grid,
    parete: Field,
    profondita: float,
    larghezza: float,
    altezza: float,
    n_canali: int,
    passo: float,
    x_inizio: float,
    x_fine: float,
    xp=np,
) -> Field:
    """Canali elicoidali CONFORMI alla parete, definiti dal campo stesso.

    Qui sta il motivo per cui si passa a geometria implicita. Un canale
    conforme in B-rep si costruisce spazzolando un profilo lungo una curva
    offset, e su una parete curva quella curva va calcolata, discretizzata e
    poi la booleana spesso fallisce. Con un campo di distanza il canale si
    scrive direttamente come CONDIZIONE sul campo:

        | parete - profondita | <= altezza/2      (a che distanza dalla parete)
        distanza_circonferenziale dal filo <= larghezza/2

    e la conformita' e' automatica, perche' `parete` E' la distanza dalla
    parete. Nessuna curva da costruire, nessuna booleana da far riuscire: il
    canale segue la parete perche' e' definito rispetto ad essa.

    `parete` deve essere il campo del VOLUME GASSOSO (negativo nel gas), cosi'
    `parete = profondita` e' la superficie a `profondita` dentro il materiale.

    NOTA ONESTA sul risultato: l'intersezione di due condizioni e' un massimo,
    e il massimo di due distanze non e' una distanza. Il campo che esce ha il
    SEGNO giusto ovunque (quindi la superficie a zero e' quella corretta e
    marching cubes la trova bene) ma sottostima la distanza vicino agli
    spigoli del canale. Va bene per generare e per fare booleane; NON va usato
    per un offset successivo di quel canale senza rigenerare il campo.
    """
    X, Y, Z = grid.coords(xp)
    theta = xp.arctan2(Z, Y)
    R = xp.sqrt(Y * Y + Z * Z)

    # fase elicoidale: a passo `passo` il canale compie un giro intero
    psi = theta - 2.0 * math.pi * X / passo
    # distanza angolare dal filo piu' vicino fra gli `n_canali` equispaziati
    passo_ang = 2.0 * math.pi / n_canali
    resto = psi - passo_ang * xp.round(psi / passo_ang)
    dist_circ = xp.abs(resto) * R

    radiale = xp.abs(parete.a - profondita) - 0.5 * altezza
    circonf = dist_circ - 0.5 * larghezza
    assiale = xp.maximum(x_inizio - X, X - x_fine)

    a = xp.maximum(xp.maximum(radiale, circonf), assiale)
    return Field(grid, xp.broadcast_to(a, grid.shape).astype(xp.float32).copy(), xp)
