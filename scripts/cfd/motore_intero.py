#!/usr/bin/env python3
"""Il motore INTERO in sezione meridiana, con dentro il campo CFD che esiste.

PERCHE' QUESTA FIGURA E NON UNO ZOOM PIU' LARGO. Il dominio CFD e' un settore
di 90 gradi che va da x = -6.6 mm a x = +12 mm: copre l'iniettore e il primo
quarto della camera, e basta. Il motore va da x = -1.5 mm a x = 49.1 mm, e la
parte che decide la spinta - convergente, gola, plug - **non e' calcolata**.
Disegnare solo il pezzo calcolato, con gli assi che finiscono a 12 mm, lascia
credere che il motore finisca li'. Disegnare tutto il motore e dire dove si
ferma il calcolo e' l'unica versione che non mente.

COME SI LEGGE. La sezione e' specchiata attorno all'asse:

    semipiano SUPERIORE  ->  piano che CONTIENE il getto
    semipiano INFERIORE  ->  piano a meta' strada fra due getti

Non sono due meta' dello stesso taglio: sono due tagli diversi disegnati
insieme, perche' la differenza fra i due E' la risposta alla domanda "i getti
si sono fusi in un anello?". Se i due semipiani sono uguali, si'. Se sopra c'e'
il combustibile e sotto no, i getti sono ancora quattro nastri separati.

Il metallo e' grigio e viene dal poligono meridiano ESATTO del sintetizzatore
(`zefiro.geometry.profile.meridian_polygon`), non da un disegno a mano: e' lo
stesso poligono che genera il solido stampato, quindi se la geometria cambia
questa figura cambia con lei.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sezione import PANNELLI, campo_nel_piano  # noqa: E402
from sezioni import leggi_vtp  # noqa: E402


def profilo_del_motore():
    """(poligono del metallo in mm, quote notevoli in mm) dal sintetizzatore."""
    from zefiro.geometry.profile import meridian_polygon
    from zefiro.sintesi import Impianto, Requisiti, progetta

    p = progetta(Requisiti(
        spinta=50.0, durata=5.0,
        impianto=Impianto(combustibile={"C3H8": 1.0}, T_bombola_min=288.15,
                          T_ossidante_iniezione=293.0,
                          T_combustibile_iniezione=283.0,
                          cd_ossidante=0.75, cd_combustibile=0.75)),
        geometria=False)
    par = p.geometria.params
    d = dict(par.derived)
    poly = meridian_polygon(d, par.plug_contour_x, par.plug_contour_r)
    x_lip = d["L_c"] + d["L_conv"]
    quote = {
        "R_c": d["R_c"] * 1e3,
        "r_cb": d["r_centerbody"] * 1e3,
        "R_lip": d["R_lip"] * 1e3,
        "L_c": d["L_c"] * 1e3,
        "x_lip": x_lip * 1e3,
        "x_gola": (x_lip + par.plug_contour_x[0]) * 1e3,
        "x_base": (x_lip + par.plug_contour_x[-1]) * 1e3,
        "t_face": d["t_face"] * 1e3,
        "t_wall": d["t_wall"] * 1e3,
    }
    return np.array(poly) * 1e3, quote


def disegna_riga(ax, dati, cmap, simm, vmax, metallo, quote, x_cfd):
    """Un campo: sopra il piano del getto, sotto il piano fra due getti."""
    vmin = -vmax if simm else 0.0
    liv = np.linspace(vmin, vmax, 41)
    im = None
    for (tri, v), segno in zip(dati, (+1.0, -1.0)):
        if segno < 0:
            #: si specchia la TRIANGOLAZIONE, non i dati: rovesciare i valori
            #: e tenere le coordinate darebbe un campo giusto su una geometria
            #: sbagliata, che e' il modo piu' elegante di dire una bugia.
            t2 = tri.__class__(tri.x, -tri.y, tri.triangles)
            t2.set_mask(tri.mask)
        else:
            t2 = tri
        im = ax.tricontourf(t2, np.clip(v, vmin, vmax), levels=liv,
                            cmap=cmap, extend="both", zorder=2)

    for segno in (+1.0, -1.0):
        ax.fill(metallo[:, 0], segno * metallo[:, 1], facecolor="0.72",
                edgecolor="0.25", linewidth=0.7, zorder=3)

    #: la frontiera del calcolo si DISEGNA. Un lettore che non la vede crede
    #: che il bianco a valle sia un risultato invece che un vuoto.
    ax.axvline(x_cfd, color="0.15", linestyle="--", linewidth=1.0, zorder=4)
    ax.text(x_cfd + 0.8, 0.0, "fine del dominio CFD\n(a valle: nessun calcolo)",
            fontsize=6.5, va="center", ha="left", zorder=5,
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.6", lw=0.5))
    ax.axhline(0.0, color="0.35", linestyle="-.", linewidth=0.6, zorder=4)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]", fontsize=8)
    ax.set_ylabel("sopra: piano del getto\nsotto: fra due getti", fontsize=7)
    ax.tick_params(labelsize=7)
    return im


def annota(ax, quote, y):
    for x, testo in ((0.0, "faccia\niniezione"),
                     (quote["L_c"], "fine camera"),
                     (quote["x_gola"], "gola"),
                     (quote["x_base"], "base plug")):
        ax.annotate(testo, xy=(x, y), xytext=(x, y * 1.30), fontsize=6.0,
                    ha="center", va="bottom", zorder=6,
                    arrowprops=dict(arrowstyle="-", lw=0.6, color="0.3"))


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", type=Path, required=True)
    ap.add_argument("--campionamento", default="superfici")
    ap.add_argument("--tempo", default=None)
    ap.add_argument("--campi", default="C3H8,vorticita,U")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    from mesh_iniettore import quote_dal_progetto

    metallo, quote = profilo_del_motore()
    q = quote_dal_progetto()
    ang = math.pi / q["n_getti"]
    piani = (("meridiano_getto", (0.0, -math.sin(ang), math.cos(ang))),
             ("meridiano_fianco", (0.0, 0.0, 1.0)))

    base = a.caso / "postProcessing" / a.campionamento
    tempi = sorted((p for p in base.iterdir() if p.is_dir()),
                   key=lambda p: float(p.name))
    t = tempi[-1] if a.tempo is None else next(p for p in tempi
                                               if p.name == a.tempo)
    letti = [leggi_vtp(t / f"{n}.vtp") for n, _ in piani]
    campi_voluti = a.campi.split(",")

    n = len(campi_voluti)
    fig, assi = plt.subplots(n, 1, figsize=(13.5, 3.05 * n), dpi=140,
                             squeeze=False)
    x_cfd = max(float(np.max(l[0][:, 0])) for l in letti) * 1e3
    R = quote["R_c"] + quote["t_wall"]
    for r, nome in enumerate(campi_voluti):
        cmap, unita, simm = PANNELLI[nome]
        dati = [campo_nel_piano(*letti[i], piani[i][1], nome)
                for i in range(2)]
        #: scala unica per i due semipiani: due scale diverse farebbero
        #: sembrare pieno anche il semipiano vuoto.
        vmax = max(float(np.percentile(np.abs(v) if simm else v, 99.0))
                   for _, v in dati)
        ax = assi[r][0]
        im = disegna_riga(ax, dati, cmap, simm, vmax, metallo, quote, x_cfd)
        ax.set_xlim(-3.0, quote["x_base"] + 3.0)
        ax.set_ylim(-1.32 * R, 1.32 * R)
        ax.set_title(unita, fontsize=9)
        fig.colorbar(im, ax=ax, fraction=0.016, pad=0.008)
        if r == 0:
            annota(ax, quote, R)

    fig.suptitle(
        f"Zefiro - motore intero in sezione, t = {float(t.name)*1e6:.1f} us   "
        f"(getti {q['n_getti']} x {q['d_getto']*1e3:.3f} mm, C = {q['C']:.2f}, "
        f"J = {q['J']:.3f})   -   il calcolo copre solo x < {x_cfd:.1f} mm",
        fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(a.out)
    print(f"scritto {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
