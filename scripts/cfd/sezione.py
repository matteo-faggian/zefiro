#!/usr/bin/env python3
"""Sezioni del motore con il campo di moto: che cosa fa davvero la geometria.

PERCHE' LA VORTICITA' E NON SOLO LA CONCENTRAZIONE. La frazione di massa dice
DOVE e' arrivato il combustibile; la vorticita' dice CHI ce l'ha portato. In un
getto trasversale il trasporto non lo fa la diffusione, lo fa una coppia di
vortici controrotanti che si forma appena il getto piega: e' quella struttura a
decidere se la miscela e' uniforme o se resta un nastro di combustibile
appiccicato a una parete. Guardare solo la concentrazione e' come giudicare un
motore dal rumore.

COME SI CALCOLA SENZA RIFARE IL CALCOLO. I piani di taglio salvano la velocita',
non la vorticita'. Su un piano meridiano - che contiene l'asse del motore -
serve pero' la sola componente AZIMUTALE,

    omega_theta = d(u_s)/dx - d(u_x)/ds

con `s` la coordinata radiale nel piano, e quella si ricava dalle sole
componenti nel piano. Non e' un ripiego: e' esattamente la componente che
descrive l'arrotolamento del getto nel piano che si sta guardando. Le altre due
non si vedrebbero comunque.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from filmato import coordinate_nel_piano, triangola  # noqa: E402
from sezioni import leggi_vtp  # noqa: E402


def _aree_doppie(tri):
    """2*area con segno di ogni triangolo (coordinate in mm)."""
    t = tri.triangles
    x, y = tri.x, tri.y
    return ((x[t[:, 1]] - x[t[:, 0]]) * (y[t[:, 2]] - y[t[:, 0]])
            - (x[t[:, 2]] - x[t[:, 0]]) * (y[t[:, 1]] - y[t[:, 0]]))


def gradiente_nodale(tri, f, a2):
    """Gradiente di `f` recuperato ai nodi, pesato sull'area dei triangoli.

    PERCHE' A MANO E NON `LinearTriInterpolator`. Le superfici di taglio
    campionate da OpenFOAM arrivano da 24 processori: hanno nodi duplicati sui
    bordi interni e qualche triangolo di area nulla. Il `TrapezoidMapTriFinder`
    di matplotlib le rifiuta ("Triangulation is invalid"), e ha ragione: e' un
    localizzatore di punti, non un derivatore. Ma per derivare non serve
    localizzare nulla. Su un triangolo a interpolazione lineare il gradiente e'
    costante e si scrive in forma chiusa dalle sole coordinate dei tre nodi:
    esattamente il gradiente di un elemento finito P1. I triangoli degeneri si
    scartano (non hanno un gradiente definito), gli altri si mediano sui nodi
    pesandoli con l'area - i triangoli grandi portano piu' informazione dei
    piccoli, ed e' la stessa media che fa un recupero alla Green-Gauss.
    """
    t = tri.triangles
    x, y = tri.x, tri.y
    buoni = np.abs(a2) > 1e-12 * np.median(np.abs(a2[np.abs(a2) > 0]))
    gx = np.zeros(len(t))
    gy = np.zeros(len(t))
    d = a2[buoni]
    f0, f1, f2 = f[t[buoni, 0]], f[t[buoni, 1]], f[t[buoni, 2]]
    x0, x1, x2 = x[t[buoni, 0]], x[t[buoni, 1]], x[t[buoni, 2]]
    y0, y1, y2 = y[t[buoni, 0]], y[t[buoni, 1]], y[t[buoni, 2]]
    gx[buoni] = ((y1 - y2) * f0 + (y2 - y0) * f1 + (y0 - y1) * f2) / d
    gy[buoni] = ((x2 - x1) * f0 + (x0 - x2) * f1 + (x1 - x0) * f2) / d
    peso = np.where(buoni, np.abs(a2), 0.0)
    sx = np.zeros(len(x))
    sy = np.zeros(len(x))
    sp = np.zeros(len(x))
    for k in range(3):
        np.add.at(sx, t[:, k], gx * peso)
        np.add.at(sy, t[:, k], gy * peso)
        np.add.at(sp, t[:, k], peso)
    vuoti = sp <= 0.0
    sp[vuoti] = 1.0
    sx[vuoti] = 0.0
    sy[vuoti] = 0.0
    return sx / sp, sy / sp


def campo_nel_piano(punti, poligoni, campi, normale, nome):
    """(triangolazione, valori) del campo richiesto sul piano di taglio."""
    import matplotlib.tri as mtri

    x, s = coordinate_nel_piano(punti, normale)
    tri = mtri.Triangulation(x * 1e3, s * 1e3, triangola(poligoni))
    a2 = _aree_doppie(tri)
    #: i triangoli degeneri (area nulla) non si mascherano per pignoleria: un
    #: triangolo di area nulla nel disegno e' una riga di colore arbitrario
    #: lunga quanto il triangolo, e nella derivata e' una divisione per zero.
    tri.set_mask(np.abs(a2) <= 1e-12 * np.median(np.abs(a2[np.abs(a2) > 0])))
    if nome == "vorticita":
        n = np.asarray(normale, dtype=float)
        n /= np.linalg.norm(n)
        if "vorticity" in campi:
            #: LA VORTICITA' DEL SOLUTORE, quando c'e'. Dalla corsa 5 il rotore
            #: viene calcolato sulle CELLE, con lo stesso schema con cui il
            #: solutore calcola tutto il resto, e poi tagliato sul piano.
            #: Ricostruirlo qui significherebbe derivare un campo gia'
            #: interpolato una volta sulla superficie di taglio, e derivare
            #: un'interpolazione amplifica il rumore dell'interpolazione.
            #:
            #: La componente che si vede su un piano meridiano e' quella lungo
            #: la NORMALE al piano: con e_s = n x e_x si ha e_x x e_s = n
            #: (perche' n e' perpendicolare all'asse), quindi
            #:     d(u_s)/dx - d(u_x)/ds = (rot U) . n
            #: cioe' esattamente la proiezione sulla normale. Non e' una
            #: approssimazione dell'altra formula: e' la stessa grandezza.
            return tri, campi["vorticity"] @ n
        #: RIPIEGO per le corse 1-4, che non avevano il campo: si ricostruisce
        #: dalle sole componenti nel piano.
        e = np.cross(n, np.array([1.0, 0.0, 0.0]))
        e /= np.linalg.norm(e)
        U = campi["U"]
        u_x = U[:, 0]
        u_s = U[:, 1:] @ e[1:]
        d_us_dx = gradiente_nodale(tri, u_s, a2)[0]
        d_ux_ds = gradiente_nodale(tri, u_x, a2)[1]
        #: le coordinate della triangolazione sono in MILLIMETRI, quindi i
        #: gradienti escono per millimetro: si riportano al secondo^-1.
        return tri, (d_us_dx - d_ux_ds) * 1e3
    v = campi[nome]
    if v.ndim > 1:
        v = np.linalg.norm(v, axis=1)
    return tri, v


#: (campo, mappa di colore, unita', simmetrico attorno allo zero)
PANNELLI = {
    "C3H8": ("inferno", "frazione di massa di GPL", False),
    "T": ("magma", "temperatura [K]", False),
    "U": ("viridis", "velocita' [m/s]", False),
    "vorticita": ("RdBu_r", "vorticita' azimutale [1/s]", True),
    #: k e nut dicono se il mescolamento che si vede e' turbolento o solo
    #: convettivo: nut/nu e' di quanto la turbolenza moltiplica la diffusione
    #: molecolare, ed e' il numero che decide se la correlazione di Holdeman -
    #: che e' una correlazione TURBOLENTA - sia applicabile a questo campo.
    "k": ("cividis", "energia cinetica turbolenta [m2/s2]", False),
    "nut": ("cividis", "viscosita' turbolenta [m2/s]", False),
}


def disegna(ax, tri, v, cmap, titolo, simmetrico, vmax=None):
    if vmax is None:
        vmax = float(np.percentile(np.abs(v) if simmetrico else v, 99.0))
    vmin = -vmax if simmetrico else 0.0
    liv = np.linspace(vmin, vmax, 41)
    c = ax.tricontourf(tri, np.clip(v, vmin, vmax), levels=liv, cmap=cmap,
                       extend="both")
    ax.set_aspect("equal")
    ax.set_title(titolo, fontsize=9)
    ax.set_xlabel("x [mm]", fontsize=8)
    ax.set_ylabel("r nel piano [mm]", fontsize=8)
    ax.tick_params(labelsize=7)
    return c


def main() -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from mesh_iniettore import quote_dal_progetto

    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", type=Path, required=True)
    ap.add_argument("--tempo", default=None, help="istante (default: l'ultimo)")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--campi", default="C3H8,vorticita,U")
    #: il campionamento di serie ("superfici") e' quello scritto DURANTE la
    #: corsa; "superfici2" e' il ricampionamento a posteriori con il piano di
    #: fianco spostato dentro al dominio (vedi _wsl/fianco.sh).
    ap.add_argument("--campionamento", default="superfici")
    a = ap.parse_args()

    d = quote_dal_progetto()
    ang = math.pi / d["n_getti"]
    piani = (("meridiano_getto", (0.0, -math.sin(ang), math.cos(ang)),
              "piano del getto"),
             ("meridiano_fianco", (0.0, 0.0, 1.0),
              "piano a meta' fra due getti"))
    base = a.caso / "postProcessing" / a.campionamento
    tempi = sorted((p for p in base.iterdir() if p.is_dir()),
                   key=lambda p: float(p.name))
    t = tempi[-1] if a.tempo is None else next(p for p in tempi if p.name == a.tempo)
    campi_voluti = a.campi.split(",")

    fig, assi = plt.subplots(len(campi_voluti), 2, figsize=(13, 3.1 * len(campi_voluti)),
                             dpi=130, squeeze=False)
    for r, nome_campo in enumerate(campi_voluti):
        cmap, unita, simm = PANNELLI[nome_campo]
        #: LA SCALA DI COLORE E' LA STESSA PER I DUE PIANI, altrimenti il
        #: confronto fra "dove c'e' il getto" e "dove non c'e'" e' truccato:
        #: due scale diverse fanno sembrare pieno anche il piano vuoto.
        dati = []
        for nome_piano, normale, _ in piani:
            f = t / f"{nome_piano}.vtp"
            punti, poligoni, campi = leggi_vtp(f)
            dati.append(campo_nel_piano(punti, poligoni, campi, normale, nome_campo))
        vmax = max(float(np.percentile(np.abs(v) if simm else v, 99.0))
                   for _, v in dati)
        for c, ((tri, v), (_, _, etichetta)) in enumerate(zip(dati, piani)):
            im = disegna(assi[r][c], tri, v, cmap,
                         f"{unita} - {etichetta}", simm, vmax)
            fig.colorbar(im, ax=assi[r][c], fraction=0.028, pad=0.02)
    fig.suptitle(
        f"Zefiro - sezione meridiana a t = {float(t.name)*1e6:.1f} us   "
        f"(getti {d['n_getti']} x {d['d_getto']*1e3:.3f} mm, C = {d['C']:.2f}, "
        f"aria {d['V_aria']:.0f} m/s, GPL {d['V_gpl']:.0f} m/s)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out = a.out or (a.caso / f"sezione_{t.name}.png")
    fig.savefig(out)
    print(f"scritto {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
