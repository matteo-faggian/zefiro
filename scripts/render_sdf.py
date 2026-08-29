#!/usr/bin/env python3
"""Guarda il pezzo: sezioni dal campo e vista 3D ombreggiata.

Le sezioni si tagliano dal CAMPO e non dalla mesh: un canale che non e' stato
scavato, o una parete che si e' chiusa dove doveva restare aperta, in una vista
esterna non si vede, in una sezione si.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


def rasterizza(v, f, larghezza, altezza, direzione, luce=(0.4, 0.5, 0.75)):
    """Z-buffer con ombreggiatura piatta. Nessuna libreria 3D: proiezione,
    ordinamento per profondita', e un array."""
    d = np.asarray(direzione, dtype=float)
    d /= np.linalg.norm(d)
    su = np.array([0.0, 0.0, 1.0])
    if abs(d @ su) > 0.95:
        su = np.array([0.0, 1.0, 0.0])
    ex = np.cross(su, d)
    ex /= np.linalg.norm(ex)
    ey = np.cross(d, ex)

    c = v.mean(axis=0)
    p = v - c
    u, w, z = p @ ex, p @ ey, p @ d
    sc = 0.9 * min(larghezza / (u.max() - u.min()), altezza / (w.max() - w.min()))
    px = ((u - u.min()) * sc + 0.5 * (larghezza - (u.max() - u.min()) * sc))
    py = ((w - w.min()) * sc + 0.5 * (altezza - (w.max() - w.min()) * sc))

    a, b, cc = f[:, 0], f[:, 1], f[:, 2]
    n = np.cross(v[b] - v[a], v[cc] - v[a])
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = n / np.where(ln > 0, ln, 1.0)
    L = np.asarray(luce, dtype=float)
    L /= np.linalg.norm(L)
    intensita = np.clip(0.22 + 0.78 * np.abs(n @ L), 0.0, 1.0)

    prof = (z[a] + z[b] + z[cc]) / 3.0
    ordine = np.argsort(-prof)              # dal piu' lontano al piu' vicino
    img = np.zeros((altezza, larghezza), dtype=float)
    X0 = np.stack([px[a], px[b], px[cc]], axis=1)
    Y0 = np.stack([py[a], py[b], py[cc]], axis=1)
    for i in ordine:
        x0, x1 = int(max(X0[i].min(), 0)), int(min(X0[i].max() + 1, larghezza))
        y0, y1 = int(max(Y0[i].min(), 0)), int(min(Y0[i].max() + 1, altezza))
        if x1 <= x0 or y1 <= y0:
            continue
        yy, xx = np.mgrid[y0:y1, x0:x1]
        ax, ay = X0[i, 0], Y0[i, 0]
        bx, by = X0[i, 1], Y0[i, 1]
        cx2, cy2 = X0[i, 2], Y0[i, 2]
        den = (by - cy2) * (ax - cx2) + (cx2 - bx) * (ay - cy2)
        if abs(den) < 1e-12:
            continue
        l1 = ((by - cy2) * (xx - cx2) + (cx2 - bx) * (yy - cy2)) / den
        l2 = ((cy2 - ay) * (xx - cx2) + (ax - cx2) * (yy - cy2)) / den
        m = (l1 >= -1e-9) & (l2 >= -1e-9) & (l1 + l2 <= 1 + 1e-9)
        if m.any():
            img[yy[m], xx[m]] = intensita[i]
    return img


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--passo", type=float, default=1.5e-4)
    ap.add_argument("--out", type=Path, default=Path("runs/motore50/vista.png"))
    a = ap.parse_args()

    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from build_50N import CAMERA, GOLA, punto_operativo
    from zefiro.sdf.core import to_numpy
    from zefiro.sdf.engine import costruisci
    from zefiro.sdf.meshing import isosurface

    op, params, l0 = punto_operativo()
    m = costruisci(params.derived, params.plug_contour_x, params.plug_contour_r,
                   [CAMERA, GOLA], a.passo, raccordo=1.5e-3)
    campo = to_numpy(m.solido.a)
    g = m.grid
    print(f"griglia {g.shape}, {g.n_voxels/1e6:.1f} Mvoxel")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(17.5, 10.0))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0], hspace=0.22, wspace=0.18)

    ax = fig.add_subplot(gs[0, :2])
    kz = g.shape[2] // 2
    sez = campo[:, :, kz].T
    ex = [g.origin[0]*1e3, (g.origin[0]+g.spacing*(g.shape[0]-1))*1e3,
          g.origin[1]*1e3, (g.origin[1]+g.spacing*(g.shape[1]-1))*1e3]
    ax.imshow((sez < 0).astype(float), origin="lower", extent=ex,
              cmap="Blues", vmin=0, vmax=1.6, aspect="equal", interpolation="nearest")
    ax.contour(np.linspace(ex[0], ex[1], g.shape[0]),
               np.linspace(ex[2], ex[3], g.shape[1]), sez, levels=[0.0],
               colors="#1a365d", linewidths=0.7)
    ax.set_title("Sezione meridiana: in blu il materiale, i vuoti sono i canali "
                 "d'acqua", fontsize=11)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")

    for k, (xq, tit) in enumerate(((0.016, "meta' camera"), (0.0365, "gola"))):
        ax = fig.add_subplot(gs[0, 2] if k == 0 else gs[1, 2])
        i = int(round((xq - g.origin[0]) / g.spacing))
        i = max(0, min(i, g.shape[0]-1))
        s = campo[i].T
        exy = [g.origin[1]*1e3, (g.origin[1]+g.spacing*(g.shape[1]-1))*1e3,
               g.origin[2]*1e3, (g.origin[2]+g.spacing*(g.shape[2]-1))*1e3]
        ax.imshow((s < 0).astype(float), origin="lower", extent=exy, cmap="Blues",
                  vmin=0, vmax=1.6, aspect="equal", interpolation="nearest")
        ax.set_title(f"Sezione a x = {xq*1e3:.1f} mm ({tit})", fontsize=10)
        ax.set_xlabel("y [mm]")
        ax.set_ylabel("z [mm]")

    print("marching cubes per la vista 3D...", flush=True)
    v, f = isosurface(m.solido)
    print(f"  {len(f)} triangoli")
    vc, fc = isosurface(m.canali)
    for k, (mesh, dirz, tit, cm) in enumerate((
            ((v, f), (0.55, 0.5, 0.67), "del pezzo, tre quarti", "bone"),
            ((vc, fc), (0.55, 0.5, 0.67), "della SOLA rete d'acqua", "copper"))):
        ax = fig.add_subplot(gs[1, k])
        img = rasterizza(mesh[0], mesh[1], 660, 540, dirz)
        ax.imshow(img, cmap=cm, origin="lower", vmin=0, vmax=1)
        ax.set_title(f"Vista {tit}", fontsize=10)
        ax.axis("off")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=125, bbox_inches="tight")
    print(f"scritto {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
