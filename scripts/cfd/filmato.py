#!/usr/bin/env python3
"""Il film: i getti che escono, si mescolano e vengono espulsi.

Due piani meridiani, non uno. Quello che CONTIENE l'asse del getto fa vedere
la penetrazione; quello a meta' strada fra due getti (il fianco del settore)
fa vedere se il combustibile arriva anche dove non e' stato sparato. Un
iniettore puo' sembrare ottimo sul primo e pessimo sul secondo: e' esattamente
il modo in cui un getto che non penetra inganna chi guarda una sola sezione.

I fotogrammi si disegnano con matplotlib su triangolazione (tricontourf) e non
con un renderer 3D: i dati sono un taglio piano, e un contorno su
triangolazione e' la rappresentazione esatta di quel taglio, senza
interpolazioni di comodo.
"""
from __future__ import annotations

import argparse
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sezioni import leggi_vtp  # noqa: E402


def triangola(poligoni):
    """Ventaglio: un poligono di n vertici diventa n-2 triangoli."""
    t = []
    for p in poligoni:
        for j in range(1, len(p) - 1):
            t.append((p[0], p[j], p[j + 1]))
    return np.asarray(t, dtype=np.int64)


def coordinate_nel_piano(punti, normale):
    """(x, s) con s = distanza dall'asse nel piano di taglio, con segno.

    Il taglio e' un piano che contiene l'asse x. La seconda coordinata utile e'
    la proiezione sul versore radiale che sta nel piano: usare y o z crudi
    darebbe un disegno schiacciato e falso quando il piano e' inclinato.
    """
    n = np.asarray(normale, dtype=float)
    n /= np.linalg.norm(n)
    #: versore radiale nel piano: perpendicolare a x e alla normale
    e = np.cross(n, np.array([1.0, 0.0, 0.0]))
    e /= np.linalg.norm(e)
    return punti[:, 0], punti[:, 1:] @ e[1:]


def fotogramma(path: Path, normale, campo: str, vmin, vmax, titolo, out: Path,
               cmap="inferno"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri

    punti, poligoni, campi = leggi_vtp(path)
    if campo not in campi or not poligoni:
        return False
    x, s = coordinate_nel_piano(punti, normale)
    tri = mtri.Triangulation(x * 1e3, s * 1e3, triangola(poligoni))
    v = campi[campo]
    if v.ndim > 1:
        v = np.linalg.norm(v, axis=1)
    fig, ax = plt.subplots(figsize=(11, 4.2), dpi=110)
    liv = np.linspace(vmin, vmax, 40)
    c = ax.tricontourf(tri, np.clip(v, vmin, vmax), levels=liv, cmap=cmap,
                       extend="both")
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("r nel piano di taglio [mm]")
    ax.set_title(titolo)
    fig.colorbar(c, ax=ax, label=campo, fraction=0.03, pad=0.02)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    return True


def main() -> int:
    from mesh_iniettore import quote_dal_progetto

    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", type=Path, required=True)
    ap.add_argument("--campo", default="C3H8")
    ap.add_argument("--piano", default="meridiano_getto",
                    choices=("meridiano_getto", "meridiano_fianco"))
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--fps", type=int, default=20)
    ap.add_argument("--vmax", type=float, default=None)
    a = ap.parse_args()

    d = quote_dal_progetto()
    ang = math.pi / d["n_getti"]
    normale = ((0.0, -math.sin(ang), math.cos(ang))
               if a.piano == "meridiano_getto" else (0.0, 0.0, 1.0))

    base = a.caso / "postProcessing" / "superfici"
    tempi = sorted((p for p in base.iterdir() if p.is_dir()),
                   key=lambda p: float(p.name))
    if not tempi:
        raise SystemExit("nessun istante campionato")

    #: la scala di colore si fissa UNA VOLTA sull'ultimo istante e non
    #: fotogramma per fotogramma: una scala che si riadatta fa sembrare
    #: uniforme anche una miscela che non lo e'.
    if a.vmax is None:
        p, poly, campi = leggi_vtp(tempi[-1] / f"{a.piano}.vtp")
        v = campi[a.campo]
        if v.ndim > 1:
            v = np.linalg.norm(v, axis=1)
        vmax = float(np.percentile(v, 99.5))
    else:
        vmax = a.vmax
    vmin = 0.0 if a.campo != "T" else 250.0
    out = a.out or (a.caso / f"film_{a.piano}_{a.campo}")
    out.mkdir(parents=True, exist_ok=True)

    n = 0
    for i, t in enumerate(tempi):
        f = t / f"{a.piano}.vtp"
        if not f.exists():
            continue
        ok = fotogramma(f, normale, a.campo, vmin, vmax,
                        f"{a.piano}   {a.campo}   t = {float(t.name)*1e6:7.2f} us",
                        out / f"f{i:05d}.png")
        n += int(ok)
    print(f"{n} fotogrammi in {out}")

    mp4 = out.with_suffix(".mp4")
    try:
        subprocess.run(
            ["ffmpeg", "-y", "-framerate", str(a.fps), "-pattern_type", "glob",
             "-i", str(out / "f*.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p",
             "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2", str(mp4)],
            check=True, capture_output=True)
        print(f"filmato: {mp4}")
    except (FileNotFoundError, subprocess.CalledProcessError) as e:
        print(f"ffmpeg non disponibile o fallito ({e}); restano i PNG")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
