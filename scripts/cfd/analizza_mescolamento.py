#!/usr/bin/env python3
"""Misura il mescolamento sui campi della CFD, e lo confronta con Holdeman.

LA DOMANDA. L'iniettore e' dimensionato su C = (S/H) sqrt(J) = 2.5, che
promette una miscela sostanzialmente uniforme entro circa due altezze di
condotto a valle dei getti. Qui si verifica se e' vero, e il numero che lo dice
e' la DISUNIFORMITA' della frazione di massa di combustibile su un piano
trasversale:

    U(x) = sqrt( < (Y - <Y>)^2 > ) / <Y>

pesata sulla portata, cioe' su rho*u*dA e non sull'area: quello che conta e'
com'e' distribuito il combustibile nel FLUSSO, non nello spazio. Un angolo
morto pieno di GPL fermo non alimenta niente e non deve contare come miscela.

U = 0 significa perfettamente miscelato; U = 1 significa che meta' del flusso
non ha combustibile per niente. La letteratura sui combustori considera
accettabile U < 0.1 all'ingresso della zona primaria.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))


def leggi_vtu(path: Path):
    """Legge un .vtu di foamToVTK senza dipendere da VTK: usa meshio se c'e',
    altrimenti il parser XML minimale qui sotto."""
    try:
        import meshio
        m = meshio.read(path)
        return m.points, {k: np.asarray(v) for k, v in m.point_data.items()}
    except ImportError:
        raise SystemExit(
            "serve meshio per leggere i file VTK:  pip install meshio")


def disuniformita(punti, dati, x0: float, spessore: float):
    """Disuniformita' pesata sulla portata in una fetta attorno a x0."""
    x = punti[:, 0]
    m = np.abs(x - x0) < 0.5 * spessore
    if m.sum() < 20:
        return None
    Y = dati["C3H8"][m]
    U = dati["U"][m]
    rho = dati.get("rho")
    w = (rho[m] if rho is not None else np.ones(m.sum())) * np.abs(U[:, 0])
    w = np.clip(w, 0.0, None)
    if w.sum() <= 0:
        return None
    Ym = np.average(Y, weights=w)
    if Ym <= 1e-12:
        return None
    var = np.average((Y - Ym) ** 2, weights=w)
    return math.sqrt(var) / Ym, Ym, m.sum()


def main() -> int:
    from mesh_iniettore import quote_dal_progetto

    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", type=Path, default=Path("runs/cfd/mescolamento"))
    ap.add_argument("--tempo", default=None, help="istante da leggere (default: l'ultimo)")
    a = ap.parse_args()

    vtk = a.caso / "VTK"
    if not vtk.exists():
        raise SystemExit(f"{vtk} non esiste: lancia prima  foamToVTK -case {a.caso}")
    interni = sorted(vtk.glob("*/internal.vtu"))
    if not interni:
        raise SystemExit(f"nessun internal.vtu in {vtk}")
    scelto = interni[-1] if a.tempo is None else next(
        p for p in interni if a.tempo in p.parent.name)
    print(f"istante letto: {scelto.parent.name}")

    punti, dati = leggi_vtu(scelto)
    mancanti = [k for k in ("C3H8", "U") if k not in dati]
    if mancanti:
        raise SystemExit(f"campi mancanti nel file: {mancanti}")

    d = quote_dal_progetto()
    Lm = d["L_mescolamento"]
    xg = d["x_getti"]
    print(f"\ngetti a x = {xg*1e3:.2f} mm; la correlazione promette miscela")
    print(f"uniforme entro {Lm*1e3:.2f} mm, cioe' a x = {(xg+Lm)*1e3:.2f} mm\n")
    print(f"  {'x [mm]':>9} {'x-x_getti':>11} {'in unita' di L_mix':>18} "
          f"{'Y medio':>10} {'disuniformita':>14}")
    stazioni = [xg + f * Lm for f in (0.25, 0.5, 0.75, 1.0, 1.5, 2.0)]
    stazioni += [0.0, 2.0e-3, 5.0e-3, 9.0e-3]
    for x0 in sorted(set(stazioni)):
        r = disuniformita(punti, dati, x0, 3.0e-4)
        if r is None:
            continue
        U, Ym, n = r
        nota = ""
        if abs(x0 - (xg + Lm)) < 1e-6:
            nota = "  <-- fine della lunghezza di mescolamento"
        if abs(x0) < 1e-9:
            nota = "  <-- sbocco in camera"
        print(f"  {x0*1e3:9.2f} {(x0-xg)*1e3:11.2f} {(x0-xg)/Lm:18.2f} "
              f"{Ym:10.4f} {U:14.3f}{nota}")

    r = disuniformita(punti, dati, xg + Lm, 3.0e-4)
    if r is not None:
        U = r[0]
        print(f"\nVERDETTO: disuniformita' = {U:.3f} alla fine della lunghezza")
        print("di mescolamento promessa dalla correlazione.")
        if U < 0.10:
            print("  La correlazione REGGE: il combustibile e' distribuito.")
        elif U < 0.25:
            print("  Mescolamento parziale: la correlazione e' ottimista, ma il")
            print("  combustibile e' comunque distribuito prima dello sbocco.")
        else:
            print("  La correlazione NON regge su questa geometria: i getti non")
            print("  penetrano come previsto e l'iniettore va ridisegnato.")
            print("  E' il motivo per cui questo caso esiste.")
    return 0


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
