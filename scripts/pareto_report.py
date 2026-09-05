#!/usr/bin/env python3
"""Legge un fronte di Pareto e lo rende UTILIZZABILE.

Un fronte e' un elenco di compromessi, non una risposta. Questo script prova a
rispondere alle tre domande che un progettista fa davvero:

  1. quali sono gli estremi, e quanto costa passare da un estremo all'altro?
  2. c'e' un GINOCCHIO, cioe' un punto oltre il quale migliorare un obiettivo
     costa sproporzionatamente sugli altri?
  3. QUALI PARAMETRI muovono il compromesso, e quali sono gia' decisi (cioe'
     saturi allo stesso valore su tutto il fronte)? I secondi non sono
     compromessi: sono conclusioni di progetto.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from zefiro.geometry.parameters import DESIGN_BOUNDS      # noqa: E402

#: Come si legge ciascun obiettivo: (etichetta, unita', fattore, verso).
#: `verso` = +1 se piu' grande e' meglio DOPO la conversione di segno.
LETTURA = {
    "neg_thrust":    ("spinta", "N", -1.0, +1),
    "neg_isp_total": ("Isp", "s", -1.0, +1),
    "neg_isp_fuel":  ("Isp sul GPL", "s", -1.0, +1),
    "neg_damkohler": ("Damkohler", "-", -1.0, +1),
    "q_throat":      ("flusso in gola", "MW/m2", 1.0e-6, -1),
    "wall_volume":   ("volume del pezzo", "cm3", 1.0e6, -1),
}


def carica(paths: list[Path]):
    X, F, nomi_f, nomi_x = [], [], None, None
    for p in paths:
        d = json.loads(p.read_text())
        if nomi_f is None:
            nomi_f = tuple(d["objective_names"])
            nomi_x = tuple(d["variable_names"])
        elif tuple(d["objective_names"]) != nomi_f:
            raise SystemExit(f"{p} ha obiettivi diversi: non si uniscono fronti diversi")
        X.append(np.array(d["X"], dtype=float))
        F.append(np.array(d["F"], dtype=float))
    return np.vstack(X), np.vstack(F), nomi_f, nomi_x


def non_dominati(F: np.ndarray) -> np.ndarray:
    """Filtro di dominanza esplicito: unendo piu' semi si ottengono punti che
    un altro seme domina, e tenerli falserebbe gli estremi."""
    n = len(F)
    tieni = np.ones(n, dtype=bool)
    for i in range(n):
        if not tieni[i]:
            continue
        dom = np.all(F <= F[i], axis=1) & np.any(F < F[i], axis=1)
        if dom.any():
            tieni[i] = False
    return tieni


def _spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    ra -= ra.mean()
    rb -= rb.mean()
    d = np.linalg.norm(ra) * np.linalg.norm(rb)
    return float(ra @ rb / d) if d else 0.0


def fmt(nome: str, valore: float) -> str:
    et, un, k, _ = LETTURA.get(nome, (nome, "", 1.0, -1))
    return f"{valore*k:>10.4g} {un}"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("fronti", nargs="+", type=Path)
    ap.add_argument("--plot", type=Path, default=None)
    ap.add_argument("--saturi-entro", type=float, default=0.05,
                    help="frazione del bound entro cui un parametro e' 'saturo'")
    a = ap.parse_args()

    X, F, nomi_f, nomi_x = carica(a.fronti)
    m = non_dominati(F)
    X, F = X[m], F[m]
    print(f"fronte unito: {len(F)} punti non dominati su {len(m)} "
          f"({len(a.fronti)} file)\n")

    # --- 1. estremi -------------------------------------------------------- #
    print("ESTREMI  (ogni riga: il migliore in quell'obiettivo, e cosa costa)")
    intest = "".join(f"{LETTURA.get(n,(n,))[0][:16]:>18}" for n in nomi_f)
    print(f"{'ottimizzando':<18}{intest}")
    idx_best = {}
    for j, n in enumerate(nomi_f):
        i = int(np.argmin(F[:, j]))
        idx_best[n] = i
        et = LETTURA.get(n, (n,))[0]
        print(f"{et[:17]:<18}" + "".join(fmt(nn, F[i, jj]).rjust(18)
                                         for jj, nn in enumerate(nomi_f)))

    # --- 2. ginocchio ------------------------------------------------------ #
    # Normalizzazione min-max sul fronte, poi il punto piu' vicino all'ideale
    # in norma euclidea. E' il criterio piu' semplice che si possa dichiarare;
    # non e' "il" punto giusto, e' un punto di partenza ragionevole.
    ideale, nadir = F.min(axis=0), F.max(axis=0)
    span = np.where(nadir > ideale, nadir - ideale, 1.0)
    Fn = (F - ideale) / span
    i_kn = int(np.argmin(np.linalg.norm(Fn, axis=1)))
    print("\nGINOCCHIO (piu' vicino all'ideale in norma normalizzata)")
    for j, n in enumerate(nomi_f):
        et = LETTURA.get(n, (n,))[0]
        print(f"  {et:<20}{fmt(n, F[i_kn, j])}   (a {Fn[i_kn,j]:.0%} del percorso "
              f"dall'ideale al nadir)")
    print("  vettore di progetto:")
    for k, nx in enumerate(nomi_x):
        lo, hi = DESIGN_BOUNDS[nx]
        print(f"    {nx:<16}{X[i_kn,k]:12.5g}   [{lo:.4g}, {hi:.4g}]"
              f"  {'  <- al bordo' if min(X[i_kn,k]-lo, hi-X[i_kn,k]) < 0.02*(hi-lo) else ''}")

    # --- 3. quali parametri sono decisi e quali sono compromesso ----------- #
    print("\nPARAMETRI: chi e' gia' DECISO e chi resta un compromesso")
    print(f"{'parametro':<18}{'min':>12}{'max':>12}{'escursione':>12}   stato")
    for k, nx in enumerate(nomi_x):
        lo, hi = DESIGN_BOUNDS[nx]
        v = X[:, k]
        frazione = (v.max() - v.min()) / (hi - lo)
        if frazione < a.saturi_entro:
            centro = 0.5 * (v.min() + v.max())
            dove = ("al MINIMO del box" if centro - lo < 0.1 * (hi - lo) else
                    "al MASSIMO del box" if hi - centro < 0.1 * (hi - lo) else
                    "su un valore interno")
            stato = f"DECISO, {dove}"
        else:
            stato = "compromesso"
        print(f"{nx:<18}{v.min():12.5g}{v.max():12.5g}{frazione:11.1%}   {stato}")
    print("\n  Un parametro 'DECISO' non e' un compromesso: su tutto il fronte")
    print("  l'ottimizzatore gli da' lo stesso valore, quindi e' una conclusione.")
    print("  Se sta al bordo del box, il box e' stretto: il vero ottimo e' fuori,")
    print("  e il bound va rimesso in discussione invece che subito.")

    # --- 4. quanto costa DAVVERO ciascun obiettivo -------------------------- #
    # La correlazione sul box dice se vale la pena mettere due metriche nel
    # fronte. Ma il compromesso VERO si legge solo sul fronte, e puo' essere
    # tutt'altra cosa: due obiettivi indipendenti sul box possono risultare
    # legatissimi una volta convergiti, o - come qui - uno dei tre puo'
    # rivelarsi quasi gratuito. E' la differenza fra "ha senso ottimizzarlo" e
    # "che prezzo ha".
    print("\nCORRELAZIONE SUL FRONTE  (da confrontare con quella sul box)")
    for i in range(len(nomi_f)):
        for j in range(i + 1, len(nomi_f)):
            ei = LETTURA.get(nomi_f[i], (nomi_f[i],))[0]
            ej = LETTURA.get(nomi_f[j], (nomi_f[j],))[0]
            print(f"  {ei:>20} <-> {ej:<20} rho = {_spearman(F[:, i], F[:, j]):+.3f}")
    print("  Attenzione: sul FRONTE una correlazione vicina a -1 e' normale")
    print("  (e' cio' che un fronte e'). Quella che conta per la scelta degli")
    print("  obiettivi e' la correlazione sul BOX: vedi objective_screening.py.")

    if len(nomi_f) >= 3:
        print("\nIL PREZZO DI OGNI OBIETTIVO  (a parita' degli altri, per quanto possibile)")
        for j, n in enumerate(nomi_f):
            altri = [t for t in range(len(nomi_f)) if t != j]
            rif = altri[0]
            et, un, k, _ = LETTURA.get(n, (n, "", 1.0, -1))
            v = F[:, j] * k
            # bande strette sul primo degli altri obiettivi: dentro una banda,
            # quegli obiettivi sono quasi costanti e resta il solo effetto di n
            bordi = np.quantile(F[:, rif], np.linspace(0, 1, 5))
            perdite = []
            for b0, b1 in zip(bordi[:-1], bordi[1:]):
                m = (F[:, rif] >= b0) & (F[:, rif] <= b1)
                if m.sum() < 8:
                    continue
                sec = [t for t in altri if t != rif][0]
                vv, ss = v[m], F[m, sec]
                i_lo, i_hi = int(np.argmin(vv)), int(np.argmax(vv))
                base = abs(ss[i_hi]) or 1.0
                perdite.append((ss[i_lo] - ss[i_hi]) / base)
            if perdite:
                pk = max(abs(x) for x in perdite)
                sec_n = LETTURA.get(nomi_f[[t for t in altri if t != rif][0]], ("",))[0]
                print(f"  portare {et:<18} da {v.max():10.4g} a {v.min():10.4g} {un:<5}"
                      f" costa al piu' {pk:6.2%} su {sec_n}")
        print("  Un costo di pochi punti percentuali su un'escursione di un ordine")
        print("  di grandezza NON e' un compromesso: e' una decisione gia' presa.")

    if a.plot:
        _plot(X, F, nomi_f, nomi_x, i_kn, a.plot)
        print(f"\nscritto {a.plot}")
    return 0


def _plot(X, F, nomi_f, nomi_x, i_kn, out: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    nf = len(nomi_f)
    coppie = [(i, j) for i in range(nf) for j in range(i + 1, nf)]
    fig = plt.figure(figsize=(5.2 * max(len(coppie), 1), 8.6))
    gs = fig.add_gridspec(2, max(len(coppie), 1), height_ratios=[1.0, 0.85],
                          hspace=0.36, wspace=0.52)

    # terzo obiettivo come colore: un fronte 3D letto su piani 2D perde
    # informazione, e il colore la restituisce senza una proiezione prospettica
    for k, (i, j) in enumerate(coppie):
        ax = fig.add_subplot(gs[0, k])
        altro = [t for t in range(nf) if t not in (i, j)]
        ei, ui, ki, _ = LETTURA.get(nomi_f[i], (nomi_f[i], "", 1.0, -1))
        ej, uj, kj, _ = LETTURA.get(nomi_f[j], (nomi_f[j], "", 1.0, -1))
        if altro:
            t = altro[0]
            et, ut, kt, _ = LETTURA.get(nomi_f[t], (nomi_f[t], "", 1.0, -1))
            sc = ax.scatter(F[:, i] * ki, F[:, j] * kj, c=F[:, t] * kt,
                            s=22, cmap="viridis", edgecolor="none")
            cb = fig.colorbar(sc, ax=ax, pad=0.03, fraction=0.045)
            cb.set_label(f"{et} [{ut}]", fontsize=9)
        else:
            ax.scatter(F[:, i] * ki, F[:, j] * kj, s=22, color="#2b6cb0")
        ax.scatter([F[i_kn, i] * ki], [F[i_kn, j] * kj], s=170, marker="*",
                   color="#c53030", zorder=5, label="ginocchio")
        ax.set_xlabel(f"{ei} [{ui}]")
        ax.set_ylabel(f"{ej} [{uj}]")
        ax.grid(alpha=0.25)
        ax.legend(loc="best", fontsize=9, frameon=False)

    # coordinate parallele sui parametri: mostra a colpo d'occhio chi e' saturo
    ax = fig.add_subplot(gs[1, :])
    lo = np.array([DESIGN_BOUNDS[n][0] for n in nomi_x])
    hi = np.array([DESIGN_BOUNDS[n][1] for n in nomi_x])
    Xn = (X - lo) / (hi - lo)
    colore = (F[:, 0] - F[:, 0].min()) / max(float(np.ptp(F[:, 0])), 1e-30)
    cmap = matplotlib.colormaps["viridis"]
    for r in range(len(Xn)):
        ax.plot(range(len(nomi_x)), Xn[r], color=cmap(colore[r]), alpha=0.28, lw=0.8)
    ax.plot(range(len(nomi_x)), Xn[i_kn], color="#c53030", lw=2.4, label="ginocchio")
    ax.set_xticks(range(len(nomi_x)))
    ax.set_xticklabels(nomi_x, rotation=38, ha="right", fontsize=9)
    ax.set_ylim(-0.04, 1.04)
    ax.set_ylabel("posizione nel box  (0 = minimo, 1 = massimo)")
    ax.set_title("Parametri sul fronte: una linea piatta a 0 o a 1 e' un parametro "
                 "saturo, cioe' deciso", fontsize=10)
    ax.grid(alpha=0.25, axis="y")
    ax.legend(loc="upper right", fontsize=9, frameon=False)

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=140, bbox_inches="tight")


if __name__ == "__main__":
    raise SystemExit(main())
