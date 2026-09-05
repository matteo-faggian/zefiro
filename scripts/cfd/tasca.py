#!/usr/bin/env python3
"""Che cosa c'e' dentro la tasca cieca davanti ai getti.

LA TASCA. `mesh_iniettore.costruisci` crea il cilindro del getto lungo
(R_anello - r_bore) + 2 mm e lo FONDE con il dominio fluido. Il "+2 mm" e' il
trucco standard per evitare facce coincidenti in un'operazione booleana - ed e'
giusto per un TAGLIO, dove l'eccedenza sparisce nel solido. Qui pero' e' una
FUSIONE: l'eccedenza resta, e diventa fluido. Il risultato e' una tasca cieca
di diametro d_getto e profonda 2 mm scavata nel metallo, esattamente di fronte
a ogni getto, che nel pezzo vero non c'e'.

Questo script non discute se sia grave: lo misura. Riporta che cosa contiene la
tasca (combustibile, velocita', temperatura) e lo confronta con il condotto
sano alla stessa stazione, cosi' la domanda "conta o non conta" ha un numero.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from filmato import coordinate_nel_piano  # noqa: E402
from sezioni import leggi_vtp  # noqa: E402


def statistiche(campi, dentro, nome):
    n = int(dentro.sum())
    if n == 0:
        print(f"  {nome:<28} nessun punto")
        return
    Y = campi["C3H8"][dentro]
    U = np.linalg.norm(campi["U"][dentro], axis=1)
    T = campi["T"][dentro]
    print(f"  {nome:<28} {n:5d} punti   Y_GPL {Y.mean():7.4f} (max {Y.max():6.4f})"
          f"   |U| {U.mean():7.2f} m/s (max {U.max():6.1f})   T {T.mean():6.1f} K")


def main() -> int:
    from mesh_iniettore import quote_dal_progetto

    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", type=Path, required=True)
    ap.add_argument("--campionamento", default="superfici2")
    a = ap.parse_args()

    d = quote_dal_progetto()
    ang = math.pi / d["n_getti"]
    base = a.caso / "postProcessing" / a.campionamento
    t = sorted((p for p in base.iterdir() if p.is_dir()),
               key=lambda p: float(p.name))[-1]
    punti, _, campi = leggi_vtp(t / "meridiano_getto.vtp")
    x, s = coordinate_nel_piano(punti, (0.0, -math.sin(ang), math.cos(ang)))

    R_e, R_i = d["R_anello"], d["R_getti"]
    xg, dg = d["x_getti"], d["d_getto"]
    print(f"istante {float(t.name)*1e6:.1f} us\n")
    print(f"tasca:   r da {R_e*1e3:.3f} a {(R_e+2.0e-3)*1e3:.3f} mm, "
          f"x da {(xg-dg/2)*1e3:.3f} a {(xg+dg/2)*1e3:.3f} mm")
    print(f"anello:  r da {R_i*1e3:.3f} a {R_e*1e3:.3f} mm\n")

    tasca = (s > R_e) & (s < R_e + 2.0e-3) & (np.abs(x - xg) < dg)
    #: il condotto alla STESSA stazione assiale: e' il confronto che conta,
    #: perche' confrontare la tasca con la media di tutto il dominio direbbe
    #: solo che il dominio e' grande.
    condotto = (s > R_i) & (s < R_e) & (np.abs(x - xg) < dg)
    #: e il condotto una lunghezza di mescolamento a valle, dove il getto
    #: dovrebbe aver attraversato: serve a dire se il getto ARRIVA al raggio
    #: della tasca oppure viene piegato prima.
    valle = (s > R_e - 0.4e-3) & (s < R_e) & (np.abs(x - (xg + d["L_mescolamento"])) < dg)

    statistiche(campi, tasca, "TASCA (non esiste nel pezzo)")
    statistiche(campi, condotto, "anello, stessa stazione")
    statistiche(campi, valle, "anello, parete esterna a valle")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
