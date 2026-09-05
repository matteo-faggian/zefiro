#!/usr/bin/env python3
"""Sistema i tipi di bordo dopo gmshToFoam.

`gmshToFoam` importa tutti i gruppi fisici come `patch` generiche. I due fianchi
del settore devono diventare `symmetryPlane` e le pareti `wall`, altrimenti:

  * con i fianchi come `patch`, il settore non e' un settore: sono due aperture
    e il fluido esce di lato. Il caso gira e il risultato e' privo di senso;
  * con le pareti come `patch`, le funzioni di parete della turbolenza non si
    attivano e lo strato limite non c'e'.

E' un passaggio obbligato che si dimentica facilmente, per questo e' uno script
e non una riga di istruzioni.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

SIMMETRIE = ("simmetria_1", "simmetria_2")
PARETI = ("pareti",)


def sistema(caso: Path) -> dict[str, str]:
    f = caso / "constant" / "polyMesh" / "boundary"
    if not f.exists():
        raise SystemExit(f"{f} non esiste: prima lancia gmshToFoam")
    s = f.read_text(encoding="utf-8")
    for n in SIMMETRIE:
        s = re.sub(rf"({n}\s*\{{\s*)type\s+\w+;", r"\1type            symmetryPlane;", s)
    for n in PARETI:
        s = re.sub(rf"({n}\s*\{{\s*)type\s+\w+;", r"\1type            wall;", s)
    f.write_text(s, encoding="utf-8")
    return dict(re.findall(r"(\w+)\s*\{\s*type\s+(\w+);", s))


if __name__ == "__main__":
    caso = Path(sys.argv[1] if len(sys.argv) > 1 else "runs/cfd/mescolamento")
    for nome, tipo in sistema(caso).items():
        print(f"  {nome:16s} {tipo}")
