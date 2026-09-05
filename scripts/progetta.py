#!/usr/bin/env python3
"""Dai requisiti al motore verificato.

    python scripts/progetta.py --spinta 50 --durata 5 --t-bombola 15

Stampa il progetto con la provenienza di ogni quota, oppure il rifiuto con
l'elenco di cio' che manca. Con --stl costruisce anche la geometria, la
verifica a due risoluzioni e scrive il file.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from zefiro.sintesi import Impianto, Processo, Requisiti, Rifiuto, progetta  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--spinta", type=float, default=50.0, help="N")
    ap.add_argument("--durata", type=float, default=5.0, help="s di regime")
    ap.add_argument("--t-bombola", type=float, default=15.0,
                    help="temperatura MINIMA di esercizio della bombola [C]. "
                         "Non e' una previsione meteo: e' un requisito operativo, "
                         "e si garantisce con un bagno d'acqua")
    ap.add_argument("--combustibile", default="C3H8:1.0",
                    help="composizione molare, es. 'C3H8:0.9,C4H10:0.1'")
    ap.add_argument("--t-aria", type=float, default=293.0,
                    help="K a valle del riduttore dell'aria (DA MISURARE)")
    ap.add_argument("--t-gpl", type=float, default=283.0,
                    help="K a valle del riduttore del GPL (DA MISURARE)")
    ap.add_argument("--cd", type=float, default=0.75,
                    help="coefficiente di efflusso degli iniettori (DA MISURARE)")
    ap.add_argument("--serbatoio", type=float, default=100.0, help="litri")
    ap.add_argument("--p-serbatoio", type=float, default=10.0, help="bar assoluti")
    ap.add_argument("--parete-minima", type=float, default=0.5, help="mm")
    ap.add_argument("--foro-minimo", type=float, default=0.6, help="mm")
    ap.add_argument("--stl", type=Path, default=None,
                    help="costruisci e verifica la geometria, e scrivi lo STL qui")
    ap.add_argument("--passo", type=float, default=1.8e-4, help="passo griglia [m]")
    a = ap.parse_args()

    comp = {}
    for pezzo in a.combustibile.split(","):
        k, v = pezzo.split(":")
        comp[k.strip()] = float(v)

    req = Requisiti(
        spinta=a.spinta, durata=a.durata,
        impianto=Impianto(
            combustibile=comp, T_bombola_min=a.t_bombola + 273.15,
            volume_serbatoio=a.serbatoio * 1e-3, p_serbatoio_max=a.p_serbatoio * 1e5,
            T_ossidante_iniezione=a.t_aria, T_combustibile_iniezione=a.t_gpl,
            cd_ossidante=a.cd, cd_combustibile=a.cd),
        processo=Processo(parete_minima=a.parete_minima * 1e-3,
                          foro_minimo=a.foro_minimo * 1e-3))

    t0 = time.perf_counter()
    p = progetta(req, geometria=a.stl is not None, passo=a.passo)
    print(p.referto())
    print(f"\n[{time.perf_counter()-t0:.1f} s]")

    if isinstance(p, Rifiuto):
        return 1
    if a.stl is not None and p.consegnabile:
        from zefiro.sdf.meshing import isosurface, write_stl
        from zefiro.sintesi.architetture import aerospike_gas_gas as arch
        m = arch.costruisci_geometria(p.geometria, a.passo)
        v, f = isosurface(m.solido)
        a.stl.parent.mkdir(parents=True, exist_ok=True)
        sha = write_stl(v, f, a.stl)
        print(f"scritto {a.stl}  {len(f)} triangoli  sha {sha[:12]}")
    elif a.stl is not None:
        print("\nSTL NON scritto: il progetto non ha passato tutte le verifiche.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
