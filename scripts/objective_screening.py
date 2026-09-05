#!/usr/bin/env python3
"""Quali obiettivi mettere in un Pareto, DIMOSTRATO invece che assunto.

Il rischio che questo script esclude e' preciso: due obiettivi quasi collineari
producono un fronte di Pareto degenere (una scheggia invece di una superficie),
e NSGA-II ci spreca una dimensione perdendo pressione selettiva. Prima di
scegliere gli obiettivi si campiona il box di progetto e si guarda la matrice
di correlazione.

Si usa la correlazione di **rango** (Spearman), non quella lineare: la
dominanza di Pareto e' una relazione d'ordine, quindi cio' che conta e' se due
metriche ordinano i candidati allo stesso modo, non se sono legate da una retta.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from zefiro.doe import latin_hypercube                       # noqa: E402
from zefiro.geometry.parameters import (                     # noqa: E402
    DESIGN_BOUNDS,
    INTEGER_PARAMETERS,
    derive,
)
from zefiro.l0.chemistry import BlowoutTable                 # noqa: E402
from zefiro.l0.mixture import MixtureModel                   # noqa: E402
from zefiro.schemas import DesignVector, FuelSpec, OperatingPoint  # noqa: E402

# --- punto operativo dello STUDIO ------------------------------------------ #
# Numeri di impianto reali (config/operating_point.yaml). Le due temperature
# sono ancora TODO G: qui sono PARAMETRI DELLO STUDIO, non misure, e il loro
# effetto e' quantificato dall'opzione --temperature-sweep.
STUDY_T_AIR = 293.0
STUDY_T_FUEL = 283.0
STUDY_CD = 0.75      # TODO I: Cd d'iniettore. Entra SOLO nei vincoli di Dp.


def study_operating_point(T_air: float, T_fuel: float) -> OperatingPoint:
    return OperatingPoint(
        p_amb=101325.0,
        p_air_supply=6.0e5,
        p_fuel_supply=8.0e5,
        fuel=FuelSpec(composition={"C3H8": 1.0}, phase_at_injection="gas",
                      thermo_source="gri30.yaml"),
        T_air_in=T_air,
        T_fuel_in=T_fuel,
        mdot_air_max=0.0718,
        cd_injector_ox=STUDY_CD,
        cd_injector_fuel=STUDY_CD,
    )


#: Metriche candidate a diventare obiettivi. Il segno e' quello della
#: MINIMIZZAZIONE, cosi' la correlazione si legge senza correggere mentalmente.
CANDIDATES = (
    "neg_thrust", "neg_isp_total", "neg_isp_fuel", "q_throat",
    "wall_volume", "neg_tau_res", "neg_damkohler", "T_adiabatic",
)


def evaluate(x_values, op, table):
    x = DesignVector(values=dict(x_values), bounds=DESIGN_BOUNDS)
    params, l0 = derive(x, op)
    d = params.derived
    tau_chem = table(x_values["p_c"], x_values["phi_core"])
    da = (l0.tau_res / tau_chem) if l0.tau_res else None
    return {
        "neg_thrust": -l0.thrust,
        "neg_isp_total": -l0.Isp_s,
        "neg_isp_fuel": -l0.Isp_fuel_s,
        "q_throat": l0.q_throat,
        "wall_volume": l0.wall_volume,
        "neg_tau_res": -l0.tau_res,
        "neg_damkohler": -da,
        "T_adiabatic": l0.T_ad,
        "_tau_res": l0.tau_res,
        "_tau_chem": tau_chem,
        "_damkohler": da,
        "_L_star": d["L_star"],
        "_eps": l0.epsilon,
    }


def spearman(M: np.ndarray) -> np.ndarray:
    R = np.apply_along_axis(lambda c: np.argsort(np.argsort(c)).astype(float), 0, M)
    R -= R.mean(axis=0)
    R /= np.linalg.norm(R, axis=0)
    return R.T @ R


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", type=int, default=400)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", type=Path, default=Path("runs/objective_screening.json"))
    a = ap.parse_args()

    op = study_operating_point(STUDY_T_AIR, STUDY_T_FUEL)
    model = MixtureModel.from_fuel(op.fuel)

    print("tabulazione del tempo chimico (PSR blowout)...", flush=True)
    table = BlowoutTable.build(
        model,
        np.linspace(DESIGN_BOUNDS["p_c"][0], DESIGN_BOUNDS["p_c"][1], 4),
        np.linspace(DESIGN_BOUNDS["phi_core"][0], DESIGN_BOUNDS["phi_core"][1], 5),
        STUDY_T_AIR, STUDY_T_FUEL,
    )
    tt = np.exp(np.array(table.log_tau))
    print(f"  tau_chem su griglia: {tt.min()*1e6:.1f} - {tt.max()*1e6:.1f} us")

    X = latin_hypercube(a.n, DESIGN_BOUNDS, seed=a.seed,
                        integer_names=INTEGER_PARAMETERS)
    rows, scartati = [], []
    for i, xv in enumerate(X):
        try:
            rows.append(evaluate(xv, op, table))
        except Exception as e:                                   # noqa: BLE001
            scartati.append(f"{type(e).__name__}: {e}")
    print(f"valutati {len(rows)}/{a.n}  (scartati {len(scartati)})")
    if scartati:
        from collections import Counter
        for msg, k in Counter(s.split(":")[0] for s in scartati).most_common():
            print(f"  {k:4d} x {msg}")

    M = np.array([[r[c] for c in CANDIDATES] for r in rows], dtype=float)
    C = spearman(M)

    w = max(len(c) for c in CANDIDATES) + 1
    print("\nCORRELAZIONE DI RANGO (Spearman) fra i candidati obiettivo")
    print(" " * w + "".join(f"{c[:7]:>8}" for c in CANDIDATES))
    for i, c in enumerate(CANDIDATES):
        print(f"{c:<{w}}" + "".join(f"{C[i, j]:8.3f}" for j in range(len(CANDIDATES))))

    print("\nCOPPIE QUASI DEGENERI (|rho| > 0.90): un solo obiettivo delle due basta")
    degeneri = []
    for i in range(len(CANDIDATES)):
        for j in range(i + 1, len(CANDIDATES)):
            if abs(C[i, j]) > 0.90:
                degeneri.append((CANDIDATES[i], CANDIDATES[j], float(C[i, j])))
                print(f"  {CANDIDATES[i]:>14} <-> {CANDIDATES[j]:<14} rho = {C[i,j]:+.4f}")
    if not degeneri:
        print("  nessuna")

    da = np.array([r["_damkohler"] for r in rows])
    tr = np.array([r["_tau_res"] for r in rows])
    print(f"\nDAMKOHLER sul box:  min {da.min():.1f}   mediana {np.median(da):.1f}"
          f"   max {da.max():.1f}")
    print(f"  tau_res  {tr.min()*1e3:.2f} - {tr.max()*1e3:.2f} ms")
    print("  -> la cinetica " + ("NON e' vincolante" if da.min() > 10
          else "PUO' essere vincolante") + f" (Da_min = {da.min():.1f})")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    a.out.write_text(json.dumps({
        "n": a.n, "seed": a.seed, "n_valid": len(rows),
        "candidates": list(CANDIDATES),
        "spearman": C.tolist(),
        "degenerate_pairs": degeneri,
        "damkohler_min": float(da.min()), "damkohler_max": float(da.max()),
        "tau_res_s": [float(tr.min()), float(tr.max())],
        "T_air_in": STUDY_T_AIR, "T_fuel_in": STUDY_T_FUEL,
    }, indent=2))
    print(f"\nscritto {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
