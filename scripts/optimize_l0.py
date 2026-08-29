#!/usr/bin/env python3
"""Fronte di Pareto su L0 con NSGA-II. E' la run lunga: ore, non minuti.

    python scripts/optimize_l0.py --pop 200 --gen 300 --seed 1 --jobs 23

Cosa fa, in ordine:
  1. tabula il tempo chimico (blowout di un PSR) sulla griglia (p_c, phi);
  2. calcola quanto flusso termico il raffreddamento DICHIARATO riesce a
     estrarre, e lo usa come vincolo;
  3. gira NSGA-II sui 12 parametri liberi;
  4. scrive il fronte in JSON e le run nel database.

Le scelte discutibili sono tutte in testa al file, dichiarate.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from zefiro.geometry.parameters import DESIGN_BOUNDS                  # noqa: E402
from zefiro.l0.chemistry import BlowoutTable                          # noqa: E402
from zefiro.l0.mixture import MixtureModel                            # noqa: E402
from zefiro.opt.driver import (                                       # noqa: E402
    OBJECTIVE_SETS,
    OptimizationSetup,
    default_n_jobs,
    front_quality,
    run_nsga2,
)
from zefiro.schemas import FuelSpec, OperatingPoint                   # noqa: E402

# --- CONFIGURAZIONE DEL RAFFREDDAMENTO [DICHIARATA] ------------------------- #
# Il vincolo termico non e' "q sotto un numero tondo": e' "q sotto cio' che
# QUESTO canale, con QUESTA acqua, riesce a portare via senza far bollire".
# Cambiare il canale cambia il vincolo, ed e' giusto cosi'.
COOLING = {
    "hydraulic_diameter": 0.5e-3,   # m   canale rettangolare stretto
    "velocity": 12.0,               # m/s
    "pressure": 4.0e5,              # Pa  pressione di rete
    "T_inlet": 288.0,               # K   [TODO H] temperatura dell'acqua di rete
    "subcooling_margin": 20.0,      # K   margine alla saturazione
}

# TODO G: temperature a valle dei riduttori. Sono PARAMETRI DI STUDIO.
STUDY_T_AIR = 293.0
STUDY_T_FUEL = 283.0
STUDY_CD = 0.75                     # TODO I


def q_removable() -> tuple[float, dict]:
    """Flusso estraibile [W/m^2] dal criterio di NON EBOLLIZIONE.

    Criterio scelto: la parete lato acqua deve restare sotto la saturazione con
    un margine. E' piu' restrittivo del flusso critico e non richiede nessuna
    proprieta' del 316L, che e' ancora un TODO (J): dipende solo da acqua e
    geometria del canale, cioe' da cose note.
    """
    from zefiro.cooling import dittus_boelter_h, reynolds, saturation_temperature

    c = COOLING
    h = dittus_boelter_h(c["hydraulic_diameter"], c["velocity"], c["T_inlet"],
                         c["pressure"])
    re = reynolds(c["hydraulic_diameter"], c["velocity"], c["T_inlet"], c["pressure"])
    t_sat = saturation_temperature(c["pressure"])
    dt = t_sat - c["T_inlet"] - c["subcooling_margin"]
    return h * dt, {"h": h, "Re": re, "T_sat": t_sat, "dT_utile": dt}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pop", type=int, default=200)
    ap.add_argument("--gen", type=int, default=300)
    ap.add_argument("--seed", type=int, required=True,
                    help="primo seme; con --seeds N si usano seed, seed+1, ...")
    ap.add_argument("--seeds", type=int, default=1,
                    help="quante ripetizioni indipendenti (verifica di robustezza)")
    ap.add_argument("--jobs", type=int, default=0, help="0 = automatico")
    ap.add_argument("--objectives", default="zefiro", choices=sorted(OBJECTIVE_SETS))
    ap.add_argument("--out", type=Path, default=Path("runs/pareto"))
    ap.add_argument("--no-thermal-constraint", action="store_true")
    a = ap.parse_args()

    op = OperatingPoint(
        p_amb=101325.0, p_air_supply=6.0e5, p_fuel_supply=8.0e5,
        fuel=FuelSpec(composition={"C3H8": 1.0}, phase_at_injection="gas",
                      thermo_source="gri30.yaml"),
        T_air_in=STUDY_T_AIR, T_fuel_in=STUDY_T_FUEL, mdot_air_max=0.0718,
        cd_injector_ox=STUDY_CD, cd_injector_fuel=STUDY_CD,
    )

    q_max, cool = q_removable()
    print(f"raffreddamento: h = {cool['h']/1e3:.1f} kW/(m2 K)  Re = {cool['Re']:.0f}"
          f"  T_sat = {cool['T_sat']-273.15:.1f} C")
    if cool["Re"] < 1.0e4:
        print(f"  ATTENZIONE: Re = {cool['Re']:.0f} < 1e4, fuori dal campo di validita'")
        print("  di Dittus-Boelter. h e' OTTIMISTICO, quindi q_max lo e' altrettanto.")
    print(f"  q estraibile = {q_max/1e6:.2f} MW/m2"
          + ("  [VINCOLO DISATTIVATO]" if a.no_thermal_constraint else ""))

    print("tabulazione del tempo chimico...", flush=True)
    t0 = time.perf_counter()
    table = BlowoutTable.build(
        MixtureModel.from_fuel(op.fuel),
        np.linspace(*DESIGN_BOUNDS["p_c"], 5),
        np.linspace(*DESIGN_BOUNDS["phi_core"], 6),
        STUDY_T_AIR, STUDY_T_FUEL,
    )
    tt = np.exp(np.array(table.log_tau))
    print(f"  tau_chem {tt.min()*1e6:.1f} - {tt.max()*1e6:.1f} us "
          f"({time.perf_counter()-t0:.0f} s)")

    setup = OptimizationSetup(
        op=op, objective_set=a.objectives, tau_chem=table,
        q_removable=None if a.no_thermal_constraint else q_max,
    )
    jobs = a.jobs or default_n_jobs()
    print(f"NSGA-II  obiettivi {setup.objective_names}")
    print(f"  pop {a.pop} x gen {a.gen} = {a.pop*a.gen} valutazioni su {jobs} processi")

    semi = [a.seed + i for i in range(a.seeds)]
    runs = []
    for s_i in semi:
        print(f"\n--- seme {s_i} ({semi.index(s_i)+1}/{len(semi)}) ---")
        r = run_nsga2(setup, pop_size=a.pop, n_gen=a.gen, seed=s_i, n_jobs=jobs,
                      verbose=(len(semi) == 1),
                      checkpoint=a.out / f"pareto_seed{s_i}.json")
        print(f"  {len(r.X):4d} punti, {r.n_eval} valutazioni, "
              f"{r.wall_time_s/60:.1f} min ({r.wall_time_s/max(r.n_eval,1)*1e3:.1f} ms/val)")
        if r.n_failed:
            print(f"  valutazioni fallite: {r.n_failed} ({r.n_failed/r.n_eval:.1%})")
            for k, v in sorted(r.failures.items(), key=lambda kv: -kv[1])[:5]:
                print(f"    {v:6d} x {k}")
        runs.append(r)

    res = runs[0]
    qual = None
    if len(runs) > 1:
        qual = front_quality(runs)
        print(f"\nROBUSTEZZA AL SEME  (ipervolume normalizzato, {len(runs)} semi)")
        for s_i, hv, n in zip(qual["seeds"], qual["hypervolume"], qual["n_points"]):
            print(f"  seme {s_i:5d}: HV = {hv:.4f}   {n:4d} punti")
        print(f"  dispersione = {qual['hv_spread']:.2%}", end="  ")
        print("-> il fronte NON dipende dal seme" if qual["hv_spread"] < 0.05
              else "-> TROPPO ALTA: budget insufficiente, le conclusioni non reggono")

    print("\nESTREMI DEL FRONTE (unione dei semi)")
    F = np.vstack([r.F for r in runs])
    for j, n in enumerate(res.objective_names):
        seg = -1.0 if n.startswith("neg_") else 1.0
        v = seg * F[:, j]
        print(f"  {n:<16} da {v.min():12.4g} a {v.max():12.4g}")

    meta = a.out / f"pareto_seed{a.seed}_meta.json"
    meta.parent.mkdir(parents=True, exist_ok=True)
    meta.write_text(json.dumps({
        "cooling": COOLING, "q_removable": q_max, "cooling_derived": cool,
        "T_air_in": STUDY_T_AIR, "T_fuel_in": STUDY_T_FUEL, "cd": STUDY_CD,
        "objective_set": a.objectives,
        "tau_chem_range_s": [float(tt.min()), float(tt.max())],
        "seeds": semi, "pop": a.pop, "gen": a.gen,
        "front_quality": qual,
    }, indent=2))
    print(f"\nscritti {a.out}/pareto_seed{a.seed}.json  e  {meta.name}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
