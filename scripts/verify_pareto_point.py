#!/usr/bin/env python3
"""Verifica INDIPENDENTE di un punto del fronte.

Indipendente vuol dire: si riparte dal vettore x scritto nel JSON, si rifa'
tutto da capo - termochimica, geometria, CAD vero con OCCT, vincoli - senza
riusare nulla di quello che l'ottimizzatore aveva in memoria. Se l'ottimizzatore
avesse un difetto (un vincolo con il segno sbagliato, un obiettivo letto dalla
colonna sbagliata), questo controllo lo vedrebbe; rileggere i suoi numeri no.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from zefiro.geometry.build import BuildOptions, build_and_export      # noqa: E402
from zefiro.geometry.parameters import (                              # noqa: E402
    DESIGN_BOUNDS,
    check_manufacturability,
    derive,
)
from zefiro.l0.chemistry import BlowoutTable                          # noqa: E402
from zefiro.l0.mixture import MixtureModel                            # noqa: E402
from zefiro.opt.objectives import objectives_l0                       # noqa: E402
from zefiro.schemas import DesignVector                               # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from optimize_l0 import (                                             # noqa: E402
    STUDY_CD,
    STUDY_T_AIR,
    STUDY_T_FUEL,
    q_removable,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("fronte", type=Path)
    ap.add_argument("--indice", type=int, default=None,
                    help="indice del punto; se assente si prende il ginocchio")
    ap.add_argument("--cad", type=Path, default=None, help="cartella per STEP/STL")
    a = ap.parse_args()

    d = json.loads(a.fronte.read_text())
    X = np.array(d["X"], dtype=float)
    F = np.array(d["F"], dtype=float)
    nomi_x = d["variable_names"]

    if a.indice is None:
        ideale, nadir = F.min(axis=0), F.max(axis=0)
        span = np.where(nadir > ideale, nadir - ideale, 1.0)
        i = int(np.argmin(np.linalg.norm((F - ideale) / span, axis=1)))
        print(f"punto: ginocchio, indice {i} di {len(F)}")
    else:
        i = a.indice
        print(f"punto: indice {i} di {len(F)}")

    from zefiro.schemas import FuelSpec, OperatingPoint
    op = OperatingPoint(
        p_amb=101325.0, p_air_supply=6.0e5, p_fuel_supply=8.0e5,
        fuel=FuelSpec(composition={"C3H8": 1.0}, phase_at_injection="gas",
                      thermo_source="gri30.yaml"),
        T_air_in=STUDY_T_AIR, T_fuel_in=STUDY_T_FUEL, mdot_air_max=0.0718,
        cd_injector_ox=STUDY_CD, cd_injector_fuel=STUDY_CD,
    )
    valori = {n: float(v) for n, v in zip(nomi_x, X[i])}
    valori["N_inj"] = float(round(valori["N_inj"]))
    x = DesignVector(values=valori, bounds=DESIGN_BOUNDS)

    params, l0 = derive(x, op)
    tab = BlowoutTable.build(MixtureModel.from_fuel(op.fuel),
                             np.linspace(*DESIGN_BOUNDS["p_c"], 5),
                             np.linspace(*DESIGN_BOUNDS["phi_core"], 6),
                             STUDY_T_AIR, STUDY_T_FUEL)
    q_max, _ = q_removable()
    o = objectives_l0("verifica", l0, op, valori["p_c"], geometry=None,
                      derived=params.derived,
                      tau_chem=tab(valori["p_c"], valori["phi_core"]),
                      q_removable=q_max)

    print("\n1. GLI OBIETTIVI COINCIDONO CON QUELLI SCRITTI NEL FRONTE?")
    ok = True
    for j, n in enumerate(d["objective_names"]):
        atteso, ottenuto = F[i, j], o.f[n]
        rel = abs(ottenuto - atteso) / max(abs(atteso), 1e-30)
        buono = rel < 1e-9
        ok &= buono
        print(f"   {n:<16}{atteso:14.7g} vs {ottenuto:14.7g}   "
              f"{'ok' if buono else f'DIVERGE ({rel:.2e})'}")

    print("\n2. I VINCOLI SONO RISPETTATI?  (g <= 0 e' fattibile)")
    for n, v in sorted(o.g.items()):
        stato = "ok" if v <= 0 else "VIOLATO"
        attivo = "  <- attivo" if -0.02 < v <= 0 else ""
        print(f"   {n:<24}{v:+12.5f}   {stato}{attivo}")
        ok &= v <= 0
    non_valutati = [n for n in ("min_feature", "watertight")
                    if n not in o.g]
    if non_valutati:
        print(f"   non valutabili a L0 (dato mancante): {non_valutati}")

    print("\n3. FISICA DEL PUNTO")
    dd = params.derived
    print(f"   spinta            {l0.thrust:8.2f} N        Isp {l0.Isp_s:7.1f} s"
          f"   Isp/GPL {l0.Isp_fuel_s:7.0f} s")
    print(f"   p_c               {valori['p_c']/1e5:8.3f} bar      T_ad {l0.T_ad:7.1f} K"
          f"   eps {l0.epsilon:6.3f}")
    print(f"   portate           aria {l0.mdot_air*1e3:6.1f} g/s   GPL "
          f"{(l0.mdot_fuel_core+l0.mdot_fuel_film)*1e3:5.2f} g/s")
    print(f"   q di gola         {l0.q_throat/1e6:8.3f} MW/m2    su un estraibile di "
          f"{q_max/1e6:.2f} MW/m2")
    print(f"   Damkohler         {l0.tau_res/tab(valori['p_c'], valori['phi_core']):8.1f}"
          f"        (tau_res {l0.tau_res*1e3:.2f} ms)")
    print(f"   volume            {l0.wall_volume*1e6:8.2f} cm3      massa in 316L "
          f"{l0.wall_volume*7990:6.3f} kg")
    print(f"   camera            D {2*dd['R_c']*1e3:6.2f} mm  L {dd['L_c']*1e3:6.2f} mm"
          f"   L* {dd['L_star']:6.3f} m")
    print(f"   iniezione         {int(dd['N_inj'])} fori   d_ox {dd['d_ox']*1e3:5.3f} mm"
          f"   d_gpl {dd['d_fuel']*1e3:5.3f} mm   J {dd['momentum_flux_ratio_J']:.3f}")
    print(f"   Dp iniettori      aria {dd['dp_inj_ox']/1e5:6.4f} bar"
          f"   GPL {dd['dp_inj_fuel']/1e5:6.4f} bar")

    print("\n4. FABBRICABILITA' (SLM)")
    note = check_manufacturability(params, None)
    for n in note:
        print(f"   {n}")
    if not note:
        print("   nessuna nota")

    if a.cad:
        print("\n5. IL SOLIDO SI COSTRUISCE DAVVERO?")
        art = build_and_export(params, a.cad, "pareto-ginocchio",
                               BuildOptions(include_injection=True))
        v_formula = l0.wall_volume
        print(f"   STEP e STL scritti in {a.cad}")
        print(f"   solido chiuso (BRep)   {art.is_valid_brep}")
        print(f"   mesh a tenuta          {art.is_watertight_mesh}")
        print(f"   volume CAD {art.volume*1e6:.3f} cm3 contro formula "
              f"{v_formula*1e6:.3f} cm3  (differenza = i fori d'iniezione: "
              f"{(v_formula-art.volume)/v_formula:+.2%})")
        ok &= bool(art.is_valid_brep and art.is_watertight_mesh)

    print("\n" + ("VERIFICA SUPERATA" if ok else "VERIFICA FALLITA"))
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
