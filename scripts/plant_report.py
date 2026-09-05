#!/usr/bin/env python3
"""Analisi dell'impianto di alimentazione: che motore permette il tuo banco.

Non usa alcun dato inventato: i limiti del compressore escono dalla potenza
all'albero per via termodinamica, il serbatoio da p, V, T, e le proprieta' del
GPL dalle equazioni di stato di riferimento (CoolProp).

Uso:
    python scripts/plant_report.py
    python scripts/plant_report.py --fad 300 --p-min 6 --bottle-T 20
"""
from __future__ import annotations

import argparse
import math

from zefiro.feed import (
    HP_TO_W,
    autorefrigeration,
    burst_duration,
    composition_from_pressure,
    compressor_bounds,
    mdot_from_fad,
    saturation_pressure,
    tank_blowdown,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--power-hp", type=float, default=3.0)
    ap.add_argument("--tank-litres", type=float, default=100.0)
    ap.add_argument("--p-max", type=float, default=10.0, help="bar")
    ap.add_argument("--p-min", type=float, default=6.0, help="bar")
    ap.add_argument("--T-amb", type=float, default=20.0, help="C")
    ap.add_argument("--fad", type=float, default=None, help="l/min di aria RESA")
    ap.add_argument("--bottle-p", type=float, default=8.0, help="bar assoluti misurati")
    ap.add_argument("--bottle-T", type=float, default=None, help="C misurati")
    args = ap.parse_args()

    T = args.T_amb + 273.15

    print("=" * 74)
    print("COMPRESSORE — limiti termodinamici dalla sola potenza")
    print("=" * 74)
    b = compressor_bounds(args.power_hp * HP_TO_W, 101325.0, args.p_max * 1e5, T)
    print(f"  potenza albero            {b.shaft_power:8.0f} W")
    print(f"  lavoro isotermo ideale    {b.w_isothermal/1e3:8.1f} kJ/kg"
          f"  ->  {b.mdot_isothermal*1e3:6.2f} g/s   (LIMITE SUPERIORE assoluto)")
    print(f"  lavoro adiabatico 1 stadio{b.w_adiabatic/1e3:8.1f} kJ/kg"
          f"  ->  {b.mdot_adiabatic*1e3:6.2f} g/s   (compressore ideale non raffreddato)")
    if args.fad is not None:
        m = mdot_from_fad(args.fad, T)
        print(f"  da targa, FAD {args.fad:.0f} l/min  ->  {m*1e3:6.2f} g/s   <-- USA QUESTO")
        mdot_cont = m
    else:
        print("  targa non fornita: passa --fad <l/min di aria resa>.")
        mdot_cont = None

    print()
    print("=" * 74)
    print(f"SERBATOIO {args.tank_litres:.0f} L — funzionamento a raffica")
    print("=" * 74)
    bd = tank_blowdown(args.tank_litres * 1e-3, args.p_max * 1e5, args.p_min * 1e5, T)
    print(f"  massa iniziale a {args.p_max:.0f} bar   {bd.mass_initial*1e3:7.0f} g")
    print(f"  utilizzabile fino a {args.p_min:.0f} bar {bd.mass_usable_adiabatic*1e3:7.0f} g"
          f"  (adiabatico, T finale {bd.T_final_adiabatic-273.15:.0f} C)")
    print(f"                          {bd.mass_usable_isothermal*1e3:7.0f} g  (isotermo, ottimistico)")
    print()
    rec = mdot_cont or 0.0
    print(f"  {'portata':>10} {'durata (adiab)':>16} {'durata (isot)':>15}")
    for mdot in (0.010, 0.020, 0.050, 0.100, 0.150, 0.200):
        ta, ti = burst_duration(bd, mdot, mdot_recharge=rec)
        fa = "illimitata" if math.isinf(ta) else f"{ta:.1f} s"
        fi = "illimitata" if math.isinf(ti) else f"{ti:.1f} s"
        print(f"  {mdot*1e3:8.0f} g/s {fa:>16} {fi:>15}")

    print()
    print("=" * 74)
    print("BOMBOLA GPL")
    print("=" * 74)
    print(f"  p_sat propano puro a {args.T_amb:.0f} C   "
          f"{saturation_pressure({'C3H8': 1.0}, T)/1e5:5.2f} bar")
    if args.bottle_T is not None:
        Tb = args.bottle_T + 273.15
        x = composition_from_pressure(args.bottle_p * 1e5, Tb)
        print(f"  misura: {args.bottle_p:.2f} bar a {args.bottle_T:.0f} C"
              f"  ->  frazione molare di propano = {x:.3f}")
        if x > 1.0:
            print("    x > 1: la misura non e' compatibile con una binaria propano/butano.")
            print("    Verifica che il manometro sia ASSOLUTO e che la bombola fosse in")
            print("    equilibrio termico. Oppure la temperatura reale e' piu' alta.")
        elif x < 0.0:
            print("    x < 0: pressione troppo bassa anche per butano puro.")
    else:
        print("  Misura pressione E temperatura insieme e rilancia con --bottle-T <C>.")
        for TC in (10, 15, 20, 25, 30):
            x = composition_from_pressure(args.bottle_p * 1e5, TC + 273.15)
            flag = "  (impossibile)" if x > 1.0 else ""
            print(f"    se fossero {TC:2d} C -> x_propano = {x:6.3f}{flag}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
