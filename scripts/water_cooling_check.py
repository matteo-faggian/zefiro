#!/usr/bin/env python3
"""Dimensionamento del raffreddamento ad acqua del banco.

Uso:
    python scripts/water_cooling_check.py
    python scripts/water_cooling_check.py --p-water 5 --T-water 12 --gap 0.5 --velocity 12
"""
from __future__ import annotations

import argparse
import math

from zefiro.cooling import (
    channel_pressure_drop,
    coolant_mass_flow,
    dittus_boelter_h,
    reynolds,
    saturation_temperature,
    steady_wall_temperatures,
    zuber_chf,
)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--q-throat", type=float, default=3.81e6, help="W/m2 in gola")
    ap.add_argument("--q-total", type=float, default=15.0e3, help="W totali")
    ap.add_argument("--p-water", type=float, default=4.0, help="bar al collettore")
    ap.add_argument("--T-water", type=float, default=15.0, help="C in ingresso")
    ap.add_argument("--gap", type=float, default=0.5, help="mm di luce del canale")
    ap.add_argument("--velocity", type=float, default=10.0, help="m/s")
    ap.add_argument("--wall", type=float, default=1.5, help="mm di parete in gola")
    ap.add_argument("--k-wall", type=float, default=None,
                    help="W/(m K). Se omesso il calcolo si ferma: e' il TODO n.J.")
    ap.add_argument("--channel-length", type=float, default=40.0, help="mm")
    ap.add_argument("--throat-radius", type=float, default=8.0, help="mm")
    args = ap.parse_args()

    if args.k_wall is None:
        print("ATTENZIONE: conducibilita' della parete non fornita (TODO n.J).")
        print("Uso 16 W/(m K), valore da MANUALE per 316L laminato, che NON e' il")
        print("316L da SLM. Passa --k-wall quando avrai il dato vero.\n")
        k_wall = 16.0
    else:
        k_wall = args.k_wall

    p = args.p_water * 1e5
    T = args.T_water + 273.15
    gap = args.gap * 1e-3
    D_h = 2.0 * gap

    h = dittus_boelter_h(D_h, args.velocity, T, p)
    w = steady_wall_temperatures(args.q_throat, args.wall * 1e-3, k_wall, h, T, p)
    dp = channel_pressure_drop(D_h, args.channel_length * 1e-3, args.velocity, T, p)
    re = reynolds(D_h, args.velocity, T, p)
    portata = 2.0 * math.pi * args.throat_radius * 1e-3 * gap * args.velocity * 1000.0 * 60.0

    print(f"ACQUA  {args.p_water:.1f} bar, {args.T_water:.0f} C   "
          f"T_sat = {saturation_temperature(p) - 273.15:.1f} C")
    print(f"CANALE gap {args.gap:.2f} mm, V {args.velocity:.1f} m/s, D_h {D_h*1e3:.2f} mm")
    print(f"       Re = {re:.0f}  {'(turbolento, Dittus-Boelter valida)' if re > 1e4 else '(SOTTO 10^4: la correlazione non vale)'}")
    print(f"       h = {h:.0f} W/(m2 K)")
    print(f"       perdita di carico nel canale = {dp/1e5:.2f} bar (solo distribuita)")
    print(f"       portata anulare = {portata:.1f} L/min")
    print()
    print(f"PARETE spessore {args.wall:.1f} mm, k = {k_wall:.0f} W/(m K)")
    print(f"       caduta nel film liquido   {w.film_drop:6.0f} K")
    print(f"       caduta per conduzione     {w.conduction_drop:6.0f} K")
    print(f"       T lato acqua = {w.T_coolant_side-273.15:.0f} C   "
          f"margine all'ebollizione {w.margin_to_boiling:.0f} K "
          f"{'OK' if w.margin_to_boiling > 0 else '<-- BOLLE'}")
    print(f"       T lato gas   = {w.T_gas_side-273.15:.0f} C")
    print()
    m = coolant_mass_flow(args.q_total, 20.0, T, p)
    print(f"SMALTIMENTO {args.q_total/1e3:.0f} kW: {m*60:.1f} L/min con salto di 20 K")
    chf = zuber_chf(p)
    print(f"\nCHF di Zuber (pool boiling saturo) = {chf/1e6:.2f} MW/m2")
    print(f"Il flusso di gola e' {args.q_throat/chf:.1f} volte questo valore.")
    print("Il progetto regge SOLO perche' il liquido resta sottoraffreddato e in")
    print("convezione forzata. Se la portata cala, si innesca l'ebollizione a film,")
    print("lo scambio crolla e la parete brucia in pochi secondi.")
    print(">>> Interblocco sulla portata d'acqua che chiude il GPL: NON opzionale.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
