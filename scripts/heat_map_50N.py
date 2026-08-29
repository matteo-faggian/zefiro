#!/usr/bin/env python3
"""Dove va il calore, lungo tutto il percorso del gas.

Serve a decidere DOVE mettere i canali di raffreddamento invece di supporlo.
Su un aerospike a espansione esterna la risposta non e' ovvia: la gola non e'
sulla parete di camera - e' sul plug, con il labbro come sola frontiera
esterna. Se il carico e' li', i canali in camera raffreddano dove non serve,
che e' esattamente l'errore gia' fatto una volta col film cooling in testa.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from zefiro.geometry.parameters import default_design_vector, derive   # noqa: E402
from zefiro.l0.cycle import evaluate_l0                                # noqa: E402
from zefiro.l0.mixture import MixtureModel                             # noqa: E402
from zefiro.schemas import FuelSpec, OperatingPoint                    # noqa: E402
from zefiro.thermal import adiabatic_wall_temperature, bartz_h_gas     # noqa: E402

P_SUP = 8.0e5
P_C = P_SUP / 1.15
T_WALL = 800.0


def punto_50N():
    fuel = FuelSpec(composition={"C3H8": 1.0}, phase_at_injection="gas",
                    thermo_source="gri30.yaml")
    def op(m):
        return OperatingPoint(p_amb=101325.0, p_air_supply=P_SUP, p_fuel_supply=8.0e5,
                              fuel=fuel, T_air_in=293.0, T_fuel_in=283.0,
                              mdot_air_max=m, cd_injector_ox=0.75, cd_injector_fuel=0.75)
    m = brentq(lambda m: evaluate_l0(op(m), p_c=P_C, phi_core=0.9, f_film=0.0,
                                     thermal_severity=False).thrust - 50.0,
               1e-3, 0.2, xtol=1e-10)
    o = op(m)
    x = default_design_vector(p_c=P_C, phi_core=0.9, f_film=0.0, N_inj=12.0,
                              Dc_over_Dt=2.6, Lc_over_Dc=1.4, plug_trunc=0.8,
                              t_wall=0.0012, d_ox_ratio=0.05, fuel_vel_ratio=1.1,
                              conv_half_angle=0.65)
    params, l0 = derive(x, o)
    return o, x, params, l0


def mappa_di_calore(op, params, l0, n=140):
    """q(x) lungo la parete di camera e lungo il plug.

    L'area di passaggio locale e' quella ANULARE fra parete e corpo centrale,
    non pi r^2: e' la differenza fra una camera anulare e un tubo, e sbagliarla
    sposta il rapporto di aree quindi h.
    """

    d = params.derived
    model = MixtureModel.from_fuel(op.fuel)
    g = model.gas
    g.TPX = l0.T_ad, P_C, l0.X_eq
    g.equilibrate("TP")
    mu, cp = g.viscosity, g.cp_mass
    Pr = mu * cp / g.thermal_conductivity
    gamma, A_t = l0.gamma_c, l0.A_t
    D_t = 2.0 * math.sqrt(A_t / math.pi)

    R_c, R_lip, L_c, L_conv = d["R_c"], d["R_lip"], d["L_c"], d["L_conv"]
    r_cb = d["r_centerbody"]
    x_lip = L_c + L_conv
    cx = np.asarray(params.plug_contour_x)
    cr = np.asarray(params.plug_contour_r)

    stazioni = []                       # (x, r_parete, A_locale, zona)
    for xx in np.linspace(0.0, L_c, n // 3):
        stazioni.append((xx, R_c, math.pi * (R_c**2 - r_cb**2), "camera"))
    for xx in np.linspace(L_c, x_lip, n // 3):
        t = (xx - L_c) / max(L_conv, 1e-12)
        r_out = R_c + t * (R_lip - R_c)
        stazioni.append((xx, r_out, math.pi * (r_out**2 - r_cb**2), "convergente"))
    for i in range(len(cx)):
        # sul plug l'area di passaggio e' quella normale al flusso: si usa il
        # rapporto d'aree isentropico gia' coerente col contorno di Angelino
        r = cr[i]
        A = math.pi * (R_lip**2 - r**2) if r < R_lip else A_t
        stazioni.append((x_lip + cx[i], r, max(A, A_t), "plug"))

    righe = []
    for xx, rw, A, zona in stazioni:
        eps = max(A / A_t, 1.0)
        M = mach_da_area(eps, gamma, supersonico=(zona == "plug"))
        T_aw = adiabatic_wall_temperature(l0.T_ad, M, gamma, Pr)
        h = bartz_h_gas(D_t, P_C, l0.c_star, mu, cp, Pr, gamma, eps, M,
                        T_WALL, l0.T_ad)
        righe.append((xx, rw, zona, eps, M, T_aw, h, h * (T_aw - T_WALL)))
    return righe


def mach_da_area(eps, gamma, supersonico):
    """Mach dal rapporto d'aree isentropico. Due radici: si sceglie il ramo."""
    if eps <= 1.0 + 1e-12:
        return 1.0
    g = gamma
    def f(M):
        return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** (
            (g + 1.0) / (2.0 * (g - 1.0))) - eps
    return brentq(f, 1.0000001, 12.0) if supersonico else brentq(f, 1e-4, 0.9999999)


def main() -> int:
    op, x, params, l0 = punto_50N()
    print(f"MOTORE 50 N   p_c {P_C/1e5:.2f} bar   spinta {l0.thrust:.2f} N   "
          f"Isp {l0.Isp_s:.1f} s   T_ad {l0.T_ad:.0f} K")
    righe = mappa_di_calore(op, params, l0)
    q = np.array([r[7] for r in righe])
    print(f"\n{'zona':<13}{'x [mm]':>9}{'r [mm]':>9}{'A/A_t':>8}{'M':>7}"
          f"{'T_aw [K]':>10}{'h [kW/m2K]':>12}{'q [MW/m2]':>11}")
    for zona in ("camera", "convergente", "plug"):
        sel = [r for r in righe if r[2] == zona]
        for r in (sel[0], sel[len(sel) // 2], sel[-1]):
            print(f"{zona:<13}{r[0]*1e3:9.2f}{r[1]*1e3:9.2f}{r[3]:8.2f}{r[4]:7.3f}"
                  f"{r[5]:10.0f}{r[6]/1e3:12.1f}{r[7]/1e6:11.3f}")
    print()
    for zona in ("camera", "convergente", "plug"):
        sel = np.array([r[7] for r in righe if r[2] == zona])
        print(f"  {zona:<13} q da {sel.min()/1e6:6.3f} a {sel.max()/1e6:6.3f} MW/m2")
    i = int(np.argmax(q))
    print(f"\nMASSIMO: {q[i]/1e6:.3f} MW/m2 in zona '{righe[i][2]}' a "
          f"x = {righe[i][0]*1e3:.2f} mm, r = {righe[i][1]*1e3:.2f} mm")

    # potenza totale, per zona: e' quella che l'acqua deve portare via
    print("\nPOTENZA TERMICA PER ZONA (integrale di q sulla superficie bagnata)")
    tot = 0.0
    for zona in ("camera", "convergente", "plug"):
        sel = [r for r in righe if r[2] == zona]
        P = 0.0
        for a, b in zip(sel[:-1], sel[1:]):
            slant = math.hypot(b[0] - a[0], b[1] - a[1])
            P += 0.5 * (a[7] + b[7]) * math.pi * (a[1] + b[1]) * slant
        tot += P
        print(f"  {zona:<13}{P/1e3:8.2f} kW")
    print(f"  {'TOTALE':<13}{tot/1e3:8.2f} kW")
    dT = 30.0
    print(f"\n  portata d'acqua per un salto di {dT:.0f} K: "
          f"{tot/(4180.0*dT)*1e3:.0f} g/s = {tot/(4180.0*dT)*60:.1f} l/min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
