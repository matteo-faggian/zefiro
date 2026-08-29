#!/usr/bin/env python3
"""Costruisce il motore da 50 N in geometria implicita ed esporta lo STL.

Tutti i numeri del raffreddamento vengono da `heat_map_50N.py` e dal
dimensionamento idraulico: sono riportati qui accanto a ciascun campo perche'
si possa risalire al motivo senza cercarlo altrove.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from scipy.optimize import brentq

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from zefiro.geometry.parameters import default_design_vector, derive       # noqa: E402
from zefiro.l0.cycle import evaluate_l0                                    # noqa: E402
from zefiro.schemas import FuelSpec, OperatingPoint                        # noqa: E402
from zefiro.sdf.engine import CircuitoRaffreddamento, costruisci           # noqa: E402
from zefiro.sdf.meshing import (                                          # noqa: E402
    is_watertight, isosurface, mesh_area, mesh_volume, write_stl,
)

P_SUP = 8.0e5
P_C = P_SUP / 1.15          # entrambi i vincoli di Dp si chiudono qui

#: Circuito CAMERA. q = 0.99 MW/m2 su grande area: 2.29 kW, il 62 % del totale.
#: 20 canali perche' sotto quel numero il setto fra canali supera i 120 K di
#: salto (l'aletta conduce lateralmente: dT = q L^2 / 2 k t). 1.0 mm e 3.5 m/s
#: danno h = 12.6 kW/m2K, che bastano: la parete lato acqua sta a 93 C con
#: 50 K di margine all'ebollizione, e la perdita di carico e' 0.1 bar.
CAMERA = CircuitoRaffreddamento(
    nome="camera", n_canali=20, lato=1.0e-3,
    parete_calda=1.2e-3, parete_fredda=1.0e-3,
    #: I due rami NON si sovrappongono. Nella prima versione si accavallavano
    #: fra 29 e 32 mm, e l'incrocio dei due reticoli elicoidali produceva
    #: pareti piu' sottili del passo della griglia: la superficie usciva aperta
    #: (56 spigoli con un solo triangolo). Fra i due resta 2 mm di pieno.
    x_inizio=0.0015, x_fine=0.0265, passo_elica=0.090, velocita=3.5,
)

#: Circuito GOLA + LABBRO. q = 4.65 MW/m2 ma su soli 18 mm di percorso.
#: In PARALLELO al primo, non in serie: e' la ragione per cui la perdita di
#: carico scende da 2.5 a 0.77 bar. Canali piu' piccoli e veloci per h = 42.8
#: kW/m2K, parete calda 0.8 mm per abbassare il salto di conduzione.
GOLA = CircuitoRaffreddamento(
    nome="gola", n_canali=18, lato=0.6e-3,
    parete_calda=0.8e-3, parete_fredda=1.8e-3,
    #: x_fine si ferma 0.5 mm PRIMA del labbro (x_lip = 39.7 mm): oltre, il
    #: mantello non esiste piu' e i canali finirebbero nel vuoto.
    x_inizio=0.0305, x_fine=0.0392, passo_elica=0.035, velocita=12.0,
    #: A U: l'acqua scende al labbro, gira, e torna su uno strato piu' esterno.
    #: Entrambi gli attacchi restano cosi' sulla parte cilindrica, dove un foro
    #: radiale ha senso. Verso il labbro la parete e' conica e nessuna
    #: profondita' di foro raggiunge il collettore senza bucare il gas.
    #: 0.8 + 0.6 + 0.4 + 0.6 + 0.8 = 3.2 mm: sta esattamente nel mantello.
    ritorno=True, setto_ritorno=0.4e-3,
)


def punto_operativo():
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
    return (o, *derive(x, o))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--passo", type=float, default=1.2e-4, help="passo griglia [m]")
    ap.add_argument("--raccordo", type=float, default=1.5e-3)
    ap.add_argument("--out", type=Path, default=Path("runs/motore50"))
    ap.add_argument("--backend", default="numpy")
    a = ap.parse_args()

    from zefiro.sdf.core import backend
    xp = backend(a.backend)

    op, params, l0 = punto_operativo()
    d = params.derived
    print(f"MOTORE 50 N   p_c {P_C/1e5:.2f} bar   spinta {l0.thrust:.2f} N   "
          f"Isp {l0.Isp_s:.1f} s   eps {l0.epsilon:.3f}")
    print(f"  camera D {2*d['R_c']*1e3:.2f} mm  L {d['L_c']*1e3:.1f} mm   "
          f"labbro R {d['R_lip']*1e3:.2f} mm   gola {l0.A_t*1e6:.1f} mm2")

    t0 = time.perf_counter()
    m = costruisci(d, params.plug_contour_x, params.plug_contour_r,
                   [CAMERA, GOLA], a.passo, raccordo=a.raccordo, xp=xp)
    print(f"\ngriglia {m.grid.shape} = {m.grid.n_voxels/1e6:.1f} Mvoxel, "
          f"{m.grid.memoria_mb():.0f} MB per campo, passo {a.passo*1e3:.3f} mm "
          f"({time.perf_counter()-t0:.1f} s)")
    for n in m.note:
        print(f"  ATTENZIONE: {n}")

    v_solido = m.solido.volume()
    v_canali = m.canali.volume()
    print(f"\nVOLUMI  solido {v_solido*1e6:8.3f} cm3   massa 316L "
          f"{v_solido*7990*1e3:6.1f} g")
    print(f"        canali {v_canali*1e6:8.3f} cm3   (cavita' del refrigerante)")

    print("\nesportazione...", flush=True)
    v, f = isosurface(m.solido)
    chiuso = is_watertight(f)
    v_mesh = mesh_volume(v, f)
    a.out.mkdir(parents=True, exist_ok=True)
    sha = write_stl(v, f, a.out / "motore50N.stl")
    print(f"  {len(f)} triangoli   chiuso: {chiuso}")
    print(f"  volume dalla mesh {v_mesh*1e6:.3f} cm3 contro {v_solido*1e6:.3f} "
          f"dal campo  (scarto {abs(v_mesh-v_solido)/v_solido:.2%})")
    print(f"  superficie {mesh_area(v, f)*1e4:.1f} cm2")
    print(f"  scritto {a.out/'motore50N.stl'}  sha {sha[:12]}")

    vc, fc = isosurface(m.canali)
    write_stl(vc, fc, a.out / "canali50N.stl")
    print(f"  scritto {a.out/'canali50N.stl'} ({len(fc)} triangoli): la sola rete "
          "di raffreddamento, per guardarla da sola")

    print("\nRAFFREDDAMENTO (dal dimensionamento, vedi heat_map_50N.py)")
    tot = 0.0
    for c in m.circuiti:
        Q = c.n_canali * c.lato**2 * c.velocita
        tot += Q
        print(f"  {c.nome:<8}{c.n_canali:3d} canali {c.lato*1e3:.1f} mm a "
              f"{c.velocita:4.1f} m/s -> {Q*6e4:5.2f} l/min   "
              f"parete calda {c.parete_calda*1e3:.1f} mm")
    print(f"  {'totale':<8}{tot*6e4:33.2f} l/min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
