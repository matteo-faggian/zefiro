#!/usr/bin/env python3
"""Dimensionamento dell'iniettore del motore da 50 N, e il compromesso che c'e'.

Non stampa "il" risultato: stampa la mappa delle scelte possibili, perche' qui
l'ottimo di mescolamento e il minimo realizzabile in stampa non coincidono e
qualcuno deve decidere dove stare in mezzo, sapendo il prezzo.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from build_50N import punto_operativo                                  # noqa: E402
from zefiro.injector import (                                          # noqa: E402
    C_OTTIMO, D_FORO_MINIMO_STAMPATO, progetta,
)

MW_ARIA, GAMMA_ARIA = 28.9649, 1.400
MW_GPL, GAMMA_GPL = 44.0956, 1.130      # propano; gamma da CoolProp a 283 K, 6.4 bar
CD = 0.75


def main() -> int:
    op, params, l0 = punto_operativo()
    d = params.derived
    p_c = params.free["p_c"]
    print(f"MOTORE 50 N   p_c {p_c/1e5:.3f} bar   aria {l0.mdot_air*1e3:.2f} g/s   "
          f"GPL {l0.mdot_fuel_core*1e3:.3f} g/s   rapporto di massa "
          f"{l0.mdot_air/l0.mdot_fuel_core:.1f}:1")
    print(f"  alimentazione: aria {op.p_air_supply/1e5:.3f} bar, "
          f"GPL {op.p_fuel_supply/1e5:.3f} bar (= p_sat a 15 C)")
    print(f"  camera: raggio {d['R_c']*1e3:.2f} mm, corpo centrale "
          f"{d['r_centerbody']*1e3:.2f} mm, lunghezza {d['L_c']*1e3:.1f} mm")

    def prog(R, n=None, disp="un_lato"):
        return progetta(
            mdot_aria=l0.mdot_air, mdot_gpl=l0.mdot_fuel_core, p_c=p_c,
            p_aria_monte=op.p_air_supply, p_gpl_monte=op.p_fuel_supply,
            T_aria=float(op.T_air_in), T_gpl=float(op.T_fuel_in),
            MW_aria=MW_ARIA, MW_gpl=MW_GPL, gamma_aria=GAMMA_ARIA,
            gamma_gpl=GAMMA_GPL, cd=CD, R_iniezione=R, n_getti=n,
            disposizione=disp)

    base = prog(d["r_centerbody"])
    print("\nCIO' CHE NON E' NEGOZIABILE (lo fissano i salti di pressione)")
    print(f"  aria    V {base.V_aria:6.1f} m/s   rho {base.rho_aria:5.2f}   "
          f"area {base.area_aria*1e6:6.2f} mm2   Dp {base.dp_aria/1e5:.3f} bar")
    print(f"  GPL     V {base.V_getto:6.1f} m/s   rho {base.rho_getto:5.2f}   "
          f"area {math.pi/4*base.n_getti*base.d_getto**2*1e6:6.3f} mm2   "
          f"Dp {base.dp_gpl/1e5:.3f} bar")
    print(f"  J = {base.J:.3f}   e' un RISULTATO: salti uguali in frazione di "
          "p_c e densita' vicine danno J vicino a 1 da soli.")

    print("\nDOVE INIETTARE: il raggio decide quanti getti stanno sulla "
          "circonferenza,\ne quindi il loro diametro a portata fissa.")
    print(f"  {'R_iniez':>8} {'H':>7} {'S':>7} {'n':>4} {'d getto':>9} "
          f"{'C':>6} {'L_mix':>7}")
    print(f"  {'[mm]':>8} {'[mm]':>7} {'[mm]':>7} {'':>4} {'[mm]':>9} "
          f"{'':>6} {'[mm]':>7}")
    for R in (d["r_centerbody"], 0.005, 0.006, 0.008, 0.010, d["R_c"]):
        i = prog(R)
        segno = "  <-- stampabile" if i.stampabile else ""
        print(f"  {R*1e3:8.2f} {i.altezza_anello*1e3:7.3f} {i.passo*1e3:7.3f} "
              f"{i.n_getti:4d} {i.d_getto*1e3:9.3f} {i.C:6.2f} "
              f"{i.lunghezza_mescolamento*1e3:7.2f}{segno}")

    print(f"\nCOMPROMESSO a raggio fisso (dal corpo centrale, R = "
          f"{d['r_centerbody']*1e3:.2f} mm):\n  meno getti = fori piu' grandi "
          "ma C piu' grande dell'ottimo, cioe' getti che\n  attraversano la "
          "corona e vanno a sbattere sulla parete opposta.")
    print(f"  {'n':>4} {'d getto':>9} {'C':>6} {'C/C*':>7} {'penetraz.':>10}")
    for n in (6, 8, 10, 12, 14, 16, 18, 20, 24):
        i = prog(d["r_centerbody"], n=n)
        pen = i.C / i.C_obiettivo
        nota = ""
        if i.d_getto < D_FORO_MINIMO_STAMPATO:
            nota = " foro sotto soglia"
        elif pen > 1.6:
            nota = " sovra-penetrazione"
        elif pen < 0.7:
            nota = " getti attaccati alla parete"
        print(f"  {n:4d} {i.d_getto*1e3:9.3f} {i.C:6.2f} {pen:7.2f} "
              f"{'ok' if 0.7 <= pen <= 1.6 else '--':>10}{nota}")

    print("\nSTRIZIONE DEL CORPO CENTRALE. Il diametro dei getti vale sempre")
    print("  d = 0.386 H  a portata e salti fissati, e H = A_aria / (2 pi R):")
    print("  RIDURRE il raggio d'iniezione ALZA H e quindi il diametro dei fori.")
    print("  Il corpo centrale puo' essere strizzato localmente: e' il solo modo")
    print("  di guadagnare diametro senza toccare nessuna pressione.")
    print(f"  {'R':>6} {'H':>7} {'n':>4} {'d':>7} {'C':>7} {'|C/C*-1|':>9}")
    #: REGOLA DI SCELTA, dichiarata prima di guardare i numeri.
    #:
    #: Fra i punti che la correlazione giudica equivalenti non si sceglie il
    #: minimo di |C/C*-1|: la stessa correlazione di Holdeman ha una dispersione
    #: molto piu' larga dell'1 %, quindi 2.518 e 2.522 sono lo stesso numero e
    #: sceglierne uno perche' e' 0.4 % piu' vicino e' rumore travestito da
    #: ottimizzazione. Si decide invece su due criteri fisici:
    #:
    #:  1. MARGINE SU UN DATO NON MISURATO. Il diametro minimo stampabile e'
    #:     ancora un TODO: chiedere alla macchina esattamente il minimo che
    #:     credo sia il minimo e' progettare sul bordo di cio' che non so.
    #:     Serve almeno il 10 % di margine.
    #:  2. ROBUSTEZZA A UN FORO OTTURATO. Con 4 getti un foro sporco toglie il
    #:     25 % del combustibile e lascia 90 gradi di camera senza GPL; con 6 ne
    #:     toglie il 17 % su 60 gradi. Piu' getti = guasto piu' benigno.
    #:
    #: Quindi: il MASSIMO numero di getti che conserva il 10 % di margine sul
    #: diametro minimo, e fra quelli il raggio con C piu' vicino all'ottimo.
    MARGINE_DIAMETRO = 1.10
    candidati = []
    R = 0.0028
    while R <= 0.00425:
        i = prog(R)
        ok = i.d_getto >= MARGINE_DIAMETRO * D_FORO_MINIMO_STAMPATO
        if ok:
            candidati.append((R, i))
        print(f"  {R*1e3:6.2f} {i.altezza_anello*1e3:7.3f} {i.n_getti:4d} "
              f"{i.d_getto*1e3:7.3f} {i.C:7.3f} {abs(i.scarto_C):9.1%}"
              f"{'' if ok else '  margine insufficiente'}")
        R += 0.00005
    if not candidati:
        raise SystemExit("nessun raggio soddisfa il margine sul diametro minimo")
    n_max = max(c[1].n_getti for c in candidati)
    migliore = min((c for c in candidati if c[1].n_getti == n_max),
                   key=lambda c: abs(c[1].scarto_C))
    R, i = migliore
    print(f"\nSCELTA: R = {R*1e3:.2f} mm, {i.n_getti} getti da "
          f"{i.d_getto*1e3:.3f} mm, C = {i.C:.3f} contro {i.C_obiettivo} "
          f"({i.scarto_C:+.1%})")
    print(f"  altezza del condotto d'aria H = {i.altezza_anello*1e3:.3f} mm "
          f"(da r = {R*1e3:.2f} a {(R+i.altezza_anello)*1e3:.2f} mm)")
    print(f"  il tratto a sezione costante deve essere lungo almeno "
          f"{i.lunghezza_mescolamento*1e3:.2f} mm prima dello sbocco in camera")

    print("\nTEMPI. Il confronto che dice se il progetto ha senso:")
    rho_c = (p_c * l0.MW_c) / (8314.462618 * l0.T_ad)
    A_ann = math.pi * (d["R_c"]**2 - d["r_centerbody"]**2)
    u_c = (l0.mdot_air + l0.mdot_fuel_core) / (rho_c * A_ann)
    t_perm = d["L_c"] / u_c
    t_mix = base.lunghezza_mescolamento / base.V_aria
    print(f"  permanenza in camera   {t_perm*1e3:8.3f} ms  "
          f"(L_c {d['L_c']*1e3:.1f} mm a {u_c:.1f} m/s)")
    print(f"  mescolamento (x/H = 2) {t_mix*1e3:8.3f} ms")
    print(f"  chimica (PSR)            0.03-0.10 ms  (zefiro.l0.chemistry)")
    print(f"  -> il mescolamento consuma il {t_mix/t_perm:.1%} della permanenza")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
