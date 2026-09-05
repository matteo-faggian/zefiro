#!/usr/bin/env python3
"""Il campo e' turbolento, e la turbolenza mescola? Due domande diverse.

PERCHE' NON BASTA DIRE "IL FLUSSO E' TURBOLENTO". Il numero di Reynolds dice
che il moto NON e' laminare; non dice che la turbolenza faccia il lavoro che
gli si sta chiedendo. La correlazione di Holdeman - C = (S/H) sqrt(J), su cui
e' dimensionato tutto l'iniettore - e' una correlazione turbolenta: presuppone
che il getto trasversale generi la coppia di vortici controrotanti e che quei
vortici, nel tempo che il fluido passa nel condotto, uniformino la miscela.
Questo modulo misura i tre numeri che dicono se quella presupposizione regge:

1. Re, per sapere se c'e' turbolenza da usare.

2. nu_t / nu, cioe' di quanto la turbolenza moltiplica la diffusione
   molecolare. E' il guadagno di mescolamento, e senza di lui il numero 3 non
   si puo' nemmeno calcolare.

3. IL CONFRONTO CHE DECIDE: tempo di diffusione turbolenta attraverso il vuoto
   fra due getti, contro tempo di permanenza.

       t_diff = L^2 / nu_t        (L = mezzo passo azimutale fra due getti)
       t_perm = x / U             (quanto dura il viaggio fino alla misura)

   Se t_diff >> t_perm la turbolenza NON fa in tempo a chiudere il vuoto, e
   allora la disuniformita' misurata non e' un difetto del calcolo: e' quello
   che il motore fa. E' esattamente il conto che, sulla corsa 4, dava 2441-7372
   us contro i 29 disponibili - cioe' due ordini di grandezza di scarto.

4. RISOLTA CONTRO MODELLATA. In URANS convivono due turbolenze: k, che il
   modello porta, e le fluttuazioni che il campo medio compie nel tempo, che il
   modello non sa di avere. Se la seconda e' trascurabile il calcolo e' di
   fatto stazionario e un'istantanea vale quanto una media; se non lo e', ogni
   numero letto su un'istantanea e' un numero letto a caso nel ciclo. Serve
   `UPrime2Mean`, che esiste solo dopo che `fieldAverage` ha accumulato: si
   analizza se c'e', altrimenti lo si dice.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from sezioni import leggi_vtp  # noqa: E402

#: viscosita' cinematica dell'aria alle condizioni di camera, da Sutherland.
#: Non e' un valore tabulato a memoria: si ricalcola sotto da mu(T)/rho.
SUTHERLAND_AS = 1.4792e-06
SUTHERLAND_TS = 116.0


def nu_molecolare(T: float, rho: float) -> float:
    mu = SUTHERLAND_AS * T ** 1.5 / (T + SUTHERLAND_TS)
    return mu / rho


def pesa(v: np.ndarray, w: np.ndarray) -> float:
    """Media pesata, con guardia: un peso nullo ovunque non e' una media."""
    s = w.sum()
    return float((v * w).sum() / s) if s > 0 else float("nan")


def analizza(caso: Path, campionamento: str, tempo: str | None) -> int:
    from mesh_iniettore import quote_dal_progetto

    d = quote_dal_progetto()
    base = caso / "postProcessing" / campionamento
    tempi = sorted((p for p in base.iterdir() if p.is_dir()),
                   key=lambda p: float(p.name))
    t = tempi[-1] if tempo is None else next(p for p in tempi if p.name == tempo)
    ang = math.pi / d["n_getti"]
    piani = (("meridiano_getto", "piano del getto"),
             ("meridiano_fianco", "fra due getti"))

    print(f"istante {float(t.name)*1e6:.1f} us   "
          f"(campionamento '{campionamento}')")
    print()

    #: --- 1. Reynolds, dai dati di progetto ------------------------------- #
    nu_a = nu_molecolare(d["T_aria"], d["rho_aria"])
    nu_g = nu_molecolare(d["T_gpl"], d["rho_gpl"])
    Re_anello = d["V_aria"] * 2.0 * d["H"] / nu_a
    Re_getto = d["V_gpl"] * d["d_getto"] / nu_g
    print("--- 1. c'e' turbolenza da usare? ---")
    print(f"  nu aria {nu_a:.3e} m2/s   nu GPL {nu_g:.3e} m2/s")
    print(f"  Re condotto (2H = {2*d['H']*1e3:.2f} mm) : {Re_anello:8.0f}")
    print(f"  Re getto    (d  = {d['d_getto']*1e3:.3f} mm): {Re_getto:8.0f}")
    print(f"  -> {'turbolento' if Re_anello > 4000 else 'NON pienamente turbolento'}"
          " nel condotto, "
          f"{'turbolento' if Re_getto > 4000 else 'transizionale'} nel getto")
    print()

    #: --- 2 e 3. dai campi campionati ------------------------------------- #
    print("--- 2. quanto mescola la turbolenza calcolata ---")
    print(f"{'piano':>16} {'k medio':>10} {'nut medio':>11} {'nut/nu':>9} "
          f"{'nut max':>11}")
    nut_medi = {}
    for nome, etichetta in piani:
        f = t / f"{nome}.vtp"
        if not f.exists():
            print(f"{etichetta:>16}   piano assente")
            continue
        punti, _, campi = leggi_vtp(f)
        if "nut" not in campi or "k" not in campi:
            print(f"{etichetta:>16}   il campionamento non porta k/nut: "
                  "e' una corsa precedente alla 5.")
            return 2
        #: si pesa sulla densita': una media aritmetica su nodi non uniformi
        #: peserebbe di piu' dove la mesh e' fitta, cioe' proprio attorno al
        #: getto, e restituirebbe la turbolenza del getto spacciata per quella
        #: del condotto.
        w = campi.get("rho", np.ones(len(punti)))
        k, nut = campi["k"], campi["nut"]
        nut_medi[nome] = pesa(nut, w)
        print(f"{etichetta:>16} {pesa(k, w):10.2f} {nut_medi[nome]:11.3e} "
              f"{nut_medi[nome]/nu_a:9.1f} {float(nut.max()):11.3e}")
    print()

    #: --- 3. il confronto che decide -------------------------------------- #
    print("--- 3. la turbolenza fa in tempo a chiudere il vuoto azimutale? ---")
    #: mezzo passo azimutale al raggio d'iniezione: e' la distanza che il
    #: combustibile deve percorrere DI LATO per riempire il vuoto fra due getti.
    S = 2.0 * math.pi * d["R_getti"] / d["n_getti"]
    L = 0.5 * S
    x_mis = -d["x_getti"]          # dal getto al piano di misura x = 0
    t_perm = x_mis / d["V_aria"]
    print(f"  passo azimutale S = {S*1e3:.3f} mm, mezzo passo L = {L*1e3:.3f} mm")
    print(f"  tratto getto -> misura = {x_mis*1e3:.3f} mm a {d['V_aria']:.0f} m/s")
    print(f"  tempo di permanenza t_perm = {t_perm*1e6:.1f} us")
    for nome, etichetta in piani:
        if nome not in nut_medi or not np.isfinite(nut_medi[nome]):
            continue
        t_diff = L * L / nut_medi[nome]
        print(f"  {etichetta:>16}: t_diff = {t_diff*1e6:9.1f} us  "
              f"-> t_diff/t_perm = {t_diff/t_perm:7.1f}")
    print()
    print("  Se il rapporto e' >> 1 la diffusione turbolenta NON chiude il")
    print("  vuoto fra i getti nel tratto disponibile: la disuniformita'")
    print("  misurata e' del motore, non del calcolo.")
    print()

    #: --- 4. risolta contro modellata ------------------------------------- #
    print("--- 4. il campo URANS e' fermo o sta oscillando? ---")
    fatto = False
    for nome, etichetta in piani:
        f = t / f"{nome}.vtp"
        if not f.exists():
            continue
        _, _, campi = leggi_vtp(f)
        if "UPrime2Mean" not in campi:
            continue
        fatto = True
        w = campi.get("rho", None)
        p2 = campi["UPrime2Mean"]
        #: la traccia del tensore e' 2k_risolta: le prime tre colonne di
        #: symmTensor sono xx, xy, xz -> la traccia e' 0, 3, 5.
        k_ris = 0.5 * (p2[:, 0] + p2[:, 3] + p2[:, 5])
        k_mod = campi["k"]
        w = w if w is not None else np.ones(len(k_mod))
        km, kr = pesa(k_mod, w), pesa(k_ris, w)
        print(f"  {etichetta:>16}: k modellata {km:8.2f}  risolta {kr:8.2f}  "
              f"risolta/totale {kr/(kr+km)*100:5.1f} %")
    if not fatto:
        print("  `UPrime2Mean` non e' fra i campi campionati: le medie non")
        print("  sono ancora state accumulate, oppure i piani non le portano.")
        print("  Si ricampiona a posteriori con `postProcess -func superfici`")
        print("  aggiungendo i campi medi all'elenco.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", type=Path, required=True)
    ap.add_argument("--campionamento", default="superfici")
    ap.add_argument("--tempo", default=None)
    a = ap.parse_args()
    return analizza(a.caso, a.campionamento, a.tempo)


if __name__ == "__main__":
    raise SystemExit(main())
