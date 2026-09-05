#!/usr/bin/env python3
"""Le pressioni si sono assestate? Con una curva, non con un'impressione.

PERCHE' QUESTO SCRIPT ESISTE. Alla domanda "la CFD gira con le pressioni
giuste?" la corsa 4 non poteva rispondere: le sonde di bordo erano scritte con
`adjustableRunTime` e in 250 microsecondi avevano prodotto DUE campioni, a 20 e
50 us. Due punti non sono una curva. Si sapeva che al termine la portata
uscente valeva 0.0111 kg/s contro 0.0086 entranti - cioe' che il dominio stava
ancora scaricando - ma non si sapeva ne' quanto in fretta stesse finendo ne'
quanto mancasse.

CHE COSA SI GUARDA, E IN CHE ORDINE.

1. LA MASSA CONTENUTA NEL DOMINIO. E' la sonda vera. Il transitorio che si sta
   aspettando E' un riempimento: il dominio parte fermo, gli ingressi ci
   spingono dentro portata, e finche' la massa contenuta cresce (o cala) il
   campo non e' quello di regime. La derivata di questa curva va a zero solo a
   riempimento finito, e non e' una differenza fra numeri grandi.

2. LO SBILANCIO DELLE PORTATE. E' la stessa informazione derivata, ed e' utile
   perche' e' adimensionale e confrontabile con una tolleranza. Ma vicino alla
   convergenza e' la differenza fra due numeri grandi e quasi uguali, quindi e'
   la piu' rumorosa delle due: si legge dopo, non prima.

3. LE PRESSIONI AI BORDI. Sono il risultato che interessa - il salto
   d'iniezione - ma sono anche le piu' lente a sembrare ferme, perche' portano
   sopra il transitorio anche l'oscillazione acustica del dominio. Si guarda la
   DERIVA (media su finestre successive), non il valore istantaneo.

IL CRITERIO SI DICHIARA PRIMA DI GUARDARE I NUMERI, altrimenti si finisce a
sceglierlo per farlo passare:
    massa contenuta stazionaria entro l'1 % su una finestra di 200 us,
    sbilancio di portata sotto il 2 %,
    deriva della pressione di camera sotto l'1 % fra le ultime due finestre.
"""
from __future__ import annotations

import argparse
import pathlib

import numpy as np

#: soglie, dichiarate qui e non sparse nel codice
DERIVA_MASSA = 0.01
SBILANCIO = 0.02
DERIVA_P = 0.01
FINESTRA = 2.0e-4      # s, l'ampiezza su cui si giudica "fermo"


def serie(base: pathlib.Path, nome: str, colonna: int) -> tuple:
    """(t, valore) uniti fra tutti i riavvii, senza doppioni e ordinati.

    Ogni riavvio del solutore crea una cartella `postProcessing/<sonda>/<t0>/`
    nuova, e le finestre si SOVRAPPONGONO quando si riparte da una scrittura
    precedente all'ultimo passo calcolato. Prendere i file in ordine
    alfabetico e concatenarli produrrebbe una serie che torna indietro nel
    tempo; qui si indicizza per istante, e l'ultimo che scrive vince.
    """
    righe: dict[float, float] = {}
    d = base / nome
    if not d.is_dir():
        return np.array([]), np.array([])
    for f in sorted(d.rglob("*.dat")):
        for l in f.read_text().splitlines():
            if l.startswith("#"):
                continue
            c = l.split()
            if len(c) > colonna:
                try:
                    righe[float(c[0])] = float(c[colonna])
                except ValueError:
                    pass
    t = np.array(sorted(righe))
    return t, np.array([righe[x] for x in t])


def finestra(t, v, t0, t1):
    m = (t >= t0) & (t <= t1)
    return v[m] if m.any() else np.array([])


def giudica(base: pathlib.Path) -> int:
    t_m, m = serie(base, "massa_nel_dominio", 1)
    t_pa, pa = serie(base, "p_ingresso_aria", 1)
    t_pg, pg = serie(base, "p_ingresso_gpl", 1)
    t_pu, pu = serie(base, "p_uscita", 1)
    t_ma, ma = serie(base, "portata_ingresso_aria", 1)
    t_mg, mg = serie(base, "portata_ingresso_gpl", 1)
    t_mu, mu = serie(base, "portata_uscita", 1)

    if t_m.size < 10:
        print(f"solo {t_m.size} campioni di massa: troppo presto per giudicare.")
        return 2
    T = t_m[-1]
    print(f"corsa arrivata a {T*1e6:.1f} us   ({t_m.size} campioni di massa, "
          f"{t_pu.size} di pressione)")
    print()

    #: ---- 1. massa contenuta -------------------------------------------- #
    print("--- massa contenuta nel dominio (la sonda del riempimento) ---")
    print(f"{'finestra [us]':>18} {'massa [ug]':>12} {'deriva':>10}")
    bordi = np.arange(0.0, T + FINESTRA, FINESTRA)
    medie = []
    for a, b in zip(bordi[:-1], bordi[1:]):
        w = finestra(t_m, m, a, b)
        if w.size:
            medie.append((0.5 * (a + b), w.mean()))
    for i, (tc, v) in enumerate(medie):
        d = "" if i == 0 else f"{(v - medie[i-1][1]) / medie[i-1][1] * 100:+9.2f}%"
        print(f"{tc*1e6:18.0f} {v*1e9:12.3f} {d:>10}")
    ok_massa = (len(medie) >= 2 and
                abs(medie[-1][1] - medie[-2][1]) / medie[-2][1] < DERIVA_MASSA)
    print(f"  -> massa {'FERMA' if ok_massa else 'ANCORA IN MOVIMENTO'} "
          f"(soglia {DERIVA_MASSA*100:.0f} % su {FINESTRA*1e6:.0f} us)")
    print()

    #: ---- 2. sbilancio di portata ---------------------------------------- #
    print("--- bilancio di massa ai bordi ---")
    ok_bil = False
    if t_mu.size and t_ma.size:
        com = np.intersect1d(np.intersect1d(t_ma, t_mg), t_mu)
        if com.size:
            ia = {x: i for i, x in enumerate(t_ma)}
            ig = {x: i for i, x in enumerate(t_mg)}
            iu = {x: i for i, x in enumerate(t_mu)}
            ent = np.array([abs(ma[ia[x]]) + abs(mg[ig[x]]) for x in com])
            usc = np.array([abs(mu[iu[x]]) for x in com])
            sb = (usc - ent) / ent
            tail = com >= T - FINESTRA
            s = sb[tail] if tail.any() else sb[-20:]
            print(f"  entrante  {ent[-1]*1e3:8.4f} g/s")
            print(f"  uscente   {usc[-1]*1e3:8.4f} g/s")
            print(f"  sbilancio ultimo {sb[-1]*100:+7.2f} %,  "
                  f"medio sull'ultima finestra {s.mean()*100:+7.2f} %")
            ok_bil = abs(s.mean()) < SBILANCIO
    print(f"  -> bilancio {'CHIUSO' if ok_bil else 'APERTO'} "
          f"(soglia {SBILANCIO*100:.0f} %)")
    print()

    #: ---- 3. pressioni ---------------------------------------------------- #
    print("--- pressioni ai bordi [bar] ---")
    print(f"{'finestra [us]':>18} {'aria':>9} {'GPL':>9} {'camera':>9} "
          f"{'salto aria':>11}")
    derive = []
    for a, b in zip(bordi[:-1], bordi[1:]):
        wa, wg, wu = (finestra(t_pa, pa, a, b), finestra(t_pg, pg, a, b),
                      finestra(t_pu, pu, a, b))
        if wa.size and wu.size:
            va, vu = wa.mean(), wu.mean()
            vg = wg.mean() if wg.size else float("nan")
            derive.append(vu)
            print(f"{0.5*(a+b)*1e6:18.0f} {va/1e5:9.4f} {vg/1e5:9.4f} "
                  f"{vu/1e5:9.4f} {(va-vu)/1e5:11.4f}")
    ok_p = (len(derive) >= 2 and
            abs(derive[-1] - derive[-2]) / derive[-2] < DERIVA_P)
    print(f"  -> pressione di camera {'FERMA' if ok_p else 'ANCORA IN DERIVA'} "
          f"(soglia {DERIVA_P*100:.0f} %)")
    print()

    tutto = ok_massa and ok_bil and ok_p
    print("=" * 62)
    if tutto:
        print("VERDETTO: pressioni assestate. Il campo si puo' leggere come")
        print("          regime, e le medie temporali hanno senso.")
    else:
        print("VERDETTO: NON assestato. Le medie e i numeri sul mescolamento")
        print("          riguardano un motore ancora in riempimento.")
        print("          Alzare endTime in system/controlDict (e' rileggibile")
        print("          a caldo) e rimisurare.")
    print("=" * 62)
    return 0 if tutto else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--caso", type=pathlib.Path, required=True)
    a = ap.parse_args()
    return giudica(a.caso / "postProcessing")


if __name__ == "__main__":
    raise SystemExit(main())
