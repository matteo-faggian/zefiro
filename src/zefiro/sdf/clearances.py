"""Quanto materiale resta fra due vuoti che non devono toccarsi.

E' la domanda che in SLM decide se il pezzo si stampa o perde, e nessuno degli
altri controlli la fa. Chiuso, in un blocco solo, circuito sigillato: tutto
vero anche con un setto da 0.15 mm fra il canale dell'acqua e la camera. Quel
setto in stampa esce poroso, e alla prima accensione l'acqua entra nel gas.

Metodo. Per due insiemi disgiunti A e B, la distanza minima fra loro e'

    min su tutto lo spazio di ( d_A(x) + d_B(x) )

dove d_A e d_B sono le distanze euclidee da A e da B: il minimo si realizza sul
segmento che li congiunge. Si usano le trasformate di distanza esatte di
`scipy.ndimage` sulle maschere booleane, non i campi SDF approssimati: i campi
dei canali nascono da massimi di condizioni e sotto stimano la distanza vicino
agli spigoli, quindi darebbero una risposta prudente ma sbagliata.

La risoluzione della griglia limita la precisione a circa un passo: il rapporto
misura una tolleranza esplicita invece di far finta di avere l'esattezza.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


#: Spessore minimo di parete progettabile in SLM su 316L [m].
#:
#: NON e' un dato verificato: e' il TODO 7. Il valore tipico dichiarato dai
#: costruttori per una parete VERTICALE e' 0.3-0.4 mm; una parete inclinata o
#: che separa due canali in pressione ha bisogno di piu'. 0.5 mm e' la scelta
#: di progetto qui, dichiarata e cambiabile in un posto solo.
MIN_WALL_SLM = 5.0e-4

#: Sotto questo valore il setto non e' una parete sottile: e' un errore di
#: modellazione, due cavita' che si compenetrano.
COMPENETRAZIONE = 1.0e-5


@dataclass(frozen=True)
class Distanza:
    a: str
    b: str
    minima: float          # m
    tolleranza: float      # m, il passo della griglia
    richiesta: float       # m

    @property
    def esito(self) -> str:
        if self.minima <= COMPENETRAZIONE:
            return "COMPENETRANO"
        if self.minima + self.tolleranza < self.richiesta:
            return "SOTTO IL MINIMO"
        if self.minima - self.tolleranza < self.richiesta:
            return "al limite"
        return "ok"


def distanza_minima(maschera_a: np.ndarray, maschera_b: np.ndarray,
                    passo: float) -> float:
    """Distanza euclidea minima fra due insiemi di voxel."""
    from scipy import ndimage

    if not maschera_a.any() or not maschera_b.any():
        return float("inf")
    if (maschera_a & maschera_b).any():
        return 0.0
    da = ndimage.distance_transform_edt(~maschera_a, sampling=passo)
    db = ndimage.distance_transform_edt(~maschera_b, sampling=passo)
    return float((da + db).min())


def _distanze_precalcolate(maschere, passo):
    """Una trasformata di distanza per PARTE, non per coppia.

    Con 18 parti le coppie sono 153, e calcolarne due trasformate ciascuna vuol
    dire 306 passate su venti milioni di voxel: minuti. Le trasformate per
    parte sono 18, e il minimo di ogni coppia si ottiene sommandole. Stesso
    risultato, un ordine di grandezza in meno.
    """
    from scipy import ndimage

    return {k: ndimage.distance_transform_edt(~v, sampling=passo)
            for k, v in maschere.items() if v.any()}


def rapporto_spessori(motore, coppie=None, minimo: float = MIN_WALL_SLM):
    """Spessore di parete fra tutte le coppie di vuoti che contano.

    `coppie` e' una lista di (nome_a, nome_b, spessore_richiesto). Se assente
    si generano tutte le coppie fra le parti, con `minimo` come requisito: e'
    il caso in cui non si sa ancora cosa ci si aspetta e si vuole solo vedere
    se qualcosa si tocca.
    """
    from zefiro.sdf.core import to_numpy

    passo = motore.grid.spacing
    maschere = {k: (to_numpy(v.a) < 0) for k, v in motore.parti.items()}
    nomi = sorted(maschere)
    if coppie is None:
        coppie = [(a, b, minimo) for i, a in enumerate(nomi) for b in nomi[i + 1:]]

    fuori = _maschera_esterno(motore)
    if fuori is not None:
        maschere["esterno"] = fuori

    dist = _distanze_precalcolate(maschere, passo)
    esiti = []
    for a, b, richiesta in coppie:
        if a not in dist or b not in dist:
            continue
        if (maschere[a] & maschere[b]).any():
            d = 0.0
        else:
            d = float((dist[a] + dist[b]).min())
        esiti.append(Distanza(a, b, d, passo, richiesta))
    return esiti


def _maschera_esterno(motore):
    """L'ambiente attorno al pezzo, come insieme di voxel.

    Serve a misurare lo spessore della parete FREDDA - quella fra canale ed
    esterno - che e' la parete che tiene l'acqua in pressione e la cui rottura
    e' il modo piu' banale di far fallire un pezzo raffreddato.
    """
    from scipy import ndimage

    from zefiro.sdf.core import to_numpy

    non_materiale = to_numpy(motore.solido.a) >= 0
    etichette, n = ndimage.label(non_materiale)
    if n == 0:
        return None
    esterno = etichette[0, 0, 0]
    if esterno == 0:
        return None
    # si toglie cio' che comunica con l'esterno attraverso i fori e gli
    # attacchi: quello e' il vuoto del circuito, non l'ambiente
    fuori = etichette == esterno
    for nome, campo in motore.parti.items():
        if nome == "gas":
            continue
        fuori &= ~(to_numpy(campo.a) < 0)
    return fuori
