"""Piani sperimentali deterministici.

Latin Hypercube e non campionamento casuale: a parita' di numero di punti, LHS
garantisce che ogni dimensione sia coperta uniformemente in tutti i suoi strati,
mentre il casuale puro lascia buchi e ammassi. Con 12 dimensioni e poche
centinaia di punti la differenza e' sostanziale.

**Il seme e' obbligatorio, non opzionale.** Un DOE non riproducibile non puo'
entrare nel database delle run (docs/architettura.md sezione 5), quindi non
esiste un default.
"""
from __future__ import annotations

from typing import Iterable, Mapping, Sequence

import numpy as np


def latin_hypercube(
    n_samples: int,
    bounds: Mapping[str, tuple[float, float]],
    seed: int,
    integer_names: Iterable[str] = (),
) -> list[dict[str, float]]:
    """Campione LHS di `n_samples` punti dentro `bounds`.

    Costruzione: ogni dimensione e' divisa in `n_samples` strati di uguale
    probabilita'; si prende un punto a caso in ciascuno strato e si permutano
    gli strati indipendentemente per ogni dimensione. Ne segue che la
    proiezione su OGNI singola dimensione ha esattamente un punto per strato,
    che e' la proprieta' che rende LHS migliore del casuale puro.

    Le dimensioni in `integer_names` vengono arrotondate DOPO il campionamento.
    L'arrotondamento puo' produrre duplicati sulla singola coordinata: e'
    corretto, perche' quella variabile ha davvero meno livelli dello spazio
    continuo.
    """
    if n_samples < 1:
        raise ValueError("n_samples deve essere >= 1")
    names = sorted(bounds)                       # ordine stabile = riproducibilita'
    rng = np.random.Generator(np.random.PCG64(seed))
    d = len(names)

    # un punto per strato, poi permutazione indipendente per dimensione
    strata = (np.arange(n_samples)[:, None] + rng.random((n_samples, d))) / n_samples
    for j in range(d):
        strata[:, j] = rng.permutation(strata[:, j])

    integer_names = set(integer_names)
    out: list[dict[str, float]] = []
    for i in range(n_samples):
        point: dict[str, float] = {}
        for j, name in enumerate(names):
            lo, hi = bounds[name]
            v = lo + strata[i, j] * (hi - lo)
            if name in integer_names:
                v = float(min(max(round(v), np.ceil(lo)), np.floor(hi)))
            point[name] = float(v)
        out.append(point)
    return out


def projection_is_stratified(
    samples: Sequence[Mapping[str, float]],
    bounds: Mapping[str, tuple[float, float]],
    name: str,
) -> bool:
    """Vero se la proiezione su `name` ha esattamente un punto per strato.

    E' la proprieta' che definisce un Latin Hypercube: se salta, il campione
    non e' LHS e va indagato invece di essere usato.
    """
    n = len(samples)
    lo, hi = bounds[name]
    idx = [
        min(int((s[name] - lo) / (hi - lo) * n), n - 1)
        for s in samples
    ]
    return sorted(idx) == list(range(n))
