"""Geometria implicita: campi di distanza con segno su griglia regolare.

PERCHE' NON B-REP
-----------------
`geometry/build.py` costruisce il pezzo come rivoluzione di un poligono in
OCCT. E' esatto e veloce, e non potra' MAI produrre canali di raffreddamento
conformi, collettori raccordati o reticoli: non e' una questione di aggiungere
parametri, e' la rappresentazione sbagliata per quelle forme.

Un campo di distanza con segno (SDF) da' invece, per costruzione:
  * unioni e differenze robuste, senza booleane che falliscono;
  * raccordi (`smooth_union`), che in B-rep sono l'operazione piu' fragile;
  * offset e gusci esatti: `sdf - t` E' la superficie offset di t;
  * canali che seguono una parete curva, che e' esattamente cio' che serve qui.

Il prezzo e' la discretizzazione: la geometria vive su una griglia, e
l'accuratezza e' dell'ordine del passo. Per questo il modulo misura il proprio
errore invece di dichiararlo (vedi `tests/test_sdf.py`).

CPU E GPU
---------
Le operazioni sono scritte su array; `xp` e' numpy oppure cupy. Su un pezzo da
40 x 40 x 70 mm a passo 0.1 mm sono 1.1e8 voxel, cioe' 450 MB per campo in
float32: gestibile in CPU ma lento, immediato su GPU.

ATTENZIONE, DICHIARATA: il percorso CUDA e' scritto ma NON e' stato eseguito
da chi ha scritto questo modulo, perche' l'ambiente di sviluppo non ha una GPU.
Il percorso CPU e' verificato dai test. Chi ha la GPU deve far girare
`tests/test_sdf.py` con ZEFIRO_SDF_BACKEND=cupy prima di fidarsi dei risultati.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Any

import numpy as np


def backend(nome: str | None = None):
    """numpy o cupy. `nome` esplicito, altrimenti ZEFIRO_SDF_BACKEND, altrimenti numpy."""
    scelto = nome or os.environ.get("ZEFIRO_SDF_BACKEND", "numpy")
    if scelto == "numpy":
        return np
    if scelto == "cupy":
        import cupy
        return cupy
    raise ValueError(f"backend {scelto!r} sconosciuto: usa 'numpy' o 'cupy'")


def to_numpy(a):
    """Riporta un array in numpy, qualunque sia il backend."""
    return a.get() if hasattr(a, "get") else np.asarray(a)


@dataclass(frozen=True)
class Grid:
    """Griglia regolare cartesiana. Le coordinate sono in METRI, come tutto
    il resto del progetto: la conversione a mm avviene solo all'export."""

    origin: tuple[float, float, float]
    spacing: float
    shape: tuple[int, int, int]

    @classmethod
    def bounding(cls, lo, hi, spacing: float, margin: float = 0.0) -> "Grid":
        """Griglia che contiene la scatola [lo, hi] piu' un margine.

        Il margine non e' cortesia: marching cubes ha bisogno che il campo
        cambi segno DENTRO la griglia, quindi il solido non deve toccare il
        bordo, altrimenti la superficie viene tagliata e la mesh esce aperta.
        """
        lo = np.asarray(lo, dtype=float) - margin
        hi = np.asarray(hi, dtype=float) + margin
        n = np.maximum(np.ceil((hi - lo) / spacing).astype(int) + 1, 2)
        return cls(tuple(lo.tolist()), float(spacing), tuple(n.tolist()))

    @property
    def n_voxels(self) -> int:
        return int(np.prod(self.shape))

    def memoria_mb(self, n_campi: int = 1) -> float:
        return self.n_voxels * 4 * n_campi / 1024**2

    def axes(self, xp=np):
        """Le tre coordinate come array 1D."""
        return tuple(
            self.origin[k] + self.spacing * xp.arange(self.shape[k], dtype=xp.float32)
            for k in range(3)
        )

    def coords(self, xp=np):
        """Le tre coordinate come array da BROADCAST (non meshgrid pieni).

        La differenza conta: tre meshgrid pieni su 1e8 voxel sono 1.2 GB, tre
        array da broadcast sono qualche kilobyte. Numpy espande al volo.
        """
        ax, ay, az = self.axes(xp)
        return (ax[:, None, None], ay[None, :, None], az[None, None, :])


class Field:
    """Un campo di distanza con segno campionato su una `Grid`.

    Convenzione: **negativo dentro il materiale**, positivo fuori, zero sulla
    superficie. E' la convenzione usuale e quella che rende `f - t` l'offset
    verso l'esterno di `t`.
    """

    __slots__ = ("grid", "a", "xp")

    def __init__(self, grid: Grid, array: Any, xp=None) -> None:
        self.grid = grid
        self.a = array
        self.xp = xp if xp is not None else (np if isinstance(array, np.ndarray)
                                             else backend("cupy"))
        if tuple(array.shape) != tuple(grid.shape):
            raise ValueError(f"array {array.shape} incompatibile con griglia {grid.shape}")

    # -- combinazione ------------------------------------------------------ #
    def union(self, other: "Field") -> "Field":
        return Field(self.grid, self.xp.minimum(self.a, other.a), self.xp)

    def intersection(self, other: "Field") -> "Field":
        return Field(self.grid, self.xp.maximum(self.a, other.a), self.xp)

    def difference(self, other: "Field") -> "Field":
        """Toglie `other` da `self`."""
        return Field(self.grid, self.xp.maximum(self.a, -other.a), self.xp)

    def smooth_union(self, other: "Field", raggio: float) -> "Field":
        """Unione raccordata di raggio `raggio` (polinomiale quadratica).

        E' l'operazione che da' alle forme implicite il loro aspetto: il
        raccordo non e' una feature da costruire, e' il modo in cui due corpi
        si uniscono. In B-rep un raccordo su un'intersezione complessa e'
        l'operazione che fallisce piu' spesso.

        La forma polinomiale (invece dell'esponenziale o del `smin` classico)
        e' scelta perche' e' esattamente la distanza corretta fuori dalla zona
        di raccordo, quindi non degrada il campo lontano dal giunto.
        """
        if raggio <= 0.0:
            return self.union(other)
        xp = self.xp
        h = xp.clip(0.5 + 0.5 * (other.a - self.a) / raggio, 0.0, 1.0)
        mix = other.a * (1 - h) + self.a * h - raggio * h * (1.0 - h)
        return Field(self.grid, mix, xp)

    def smooth_difference(self, other: "Field", raggio: float) -> "Field":
        if raggio <= 0.0:
            return self.difference(other)
        xp = self.xp
        h = xp.clip(0.5 - 0.5 * (other.a + self.a) / raggio, 0.0, 1.0)
        mix = self.a * (1 - h) + (-other.a) * h + raggio * h * (1.0 - h)
        return Field(self.grid, mix, xp)

    # -- forma ------------------------------------------------------------- #
    def offset(self, t: float) -> "Field":
        """Superficie spostata di `t` verso l'esterno. Esatta per costruzione:
        e' il significato stesso di un campo di distanza."""
        return Field(self.grid, self.a - t, self.xp)

    def shell(self, spessore: float, verso: str = "interno") -> "Field":
        """Guscio di spessore costante a partire dalla superficie."""
        if verso == "interno":
            return self.offset(0.0).intersection(
                Field(self.grid, -(self.a + spessore), self.xp))
        if verso == "esterno":
            return self.offset(spessore).intersection(
                Field(self.grid, -self.a, self.xp))
        raise ValueError("verso: 'interno' o 'esterno'")

    # -- misura ------------------------------------------------------------- #
    def volume(self) -> float:
        """Volume del solido [m^3], con correzione sub-voxel.

        Contare i voxel negativi darebbe un errore di ordine h (il volume del
        guscio di voxel a cavallo della superficie). Si usa invece una
        frazione lineare basata sul valore del campo: ogni voxel contribuisce
        `clip(0.5 - phi/h, 0, 1)`, che e' esatto al primo ordine per una
        superficie piana e riduce l'errore a ordine h^2. Il test lo misura.
        """
        h = self.grid.spacing
        frazione = self.xp.clip(0.5 - self.a / h, 0.0, 1.0)
        return float(to_numpy(frazione.sum())) * h**3

    def sample(self, punti) -> np.ndarray:
        """Valore del campo in punti arbitrari, per interpolazione trilineare.

        Serve a verificare il campo campionato CONTRO la formula analitica in
        posizioni che non cadono sui nodi: e' il modo di misurare l'errore di
        discretizzazione invece di supporlo.
        """
        a = to_numpy(self.a)
        g = self.grid
        p = (np.asarray(punti, dtype=float) - np.asarray(g.origin)) / g.spacing
        i0 = np.floor(p).astype(int)
        f = p - i0
        i0 = np.clip(i0, 0, np.array(g.shape) - 2)
        out = np.zeros(len(p))
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    w = ((f[:, 0] if dx else 1 - f[:, 0])
                         * (f[:, 1] if dy else 1 - f[:, 1])
                         * (f[:, 2] if dz else 1 - f[:, 2]))
                    out += w * a[i0[:, 0] + dx, i0[:, 1] + dy, i0[:, 2] + dz]
        return out

    def tocca_il_bordo(self, margine: int = 1) -> bool:
        """True se il solido arriva al bordo della griglia: in quel caso la
        superficie viene tagliata e la mesh esce APERTA. Meglio saperlo prima
        di esportare che dopo aver mandato in stampa un STL non chiuso."""
        a = to_numpy(self.a)
        m = margine
        facce = (a[:m], a[-m:], a[:, :m], a[:, -m:], a[:, :, :m], a[:, :, -m:])
        return bool(any(float(f.min()) <= 0.0 for f in facce))
