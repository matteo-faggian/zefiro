"""Dal campo alla superficie: marching cubes, verifica, export STL."""
from __future__ import annotations

import hashlib
import struct
from pathlib import Path

import numpy as np

from zefiro.sdf.core import Field, to_numpy


def isosurface(campo: Field, livello: float = 0.0):
    """Triangolazione della superficie a `livello`. Ritorna (vertici, facce).

    I vertici escono in METRI nel sistema della griglia. Marching cubes lavora
    per interpolazione lineare lungo gli spigoli, quindi l'errore di posizione
    e' di ordine h^2 per una superficie liscia, non h: e' il motivo per cui una
    griglia da 0.1 mm da' una superficie molto meglio di 0.1 mm.
    """
    from skimage import measure

    a = to_numpy(campo.a)
    if a.min() > livello or a.max() < livello:
        raise ValueError(
            f"il campo non attraversa il livello {livello}: intervallo "
            f"[{a.min():.4g}, {a.max():.4g}]. Griglia sbagliata o solido vuoto."
        )
    v, f, _, _ = measure.marching_cubes(a, level=livello,
                                        spacing=(campo.grid.spacing,) * 3)
    return v + np.asarray(campo.grid.origin), f


def mesh_volume(vertici: np.ndarray, facce: np.ndarray) -> float:
    """Volume racchiuso [m^3], per il teorema della divergenza.

    Verifica INDIPENDENTE dal conteggio dei voxel: se le due misure
    concordano, l'errore non e' in nessuna delle due.
    """
    a, b, c = vertici[facce[:, 0]], vertici[facce[:, 1]], vertici[facce[:, 2]]
    return float(np.abs(np.einsum("ij,ij->i", a, np.cross(b, c)).sum()) / 6.0)


def mesh_area(vertici: np.ndarray, facce: np.ndarray) -> float:
    a, b, c = vertici[facce[:, 0]], vertici[facce[:, 1]], vertici[facce[:, 2]]
    return float(np.linalg.norm(np.cross(b - a, c - a), axis=1).sum() / 2.0)


def is_watertight(facce: np.ndarray) -> bool:
    """Ogni spigolo appartiene esattamente a due triangoli.

    Controllo topologico, indipendente da marching cubes: se l'algoritmo
    sbagliasse un caso, o se il solido toccasse il bordo della griglia, qui si
    vedrebbe. Uno STL non chiuso non e' stampabile.
    """
    spigoli = np.vstack([facce[:, [0, 1]], facce[:, [1, 2]], facce[:, [2, 0]]])
    spigoli = np.sort(spigoli, axis=1)
    _, conteggi = np.unique(spigoli, axis=0, return_counts=True)
    return bool((conteggi == 2).all())


def write_stl(vertici: np.ndarray, facce: np.ndarray, path: Path,
              scala: float = 1000.0) -> str:
    """STL binario in MILLIMETRI (`scala` 1000 da metri). Ritorna lo sha256.

    Il file e' deterministico a parita' di ingresso: l'intestazione e' fissa e
    non contiene data, cosi' due run identiche danno lo stesso sha256 e il
    versionamento delle run (docs sezione 5) resta valido.
    """
    v = np.asarray(vertici, dtype=np.float64) * scala
    f = np.asarray(facce, dtype=np.int64)
    a, b, c = v[f[:, 0]], v[f[:, 1]], v[f[:, 2]]
    n = np.cross(b - a, c - a)
    ln = np.linalg.norm(n, axis=1, keepdims=True)
    n = np.divide(n, np.where(ln > 0, ln, 1.0))

    dati = bytearray()
    dati += b"zefiro sdf mesh".ljust(80, b"\0")
    dati += struct.pack("<I", len(f))
    blocco = np.zeros((len(f), 12), dtype="<f4")
    blocco[:, 0:3] = n
    blocco[:, 3:6] = a
    blocco[:, 6:9] = b
    blocco[:, 9:12] = c

    grezzo = blocco.tobytes()
    for i in range(len(f)):
        dati += grezzo[i * 48:(i + 1) * 48]
        dati += b"\0\0"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(bytes(dati))
    return hashlib.sha256(bytes(dati)).hexdigest()


def connected_components(vertici: np.ndarray, facce: np.ndarray) -> list[int]:
    """Dimensione (in triangoli) di ogni componente connessa, decrescente.

    Serve perche' una sezione meridiana INGANNA: taglia dove capita, e due
    parti collegate girando attorno all'asse sembrano staccate. L'unica
    risposta affidabile e' topologica. Un pezzo in due componenti non si
    stampa in un colpo solo, e va saputo prima di mandarlo in macchina.
    """
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components as cc

    n = len(vertici)
    righe = np.concatenate([facce[:, 0], facce[:, 1], facce[:, 2]])
    colonne = np.concatenate([facce[:, 1], facce[:, 2], facce[:, 0]])
    g = sp.coo_matrix((np.ones(len(righe)), (righe, colonne)), shape=(n, n))
    n_comp, etichette = cc(g, directed=False)
    per_faccia = etichette[facce[:, 0]]
    return sorted(np.bincount(per_faccia, minlength=n_comp).tolist(), reverse=True)


def cavities_are_closed(campo, canali) -> tuple[bool, float]:
    """La rete di raffreddamento e' un volume CHIUSO dentro il pezzo?

    Ritorna (chiusa, frazione_di_canale_che_sbuca). Un canale che sbuca nella
    camera scarica acqua nel gas: e' il modo piu' rapido di spegnere il motore
    e riempire d'acqua l'impianto. Si controlla che il canale non tocchi la
    cavita' gassosa, non che "sembri" dentro.
    """
    from zefiro.sdf.core import to_numpy

    c = to_numpy(canali.a)
    g = to_numpy(campo.a)
    dentro_canale = c < 0
    if not dentro_canale.any():
        return True, 0.0
    sbuca = dentro_canale & (g < 0)
    return (not sbuca.any()), float(sbuca.sum() / dentro_canale.sum())
