"""Verifica topologica della pelle: buchi indesiderati, e come accorgersene.

"E' chiusa?" non basta. Una superficie puo' essere perfettamente chiusa e
avere un buco passante che nessuno ha voluto: un canale che ha sfondato in
camera, una porta che ha attraversato il pezzo, due cavita' che si sono unite.
La superficie resta chiusa, orientabile e senza bordi - e sbagliata.

LO STRUMENTO GIUSTO E' LA CARATTERISTICA DI EULERO.
Per una superficie chiusa e orientabile,

    chi = V - E + F = 2 (1 - g)

dove `g` e' il GENERE, cioe' il numero di manici: quante volte la superficie e'
bucata da parte a parte. Una sfera ha g = 0, una ciambella g = 1, un pezzo con
n fori passanti indipendenti ha g = n.

Il genere si puo' PREVEDERE dal progetto: si conta quanti passaggi passanti si
sono voluti, e si confronta. Se il genere misurato e' piu' alto di quello
previsto, **esiste un buco che nessuno ha disegnato**, e il numero dice quanti.
E' un controllo intero, senza tolleranze: o torna o non torna.

Gli altri controlli qui sono le precondizioni perche' la formula di Eulero
significhi qualcosa: bordi liberi, spigoli non-manifold, vertici non-manifold,
orientazione coerente, triangoli degeneri.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass(frozen=True)
class RapportoTopologia:
    n_vertici: int
    n_spigoli: int
    n_facce: int
    n_componenti: int
    euler: int
    genere: float
    spigoli_di_bordo: int
    spigoli_non_manifold: int
    vertici_non_manifold: int
    facce_degeneri: int
    orientazione_coerente: bool
    genere_atteso: int | None = None
    note: tuple[str, ...] = field(default_factory=tuple)

    @property
    def sana(self) -> bool:
        """Topologia corretta. I triangoli di area nulla NON entrano qui: non
        cambiano V-E+F ne' aprono la superficie, sono un difetto di qualita'
        del file (che alcuni slicer rifiutano), non di geometria. Stanno in
        `pulita` e nelle note."""
        base = (self.spigoli_di_bordo == 0 and self.spigoli_non_manifold == 0
                and self.vertici_non_manifold == 0
                and self.orientazione_coerente)
        if self.genere_atteso is None:
            return base
        return base and abs(self.genere - self.genere_atteso) < 0.5

    @property
    def pulita(self) -> bool:
        """Sana E senza triangoli degeneri: cioe' un file che si puo' mandare
        in stampa senza che nessuno storca il naso."""
        return self.sana and self.facce_degeneri == 0


def analizza(vertici: np.ndarray, facce: np.ndarray,
             genere_atteso: int | None = None,
             area_minima: float = 1.0e-14) -> RapportoTopologia:
    """Tutti i controlli di integrita' della pelle, in una passata."""
    V = int(len(vertici))
    F = int(len(facce))

    spigoli = np.sort(np.vstack([facce[:, [0, 1]], facce[:, [1, 2]],
                                 facce[:, [2, 0]]]), axis=1)
    unici, conteggi = np.unique(spigoli, axis=0, return_counts=True)
    E = int(len(unici))
    bordo = int((conteggi == 1).sum())
    non_manifold_e = int((conteggi > 2).sum())

    # --- vertici non-manifold: due lembi che si toccano in un punto solo ---- #
    # Sono il difetto piu' insidioso, perche' la superficie resta chiusa e con
    # tutti gli spigoli a due facce: solo il ventaglio di triangoli attorno al
    # vertice si spezza in due. Uno slicer li accetta e stampa una clessidra.
    nm_v = _vertici_non_manifold(facce, V)

    # --- triangoli degeneri ------------------------------------------------- #
    a, b, c = vertici[facce[:, 0]], vertici[facce[:, 1]], vertici[facce[:, 2]]
    aree = np.linalg.norm(np.cross(b - a, c - a), axis=1) / 2.0
    degeneri = int((aree <= area_minima).sum())

    # --- orientazione: ogni spigolo interno percorso una volta per verso ---- #
    orientati = np.vstack([facce[:, [0, 1]], facce[:, [1, 2]], facce[:, [2, 0]]])
    coerente = _orientazione_coerente(orientati)

    # --- componenti connesse ------------------------------------------------ #
    n_comp = _componenti(facce, V)

    chi = V - E + F
    # chi = 2 (n_componenti - genere_totale)  =>  g = n_comp - chi/2
    genere = n_comp - chi / 2.0

    note: list[str] = []
    if bordo:
        note.append(f"{bordo} spigoli di bordo: la superficie e' APERTA, "
                    "il genere calcolato non ha significato")
    if non_manifold_e:
        note.append(f"{non_manifold_e} spigoli condivisi da piu' di due facce")
    if nm_v:
        note.append(f"{nm_v} vertici non-manifold: due lembi si toccano in un "
                    "punto solo, e uno slicer ci stampa una clessidra")
    if degeneri:
        note.append(f"{degeneri} triangoli di area nulla")
    if not coerente:
        note.append("orientazione incoerente: alcune facce hanno la normale "
                    "invertita rispetto alle vicine")
    if genere_atteso is not None and abs(genere - genere_atteso) >= 0.5:
        extra = genere - genere_atteso
        note.append(
            f"GENERE {genere:.0f} contro {genere_atteso} previsto: ci sono "
            f"{abs(extra):.0f} passaggi passanti "
            f"{'IN PIU' + chr(39) + ' di quelli disegnati' if extra > 0 else 'IN MENO'}"
        )

    return RapportoTopologia(
        n_vertici=V, n_spigoli=E, n_facce=F, n_componenti=n_comp,
        euler=int(chi), genere=float(genere), spigoli_di_bordo=bordo,
        spigoli_non_manifold=non_manifold_e, vertici_non_manifold=nm_v,
        facce_degeneri=degeneri, orientazione_coerente=coerente,
        genere_atteso=genere_atteso, note=tuple(note),
    )


def _componenti(facce: np.ndarray, n_vertici: int) -> int:
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components

    r = np.concatenate([facce[:, 0], facce[:, 1], facce[:, 2]])
    c = np.concatenate([facce[:, 1], facce[:, 2], facce[:, 0]])
    g = sp.coo_matrix((np.ones(len(r)), (r, c)), shape=(n_vertici, n_vertici))
    n, etichette = connected_components(g, directed=False)
    # si contano solo le componenti che hanno almeno una faccia: i vertici
    # isolati non fanno superficie
    usate = np.unique(etichette[facce[:, 0]])
    return int(len(usate))


def _orientazione_coerente(spigoli_orientati: np.ndarray) -> bool:
    """Ogni spigolo interno deve comparire una volta in un verso e una
    nell'altro. Se compare due volte nello stesso verso, le due facce hanno
    normali opposte e il solido e' rovesciato su se stesso."""
    a = spigoli_orientati
    chiavi = a[:, 0].astype(np.int64) * (a.max() + 1) + a[:, 1]
    _, conte = np.unique(chiavi, return_counts=True)
    return bool((conte == 1).all())


def _vertici_non_manifold(facce: np.ndarray, n_vertici: int) -> int:
    """Un vertice e' manifold se i triangoli che lo contengono formano UN solo
    ventaglio. Si costruisce, per ogni vertice, il grafo degli spigoli opposti
    e si contano le componenti: piu' di una vuol dire due lembi che si toccano
    in quel punto.
    """
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components

    opposti = {}
    for i in range(3):
        v = facce[:, i]
        u = facce[:, (i + 1) % 3]
        w = facce[:, (i + 2) % 3]
        for vv, uu, ww in zip(v, u, w):
            opposti.setdefault(int(vv), []).append((int(uu), int(ww)))

    cattivi = 0
    for v, archi in opposti.items():
        nodi = sorted({x for arco in archi for x in arco})
        if len(nodi) < 3:
            continue
        indice = {n: k for k, n in enumerate(nodi)}
        r = [indice[u] for u, _ in archi]
        c = [indice[w] for _, w in archi]
        g = sp.coo_matrix((np.ones(len(r)), (r, c)), shape=(len(nodi), len(nodi)))
        n, _ = connected_components(g, directed=False)
        if n > 1:
            cattivi += 1
    return cattivi


def salda_e_pulisci(vertici: np.ndarray, facce: np.ndarray,
                    tolleranza: float = 1.0e-9):
    """Salda i vertici coincidenti e toglie i triangoli degeneri.

    Marching cubes produce triangoli di area nulla ogni volta che
    l'isosuperficie passa esattamente per un nodo della griglia - cosa che
    succede sistematicamente sulle facce piane allineate alla griglia. Non
    cambiano la topologia (la caratteristica di Eulero resta giusta), ma alcuni
    slicer ci si strozzano, e uno STL con triangoli nulli e' un file che il
    fornitore puo' rifiutare.

    La riparazione e' quella standard: prima si SALDANO i vertici coincidenti,
    poi si tolgono le facce che hanno due indici uguali. In quest'ordine, e
    solo in quest'ordine, l'operazione non apre buchi: una faccia con due
    vertici nello stesso punto contribuisce uno spigolo percorso avanti e
    indietro, che togliendola sparisce senza lasciare bordo.
    """
    q = np.round(vertici / tolleranza).astype(np.int64)
    _, primo, inverso = np.unique(q, axis=0, return_index=True,
                                  return_inverse=True)
    nuovi = vertici[primo]
    f = inverso[facce]
    buone = (f[:, 0] != f[:, 1]) & (f[:, 1] != f[:, 2]) & (f[:, 2] != f[:, 0])
    f = f[buone]
    usati, f = np.unique(f, return_inverse=True)
    return nuovi[usati], f.reshape(-1, 3)
