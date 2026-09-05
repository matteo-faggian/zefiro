"""Prova di tenuta NUMERICA: quali vuoti comunicano con quali.

Perche' serve un controllo in piu' rispetto a quelli che gia' ci sono.

`is_watertight` dice che la superficie non ha bordi liberi. `topology.analizza`
dice quanti manici ha, e se sono piu' del previsto dice che un buco esiste. Ma
nessuno dei due dice **quale buco**: non sanno che l'acqua non deve toccare il
gas, che il GPL non deve toccare l'aria se non attraverso i fori d'iniezione, e
che niente di tutto cio' deve toccare l'esterno. Un motore con un manico in
piu' e' sbagliato; un motore con il numero giusto di manici ma con l'acqua che
comunica col gas e' sbagliato in modo molto peggiore, e la formula di Eulero da'
lo stesso numero nei due casi.

L'ANALOGO FISICO. In collaudo si tappano tutti gli attacchi tranne uno, si mette
in pressione un circuito alla volta e si guarda se il manometro scende. Questo
modulo fa esattamente quello, sul reticolo di voxel: mette un "gas traccia" in
un punto del circuito, lo lascia diffondere in tutto il vuoto raggiungibile, e
guarda dove arriva.

L'ALGORITMO. Sul complemento del solido si etichettano le componenti connesse
(scipy.ndimage.label, connettivita' a 6 facce: due voxel comunicano solo se
condividono una faccia, mai se si toccano per spigolo o vertice - un passaggio
spesso un voxel in diagonale non e' un passaggio, e' un artefatto di
discretizzazione). Poi:

  * ogni dominio dichiarato deve cadere in UNA componente;
  * due domini che non devono comunicare devono cadere in componenti DIVERSE;
  * due domini che DEVONO comunicare (aria e GPL attraverso i fori) devono
    cadere nella stessa;
  * ogni componente deve essere reclamata da almeno un dominio. Una componente
    che nessuno rivendica e' un vuoto chiuso che nessuno ha disegnato: polvere
    intrappolata se e' interna, un errore di modello altrimenti.

E' un controllo INTERO: o le etichette coincidono o no, senza tolleranze.

LIMITE, dichiarato. Il metodo vede solo cio' che la griglia risolve. Una cricca
piu' sottile del passo non c'e'; una parete piu' sottile del passo puo'
sparire e produrre un falso allarme. Per questo `verifica_tenuta` va rilanciata
almeno a due risoluzioni: **il risultato conta solo se non cambia**. La
funzione `convergenza` fa proprio questo confronto.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage

from zefiro.sdf.core import to_numpy

#: Connettivita' a 6 facce. Vedi il docstring: la connettivita' a 26 farebbe
#: comunicare due canali che si sfiorano per uno spigolo, cioe' una parete
#: perfettamente presente ma spessa meno del passo della griglia.
_STRUTTURA = ndimage.generate_binary_structure(3, 1)


@dataclass(frozen=True)
class Sonda:
    """Un punto che deve stare dentro un certo vuoto.

    `nome` e' il dominio fluido (per es. "acqua_camera"); `punto` sono le
    coordinate in metri. Il punto va scelto in modo che sia certamente dentro
    quel vuoto anche se la geometria cambia un poco: al centro di un canale,
    non contro una parete.
    """
    nome: str
    punto: tuple[float, float, float]


@dataclass(frozen=True)
class RapportoTenuta:
    etichette: dict[str, int]              # dominio -> etichetta della componente
    volumi: dict[int, float]               # etichetta -> volume [m3]
    comunicazioni: list[tuple[str, str]]   # coppie che condividono la componente
    non_rivendicate: list[int]             # componenti che nessun dominio occupa
    volume_non_rivendicato: float
    sonde_fuori: list[str]                 # sonde finite nel SOLIDO
    n_componenti: int
    passo: float
    attese: frozenset[tuple[str, str]] = field(default_factory=frozenset)

    @property
    def inattese(self) -> list[tuple[str, str]]:
        return [c for c in self.comunicazioni if c not in self.attese]

    @property
    def mancate(self) -> list[tuple[str, str]]:
        """Comunicazioni VOLUTE che non ci sono. Un motore in cui il GPL non
        raggiunge la camera e' rotto quanto uno che perde."""
        return [c for c in sorted(self.attese) if c not in self.comunicazioni]

    @property
    def stagno(self) -> bool:
        return (not self.inattese and not self.mancate and not self.sonde_fuori
                and not self.non_rivendicate)

    def riassunto(self) -> str:
        r = [f"componenti di vuoto: {self.n_componenti}   passo {self.passo*1e3:.3f} mm"]
        for nome, et in sorted(self.etichette.items()):
            r.append(f"  {nome:<26} componente {et:>3}  "
                     f"{self.volumi.get(et, 0.0)*1e6:8.3f} cm3")
        for a, b in self.inattese:
            r.append(f"  PERDITA   {a} <-> {b}")
        for a, b in self.mancate:
            r.append(f"  INTERROTTO  {a} non comunica con {b}")
        for s in self.sonde_fuori:
            r.append(f"  SONDA NEL PIENO  {s}")
        if self.non_rivendicate:
            r.append(f"  vuoti non rivendicati: {len(self.non_rivendicate)}, "
                     f"{self.volume_non_rivendicato*1e9:.1f} mm3")
        return "\n".join(r)


def punto_piu_interno(campo, grid, dentro=None) -> tuple[float, float, float]:
    """Il punto piu' lontano da ogni parete DENTRO un vuoto.

    Serve a piazzare le sonde. Scriverne a mano le coordinate significa
    ricalcolarle ogni volta che la geometria cambia di un decimo, e sbagliarle
    in silenzio: una sonda finita nel materiale fa fallire la prova per il
    motivo sbagliato, una sonda finita nel vuoto ACCANTO la fa passare per il
    motivo sbagliato, che e' molto peggio.

    ATTENZIONE AL SEGNO, ed e' facile sbagliarlo. La convenzione del progetto e'
    "negativo dentro il materiale". Un campo che descrive una CAVITA' (un
    canale, la camera, il condotto d'aria) e' quindi negativo DENTRO la
    cavita': il punto piu' interno e' l'argomento del MINIMO, non del massimo.
    Con l'argmax si prende il punto piu' lontano dalla cavita', cioe' uno
    spigolo della griglia in aria libera - e ogni sonda finisce nella stessa
    componente, quella dell'esterno, facendo risultare "tutto comunica con
    tutto". E' successo, e la prova sembrava funzionare.

    SECONDA INSIDIA. Un campo costruito per intersezioni non e' un campo di
    distanza valido lontano dalla superficie: puo' essere molto negativo in
    regioni che non appartengono affatto alla cavita'. Per questo si passa
    `dentro`, una maschera booleana che dice dove il punto ha diritto di
    stare - tipicamente "nel vuoto del pezzo finito e dentro l'involucro".
    Senza, la sonda del ramo di gola finiva in aria libera e il circuito
    risultava comunicante con l'esterno.
    """
    a = to_numpy(campo.a)
    if dentro is None:
        i = np.unravel_index(int(np.argmin(a)), a.shape)
        valore = a[i]
    else:
        m = np.asarray(dentro, dtype=bool)
        if not m.any():
            raise ValueError("la maschera non contiene nessun voxel.")
        b = np.where(m, a, np.inf)
        i = np.unravel_index(int(np.argmin(b)), b.shape)
        valore = a[i]
    if valore >= 0.0:
        raise ValueError(
            "il campo non ha nessun punto interno: il vuoto non esiste a questa "
            "risoluzione. Non e' una sonda mal messa, e' una geometria che la "
            "griglia non risolve.")
    return tuple(float(o + grid.spacing * k) for o, k in zip(grid.origin, i))


def _indice(grid, punto):
    i = [int(round((p - o) / grid.spacing)) for p, o in zip(punto, grid.origin)]
    for k, (v, n) in enumerate(zip(i, grid.shape)):
        if not 0 <= v < n:
            raise ValueError(
                f"sonda fuori dalla griglia sull'asse {k}: indice {v} su {n}. "
                "La griglia non contiene il punto, quindi la prova non "
                "significa niente: allarga il dominio o correggi la sonda."
            )
    return tuple(i)


def verifica_tenuta(
    solido,
    grid,
    sonde: list[Sonda],
    attese: set[tuple[str, str]] | None = None,
) -> RapportoTenuta:
    """Etichetta il vuoto e dice chi comunica con chi.

    `solido` e' il campo di distanza del pezzo (negativo dentro il materiale).
    `sonde` deve contenere anche l'ESTERNO: senza, un foro nella pelle non si
    distingue da un canale regolare.

    `attese` sono le coppie che DEVONO comunicare, come tuple ordinate
    alfabeticamente. Nel motore sono l'aria e il GPL attraverso i fori
    d'iniezione (che e' il modo in cui il motore funziona), e nient'altro.
    """
    attese = frozenset(tuple(sorted(c)) for c in (attese or set()))
    vuoto = to_numpy(solido.a) > 0.0
    etichettato, n = ndimage.label(vuoto, structure=_STRUTTURA)

    dv = grid.spacing ** 3
    conteggi = np.bincount(etichettato.ravel(), minlength=n + 1)
    volumi = {int(k): float(conteggi[k]) * dv for k in range(1, n + 1)}

    etichette: dict[str, int] = {}
    fuori: list[str] = []
    for s in sonde:
        e = int(etichettato[_indice(grid, s.punto)])
        if e == 0:
            fuori.append(s.nome)
        else:
            etichette[s.nome] = e

    comunicazioni = []
    nomi = sorted(etichette)
    for i, a in enumerate(nomi):
        for b in nomi[i + 1:]:
            if etichette[a] == etichette[b]:
                comunicazioni.append((a, b))

    rivendicate = set(etichette.values())
    non_riv = [k for k in volumi if k not in rivendicate]
    return RapportoTenuta(
        etichette=etichette,
        volumi=volumi,
        comunicazioni=comunicazioni,
        non_rivendicate=sorted(non_riv, key=lambda k: -volumi[k]),
        volume_non_rivendicato=sum(volumi[k] for k in non_riv),
        sonde_fuori=fuori,
        n_componenti=n,
        passo=grid.spacing,
        attese=attese,
    )


def confronta(a: RapportoTenuta, b: RapportoTenuta) -> list[str]:
    """Differenze fra due prove alla stessa geometria e risoluzione diversa.

    Il risultato di una prova di tenuta su griglia vale solo se NON dipende dal
    passo. Se cambia, non si e' scoperto un difetto: si e' scoperto che la
    griglia non risolve la geometria, il che e' a sua volta un'informazione
    (quel dettaglio e' al limite anche per la macchina) ma non e' un verdetto.
    """
    d = []
    if a.n_componenti != b.n_componenti:
        d.append(f"numero di vuoti: {a.n_componenti} a {a.passo*1e3:.3f} mm, "
                 f"{b.n_componenti} a {b.passo*1e3:.3f} mm")
    ca, cb = set(a.comunicazioni), set(b.comunicazioni)
    for c in sorted(ca ^ cb):
        dove = f"{a.passo*1e3:.3f}" if c in ca else f"{b.passo*1e3:.3f}"
        d.append(f"comunicazione {c[0]} <-> {c[1]} presente solo a {dove} mm")
    if set(a.etichette) != set(b.etichette):
        d.append("insiemi di sonde diversi: le due prove non sono confrontabili")
    return d
