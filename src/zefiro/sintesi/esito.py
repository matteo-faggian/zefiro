"""Che cosa restituisce il modello: un progetto, oppure un rifiuto motivato.

Il rifiuto e' parte del contratto e non un caso d'errore. Un sintetizzatore che
restituisce sempre una geometria e' un sintetizzatore che inventa i dati che
non ha, e la geometria che consegna e' credibile e falsa. Su un modello che
gira dodici ore senza nessuno che guardi, questa e' l'unica difesa che resta.

Un progetto non e' uno STL: e' uno STL PIU' il suo inviluppo di validita' e il
referto di tutte le verifiche che ha passato. "50 N, Isp 143.7 s, valido per
bombola sopra i 15 C, raffica 7.0 s, con questi otto dati ancora da misurare"
e' un progetto. Il solo STL e' un'immagine.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any

from zefiro.sintesi.provenienza import Registro


class Esito(Enum):
    PASSATA = "passata"
    #: La verifica non ha trovato difetti ma non e' conclusiva: tipicamente
    #: perche' la griglia non risolve il dettaglio che doveva controllare.
    #: NON e' un successo, ed e' importante che non lo sembri.
    NON_CONCLUSIVA = "non conclusiva"
    FALLITA = "fallita"


@dataclass(frozen=True)
class Verifica:
    """Una singola verifica, con il suo esito e i numeri che l'hanno prodotto."""
    nome: str
    esito: Esito
    #: Che cosa e' stato controllato, in una riga, in modo che chi legge il
    #: referto capisca che cosa NON e' stato controllato.
    domanda: str
    dettaglio: str = ""
    valori: dict[str, Any] = field(default_factory=dict)
    #: Se la verifica e' stata rifatta a piu' risoluzioni: True solo se il
    #: verdetto e' lo stesso. Un verdetto che cambia col passo non e' un
    #: verdetto.
    convergente: bool | None = None

    @property
    def bloccante(self) -> bool:
        return self.esito is not Esito.PASSATA


@dataclass(frozen=True)
class Inviluppo:
    """Dove il progetto vale, e che cosa succede fuori.

    Non e' documentazione: e' il risultato. Un motore progettato per una
    bombola a 15 C e acceso a 10 C non e' "leggermente diverso", e' un motore
    a cui il combustibile non entra piu' in camera.
    """
    T_bombola_minima: float          # K
    durata_raffica: float            # s
    spinta: float                    # N
    isp: float                       # s
    p_camera: float                  # Pa
    epsilon: float
    #: Frase per frase, che cosa degrada e in che verso uscendo dall'inviluppo.
    fuori_inviluppo: tuple[str, ...] = ()


@dataclass
class Progetto:
    """Un progetto consegnabile."""
    architettura: str
    requisiti: Any
    registro: Registro
    inviluppo: Inviluppo
    verifiche: list[Verifica] = field(default_factory=list)
    #: Oggetti pesanti, presenti solo se la geometria e' stata costruita.
    geometria: Any = None
    percorso_stl: Any = None
    #: Dati che mancano e su cui il progetto sta comunque in piedi, con l'effetto
    #: che avranno quando arriveranno.
    dati_mancanti: list[str] = field(default_factory=list)

    @property
    def consegnabile(self) -> bool:
        """Nessuna verifica bloccante e nessuna quota MANCANTE.

        I `dati_mancanti` dell'impianto NON bloccano: sono misure che l'utente
        deve fare e il cui effetto e' dichiarato. Le quote MANCANTI del
        registro si': quelle sono buchi dentro il progetto stesso.
        """
        return (not any(v.bloccante for v in self.verifiche)
                and not self.registro.mancanti)

    def referto(self) -> str:
        r = [f"PROGETTO  architettura '{self.architettura}'",
             f"  spinta {self.inviluppo.spinta:.2f} N   "
             f"Isp {self.inviluppo.isp:.1f} s   "
             f"p_c {self.inviluppo.p_camera/1e5:.3f} bar   "
             f"eps {self.inviluppo.epsilon:.3f}",
             f"  valido per bombola >= {self.inviluppo.T_bombola_minima-273.15:.1f} C, "
             f"raffica {self.inviluppo.durata_raffica:.2f} s"]
        for f in self.inviluppo.fuori_inviluppo:
            r.append(f"    fuori inviluppo: {f}")
        if self.verifiche:
            r.append("\nVERIFICHE")
            for v in self.verifiche:
                r.append(f"  [{v.esito.value:<15}] {v.nome:<30} {v.dettaglio}")
        if self.dati_mancanti:
            r.append("\nDATI ANCORA DA MISURARE (non bloccanti, effetto dichiarato)")
            for d in self.dati_mancanti:
                r.append(f"  - {d}")
        r.append(self.registro.referto())
        r.append(f"\nCONSEGNABILE: {self.consegnabile}")
        return "\n".join(r)


@dataclass
class Rifiuto:
    """Il modello non produce niente, e dice perche'.

    `dati_mancanti` sono elencati TUTTI insieme: uno alla volta significa far
    ripartire l'utente in laboratorio cinque volte.
    """
    motivo: str
    dati_mancanti: list[str] = field(default_factory=list)
    vincoli_violati: list[str] = field(default_factory=list)
    architetture_provate: list[tuple[str, str]] = field(default_factory=list)

    @property
    def consegnabile(self) -> bool:
        return False

    def referto(self) -> str:
        r = [f"RIFIUTO: {self.motivo}"]
        if self.dati_mancanti:
            r.append("\nDATI MANCANTI (vanno misurati, non stimati)")
            r += [f"  - {d}" for d in self.dati_mancanti]
        if self.vincoli_violati:
            r.append("\nVINCOLI VIOLATI")
            r += [f"  - {v}" for v in self.vincoli_violati]
        if self.architetture_provate:
            r.append("\nARCHITETTURE PROVATE E PERCHE' NON VANNO")
            r += [f"  - {n}: {p}" for n, p in self.architetture_provate]
        return "\n".join(r)
