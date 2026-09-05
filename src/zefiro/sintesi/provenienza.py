"""Ogni numero del progetto porta con se' da dove viene.

PERCHE' QUESTO MODULO ESISTE, ed e' la cosa piu' importante di tutta la
sintesi. Un modello che produce una geometria produce qualche centinaio di
quote. Alcune sono conseguenze inevitabili della fisica: se la spinta e' 50 N e
la c* e' quella, l'area di gola e' quella e non c'e' niente da discutere.
Altre sono decisioni prese da qualcuno, con una motivazione, che potrebbero
essere diverse. Altre ancora sono dati che nessuno ha misurato e che stanno li'
come segnaposto.

**Guardando lo STL le tre cose sono indistinguibili.** E' esattamente cosi' che
si consegna un pezzo bellissimo e sbagliato: il numero inventato ha lo stesso
aspetto di quello dimostrato.

Qui ogni quota e' una `Scelta` che porta con se' la propria origine e il
proprio motivo. Il referto finale puo' quindi dire, e deve dire, quante quote
sono DERIVATE e quante sono ancora DECISE a mano: e' la misura onesta di
quanto il modello sia un modello e non un progetto travestito da modello.

Non e' burocrazia. E' l'unico modo di rispondere alla domanda "perche' questo
canale e' largo 1.0 mm" sei mesi dopo, quando chi l'ha deciso non se lo ricorda
piu' - e di accorgersi che una quota era giusta per una versione precedente
della geometria, che e' il difetto che in questo progetto e' gia' comparso tre
volte (il raccordo da 1.5 mm, il passo dell'elica, le quote dei circuiti).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any, Iterator


class Origine(Enum):
    """Da dove viene un numero. L'ordine e' quello di affidabilita' crescente
    al contrario: `DERIVATA` non si discute, `MANCANTE` blocca tutto."""

    #: Conseguenza di una legge fisica o di una relazione geometrica, dati i
    #: requisiti. Non si puo' cambiare senza cambiare i requisiti.
    DERIVATA = "derivata"
    #: Valore pubblicato, di norma o di letteratura: una tabella BSPP, un
    #: coefficiente di Holdeman, una tensione di vapore da CoolProp.
    #: Non e' negoziabile ma ha una fonte, e la fonte va scritta.
    NORMATIVA = "normativa"
    #: Misurato sul campo o su una macchina reale.
    MISURATA = "misurata"
    #: Decisa da un umano con una motivazione. E' il numero che un domani
    #: qualcuno rimettera' in discussione, e deve poterlo fare leggendo il
    #: motivo, non indovinandolo.
    DECISA = "decisa"
    #: Segnaposto: nessuno l'ha misurata e il progetto ci sta sopra lo stesso.
    #: Un progetto con anche una sola quota MANCANTE non e' consegnabile.
    MANCANTE = "mancante"


@dataclass(frozen=True)
class Scelta:
    """Un numero con la sua provenienza.

    `valore` e' in unita' SI, sempre, come tutto il resto del progetto.
    `motivo` non e' una descrizione: e' la ragione per cui il numero e' QUEL
    numero. "larghezza del canale" non e' un motivo; "1.0 mm perche' sotto 20
    canali il setto fra canali supera i 120 K di salto" lo e'.
    """
    nome: str
    valore: float
    unita: str
    origine: Origine
    motivo: str
    #: Da dove si puo' rifare il conto: funzione, documento, norma.
    fonte: str = ""
    #: Alternative scartate, con il loro costo. Serve a non rifare due volte la
    #: stessa analisi e a far capire subito quanto margine c'e'.
    scartate: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if not self.motivo.strip():
            raise ValueError(
                f"la quota '{self.nome}' non ha motivo. Una quota senza motivo e' "
                "un numero inventato con un'etichetta sopra: e' proprio il caso "
                "che questo modulo esiste per rendere impossibile.")
        if self.origine is Origine.MANCANTE and self.valore == self.valore:
            # NaN e' l'unico valore ammesso per una quota mancante: qualunque
            # numero finito verrebbe usato dai conti a valle come se fosse vero.
            if self.valore is not None and not _e_nan(self.valore):
                raise ValueError(
                    f"la quota '{self.nome}' e' dichiarata MANCANTE ma ha un valore "
                    f"finito ({self.valore}). Un segnaposto numerico si propaga nei "
                    "conti e produce un risultato credibile e falso: usare float('nan').")

    def __float__(self) -> float:
        return float(self.valore)

    def formatta(self, fattore: float = 1.0, cifre: int = 3) -> str:
        return f"{self.valore * fattore:.{cifre}f} {self.unita}"


def _e_nan(x: Any) -> bool:
    return isinstance(x, float) and x != x


@dataclass
class Registro:
    """L'insieme delle quote di un progetto, con il conto delle origini.

    Il conto e' il punto: `maturita()` dice quale frazione del progetto e'
    conseguenza della fisica e quale e' ancora giudizio umano. Non c'e' un
    valore giusto - un progetto con zero quote DECISE probabilmente non sta
    dicendo la verita' - ma c'e' un modo giusto di guardarlo, ed e' vederlo.
    """
    quote: dict[str, Scelta] = field(default_factory=dict)

    def aggiungi(self, s: Scelta) -> Scelta:
        if s.nome in self.quote and self.quote[s.nome] != s:
            raise ValueError(
                f"la quota '{s.nome}' e' gia' nel registro con un altro valore "
                f"({self.quote[s.nome].valore} contro {s.valore}). Due regole che "
                "producono lo stesso nome sono due regole in conflitto, non un "
                "aggiornamento.")
        self.quote[s.nome] = s
        return s

    def __iter__(self) -> Iterator[Scelta]:
        return iter(self.quote.values())

    def __len__(self) -> int:
        return len(self.quote)

    def __getitem__(self, nome: str) -> Scelta:
        return self.quote[nome]

    def valore(self, nome: str) -> float:
        return float(self.quote[nome].valore)

    def per_origine(self) -> dict[Origine, list[Scelta]]:
        d: dict[Origine, list[Scelta]] = {o: [] for o in Origine}
        for s in self.quote.values():
            d[s.origine].append(s)
        return d

    @property
    def mancanti(self) -> list[Scelta]:
        return [s for s in self.quote.values() if s.origine is Origine.MANCANTE]

    def maturita(self) -> dict[str, float]:
        """Frazione delle quote per origine. E' un indicatore, non un voto."""
        n = max(len(self.quote), 1)
        return {o.value: len(v) / n for o, v in self.per_origine().items()}

    def referto(self, fattori: dict[str, float] | None = None) -> str:
        fattori = fattori or {}
        righe = []
        for o in Origine:
            gruppo = sorted(self.per_origine()[o], key=lambda s: s.nome)
            if not gruppo:
                continue
            righe.append(f"\n{o.value.upper()}  ({len(gruppo)} quote)")
            for s in gruppo:
                f = fattori.get(s.unita, 1.0)
                righe.append(f"  {s.nome:<34}{s.valore * f:12.4f} {s.unita:<6} {s.motivo}")
        m = self.maturita()
        righe.append(
            f"\nMATURITA': {m['derivata']:.0%} derivate, {m['decisa']:.0%} decise, "
            f"{m['normativa']:.0%} da norma, {m['misurata']:.0%} misurate, "
            f"{m['mancante']:.0%} MANCANTI")
        return "\n".join(righe)
