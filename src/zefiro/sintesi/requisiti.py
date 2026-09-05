"""Che cosa si chiede al modello.

Tre gruppi, e la divisione non e' cosmetica:

  * la MISSIONE, cioe' che cosa deve fare il motore (spinta, durata, ambiente);
  * l'IMPIANTO, cioe' con che cosa lo si alimenta - e qui sta la lezione della
    revisione 0.4.0: **la pressione dei propellenti non e' un requisito**. La
    pressione del combustibile e' p_sat(T) e la fissa la termodinamica; quella
    dell'ossidante la fissa il serbatoio e il punto in cui lo si smette di
    scaricare. Chiedere all'utente "a che pressione hai il combustibile" invita
    a scrivere un numero che non esiste, e in questo progetto e' successo;
  * il PROCESSO, cioe' come lo si costruisce. Non e' un vincolo a valle: oggi
    abbiamo visto la taglia di un raccordo allungare la testa d'iniezione di 16
    millimetri. La fabbricazione risale la catena fino alla termodinamica.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Mapping

from zefiro.schemas import MissingDatum


@dataclass(frozen=True)
class Impianto:
    """Con che cosa si alimenta il motore.

    NOTA SULLE PRESSIONI. Qui non si scrive nessuna pressione di combustibile:
    si scrive la temperatura MINIMA a cui la bombola verra' usata, e la
    pressione la calcola `feed.supply_pressures` come p_sat a quella
    temperatura. E' l'unico modo di non poter esprimere una bombola che non
    esiste - il file di configurazione conteneva 8.00 bar, che sono il propano
    puro a 18.3 C, e nessuno se n'era accorto per mesi.
    """
    #: Composizione MOLARE del combustibile (es. {"C3H8": 1.0}).
    combustibile: Mapping[str, float]
    #: Temperatura minima di esercizio della bombola [K]. E' un requisito
    #: OPERATIVO ("sotto questa temperatura non si accende"), non una previsione
    #: meteo, e si garantisce con un bagno d'acqua.
    T_bombola_min: float
    #: Ossidante. Per ora solo aria: il modello e' gas-gas.
    ossidante: str = "aria"
    #: Serbatoio dell'ossidante.
    volume_serbatoio: float = 0.100          # m3
    p_serbatoio_max: float = 10.0e5          # Pa assoluti
    T_serbatoio: float = 293.15              # K
    #: Temperature a VALLE dei riduttori. None = da misurare: l'espansione
    #: attraverso un riduttore e' quasi isentalpica e raffredda il gas, quindi
    #: NON sono la temperatura ambiente.
    T_ossidante_iniezione: float | None = None
    T_combustibile_iniezione: float | None = None
    #: Coefficienti di efflusso degli iniettori. None = da misurare.
    cd_ossidante: float | None = None
    cd_combustibile: float | None = None
    #: Portata massima del riduttore del combustibile [kg/s]. None = da misurare.
    portata_max_riduttore: float | None = None
    #: Meccanismo cinetico per Cantera.
    meccanismo: str = "gri30.yaml"

    def dati_mancanti(self) -> list[str]:
        """Elenco dei dati che mancano, TUTTI insieme.

        Uno alla volta significa far ripartire l'utente in laboratorio cinque
        volte. Il modello li elenca tutti al primo colpo.
        """
        m = []
        if self.T_ossidante_iniezione is None:
            m.append("T_ossidante_iniezione: temperatura dell'aria a valle del "
                     "riduttore, da misurare con un termometro sul tubo")
        if self.T_combustibile_iniezione is None:
            m.append("T_combustibile_iniezione: temperatura del GPL a valle del "
                     "riduttore")
        if self.cd_ossidante is None or self.cd_combustibile is None:
            m.append("cd_ossidante / cd_combustibile: coefficienti di efflusso "
                     "degli iniettori, da prova di soffiaggio a freddo")
        return m


@dataclass(frozen=True)
class Processo:
    """Come si costruisce. Vincoli di fabbricazione, con la loro provenienza.

    Il campo `misurati` dice quali di questi numeri vengono da una prova sulla
    macchina che stampera' davvero il pezzo e quali sono soglie di progetto
    conservative. La differenza conta: progettare esattamente sul minimo che si
    CREDE di conoscere e' progettare sul bordo di cio' che non si sa, ed e' il
    motivo per cui l'iniettore di Zefiro tiene il 10 % di margine sul diametro
    minimo dei fori.
    """
    nome: str = "SLM 316L"
    #: Parete minima che si stampa densa e a tenuta [m].
    parete_minima: float = 5.0e-4
    #: Foro passante minimo ottenibile senza ripresa meccanica [m].
    foro_minimo: float = 6.0e-4
    #: Angolo di autosostentamento del processo [gradi dal piano].
    angolo_autosostentamento: float = 45.0
    #: Margine di progetto sull'angolo: si progetta CON margine sulla soglia di
    #: processo, non uguale ad essa.
    angolo_progetto: float = 52.0
    #: Sottomisura con cui si stampa un foro che verra' forato e maschiato.
    sottomisura_foro: float = 1.2e-3
    #: Densita' del materiale [kg/m3], per la massa.
    densita: float = 7990.0
    #: Quali dei campi sopra vengono da una misura sulla macchina reale.
    misurati: frozenset[str] = frozenset()

    def margine_foro(self, fattore: float = 1.10) -> float:
        """Diametro minimo che si accetta in progetto.

        Se `foro_minimo` non e' misurato si tiene un margine, perche' quel
        numero e' una credenza e non un dato.
        """
        if "foro_minimo" in self.misurati:
            return self.foro_minimo
        return fattore * self.foro_minimo


@dataclass(frozen=True)
class Requisiti:
    """La missione, l'impianto e il processo."""
    spinta: float                  # N, in regime stazionario
    durata: float                  # s, di regime stazionario
    impianto: Impianto
    processo: Processo = field(default_factory=Processo)
    p_ambiente: float = 101325.0   # Pa
    #: Rapporto di equivalenza di progetto in camera. 0.9 = leggermente magro:
    #: la temperatura adiabatica e' quasi al massimo ma si evita di lasciare
    #: combustibile incombusto, che su un banco a cielo aperto e' fuoco che
    #: esce dall'ugello.
    phi: float = 0.9
    #: Frazione minima di p_c richiesta come salto d'iniezione, per lato.
    #: Sotto questa soglia l'iniettore non disaccoppia piu' l'alimentazione
    #: dalla camera e il sistema puo' entrare in chugging.
    frazione_dp_iniezione: float = 0.15
    #: Rapporto d'area minimo perche' l'ugello sia un ugello.
    epsilon_minimo: float = 1.20

    def __post_init__(self) -> None:
        if self.spinta <= 0.0:
            raise MissingDatum("spinta deve essere positiva.")
        if self.durata <= 0.0:
            raise MissingDatum("durata deve essere positiva.")
        if not 0.3 <= self.phi <= 1.5:
            raise MissingDatum(
                f"phi = {self.phi}: fuori da un intervallo in cui una fiamma "
                "propano-aria sta accesa in modo stabile.")
        if self.impianto.ossidante != "aria":
            raise NotImplementedError(
                f"ossidante {self.impianto.ossidante!r}: il modello e' gas-gas e "
                "conosce solo l'aria. Un ossidante liquido richiede un modello "
                "di flash all'iniettore che non esiste (TODO n.3).")

    def dati_mancanti(self) -> list[str]:
        return self.impianto.dati_mancanti()
