"""Sintesi: dai requisiti alla geometria, con la provenienza di ogni quota.

L'interfaccia e' una sola funzione:

    progetta(requisiti) -> Progetto | Rifiuto

e il fatto che possa restituire un RIFIUTO e' parte del contratto, non un caso
d'errore. Un modello che restituisce sempre una geometria e' un modello che
inventa i dati che non ha.
"""
from zefiro.sintesi.esito import Progetto, Rifiuto, Verifica
from zefiro.sintesi.provenienza import Origine, Registro, Scelta
from zefiro.sintesi.requisiti import Impianto, Processo, Requisiti
from zefiro.sintesi.sintetizzatore import progetta

__all__ = ["progetta", "Progetto", "Rifiuto", "Verifica", "Origine", "Registro", "Scelta",
           "Impianto", "Processo", "Requisiti"]
