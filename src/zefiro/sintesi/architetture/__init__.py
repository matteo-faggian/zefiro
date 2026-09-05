"""Le famiglie architetturali che il sintetizzatore sa costruire.

Finche' ce n'e' UNA, questo registro e' una promessa e non un'astrazione: si
scopre che cosa era davvero generale solo aggiungendo la seconda. Il fatto e'
scritto qui perche' non venga dimenticato.
"""
from zefiro.sintesi.architetture import aerospike_gas_gas

ARCHITETTURE = (aerospike_gas_gas,)

__all__ = ["ARCHITETTURE", "aerospike_gas_gas"]
