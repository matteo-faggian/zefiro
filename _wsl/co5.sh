#!/bin/bash
# Alza maxCo a caldo. `runTimeModifiable yes` fa rileggere il controlDict al
# solutore senza riavviarlo, quindi non si perde niente.
#
# PERCHE' SI PUO'. Al momento della modifica, a 55 us: residui iniziali sotto
# 1.2e-3 su tutte le equazioni, tutte risolte in UNA iterazione; errore di
# continuita' locale 1e-10 e cumulato 1.6e-6 fermo; limitatore di temperatura
# mai intervenuto (260-302 K contro una finestra 250-500). Il margine c'e'.
#
# PERCHE' 0.8 E NON DI PIU'. `nOuterCorrectors 1` significa che il PIMPLE sta
# girando in modalita' PISO, e il PISO e' stabile fino a Courant ~1. Oltre
# servirebbero correttori esterni, cioe' cambiare numerica a meta' corsa - la
# cosa che ha fatto piantare tre volte la corsa 3. 0.8 sta sotto il limite e
# raddoppia il passo, che e' quello che serve.
C=/root/zefiro/runs/cfd/run5/caso/system/controlDict
NUOVO="${1:-0.8}"
echo "prima : $(grep 'maxCo' $C)"
sed -i "s/^maxCo .*/maxCo               $NUOVO;/" $C
echo "dopo  : $(grep 'maxCo' $C)"
