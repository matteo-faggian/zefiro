#!/bin/bash
C=/root/zefiro/runs/cfd/run2/caso
echo "--- lancio.log completo ---"
cat /root/zefiro/runs/cfd/run2/lancio.log
echo "--- ultime 25 righe del solver ---"
tail -25 "$C/log.reactingFoam"
echo "--- andamento del passo temporale ---"
grep "^deltaT" "$C/log.reactingFoam" | awk 'NR%4000==1'
echo "--- Courant max recenti ---"
grep "Courant Number" "$C/log.reactingFoam" | tail -3
echo "--- min/max T recenti ---"
grep "min/max(T)" "$C/log.reactingFoam" | tail -3
