#!/bin/bash
# Fermata PULITA del solutore in stallo: `stopAt writeNow` gli fa scrivere
# l'istante corrente e uscire con "End", cosi' il supervisore vede la corsa
# come conclusa e non la rilancia. Ammazzare mpirun invece lascerebbe il log
# senza "End" e il supervisore ripartirebbe subito da capo.
C=/root/zefiro/runs/cfd/run3/caso
cp "$C/system/controlDict" "$C/system/controlDict.prima_della_fermata"
sed -i 's/^stopAt .*/stopAt              writeNow;/' "$C/system/controlDict"
grep -n "^stopAt" "$C/system/controlDict"
for i in $(seq 1 60); do
  pgrep -f reactingFoam > /dev/null || break
  sleep 2
done
echo "solutore: $(pgrep -c -f reactingFoam) processi"
tail -4 "$C/log.reactingFoam"
echo "--- istanti di volume scritti ---"
ls "$C/processor0" | sort -g
