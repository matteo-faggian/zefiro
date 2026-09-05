#!/bin/bash
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
C=/root/zefiro/runs/cfd/run3/caso
cd "$C"
for T in 5e-05 5.25e-05 5.5e-05 5.62272055e-05; do
  echo "########## t = $T"
  mpirun -np 24 postProcess -parallel -time "$T" \
      -func 'fieldMinMax(fields=(U p T k),location=true)' 2>&1 \
    | grep -A2 -iE "^ *(min|max)\(" | head -60
done
