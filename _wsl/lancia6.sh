#!/bin/bash
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
C=/root/zefiro/runs/cfd/run6/caso
[ -d "$C/processor0" ] || { echo "caso non decomposto"; exit 2; }
cd "$C"
for t in $(seq 1 12); do
  echo "=== avvio $t === $(date)"
  mpirun -np 24 reactingFoam -parallel >> log.reactingFoam 2>&1
  echo "uscita $?  $(grep '^Time = ' log.reactingFoam | tail -1)  $(date)"
  if grep -q "^End$" log.reactingFoam; then echo FINITO; break; fi
  sleep 5
done
