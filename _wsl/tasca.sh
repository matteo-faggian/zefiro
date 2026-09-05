#!/bin/bash
Z=/root/zefiro
export PYTHONPATH=$Z/src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
for f in tasca.py sezioni.py filmato.py mesh_iniettore.py; do
  tr -d '\r' < /mnt/i/AA_ENGINE/scripts/cfd/$f > $Z/scripts/cfd/$f
done
/root/zef/bin/python $Z/scripts/cfd/tasca.py --caso $Z/runs/cfd/run4/caso 2>&1 | tail -20
