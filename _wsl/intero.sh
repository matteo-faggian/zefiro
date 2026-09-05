#!/bin/bash
Z=/root/zefiro
export PYTHONPATH=$Z/src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
for f in sezione.py sezioni.py filmato.py motore_intero.py mesh_iniettore.py; do
  tr -d '\r' < /mnt/i/AA_ENGINE/scripts/cfd/$f > $Z/scripts/cfd/$f
done
C=$Z/runs/cfd/run3/caso
/root/zef/bin/python $Z/scripts/cfd/motore_intero.py \
   --caso "$C" --campionamento superfici2 --campi C3H8,vorticita,U \
   --out /mnt/i/AA_ENGINE/runs_cfd/motore_intero.png
