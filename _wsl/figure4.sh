#!/bin/bash
# Le figure della corsa 4, all'ultimo istante calcolato (250 us).
Z=/root/zefiro
export PYTHONPATH=$Z/src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
for f in sezione.py sezioni.py filmato.py motore_intero.py mesh_iniettore.py; do
  tr -d '\r' < /mnt/i/AA_ENGINE/scripts/cfd/$f > $Z/scripts/cfd/$f
done
C=$Z/runs/cfd/run4/caso
O=/mnt/i/AA_ENGINE/runs_cfd
/root/zef/bin/python $Z/scripts/cfd/sezione.py --caso "$C" \
   --campionamento superfici2 --campi C3H8,vorticita,U,T --out $O/run4_sezione.png
/root/zef/bin/python $Z/scripts/cfd/motore_intero.py --caso "$C" \
   --campionamento superfici2 --campi C3H8,vorticita,U --out $O/run4_motore.png
ls -la $O/run4_*.png
