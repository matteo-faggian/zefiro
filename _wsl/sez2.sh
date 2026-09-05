#!/bin/bash
Z=/root/zefiro
export PYTHONPATH=$Z/src
for f in sezione.py sezioni.py filmato.py; do
  tr -d '\r' < /mnt/i/AA_ENGINE/scripts/cfd/$f > $Z/scripts/cfd/$f
done
C=$Z/runs/cfd/run3/caso
/root/zef/bin/python $Z/scripts/cfd/sezione.py \
   --caso "$C" --campionamento superfici2 --campi C3H8,vorticita,U,T \
   --out /mnt/i/AA_ENGINE/runs_cfd/sezione_finale.png || exit 4
ls -la /mnt/i/AA_ENGINE/runs_cfd/
