#!/bin/bash
Z=/root/zefiro
export PYTHONPATH=$Z/src
for f in sezione.py sezioni.py filmato.py; do
  tr -d '\r' < /mnt/i/AA_ENGINE/scripts/cfd/$f > $Z/scripts/cfd/$f
done
/root/zef/bin/python -c "import matplotlib; print('matplotlib', matplotlib.__version__)" || exit 3
C=$Z/runs/cfd/run3/caso
mkdir -p /mnt/i/AA_ENGINE/runs_cfd
for T in "$@"; do
  /root/zef/bin/python $Z/scripts/cfd/sezione.py --caso "$C" --tempo "$T" \
     --campi C3H8,vorticita,U --out /mnt/i/AA_ENGINE/runs_cfd/sezione_$T.png || exit 4
done
ls -la /mnt/i/AA_ENGINE/runs_cfd/
