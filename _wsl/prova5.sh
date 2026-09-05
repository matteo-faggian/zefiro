#!/bin/bash
# Sincronizza, prova a generare SOLO il caso (senza mesh) e controlla che i
# blocchi nuovi ci siano davvero nel controlDict prodotto.
set -e
bash /mnt/i/AA_ENGINE/_wsl/run.sh sync5 >/dev/null 2>&1 || true
cd /root/zefiro
export PYTHONPATH=/root/zefiro/src
/root/zef/bin/python -m pip install -q pytest 2>&1 | tail -2 || echo "(pip non disponibile)"
rm -rf /tmp/prova5
/root/zef/bin/python scripts/cfd/caso_mescolamento.py \
    --caso /tmp/prova5 --solo-caso --tempo 1.0e-3 --scrittura 2.0e-5 \
    --campionamento 2.0e-6 --proc 24 --binario \
    --passo-sonde 20 --media-da 5.0e-4 2>&1 | tail -12
echo
echo "=== blocchi nel controlDict ==="
grep -nE '^\s{4}[a-z_]+$|writeControl|writeInterval|executeControl|timeStart|endTime|maxCo' \
     /tmp/prova5/system/controlDict | grep -v '^\s*//'
