#!/bin/bash
# Costruisce la geometria e scrive lo STL sulla macchina dell'utente: qui la
# memoria c'e' (23 GB contro 8), e il passo di confronto piu' fine ci sta.
set +e
mkdir -p /root/zefiro/src /root/zefiro/scripts /root/zefiro/config
cp -r /mnt/i/AA_ENGINE/src/. /root/zefiro/src/
cp -r /mnt/i/AA_ENGINE/scripts/. /root/zefiro/scripts/
cp -r /mnt/i/AA_ENGINE/config/. /root/zefiro/config/ 2>/dev/null
cd /root/zefiro
nice -n 10 /root/zef/bin/python scripts/progetta.py --spinta 50 --durata 5 \
    --t-bombola 15 --passo 2.5e-4 --stl /root/zefiro/runs/zefiro_50N.stl \
    > /root/zefiro/runs/progetta.log 2>&1
echo "uscita: $?"
tail -3 /root/zefiro/runs/progetta.log
