#!/bin/bash
export PYTHONPATH=/root/zefiro/src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cp /mnt/i/AA_ENGINE/scripts/cfd/*.py /root/zefiro/scripts/cfd/ 2>/dev/null
sed -i 's/\r$//' /root/zefiro/scripts/cfd/*.py
cd /root/zefiro
nice -n 19 /root/zef/bin/python scripts/cfd/sezioni.py \
    --caso runs/cfd/run4/caso --storia 2>&1 | tail -45
