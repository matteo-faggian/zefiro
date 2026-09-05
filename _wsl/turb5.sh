#!/bin/bash
cp /mnt/i/AA_ENGINE/scripts/cfd/turbolenza.py /root/zefiro/scripts/cfd/
cp /mnt/i/AA_ENGINE/scripts/cfd/sezione.py /root/zefiro/scripts/cfd/
sed -i 's/\r$//' /root/zefiro/scripts/cfd/turbolenza.py /root/zefiro/scripts/cfd/sezione.py
export PYTHONPATH=/root/zefiro/src:/root/zefiro/scripts/cfd
/root/zef/bin/python /root/zefiro/scripts/cfd/turbolenza.py \
    --caso "${1:-/root/zefiro/runs/cfd/run5/caso}"
