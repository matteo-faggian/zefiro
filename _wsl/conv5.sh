#!/bin/bash
cp /mnt/i/AA_ENGINE/scripts/cfd/convergenza.py /root/zefiro/scripts/cfd/
sed -i 's/\r$//' /root/zefiro/scripts/cfd/convergenza.py
/root/zef/bin/python /root/zefiro/scripts/cfd/convergenza.py \
    --caso "${1:-/root/zefiro/runs/cfd/run5/caso}"
