#!/bin/bash
cp /mnt/i/AA_ENGINE/scripts/cfd/swirl.py /root/zefiro/scripts/cfd/
sed -i 's/\r$//' /root/zefiro/scripts/cfd/swirl.py
export PYTHONPATH=/root/zefiro/src:/root/zefiro/scripts/cfd
/root/zef/bin/python /root/zefiro/scripts/cfd/swirl.py
