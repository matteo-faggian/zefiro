#!/bin/bash
export PYTHONPATH=/root/zefiro/src
cd /root/zefiro
/root/zef/bin/python scripts/cfd/sezioni.py --caso runs/cfd/run3/caso --storia 2>&1 | tail -40
