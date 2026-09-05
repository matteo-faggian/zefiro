#!/bin/bash
set +e
source /root/of/etc/bashrc > /root/of_env.log 2>&1
set -e
pkill -f reactingFoam; sleep 1
export PYTHONPATH=/root/zefiro/src
cd /root/zefiro
rm -rf runs/cfd/run1
mkdir -p runs/cfd/run1
echo "=== mesh fine + caso ==="
/root/zef/bin/python scripts/cfd/caso_mescolamento.py --caso runs/cfd/run1/caso \
    --fine 4.0e-5 --grossa 2.5e-4 --tempo 2.0e-4 --scrittura 2.0e-5 \
    --campionamento 1.0e-6 --proc 24 --binario 2>&1 | tail -6
cd /root/zefiro/runs/cfd/run1/caso
echo "=== gmshToFoam ==="
gmshToFoam ../iniettore.msh > log.gmshToFoam 2>&1 || { tail -20 log.gmshToFoam; exit 1; }
cd /root/zefiro
/root/zef/bin/python scripts/cfd/patch_bordi.py runs/cfd/run1/caso
cd /root/zefiro/runs/cfd/run1/caso
echo "=== checkMesh ==="
checkMesh > log.checkMesh 2>&1 || true
grep -E "cells:|Mesh OK|non-orthogonality|skewness|\*\*\*" log.checkMesh | head -20
