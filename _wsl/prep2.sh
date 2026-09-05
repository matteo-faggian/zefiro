#!/bin/bash
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
set -e
export PYTHONPATH=/root/zefiro/src
cd /root/zefiro
rm -rf runs/cfd/run2; mkdir -p runs/cfd/run2
/root/zef/bin/python scripts/cfd/caso_mescolamento.py --caso runs/cfd/run2/caso \
    --fine 7.0e-5 --grossa 3.5e-4 --tempo 1.2e-4 --scrittura 2.0e-5 \
    --campionamento 6.0e-7 --proc 24 --binario 2>&1 | grep -E "^mesh|^getti|^aria|^lunghezza"
cd /root/zefiro/runs/cfd/run2/caso
gmshToFoam ../iniettore.msh > log.gmshToFoam 2>&1
cd /root/zefiro
/root/zef/bin/python scripts/cfd/patch_bordi.py runs/cfd/run2/caso > /dev/null
cd /root/zefiro/runs/cfd/run2/caso
checkMesh > log.checkMesh 2>&1
grep -E "cells:|Mesh OK|non-orthogonality Max|Max skewness|\*\*\*" log.checkMesh
