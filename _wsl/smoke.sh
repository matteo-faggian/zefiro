#!/bin/bash
set +e
source /root/of/etc/bashrc > /root/of_env.log 2>&1
set -e
export PYTHONPATH=/root/zefiro/src
cd /root/zefiro
rm -rf runs/cfd/smoke
echo "=== 1-2. mesh grossa + caso ==="
/root/zef/bin/python scripts/cfd/caso_mescolamento.py --caso runs/cfd/smoke \
    --fine 1.2e-4 --grossa 6e-4 --tempo 2.0e-5 --scrittura 1.0e-5 --proc 4 2>&1 | tail -6
cd /root/zefiro/runs/cfd/smoke
echo "=== 3. gmshToFoam ==="
gmshToFoam ../iniettore.msh > log.gmshToFoam 2>&1 || { tail -20 log.gmshToFoam; exit 1; }
cd /root/zefiro
/root/zef/bin/python scripts/cfd/patch_bordi.py runs/cfd/smoke
cd /root/zefiro/runs/cfd/smoke
echo "=== 4. checkMesh ==="
checkMesh > log.checkMesh 2>&1 || true
tail -18 log.checkMesh
echo "=== 5. reactingFoam (breve) ==="
timeout 900 reactingFoam > log.reactingFoam 2>&1 || true
tail -25 log.reactingFoam
