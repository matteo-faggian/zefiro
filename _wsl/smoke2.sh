#!/bin/bash
set +e
source /root/of/etc/bashrc > /root/of_env.log 2>&1
set -e
cp /mnt/i/AA_ENGINE/scripts/cfd/caso_mescolamento.py /root/zefiro/scripts/cfd/
md5sum /root/zefiro/scripts/cfd/caso_mescolamento.py
export PYTHONPATH=/root/zefiro/src
cd /root/zefiro
rm -rf runs/cfd/smoke
/root/zef/bin/python scripts/cfd/caso_mescolamento.py --caso runs/cfd/smoke \
    --fine 1.2e-4 --grossa 6e-4 --tempo 1.5e-5 --scrittura 5.0e-6 \
    --campionamento 1.0e-6 --proc 4 2>&1 | tail -4
cd /root/zefiro/runs/cfd/smoke
gmshToFoam ../iniettore.msh > log.gmshToFoam 2>&1 || { tail -20 log.gmshToFoam; exit 1; }
cd /root/zefiro
/root/zef/bin/python scripts/cfd/patch_bordi.py runs/cfd/smoke > /dev/null
cd /root/zefiro/runs/cfd/smoke
echo "=== reactingFoam con i functionObjects ==="
timeout 900 reactingFoam > log.reactingFoam 2>&1
echo "uscita: $?"
tail -12 log.reactingFoam
echo "=== superfici prodotte ==="
find postProcessing -name '*.vtk' -o -name '*.vtp' | head -20
find postProcessing -type d | head -20
