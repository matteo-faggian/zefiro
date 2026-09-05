#!/bin/bash
# CORSA 5 - preparazione. Mesh RIFATTA da zero, non ereditata.
#
# Perche' da zero: le corse 3 e 4 giravano sulla stessa mesh, che conteneva
# una tasca cieca di 0.731 x 2.000 mm nel metallo davanti a ogni getto - un
# residuo di un `+2 mm` di sfondamento pensato per un taglio booleano e finito
# dentro una fusione. E' stata misurata (aria ferma, Y_GPL = 0 su 927 punti) e
# non aveva contaminato il verdetto, ma e' geometria che il motore non ha.
#
# Raffinamento IDENTICO a quello delle corse 3 e 4 (fine 7.0e-5, grossa
# 3.5e-4): fra la 4 e la 5 deve cambiare UNA cosa sola, altrimenti la
# differenza fra i due risultati non e' attribuibile a niente.
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
bash /mnt/i/AA_ENGINE/_wsl/run.sh sync5 > /dev/null 2>&1
export PYTHONPATH=/root/zefiro/src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd /root/zefiro
rm -rf runs/cfd/run5; mkdir -p runs/cfd/run5
echo "=== mesh + caso === $(date)"
/root/zef/bin/python scripts/cfd/caso_mescolamento.py \
    --caso runs/cfd/run5/caso --fine 7.0e-5 --grossa 3.5e-4 \
    --tempo 1.0e-3 --scrittura 2.0e-5 --campionamento 2.0e-6 \
    --proc 24 --binario --passo-sonde 20 --media-da 5.0e-4 2>&1 \
  | grep -E "^mesh|^getti|^aria|^portate|^velocita|^J con|^lunghezza|^caso"
C=/root/zefiro/runs/cfd/run5/caso
[ -s runs/cfd/run5/iniettore.msh ] || { echo "MESH ASSENTE"; exit 2; }
ls -la runs/cfd/run5/iniettore.msh
cd "$C"
echo "=== gmshToFoam === $(date)"
gmshToFoam ../iniettore.msh > log.gmshToFoam 2>&1
cd /root/zefiro
/root/zef/bin/python scripts/cfd/patch_bordi.py runs/cfd/run5/caso > /dev/null
cd "$C"
echo "=== checkMesh === $(date)"
checkMesh > log.checkMesh 2>&1
grep -E "Mesh OK|\*\*\*|cells:|max.*non-orth|Max skewness" log.checkMesh
echo "=== decomposePar === $(date)"
decomposePar -force > log.decomposePar 2>&1
tail -3 log.decomposePar
echo "=== PRONTO === $(date)"
