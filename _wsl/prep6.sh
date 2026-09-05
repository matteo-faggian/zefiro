#!/bin/bash
# CORSA 6. Cambia UNA cosa rispetto alla 5, in due pezzi che stanno insieme:
#   - gli ingressi salgono in rampa in 100 us invece che a gradino;
#   - lo sbocco ancora la pressione media a p_c (fixedMean) invece di
#     lasciarla galleggiare (waveTransmissive con lInf 2.7 volte il dominio).
# La mesh e' LA STESSA della corsa 5 (verificata: 6 bordi, ingresso_gpl 219
# facce, aree entro lo 0.6 % dell'analitico). Rifarla introdurrebbe una
# variabile in piu' fra due corse che devono essere confrontabili.
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
bash /mnt/i/AA_ENGINE/_wsl/run.sh sync5 > /dev/null 2>&1
export PYTHONPATH=/root/zefiro/src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
cd /root/zefiro
M=/root/zefiro/runs/cfd/run5/iniettore.msh
[ -s "$M" ] || { echo "manca la mesh della corsa 5"; exit 2; }
rm -rf runs/cfd/run6; mkdir -p runs/cfd/run6
cp "$M" runs/cfd/run6/iniettore.msh
echo "=== caso === $(date)"
/root/zef/bin/python scripts/cfd/caso_mescolamento.py \
    --caso runs/cfd/run6/caso --solo-caso \
    --tempo 1.0e-3 --scrittura 2.0e-5 --campionamento 2.0e-6 \
    --proc 24 --binario --passo-sonde 20 --media-da 6.0e-4 2>&1 \
  | grep -E "^getti|^aria|^portate|^velocita|^J con|^caso"
C=/root/zefiro/runs/cfd/run6/caso
echo "=== rampa e sbocco, come sono finiti nel caso ==="
sed -n '/ingresso_aria/,/}/p' $C/0/U | grep -A5 massFlowRate
sed -n '/uscita/,/}/p' $C/0/p | head -6
cd "$C"
echo "=== gmshToFoam === $(date)"
gmshToFoam ../iniettore.msh > log.gmshToFoam 2>&1
cd /root/zefiro
/root/zef/bin/python scripts/cfd/patch_bordi.py runs/cfd/run6/caso > /dev/null
cd "$C"
grep -c . constant/polyMesh/boundary > /dev/null
echo "bordi: $(sed -n '19p' constant/polyMesh/boundary)"
checkMesh > log.checkMesh 2>&1
grep -E "Mesh OK|\*\*\*|cells:" log.checkMesh
decomposePar -force > log.decomposePar 2>&1
tail -2 log.decomposePar
echo "=== PRONTO === $(date)"
