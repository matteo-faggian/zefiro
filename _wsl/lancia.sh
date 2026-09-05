#!/bin/bash
# Supervisore: rilancia il solutore se muore per cause esterne (la macchina
# virtuale WSL si e' gia' spenta due volte sotto di lui). Con
# `startFrom latestTime` ogni ripartenza riprende dall'ultimo istante scritto,
# quindi il costo di un'interruzione e' al massimo un intervallo di scrittura.
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
cp /mnt/i/AA_ENGINE/scripts/cfd/*.py /root/zefiro/scripts/cfd/
export PYTHONPATH=/root/zefiro/src
C=/root/zefiro/runs/cfd/run3/caso
if [ ! -d "$C/processor0" ]; then
  cd /root/zefiro
  rm -rf runs/cfd/run3; mkdir -p runs/cfd/run3
  /root/zef/bin/python scripts/cfd/caso_mescolamento.py --caso runs/cfd/run3/caso \
      --fine 7.0e-5 --grossa 3.5e-4 --tempo 1.6e-4 --scrittura 1.0e-5 \
      --campionamento 5.0e-7 --proc 24 --binario 2>&1 | grep -E "^mesh|^getti"
  cd "$C"
  gmshToFoam ../iniettore.msh > log.gmshToFoam 2>&1
  cd /root/zefiro
  /root/zef/bin/python scripts/cfd/patch_bordi.py runs/cfd/run3/caso > /dev/null
  cd "$C"
  checkMesh > log.checkMesh 2>&1
  decomposePar -force > log.decomposePar 2>&1
fi
cd "$C"
for tentativo in 1 2 3 4 5 6 7 8 9 10; do
  echo "=== avvio $tentativo === $(date)"
  mpirun -np 24 reactingFoam -parallel >> log.reactingFoam 2>&1
  rc=$?
  ultimo=$(grep "^Time = " log.reactingFoam | tail -1)
  echo "uscita $rc  $ultimo  $(date)"
  #: se e' arrivato in fondo, `Time` ha raggiunto endTime e non serve altro
  if grep -q "^End$" log.reactingFoam; then echo FINITO; break; fi
  sleep 5
done
