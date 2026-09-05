#!/bin/bash
# CORSA 5 - esecuzione. 24 processi, endTime 1.0e-3 s.
#
# PERCHE' 1000 MICROSECONDI E NON 250. La corsa 4 si fermava a 250 us con le
# pressioni ancora in movimento (7.32 -> 6.75 bar sull'alimentazione aria) e
# la portata uscente diversa dalla entrante: il dominio stava ancora
# scaricando la massa presa in piu' durante il riempimento. Il tempo con cui
# quel transitorio si spegne NON e' il tempo di transito della camera
# (1264 us): e' il tempo di rilassamento della condizione allo sbocco,
#     lInf / c = 0.05 m / 340 m/s = 147 us,
# che e' anche quello che si legge nel decadimento misurato della corsa 4
# (tau ~ 140 us fra 50 e 250 us). Mille microsecondi sono ~7 tempi di
# rilassamento: l'eccesso di pressione cala di e^-7, cioe' sotto il per mille.
#
# NON e' una previsione su cui fidarsi: `runTimeModifiable yes` e le sonde a
# 20 passi permettono di GUARDARE la curva mentre gira e di alzare endTime se
# non e' piatta. Il criterio si dichiara prima: massa contenuta nel dominio
# stazionaria entro l'1 % su 200 us, e sbilancio di portata sotto il 2 %.
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
C=/root/zefiro/runs/cfd/run5/caso
[ -d "$C/processor0" ] || { echo "caso non decomposto: lancia prima prep5"; exit 2; }
cd "$C"
# il ciclo di riavvio non e' pigrizia: `startFrom latestTime` riprende
# dall'ultima scrittura, e nelle corse precedenti la macchina virtuale si e'
# spenta sotto il solutore piu' di una volta. Dieci tentativi, e si ferma da
# solo quando trova "End".
for t in $(seq 1 12); do
  echo "=== avvio $t === $(date)"
  mpirun -np 24 reactingFoam -parallel >> log.reactingFoam 2>&1
  echo "uscita $?  $(grep '^Time = ' log.reactingFoam | tail -1)  $(date)"
  if grep -q "^End$" log.reactingFoam; then echo FINITO; break; fi
  sleep 5
done
