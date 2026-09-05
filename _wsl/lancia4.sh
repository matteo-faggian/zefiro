#!/bin/bash
# Corsa 4. Differenze rispetto alla 3, tutte nel caso, non qui:
#   - ingressi a PORTATA imposta (rompe l'anello p->rho->mdot->p che ha
#     fermato le corse 1, 2 e 3 sempre allo stesso labbro di foro);
#   - pMin/pMax nel PIMPLE al posto del limitatore di velocita';
#   - due correttori non-ortogonali (la mesh e' tetraedrica, 57.9 gradi).
#
# endTime = 2.5e-4 s e NON 1.3e-4. Il dominio ha due tempi caratteristici che
# differiscono di 20 volte:
#     anello  (-6.6 -> 0 mm) a 120 m/s  ->    55 us per transito
#     camera  ( 0 -> 12 mm)  a  10.6 m/s -> 1132 us per transito
# 2.5e-4 s sono ~4.5 transiti dell'anello: bastano per i piani di misura del
# mescolamento, che stanno tutti fra -2.5 e +1.7 mm. NON bastano per la camera,
# e non e' un problema di pazienza: quattro transiti di camera sarebbero 4.5 ms,
# cioe' giorni di calcolo. Quello che la camera fa a regime va chiesto a un
# calcolo STAZIONARIO, non a questo.
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
cp /mnt/i/AA_ENGINE/scripts/cfd/*.py /root/zefiro/scripts/cfd/
sed -i 's/\r$//' /root/zefiro/scripts/cfd/*.py
export PYTHONPATH=/root/zefiro/src
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
C=/root/zefiro/runs/cfd/run4/caso
if [ ! -d "$C/processor0" ]; then
  cd /root/zefiro
  rm -rf runs/cfd/run4; mkdir -p runs/cfd/run4
  #: --solo-caso: la mesh di run3 e' gia' stata verificata (checkMesh OK,
  #: 1.15 M celle, non-ortogonalita' 57.9). Rifarla non aggiungerebbe niente e
  #: introdurrebbe una variabile in piu' fra le due corse.
  #: la mesh e' quella di run3, per nome esplicito: vedi prep4.sh.
  M=/root/zefiro/runs/cfd/run3/iniettore.msh
  [ -s "$M" ] || { echo "manca la mesh di run3"; exit 2; }
  cp "$M" runs/cfd/run4/iniettore.msh
  /root/zef/bin/python scripts/cfd/caso_mescolamento.py --caso runs/cfd/run4/caso \
      --solo-caso --fine 7.0e-5 --grossa 3.5e-4 --tempo 2.5e-4 \
      --scrittura 1.0e-5 --campionamento 5.0e-7 --proc 24 --binario 2>&1 \
    | grep -E "^mesh|^getti|^aria|^portate|^velocita|^J con|^lunghezza"
  cd "$C"
  gmshToFoam ../iniettore.msh > log.gmshToFoam 2>&1
  cd /root/zefiro
  /root/zef/bin/python scripts/cfd/patch_bordi.py runs/cfd/run4/caso > /dev/null
  cd "$C"
  checkMesh > log.checkMesh 2>&1
  grep -E "Mesh OK|\*\*\*" log.checkMesh
  decomposePar -force > log.decomposePar 2>&1
fi
cd "$C"
for t in 1 2 3 4 5 6 7 8 9 10; do
  echo "=== avvio $t === $(date)"
  mpirun -np 24 reactingFoam -parallel >> log.reactingFoam 2>&1
  echo "uscita $?  $(grep '^Time = ' log.reactingFoam | tail -1)  $(date)"
  if grep -q "^End$" log.reactingFoam; then echo FINITO; break; fi
  sleep 5
done
