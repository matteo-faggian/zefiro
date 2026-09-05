#!/bin/bash
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
pkill -f reactingFoam; pkill -f _lancia.sh; sleep 2
cd /root/zefiro/runs/cfd/smoke
foamDictionary -entry endTime -set 1.2e-06 system/controlDict > /dev/null 2>&1
foamDictionary -entry functions -remove system/controlDict > /dev/null 2>&1
misura () {   # $1 = etichetta, resto = comando
  local eti="$1"; shift
  rm -rf 0.* [1-9]* processor*/[1-9]* processor*/0.* 2>/dev/null
  local t0=$(date +%s.%N)
  "$@" > log.bench 2>&1
  local t1=$(date +%s.%N)
  local n=$(grep -c "^Time = " log.bench)
  echo "$eti: $n passi in $(echo "$t1 - $t0" | bc) s  -> $(echo "scale=3; ($t1 - $t0)/$n" | bc) s/passo"
}
misura "seriale        " reactingFoam
for np in 4 8 16 24; do
  decomposePar -force > /dev/null 2>&1
  sed -i "s/numberOfSubdomains.*/numberOfSubdomains  $np;/" system/decomposeParDict
  decomposePar -force > /dev/null 2>&1
  misura "mpi $np rank    " mpirun -np $np reactingFoam -parallel
done
export FI_PROVIDER=shm
misura "mpi 24 shm      " mpirun -np 24 reactingFoam -parallel
