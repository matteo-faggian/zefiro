#!/bin/bash
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
echo "FOAM_MPI=$FOAM_MPI  WM_MPLIB=$WM_MPLIB"
echo "LD_LIBRARY_PATH=$LD_LIBRARY_PATH" | tr ':' '\n' | head -8
cd /root/zefiro/runs/cfd/smoke
mpirun -np 4 reactingFoam -parallel > log.par 2>&1
echo "uscita: $?"
tail -8 log.par
