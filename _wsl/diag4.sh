#!/bin/bash
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
echo "FOAM_LIBBIN=$FOAM_LIBBIN  FOAM_MPI=$FOAM_MPI"
ls -d /root/of/lib/*/ | head -20
echo "--- creo il collegamento mancante ---"
ln -sfn /root/of/lib/mpich-3.3 /root/of/lib/sys-mpich
ls -la /root/of/lib/sys-mpich
cd /root/zefiro/runs/cfd/smoke
mpirun -np 4 reactingFoam -parallel > log.par 2>&1
echo "uscita: $?"
tail -8 log.par
