#!/bin/bash
ls /root/of/etc/conda/activate.d/ 2>/dev/null && cat /root/of/etc/conda/activate.d/* 2>/dev/null | head -60
echo "=== prova con mpich ==="
set +e
source /root/of/etc/bashrc > /dev/null 2>&1
export FOAM_MPI=mpich-3.3
export LD_LIBRARY_PATH=/root/of/lib/mpich-3.3:$LD_LIBRARY_PATH
cd /root/zefiro/runs/cfd/smoke
mpirun -np 4 reactingFoam -parallel > log.par 2>&1
echo "uscita: $?"
tail -8 log.par
