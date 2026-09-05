#!/bin/bash
set +e
source /root/of/etc/bashrc > /root/of_env.log 2>&1
cd /root/zefiro/runs/cfd/smoke
rm -rf processor*
foamDictionary -entry endTime -set 3e-07 system/controlDict > /dev/null 2>&1
decomposePar -force > log.decomposePar 2>&1
tail -3 log.decomposePar
mpirun -np 4 reactingFoam -parallel > log.par 2>&1
echo "uscita mpirun: $?"
tail -6 log.par
