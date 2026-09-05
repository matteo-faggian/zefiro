#!/bin/bash
set +e
source /root/of/etc/bashrc > /root/of_env.log 2>&1
echo "WM_MPLIB=$WM_MPLIB  FOAM_MPI=$FOAM_MPI"
which mpirun mpiexec
mpirun --version 2>&1 | head -3
echo "--- log parallelo ---"
head -40 /root/zefiro/runs/cfd/smoke/log.par
echo "--- libPstream ---"
ls -la /root/of/lib/*/libPstream* /root/of/lib/libPstream* 2>/dev/null
find /root/of -name 'libPstream*' 2>/dev/null | head
