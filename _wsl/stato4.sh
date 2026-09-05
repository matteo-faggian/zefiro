#!/bin/bash
C=/root/zefiro/runs/cfd/run4/caso
echo "--- lancio ---"; tail -14 /root/zefiro/runs/cfd/lancio4.log
echo "--- avanzamento ---"
echo "passi: $(grep -c '^Time = ' $C/log.reactingFoam 2>/dev/null)"
grep "^Time = " "$C/log.reactingFoam" 2>/dev/null | tail -1
grep "^deltaT = " "$C/log.reactingFoam" 2>/dev/null | tail -1
grep "^Courant Number" "$C/log.reactingFoam" 2>/dev/null | tail -1
grep "ExecutionTime" "$C/log.reactingFoam" 2>/dev/null | tail -1
echo "fotogrammi: $(ls $C/postProcessing/superfici 2>/dev/null | wc -l)"
echo "processi: $(pgrep -c -f reactingFoam)"
