#!/bin/bash
C=${1:-/root/zefiro/runs/cfd/smoke}
echo "--- ultime righe log ---"
tail -30 "$C"/log.reactingFoam 2>/dev/null
echo "--- tempi scritti ---"
ls -d "$C"/[0-9]* 2>/dev/null | sort -g | tail -5
echo "--- proc ---"
ps -eo pid,etime,pcpu,cmd | grep -E 'reactingFoam' | grep -v grep | head
