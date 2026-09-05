#!/bin/bash
C=/root/zefiro/runs/cfd/run3/caso
echo "--- lancio ---"; tail -6 /root/zefiro/runs/cfd/run3/lancio3.log
echo "--- avanzamento ---"
grep -c "^Time = " "$C/log.reactingFoam" 2>/dev/null
tail -3 "$C/log.reactingFoam" 2>/dev/null
grep "^Time = " "$C/log.reactingFoam" 2>/dev/null | tail -1
grep "ExecutionTime" "$C/log.reactingFoam" 2>/dev/null | tail -1
echo "--- fotogrammi ---"
ls "$C/postProcessing/superfici" 2>/dev/null | wc -l
echo "--- memoria / cpu ---"
free -g | head -2
ps -eo pcpu,rss,cmd --sort=-pcpu | head -4
