#!/bin/bash
C=/root/zefiro/runs/cfd/run3/caso
echo "--- coda del log del solutore ---"
tail -6 "$C/log.reactingFoam"
echo "--- errori ---"
grep -iE "FOAM FATAL|error|signal|Aborted|terminate" "$C/log.reactingFoam" | tail -5
echo "--- log del lancio ---"
ls -la /root/zefiro/runs/cfd/lancio3.log 2>/dev/null && cat /root/zefiro/runs/cfd/lancio3.log
echo "--- uptime della WSL ---"
uptime
echo "--- dmesg recente ---"
dmesg 2>/dev/null | tail -8
