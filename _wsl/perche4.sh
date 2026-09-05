#!/bin/bash
C=/root/zefiro/runs/cfd/run4/caso
echo "=== inizio del lancio ==="; head -30 /root/zefiro/runs/cfd/lancio4.log
echo "=== log solutore ==="; tail -35 "$C/log.reactingFoam" 2>/dev/null || echo "nessun log"
echo "=== esistono i processor? ==="; ls -d "$C"/processor* 2>/dev/null | head -3
