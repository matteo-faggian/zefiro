#!/bin/bash
# Rilancia il supervisore SOLO se non sta girando niente. Idempotente: se il
# solutore e' vivo non tocca nulla.
if pgrep -f reactingFoam > /dev/null; then
  echo "vivo: $(grep '^Time = ' /root/zefiro/runs/cfd/run3/caso/log.reactingFoam 2>/dev/null | tail -1)"
  exit 0
fi
if pgrep -f _lancia.sh > /dev/null; then
  echo "supervisore vivo, solutore in avvio"
  exit 0
fi
if grep -q "^End$" /root/zefiro/runs/cfd/run3/caso/log.reactingFoam 2>/dev/null; then
  echo "corsa finita"
  exit 0
fi
echo "RILANCIO $(date)"
tr -d '\r' < /mnt/i/AA_ENGINE/_wsl/lancia.sh > /root/_lancia.sh
chmod +x /root/_lancia.sh
setsid nohup /root/_lancia.sh >> /root/zefiro/runs/cfd/lancio3.log 2>&1 < /dev/null &
