#!/bin/bash
# Perche' la corsa 5 si e' fermata: crash numerico o macchina caduta sotto?
C=/root/zefiro/runs/cfd/run5/caso
L=$C/log.reactingFoam
echo "=== il lanciatore che cosa dice ==="
cat /root/lancia5.out 2>/dev/null | tail -30
echo
echo "=== ultime righe del log del solutore ==="
tail -25 "$L"
echo
echo "=== c'e' un 'End'? ==="
grep -c '^End$' "$L"
echo
echo "=== segni di divergenza ==="
grep -nE 'FOAM FATAL|Floating point|diverg|bounding|Maximum number of iterations' "$L" | tail -15
echo
echo "=== il limitatore di pressione e' intervenuto? ==="
grep -c 'pMin\|pMax' "$L"
echo
echo "=== uptime della macchina virtuale ==="
uptime
