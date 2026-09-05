#!/bin/bash
# Residui e sanita' numerica: e' il margine che decide se si puo' alzare maxCo.
L=/root/zefiro/runs/cfd/run5/caso/log.reactingFoam
echo "=== residui iniziali dell'ultimo passo ==="
tail -60 "$L" | grep -E 'Solving for (Ux|Uy|Uz|C3H8|h|p|k|omega)' | \
  awk '{for(i=1;i<=NF;i++) if($i=="residual"){print $1" "$2" "$3" "$4" iniz="$(i+2)" iter="$NF; break}}' | tail -14
echo
echo "=== errore di continuita' (deve restare piccolo e non crescere) ==="
grep 'cumulative =' "$L" | tail -3
echo
echo "=== il limitatore di temperatura interviene? ==="
grep 'LimitedCells' "$L" | tail -2
echo
echo "=== min/max T e pressione ==="
grep 'min/max(T)' "$L" | tail -2
echo
echo "=== maxCo attuale ==="
grep 'maxCo' /root/zefiro/runs/cfd/run5/caso/system/controlDict
