#!/bin/bash
# Impostazioni numeriche effettive del run 4 e costo per microsecondo simulato.
C=/root/zefiro/runs/cfd/run4/caso
echo "=== controlDict (senza functions) ==="
sed -n '1,60p' $C/system/controlDict | grep -v '^//' | grep -v '^ *$'
echo
echo "=== PIMPLE ==="
sed -n '/PIMPLE/,/^}/p' $C/system/fvSolution
echo
echo "=== ddt e div ==="
sed -n '/ddtSchemes/,/^}/p' $C/system/fvSchemes
echo
echo "=== costo ==="
grep -c '^Time = ' $C/log.reactingFoam
grep 'ExecutionTime' $C/log.reactingFoam | tail -1
grep '^Time = ' $C/log.reactingFoam | tail -1
echo "--- Courant ---"
grep 'Courant Number' $C/log.reactingFoam | tail -3
echo "--- deltaT ---"
grep 'deltaT = ' $C/log.reactingFoam | tail -2
echo "--- celle ---"
grep -E 'cells:|faces:' $C/log.checkMesh | head -4
