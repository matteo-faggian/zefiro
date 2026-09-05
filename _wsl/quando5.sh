#!/bin/bash
# QUANDO e' cominciata la degenerazione? Si cerca il PRIMO segno di ciascun
# sintomo e si legge il tempo simulato a quel punto. Senza questo, attribuire
# la causa e' un'opinione.
L=/root/zefiro/runs/cfd/run5/caso/log.reactingFoam
echo "=== primo 'bounding k' ==="
awk '/^Time = /{t=$3} /bounding k/{print "  riga "NR"  t = "t*1e6" us  "$0; exit}' "$L"
echo "=== primo 'bounding omega' ==="
awk '/^Time = /{t=$3} /bounding omega/{print "  riga "NR"  t = "t*1e6" us  "$0; exit}' "$L"
echo "=== prima cella limitata in temperatura ==="
awk '/^Time = /{t=$3} /LimitedCells=[1-9]/{print "  riga "NR"  t = "t*1e6" us  "$0; exit}' "$L"
echo "=== primo intervento di pressureControl ==="
awk '/^Time = /{t=$3} /pressureControl/{print "  riga "NR"  t = "t*1e6" us  "$0; exit}' "$L"
echo
echo "=== quanto vale maxCo e da quando (la modifica e' delle 17:14 circa) ==="
grep 'maxCo' /root/zefiro/runs/cfd/run5/caso/system/controlDict
echo
echo "=== istanti ancora su disco (purgeWrite 4) ==="
ls -d /root/zefiro/runs/cfd/run5/caso/processor0/*/ | grep -v constant | sort -g | tr '\n' ' '
echo
echo "=== la corsa 4, per confronto: aveva mai avuto questi sintomi? ==="
grep -c 'bounding k' /root/zefiro/runs/cfd/run4/caso/log.reactingFoam
grep -c 'LimitedCells=[1-9]' /root/zefiro/runs/cfd/run4/caso/log.reactingFoam
