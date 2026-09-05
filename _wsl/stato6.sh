#!/bin/bash
# Stato della corsa 6 PIU' la guardia sui sintomi che hanno ucciso la 5.
C=/root/zefiro/runs/cfd/run6/caso
L=$C/log.reactingFoam
t=$(grep '^Time = ' $L | tail -1 | awk '{print $3}')
echo "t = $(awk -v t=$t 'BEGIN{printf "%.1f", t*1e6}') us   processi: $(pgrep -c -f reactingFoam)"
grep 'Courant Number' $L | tail -1
grep 'deltaT = ' $L | tail -1
echo
echo "--- GUARDIA: i sintomi della corsa 5 ---"
nk=$(grep -c 'bounding k' $L)
nT=$(grep -c 'LimitedCells=[1-9]' $L)
echo "  bounding k          : $nk   $([ $nk -gt 0 ] && echo '<-- ALLARME' )"
echo "  celle limitate in T : $nT   $([ $nT -gt 0 ] && echo '<-- ALLARME' )"
grep 'min/max(T)' $L | tail -1
grep 'cumulative =' $L | tail -1
echo
echo "--- pressioni [bar] ---"
for s in p_ingresso_aria p_ingresso_gpl p_uscita; do
  f=$(find $C/postProcessing/$s -name 'surfaceFieldValue.dat' 2>/dev/null | sort | tail -1)
  [ -n "$f" ] && echo "  $s $(grep -v '^#' "$f" | tail -1 | awk '{printf "%8.4f", $2/1e5}')"
done
echo "--- portate [g/s] ---"
for s in portata_ingresso_aria portata_ingresso_gpl portata_uscita; do
  f=$(find $C/postProcessing/$s -name 'surfaceFieldValue.dat' 2>/dev/null | sort | tail -1)
  [ -n "$f" ] && echo "  $s $(grep -v '^#' "$f" | tail -1 | awk '{printf "%9.4f", $2*1e3}')"
done
