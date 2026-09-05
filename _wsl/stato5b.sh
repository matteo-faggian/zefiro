#!/bin/bash
C=/root/zefiro/runs/cfd/run5/caso
L=$C/log.reactingFoam
t=$(grep '^Time = ' $L | tail -1 | awk '{print $3}')
e=$(grep 'ExecutionTime' $L | tail -1 | awk '{print $3}')
echo "tempo simulato : $t s"
echo "tempo di calcolo: $e s"
awk -v t="$t" -v e="$e" 'BEGIN{
  if (t>0) {
    r = e/(t*1e6);
    printf "costo          : %.1f s per microsecondo simulato\n", r;
    printf "per arrivare a 1000 us: %.1f ore in tutto, ne mancano %.1f\n",
           r*1000/3600, r*(1000-t*1e6)/3600;
  }}'
echo
grep 'Courant Number' $L | tail -1
grep 'deltaT = ' $L | tail -1
echo "processi vivi: $(pgrep -c -f reactingFoam)"
echo
echo "--- ultime pressioni misurate [bar] ---"
for s in p_ingresso_aria p_ingresso_gpl p_uscita; do
  f=$(find $C/postProcessing/$s -name 'surfaceFieldValue.dat' | sort | tail -1)
  [ -n "$f" ] && echo "$s $(grep -v '^#' "$f" | tail -1 | awk '{printf "%8.4f bar  %6.1f K", $2/1e5, $3}')"
done
echo "--- portate [g/s] ---"
for s in portata_ingresso_aria portata_ingresso_gpl portata_uscita; do
  f=$(find $C/postProcessing/$s -name 'surfaceFieldValue.dat' | sort | tail -1)
  [ -n "$f" ] && echo "$s $(grep -v '^#' "$f" | tail -1 | awk '{printf "%9.4f", $2*1e3}')"
done
