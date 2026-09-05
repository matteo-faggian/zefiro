#!/bin/bash
# y+ : le funzioni di parete sono lecite su questa mesh?
# k-omega SST con wall function vuole y+ fuori dalla zona cuscinetto
# (5 < y+ < 30 e' la terra di nessuno), oppure y+ < 1 se si vuole risolvere
# il sottostrato viscoso. Se cade in mezzo, l'attrito a parete non e' ne'
# risolto ne' modellato - e va detto, non nascosto.
L=/root/zefiro/runs/cfd/run5/caso/log.reactingFoam
echo "=== y+ scritto dal solutore ==="
grep -i 'yPlus' "$L" | tail -8
echo
echo "=== istanti di volume scritti ==="
ls -d /root/zefiro/runs/cfd/run5/caso/processor0/*/ 2>/dev/null | tail -6
echo
echo "=== campi presenti nell'ultimo istante ==="
ultimo=$(ls -d /root/zefiro/runs/cfd/run5/caso/processor0/*/ | grep -v constant | sort -g | tail -1)
echo "$ultimo"
ls "$ultimo" | tr '\n' ' '
echo
