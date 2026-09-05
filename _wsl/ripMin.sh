#!/bin/bash
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
C=/root/zefiro/runs/cfd/run4/caso
F=$C/system/fvSolution
sed -i 's/^\( *pMin *\).*$/\1                1.0e5;/; s/^\( *pMax *\).*$/\1                5.0e6;/' "$F"
sed -n '/^PIMPLE/,/^}/p' "$F"
