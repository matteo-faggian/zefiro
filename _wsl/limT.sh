#!/bin/bash
L=/root/zefiro/runs/cfd/run3/caso/log.reactingFoam
echo "--- ultime righe del limite di temperatura ---"
grep limitTemperature "$L" | tail -4
echo "--- righe totali / righe con celle tagliate ---"
tot=$(grep -c limitTemperature "$L")
att=$(grep limitTemperature "$L" | grep -vc "LimitedCells=0,")
echo "$att su $tot"
