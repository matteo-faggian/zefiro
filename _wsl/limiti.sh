#!/bin/bash
C=/root/zefiro/runs/cfd/run3/caso
echo "--- attivita' dei limiti (ultime) ---"
grep -E "limitTemperature|limitVelocity" "$C/log.reactingFoam" | tail -4
echo "--- quante volte hanno tagliato qualcosa ---"
grep -E "limitTemperature|limitVelocity" "$C/log.reactingFoam" | grep -v "Limited 0 (0%)" | wc -l
echo "--- estremi di T ---"
grep "min/max(T)" "$C/log.reactingFoam" | tail -2
grep "min/max(T)" "$C/log.reactingFoam" | awk 'NR%400==1' | tail -6
echo "--- passo temporale ---"
grep "^deltaT" "$C/log.reactingFoam" | tail -1
