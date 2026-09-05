#!/bin/bash
C=/root/zefiro/runs/cfd/run3/caso
echo "--- deltaT ultimi ---"
grep "^deltaT = " "$C/log.reactingFoam" | tail -5
echo "--- deltaT ogni 200 passi ---"
grep "^deltaT = " "$C/log.reactingFoam" | awk 'NR%200==1'
echo "--- Courant ultimi ---"
grep "^Courant Number" "$C/log.reactingFoam" | tail -3
echo "--- superfici disponibili ---"
ls "$C/postProcessing/superfici" | sort -g | tail -4
echo "--- contenuto ultimo ---"
L=$(ls "$C/postProcessing/superfici" | sort -g | tail -1)
ls "$C/postProcessing/superfici/$L"
