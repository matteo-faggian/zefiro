#!/bin/bash
C=/root/zefiro/runs/cfd/run3/caso
echo "=== 0/U ==="; sed -n '/boundaryField/,$p' "$C/0/U" | head -40
echo "=== 0/T ingressi ==="; sed -n '/boundaryField/,$p' "$C/0/T" | head -22
echo "=== quote (tutte) ==="; cat "$C/quote.txt"
