#!/bin/bash
# Dove va il tempo: iterazioni dei solutori lineari per passo.
L=/root/zefiro/runs/cfd/run5/caso/log.reactingFoam
echo "=== ultimo passo, per esteso ==="
awk '/^Time = /{n++} n>=1{buf=buf"\n"$0} /^ExecutionTime/{if(n>=1){print buf; buf=""; n=0}}' \
    "$L" | tail -45
echo
echo "=== ritmo istantaneo sugli ultimi passi ==="
grep -E '^Time = |^ExecutionTime' "$L" | tail -20 | \
awk '/^Time/{t=$3} /^ExecutionTime/{e=$3; if(tp>0) printf "  t=%10.4e us  dt_calcolo=%6.2f s\n", t*1e6, e-ep; tp=t; ep=e}'
