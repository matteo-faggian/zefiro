#!/bin/bash
# Ritmo ISTANTANEO, non la media dall'inizio: dopo un cambio di maxCo la media
# cumulata resta zavorrata dalla parte lenta e non dice piu' niente di utile.
L=/root/zefiro/runs/cfd/run5/caso/log.reactingFoam
grep -E '^Time = |^ExecutionTime' "$L" | tail -60 | \
awk '/^Time/{t=$3} /^ExecutionTime/{e=$3;
       if(tp>0){dt+=t-tp; de+=e-ep; n++} tp=t; ep=e}
     END{ if(n>0){ r=de/(dt*1e6);
       printf "su %d passi: %.1f s per microsecondo simulato\n", n, r;
       printf "tempo simulato ora: %.1f us\n", t*1e6;
       printf "mancano %.1f ore per arrivare a 1000 us\n", r*(1000-t*1e6)/3600;
     }}'
