#!/bin/bash
C=/root/zefiro/runs/cfd/run2/caso
awk '
/^Time = /       {t=$3}
/min\/max\(T\)/  {gsub(",","",$3); tmin=$3+0; tmax=$4+0;
                  if (tmin<260 && !visto) {print "prima T<260 K a t =", t, " Tmin =", tmin, " Tmax =", tmax; visto=1}
                  if (tmax>400 && !vistoM) {print "prima T>400 K a t =", t, " Tmax =", tmax; vistoM=1}}
/^deltaT = /     {dt=$3+0; if (dt<1e-10 && !vistoD) {print "prima deltaT<1e-10 a t =", t, " dt =", dt; vistoD=1}}
' "$C/log.reactingFoam"
echo "--- estremi di T ogni 6000 righe ---"
grep "min/max(T)" "$C/log.reactingFoam" | awk 'NR%6000==1'
