#!/bin/bash
C=/root/zefiro/runs/cfd/run3/caso
ls "$C/processor0/5.62272055e-05"
echo "--- checkMesh: qualita' ---"
grep -iE "volume|skew|non-ortho|Mesh OK|\*\*\*" "$C/log.checkMesh" | head -30
