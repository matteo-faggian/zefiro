#!/bin/bash
C=/root/zefiro/runs/cfd/run4/caso
echo "=== boundary del caso ricomposto ==="
sed -n '/^[0-9]/,$p' "$C/constant/polyMesh/boundary" | head -60
echo "=== boundary del processor0 ==="
grep -A4 -E "^\s+(uscita|ingresso_aria|pareti)$" "$C/processor0/constant/polyMesh/boundary" | head -30
echo "=== patch_bordi.py ==="
cat /mnt/i/AA_ENGINE/scripts/cfd/patch_bordi.py
