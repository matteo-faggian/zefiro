#!/bin/bash
# Primi controlli sulla corsa 5, quelli che si fanno PRIMA di lasciarla andare
# per ore: le aree delle patch sono quelle analitiche? le sonde scrivono?
C=/root/zefiro/runs/cfd/run5/caso
echo "=== avanzamento ==="
grep '^Time = ' $C/log.reactingFoam | tail -1
grep 'ExecutionTime' $C/log.reactingFoam | tail -1
echo
echo "=== aree delle patch misurate dal solutore ==="
for s in p_ingresso_aria p_ingresso_gpl p_uscita; do
  f=$(find $C/postProcessing/$s -name 'surfaceFieldValue.dat' 2>/dev/null | head -1)
  [ -n "$f" ] && echo "$s: $(grep -m1 Area "$f")"
done
echo
echo "=== confronto con l'area analitica ==="
/root/zef/bin/python - <<'PY'
import math, sys
sys.path.insert(0, "/root/zefiro/scripts/cfd")
from mesh_iniettore import quote_dal_progetto
d = quote_dal_progetto(); n = d["n_getti"]
print(f"  ingresso_aria  {math.pi*(d['R_anello']**2-d['R_getti']**2)/n:.6e} m2")
print(f"  ingresso_gpl   {math.pi/4*d['d_getto']**2:.6e} m2")
print(f"  uscita         {math.pi*(d['R_c']**2-d['r_centerbody']**2)/n:.6e} m2")
PY
echo
echo "=== le sonde scrivono davvero? (il difetto della corsa 4) ==="
for s in massa_nel_dominio p_uscita portata_uscita; do
  f=$(find $C/postProcessing/$s -name '*.dat' 2>/dev/null | head -1)
  [ -n "$f" ] && echo "$s: $(grep -vc '^#' "$f") campioni"
done
echo
echo "=== y+ (validita' delle funzioni di parete) ==="
grep -A2 'yPlus' $C/log.reactingFoam | grep -E 'min|max|average' | tail -3
