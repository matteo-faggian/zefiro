#!/bin/bash
# Misura pressione e portata sui tre contorni del caso indicato, sull'ultimo
# istante di volume scritto. Risponde con un numero alla domanda "gira con le
# pressioni giuste, da ogni lato?": sugli ingressi non si impone nessuna
# pressione (si impone la portata), quindi quella che si sviluppa e' un
# RISULTATO, e va confrontata con i 7.315 bar del dimensionamento.
set +e
source /mnt/i/AA_ENGINE/_wsl/ambiente.sh
C=${1:-/root/zefiro/runs/cfd/run4/caso}
cd "$C" || exit 2
T="${2:-$(ls -d processor0/[0-9]* 2>/dev/null | xargs -n1 basename | sort -g | tail -1)}"
echo "istante: $T s"
cat > system/misura_contorni <<'DICT'
FoamFile
{ version 2.0; format ascii; class dictionary; location "system";
  object misura_contorni; }
functions
{
DICT
for P in ingresso_aria ingresso_gpl uscita; do
cat >> system/misura_contorni <<DICT
    p_$P
    {
        type surfaceFieldValue; libs (fieldFunctionObjects);
        regionType patch; name $P; operation areaAverage;
        fields (p T); writeFields no; log yes;
    }
    portata_$P
    {
        type surfaceFieldValue; libs (fieldFunctionObjects);
        regionType patch; name $P; operation sum;
        fields (phi); writeFields no; log yes;
    }
DICT
done
echo "}" >> system/misura_contorni
mpirun -np 24 postProcess -parallel -time "$T" -fields "(p T U phi)" \
    -dict system/misura_contorni > log.contorni 2>&1
grep -E "areaAverage|sum\(" log.contorni | grep -vE "^\s*$"
