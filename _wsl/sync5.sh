#!/bin/bash
# Allinea la copia di lavoro WSL al master su I:\AA_ENGINE.
# Il master e' il repo git; /root/zefiro e' solo la copia su disco veloce da
# cui gira OpenFOAM. I risultati (runs/) NON si toccano.
set -e
M=/mnt/i/AA_ENGINE
W=/root/zefiro
for d in src scripts tests config; do
  rsync -a --delete --exclude '__pycache__' "$M/$d/" "$W/$d/"
done
cp "$M/pyproject.toml" "$W/pyproject.toml"
# i .sh e i .py arrivano da Windows: via i CR, altrimenti python e bash si
# lamentano di caratteri invisibili.
find "$W/scripts" "$W/src" "$W/tests" -type f \( -name '*.py' -o -name '*.sh' \) \
     -exec sed -i 's/\r$//' {} +
echo "=== dopo la sincronizzazione ==="
grep -n 'SCONFINAMENTO_GETTO\|def tratto_del_getto' "$W/scripts/cfd/mesh_iniettore.py"
echo "---"
grep -n 'scarto_fianco' "$W/scripts/cfd/caso_mescolamento.py"
echo "---"
ls -la "$W/scripts/cfd/" | grep -E 'mesh_iniettore|caso_mescol'
