#!/bin/bash
# I tre sospetti, riga per riga: falsi positivi o roba vera?
M=/mnt/i/AA_ENGINE
for f in scripts/cfd/caso_mescolamento.py scripts/cfd/LANCIA_SU_WSL.md; do
  echo "=== $f ==="
  grep -nIE '(api[_-]?key|secret|password|passwd|token|BEGIN [A-Z ]*PRIVATE KEY|ghp_|sk-)' \
       "$M/$f" | head -8
  echo
done
echo "=== e nei file di configurazione? ==="
ls $M/config/
grep -rIn --include='*.yaml' --include='*.yml' --include='*.json' --include='*.toml' \
     -E '(key|secret|token|password)' $M/config $M/pyproject.toml $M/environment.yml 2>/dev/null | head -10
echo "  (vuoto = pulito)"
