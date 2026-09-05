#!/bin/bash
# Cosa c'e' gia' in casa, prima di decidere come pubblicare.
M=/mnt/i/AA_ENGINE
echo "=== git e gh disponibili? ==="
git --version 2>/dev/null || echo "  git ASSENTE in WSL"
gh --version 2>/dev/null | head -1 || echo "  gh (GitHub CLI) ASSENTE in WSL"
echo
echo "=== identita' git gia' configurata? ==="
git config --global user.name 2>/dev/null || echo "  user.name non impostato"
git config --global user.email 2>/dev/null || echo "  user.email non impostato"
echo
echo "=== gh gia' autenticato? ==="
gh auth status 2>&1 | head -5
echo
echo "=== CACCIA AI SEGRETI (prima di qualunque push) ==="
grep -rIl --exclude-dir=.git --exclude-dir=__pycache__ --exclude-dir='(eliminare)' \
     -E '(api[_-]?key|secret|password|passwd|token|BEGIN [A-Z ]*PRIVATE KEY|ghp_|sk-)' \
     "$M" 2>/dev/null | head -20
echo "  (nessuna riga sopra = nessun sospetto)"
echo
echo "=== quanto peserebbe il repo con il .gitignore attuale ==="
cd "$M" 2>/dev/null || exit 1
echo "  file che NON sono ignorati, per estensione:"
find . -type f -not -path './.git/*' -not -path './runs/*' \
       -not -path './(eliminare)/*' -not -path '*__pycache__*' \
       -not -name '*.stl' -not -name '*.msh' -not -name '*.step' \
     | sed 's/.*\.//' | sort | uniq -c | sort -rn | head -12
echo "  peso totale di quei file:"
find . -type f -not -path './.git/*' -not -path './runs/*' \
       -not -path './(eliminare)/*' -not -path '*__pycache__*' \
       -not -name '*.stl' -not -name '*.msh' -not -name '*.step' \
     -printf '%s\n' | awk '{s+=$1} END {printf "  %.2f MB su %d file\n", s/1048576, NR}'
