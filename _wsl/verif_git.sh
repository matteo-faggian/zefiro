#!/bin/bash
# La storia ricostruita e' credibile? Tre prove.
R=/root/zefiro_git
cd $R || exit 1
echo "=== 1. i commit sono INCREMENTALI o riscrivono tutto? ==="
echo "   (se ogni passo tocca ~tutti i file, la normalizzazione della radice"
echo "    ha fallito e la storia non vale niente)"
prec=""
for c in $(git log --reverse --format=%h); do
  if [ -n "$prec" ]; then
    s=$(git diff --shortstat $prec $c)
    printf "  %s -> %s  %s\n" "$prec" "$c" "${s:-nessuna differenza}"
  fi
  prec=$c
done
echo
echo "=== 2. file piu' grossi finiti nel repo ==="
git rev-list --objects --all \
 | git cat-file --batch-check='%(objecttype) %(objectname) %(objectsize) %(rest)' \
 | awk '$1=="blob"{print $3, $4}' | sort -rn | head -5 \
 | awk '{printf "  %8.2f KB  %s\n", $1/1024, $2}'
echo
echo "=== 3. integrita' e coerenza con la cartella di lavoro ==="
git fsck --no-progress 2>&1 | head -3
echo "  file nel commit finale: $(git ls-tree -r HEAD --name-only | wc -l)"
echo "  differenze fra HEAD e I:\\AA_ENGINE (escluse le esclusioni):"
diff -rq --exclude='.git' --exclude='__pycache__' --exclude='runs' \
     --exclude='(eliminare)' --exclude='.pytest_cache' --exclude='.ruff_cache' \
     $R /mnt/i/AA_ENGINE 2>/dev/null \
  | grep -vE '\.(stl|msh|step|parquet|db)$' | head -10
echo "  (nessuna riga sopra = il commit finale E' lo stato attuale)"
