#!/bin/bash
# I fine-riga. Il master sta su un disco Windows, il calcolo gira in WSL, e gli
# zip venivano da chissa' dove: se meta' repo e' CRLF e meta' LF, ogni futuro
# commit fatto dall'altro lato risultera' una riscrittura totale del file e i
# diff diventeranno illeggibili. Meglio saperlo adesso.
R=/root/zefiro_git
echo "=== quanti file hanno CRLF nel commit finale? ==="
crlf=0; lf=0
while IFS= read -r f; do
  if grep -qU $'\r' "$R/$f" 2>/dev/null; then crlf=$((crlf+1)); else lf=$((lf+1)); fi
done < <(cd $R && git ls-tree -r HEAD --name-only | grep -E '\.(py|sh|md|yaml|yml|toml|txt|js|html)$')
echo "  CRLF: $crlf    LF: $lf"
echo
echo "=== e negli zip di agosto? (primo commit) ==="
cd $R && git stash list >/dev/null 2>&1
git show 84e85c7:src/zefiro/l0/cycle.py 2>/dev/null | head -c 200 | od -c | grep -c '\\r' \
  && echo "  il primo commit ha CRLF" || echo "  il primo commit e' LF"
echo
echo "=== esiste gia' un .gitattributes? ==="
ls -la $R/.gitattributes 2>/dev/null || echo "  no"
