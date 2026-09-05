#!/bin/bash
M=/mnt/i/AA_ENGINE
echo "=== esiste .git nel master? com'e' fatto? ==="
ls -la "$M/.git" 2>/dev/null | head -12
echo
echo "=== git lo riconosce? ==="
cd "$M" && git rev-parse --git-dir 2>&1 | head -3
echo
echo "=== HEAD e config ci sono? ==="
head -1 "$M/.git/HEAD" 2>/dev/null || echo "  HEAD assente"
ls "$M/.git/refs/heads/" 2>/dev/null || echo "  refs/heads assente"
echo
echo "=== il repo di partenza sta bene? ==="
cd /root/zefiro_git && git log --oneline | head -2 && git status --short | head -3
