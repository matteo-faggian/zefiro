#!/bin/bash
# Sposta in "(eliminare)" quello che e' sicuramente ridondante.
# NON cancella niente: lo spostamento e' reversibile, e la struttura di
# origine si conserva dentro la cartella cosi' si sa sempre da dove veniva.
set -e
M=/mnt/i/AA_ENGINE
E="$M/(eliminare)"
mkdir -p "$E"

sposta() {   # $1 = percorso relativo a M
  local src="$M/$1"
  [ -e "$src" ] || { echo "  (assente) $1"; return 0; }
  local dst="$E/$(dirname "$1")"
  mkdir -p "$dst"
  mv "$src" "$dst/"
  echo "  spostato: $1"
}

echo "=== 1. duplicato esatto dello STL (stesso md5 dell'altro) ==="
sposta "runs/prova_mat"

echo "=== 2. cartella gia' marcata da te ==="
sposta "_da_cancellare"

echo "=== 3. cache rigenerate da sole ==="
sposta ".pytest_cache"
sposta ".ruff_cache"

echo "=== 4. script con il nome scritto in cirillico (la 'o' di 'fianco' e'"
echo "       U+043E: non e' invocabile come 'rifianco', e' nato morto) ==="
sposta "_wsl/rifiancо.sh"

echo
echo "=== risultato ==="
du -sh "$E"
echo "--- cosa resta nel progetto ---"
du -sh $M/* 2>/dev/null | sort -h | tail -8
