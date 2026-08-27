#!/usr/bin/env bash
# Prepara l'ambiente da zero e verifica che funzioni.
# Idempotente: rilanciarlo non rompe niente.
#
#   bash scripts/bootstrap.sh
set -euo pipefail

cd "$(dirname "$0")/.."
REPO="$(pwd)"
echo "repo: $REPO"

if [[ ! -f pyproject.toml ]]; then
  echo "ERRORE: non sembra la radice del repo (manca pyproject.toml)." >&2
  exit 1
fi

# --- 1. ambiente isolato ---------------------------------------------------
if [[ ! -d .venv ]]; then
  echo "==> creo l'ambiente virtuale"
  python3 -m venv .venv
fi
# shellcheck disable=SC1091
source .venv/bin/activate
python -m pip install --quiet --upgrade pip

# --- 2. dipendenze ---------------------------------------------------------
echo "==> installo (la prima volta scarica ~500 MB di OCCT, abbi pazienza)"
pip install --quiet --timeout 300 --retries 5 -e ".[dev]"

# --- 3. cartella delle run fuori da drvfs ----------------------------------
if [[ -z "${ZEFIRO_RUNS:-}" ]]; then
  RUNS="$HOME/zefiro-runs"
  mkdir -p "$RUNS"
  export ZEFIRO_RUNS="$RUNS"
  if ! grep -q 'ZEFIRO_RUNS' "$HOME/.bashrc" 2>/dev/null; then
    echo "export ZEFIRO_RUNS=\"$RUNS\"" >> "$HOME/.bashrc"
    echo "==> aggiunto ZEFIRO_RUNS a ~/.bashrc  ($RUNS)"
  fi
fi
echo "    ZEFIRO_RUNS=$ZEFIRO_RUNS"

# --- 4. git: serve per l'identita' delle run -------------------------------
if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "==> inizializzo il repository git"
  echo "    (il run_id include la revisione del codice: senza git non"
  echo "     distinguerebbe due versioni diverse del software)"
  git init -q
  git add -A
  git commit -qm "stato iniziale" || true
fi

# --- 5. verifica -----------------------------------------------------------
echo
python scripts/doctor.py
echo
echo "==> test"
pytest -q

cat <<'MSG'

Fatto. Ricorda di attivare l'ambiente a ogni nuova sessione:

    source .venv/bin/activate

Primo giro:

    python scripts/plant_report.py --fad 300 --bottle-T 20
    zefiro-l0
    zefiro-geometry
MSG
