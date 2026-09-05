#!/bin/bash
# Installa OpenFOAM in WSL. Da lanciare con sudo, UNA VOLTA SOLA.
#
#   sudo bash /mnt/i/AA_ENGINE/scripts/cfd/installa_openfoam.sh
#
# PERCHE' NON LO SCRIPT UFFICIALE. `add-debian-repo.sh` di openfoam.com ricava
# il nome in codice della distribuzione da /etc/os-release, e su Ubuntu 26.04
# trova "resolute", per cui il repository non ha pacchetti. Qui si aggancia
# esplicitamente il ramo "noble" (Ubuntu 24.04 LTS): i binari compilati per
# noble girano senza problemi su una distribuzione piu' recente, perche' la
# glibc e' compatibile all'indietro. Il contrario non sarebbe vero.
set -euo pipefail

VERSIONE="${1:-openfoam2606-default}"
RAMO="noble"

if [ "$(id -u)" -ne 0 ]; then
    echo "Va lanciato con sudo:  sudo bash $0" >&2
    exit 1
fi

echo "== dipendenze minime"
apt-get update -qq
apt-get install -y --no-install-recommends curl gnupg ca-certificates

echo "== chiave del repository OpenFOAM"
curl -fsSL https://dl.openfoam.com/pubkey.gpg \
    | gpg --dearmor > /etc/apt/trusted.gpg.d/openfoam.gpg

echo "== sorgente apt, agganciata al ramo $RAMO"
cat > /etc/apt/sources.list.d/openfoam.list <<LISTA
# OpenFOAM (www.openfoam.com), ramo $RAMO scelto a mano: vedi installa_openfoam.sh
deb [arch=amd64] https://dl.openfoam.com/repos/deb $RAMO main
LISTA

echo "== installazione di $VERSIONE (qualche minuto, sono circa 2 GB)"
apt-get update -qq
apt-get install -y "$VERSIONE"

BASHRC="$(ls -d /usr/lib/openfoam/openfoam*/etc/bashrc 2>/dev/null | tail -1 || true)"
if [ -z "$BASHRC" ]; then
    echo "installato, ma non trovo il file etc/bashrc: controlla /usr/lib/openfoam" >&2
    exit 1
fi

echo
echo "== verifica"
# shellcheck disable=SC1090
source "$BASHRC"
for u in blockMesh checkMesh gmshToFoam decomposePar reactingFoam foamToVTK; do
    if command -v "$u" >/dev/null; then echo "  ok      $u"; else echo "  MANCA   $u"; fi
done

echo
echo "FATTO. Per usarlo in una shell nuova:"
echo "  source $BASHRC"
echo
echo "Riga da aggiungere al ~/.bashrc se vuoi che sia sempre pronto:"
echo "  echo 'source $BASHRC' >> ~/.bashrc"
