#!/bin/bash
# dispatcher: converte CRLF ed esegue lo script richiesto
n="$1"; shift
tr -d '\r' < "/mnt/i/AA_ENGINE/_wsl/$n.sh" > "/root/_$n.sh"
chmod +x "/root/_$n.sh"
exec bash "/root/_$n.sh" "$@"
