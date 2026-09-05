#!/bin/bash
# Che cosa c'e' davvero dentro gli zip, prima di ricostruirci una storia sopra.
Z="/mnt/i/AA_ENGINE/(eliminare)/_da_cancellare"
ls -la "$Z"/*.zip 2>/dev/null | wc -l
echo
echo "=== zip in ordine di data, con radice e numero di file ==="
for f in "$Z"/*.zip; do
  [ -e "$f" ] || continue
  d=$(stat -c '%y' "$f" | cut -c1-16)
  n=$(unzip -l "$f" 2>/dev/null | tail -1 | awk '{print $2}')
  radice=$(unzip -l "$f" 2>/dev/null | awk 'NR==4{print $4}')
  printf "%s  %-28s  %4s file  radice: %s\n" "$d" "$(basename "$f")" "$n" "$radice"
done | sort
echo
echo "=== e la cartella 'zefiro' scompattata che c'e' accanto? ==="
ls "$Z/zefiro" 2>/dev/null | head
echo
echo "=== esempio: che cosa contiene il piu' recente ==="
ultimo=$(ls -t "$Z"/*.zip 2>/dev/null | head -1)
echo "  $ultimo"
unzip -l "$ultimo" 2>/dev/null | head -20
