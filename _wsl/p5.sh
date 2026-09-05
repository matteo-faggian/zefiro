#!/bin/bash
# Storia delle pressioni e delle portate ai bordi del run 4, letta dai file
# grezzi di surfaceFieldValue. Serve a capire SE e QUANTO le pressioni siano
# ancora in transitorio, prima di decidere come impostare il run 5.
P=/root/zefiro/runs/cfd/run4/caso/postProcessing
for d in p_ingresso_aria p_ingresso_gpl p_uscita portata_ingresso_aria portata_ingresso_gpl portata_uscita; do
  echo "### $d"
  find $P/$d -name '*.dat' | sort | head -5
done
echo
echo "=== esempio di intestazione ==="
f=$(find $P/p_uscita -name '*.dat' | sort | head -1)
head -6 "$f"
echo "..."
tail -3 "$f"
