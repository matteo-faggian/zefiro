#!/bin/bash
C=/root/zefiro/runs/cfd/smoke
echo "--- alberatura postProcessing ---"
find "$C/postProcessing" -maxdepth 3 -type d 2>/dev/null | head -20
echo "--- file ---"
find "$C/postProcessing" -type f 2>/dev/null | head -20
echo "--- conteggio ---"
find "$C/postProcessing" -type f 2>/dev/null | wc -l
