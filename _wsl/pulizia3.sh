#!/bin/bash
M=/mnt/i/AA_ENGINE
echo "=== gli STL: date, dimensioni, sono lo stesso file? ==="
ls -la $M/runs/mio/motore.stl $M/runs/prova_mat/zefiro50N.stl
echo "md5 (ci mette un po', sono 80 MB l'uno):"
md5sum $M/runs/mio/motore.stl $M/runs/prova_mat/zefiro50N.stl
echo
echo "=== si possono rigenerare? esiste l'esportatore? ==="
grep -n 'def build_and_export' $M/src/zefiro/geometry/*.py
echo "--- chi lo chiama ---"
grep -rln 'build_and_export' $M/scripts $M/tests | tr '\n' ' '; echo
echo
echo "=== il file dal nome storto ==="
ls $M/_wsl | grep -n 'rifianc'
