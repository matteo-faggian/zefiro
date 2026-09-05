#!/bin/bash
M=/mnt/i/AA_ENGINE
echo "=== runs/mio ==="
du -sh $M/runs/mio/* 2>/dev/null | sort -h | tail -12
echo
echo "=== runs/prova_mat ==="
du -sh $M/runs/prova_mat/* 2>/dev/null | sort -h | tail -12
echo
echo "=== i due sono lo stesso contenuto? ==="
ls $M/runs/mio | tr '\n' ' '; echo
ls $M/runs/prova_mat | tr '\n' ' '; echo
echo
echo "=== scripts ==="
ls $M/scripts | tr '\n' ' '; echo
echo
echo "=== scripts/cfd ==="
ls $M/scripts/cfd | tr '\n' ' '; echo
echo
echo "=== _wsl: quali script sono ancora usati ==="
ls $M/_wsl | wc -l
echo
echo "=== .gitignore ==="
cat $M/.gitignore 2>/dev/null
