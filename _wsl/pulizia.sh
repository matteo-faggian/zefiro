#!/bin/bash
# Inventario della cartella di progetto: cosa c'e' e quanto pesa.
# NON cancella niente: serve solo a decidere.
M=/mnt/i/AA_ENGINE
echo "=== peso delle cartelle di primo livello ==="
du -sh $M/* $M/.[a-z]* 2>/dev/null | sort -h
echo
echo "=== _da_cancellare ==="
ls -la $M/_da_cancellare 2>/dev/null | head -30
echo
echo "=== runs (primo livello) ==="
du -sh $M/runs/* 2>/dev/null | sort -h | tail -20
echo
echo "=== runs_cfd (primo livello) ==="
du -sh $M/runs_cfd/* 2>/dev/null | sort -h | tail -20
echo
echo "=== docs ==="
ls -la $M/docs 2>/dev/null | head -30
