#!/bin/bash
cd /mnt/i/AA_ENGINE || exit 2
git status --short | head -20
echo "=== branch ==="
git rev-parse --abbrev-ref HEAD
