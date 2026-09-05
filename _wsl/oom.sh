#!/bin/bash
echo "--- OOM nel log del kernel ---"
dmesg 2>/dev/null | grep -iE "out of memory|killed process|oom" | tail -8
echo "--- memoria ---"
free -g | head -2
