#!/bin/bash
C=/root/zefiro/runs/cfd/run3/caso
ls "$C/postProcessing/"
echo "=== log ==="
tail -25 "$C/log.fianco"
