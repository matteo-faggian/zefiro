#!/bin/bash
C=/root/zefiro/runs/cfd/run2/caso
sed -n '/Checking geometry/,$p' "$C/log.checkMesh"
