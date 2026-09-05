#!/bin/bash
C=/root/zefiro/runs/cfd/run3/caso
grep -inE "unknown|not found|error|cannot|CourantNo|cellCentres|available" "$C/log.post" | head -40
