#!/bin/bash
pkill -f reactingFoam; pkill -f _lancia.sh
tr -d '\r' < /mnt/i/AA_ENGINE/_wsl/lancia.sh > /root/_lancia.sh
chmod +x /root/_lancia.sh
setsid nohup /root/_lancia.sh >> /root/zefiro/runs/cfd/lancio3.log 2>&1 < /dev/null &
sleep 6
tail -3 /root/zefiro/runs/cfd/lancio3.log 2>/dev/null
ps -eo pid,etime,cmd | grep -E '_lancia|python|gmsh' | grep -v grep | head -3
