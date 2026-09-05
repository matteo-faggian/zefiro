#!/bin/bash
pkill -f reactingFoam; pkill -f _lancia.sh; sleep 2
rm -rf /root/zefiro/runs/cfd/run3
rm -f /root/zefiro/runs/cfd/lancio3.log
tr -d '\r' < /mnt/i/AA_ENGINE/_wsl/lancia.sh > /root/_lancia.sh
chmod +x /root/_lancia.sh
setsid nohup /root/_lancia.sh >> /root/zefiro/runs/cfd/lancio3.log 2>&1 < /dev/null &
sleep 8
ps -eo pid,etime,cmd | grep -E '_lancia|python' | grep -v grep | head -2
