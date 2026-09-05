#!/bin/bash
tr -d '\r' < /mnt/i/AA_ENGINE/_wsl/geom.sh > /root/_geom.sh
chmod +x /root/_geom.sh
setsid nohup /root/_geom.sh > /root/zefiro/runs/geom.log 2>&1 < /dev/null &
sleep 5
ps -eo pid,etime,cmd | grep -E '[_]geom|[p]rogetta' | head -3
