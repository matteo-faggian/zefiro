#!/bin/bash
set +e
cp /mnt/i/AA_ENGINE/scripts/cfd/sezioni.py /mnt/i/AA_ENGINE/scripts/cfd/filmato.py /root/zefiro/scripts/cfd/
export PYTHONPATH=/root/zefiro/src
cd /root/zefiro
echo "=== sezioni ==="
/root/zef/bin/python scripts/cfd/sezioni.py --caso runs/cfd/run2/caso 2>&1 | tail -25
echo "=== un fotogramma ==="
/root/zef/bin/python scripts/cfd/filmato.py --caso runs/cfd/run2/caso --piano meridiano_getto 2>&1 | tail -5
ls runs/cfd/run2/caso/film_meridiano_getto_C3H8/ 2>/dev/null | tail -3
