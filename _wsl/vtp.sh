#!/bin/bash
F=$(find /root/zefiro/runs/cfd/smoke/postProcessing -name 'sez_02.vtp' | head -1)
echo "file: $F"
head -c 400 "$F"; echo
/root/zef/bin/python - "$F" <<'PY'
import sys
try:
    import meshio
    m = meshio.read(sys.argv[1])
    print("meshio ok:", m)
except Exception as e:
    print("meshio KO:", type(e).__name__, e)
PY
