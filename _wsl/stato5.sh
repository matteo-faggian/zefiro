#!/bin/bash
# Stato del codice prima del run 5: la correzione della tasca c'e'? i test passano?
# e quali sono i numeri di progetto che servono a dimensionare il run.
cd /root/zefiro
echo "=== correzione tasca in mesh_iniettore.py ==="
grep -n 'SCONFINAMENTO_GETTO\|def tratto_del_getto' scripts/cfd/mesh_iniettore.py
echo
echo "=== piano meridiano_fianco nel generatore ==="
grep -n 'scarto_fianco' scripts/cfd/caso_mescolamento.py
echo
echo "=== test ==="
/root/zef/bin/python -m pytest -q 2>&1 | tail -5
echo
echo "=== numeri di progetto ==="
/root/zef/bin/python - <<'PY'
import sys; sys.path.insert(0, "scripts/cfd")
from mesh_iniettore import quote_dal_progetto
import math
d = quote_dal_progetto()
for k in sorted(d):
    v = d[k]
    print(f"  {k:22s} {v!r}" if not isinstance(v, float) else f"  {k:22s} {v:.6g}")
PY
