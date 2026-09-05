#!/bin/bash
# Verifica che il caso 5 sia davvero pronto PRIMA di impegnare la macchina
# per ore: decomposizione, tipi dei bordi, qualita' della mesh.
C=/root/zefiro/runs/cfd/run5/caso
echo "processori decomposti: $(ls -d $C/processor* 2>/dev/null | wc -l)"
echo
echo "=== tipi dei bordi (i due fianchi DEVONO essere symmetryPlane) ==="
awk '/^[a-z_]+$/{n=$1} /type/{print n" -> "$2} /nFaces/{print "    facce "$2}' \
    $C/constant/polyMesh/boundary | head -40
echo
echo "=== qualita' della mesh ==="
grep -E 'cells:|non-orthogonality|skewness|Mesh OK|\*\*\*' $C/log.checkMesh
echo
echo "=== confronto con la corsa 4 ==="
grep -E '^    cells:' /root/zefiro/runs/cfd/run4/caso/log.checkMesh
echo
echo "=== il getto non tocca piu' il metallo esterno? ==="
grep -n 'r0, lung = tratto_del_getto' /root/zefiro/scripts/cfd/mesh_iniettore.py
/root/zef/bin/python - <<'PY'
import sys; sys.path.insert(0, "/root/zefiro/scripts/cfd")
from mesh_iniettore import quote_dal_progetto, tratto_del_getto
d = quote_dal_progetto()
r0, lung = tratto_del_getto(d["r_bore"], d["R_getti"], d["R_anello"])
print(f"  cilindro del getto: da r = {r0*1e3:.3f} a {(r0+lung)*1e3:.3f} mm")
print(f"  parete esterna dell'anello a r = {d['R_anello']*1e3:.3f} mm")
print(f"  margine: {(d['R_anello']-(r0+lung))*1e3:.3f} mm  "
      f"({'OK, nessuna tasca' if r0+lung < d['R_anello'] else 'TASCA!'})")
PY
