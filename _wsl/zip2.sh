#!/bin/bash
/root/zef/bin/python - <<'PY'
import zipfile, pathlib, datetime
Z = pathlib.Path("/mnt/i/AA_ENGINE/(eliminare)/_da_cancellare")
zips = sorted(Z.glob("*.zip"), key=lambda p: p.stat().st_mtime)
print(f"{'data':16} {'nome':24} {'file':>5} {'radici'}")
info = []
for z in zips:
    t = datetime.datetime.fromtimestamp(z.stat().st_mtime)
    try:
        with zipfile.ZipFile(z) as f:
            nomi = [n for n in f.namelist() if not n.endswith("/")]
    except Exception as e:
        print(f"  {z.name}: NON LEGGIBILE ({e})"); continue
    radici = sorted({n.split("/")[0] for n in nomi})
    info.append((t, z.name, nomi, radici))
    print(f"{t:%Y-%m-%d %H:%M} {z.name:24} {len(nomi):5d} {radici[:4]}")

print()
print("=== istantanee complete o consegne parziali? ===")
print("Si guarda se un file che c'e' SEMPRE nel progetto (il modello L0)")
print("compare in ognuno degli zip. Se manca in alcuni, sono consegne")
print("parziali e ricostruire una storia sostituendo l'albero farebbe")
print("SPARIRE dei file a ogni commit: sarebbe una storia falsa.")
for t, nome, nomi, _ in info:
    ha_l0 = any("l0/cycle.py" in n for n in nomi)
    ha_src = any(n.startswith(("src/", "zefiro/src/")) for n in nomi)
    print(f"  {nome:24} cycle.py {'SI' if ha_l0 else 'no':>3}   "
          f"src/ {'SI' if ha_src else 'no':>3}   {len(nomi):4d} file")
PY
