#!/usr/bin/env python3
"""Scarica il meccanismo cinetico San Diego e ne registra lo SHA256.

Perche' non e' vendorizzato nel repo: il file e' un DATO di terze parti con una
sua licenza e una sua versione. Metterlo in Git significherebbe che fra sei mesi
non sapresti piu' quale versione hai usato. Cosi' invece il nome del file porta
il proprio hash e `config/mechanisms.lock` registra fonte, data e hash.

Uso:
    python scripts/fetch_mechanism.py --out config/mechanisms/

Se la rete non raggiunge il sito UCSD, scarica il file a mano dalla pagina
ufficiale del gruppo Combustion Research di UC San Diego e passa
--from-file <percorso>. Lo script fara' comunque conversione, hashing e lock.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

# Pagina ufficiale del meccanismo (da cui prendere i file Chemkin):
SOURCE_PAGE = "https://web.eng.ucsd.edu/mae/groups/combustion/mechanism.html"


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", type=Path, default=Path("config/mechanisms"))
    ap.add_argument("--url", type=str, default=None,
                    help="URL diretto del file .yaml o .inp gia' convertito")
    ap.add_argument("--from-file", type=Path, default=None,
                    help="usa un file locale gia' scaricato a mano")
    ap.add_argument("--name", type=str, default="sandiego")
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    tmp = args.out / f"{args.name}.raw"

    if args.from_file is not None:
        shutil.copy(args.from_file, tmp)
        origin = str(args.from_file)
    elif args.url is not None:
        print(f"scarico {args.url}")
        urllib.request.urlopen(args.url, timeout=60)  # noqa: S310
        with urllib.request.urlopen(args.url, timeout=60) as r, tmp.open("wb") as fh:  # noqa: S310
            shutil.copyfileobj(r, fh)
        origin = args.url
    else:
        print(
            "Nessuna sorgente indicata.\n"
            f"Scarica i file Chemkin del meccanismo da:\n  {SOURCE_PAGE}\n"
            "poi rilancia con --from-file <percorso>, oppure passa --url.",
            file=sys.stderr,
        )
        return 2

    if tmp.suffix != ".yaml" and not tmp.read_bytes()[:64].lstrip().startswith(b"description"):
        # formato Chemkin: conversione con ck2yaml (fornito da Cantera)
        converted = args.out / f"{args.name}.yaml"
        subprocess.run(
            [sys.executable, "-m", "cantera.ck2yaml",
             f"--input={tmp}", f"--output={converted}", "--permissive"],
            check=True,
        )
    else:
        converted = args.out / f"{args.name}.yaml"
        shutil.copy(tmp, converted)

    digest = sha256_of(converted)
    final = args.out / f"{args.name}_{digest[:8]}.yaml"
    converted.rename(final)
    tmp.unlink(missing_ok=True)

    lock = Path("config/mechanisms.lock")
    entry = {
        "name": args.name,
        "file": final.name,
        "sha256": digest,
        "origin": origin,
        "source_page": SOURCE_PAGE,
        "fetched_utc": datetime.now(timezone.utc).isoformat(),
    }
    data = json.loads(lock.read_text()) if lock.exists() else {}
    data[args.name] = entry
    lock.write_text(json.dumps(data, indent=2), encoding="utf-8")

    print(f"OK -> {final}")
    print(f"sha256 {digest}")
    print("Ora imposta in config/operating_point.yaml:")
    print(f"  fuel.thermo_source: {final}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
