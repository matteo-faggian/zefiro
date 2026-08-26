"""Identita' deterministica di una run (docs/architettura.md sezione 5.1).

run_id = blake2b_16( canonical_json({schema, design, operating, code_rev}) )

Proprieta' volute:
  * stesso input -> stesso id su qualunque macchina;
  * cambiare il codice cambia gli id, quindi non si confrontano per sbaglio
    risultati prodotti da versioni diverse;
  * working tree sporco -> suffisso "-dirty", e quelle run vanno in una tabella
    separata. Servono a esplorare, non a concludere.
"""
from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path
from typing import Any

from zefiro.schemas import SCHEMA_VERSION, DesignVector, OperatingPoint


def _canonical(obj: Any) -> Any:
    """JSON canonico: chiavi ordinate, float via repr() (round-trip esatto)."""
    if isinstance(obj, dict):
        return {k: _canonical(obj[k]) for k in sorted(obj, key=str)}
    if isinstance(obj, (list, tuple)):
        return [_canonical(v) for v in obj]
    if isinstance(obj, float):
        return repr(obj)
    if isinstance(obj, Path):
        return str(obj)
    return obj


def git_revision(repo: Path | None = None) -> str:
    """`git describe --always --dirty`, oppure 'nogit' se non e' un repo."""
    repo = repo or Path(__file__).resolve().parents[3]
    try:
        out = subprocess.run(
            ["git", "describe", "--always", "--dirty", "--tags"],
            cwd=repo, capture_output=True, text=True, timeout=10, check=True,
        )
        return out.stdout.strip() or "nogit"
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return "nogit"


def make_run_id(x: DesignVector, op: OperatingPoint, code_rev: str | None = None) -> str:
    rev = code_rev if code_rev is not None else git_revision()
    payload = _canonical({
        "schema_version": SCHEMA_VERSION,
        "design": dict(x.values),
        "operating": op.to_dict(),
        "code_rev": rev,
    })
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    rid = hashlib.blake2b(blob, digest_size=16).hexdigest()
    return f"{rid}-dirty" if rev.endswith("-dirty") else rid


def env_fingerprint() -> dict[str, str]:
    """Cio' che rende falsificabile l'affermazione 'riproducibile'."""
    import cantera as ct

    fp: dict[str, str] = {
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "machine": platform.machine(),
        "cantera": ct.__version__,
        "code_rev": git_revision(),
    }
    for mod in ("numpy", "build123d", "gmsh"):
        try:
            fp[mod] = __import__(mod).__version__
        except Exception:
            fp[mod] = "assente"
    return fp
