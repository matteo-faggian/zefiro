"""Accoppiamento CFD -> FEM, one-way. STUB."""
from __future__ import annotations

from pathlib import Path


def map_wall_loads(
    source: Path, target_mesh: Path, out: Path, tolerance: float = 1.0e-3
) -> float:
    """Trasferisce p e q dalla mesh fluida a quella strutturale.

    Ritorna l'errore RELATIVO di conservazione dell'integrale:

        err = |INT q dA (CFD) - INT q dA (FEM)| / |INT q dA (CFD)|

    e l'implementazione deve FALLIRE se err > tolerance. Motivo: se
    l'integrale non si conserva, il FEM sta risolvendo un problema termico
    diverso da quello posto dalla CFD, e i margini strutturali che ne escono
    non significano nulla. E' l'unico controllo che rende l'accoppiamento
    verificabile invece che sperato.
    """
    raise NotImplementedError("Fase 4 della ROADMAP.")
