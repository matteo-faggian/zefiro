"""RANS reattivo con OpenFOAM. STUB."""
from __future__ import annotations

from pathlib import Path

from zefiro.schemas import CFDResult, L0Result, MeshArtifact, OperatingPoint


def run_reacting_foam(
    mesh: MeshArtifact,
    op: OperatingPoint,
    l0: L0Result,
    case_dir: Path,
    n_procs: int,
    max_iterations: int,
) -> CFDResult:
    """Lancia il caso reattivo e ne estrae p e q a parete.

    `l0` non e' decorativo: fornisce la condizione iniziale (T_ad, composizione
    di equilibrio) da cui partire. Inizializzare a T ambiente un caso reattivo
    a queste scale porta quasi sempre allo spegnimento numerico.
    """
    raise NotImplementedError("Fase 3 della ROADMAP.")
