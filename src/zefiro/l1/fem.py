"""FEM termo-meccanico con CalculiX. STUB."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from zefiro.schemas import CFDResult, FEMResult, GeometryArtifact


def run_calculix(
    geometry: GeometryArtifact,
    loads: CFDResult,
    material: dict[str, Any],
    case_dir: Path,
    transient: bool,
    burn_time: float,
) -> FEMResult:
    """Analisi termica (transitoria) seguita da analisi meccanica.

    `material` deve contenere le proprieta' del 316L SLM: se sono `null`
    (TODO n.6) l'implementazione dovra' sollevare MissingDatum, non usare
    valori da manuale del 316L laminato, che e' un materiale diverso.
    """
    raise NotImplementedError("Fase 4 della ROADMAP.")
