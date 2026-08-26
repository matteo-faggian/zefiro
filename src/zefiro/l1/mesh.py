"""Generazione mesh con Gmsh (API Python). STUB."""
from __future__ import annotations

from pathlib import Path

from zefiro.schemas import CANONICAL_BOUNDARIES, GeometryArtifact, MeshArtifact


def generate_mesh(
    geometry: GeometryArtifact,
    out_dir: Path,
    target_y_plus: float,
    n_boundary_layers: int,
    base_size: float,
) -> MeshArtifact:
    """Mesh del VOLUME FLUIDO a partire dallo STEP del solido.

    Nota di metodo: il dominio fluido e' il complemento del solido dentro un
    contenitore che include la regione di scarico. Non si mesha il solido: si
    mesha il negativo. Il FEM userA' invece il solido, con una mesh distinta.

    Precondizioni verificate dall'implementazione:
      * `geometry.is_watertight_mesh` e `geometry.is_valid_brep` entrambi veri;
      * i nomi delle physical groups prodotte coincidono ESATTAMENTE con
        `CANONICAL_BOUNDARIES` (un refuso qui = condizione al contorno
        sbagliata in silenzio).
    """
    raise NotImplementedError(
        "Fase 3 della ROADMAP. Contratto gia' fissato: MeshArtifact, "
        f"boundary attese = {CANONICAL_BOUNDARIES}"
    )
