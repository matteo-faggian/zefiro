"""Assemblaggio di obiettivi e vincoli. STUB parzialmente implementato.

Convenzione RIGIDA (docs/architettura.md sezione 4.7): minimizzare `f`,
fattibile se `g <= 0`. L'inversione di segno la fa QUESTO modulo, mai
l'ottimizzatore.
"""
from __future__ import annotations

from typing import Any

from zefiro.schemas import FEMResult, GeometryArtifact, L0Result, Objectives


def objectives_l0(
    run_id: str, l0: L0Result, geometry: GeometryArtifact | None = None
) -> Objectives:
    """Obiettivi e vincoli valutabili gia' a L0 (millisecondi).

    Serve a scartare il grosso dei candidati prima di spendere una run L1.
    """
    f = {
        "neg_thrust": -l0.thrust,          # massimizzare spinta
        "neg_isp": -l0.Isp_s,              # massimizzare Isp di sistema
    }
    g = {
        # il combustibile deve poter entrare in camera con un Dp del 20 %
        "fuel_dp_margin": l0.phi_core * 0.0,   # placeholder: richiede Cd (TODO n.5)
    }
    if geometry is not None:
        f["mass"] = geometry.volume            # minimizzare massa (a densita' fissa)
        g["watertight"] = 0.0 if geometry.is_watertight_mesh else 1.0
    return Objectives(
        run_id=run_id, fidelity="L0", f=f, g=g,
        source={k: "l0.cycle" for k in f},
    )


def objectives_l1(
    run_id: str, l0: L0Result, fem: FEMResult, extra: dict[str, Any] | None = None
) -> Objectives:
    """Obiettivi completi, con i margini strutturali e termici dal FEM."""
    raise NotImplementedError("Fase 5 della ROADMAP.")
