#!/usr/bin/env python3
"""Guarda la mesh: sezione meridiana, faccia d'iniezione, qualita'.

Una mesh si giudica guardandola. I numeri di `checkMesh` dicono se e'
accettabile, non se e' la mesh che volevi: un dominio costruito male, un foro
finito nel posto sbagliato, un raffinamento che non copre la zona giusta
passano tutti i controlli numerici.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


def leggi(msh: Path):
    import gmsh

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.merge(str(msh))
        tag, coord, _ = gmsh.model.mesh.getNodes()
        xyz = np.array(coord, dtype=float).reshape(-1, 3)
        idx = np.zeros(int(np.max(tag)) + 1, dtype=np.int64)
        idx[np.array(tag, dtype=np.int64)] = np.arange(len(tag))
        tipi, _, nodi = gmsh.model.mesh.getElements(3)
        tet = np.array(nodi[list(tipi).index(4)], dtype=np.int64).reshape(-1, 4)
        patch = {}
        for dim, ptag in gmsh.model.getPhysicalGroups(2):
            nome = gmsh.model.getPhysicalName(dim, ptag)
            tri = []
            for sup in gmsh.model.getEntitiesForPhysicalGroup(dim, ptag):
                tp, _, nd = gmsh.model.mesh.getElements(2, sup)
                if 2 in tp:
                    tri.append(np.array(nd[list(tp).index(2)],
                                        dtype=np.int64).reshape(-1, 3))
            if tri:
                patch[nome] = idx[np.vstack(tri)]
        return xyz, idx[tet], patch
    finally:
        gmsh.finalize()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("msh", type=Path)
    ap.add_argument("--out", type=Path, default=Path("runs/mesh/mesh.png"))
    a = ap.parse_args()

    xyz, tet, patch = leggi(a.msh)
    print(f"{len(tet)} tetraedri, {len(xyz)} nodi, {len(patch)} patch")
    for nome, tri in sorted(patch.items()):
        P = xyz[tri]
        area = float(np.linalg.norm(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]),
                                    axis=1).sum() / 2.0)
        print(f"  {nome:<20}{len(tri):7d} triangoli   {area*1e6:10.2f} mm2")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PolyCollection

    fig, axes = plt.subplots(1, 3, figsize=(19, 6.2))

    # --- 1. traccia della mesh sul piano periodico ------------------------- #
    ax = axes[0]
    tri = patch.get("periodic_a")
    if tri is not None:
        P = xyz[tri]
        seg = []
        for i, j in ((0, 1), (1, 2), (2, 0)):
            seg.extend(zip(zip(P[:, i, 0] * 1e3, np.hypot(P[:, i, 1], P[:, i, 2]) * 1e3),
                           zip(P[:, j, 0] * 1e3, np.hypot(P[:, j, 1], P[:, j, 2]) * 1e3)))
        ax.add_collection(LineCollection(seg, linewidths=0.25, colors="#2b6cb0"))
        ax.autoscale()
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_title("Sezione meridiana (faccia periodica)", fontsize=11)

    # --- 2. zoom sulla camera, con le pareti evidenziate ------------------- #
    ax = axes[1]
    if tri is not None:
        ax.add_collection(LineCollection(seg, linewidths=0.4, colors="#a0aec0"))
    colori = {"wall_chamber": "#c53030", "wall_convergent": "#dd6b20",
              "wall_plug": "#2f855a", "wall_faceplate": "#6b46c1",
              "inlet_air": "#2b6cb0", "inlet_fuel_core": "#000000",
              "inlet_fuel_film": "#b7791f"}
    for nome, col in colori.items():
        if nome not in patch:
            continue
        P = xyz[patch[nome]]
        cx = P[:, :, 0].mean(axis=1) * 1e3
        cr = np.hypot(P[:, :, 1], P[:, :, 2]).mean(axis=1) * 1e3
        ax.scatter(cx, cr, s=2.5, color=col, label=nome)
    lim_x = xyz[:, 0].max() * 1e3
    ax.set_xlim(-2, 0.28 * lim_x)
    ax.set_ylim(0, np.hypot(xyz[:, 1], xyz[:, 2]).max() * 1e3 * 0.45)
    ax.set_aspect("equal")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_title("Camera e frontiere nominate", fontsize=11)
    ax.legend(fontsize=8, markerscale=3, loc="upper right", frameon=False)

    # --- 3. la faccia d'iniezione vista di fronte -------------------------- #
    ax = axes[2]
    for nome, col, z in (("wall_faceplate", "#e2e8f0", 1), ("inlet_air", "#2b6cb0", 2),
                         ("inlet_fuel_core", "#c53030", 3),
                         ("inlet_fuel_film", "#b7791f", 3)):
        if nome not in patch:
            continue
        P = xyz[patch[nome]][:, :, 1:] * 1e3
        ax.add_collection(PolyCollection(P, facecolors=col, edgecolors="#4a5568",
                                         linewidths=0.15, zorder=z, label=nome))
    ax.autoscale()
    ax.set_aspect("equal")
    ax.set_xlabel("y [mm]")
    ax.set_ylabel("z [mm]")
    ax.set_title("Faccia d'iniezione: i fori sono imprintati", fontsize=11)
    ax.legend(fontsize=8, loc="upper right", frameon=False)

    a.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(a.out, dpi=135, bbox_inches="tight")
    print(f"\nscritto {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
