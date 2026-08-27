"""Problemi strutturali: uno di verifica, uno vero.

`mbb_beam` non serve a Zefiro: serve a sapere se l'ottimizzatore funziona. E' il
caso di riferimento della letteratura sull'ottimizzazione topologica, e la
struttura che ne esce e' nota — un traliccio con un corrente superiore in
compressione, uno inferiore in trazione e diagonali fra i due. Se non esce
quella, il codice ha un problema, e conviene scoprirlo li' invece che sulla
geometria vera.

`chamber_jacket` e' il problema di Zefiro: la sezione meridiana del mantello,
con la pressione interna e il gradiente termico nella parete, e la domanda
"dove serve materiale". La risposta e' la nervatura, e non la decide chi scrive
il codice.
"""
from __future__ import annotations

import numpy as np

from zefiro.topopt.fem2d import AXISYMMETRIC, PLANE_STRESS, Grid
from zefiro.topopt.simp import Material, Problem


def mbb_beam(nx: int = 120, ny: int = 40, E: float = 1.0, nu: float = 0.3) -> Problem:
    """Trave MBB, mezza per simmetria. Caso di verifica dell'ottimizzatore.

    Carico verticale nell'angolo in alto a sinistra, simmetria sul bordo
    sinistro (spostamento orizzontale impedito), appoggio verticale nell'angolo
    in basso a destra. Unita' arbitrarie: qui interessa la FORMA, non i newton.
    """
    g = Grid(nx=nx, ny=ny, dx=1.0, dy=1.0, mode=PLANE_STRESS)
    f = np.zeros(g.n_dof)
    f[2 * g.node(0, ny) + 1] = -1.0
    fixed = np.array([2 * g.node(0, j) for j in range(ny + 1)]
                     + [2 * g.node(nx, 0) + 1])
    return Problem(grid=g, material=Material(E=E, nu=nu), fixed_dofs=fixed,
                   force=f, name="trave MBB (verifica)")


def chamber_jacket(
    r_inner: float,
    length: float,
    thickness_envelope: float,
    liner_thickness: float,
    p_chamber: float,
    delta_T_wall: float,
    material: Material,
    nr: int = 40,
    nz: int = 120,
) -> Problem:
    """Mantello della camera in sezione meridiana, assialsimmetrico.

    Dominio di progetto: l'anello fra il raggio interno e l'inviluppo massimo
    che accetti di occupare. Il LINER (lo strato piu' interno) e' congelato a
    pieno: e' il confine di pressione e la superficie bagnata, non e' negoziabile
    — l'ottimizzatore non puo' aprirci un buco. Tutto il resto e' libero.

    Carichi, entrambi presenti perche' agiscono insieme:
      * pressione interna sulla faccia del liner;
      * gradiente termico nello spessore, applicato come dilatazione impedita.
        E' quello che domina: a 5 bar la pressione produce megapascal, il
        gradiente ne produce centinaia.

    Vincoli: appoggio assiale alle due estremita' (la camera e' flangiata alla
    testa e trattenuta all'ugello). E' un'idealizzazione, e va rifatta quando il
    montaggio reale sara' deciso.
    """
    dr = thickness_envelope / nr
    dz = length / nz
    g = Grid(nx=nr, ny=nz, dx=dr, dy=dz, mode=AXISYMMETRIC, x0=r_inner)

    n_liner = max(1, int(round(liner_thickness / dr)))
    idx = np.arange(g.n_elem)
    i_col = idx // nz                        # indice radiale
    solid = idx[i_col < n_liner]             # liner congelato a pieno

    # pressione sulla faccia interna: forza = p * 2 pi r dz, ripartita sui nodi
    f = np.zeros(g.n_dof)
    for j in range(nz):
        for node, w in ((g.node(0, j), 0.5), (g.node(0, j + 1), 0.5)):
            f[2 * node] += p_chamber * 2.0 * np.pi * r_inner * dz * w

    # gradiente termico: lineare nello spessore, massimo sul liner
    dT = np.zeros(g.n_elem)
    for i in range(nr):
        frazione = 1.0 - i / max(nr - 1, 1)
        dT[i * nz:(i + 1) * nz] = delta_T_wall * frazione

    fixed = np.array(
        [2 * g.node(i, 0) + 1 for i in range(nr + 1)]
        + [2 * g.node(i, nz) + 1 for i in range(nr + 1)]
        + [2 * g.node(0, 0)]                       # blocca il moto rigido radiale
    )
    return Problem(grid=g, material=material, fixed_dofs=fixed, force=f,
                   delta_T=dT, solid=solid, name="mantello di camera")


def to_ascii(xi: np.ndarray, nx: int, ny: int, cols: int = 100,
             rows: int = 26, transpose: bool = False) -> str:
    """Anteprima testuale del risultato, per guardarlo senza aprire nulla."""
    X = xi.reshape(nx, ny)
    if transpose:
        X = X.T
    h, w = X.shape
    ci = np.linspace(0, w - 1, min(cols, w)).astype(int)
    ri = np.linspace(h - 1, 0, min(rows, h)).astype(int)
    palette = " .:-=+*#%@"
    out = []
    for r in ri:
        out.append("".join(palette[min(int(X[r, c] * len(palette)), len(palette) - 1)]
                           for c in ci))
    return "\n".join(out)
