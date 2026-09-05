"""FEM lineare elastico su griglia strutturata Q4: piano e assialsimmetrico.

E' il motore dell'ottimizzazione topologica. Deve essere ESATTO, perche' un
errore qui non produce un risultato sbagliato in modo evidente: produce una
struttura ottimizzata alla perfezione per un problema che non esiste.

Scelte, e i loro motivi:

* **La matrice di rigidezza dell'elemento e' INTEGRATA, non copiata.** Le
  versioni tabellate che girano in letteratura valgono per un quadrato unitario
  in tensione piana con nu = 0.3; qui servono elementi rettangolari, tensione e
  deformazione piana, e assialsimmetrico. Integrare con Gauss 2x2 costa
  microsecondi una volta sola e vale per tutti i casi.
* **Assialsimmetrico incluso.** La camera di Zefiro e' un solido di rivoluzione:
  la sezione meridiana con carico di pressione e gradiente termico e' il
  problema strutturale vero, non un'astrazione 2D.
* **Verifica con il patch test.** E' il controllo standard per un elemento
  finito: un campo di deformazione costante deve essere riprodotto ESATTAMENTE,
  a meno dell'errore di macchina. Un elemento che non lo passa non converge.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

PLANE_STRESS = "plane_stress"
PLANE_STRAIN = "plane_strain"
AXISYMMETRIC = "axisymmetric"

#: Punti e pesi di Gauss 2x2 (esatti per polinomi fino al terzo grado).
_G = 1.0 / np.sqrt(3.0)
GAUSS = [(-_G, -_G), (_G, -_G), (_G, _G), (-_G, _G)]


def constitutive(E: float, nu: float, mode: str) -> np.ndarray:
    """Matrice costitutiva D.

    3x3 nei casi piani (eps_x, eps_y, gamma_xy), 4x4 in assialsimmetrico, dove
    si aggiunge la deformazione circonferenziale eps_theta = u_r / r. E' proprio
    quel termine a rendere il problema assialsimmetrico diverso da uno piano:
    un cilindro in pressione resiste in circonferenziale, non in radiale.
    """
    if mode == PLANE_STRESS:
        f = E / (1.0 - nu * nu)
        return f * np.array([[1.0, nu, 0.0],
                             [nu, 1.0, 0.0],
                             [0.0, 0.0, (1.0 - nu) / 2.0]])
    if mode == PLANE_STRAIN:
        f = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
        return f * np.array([[1.0 - nu, nu, 0.0],
                             [nu, 1.0 - nu, 0.0],
                             [0.0, 0.0, (1.0 - 2.0 * nu) / 2.0]])
    if mode == AXISYMMETRIC:
        f = E / ((1.0 + nu) * (1.0 - 2.0 * nu))
        return f * np.array([[1.0 - nu, nu, 0.0, nu],
                             [nu, 1.0 - nu, 0.0, nu],
                             [0.0, 0.0, (1.0 - 2.0 * nu) / 2.0, 0.0],
                             [nu, nu, 0.0, 1.0 - nu]])
    raise ValueError(f"modo non riconosciuto: {mode!r}")


def _shape(xi: float, eta: float) -> tuple[np.ndarray, np.ndarray]:
    """Funzioni di forma Q4 bilineari e loro derivate in coordinate naturali.

    Ordine dei nodi, antiorario: (-1,-1), (+1,-1), (+1,+1), (-1,+1).
    """
    N = 0.25 * np.array([(1 - xi) * (1 - eta), (1 + xi) * (1 - eta),
                         (1 + xi) * (1 + eta), (1 - xi) * (1 + eta)])
    dN = 0.25 * np.array([
        [-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)],   # d/dxi
        [-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)],       # d/deta
    ])
    return N, dN


def element_stiffness(dx: float, dy: float, D: np.ndarray, mode: str,
                      r_center: float = 0.0) -> np.ndarray:
    """Rigidezza 8x8 di un elemento rettangolare, per integrazione di Gauss.

    `r_center` e' il raggio del centro dell'elemento e serve solo in
    assialsimmetrico, dove il volume dell'anello e' 2*pi*r*dA: gli elementi
    lontani dall'asse pesano di piu'.
    """
    J = np.array([[dx / 2.0, 0.0], [0.0, dy / 2.0]])
    detJ = dx * dy / 4.0
    Jinv = np.array([[2.0 / dx, 0.0], [0.0, 2.0 / dy]])
    K = np.zeros((8, 8))
    for xi, eta in GAUSS:
        N, dN = _shape(xi, eta)
        dNxy = Jinv @ dN                       # derivate in coordinate fisiche
        if mode == AXISYMMETRIC:
            r = r_center + 0.5 * dx * xi
            if r <= 0.0:
                r = 1.0e-9                     # elementi sull'asse: evita 1/0
            B = np.zeros((4, 8))
            B[0, 0::2] = dNxy[0]               # eps_r  = du/dr
            B[1, 1::2] = dNxy[1]               # eps_z  = dw/dz
            B[2, 0::2] = dNxy[1]               # gamma  = du/dz + dw/dr
            B[2, 1::2] = dNxy[0]
            B[3, 0::2] = N / r                 # eps_th = u/r
            w = detJ * 2.0 * np.pi * r
        else:
            B = np.zeros((3, 8))
            B[0, 0::2] = dNxy[0]
            B[1, 1::2] = dNxy[1]
            B[2, 0::2] = dNxy[1]
            B[2, 1::2] = dNxy[0]
            w = detJ
        K += w * (B.T @ D @ B)
    return K


def element_thermal_load(dx: float, dy: float, D: np.ndarray, mode: str,
                         alpha: float, delta_T: float,
                         r_center: float = 0.0) -> np.ndarray:
    """Vettore dei carichi termici dell'elemento: f = INT B^T D eps_0 dV.

    `eps_0 = alpha * dT` e' la deformazione libera impedita. E' il carico che
    domina in Zefiro: a 5.22 bar la pressione produce 6 MPa, un gradiente di
    100 K in parete impedita ne produce circa 460.
    """
    Jinv = np.array([[2.0 / dx, 0.0], [0.0, 2.0 / dy]])
    detJ = dx * dy / 4.0
    n = 4 if mode == AXISYMMETRIC else 3
    eps0 = np.zeros(n)
    eps0[:2] = alpha * delta_T
    if mode == AXISYMMETRIC:
        eps0[3] = alpha * delta_T
    f = np.zeros(8)
    for xi, eta in GAUSS:
        N, dN = _shape(xi, eta)
        dNxy = Jinv @ dN
        if mode == AXISYMMETRIC:
            r = max(r_center + 0.5 * dx * xi, 1.0e-9)
            B = np.zeros((4, 8))
            B[0, 0::2] = dNxy[0]; B[1, 1::2] = dNxy[1]
            B[2, 0::2] = dNxy[1]; B[2, 1::2] = dNxy[0]
            B[3, 0::2] = N / r
            w = detJ * 2.0 * np.pi * r
        else:
            B = np.zeros((3, 8))
            B[0, 0::2] = dNxy[0]; B[1, 1::2] = dNxy[1]
            B[2, 0::2] = dNxy[1]; B[2, 1::2] = dNxy[0]
            w = detJ
        f += w * (B.T @ D @ eps0)
    return f


@dataclass
class Grid:
    """Griglia strutturata di elementi rettangolari.

    Numerazione: elemento (i, j) con i lungo x (colonne) e j lungo y (righe);
    indice lineare `e = i * ny + j`, cioe' colonna per colonna. E' la stessa
    convenzione dei codici classici di ottimizzazione topologica, e serve a
    poter confrontare i risultati con la letteratura senza rimescolare indici.
    """
    nx: int
    ny: int
    dx: float
    dy: float
    mode: str = PLANE_STRESS
    x0: float = 0.0     # coordinata del bordo sinistro (raggio, se assialsim.)

    @property
    def n_elem(self) -> int:
        return self.nx * self.ny

    @property
    def n_node(self) -> int:
        return (self.nx + 1) * (self.ny + 1)

    @property
    def n_dof(self) -> int:
        return 2 * self.n_node

    def node(self, i: int, j: int) -> int:
        return i * (self.ny + 1) + j

    def element_dofs(self) -> np.ndarray:
        """(n_elem, 8) gradi di liberta' di ogni elemento."""
        out = np.empty((self.n_elem, 8), dtype=np.int64)
        for i in range(self.nx):
            for j in range(self.ny):
                n0 = self.node(i, j)
                n1 = self.node(i + 1, j)
                n2 = self.node(i + 1, j + 1)
                n3 = self.node(i, j + 1)
                out[i * self.ny + j] = [2 * n0, 2 * n0 + 1, 2 * n1, 2 * n1 + 1,
                                        2 * n2, 2 * n2 + 1, 2 * n3, 2 * n3 + 1]
        return out

    def element_centers(self) -> np.ndarray:
        """(n_elem, 2) coordinate dei centri."""
        i = np.repeat(np.arange(self.nx), self.ny)
        j = np.tile(np.arange(self.ny), self.nx)
        return np.column_stack([self.x0 + (i + 0.5) * self.dx, (j + 0.5) * self.dy])

    def element_stiffness_set(self, D: np.ndarray) -> np.ndarray:
        """(n_elem, 8, 8). In assialsimmetrico varia con il raggio, altrimenti
        e' la stessa matrice per tutti e la si calcola una volta sola."""
        if self.mode != AXISYMMETRIC:
            ke = element_stiffness(self.dx, self.dy, D, self.mode)
            return np.broadcast_to(ke, (self.n_elem, 8, 8)).copy()
        out = np.empty((self.n_elem, 8, 8))
        for i in range(self.nx):
            rc = self.x0 + (i + 0.5) * self.dx
            ke = element_stiffness(self.dx, self.dy, D, self.mode, r_center=rc)
            out[i * self.ny:(i + 1) * self.ny] = ke
        return out

    def element_thermal_set(self, D: np.ndarray, alpha: float,
                            delta_T: np.ndarray | float) -> np.ndarray:
        """(n_elem, 8) carichi termici, con dT eventualmente diverso per elemento."""
        dT = np.full(self.n_elem, float(delta_T)) if np.isscalar(delta_T) else np.asarray(delta_T)
        out = np.empty((self.n_elem, 8))
        for i in range(self.nx):
            rc = self.x0 + (i + 0.5) * self.dx
            base = element_thermal_load(self.dx, self.dy, D, self.mode, alpha, 1.0,
                                        r_center=rc)
            sl = slice(i * self.ny, (i + 1) * self.ny)
            out[sl] = base[None, :] * dT[sl, None]
        return out


def assemble(grid: Grid, ke_set: np.ndarray, density_scale: np.ndarray,
             edofs: np.ndarray) -> sp.csc_matrix:
    """Assembla K = SOMMA_e  scale_e * KE_e, in formato sparso."""
    n = ke_set.shape[0]
    rows = np.repeat(edofs, 8, axis=1).ravel()
    cols = np.tile(edofs, (1, 8)).ravel()
    vals = (ke_set * density_scale[:, None, None]).reshape(n, 64).ravel()
    K = sp.coo_matrix((vals, (rows, cols)), shape=(grid.n_dof, grid.n_dof))
    return K.tocsc()


def solve(K: sp.csc_matrix, f: np.ndarray, fixed: np.ndarray) -> np.ndarray:
    """Risolve K u = f con i gradi di liberta' vincolati eliminati."""
    free = np.setdiff1d(np.arange(K.shape[0]), fixed, assume_unique=False)
    u = np.zeros(K.shape[0])
    Kff = K[free, :][:, free]
    u[free] = spla.spsolve(Kff.tocsc(), f[free])
    return u
