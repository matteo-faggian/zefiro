"""Obiettivo basato sulla TENSIONE, non sulla compliance.

Perche' serve, e perche' la compliance non basta. Con un carico meccanico fisso,
minimizzare la compliance e' un buon surrogato di "non rompersi": la struttura
piu' rigida e' anche quella meno sollecitata, a parita' di materiale. Con un
carico TERMICO il ragionamento salta: il carico e' la dilatazione impedita,
quindi e' proporzionale al materiale presente. Aggiungere materia aggiunge
rigidezza *e* carico, e la compliance puo' benissimo salire. Su Zefiro succede
davvero — provato: minimizzando la compliance sotto gradiente termico
l'ottimizzatore peggiora l'obiettivo e svuota il dominio.

L'obiettivo giusto e' la tensione. Qui si minimizza la norma-P della tensione di
von Mises normalizzata all'ammissibile, che e' un'approssimazione derivabile del
massimo:

    J = [ SOMMA_e ( sigma_vm,e / sigma_amm )^P ]^(1/P)

Due accorgimenti obbligatori, entrambi ben noti e entrambi facili da dimenticare:

* **Rilassamento della tensione.** Un elemento con densita' che tende a zero ha
  tensione che NON tende a zero (il rapporto forza/area resta finito), e questo
  crea minimi locali dove l'ottimizzatore non riesce a svuotare un elemento
  perche' la sua tensione esplode. Si usa `sigma_rilassata = xi^eta * sigma`,
  con eta < 1, che rende il vincolo inattivo nel vuoto.
* **La tensione va calcolata sulla deformazione ELASTICA**, cioe' sottraendo la
  dilatazione libera: sigma = D (B u - eps_termica). Usare B u e basta significa
  attribuire al materiale una tensione che non prova.

La derivata e' ottenuta per variabile aggiunta e verificata contro differenze
finite nei test: e' l'unico modo di sapere che e' giusta.
"""
from __future__ import annotations

import numpy as np

from zefiro.topopt.fem2d import AXISYMMETRIC, Grid, assemble, solve

#: Esponente di rilassamento della tensione nel vuoto.
STRESS_RELAX_ETA = 0.5
#: Esponente di aggregazione. Alto approssima meglio il massimo ma peggiora il
#: condizionamento; 8-12 e' l'intervallo praticabile.
PNORM = 10.0


def _B_centroid(grid: Grid, i: int) -> np.ndarray:
    """Matrice B valutata al centro dell'elemento della colonna radiale `i`."""
    dNxy = np.array([[-1.0, 1.0, 1.0, -1.0], [-1.0, -1.0, 1.0, 1.0]]) * 0.25
    dNxy = np.array([[2.0 / grid.dx, 0.0], [0.0, 2.0 / grid.dy]]) @ dNxy
    if grid.mode == AXISYMMETRIC:
        r = max(grid.x0 + (i + 0.5) * grid.dx, 1e-9)
        B = np.zeros((4, 8))
        B[0, 0::2] = dNxy[0]; B[1, 1::2] = dNxy[1]
        B[2, 0::2] = dNxy[1]; B[2, 1::2] = dNxy[0]
        B[3, 0::2] = 0.25 / r
        return B
    B = np.zeros((3, 8))
    B[0, 0::2] = dNxy[0]; B[1, 1::2] = dNxy[1]
    B[2, 0::2] = dNxy[1]; B[2, 1::2] = dNxy[0]
    return B


def _vm_matrix(n: int) -> np.ndarray:
    """Matrice V tale che sigma_vm^2 = sigma' V sigma."""
    if n == 3:
        return np.array([[1.0, -0.5, 0.0], [-0.5, 1.0, 0.0], [0.0, 0.0, 3.0]])
    return np.array([[1.0, -0.5, 0.0, -0.5],
                     [-0.5, 1.0, 0.0, -0.5],
                     [0.0, 0.0, 3.0, 0.0],
                     [-0.5, -0.5, 0.0, 1.0]])


def build_B_set(grid: Grid) -> np.ndarray:
    """(n_elem, ncomp, 8): la B di ogni elemento, al centro."""
    ncomp = 4 if grid.mode == AXISYMMETRIC else 3
    out = np.empty((grid.n_elem, ncomp, 8))
    for i in range(grid.nx):
        B = _B_centroid(grid, i)
        out[i * grid.ny:(i + 1) * grid.ny] = B
    return out


def stress_objective(prob, xi, opt, ke_set, fth_set, edofs, D, B_set,
                     sigma_allow: float, eps0_set: np.ndarray | None = None):
    """(J, dJ/dxi, sigma_vm) — norma-P della tensione, con la sua derivata.

    `eps0_set` (n_elem, ncomp) e' la deformazione termica libera di ogni
    elemento, da sottrarre per ottenere la tensione vera.
    """
    m = prob.material
    # Fondo di densita': il rilassamento usa xi^(eta-1) con eta < 1, che diverge
    # in zero. Il fondo NON e' una toppa: e' la stessa idea di E_min, cioe' che
    # il "vuoto" numerico resta un materiale debolissimo ma esistente.
    xi = np.maximum(np.asarray(xi, dtype=float), 1e-6)
    E_scale = m.E_min_ratio + xi ** opt.penal * (1.0 - m.E_min_ratio)
    K = assemble(prob.grid, ke_set, E_scale, edofs)
    f = prob.force.copy()
    thermal = fth_set is not None
    if thermal:
        t_scale = xi ** opt.thermal_penal
        np.add.at(f, edofs.ravel(), (fth_set * t_scale[:, None]).ravel())
    u = solve(K, f, prob.fixed_dofs)
    ue = u[edofs]

    # tensione elastica: sigma = D (B u - eps_termica)
    strain = np.einsum("ecd,ed->ec", B_set, ue)
    if eps0_set is not None:
        strain = strain - eps0_set * (xi ** opt.thermal_penal)[:, None]
    sigma = strain @ D.T
    V = _vm_matrix(sigma.shape[1])
    vm2 = np.einsum("ec,cd,ed->e", sigma, V, sigma)
    vm = np.sqrt(np.maximum(vm2, 1e-30))
    relax = xi ** STRESS_RELAX_ETA
    s = relax * vm / sigma_allow

    Ssum = float(np.sum(s ** PNORM))
    J = Ssum ** (1.0 / PNORM)

    # dJ/ds_e
    dJds = (Ssum ** (1.0 / PNORM - 1.0)) * s ** (PNORM - 1.0)
    # ds/dsigma  (attraverso vm)
    dsdvm = relax / sigma_allow
    dvmdsigma = (sigma @ V.T + sigma @ V) / (2.0 * vm[:, None])
    w = (dJds * dsdvm)[:, None] * dvmdsigma            # (n_elem, ncomp)

    # --- termine esplicito in xi (rilassamento + carico termico nella tensione)
    dJdxi = dJds * (STRESS_RELAX_ETA * xi ** (STRESS_RELAX_ETA - 1.0)) * vm / sigma_allow
    if eps0_set is not None:
        d_eps = -eps0_set * (opt.thermal_penal * xi ** (opt.thermal_penal - 1.0))[:, None]
        dJdxi = dJdxi + np.einsum("ec,cd,ed->e", w, D.T, d_eps)

    # --- variabile aggiunta: K lam = dJ/du
    dJdu = np.zeros_like(u)
    contrib = np.einsum("ec,cd,edk->ek", w, D.T, B_set)   # (n_elem, 8)
    np.add.at(dJdu, edofs.ravel(), contrib.ravel())
    lam = solve(K, dJdu, prob.fixed_dofs)
    lame = lam[edofs]

    dEdxi = opt.penal * xi ** (opt.penal - 1.0) * (1.0 - m.E_min_ratio)
    term = -dEdxi * np.einsum("ei,eij,ej->e", lame, ke_set, ue)
    if thermal:
        dfdxi = opt.thermal_penal * xi ** (opt.thermal_penal - 1.0)
        term = term + dfdxi * np.einsum("ei,ei->e", lame, fth_set)
    return J, dJdxi + term, vm


def thermal_free_strain(grid: Grid, alpha: float, delta_T) -> np.ndarray:
    """(n_elem, ncomp) deformazione libera alpha*dT, isotropa."""
    ncomp = 4 if grid.mode == AXISYMMETRIC else 3
    dT = np.full(grid.n_elem, float(delta_T)) if np.isscalar(delta_T) else np.asarray(delta_T)
    e = np.zeros((grid.n_elem, ncomp))
    e[:, 0] = e[:, 1] = alpha * dT
    if ncomp == 4:
        e[:, 3] = alpha * dT
    return e
