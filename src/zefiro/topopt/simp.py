"""Ottimizzazione topologica per densita' (SIMP) con criterio di ottimalita'.

La topologia e' il RISULTATO, non un ingresso: si parte da un dominio pieno di
materiale a densita' uniforme e l'ottimizzatore decide dove serve, aprendo buchi
e facendo comparire nervature dove il carico le richiede. E' cio' che distingue
il design generativo da quello parametrico, dove la forma la decide chi scrive
il codice e l'ottimizzatore muove solo delle quote.

Il punto delicato: **i carichi termici dipendono dal progetto**. Con la
pressione, il carico e' fisso e la compliance e' autoaggiunta: dc/dx =
-u' (dK/dx) u, un solo termine. Con la dilatazione impedita il carico stesso e'
proporzionale alla materia presente, e la derivata corretta e'

    dc/dx = 2 u' (df/dx) - u' (dK/dx) u

Sbagliare quel fattore 2 (o dimenticare il primo termine del tutto) produce un
ottimizzatore che converge tranquillamente verso qualcosa che non e' un ottimo,
senza dare alcun segnale. Il test `test_topopt.py` confronta la sensibilita' con
differenze finite proprio per questo.
"""
from __future__ import annotations

import json
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
import scipy.sparse as sp

from zefiro.topopt import filters as flt
from zefiro.topopt.fem2d import Grid, assemble, constitutive, solve


@dataclass
class Material:
    """Proprieta' del materiale. Nessun default: sono il TODO n.J."""
    E: float
    nu: float
    alpha: float = 0.0          # 1/K, dilatazione termica
    E_min_ratio: float = 1e-9   # rigidezza residua del vuoto, per non singolare


@dataclass
class Problem:
    """Il problema strutturale: dominio, carichi, vincoli, zone congelate."""
    grid: Grid
    material: Material
    fixed_dofs: np.ndarray
    force: np.ndarray                                  # carico indipendente dal progetto
    delta_T: np.ndarray | float = 0.0                  # per elemento, K
    void: np.ndarray | None = None                     # elementi forzati a vuoto
    solid: np.ndarray | None = None                    # elementi forzati a pieno
    name: str = "problema"


@dataclass
class Options:
    volume_fraction: float = 0.4
    penal: float = 3.0
    #: Esponente di penalizzazione del CARICO termico. Tenuto a 1 e non uguale a
    #: `penal`: con lo stesso esponente il carico svanisce alla stessa velocita'
    #: della rigidezza e le densita' intermedie diventano artificialmente
    #: convenienti, perche' portano carico quasi nullo a costo di volume ridotto.
    thermal_penal: float = 1.0
    rmin: float = 2.0                  # unita' FISICHE, non elementi
    move: float = 0.2
    max_iter: int = 200
    tol: float = 1.0e-3
    beta: float = 1.0                  # nitidezza iniziale di Heaviside
    beta_max: float = 16.0
    beta_every: int = 40               # raddoppia beta ogni N iterazioni
    printable: bool = True             # applica il filtro di stampabilita' SLM
    checkpoint_every: int = 25
    outer_per_beta: int = 6            # cicli del Lagrangiano per livello di beta
    inner_maxiter: int = 30            # iterazioni L-BFGS-B per ciclo
    #: "compliance" oppure "stress". Con carico TERMICO la compliance non e' un
    #: surrogato valido di resistenza (aggiungere materia aggiunge anche carico):
    #: in quel caso serve "stress". Vedi il docstring di topopt/stress.py.
    objective: str = "compliance"
    sigma_allow: float | None = None   # Pa, obbligatorio se objective = "stress"
    #: "auto" sceglie OC per la compliance e il Lagrangiano aumentato per la
    #: tensione. Forzare OC sulla tensione fa divergere l'ottimizzazione: la
    #: sensibilita' cambia segno e l'aggiornamento OC non e' piu' valido.
    solver: str = "auto"


@dataclass
class Result:
    x: np.ndarray                  # progetto (blueprint)
    xi: np.ndarray                 # densita' fisica, cioe' cio' che si stampa
    compliance: float              # valore dell'obiettivo scelto
    sigma_vm_max: float = float("nan")
    history: list[dict] = field(default_factory=list)
    iterations: int = 0
    seconds: float = 0.0
    overhang_violation: float = 0.0


def _physical(x: np.ndarray, H: sp.csr_matrix, opt: Options, beta: float,
              nx: int, ny: int):
    """x -> filtro di densita' -> Heaviside -> stampabilita' -> densita' fisica."""
    xf = H @ x
    xh = flt.heaviside(xf, beta)
    xi = flt.printable(xh, nx, ny) if opt.printable else xh
    return xf, xh, xi


def _chain(dfdxi: np.ndarray, xf: np.ndarray, xh: np.ndarray,
           H: sp.csr_matrix, opt: Options, beta: float, nx: int, ny: int):
    """Riporta una derivata da xi fino a x, attraversando i tre filtri."""
    d = flt.printable_chain(xh, nx, ny, dfdxi) if opt.printable else dfdxi
    d = d * flt.heaviside_derivative(xf, beta)
    return H.T @ d


def compliance_and_sensitivity(prob: Problem, xi: np.ndarray, opt: Options,
                               ke_set: np.ndarray, fth_set: np.ndarray,
                               edofs: np.ndarray):
    """(compliance, dc/dxi, u). La derivata include il termine di carico termico."""
    m = prob.material
    E_scale = m.E_min_ratio + xi ** opt.penal * (1.0 - m.E_min_ratio)
    K = assemble(prob.grid, ke_set, E_scale, edofs)

    f = prob.force.copy()
    thermal = fth_set is not None
    if thermal:
        t_scale = xi ** opt.thermal_penal
        fe = fth_set * t_scale[:, None]
        np.add.at(f, edofs.ravel(), fe.ravel())

    u = solve(K, f, prob.fixed_dofs)
    ue = u[edofs]                                        # (n_elem, 8)

    # termine di rigidezza: -u' dK/dxi u
    uKu = np.einsum("ei,eij,ej->e", ue, ke_set, ue)
    dEdxi = opt.penal * xi ** (opt.penal - 1.0) * (1.0 - m.E_min_ratio)
    dc = -dEdxi * uKu

    c = float(f @ u)
    if thermal:
        # + 2 u' df/dxi  (carico dipendente dal progetto: NON autoaggiunto)
        dfdxi = opt.thermal_penal * xi ** (opt.thermal_penal - 1.0)
        dc = dc + 2.0 * dfdxi * np.einsum("ei,ei->e", ue, fth_set)
    return c, dc, u


def _oc_update(x: np.ndarray, dc: np.ndarray, dv: np.ndarray, opt: Options,
               target_volume: float, void: np.ndarray | None,
               solid: np.ndarray | None,
               physical: Callable[[np.ndarray], np.ndarray]) -> np.ndarray:
    """Aggiornamento per criterio di ottimalita', con bisezione sul moltiplicatore.

    Il vincolo si valuta sulla densita' FISICA, non sul progetto grezzo. Non e'
    un dettaglio: i filtri (densita', Heaviside, stampabilita') cambiano il
    volume, e imporre il vincolo prima di attraversarli significa consegnare un
    pezzo con una frazione di materiale diversa da quella richiesta. Costa una
    passata di filtri per iterazione di bisezione — trascurabile accanto alla
    soluzione del sistema lineare.
    """
    lo, hi = 1e-12, 1e12
    xn = x
    dcn = np.minimum(dc, -1e-14)          # deve essere negativa: piu' materia, meno obiettivo
    for _ in range(60):
        lmid = 0.5 * (lo + hi)
        step = np.sqrt(-dcn / (lmid * np.maximum(dv, 1e-14)))
        xn = np.clip(np.clip(x * step, x - opt.move, x + opt.move), 0.0, 1.0)
        if void is not None:
            xn[void] = 0.0
        if solid is not None:
            xn[solid] = 1.0
        if physical(xn).sum() > target_volume:
            lo = lmid
        else:
            hi = lmid
        if (hi - lo) / (hi + lo) < 1e-9:
            break
    return xn


def optimize(prob: Problem, opt: Options | None = None,
             x0: np.ndarray | None = None,
             checkpoint: Path | None = None,
             on_iteration: Callable[[dict], None] | None = None) -> Result:
    """Esegue l'ottimizzazione. Pensata per girare per ore, quindi salva."""
    opt = opt or Options()
    g = prob.grid
    nx, ny, n = g.nx, g.ny, g.n_elem
    D = constitutive(prob.material.E, prob.material.nu, g.mode)
    ke_set = g.element_stiffness_set(D)
    edofs = g.element_dofs()
    fth_set = None
    eps0_set = None
    if prob.material.alpha != 0.0 and np.any(np.asarray(prob.delta_T) != 0.0):
        fth_set = g.element_thermal_set(D, prob.material.alpha, prob.delta_T)

    if opt.objective == "stress":
        from zefiro.topopt.stress import build_B_set, stress_objective, thermal_free_strain
        if opt.sigma_allow is None:
            raise ValueError("objective = 'stress' richiede sigma_allow (Pa)")
        B_set = build_B_set(g)
        if fth_set is not None:
            eps0_set = thermal_free_strain(g, prob.material.alpha, prob.delta_T)

        def objective(xi):
            J, dJ, vm = stress_objective(prob, xi, opt, ke_set, fth_set, edofs,
                                         D, B_set, opt.sigma_allow, eps0_set)
            return J, dJ, float(vm.max())
    elif opt.objective == "compliance":
        def objective(xi):
            c, dc, _ = compliance_and_sensitivity(prob, xi, opt, ke_set, fth_set, edofs)
            return c, dc, float("nan")
    else:
        raise ValueError(f"obiettivo non riconosciuto: {opt.objective!r}")

    H = flt.density_filter_matrix(nx, ny, opt.rmin, g.dx, g.dy)
    target = opt.volume_fraction * n

    x = np.full(n, opt.volume_fraction) if x0 is None else x0.copy()
    if prob.void is not None:
        x[prob.void] = 0.0
    if prob.solid is not None:
        x[prob.solid] = 1.0

    beta = opt.beta
    history: list[dict] = []
    t0 = time.perf_counter()
    change = 1.0
    it = 0

    solver = opt.solver
    if solver == "auto":
        solver = "nlp" if opt.objective == "stress" else "oc"

    if solver == "nlp":
        from zefiro.topopt.nlp import NlpOptions, solve_augmented_lagrangian

        def _obj(v):
            xf, xh, xi = _physical(v, H, opt, beta, nx, ny)
            J, dJdxi, vmax = objective(xi)
            return J, _chain(dJdxi, xf, xh, H, opt, beta, nx, ny)

        def _vol(v):
            xf, xh, xi = _physical(v, H, opt, beta, nx, ny)
            return float(xi.sum()), _chain(np.ones(n), xf, xh, H, opt, beta, nx, ny)

        n_beta = max(1, int(np.ceil(np.log2(max(opt.beta_max / max(beta, 0.5), 1.0)))) + 1)
        per = max(2, opt.outer_per_beta)
        for _ in range(n_beta):
            x, hist = solve_augmented_lagrangian(
                x, _obj, _vol, target,
                NlpOptions(outer=per, inner_maxiter=opt.inner_maxiter),
                freeze_void=prob.void, freeze_solid=prob.solid,
                on_outer=on_iteration)
            for r_ in hist:
                r_["beta"] = beta
            history.extend(hist)
            if checkpoint:
                save_checkpoint(checkpoint, x, history, opt, prob.name)
            if beta >= opt.beta_max:
                break
            beta = min(max(beta, 0.5) * 2.0, opt.beta_max)

        xf, xh, xi = _physical(x, H, opt, beta, nx, ny)
        c, _, vmax = objective(xi)
        return Result(x=x, xi=xi, compliance=c, sigma_vm_max=vmax, history=history,
                      iterations=len(history), seconds=time.perf_counter() - t0,
                      overhang_violation=flt.max_overhang_violation(xi, nx, ny))

    for it in range(1, opt.max_iter + 1):
        xf, xh, xi = _physical(x, H, opt, beta, nx, ny)
        c, dcdxi, vmax = objective(xi)
        dc = _chain(dcdxi, xf, xh, H, opt, beta, nx, ny)
        dv = _chain(np.ones(n), xf, xh, H, opt, beta, nx, ny)

        if not np.all(np.isfinite(dc)):
            raise FloatingPointError(
                f"sensibilita' non finita all'iterazione {it}: {int(np.sum(~np.isfinite(dc)))} "
                "elementi. Meglio fermarsi che proseguire producendo una struttura "
                "ottimizzata su numeri privi di senso."
            )
        xnew = _oc_update(x, dc, dv, opt, target, prob.void, prob.solid,
                          lambda xx: _physical(xx, H, opt, beta, nx, ny)[2])
        change = float(np.max(np.abs(xnew - x)))
        x = xnew

        rec = {"iter": it, "objective": c, "sigma_vm_max": vmax, "change": change,
               "beta": beta, "vol": float(xi.mean()),
               "seconds": time.perf_counter() - t0}
        history.append(rec)
        if on_iteration:
            on_iteration(rec)
        if checkpoint and it % opt.checkpoint_every == 0:
            save_checkpoint(checkpoint, x, history, opt, prob.name)

        if it % opt.beta_every == 0 and beta < opt.beta_max:
            beta = min(beta * 2.0, opt.beta_max)
            change = 1.0                    # la continuazione riapre il problema
        elif change < opt.tol and beta >= opt.beta_max:
            break

    xf, xh, xi = _physical(x, H, opt, beta, nx, ny)
    c, _, vmax = objective(xi)
    res = Result(x=x, xi=xi, compliance=c, sigma_vm_max=vmax, history=history, iterations=it,
                 seconds=time.perf_counter() - t0,
                 overhang_violation=flt.max_overhang_violation(xi, nx, ny))
    if checkpoint:
        save_checkpoint(checkpoint, x, history, opt, prob.name)
    return res


def save_checkpoint(path: Path, x: np.ndarray, history: list[dict],
                    opt: Options, name: str) -> None:
    """Una run di ore deve poter essere ripresa, e guardata mentre gira."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, x=x)
    path.with_suffix(".json").write_text(json.dumps({
        "name": name, "options": opt.__dict__, "history": history,
    }, indent=1), encoding="utf-8")


def load_checkpoint(path: Path) -> tuple[np.ndarray, list[dict]]:
    path = Path(path)
    x = np.load(path)["x"]
    meta = json.loads(path.with_suffix(".json").read_text(encoding="utf-8"))
    return x, meta["history"]
