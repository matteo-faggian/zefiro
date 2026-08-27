"""Ottimizzatore generale per obiettivi non monotoni (tensione).

Perche' non basta il criterio di ottimalita'. L'aggiornamento OC classico,

    x_nuovo = x * ( -dJ/dx / (lambda * dV/dx) )^(1/2)

e' DERIVATO dalle condizioni di stazionarieta' della compliance con un solo
vincolo di risorsa, e presuppone che dJ/dx sia negativa ovunque: piu' materiale,
meno obiettivo. Per la compliance con carico fisso e' vero e OC e' imbattibile
per semplicita' e velocita'.

Per la TENSIONE e' falso. Aggiungere materiale in un punto irrigidisce la
struttura e puo' spostare il flusso delle forze in modo da ALZARE la tensione
altrove: la sensibilita' cambia segno da elemento a elemento. Troncare i valori
positivi, come fa OC quando lo si forza, non e' un'approssimazione: e'
un'euristica che diverge — provato, la tensione di picco e' salita di quattro
ordini di grandezza.

Qui si usa quindi un **Lagrangiano aumentato** sul vincolo di volume, con
L-BFGS-B all'interno:

    L(x, lam, mu) = J(x) + lam * g(x) + (mu/2) * g(x)^2 ,   g = V(x)/V* - 1

L-BFGS-B gestisce nativamente i limiti 0 <= x <= 1 e scala a centinaia di
migliaia di variabili usando i gradienti esatti che gia' abbiamo. Il ciclo
esterno aggiorna il moltiplicatore (lam += mu*g) e alza la penalita': e' lo
schema standard, senza costanti inventate.
"""
from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Callable

import numpy as np
from scipy.optimize import minimize


@dataclass
class NlpOptions:
    outer: int = 12                 # cicli del Lagrangiano aumentato
    inner_maxiter: int = 30         # iterazioni L-BFGS-B per ciclo
    mu0: float = 10.0               # penalita' iniziale
    mu_growth: float = 1.6
    mu_max: float = 5.0e3
    tol_volume: float = 2.0e-3      # |V/V* - 1| accettabile
    verbose: bool = False


def solve_augmented_lagrangian(
    x0: np.ndarray,
    objective: Callable[[np.ndarray], tuple[float, np.ndarray]],
    volume: Callable[[np.ndarray], tuple[float, np.ndarray]],
    target_volume: float,
    opt: NlpOptions | None = None,
    freeze_void: np.ndarray | None = None,
    freeze_solid: np.ndarray | None = None,
    on_outer: Callable[[dict], None] | None = None,
):
    """Minimizza `objective` con V(x) = target, x in [0, 1].

    `objective` e `volume` restituiscono (valore, gradiente) rispetto a x, cioe'
    gia' riportati attraverso i filtri.
    """
    opt = opt or NlpOptions()
    x = np.clip(x0.astype(float), 0.0, 1.0)
    lo = np.zeros_like(x)
    hi = np.ones_like(x)
    if freeze_void is not None:
        hi[freeze_void] = 0.0
    if freeze_solid is not None:
        lo[freeze_solid] = 1.0
    x = np.clip(x, lo, hi)

    lam = 0.0
    mu = opt.mu0
    history: list[dict] = []
    t0 = time.perf_counter()

    # normalizzazione dell'obiettivo: rende la penalita' confrontabile con J
    J0 = max(abs(objective(x)[0]), 1e-30)

    for k in range(opt.outer):
        def L(v: np.ndarray):
            J, dJ = objective(v)
            V, dV = volume(v)
            g = V / target_volume - 1.0
            dg = dV / target_volume
            val = J / J0 + lam * g + 0.5 * mu * g * g
            grad = dJ / J0 + (lam + mu * g) * dg
            return float(val), grad

        res = minimize(L, x, jac=True, method="L-BFGS-B",
                       bounds=np.column_stack([lo, hi]),
                       options={"maxiter": opt.inner_maxiter, "maxcor": 15})
        x = np.clip(res.x, lo, hi)

        J, _ = objective(x)
        V, _ = volume(x)
        g = V / target_volume - 1.0
        lam += mu * g
        mu = min(mu * opt.mu_growth, opt.mu_max)

        rec = {"outer": k + 1, "objective": float(J), "volume_error": float(g),
               "mu": mu, "lambda": lam, "inner_iter": int(res.nit),
               "seconds": time.perf_counter() - t0}
        history.append(rec)
        if on_outer:
            on_outer(rec)
        if abs(g) < opt.tol_volume and k >= 3:
            break

    return x, history
