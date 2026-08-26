"""Contorno del plug di un aerospike a espansione esterna (metodo di Angelino).

La derivazione completa e' in docs/architettura.md sezione 4.3. In sintesi:
il ventaglio di espansione e' centrato sul labbro (raggio R_lip), le
caratteristiche sono rette, e il contorno del plug e' la linea di corrente che
chiude il campo. Imponendo la conservazione della massa attraverso ciascuna
caratteristica si ottiene, in forma chiusa,

    r_p(M)/R_lip = sqrt( 1 - M * (A/A*)(M) * sin(alpha) / eps )
    alpha(M)     = [nu(M_e) - nu(M)] + asin(1/M)
    l            = (R_lip - r_p) / sin(alpha)
    x_p          = l * cos(alpha)

Nessun dato tabellato: entrano solo gamma e il Mach di progetto.

Origine del sistema di riferimento: il LABBRO, in (x=0, r=R_lip). Il punto di
gola del plug ha x < 0 (la gola e' inclinata di nu(M_e) rispetto al piano
radiale, quindi cade a monte del labbro): e' corretto, non un errore di segno.
"""
from __future__ import annotations

import math
from dataclasses import dataclass


def prandtl_meyer(M: float, gamma: float) -> float:
    """Angolo di Prandtl-Meyer nu(M) [rad]. nu(1) = 0."""
    if M < 1.0:
        raise ValueError(f"Prandtl-Meyer definito solo per M >= 1, ricevuto {M!r}")
    if M == 1.0:
        return 0.0
    b = math.sqrt((gamma + 1.0) / (gamma - 1.0))
    m2 = M * M - 1.0
    return b * math.atan(math.sqrt(m2) / b) - math.atan(math.sqrt(m2))


def area_ratio(M: float, gamma: float) -> float:
    """Rapporto isentropico A/A*."""
    if M <= 0.0:
        raise ValueError(f"M deve essere positivo, ricevuto {M!r}")
    g = gamma
    return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** (
        (g + 1.0) / (2.0 * (g - 1.0))
    )


def mach_from_pressure_ratio(pr: float, gamma: float) -> float:
    """M da p0/p isentropico."""
    if pr <= 1.0:
        raise ValueError(f"p0/p deve essere > 1, ricevuto {pr!r}")
    return math.sqrt(2.0 / (gamma - 1.0) * (pr ** ((gamma - 1.0) / gamma) - 1.0))


def mach_from_area_ratio(eps: float, gamma: float, tol: float = 1e-12) -> float:
    """Ramo SUPERSONICO di A/A* = eps, per bisezione (monotono per M > 1)."""
    if eps < 1.0:
        raise ValueError(f"A/A* deve essere >= 1, ricevuto {eps!r}")
    lo, hi = 1.0, 2.0
    while area_ratio(hi, gamma) < eps:
        hi *= 2.0
        if hi > 1.0e4:
            raise ValueError("Nessuna soluzione supersonica: eps troppo grande.")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if area_ratio(mid, gamma) < eps:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


@dataclass(frozen=True)
class PlugContour:
    x: tuple[float, ...]        # m, labbro in x = 0
    r: tuple[float, ...]        # m
    M: tuple[float, ...]
    R_lip: float                # m
    epsilon: float              # A_e/A_t coerente con gamma e M_e
    M_e: float
    nu_e: float                 # rad
    r_throat: float             # m, raggio del plug alla gola
    x_throat: float             # m (negativo: a monte del labbro)
    x_tip_full: float           # m, apice del plug non troncato
    truncation: float           # frazione di lunghezza mantenuta


def plug_contour(
    gamma: float,
    M_e: float,
    R_lip: float,
    n_points: int = 160,
    truncation: float = 1.0,
) -> PlugContour:
    """Genera il contorno del plug.

    `epsilon` NON e' un ingresso: e' A/A*(M_e) con questo gamma. Imporre un eps
    ricavato altrove (per esempio dall'espansione in equilibrio spostato di L0)
    renderebbe il contorno incoerente con se' stesso e il plug non chiuderebbe
    esattamente sull'asse. La differenza fra i due eps e' un indicatore
    diagnostico utile, ed e' esposta da geometry.parameters.derive().
    """
    if not (0.0 < truncation <= 1.0):
        raise ValueError(f"truncation deve stare in (0, 1], vale {truncation!r}")
    if M_e <= 1.0:
        raise ValueError(f"M_e deve essere > 1, vale {M_e!r}")

    eps = area_ratio(M_e, gamma)
    nu_e = prandtl_meyer(M_e, gamma)

    xs: list[float] = []
    rs: list[float] = []
    ms: list[float] = []
    for i in range(n_points):
        # distribuzione uniforme in nu: concentra i punti dove il flusso gira
        # di piu', cioe' vicino alla gola, che e' dove il contorno curva.
        nu = nu_e * i / (n_points - 1)
        if i == 0:
            M = 1.0
        elif i == n_points - 1:
            M = M_e
        else:
            M = _mach_from_nu(nu, gamma)
        mu = math.asin(min(1.0, 1.0 / M))
        alpha = (nu_e - nu) + mu
        sa = math.sin(alpha)
        if i == n_points - 1:
            # Punto di apice. L'algebra da' esattamente r_p = 0 (vedi il test
            # test_il_plug_chiude_esattamente_sullasse), ma l'espressione
            # generale e' li' una differenza fra numeri quasi uguali e perde
            # ~7 cifre. Si usa il valore analitico, non la sua versione mal
            # condizionata: il profilo CAD deve chiudere sull'asse in modo
            # esatto, altrimenti la rivoluzione genera una faccia parassita.
            r_p = 0.0
        else:
            arg = 1.0 - M * area_ratio(M, gamma) * sa / eps
            r_p = R_lip * math.sqrt(max(arg, 0.0))
        ell = (R_lip - r_p) / sa
        xs.append(ell * math.cos(alpha))
        rs.append(r_p)
        ms.append(M)

    x_tip_full = xs[-1]
    if truncation < 1.0:
        x_cut = xs[0] + truncation * (x_tip_full - xs[0])
        keep = [i for i, x in enumerate(xs) if x <= x_cut]
        if len(keep) < 3:
            raise ValueError("truncation troppo aggressiva: contorno degenere.")
        xs, rs, ms = [xs[i] for i in keep], [rs[i] for i in keep], [ms[i] for i in keep]

    return PlugContour(
        x=tuple(xs), r=tuple(rs), M=tuple(ms),
        R_lip=R_lip, epsilon=eps, M_e=M_e, nu_e=nu_e,
        r_throat=rs[0], x_throat=xs[0], x_tip_full=x_tip_full, truncation=truncation,
    )


def _mach_from_nu(nu: float, gamma: float, tol: float = 1e-12) -> float:
    """Inverte nu(M) per bisezione (nu e' monotona crescente in M >= 1)."""
    if nu <= 0.0:
        return 1.0
    lo, hi = 1.0, 2.0
    while prandtl_meyer(hi, gamma) < nu:
        hi *= 2.0
        if hi > 1.0e4:
            raise ValueError("nu fuori dal dominio raggiungibile.")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if prandtl_meyer(mid, gamma) < nu:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)
