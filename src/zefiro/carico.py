"""Dove va il calore lungo il percorso del gas, come funzione riusabile.

Era dentro `scripts/heat_map_50N.py`, cioe' dentro uno script: il sintetizzatore
ha bisogno della stessa mappa per dimensionare i circuiti di raffreddamento, e
una funzione che sta in uno script non e' riusabile senza importare lo script.

Su un aerospike a espansione esterna la risposta non e' ovvia: la gola non e'
sulla parete di camera, e' sul plug, con il labbro come sola frontiera esterna.
Se il carico fosse li', i canali in camera raffredderebbero dove non serve -
che e' esattamente l'errore gia' fatto una volta col film cooling in testa.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.optimize import brentq

from zefiro.l0.mixture import MixtureModel
from zefiro.thermal import adiabatic_wall_temperature, bartz_h_gas

#: Temperatura di parete assunta per valutare Bartz. Non e' la temperatura che
#: la parete AVRA': e' un valore di riferimento con cui si calcola h, e h
#: dipende poco da questa scelta (entra come (T_w/T_0)^-0.68 nel fattore di
#: correzione di Bartz). 800 K e' un compromesso fra la camera raffreddata e
#: il labbro, dove la parete corre piu' calda.
T_PARETE_RIFERIMENTO = 800.0


def mappa_di_calore(op, params, l0, p_c, n=140, T_wall=T_PARETE_RIFERIMENTO):
    """q(x) lungo la parete di camera e lungo il plug.

    L'area di passaggio locale e' quella ANULARE fra parete e corpo centrale,
    non pi r^2: e' la differenza fra una camera anulare e un tubo, e sbagliarla
    sposta il rapporto di aree quindi h.
    """

    d = params.derived
    model = MixtureModel.from_fuel(op.fuel)
    g = model.gas
    g.TPX = l0.T_ad, p_c, l0.X_eq
    g.equilibrate("TP")
    mu, cp = g.viscosity, g.cp_mass
    Pr = mu * cp / g.thermal_conductivity
    gamma, A_t = l0.gamma_c, l0.A_t
    D_t = 2.0 * math.sqrt(A_t / math.pi)

    R_c, R_lip, L_c, L_conv = d["R_c"], d["R_lip"], d["L_c"], d["L_conv"]
    r_cb = d["r_centerbody"]
    x_lip = L_c + L_conv
    cx = np.asarray(params.plug_contour_x)
    cr = np.asarray(params.plug_contour_r)

    stazioni = []                       # (x, r_parete, A_locale, zona)
    for xx in np.linspace(0.0, L_c, n // 3):
        stazioni.append((xx, R_c, math.pi * (R_c**2 - r_cb**2), "camera"))
    for xx in np.linspace(L_c, x_lip, n // 3):
        t = (xx - L_c) / max(L_conv, 1e-12)
        r_out = R_c + t * (R_lip - R_c)
        stazioni.append((xx, r_out, math.pi * (r_out**2 - r_cb**2), "convergente"))
    for i in range(len(cx)):
        # sul plug l'area di passaggio e' quella normale al flusso: si usa il
        # rapporto d'aree isentropico gia' coerente col contorno di Angelino
        r = cr[i]
        A = math.pi * (R_lip**2 - r**2) if r < R_lip else A_t
        stazioni.append((x_lip + cx[i], r, max(A, A_t), "plug"))

    righe = []
    for xx, rw, A, zona in stazioni:
        eps = max(A / A_t, 1.0)
        M = mach_da_area(eps, gamma, supersonico=(zona == "plug"))
        T_aw = adiabatic_wall_temperature(l0.T_ad, M, gamma, Pr)
        h = bartz_h_gas(D_t, p_c, l0.c_star, mu, cp, Pr, gamma, eps, M,
                        T_wall, l0.T_ad)
        righe.append((xx, rw, zona, eps, M, T_aw, h, h * (T_aw - T_wall)))
    return righe


def mach_da_area(eps, gamma, supersonico):
    """Mach dal rapporto d'aree isentropico. Due radici: si sceglie il ramo."""
    if eps <= 1.0 + 1e-12:
        return 1.0
    g = gamma
    def f(M):
        return (1.0 / M) * ((2.0 / (g + 1.0)) * (1.0 + 0.5 * (g - 1.0) * M * M)) ** (
            (g + 1.0) / (2.0 * (g - 1.0))) - eps
    return brentq(f, 1.0000001, 12.0) if supersonico else brentq(f, 1e-4, 0.9999999)




def potenza_per_zona(righe, params) -> dict[str, float]:
    """Integra q sulla superficie bagnata, zona per zona [W].

    L'integrale e' su un tronco di cono fra due stazioni consecutive:
    dA = 2 pi r_medio * lunghezza dello spigolo. Usare 2 pi r dx invece della
    lunghezza vera dello spigolo sottostima l'area del convergente del 25 %.
    """
    tot: dict[str, float] = {}
    for a, b in zip(righe[:-1], righe[1:]):
        if a[2] != b[2]:
            continue
        dx, dr = b[0] - a[0], b[1] - a[1]
        ds = math.hypot(dx, dr)
        area = math.pi * (a[1] + b[1]) * ds
        tot[a[2]] = tot.get(a[2], 0.0) + 0.5 * (a[7] + b[7]) * area
    return tot


def q_massimo_per_zona(righe) -> dict[str, float]:
    m: dict[str, float] = {}
    for r in righe:
        m[r[2]] = max(m.get(r[2], 0.0), r[7])
    return m
