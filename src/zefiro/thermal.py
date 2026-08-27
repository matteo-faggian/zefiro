"""Carico termico a parete e risposta transitoria della parete.

Serve a rispondere a una domanda precisa: con un tempo di funzionamento
FINITO, la parete puo' semplicemente assorbire il calore (pozzo termico), o
serve davvero il film cooling? La risposta cambia il motore: senza film, tutto
il GPL va alla spinta e spariscono i fori sotto il mezzo millimetro.

Due pezzi, con statuti epistemici DIVERSI e va tenuto presente:

1. `bartz_h_gas` e' una CORRELAZIONE EMPIRICA (Bartz, "A Simple Equation for
   Rapid Estimation of Rocket Nozzle Convective Heat Transfer Coefficients",
   Jet Propulsion 27(1), 1957). Le sue costanti non si derivano: sono tarate su
   dati. L'accuratezza tipica riportata in letteratura e' dell'ordine del
   +-30 %, e peggiora fuori dal campo su cui e' stata tarata (motori a
   propellente liquido, gole molto piu' grandi di questa). Va usata come
   ordine di grandezza per decidere se vale la pena fare la CFD, MAI come
   verifica.

2. `transient_slab` e' invece FISICA ESATTA: l'equazione del calore 1D
   risolta numericamente. Non contiene costanti tarate. I suoi test la
   verificano contro due soluzioni analitiche indipendenti.

Le proprieta' del gas (mu, cp, Pr) vengono da Cantera, quindi dai dati di
trasporto del meccanismo. Le proprieta' del METALLO sono argomenti
obbligatori: sono il TODO n.J e non hanno un default.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

BARTZ_REFERENCE = (
    "Bartz, D.R., 'A Simple Equation for Rapid Estimation of Rocket Nozzle "
    "Convective Heat Transfer Coefficients', Jet Propulsion 27(1), 1957"
)
#: Accuratezza tipica riportata per la correlazione di Bartz.
BARTZ_UNCERTAINTY = 0.30


def recovery_factor(prandtl: float, turbulent: bool = True) -> float:
    """Fattore di recupero: r = Pr^(1/3) turbolento, Pr^(1/2) laminare.

    E' un risultato di strato limite, non una costante tarata.
    """
    return prandtl ** (1.0 / 3.0) if turbulent else math.sqrt(prandtl)


def adiabatic_wall_temperature(
    T_stagnation: float, mach: float, gamma: float, prandtl: float
) -> float:
    """Temperatura di parete adiabatica.

        T_aw = T_statica * (1 + r (g-1)/2 M^2)
        T_statica = T_ristagno / (1 + (g-1)/2 M^2)

    E' la temperatura che "vede" la parete, non quella di ristagno: in camera
    (M ~ 0) coincidono, in gola no.
    """
    f = 1.0 + 0.5 * (gamma - 1.0) * mach * mach
    return (T_stagnation / f) * (1.0 + recovery_factor(prandtl) * 0.5 * (gamma - 1.0) * mach**2)


def bartz_h_gas(
    D_t: float,
    p_c: float,
    c_star: float,
    mu: float,
    cp: float,
    prandtl: float,
    gamma: float,
    area_ratio_local: float,
    mach_local: float,
    T_wall: float,
    T_c: float,
    curvature_radius: float | None = None,
) -> float:
    """Coefficiente di scambio lato gas [W/(m^2 K)] secondo Bartz.

        h = (0.026/D_t^0.2) (mu^0.2 cp / Pr^0.6) (p_c/c*)^0.8
            (D_t/R_curv)^0.1 (A_t/A)^0.9 sigma

    con il fattore di correzione per la variazione di proprieta' attraverso lo
    strato limite

        sigma = [0.5 (T_w/T_c)(1 + (g-1)/2 M^2) + 0.5]^-0.68
                (1 + (g-1)/2 M^2)^-0.12

    `area_ratio_local` e' A_locale/A_gola (>= 1 ovunque).
    `curvature_radius` e' il raggio di curvatura della parete alla gola; se
    None si usa 1.5 D_t, che e' il valore su cui la correlazione e' comunemente
    applicata: quel termine entra alla potenza 0.1, quindi un fattore 2 di
    errore su di esso sposta h del 7 %.

    ATTENZIONE: correlazione empirica, +-30 % tipico. Vedi il docstring del
    modulo.
    """
    if curvature_radius is None:
        curvature_radius = 1.5 * D_t
    f = 1.0 + 0.5 * (gamma - 1.0) * mach_local**2
    sigma = (0.5 * (T_wall / T_c) * f + 0.5) ** -0.68 * f ** -0.12
    return (
        (0.026 / D_t**0.2)
        * (mu**0.2 * cp / prandtl**0.6)
        * (p_c / c_star) ** 0.8
        * (D_t / curvature_radius) ** 0.1
        * (1.0 / area_ratio_local) ** 0.9
        * sigma
    )


# --------------------------------------------------------------------------- #
# Risposta transitoria della parete: fisica esatta, nessuna costante tarata
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class SlabResponse:
    t: np.ndarray            # s
    T_hot: np.ndarray        # K, faccia calda
    T_cold: np.ndarray       # K, faccia fredda (dorso)
    q_hot: np.ndarray        # W/m^2 entrante
    energy_in: float         # J/m^2 netti entrati (entrata dal gas meno uscita dal dorso)
    energy_stored: float     # J/m^2 nella parete a fine transitorio
    penetration_depth: float # m, sqrt(alpha t_end)
    biot: float


def transient_slab(
    thickness: float,
    rho: float,
    cp: float,
    k: float,
    h_gas: float,
    T_aw: float,
    t_end: float,
    T_initial: float = 293.15,
    n_nodes: int = 81,
    n_steps: int = 4000,
    h_back: float = 0.0,
    T_back: float = 293.15,
) -> SlabResponse:
    """Parete piana 1D: convezione sulla faccia calda, dorso adiabatico o raffreddato.

    Con `h_back = 0` (default) il dorso e' ADIABATICO, che e' il caso peggiore:
    qualunque perdita reale verso l'esterno puo' solo abbassare la temperatura.
    Un progetto che sopravvive con dorso adiabatico sopravvive.

    Con `h_back > 0` il dorso scambia con un refrigerante a `T_back`: e' il caso
    della parete raffreddata ad acqua, dove il transitorio si esaurisce in
    frazioni di secondo e il regime stazionario e' il caso dimensionante. Serve
    per confrontare le due strategie sullo stesso grafico.

    Perche' 1D: lo spessore e' molto minore del raggio di camera, quindi la
    curvatura e la conduzione assiale sono correzioni del secondo ordine sul
    tempo scala di pochi secondi.

    Schema: Crank-Nicolson (implicito, incondizionatamente stabile e accurato
    al secondo ordine nel tempo), risolto con un sistema tridiagonale.
    Proprieta' costanti: e' l'approssimazione principale, e per l'acciaio fra
    300 K e 1000 K k e cp variano del 30-40 %, quindi il risultato e' un ordine
    di grandezza, non una verifica.
    """
    if thickness <= 0 or rho <= 0 or cp <= 0 or k <= 0:
        raise ValueError("Proprieta' del materiale e spessore devono essere positivi.")
    alpha = k / (rho * cp)
    dx = thickness / (n_nodes - 1)
    dt = t_end / n_steps
    r = alpha * dt / (dx * dx)

    T = np.full(n_nodes, float(T_initial))
    t_hist = np.empty(n_steps + 1)
    hot = np.empty(n_steps + 1)
    cold = np.empty(n_steps + 1)
    q_hist = np.empty(n_steps + 1)
    t_hist[0], hot[0], cold[0] = 0.0, T[0], T[-1]
    q_hist[0] = h_gas * (T_aw - T[0])

    # matrice tridiagonale di Crank-Nicolson, costante nel tempo
    a = np.full(n_nodes, -r / 2.0)          # sotto-diagonale
    b = np.full(n_nodes, 1.0 + r)           # diagonale
    c = np.full(n_nodes, -r / 2.0)          # sopra-diagonale
    # nodo 0: bilancio su mezzo volume con flusso convettivo entrante
    bi_dx = h_gas * dx / k
    b[0] = 1.0 + r + r * bi_dx
    c[0] = -r
    # nodo N-1: dorso, mezzo volume. h_back = 0 -> adiabatico.
    bo_dx = h_back * dx / k
    a[-1] = -r
    b[-1] = 1.0 + r + r * bo_dx

    energy_in = 0.0
    for n in range(n_steps):
        rhs = np.empty(n_nodes)
        rhs[1:-1] = T[1:-1] + 0.5 * r * (T[2:] - 2.0 * T[1:-1] + T[:-2])
        rhs[0] = T[0] + r * (T[1] - T[0] - bi_dx * (T[0] - T_aw)) + r * bi_dx * T_aw
        rhs[-1] = T[-1] + r * (T[-2] - T[-1]) - r * bo_dx * T[-1] + 2.0 * r * bo_dx * T_back
        T_new = _thomas(a, b, c, rhs)
        q_mid = h_gas * (T_aw - 0.5 * (T[0] + T_new[0]))
        q_out = h_back * (0.5 * (T[-1] + T_new[-1]) - T_back)
        energy_in += (q_mid - q_out) * dt
        T = T_new
        t_hist[n + 1], hot[n + 1], cold[n + 1] = (n + 1) * dt, T[0], T[-1]
        q_hist[n + 1] = h_gas * (T_aw - T[0])

    w = np.full(n_nodes, dx)
    w[0] = w[-1] = 0.5 * dx
    stored = float(np.sum(rho * cp * w * (T - T_initial)))
    return SlabResponse(
        t=t_hist, T_hot=hot, T_cold=cold, q_hot=q_hist,
        energy_in=energy_in, energy_stored=stored,
        penetration_depth=math.sqrt(alpha * t_end),
        biot=h_gas * thickness / k,
    )


def _thomas(a: np.ndarray, b: np.ndarray, c: np.ndarray, d: np.ndarray) -> np.ndarray:
    """Algoritmo di Thomas per sistemi tridiagonali."""
    n = len(d)
    cp_ = np.empty(n)
    dp = np.empty(n)
    cp_[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp_[i - 1]
        cp_[i] = c[i] / m if i < n - 1 else 0.0
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = np.empty(n)
    x[-1] = dp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp_[i] * x[i + 1]
    return x


def semi_infinite_surface_temperature(
    q: float, t: float, rho: float, cp: float, k: float, T_initial: float = 293.15
) -> float:
    """Soluzione analitica: semispazio con flusso COSTANTE imposto in superficie.

        T(0,t) = T_0 + 2 q sqrt(alpha t / pi) / k

    Serve come verifica indipendente di `transient_slab` ai tempi brevi, cioe'
    finche' la profondita' di penetrazione e' molto minore dello spessore.
    """
    alpha = k / (rho * cp)
    return T_initial + 2.0 * q * math.sqrt(alpha * t / math.pi) / k


STEFAN_BOLTZMANN = 5.670374419e-8   # W/(m^2 K^4), esatta per definizione SI


def gas_radiation_flux(
    emissivity_gas: float, T_gas: float, T_wall: float, absorptivity_wall: float = 1.0
) -> float:
    """Flusso radiativo dai gas alla parete [W/m^2].

        q = alpha_w sigma (eps_g T_g^4 - eps_w T_w^4)

    Qui in forma semplificata con parete assunta grigia e assorbente.

    PERCHE' SERVE: la correlazione di Bartz e' puramente CONVETTIVA. In camera,
    dove il flusso convettivo e' basso perche' la velocita' e' bassa, il
    contributo radiativo di CO2 e H2O a oltre 2000 K e' dello stesso ordine, e
    ignorarlo sottostima la temperatura di parete. In gola invece il convettivo
    domina di quasi un ordine di grandezza e il radiativo si puo' trascurare.

    `emissivity_gas` NON e' calcolabile qui: dipende da pressioni parziali di
    CO2 e H2O e dalla lunghezza caratteristica del volume di gas. Va preso da
    correlazioni di emissivita' dei gas (Hottel) o dalla CFD. Passare un
    intervallo e guardare la sensibilita' e' l'uso corretto di questa funzione.
    """
    return absorptivity_wall * STEFAN_BOLTZMANN * (
        emissivity_gas * T_gas**4 - absorptivity_wall * T_wall**4
    )


def mean_beam_length(volume: float, surface: float) -> float:
    """Lunghezza media del raggio: L = 3.6 V / A.

    E' la lunghezza caratteristica con cui si entra nelle correlazioni di
    emissivita' dei gas. Il coefficiente 3.6 e' un risultato geometrico medio
    su forme convesse, non un parametro tarato.
    """
    return 3.6 * volume / surface
