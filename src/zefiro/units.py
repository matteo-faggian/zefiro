"""Unita' e costanti.

Regola del progetto: **dentro il codice tutto e' SI stretto**.
Questi fattori esistono solo per il confine I/O (lettura YAML, stampa report).
Se trovi una conversione altrove, e' un bug.
"""
from __future__ import annotations

# --- fattori di conversione: moltiplicare per passare A SI ---
BAR = 1.0e5          # bar   -> Pa
MM = 1.0e-3          # mm    -> m
DEG = 3.141592653589793 / 180.0   # deg -> rad
G_PER_S = 1.0e-3     # g/s   -> kg/s

# --- costanti fisiche ---
R_UNIVERSAL = 8.31446261815324   # J/(mol K)  (CODATA, esatta per definizione SI)
G0 = 9.80665                     # m/s^2      (definizione ISO 80000, usata per I_sp)

# --- composizione dell'aria secca (frazioni MOLARI) ---
# Fonte: composizione standard dell'atmosfera secca. Usata coerentemente ovunque:
# cambiare qui cambia AFR e T_ad, e questo e' voluto (un solo posto).
AIR_MOLE_FRACTIONS: dict[str, float] = {
    "O2": 0.20946,
    "N2": 0.78084,
    "AR": 0.00934,
    "CO2": 0.00036,
}


def air_composition_string(species_names: list[str]) -> str:
    """Stringa di composizione dell'aria filtrata sulle specie disponibili.

    Cantera vuole ``"O2:0.2, N2:0.78"``. Se il meccanismo non ha AR o CO2 li
    ometto e **rinormalizzo**, invece di fallire o di far finta che ci siano.
    La rinormalizzazione e' esplicita perche' cambia (di poco) AFR e T_ad.
    """
    upper = {s.upper(): s for s in species_names}
    kept = {upper[k]: v for k, v in AIR_MOLE_FRACTIONS.items() if k in upper}
    if "O2" not in {k.upper() for k in kept}:
        raise ValueError("Il meccanismo non contiene O2: non e' utilizzabile per l'aria.")
    tot = sum(kept.values())
    return ", ".join(f"{k}:{v / tot!r}" for k, v in kept.items())
