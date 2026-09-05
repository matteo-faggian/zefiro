"""Stime strutturali analitiche di primo ordine.

Non sostituiscono il FEM: servono a (a) sapere PRIMA di lanciare CalculiX quale
carico domina, e (b) avere un caso a soluzione nota contro cui verificare il
FEM stesso. Un FEM che non riproduce Lame' su un cilindro in pressione ha un
problema di modello, non di mesh, e conviene scoprirlo su questi due casi
invece che sulla geometria vera.

Tutte le formule qui dentro sono derivabili in mezza pagina. Le proprieta' del
materiale sono argomenti obbligatori (TODO n.J): non esistono default.
"""
from __future__ import annotations


def thin_wall_hoop_stress(p: float, radius: float, thickness: float) -> float:
    """Tensione circonferenziale in un cilindro a parete sottile: sigma = p r / t.

    Si ricava dall'equilibrio di mezzo cilindro: la forza di pressione sulla
    proiezione diametrale, `p * 2r * L`, e' equilibrata dalle due sezioni di
    parete, `2 * sigma * t * L`. Vale per t/r sotto circa 0.1.
    """
    if thickness <= 0.0:
        raise ValueError("Lo spessore deve essere positivo.")
    return p * radius / thickness


def lame_hoop_stress_inner(p: float, r_inner: float, r_outer: float) -> float:
    """Lame': tensione circonferenziale alla parete INTERNA, parete spessa.

        sigma_theta(r_i) = p (r_o^2 + r_i^2) / (r_o^2 - r_i^2)

    E' la soluzione esatta, senza l'ipotesi di parete sottile. Serve come caso
    di verifica del FEM (criterio di chiusura della fase 4).
    """
    if r_outer <= r_inner:
        raise ValueError("r_outer deve essere maggiore di r_inner.")
    return p * (r_outer**2 + r_inner**2) / (r_outer**2 - r_inner**2)


def fully_constrained_thermal_stress(
    youngs_modulus: float, alpha: float, delta_T: float, poisson: float
) -> float:
    """Tensione termica in una parete completamente impedita di dilatarsi:

        sigma = E alpha dT / (1 - nu)

    E' un LIMITE SUPERIORE: un vincolo reale e' sempre parziale, quindi la
    tensione vera e' minore. Va usata per capire se il termico domina sul
    meccanico, non come verifica.

    Il `(1 - nu)` viene dallo stato di deformazione piana biassiale: la parete
    e' impedita in DUE direzioni nel piano, non una.
    """
    return youngs_modulus * alpha * delta_T / (1.0 - poisson)


def biot_number(h: float, length: float, k: float) -> float:
    """Bi = h L / k. Sotto ~0.1 la parete si puo' trattare a temperatura uniforme;
    sopra, il gradiente nello spessore e' il protagonista e serve risolverlo."""
    return h * length / k
