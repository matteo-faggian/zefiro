"""Obiettivi e vincoli.

Convenzione RIGIDA (docs/architettura.md sezione 4.7):
  * `f` si MINIMIZZA sempre. Un obiettivo da massimizzare entra come -valore, e
    l'inversione la fa questo modulo, non l'ottimizzatore: cosi' cambiare
    ottimizzatore non cambia il significato dei numeri.
  * `g <= 0` e' FATTIBILE, sempre.

I vincoli sono **adimensionalizzati**. Non e' cosmesi: pymoo e BoTorch
aggregano le violazioni, e sommare un Dp in pascal con un margine strutturale
puro darebbe alla pressione un peso arbitrario di 10^5. Normalizzare rende la
somma delle violazioni una quantita' con un significato.

I nomi qui sotto sono un REGISTRO: sono le colonne del database (store/db.py).
Aggiungerne uno e' un cambio di schema e va versionato.
"""
from __future__ import annotations

from typing import Mapping

from zefiro.schemas import (
    FEMResult,
    GeometryArtifact,
    L0Result,
    Objectives,
    OperatingPoint,
)

#: Registro degli obiettivi. Ordine stabile: e' l'ordine delle colonne.
OBJECTIVE_NAMES: tuple[str, ...] = ("neg_thrust", "neg_isp_total", "wall_volume")

#: Registro dei vincoli. `g <= 0` fattibile.
CONSTRAINT_NAMES: tuple[str, ...] = (
    "fuel_dp_stability",     # il Dp iniettore GPL deve essere >= 15 % di p_c
    "fuel_supply_pressure",  # p_c + Dp_GPL deve stare sotto la pressione di bombola
    "ox_supply_pressure",    # idem lato aria
    "air_supply",            # la portata richiesta deve essere disponibile
    "min_feature",           # quote minime fabbricabili in SLM
    "watertight",            # la geometria deve essere un solido chiuso
    "margin_yield",          # (solo L1) margine strutturale
    "margin_temp",           # (solo L1) margine termico
)

#: Frazione minima di p_c che il Dp dell'iniettore combustibile deve garantire.
#: Sotto questa soglia l'iniezione si accoppia con l'acustica di camera. E' una
#: SCELTA di progetto, dichiarata qui e non sepolta in una formula.
MIN_INJECTOR_DP_FRACTION = 0.15


def objectives_l0(
    run_id: str,
    l0: L0Result,
    op: OperatingPoint,
    p_c: float,
    geometry: GeometryArtifact | None = None,
    min_feature_size: float | None = None,
    derived: Mapping[str, float] | None = None,
) -> Objectives:
    """Obiettivi e vincoli valutabili a L0, in millisecondi.

    Serve a scartare il grosso dei candidati prima di spendere una run L1.

    Un vincolo che NON puo' essere valutato (perche' manca un dato) resta
    **fuori** dal dizionario invece di essere messo a zero. Metterlo a zero
    direbbe "soddisfatto", che e' falso: la differenza fra "verificato" e
    "non verificabile" non va nascosta.
    """
    f: dict[str, float] = {
        "neg_thrust": -l0.thrust,
        "neg_isp_total": -l0.Isp_s,
    }
    g: dict[str, float] = {}
    source: dict[str, str] = {k: "l0.cycle" for k in f}

    d = dict(derived or {})

    # --- vincoli di alimentazione ---------------------------------------- #
    dp_f = d.get("dp_inj_fuel")
    if dp_f is not None:
        # Dp >= 0.15 p_c   ->   g = (0.15 p_c - Dp) / p_c
        g["fuel_dp_stability"] = (MIN_INJECTOR_DP_FRACTION * p_c - dp_f) / p_c
        # p_c + Dp <= p_bombola   ->   g = (p_c + Dp - p_supply) / p_supply
        g["fuel_supply_pressure"] = (p_c + dp_f - op.p_fuel_supply) / op.p_fuel_supply
        source["fuel_dp_stability"] = "geometry.parameters"
        source["fuel_supply_pressure"] = "geometry.parameters"
    dp_ox = d.get("dp_inj_ox")
    if dp_ox is not None:
        g["ox_supply_pressure"] = (p_c + dp_ox - op.p_air_supply) / op.p_air_supply
        source["ox_supply_pressure"] = "geometry.parameters"

    if op.mdot_air_max is not None:
        g["air_supply"] = (l0.mdot_air - op.mdot_air_max) / op.mdot_air_max
        source["air_supply"] = "operating_point"

    # --- vincoli geometrici ----------------------------------------------- #
    if geometry is not None:
        f["wall_volume"] = geometry.volume       # proxy di massa a densita' fissa
        source["wall_volume"] = "geometry.build"
        g["watertight"] = 0.0 if geometry.is_watertight_mesh else 1.0
        source["watertight"] = "geometry.build"

    if min_feature_size is not None and d:
        quote = [d[k] for k in ("d_ox", "d_fuel", "d_film", "film_land", "t_wall")
                 if d.get(k, 0.0) > 0.0]
        if quote:
            # g = (min_feature - quota_minima) / min_feature
            g["min_feature"] = (min_feature_size - min(quote)) / min_feature_size
            source["min_feature"] = "geometry.parameters"

    return Objectives(run_id=run_id, fidelity="L0", f=f, g=g, source=source)


def objectives_l1(
    run_id: str,
    l0: L0Result,
    op: OperatingPoint,
    p_c: float,
    fem: FEMResult,
    geometry: GeometryArtifact | None = None,
    min_feature_size: float | None = None,
    derived: Mapping[str, float] | None = None,
) -> Objectives:
    """Obiettivi completi, con i margini dal FEM.

    Nota sulla convenzione: `FEMResult.margin_yield` e' un MARGINE (positivo =
    buono), mentre `g` vuole il contrario. L'inversione avviene qui, una volta
    sola, invece che in ogni punto d'uso.
    """
    base = objectives_l0(run_id, l0, op, p_c, geometry, min_feature_size, derived)
    g = dict(base.g)
    g["margin_yield"] = -fem.margin_yield
    g["margin_temp"] = -fem.margin_temp
    src = dict(base.source)
    src["margin_yield"] = "l1.fem"
    src["margin_temp"] = "l1.fem"
    return Objectives(run_id=run_id, fidelity="L1", f=base.f, g=g, source=src)
