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
#:
#: Sono le metriche CALCOLATE e archiviate a ogni run. Quali di queste vengano
#: effettivamente messe in un fronte di Pareto e' una scelta separata, che vive
#: nel driver (`opt.driver.OBJECTIVE_SETS`) e non qui: archiviare e' gratis,
#: ottimizzare no.
OBJECTIVE_NAMES: tuple[str, ...] = (
    "neg_thrust",       # -F            [N]      spinta
    "neg_isp_total",    # -Isp          [s]      impulso specifico sul totale
    "neg_isp_fuel",     # -Isp_fuel     [s]      spinta per kg di GPL
    "q_throat",         #  q            [W/m^2]  severita' termica in gola
    "wall_volume",      #  V            [m^3]    proxy di massa (Pappo-Guldino)
    "neg_damkohler",    # -Da           [-]      completezza della combustione
)

#: Registro dei vincoli. `g <= 0` fattibile.
CONSTRAINT_NAMES: tuple[str, ...] = (
    "fuel_dp_stability",     # il Dp iniettore GPL deve essere >= 15 % di p_c
    "fuel_supply_pressure",  # p_c + Dp_GPL deve stare sotto la pressione di bombola
    "ox_dp_stability",       # idem, lato ARIA: e' il 94 % della portata
    "ox_supply_pressure",    # p_c + Dp_aria sotto la pressione di serbatoio
    "air_supply",            # la portata richiesta deve essere disponibile
    "min_feature",           # quote minime fabbricabili in SLM
    "watertight",            # la geometria deve essere un solido chiuso
    "throat_heat_flux",      # q di gola entro cio' che il raffreddamento estrae
    "combustion_residence",  # Damkohler minimo: la camera deve bruciare
    "margin_yield",          # (solo L1) margine strutturale
    "margin_temp",           # (solo L1) margine termico
)

#: Frazione minima di p_c che il Dp dell'iniettore combustibile deve garantire.
#: Sotto questa soglia l'iniezione si accoppia con l'acustica di camera. E' una
#: SCELTA di progetto, dichiarata qui e non sepolta in una formula.
MIN_INJECTOR_DP_FRACTION = 0.15

#: Lo stesso, lato ARIA. Costante separata perche' le due soglie sono scelte
#: indipendenti, anche se oggi coincidono.
#:
#: PERCHE' ESISTE. Nella prima versione questo vincolo NON c'era, e la prima
#: run di NSGA-II lo ha scoperto e sfruttato: ha portato `d_ox_ratio` al
#: massimo del box (fori d'aria larghissimi -> velocita' bassa -> Dp d'aria
#: quasi nullo) per potersi permettere p_c = 5.94 bar contro i 6.0 bar di
#: serbatoio. Sulla carta piu' spinta; in giardino un motore in cui la pressione
#: di camera comanda la portata d'aria invece del contrario, cioe' il classico
#: accoppiamento fra camera e alimentazione.
#:
#: L'aria e' il 94 % della portata: qui il disaccoppiamento conta PIU' che sul
#: combustibile, non meno. Con questo vincolo attivo il massimo ammesso diventa
#: p_c <= p_aria / 1.15 = 5.22 bar, che e' esattamente il valore ricavato a mano
#: in config/design_default.yaml: due strade indipendenti allo stesso numero.
MIN_OX_INJECTOR_DP_FRACTION = 0.15

#: Damkohler minimo ammesso, Da = tau_residenza / tau_blowout.
#:
#: Non e' una preferenza di progetto: e' il CONFINE DI VALIDITA' del modello L0.
#: Tutta la prestazione L0 poggia sull'ipotesi di equilibrio chimico raggiunto
#: in camera (l0/equilibrium.py). Se il gas esce prima che la cinetica abbia
#: finito, quel numero non e' ottimistico: e' falso, e ottimizzarlo produce un
#: motore che sulla carta va e in giardino no.
#:
#: Da = 1 e' la soglia di spegnimento di un reattore PERFETTAMENTE miscelato.
#: Una camera reale non lo e': ha una distribuzione dei tempi di residenza, e
#: la frazione di fluido che attraversa piu' in fretta della media e' quella
#: che decide. Il fattore 5 copre quella coda. E' una scelta, dichiarata qui.
MIN_DAMKOHLER = 5.0


def objectives_l0(
    run_id: str,
    l0: L0Result,
    op: OperatingPoint,
    p_c: float,
    geometry: GeometryArtifact | None = None,
    min_feature_size: float | None = None,
    derived: Mapping[str, float] | None = None,
    tau_chem: float | None = None,
    q_removable: float | None = None,
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
        "neg_isp_fuel": -l0.Isp_fuel_s,
    }
    g: dict[str, float] = {}
    source: dict[str, str] = {k: "l0.cycle" for k in f}

    d = dict(derived or {})

    # --- severita' termica ------------------------------------------------- #
    # q_throat resta None se Cantera non ha potuto dare le proprieta' di
    # trasporto: in quel caso NON entra fra gli obiettivi, invece di entrarci
    # come zero (che vorrebbe dire "gola fredda", il contrario del vero).
    if l0.q_throat is not None:
        f["q_throat"] = l0.q_throat
        source["q_throat"] = "l0.cycle/bartz"
        if q_removable is not None:
            # q <= q_estraibile  ->  g = (q - q_max) / q_max
            g["throat_heat_flux"] = (l0.q_throat - q_removable) / q_removable
            source["throat_heat_flux"] = "cooling"

    # --- completezza della combustione -------------------------------------- #
    # tau_res arriva dalla geometria, tau_chem dalla cinetica: servono ENTRAMBI,
    # e se manca uno dei due il Damkohler non e' zero, e' ignoto.
    if tau_chem is not None and l0.tau_res is not None and tau_chem > 0.0:
        da = l0.tau_res / tau_chem
        f["neg_damkohler"] = -da
        source["neg_damkohler"] = "l0.chemistry"
        g["combustion_residence"] = (MIN_DAMKOHLER - da) / MIN_DAMKOHLER
        source["combustion_residence"] = "l0.chemistry"

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
        g["ox_dp_stability"] = (MIN_OX_INJECTOR_DP_FRACTION * p_c - dp_ox) / p_c
        g["ox_supply_pressure"] = (p_c + dp_ox - op.p_air_supply) / op.p_air_supply
        source["ox_dp_stability"] = "geometry.parameters"
        source["ox_supply_pressure"] = "geometry.parameters"

    if op.mdot_air_max is not None:
        g["air_supply"] = (l0.mdot_air - op.mdot_air_max) / op.mdot_air_max
        source["air_supply"] = "operating_point"

    # --- vincoli geometrici ----------------------------------------------- #
    if geometry is not None:
        f["wall_volume"] = geometry.volume       # volume VERO, dal CAD
        source["wall_volume"] = "geometry.build"
        g["watertight"] = 0.0 if geometry.is_watertight_mesh else 1.0
        source["watertight"] = "geometry.build"
    elif l0.wall_volume is not None:
        # Volume in forma chiusa, senza costruire il solido: coincide con OCCT
        # a precisione di macchina (test_objectives lo misura) e costa
        # microsecondi invece di ~100 ms. `source` dice comunque quale dei due
        # e' finito qui: in un database da migliaia di run, non poterli
        # distinguere sarebbe irreparabile.
        f["wall_volume"] = l0.wall_volume
        source["wall_volume"] = "geometry.profile/esatto"

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
    tau_chem: float | None = None,
    q_removable: float | None = None,
) -> Objectives:
    """Obiettivi completi, con i margini dal FEM.

    Nota sulla convenzione: `FEMResult.margin_yield` e' un MARGINE (positivo =
    buono), mentre `g` vuole il contrario. L'inversione avviene qui, una volta
    sola, invece che in ogni punto d'uso.
    """
    base = objectives_l0(run_id, l0, op, p_c, geometry, min_feature_size,
                         derived, tau_chem, q_removable)
    g = dict(base.g)
    g["margin_yield"] = -fem.margin_yield
    g["margin_temp"] = -fem.margin_temp
    src = dict(base.source)
    src["margin_yield"] = "l1.fem"
    src["margin_temp"] = "l1.fem"
    return Objectives(run_id=run_id, fidelity="L1", f=base.f, g=g, source=src)
