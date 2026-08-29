"""Bilancio L0 completo: dai parametri e dal punto operativo a un L0Result."""
from __future__ import annotations

import math

from zefiro.l0.equilibrium import chamber_state
from zefiro.l0.mixture import MixtureModel
from zefiro.l0.nozzle import AEROSPIKE_MEANINGFUL_PR, SHIFTING, expand_to_pressure
from zefiro.schemas import L0Result, OperatingPoint
from zefiro.units import G0


def evaluate_l0(
    op: OperatingPoint,
    p_c: float,
    phi_core: float,
    f_film: float = 0.0,
    mdot_air: float | None = None,
    expansion_mode: str = SHIFTING,
    chamber_volume: float | None = None,
    thermal_severity: bool = True,
    wall_temperature: float = 800.0,
) -> L0Result:
    """Valutazione L0 di un punto di progetto.

    Parameters
    ----------
    p_c : pressione di camera [Pa]
    phi_core : rapporto di equivalenza del core (tutta l'aria, il combustibile
        al netto del film)
    f_film : frazione del combustibile totale destinata al film cooling
    mdot_air : portata d'aria [kg/s]. Se None si usa `op.mdot_air_max`, che e'
        un TODO bloccante finche' non fornito.
    chamber_volume : se noto (dalla geometria), abilita il calcolo di tau_res.
    """
    op.require("T_air_in", "T_fuel_in")
    if mdot_air is None:
        op.require("mdot_air_max")
        mdot_air = float(op.mdot_air_max)  # type: ignore[arg-type]

    if p_c >= op.p_fuel_supply:
        raise ValueError(
            f"p_c = {p_c/1e5:.2f} bar >= pressione di alimentazione GPL "
            f"({op.p_fuel_supply/1e5:.2f} bar): il combustibile non puo' entrare in camera."
        )

    model = MixtureModel.from_fuel(op.fuel)
    afr = model.afr_stoichiometric()
    mdot_f_core, mdot_f_film, phi_global = model.mass_flows(mdot_air, phi_core, f_film)
    mdot_fuel_tot = mdot_f_core + mdot_f_film
    mdot_tot = mdot_air + mdot_fuel_tot

    warnings: list[str] = []
    assumptions: list[str] = [
        "equilibrio chimico completo in camera (limite superiore di T_ad e c*)",
        "velocita' in camera trascurabile: T statica = T di ristagno",
        "camera adiabatica: nessuna perdita verso la parete",
        f"espansione isentropica 1D, modo = {expansion_mode}",
        "aerospike al punto di progetto: p_e = p_amb (espansione adattata)",
        "nessun modello di pressione di base per il plug troncato: "
        "la prestazione e' un limite superiore",
    ]

    # --- temperatura di fiamma del CORE: e' cio' che vede la fiamma --------- #
    core = chamber_state(model, p_c, op.T_air_in, op.T_fuel_in, mdot_air, mdot_f_core)
    w = model.equilibrium_validity_warning(phi_core)
    if w:
        warnings.append(w)

    # --- stato per la prestazione: tutto il combustibile mescolato e bruciato #
    if f_film > 0.0:
        assumptions.append(
            "film cooling completamente mescolato e bruciato a monte della gola: "
            "limite SUPERIORE di prestazione (il film reale brucia in parte o affatto)"
        )
        perf = chamber_state(model, p_c, op.T_air_in, op.T_fuel_in, mdot_air, mdot_fuel_tot)
        w2 = model.equilibrium_validity_warning(phi_global)
        if w2:
            warnings.append(w2)
    else:
        perf = core

    noz = expand_to_pressure(model.gas, perf, p_e=op.p_amb, p_amb=op.p_amb, mode=expansion_mode)

    A_t = mdot_tot * noz.A_t_over_mdot
    thrust = noz.C_F * p_c * A_t

    pr = p_c / op.p_amb
    if pr < AEROSPIKE_MEANINGFUL_PR:
        warnings.append(
            f"p_c/p_amb = {pr:.2f} < {AEROSPIKE_MEANINGFUL_PR}: con eps = {noz.epsilon:.2f} "
            "l'ugello e' appena supersonico e il vantaggio di compensazione di quota "
            "dell'aerospike e' trascurabile. La scelta va giustificata come studio del "
            "metodo, non come scelta di prestazione."
        )
    if p_c > 0.8 * op.p_fuel_supply:
        warnings.append(
            f"p_c = {p_c/1e5:.2f} bar lascia meno del 20 % di Dp all'iniettore GPL "
            f"(alimentazione {op.p_fuel_supply/1e5:.2f} bar): rischio di accoppiamento "
            "fra iniezione e acustica di camera."
        )

    # --- severita' termica in gola (Bartz) ---------------------------------- #
    q_throat = T_aw = None
    if thermal_severity:
        try:
            from zefiro.thermal import adiabatic_wall_temperature, bartz_h_gas

            g = model.gas
            g.TPX = perf.T, perf.p, perf.X
            g.equilibrate("TP")
            mu, cpg = g.viscosity, g.cp_mass
            Pr = mu * cpg / g.thermal_conductivity
            D_t = 2.0 * math.sqrt(A_t / math.pi)
            T_aw = adiabatic_wall_temperature(perf.T, 1.0, perf.gamma, Pr)
            h = bartz_h_gas(D_t, p_c, noz.c_star, mu, cpg, Pr, perf.gamma,
                            1.0, 1.0, wall_temperature, perf.T)
            q_throat = h * (T_aw - wall_temperature)
            assumptions.append(
                "flusso termico di gola da correlazione di Bartz (empirica, +-30 %), "
                f"a parete assunta a {wall_temperature:.0f} K"
            )
        except Exception:  # noqa: BLE001
            q_throat = T_aw = None

    tau_res = None
    if chamber_volume is not None:
        rho_c = perf.p * perf.MW / (8314.462618153242 * perf.T)
        tau_res = rho_c * chamber_volume / mdot_tot

    return L0Result(
        phi_core=phi_core,
        phi_global=phi_global,
        AFR_stoich=afr,
        mdot_air=mdot_air,
        mdot_fuel_core=mdot_f_core,
        mdot_fuel_film=mdot_f_film,
        T_ad=core.T,
        X_eq=core.X,
        gamma_c=core.gamma,
        MW_c=core.MW,
        cp_c=core.cp,
        c_star=noz.c_star,
        M_e=noz.M_e,
        epsilon=noz.epsilon,
        C_F=noz.C_F,
        Isp_s=thrust / (mdot_tot * G0),
        Isp_fuel_s=thrust / (mdot_fuel_tot * G0),
        thrust=thrust,
        A_t=A_t,
        q_throat=q_throat,
        T_wall_adiabatic=T_aw,
        tau_res=tau_res,
        tau_chem=None,
        damkohler=None,
        warnings=tuple(warnings),
        assumptions=tuple(assumptions),
    )
