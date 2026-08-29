"""Obiettivi e vincoli: convenzione dei segni e adimensionalizzazione.

La convenzione (minimizzare `f`, `g <= 0` fattibile) e' facile da rompere e
impossibile da accorgersene guardando i numeri: un segno invertito produce un
ottimizzatore che converge tranquillamente sulla soluzione peggiore.
"""
from __future__ import annotations

import dataclasses

import pytest

from zefiro.geometry.parameters import derive
from zefiro.opt.objectives import (
    CONSTRAINT_NAMES,
    MIN_INJECTOR_DP_FRACTION,
    OBJECTIVE_NAMES,
    objectives_l0,
    objectives_l1,
)
from zefiro.schemas import FEMResult


@pytest.fixture
def evaluated(operating_point, design_vector):
    params, l0 = derive(design_vector, operating_point)
    return params, l0, design_vector.values["p_c"]


def test_gli_obiettivi_da_massimizzare_entrano_col_segno_invertito(
    operating_point, evaluated
):
    params, l0, p_c = evaluated
    o = objectives_l0("r", l0, operating_point, p_c, derived=params.derived)
    assert o.f["neg_thrust"] == -l0.thrust
    assert o.f["neg_isp_total"] == -l0.Isp_s
    # spinta maggiore -> obiettivo minore: e' cio' che l'ottimizzatore cerca
    migliore = dataclasses.replace(l0, thrust=l0.thrust * 2.0)
    o2 = objectives_l0("r", migliore, operating_point, p_c, derived=params.derived)
    assert o2.f["neg_thrust"] < o.f["neg_thrust"]


def test_i_nomi_prodotti_stanno_nel_registro(operating_point, evaluated):
    """Il registro sono le colonne del database: un nome fuori registro
    verrebbe silenziosamente perso all'inserimento."""
    params, l0, p_c = evaluated
    o = objectives_l0("r", l0, operating_point, p_c, min_feature_size=0.4e-3,
                      derived=params.derived)
    assert set(o.f) <= set(OBJECTIVE_NAMES)
    assert set(o.g) <= set(CONSTRAINT_NAMES)


def test_vincolo_di_stabilita_delliniezione(operating_point, evaluated):
    """g = (0.15 p_c - Dp) / p_c: negativo quando il Dp e' sufficiente."""
    params, l0, p_c = evaluated
    dp = params.derived["dp_inj_fuel"]
    o = objectives_l0("r", l0, operating_point, p_c, derived=params.derived)
    atteso = (MIN_INJECTOR_DP_FRACTION * p_c - dp) / p_c
    assert o.g["fuel_dp_stability"] == pytest.approx(atteso, rel=1e-12)
    assert (o.g["fuel_dp_stability"] <= 0.0) == (dp >= MIN_INJECTOR_DP_FRACTION * p_c)


def test_vincolo_di_pressione_di_alimentazione(operating_point, evaluated):
    """p_c + Dp deve stare sotto la pressione di bombola."""
    params, l0, p_c = evaluated
    o = objectives_l0("r", l0, operating_point, p_c, derived=params.derived)
    fattibile = p_c + params.derived["dp_inj_fuel"] <= operating_point.p_fuel_supply
    assert (o.g["fuel_supply_pressure"] <= 0.0) == fattibile


def test_i_vincoli_sono_adimensionali(operating_point, evaluated):
    """Sommare un Dp in pascal con un margine strutturale puro darebbe alla
    pressione un peso arbitrario di 10^5 nella violazione aggregata."""
    params, l0, p_c = evaluated
    o = objectives_l0("r", l0, operating_point, p_c, min_feature_size=0.4e-3,
                      derived=params.derived)
    for name, v in o.g.items():
        assert abs(v) < 1.0e3, f"{name} = {v}: sembra dimensionale"


def test_un_vincolo_non_valutabile_resta_fuori(operating_point, evaluated):
    """Metterlo a zero direbbe 'soddisfatto', che e' falso. La differenza fra
    'verificato' e 'non verificabile' non va nascosta."""
    params, l0, p_c = evaluated
    senza_cd = dataclasses.replace(operating_point, cd_injector_fuel=None,
                                   cd_injector_ox=None)
    o = objectives_l0("r", l0, senza_cd, p_c, derived={})
    assert "fuel_dp_stability" not in o.g
    assert "min_feature" not in o.g


def test_fabbricabilita_come_vincolo(operating_point, evaluated):
    """Il foro piu' piccolo deve superare min_feature_size."""
    params, l0, p_c = evaluated
    d = params.derived
    quota_min = min(d[k] for k in ("d_ox", "d_fuel", "d_film", "film_land", "t_wall")
                    if d[k] > 0.0)
    o = objectives_l0("r", l0, operating_point, p_c, min_feature_size=0.4e-3,
                      derived=d)
    assert o.g["min_feature"] == pytest.approx((0.4e-3 - quota_min) / 0.4e-3, rel=1e-12)
    # a questa scala il vincolo e' ATTIVO: e' il risultato di progetto, non un bug
    assert o.g["min_feature"] > 0.0


def test_il_margine_del_fem_entra_col_segno_invertito(operating_point, evaluated):
    """FEMResult.margin_yield e' un margine (positivo = buono), `g` vuole il
    contrario. L'inversione avviene una volta sola, in objectives_l1."""
    params, l0, p_c = evaluated
    fem = FEMResult(run_id="r", cfd_id="c", T_wall_max=900.0, von_mises_max=1.0e8,
                    margin_yield=1.5, margin_temp=0.3, solver_version="test")
    o = objectives_l1("r", l0, operating_point, p_c, fem, derived=params.derived)
    assert o.g["margin_yield"] == -1.5
    assert o.g["margin_temp"] == -0.3
    assert o.fidelity == "L1"


def test_feasible_richiede_tutti_i_vincoli_non_positivi(operating_point, evaluated):
    params, l0, p_c = evaluated
    o = objectives_l0("r", l0, operating_point, p_c, derived=params.derived)
    assert o.feasible == all(v <= 0.0 for v in o.g.values())


# --- proxy di massa contro il CAD ------------------------------------------- #

@pytest.mark.slow
def test_volume_in_forma_chiusa_coincide_col_cad(operating_point, design_vector):
    """`objectives_l0` usa la formula chiusa quando il CAD non c'e'.

    E' quel numero che l'ottimizzatore minimizza per decine di migliaia di
    valutazioni: se divergesse dal volume vero, si ottimizzerebbe una massa
    che non esiste. Qui si costruisce davvero il solido con OCCT e si misura.

    La soglia e' 1e-9 e non una percentuale perche' la formula non e'
    un'approssimazione: il solido E' la rivoluzione di quel poligono. Una
    versione precedente stimava a parete sottile e sbagliava del 7.7 %,
    contando due volte il materiale agli spigoli concavi.
    """
    from zefiro.geometry.build import BuildOptions, build_solid
    from zefiro.geometry.parameters import derive

    params, l0 = derive(design_vector, operating_point)
    # senza i fori d'iniezione: il proxy non li modella, e confrontare un
    # solido forato con una stima non forata misurerebbe i fori, non il proxy.
    solido = build_solid(params, BuildOptions(include_injection=False))
    v_cad = solido.volume * 1.0e-9                     # mm^3 -> m^3
    errore = abs(l0.wall_volume - v_cad) / v_cad
    assert errore < 1.0e-9, (
        f"proxy {l0.wall_volume*1e6:.1f} cm^3 contro CAD "
        f"{v_cad*1e6:.1f} cm^3: errore {errore:.1%}"
    )


def test_formula_e_cad_non_si_confondono(operating_point, design_vector):
    """`source` deve dire quale dei due volumi e' finito nell'obiettivo:
    confonderli in un database di mille run e' irreparabile."""
    from zefiro.geometry.parameters import derive
    from zefiro.opt.objectives import objectives_l0

    params, l0 = derive(design_vector, operating_point)
    o = objectives_l0("r", l0, operating_point, design_vector.values["p_c"],
                      geometry=None, derived=params.derived)
    assert o.source["wall_volume"] == "geometry.profile/esatto"


def test_damkohler_assente_non_diventa_zero(operating_point, design_vector):
    """Senza tau_chem il vincolo di residenza deve restare FUORI dal dizionario:
    metterlo a zero direbbe 'verificato e soddisfatto', che e' falso."""
    from zefiro.geometry.parameters import derive
    from zefiro.opt.objectives import objectives_l0

    params, l0 = derive(design_vector, operating_point)
    o = objectives_l0("r", l0, operating_point, design_vector.values["p_c"],
                      derived=params.derived, tau_chem=None)
    assert "combustion_residence" not in o.g
    assert "neg_damkohler" not in o.f

    o2 = objectives_l0("r", l0, operating_point, design_vector.values["p_c"],
                       derived=params.derived, tau_chem=40.0e-6)
    assert "combustion_residence" in o2.g
    da = -o2.f["neg_damkohler"]
    assert da == pytest.approx(l0.tau_res / 40.0e-6)


def test_vincolo_di_dp_aria_esiste_ed_e_simmetrico(operating_point, design_vector):
    """Regressione su un difetto TROVATO DALL'OTTIMIZZATORE.

    Senza `ox_dp_stability` NSGA-II allargava i fori d'aria fino ad azzerarne
    il Dp, per guadagnare pressione di camera contro il serbatoio. Il vincolo
    lato combustibile c'era e quello lato aria no: un'asimmetria che nessuno
    aveva scritto apposta, e che l'ottimizzatore ha trovato in una run.

    Qui si verifica che i due vincoli esistano entrambi e abbiano la STESSA
    forma, cosi' l'asimmetria non puo' tornare di nascosto.
    """
    from zefiro.geometry.parameters import derive
    from zefiro.opt.objectives import (
        MIN_INJECTOR_DP_FRACTION,
        MIN_OX_INJECTOR_DP_FRACTION,
        objectives_l0,
    )

    params, l0 = derive(design_vector, operating_point)
    p_c = design_vector.values["p_c"]
    o = objectives_l0("r", l0, operating_point, p_c, derived=params.derived)

    assert "ox_dp_stability" in o.g and "fuel_dp_stability" in o.g
    d = params.derived
    assert o.g["ox_dp_stability"] == pytest.approx(
        (MIN_OX_INJECTOR_DP_FRACTION * p_c - d["dp_inj_ox"]) / p_c)
    assert o.g["fuel_dp_stability"] == pytest.approx(
        (MIN_INJECTOR_DP_FRACTION * p_c - d["dp_inj_fuel"]) / p_c)


def test_il_tetto_di_p_c_e_quello_ricavato_a_mano():
    """Con Dp_aria >= 15 % di p_c e serbatoio a 6 bar, il massimo ammesso e'
    p_c = 6/1.15 = 5.22 bar, cioe' il valore in config/design_default.yaml.
    Se qualcuno cambia la soglia senza accorgersene, questo test lo dice."""
    from zefiro.opt.objectives import MIN_OX_INJECTOR_DP_FRACTION

    p_c_max = 6.0e5 / (1.0 + MIN_OX_INJECTOR_DP_FRACTION)
    assert p_c_max / 1e5 == pytest.approx(5.22, abs=0.01)
