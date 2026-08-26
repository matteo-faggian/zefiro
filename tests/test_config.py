"""La configurazione del repo deve caricarsi, essere coerente, e i TODO
rimasti devono restare TODO."""
from __future__ import annotations

import pytest

from zefiro.config import CONFIG_DIR, load_design_vector, load_material, load_operating_point
from zefiro.feed import tank_blowdown
from zefiro.schemas import MissingDatum

BURN_TIME = 5.0
TANK_VOLUME = 0.100
TANK_P_MAX = 10.0e5


def test_i_file_di_config_esistono():
    assert (CONFIG_DIR / "operating_point.yaml").exists()
    assert (CONFIG_DIR / "design_default.yaml").exists()
    assert (CONFIG_DIR / "materials" / "aisi316l.yaml").exists()


def test_il_punto_di_progetto_deciso_e_caricato():
    op = load_operating_point()
    assert op.mdot_air_max == pytest.approx(0.0718)
    assert op.burn_time == pytest.approx(BURN_TIME)
    assert op.p_fuel_supply == pytest.approx(8.0e5)


def test_la_pressione_daria_di_config_e_la_MINIMA_non_la_massima():
    """Semantica del campo, ed e' facile sbagliarla: il riduttore perde il
    controllo alla FINE della raffica, quando il serbatoio e' al minimo.
    Dimensionare sui 10 bar iniziali significherebbe progettare per un istante
    che dura zero."""
    op = load_operating_point()
    assert op.p_air_supply == pytest.approx(6.0e5)
    assert op.p_air_supply < TANK_P_MAX


def test_la_portata_di_progetto_e_sostenibile_per_i_5_secondi():
    """Controllo incrociato fra config e fisica dell'impianto: la massa
    estraibile dal serbatoio piu' la ricarica del compressore deve bastare a
    tenere quella portata per la durata dichiarata.

    E' il test che si accorgerebbe se qualcuno cambiasse `mdot_air_max` senza
    ricalcolare la durata, o viceversa.
    """
    op = load_operating_point()
    bd = tank_blowdown(TANK_VOLUME, TANK_P_MAX, op.p_air_supply, 293.15)
    richiesta = op.mdot_air_max * BURN_TIME
    # La ricarica del compressore NON entra nel conto: e' ancora un TODO e non
    # si dimensiona su un numero non misurato. Quando ci sara', sara' margine.
    assert richiesta <= bd.mass_usable_adiabatic, (
        f"richiesti {richiesta*1e3:.0f} g, disponibili "
        f"{bd.mass_usable_adiabatic*1e3:.0f} g senza ricarica"
    )
    # e non deve avanzare massa: avanzarne significherebbe buttare via spinta
    assert richiesta > 0.95 * bd.mass_usable_adiabatic


def test_la_pressione_di_camera_rispetta_il_vincolo_lato_aria():
    """p_c + Dp_ox <= p_aria a fine raffica, con Dp_ox >= 15 % di p_c."""
    from zefiro.opt.objectives import MIN_INJECTOR_DP_FRACTION

    op = load_operating_point()
    x = load_design_vector()
    p_c = x.values["p_c"]
    assert p_c * (1.0 + MIN_INJECTOR_DP_FRACTION) <= op.p_air_supply * 1.001
    # e il lato GPL deve avere ancora margine: il collo di bottiglia si e'
    # spostato sull'aria, ed e' un fatto di progetto da non perdere
    assert p_c * (1.0 + MIN_INJECTOR_DP_FRACTION) < op.p_fuel_supply


def test_il_film_cooling_e_disattivato_ma_la_capacita_resta():
    """L'analisi termica a 5 s mostra che il film iniettato in testa protegge
    una camera che non ne ha bisogno e non arriva alla gola, che invece fonde.
    Quindi f_film = 0 nel progetto di default, ma i bounds restano aperti per
    quando il raffreddamento verra' rimesso alla gola."""
    from zefiro.geometry.parameters import DESIGN_BOUNDS

    assert load_design_vector().values["f_film"] == 0.0
    assert DESIGN_BOUNDS["f_film"][1] > 0.0


def test_design_default_rispetta_i_bounds():
    x = load_design_vector()          # il costruttore valida i bounds
    assert x.values["N_inj"] == 12.0


def test_i_todo_rimasti_sono_ancora_todo():
    """Quello che non e' stato misurato non deve essere stato riempito con
    valori plausibili."""
    op = load_operating_point()
    assert op.T_air_in is None
    assert op.T_fuel_in is None
    assert op.cd_injector_ox is None
    assert op.cd_injector_fuel is None
    with pytest.raises(MissingDatum, match="T_air_in"):
        op.require("T_air_in", "T_fuel_in")


def test_materiale_316l_non_contiene_valori_inventati():
    """TODO n.J: mettere i valori del 316L LAMINATO darebbe margini sbagliati
    con l'aria di essere giusti. Restano null finche' non li misuri.

    Nota: l'analisi termica gia' fatta usa valori PLACEHOLDER dichiarati come
    tali negli script, non presi da qui."""
    m = load_material()
    assert m["density_kg_m3"] is None
    assert m["yield_strength_Pa_Z"] is None
    assert m["thermal_conductivity_W_mK"] is None
    assert m["process"]["min_feature_size_m"] is None
    assert m["anisotropic"] is True
