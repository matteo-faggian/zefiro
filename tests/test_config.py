"""La configurazione del repo deve caricarsi, e i TODO devono restare TODO."""
from __future__ import annotations

import pytest

from zefiro.config import CONFIG_DIR, load_design_vector, load_material, load_operating_point
from zefiro.schemas import MissingDatum


def test_operating_point_si_carica_e_conserva_i_none():
    op = load_operating_point()
    assert op.p_air_supply == pytest.approx(10.0e5)
    assert op.p_fuel_supply == pytest.approx(8.0e5)
    # i TODO NON devono essere stati sostituiti da default plausibili
    assert op.mdot_air_max is None
    assert op.T_air_in is None


def test_usare_un_todo_fallisce_con_messaggio_utile():
    op = load_operating_point()
    with pytest.raises(MissingDatum, match="mdot_air_max"):
        op.require("mdot_air_max")


def test_design_default_rispetta_i_bounds():
    x = load_design_vector()          # il costruttore valida i bounds
    assert x.values["N_inj"] == 12.0
    assert x.values["p_c"] == pytest.approx(4.0e5)


def test_materiale_316l_non_contiene_valori_inventati():
    """TODO n.6: mettere i valori del 316L LAMINATO darebbe margini sbagliati
    con l'aria di essere giusti. Devono restare null finche' non li misuri."""
    m = load_material()
    assert m["density_kg_m3"] is None
    assert m["yield_strength_Pa_Z"] is None
    assert m["process"]["min_feature_size_m"] is None
    assert m["anisotropic"] is True


def test_i_file_di_config_esistono():
    assert (CONFIG_DIR / "operating_point.yaml").exists()
    assert (CONFIG_DIR / "design_default.yaml").exists()
    assert (CONFIG_DIR / "materials" / "aisi316l.yaml").exists()
