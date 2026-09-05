"""La configurazione del repo deve caricarsi, essere coerente, e i TODO
rimasti devono restare TODO."""
from __future__ import annotations

import pytest

from pathlib import Path

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
    assert op.mdot_air_max == pytest.approx(0.0470)
    assert op.burn_time == pytest.approx(BURN_TIME)
    # NON un numero scritto a mano: e' p_sat(propano, 15 C) da CoolProp.
    from zefiro.feed import saturation_pressure
    assert op.p_fuel_supply == pytest.approx(saturation_pressure({"C3H8": 1.0}, 288.15))


def test_la_pressione_daria_di_config_e_la_MINIMA_non_la_massima():
    """Semantica del campo, ed e' facile sbagliarla: il riduttore perde il
    controllo alla FINE della raffica, quando il serbatoio e' al minimo.
    Dimensionare sui 10 bar iniziali significherebbe progettare per un istante
    che dura zero."""
    op = load_operating_point()
    assert op.p_air_supply == pytest.approx(op.p_fuel_supply)
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
    # I due lati si chiudono INSIEME per costruzione: il fondo scarico del
    # serbatoio e' posto uguale a p_sat, che e' il solo valore che non spreca
    # aria ne' lascia il riduttore dell'aria a corto prima di quello del GPL.
    assert p_c * (1.0 + MIN_INJECTOR_DP_FRACTION) <= op.p_fuel_supply * 1.001


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


def test_la_cartella_delle_run_e_configurabile(monkeypatch, tmp_path):
    """Sotto WSL il repo sta su drvfs, dove l'I/O su file piccoli e numerosi e'
    5-10 volte piu' lento. Gli artefatti devono poter stare altrove."""
    from zefiro.cli import runs_root

    monkeypatch.delenv("ZEFIRO_RUNS", raising=False)
    assert runs_root() == Path("runs")
    monkeypatch.setenv("ZEFIRO_RUNS", str(tmp_path / "altrove"))
    assert runs_root() == tmp_path / "altrove"


# --------------------------------------------------------------------------- #
# La pressione del GPL e' termodinamica, non configurazione
# --------------------------------------------------------------------------- #
def _config_modificata(tmp_path, **campi):
    """Copia operating_point.yaml cambiando alcuni campi di primo livello."""
    import yaml

    d = yaml.safe_load((CONFIG_DIR / "operating_point.yaml").read_text(encoding="utf-8"))
    for k, v in campi.items():
        if k == "T_design_K":
            d["fuel"]["bottle"]["T_design_K"] = v
        else:
            d[k] = v
    q = tmp_path / "op.yaml"
    q.write_text(yaml.safe_dump(d), encoding="utf-8")
    return q


def test_una_pressione_di_bombola_impossibile_viene_rifiutata(tmp_path):
    """8 bar a 15 C non esistono: il propano PURO ne fa 7.32, e non c'e' miscela
    piu' volatile del propano puro fra quelle di una bombola da barbecue.

    E' il difetto che questo controllo esiste per intercettare: fino alla
    revisione 0.4.0 il file conteneva 8.0 bar scritti a mano, cioe' la
    tensione di vapore a 18.3 C, e l'intero punto operativo era valido solo
    sopra quella temperatura senza che nulla lo dicesse."""
    from zefiro.schemas import MissingDatum

    q = _config_modificata(tmp_path, p_fuel_supply_bar=8.0)
    with pytest.raises(MissingDatum, match="satura"):
        load_operating_point(q)


def test_la_stessa_pressione_e_accettata_se_la_bombola_e_calda(tmp_path):
    """Non e' il numero 8.0 a essere vietato: e' 8.0 A 15 GRADI. A 20 C la
    bombola satura a 8.36 bar e 8.0 diventa un dato legittimo."""
    q = _config_modificata(tmp_path, p_fuel_supply_bar=8.0, T_design_K=293.15)
    assert load_operating_point(q).p_fuel_supply == pytest.approx(8.0e5)


def test_senza_temperatura_di_progetto_non_c_e_punto_operativo(tmp_path):
    from zefiro.schemas import MissingDatum

    q = _config_modificata(tmp_path, T_design_K=None)
    with pytest.raises(MissingDatum, match="T_design_K"):
        load_operating_point(q)


def test_la_frazione_di_dp_del_loader_e_quella_dell_ottimizzatore():
    """_DP_FRACTION e' duplicata in config.py per non far dipendere il
    caricamento della configurazione dall'ottimizzatore. La duplicazione e'
    lecita solo se qualcuno verifica che i due numeri coincidano."""
    from zefiro.config import _DP_FRACTION
    from zefiro.opt.objectives import MIN_INJECTOR_DP_FRACTION

    assert _DP_FRACTION == MIN_INJECTOR_DP_FRACTION


def test_la_bombola_piu_calda_accorcia_la_raffica():
    """Il risultato controintuitivo che dimensiona il motore.

    Una bombola piu' calda alza p_c e quindi l'Isp. Ma il fondo scarico del
    serbatoio d'aria vale p_sat, quindi alzarlo toglie massa utilizzabile
    dal serbatoio: fra 10 bar e p_sat resta meno aria, e la raffica si
    accorcia. Il requisito dei 5 s cade a 18.8 C.

    Conseguenza pratica: la bombola va tenuta FRESCA, non tiepida."""
    from zefiro.feed import saturation_pressure, tank_blowdown

    durate = []
    for T in (288.15, 293.15):
        p = saturation_pressure({"C3H8": 1.0}, T)
        durate.append(tank_blowdown(0.100, TANK_P_MAX, p, 293.15).mass_usable_adiabatic)
    assert durate[1] < durate[0]
