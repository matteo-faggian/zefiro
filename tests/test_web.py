"""Interfaccia web: contratto dell'API, errori strutturati, coda dei lavori.

Il test piu' importante di questo file e' `test_la_gui_produce_lo_stesso_run_id
_della_riga_di_comando`: e' quello che rende verificabile la regola
architetturale su cui poggia tutto, cioe' che la GUI sia un CLIENT delle stesse
funzioni e non una seconda strada per far girare i calcoli.
"""
from __future__ import annotations

import time

import pytest

pytest.importorskip("fastapi")
from fastapi.testclient import TestClient  # noqa: E402

from zefiro.geometry.parameters import DESIGN_BOUNDS  # noqa: E402
from zefiro.schemas import (  # noqa: E402
    ContractViolation,
    MissingDatum,
    MissingThermoData,
)
from zefiro.web import create_app  # noqa: E402
from zefiro.web.errors import DOMAIN_CODES, ApiError, to_api_error  # noqa: E402

API = "/api/v1"
#: Valori PROVVISORI per i dati non ancora misurati: servono a far girare i
#: test, non descrivono il motore reale.
OPERATING = {"T_air_in": 293.0, "T_fuel_in": 283.0,
             "cd_injector_ox": 0.75, "cd_injector_fuel": 0.75}


@pytest.fixture(scope="module")
def client(tmp_path_factory):
    import os

    os.environ["ZEFIRO_RUNS"] = str(tmp_path_factory.mktemp("runs"))
    with TestClient(create_app()) as c:
        yield c
    os.environ.pop("ZEFIRO_RUNS", None)


@pytest.fixture(scope="module")
def base_design(client):
    return {f["name"]: f["value"] for f in client.get(f"{API}/schema").json()["design"]}


# --------------------------------------------------------------------------- #
# Contratto
# --------------------------------------------------------------------------- #
def test_health(client):
    d = client.get(f"{API}/health").json()
    assert d["status"] == "ok"
    assert isinstance(d["dati_mancanti"], list)


def test_lo_schema_espone_tutti_i_parametri_liberi(client):
    s = client.get(f"{API}/schema").json()
    nomi = {f["name"] for f in s["design"]}
    assert nomi == set(DESIGN_BOUNDS), nomi ^ set(DESIGN_BOUNDS)
    for f in s["design"]:
        lo, hi = DESIGN_BOUNDS[f["name"]]
        assert f["min"] <= f["value"] <= f["max"]
        assert f["help"], f"{f['name']} senza spiegazione"


def test_lo_schema_dichiara_quali_dati_mancano(client):
    s = client.get(f"{API}/schema").json()
    assert any(f["missing"] for f in s["operating"])
    assert s["material_is_placeholder"] is True


def test_openapi_e_generato(client):
    o = client.get("/api/openapi.json").json()
    assert o["info"]["title"] == "Zefiro"
    assert f"{API}/evaluate" in o["paths"]


def test_la_pagina_e_i_suoi_asset_sono_serviti(client):
    assert client.get("/").status_code == 200
    assert "no-cache" in client.get("/").headers.get("cache-control", "")
    js = client.get("/assets/app.js")
    assert js.status_code == 200 and len(js.content) > 5000


# --------------------------------------------------------------------------- #
# Errori: codice stabile, non stringhe libere
# --------------------------------------------------------------------------- #
def test_ogni_errore_di_dominio_ha_un_codice(client):
    """Un errore di dominio non mappato arriverebbe al client come 500, cioe'
    come se fosse colpa del server."""
    for cls in (MissingDatum, MissingThermoData, ContractViolation):
        assert cls in DOMAIN_CODES
        api = to_api_error(cls("prova"))
        assert isinstance(api, ApiError) and api.code and api.status < 500


def test_dati_mancanti_indica_quali(client, base_design):
    r = client.post(f"{API}/evaluate", json={"design": base_design, "operating": {}})
    assert r.status_code == 422
    e = r.json()["error"]
    assert e["code"] == "DATO_MANCANTE"
    assert "T_air_in" in e["details"]["fields"]


def test_fuori_dai_bounds_e_un_errore_di_contratto(client, base_design):
    bad = dict(base_design)
    bad["p_c"] = 99.0
    r = client.post(f"{API}/evaluate", json={"design": bad, "operating": OPERATING})
    assert r.status_code == 422
    assert r.json()["error"]["code"] == "CONTRATTO_VIOLATO"


def test_parametro_mancante_nel_corpo(client, base_design):
    incompleto = {k: v for k, v in base_design.items() if k != "t_wall"}
    r = client.post(f"{API}/evaluate", json={"design": incompleto, "operating": OPERATING})
    assert r.status_code == 422
    assert "t_wall" in r.json()["error"]["message"]


def test_lavoro_sconosciuto(client):
    r = client.get(f"{API}/jobs/nonesiste")
    assert r.status_code == 404 and r.json()["error"]["code"] == "LAVORO_SCONOSCIUTO"


def test_geometria_non_ancora_generata(client):
    r = client.get(f"{API}/geometry/0123456789abcdef/stl")
    assert r.status_code == 404 and r.json()["error"]["code"] == "ARTEFATTO_ASSENTE"


def test_pressione_minima_sopra_la_massima(client):
    r = client.get(f"{API}/plant", params={"p_min": 12, "p_max": 10})
    assert r.status_code == 422 and r.json()["error"]["code"] == "PARAMETRO_NON_VALIDO"


# --------------------------------------------------------------------------- #
# Valutazione
# --------------------------------------------------------------------------- #
@pytest.fixture(scope="module")
def evaluated(client, base_design):
    r = client.post(f"{API}/evaluate", json={"design": base_design, "operating": OPERATING})
    assert r.status_code == 200, r.text
    return r.json()


def test_valutazione_restituisce_cio_che_serve_al_disegno(evaluated):
    d = evaluated
    assert d["l0"]["thrust"] > 0 and d["l0"]["T_ad"] > 1500
    sec = d["section"]
    assert len(sec["profile_mm"]) > 20
    assert all(len(p) == 2 for p in sec["profile_mm"])
    assert sec["injectors"]["N"] >= 6
    assert sec["stations_mm"]["R_c"] > sec["stations_mm"]["R_lip"] > 0


def test_gli_avvisi_arrivano_al_client(evaluated):
    """Il posto peggiore per un avviso e' il log del server."""
    assert any("compensazione di quota" in w for w in evaluated["warnings"])
    assert evaluated["assumptions"]


def test_la_valutazione_e_salvata_nel_database(client, evaluated):
    assert evaluated["stored"] is True
    ids = {r["run_id"] for r in client.get(f"{API}/runs").json()["rows"]}
    assert evaluated["run_id"] in ids


def test_la_gui_produce_lo_stesso_run_id_della_riga_di_comando(client, base_design, evaluated):
    """LA REGOLA ARCHITETTURALE, resa verificabile.

    Se la GUI calcolasse per conto suo — altre unita', altri default, un
    percorso di codice diverso — produrrebbe un identificativo diverso a parita'
    di progetto, e i suoi risultati non sarebbero confrontabili con quelli della
    riga di comando ne' riproducibili. Questo test fallisce nell'istante in cui
    qualcuno introduce una scorciatoia nel backend web.
    """
    from zefiro.schemas import DesignVector
    from zefiro.store import make_run_id
    from zefiro.web.service import design_to_si, operating_from

    x = DesignVector(values=design_to_si(base_design), bounds=DESIGN_BOUNDS)
    atteso = make_run_id(x, operating_from(OPERATING))
    assert evaluated["run_id"] == atteso


def test_gli_stessi_ingressi_danno_lo_stesso_identificativo(client, base_design, evaluated):
    r = client.post(f"{API}/evaluate", json={"design": base_design, "operating": OPERATING})
    assert r.json()["run_id"] == evaluated["run_id"]


# --------------------------------------------------------------------------- #
# Termica e impianto
# --------------------------------------------------------------------------- #
def test_termica(client, base_design):
    r = client.post(f"{API}/thermal", json={"design": base_design, "operating": OPERATING})
    assert r.status_code == 200
    d = r.json()
    zone = {s["label"] for s in d["stations"]}
    assert zone == {"camera", "gola"}
    gola = next(s for s in d["stations"] if s["label"] == "gola")
    cam = next(s for s in d["stations"] if s["label"] == "camera")
    # la gola vede molto piu' flusso della camera: e' il risultato di progetto
    assert gola["q0_MW_m2"] > 5 * cam["q0_MW_m2"]
    # il raffreddamento sul dorso deve abbassare la temperatura, non alzarla
    assert gola["cooled"][-1] < gola["sink"][-1]
    assert len(d["t"]) == len(gola["sink"])
    assert d["material_is_placeholder"] is True
    assert "Bartz" in d["correlation"]


def test_impianto(client):
    d = client.get(f"{API}/plant").json()
    assert d["compressore"]["mdot_isotermo_g_s"] > d["compressore"]["mdot_adiabatico_g_s"]
    assert d["serbatoio"]["utilizzabile_adiabatico_g"] < d["serbatoio"]["utilizzabile_isotermo_g"]
    assert d["per_durata"]["mdot_g_s"] > 0


# --------------------------------------------------------------------------- #
# Coda dei lavori
# --------------------------------------------------------------------------- #
def _attendi(client, job_id, timeout=180.0):
    t0 = time.time()
    while time.time() - t0 < timeout:
        st = client.get(f"{API}/jobs/{job_id}").json()
        if st["status"] in ("done", "failed", "cancelled"):
            return st
        time.sleep(0.2)
    raise AssertionError("il lavoro non si e' concluso")


@pytest.mark.slow
def test_ciclo_di_vita_di_un_lavoro_di_geometria(client, base_design):
    r = client.post(f"{API}/jobs/geometry",
                    json={"design": base_design, "operating": OPERATING, "sector": False})
    assert r.status_code == 202
    job = r.json()
    assert job["status"] in ("pending", "running")

    st = _attendi(client, job["id"])
    assert st["status"] == "done", st.get("error")
    g = st["result"]
    assert g["watertight"] is True and g["valid_brep"] is True
    assert g["n_triangles"] > 1000
    assert g["mass_is_placeholder"] is True

    stl = client.get(g["stl_url"])
    assert stl.status_code == 200 and len(stl.content) > 10_000
    assert stl.headers["content-type"] == "model/stl"


def test_un_progetto_non_valido_e_rifiutato_subito_non_dopo(client, base_design):
    """Meglio un 422 adesso che un lavoro fallito fra dieci secondi."""
    incompleto = {k: v for k, v in base_design.items() if k != "p_c"}
    r = client.post(f"{API}/jobs/geometry",
                    json={"design": incompleto, "operating": OPERATING})
    assert r.status_code == 422


def test_sweep_richiede_un_seme(client):
    r = client.post(f"{API}/jobs/sweep", json={"n": 5})
    assert r.status_code == 422       # validazione di pydantic: seed obbligatorio


@pytest.mark.slow
def test_sweep_popola_il_database(client):
    r = client.post(f"{API}/jobs/sweep", json={"n": 12, "seed": 7, "operating": OPERATING})
    st = _attendi(client, r.json()["id"])
    assert st["status"] == "done", st.get("error")
    assert st["result"]["valutati"] + st["result"]["scartati"] == 12
    assert st["result"]["seed"] == 7


def test_annullare_un_lavoro_concluso_e_un_conflitto(client, base_design):
    r = client.post(f"{API}/jobs/sweep", json={"n": 1, "seed": 1, "operating": OPERATING})
    jid = r.json()["id"]
    _attendi(client, jid)
    assert client.delete(f"{API}/jobs/{jid}").status_code == 409
