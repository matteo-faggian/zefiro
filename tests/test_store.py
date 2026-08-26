"""Database delle run: schema, routing delle run sporche, query di riferimento."""
from __future__ import annotations

import pytest

from zefiro.geometry.parameters import derive
from zefiro.opt.objectives import objectives_l0
from zefiro.store import RunStore, read_wall_field, write_wall_field
from zefiro.store.db import schema_columns

NOW = "2026-08-26T12:00:00+00:00"


@pytest.fixture
def store(tmp_path):
    with RunStore(tmp_path / "runs.db") as s:
        yield s


@pytest.fixture
def evaluated(operating_point, design_vector):
    params, l0 = derive(design_vector, operating_point)
    obj = objectives_l0(
        "abc", l0, operating_point, design_vector.values["p_c"],
        min_feature_size=0.4e-3, derived=params.derived,
    )
    return design_vector, l0, obj


def test_nessuna_colonna_duplicata():
    """`phi_core` sta sia nel vettore di progetto sia in L0Result: una sola
    colonna, altrimenti SQLite rifiuta la CREATE TABLE."""
    names = [n for n, _ in schema_columns()]
    assert len(names) == len(set(names)), \
        [n for n in names if names.count(n) > 1]


def test_round_trip_esatto_sui_float(store, evaluated):
    x, l0, obj = evaluated
    store.insert("run0001", x, l0, NOW, "rev1", "envhash", objectives=obj)
    row = store.get("run0001")
    assert row is not None
    assert row["p_c"] == x.values["p_c"]          # REAL = float64: esatto
    assert row["thrust"] == l0.thrust
    assert row["Isp_s"] == l0.Isp_s
    assert row["c_star"] == l0.c_star


def test_le_run_sporche_finiscono_in_una_tabella_separata(store, evaluated):
    """Le run con working tree sporco servono a esplorare, non a concludere.
    La separazione e' STRUTTURALE: dipende dal run_id, non da un flag che
    qualcuno puo' dimenticare."""
    x, l0, obj = evaluated
    assert store.insert("pulita", x, l0, NOW, "v1", "e", objectives=obj) == "runs"
    assert store.insert("abc-dirty", x, l0, NOW, "v1-dirty", "e", objectives=obj) \
        == "runs_dirty"
    assert store.count("runs") == 1
    assert store.count("runs_dirty") == 1


def test_le_run_sporche_non_compaiono_nelle_query_di_sintesi(store, evaluated):
    x, l0, obj = evaluated
    store.insert("sporca-dirty", x, l0, NOW, "v-dirty", "e", objectives=obj)
    assert store.pareto_candidates() == []


def test_query_di_riferimento_della_roadmap(store, evaluated):
    """La query che la fase 2 deve saper rispondere."""
    x, l0, obj = evaluated
    store.insert("r1", x, l0, NOW, "v1", "e", objectives=obj)
    rows = store.query(
        "SELECT run_id FROM runs WHERE plug_trunc BETWEEN ? AND ? AND Isp_s > ?",
        (0.2, 1.0, 50.0),
    )
    assert [r["run_id"] for r in rows] == ["r1"]


def test_sentinella_di_coerenza_fra_design_e_l0(store, evaluated):
    """Se `derive` alterasse una variabile libera, l'inserimento deve fallire.

    E' l'unico punto del sistema che se ne accorgerebbe: senza questa
    sentinella, il database conterrebbe un `phi_core` che non corrisponde a
    quello con cui il risultato e' stato calcolato.
    """
    import dataclasses

    x, l0, obj = evaluated
    corrotto = dataclasses.replace(l0, phi_core=l0.phi_core + 0.01)
    with pytest.raises(ValueError, match="phi_core"):
        store.insert("r2", x, corrotto, NOW, "v1", "e", objectives=obj)


def test_fattibilita_registrata(store, evaluated):
    x, l0, obj = evaluated
    store.insert("r3", x, l0, NOW, "v1", "e", objectives=obj)
    assert store.get("r3")["feasible"] == int(obj.feasible)


# --------------------------------------------------------------------------- #
# Campi di parete in Parquet
# --------------------------------------------------------------------------- #
def test_round_trip_campo_di_parete(tmp_path):
    path = write_wall_field(
        tmp_path / "q_wall.parquet",
        patch_id=["wall_chamber", "wall_throat"],
        x=[0.0, 1.5], y=[0.1, 0.2], z=[0.3, 0.4], value=[1.0e6, 2.5e6],
    )
    back = read_wall_field(path)
    assert back["patch_id"] == ["wall_chamber", "wall_throat"]
    assert back["value"] == [1.0e6, 2.5e6]
    assert back["x"] == [0.0, 1.5]


def test_schema_parquet_sbagliato_viene_rifiutato(tmp_path):
    """E' l'unico formato che attraversa il confine CFD -> FEM: se lo schema
    cambia in silenzio, il FEM legge la colonna sbagliata."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    bad = tmp_path / "bad.parquet"
    pq.write_table(pa.table({"a": [1.0], "b": [2.0]}), str(bad))
    with pytest.raises(ValueError, match="Schema Parquet inatteso"):
        read_wall_field(bad)
