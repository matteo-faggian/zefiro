"""Indice SQLite delle run + campi spaziali in Parquet.

Perche' SQLite con colonne SCALARI e non un blob JSON: le query che servono
davvero sono del tipo

    SELECT run_id FROM runs
    WHERE plug_trunc BETWEEN 0.2 AND 0.4 AND Isp_s > 120 AND feasible = 1

cioe' SQL su scalari indicizzabili. Una colonna JSON renderebbe questa query
lenta e non indicizzabile, e su 10^4-10^5 run e' la differenza fra un'analisi
interattiva e una che non fai.

Perche' i campi spaziali NON stanno qui: p e q a parete sono decine di migliaia
di righe per run, non vengono mai filtrati per valore, e servono solo al FEM.
Vanno in Parquet, referenziati per path.

Le run con working tree sporco finiscono in una tabella SEPARATA (`runs_dirty`).
Servono a esplorare, non a concludere, e non devono inquinare le query di
sintesi: la separazione e' strutturale, non una convenzione da ricordare.
"""
from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any, Iterable, Sequence

from zefiro.opt.objectives import CONSTRAINT_NAMES, OBJECTIVE_NAMES
from zefiro.schemas import (
    SCHEMA_VERSION,
    ContractViolation,
    DesignVector,
    GeometryArtifact,
    L0Result,
    Objectives,
)

#: Campi scalari di L0Result che diventano colonne. Esclusi di proposito:
#: X_eq (mappa), warnings e assumptions (tuple di stringhe). Quelli vivono nel
#: JSON accanto alla run: non li filtri mai per valore.
L0_SCALAR_FIELDS: tuple[str, ...] = (
    "phi_core", "phi_global", "AFR_stoich",
    "mdot_air", "mdot_fuel_core", "mdot_fuel_film",
    "T_ad", "gamma_c", "MW_c", "cp_c",
    "c_star", "M_e", "epsilon", "C_F", "Isp_s", "Isp_fuel_s", "thrust", "A_t",
    "tau_res", "tau_chem", "damkohler",
    "q_throat", "T_wall_adiabatic", "wall_volume",
)

GEOMETRY_SCALAR_FIELDS: tuple[str, ...] = (
    "volume", "wetted_area", "n_triangles", "mesh_tolerance",
)

#: Colonne su cui si crea un indice. Scelte perche' sono quelle su cui
#: filtrerai: i parametri di progetto e le prestazioni.
_INDEXED = ("thrust", "Isp_s", "T_ad", "feasible", "fidelity")


def _design_columns() -> tuple[str, ...]:
    from zefiro.geometry.parameters import DESIGN_BOUNDS
    return tuple(sorted(DESIGN_BOUNDS))


def _l0_only_columns() -> tuple[str, ...]:
    """Campi di L0Result che non sono GIA' colonne del vettore di progetto.

    `phi_core` compare in entrambi: e' una variabile libera e L0 lo restituisce
    invariato. Una sola colonna, e all'inserimento si VERIFICA che i due valori
    coincidano - se un giorno non coincidessero, ci sarebbe un bug in `derive`
    e questa e' l'unica sentinella che lo vedrebbe.
    """
    design = set(_design_columns())
    return tuple(f for f in L0_SCALAR_FIELDS if f not in design)


def schema_columns() -> list[tuple[str, str]]:
    """(nome, tipo SQL) di tutte le colonne della tabella `runs`."""
    cols: list[tuple[str, str]] = [
        ("run_id", "TEXT PRIMARY KEY"),
        ("schema_version", "TEXT NOT NULL"),
        ("code_rev", "TEXT"),
        ("env_hash", "TEXT"),
        ("created_utc", "TEXT NOT NULL"),
        ("fidelity", "TEXT NOT NULL"),
        ("feasible", "INTEGER"),
    ]
    cols += [(c, "REAL") for c in _design_columns()]
    cols += [(c, "REAL") for c in _l0_only_columns()]
    cols += [("geom_" + c, "REAL") for c in GEOMETRY_SCALAR_FIELDS]
    cols += [("geom_is_valid_brep", "INTEGER"), ("geom_is_watertight", "INTEGER")]
    cols += [("geom_sha256_step", "TEXT"), ("geom_sha256_stl", "TEXT")]
    cols += [("f_" + n, "REAL") for n in OBJECTIVE_NAMES]
    cols += [("g_" + n, "REAL") for n in CONSTRAINT_NAMES]

    # Sentinella. I prefissi "f_" e "g_" degli obiettivi vivono nello stesso
    # spazio dei nomi dei parametri di progetto, e uno di questi si chiama gia'
    # `f_film`: basterebbe un obiettivo chiamato "film" per avere due colonne
    # con lo stesso nome, che SQLite accetta in silenzio e rende una delle due
    # irraggiungibile. Meglio esplodere qui che scoprirlo su mille run.
    visti: dict[str, int] = {}
    for n, _t in cols:
        visti[n] = visti.get(n, 0) + 1
    doppie = sorted(n for n, k in visti.items() if k > 1)
    if doppie:
        raise ContractViolation(
            f"Nomi di colonna duplicati nello schema: {doppie}. "
            "Rinomina l'obiettivo, il vincolo o il parametro di progetto."
        )
    return cols


class RunStore:
    """Indice delle run. Usare come context manager."""

    def __init__(self, path: Path | str) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(str(self.path))
        self.conn.row_factory = sqlite3.Row
        self._create()

    def __enter__(self) -> "RunStore":
        return self

    def __exit__(self, *exc: Any) -> None:
        self.close()

    def close(self) -> None:
        self.conn.commit()
        self.conn.close()

    # -- schema ------------------------------------------------------------ #
    def _create(self) -> None:
        cols = schema_columns()
        ddl = ", ".join(f'"{n}" {t}' for n, t in cols)
        cur = self.conn.cursor()
        for table in ("runs", "runs_dirty"):
            cur.execute(f'CREATE TABLE IF NOT EXISTS "{table}" ({ddl})')
        self.migrate()
        for c in _INDEXED + _design_columns():
            cur.execute(f'CREATE INDEX IF NOT EXISTS "idx_runs_{c}" ON runs ("{c}")')
        self.conn.commit()

    def migrate(self) -> dict[str, list[str]]:
        """Allinea un database esistente allo schema corrente.

        Serve perche' `OBJECTIVE_NAMES` e `CONSTRAINT_NAMES` DEFINISCONO le
        colonne: aggiungere un obiettivo cambia lo schema, e un file scritto
        prima non ha quella colonna. `CREATE TABLE IF NOT EXISTS` non se ne
        accorge e le scritture successive fallirebbero a meta' popolamento.

        Si aggiungono solo colonne (ALTER TABLE ADD COLUMN), mai si tolgono:
        le righe vecchie restano leggibili e i campi nuovi valgono NULL, che e'
        l'unica risposta onesta ("quella run non l'ha misurato") - diverso da
        zero, che significherebbe "misurato, vale zero".

        Una colonna che esiste nel file ma NON nello schema corrente e' invece
        un segnale di ALLARME, non di routine: vuol dire che questo codice e'
        piu' vecchio del database, e scriverci sopra corromperebbe dati che non
        sa produrre. In quel caso si solleva.
        """
        attesi = {n: t for n, t in schema_columns()}
        aggiunte: dict[str, list[str]] = {}
        cur = self.conn.cursor()
        for table in ("runs", "runs_dirty"):
            presenti = {r["name"] for r in cur.execute(f'PRAGMA table_info("{table}")')}
            if not presenti:
                continue
            orfane = presenti - set(attesi)
            if orfane:
                raise ContractViolation(
                    f'La tabella "{table}" di {self.path.name} ha colonne che questo '
                    f"codice non conosce: {sorted(orfane)}. Il database e' stato scritto "
                    "da una versione PIU' RECENTE; aggiorna il codice invece di "
                    "scriverci sopra."
                )
            nuove = [c for c in attesi if c not in presenti]
            for c in nuove:
                tipo = attesi[c].replace(" PRIMARY KEY", "").replace(" NOT NULL", "")
                cur.execute(f'ALTER TABLE "{table}" ADD COLUMN "{c}" {tipo}')
            if nuove:
                aggiunte[table] = nuove
        self.conn.commit()
        return aggiunte

    # -- scrittura --------------------------------------------------------- #
    def insert(
        self,
        run_id: str,
        x: DesignVector,
        l0: L0Result,
        created_utc: str,
        code_rev: str,
        env_hash: str,
        geometry: GeometryArtifact | None = None,
        objectives: Objectives | None = None,
        fidelity: str = "L0",
        replace: bool = False,
    ) -> str:
        """Inserisce una run. Ritorna il nome della tabella usata.

        Il routing su `runs_dirty` e' automatico dal suffisso di `run_id`: non
        c'e' modo di sbagliarsi dimenticando un flag.
        """
        table = "runs_dirty" if run_id.endswith("-dirty") else "runs"
        row: dict[str, Any] = {
            "run_id": run_id,
            "schema_version": SCHEMA_VERSION,
            "code_rev": code_rev,
            "env_hash": env_hash,
            "created_utc": created_utc,
            "fidelity": fidelity,
            "feasible": None,
        }
        row.update({c: float(x.values[c]) for c in _design_columns()})
        row.update({c: _num(getattr(l0, c)) for c in _l0_only_columns()})

        # sentinella: L0 non deve alterare le variabili libere che riceve
        for shared in set(_design_columns()) & set(L0_SCALAR_FIELDS):
            if abs(getattr(l0, shared) - x.values[shared]) > 1e-12:
                raise ValueError(
                    f"Incoerenza su {shared!r}: design={x.values[shared]!r}, "
                    f"L0={getattr(l0, shared)!r}. derive() ha modificato una "
                    "variabile libera."
                )

        if geometry is not None:
            row.update({"geom_" + c: _num(getattr(geometry, c))
                        for c in GEOMETRY_SCALAR_FIELDS})
            row["geom_is_valid_brep"] = int(geometry.is_valid_brep)
            row["geom_is_watertight"] = int(geometry.is_watertight_mesh)
            row["geom_sha256_step"] = geometry.sha256_step
            row["geom_sha256_stl"] = geometry.sha256_stl

        if objectives is not None:
            row["feasible"] = int(objectives.feasible)
            row["fidelity"] = objectives.fidelity
            row.update({"f_" + k: float(v) for k, v in objectives.f.items()
                        if "f_" + k in dict(schema_columns())})
            row.update({"g_" + k: float(v) for k, v in objectives.g.items()
                        if "g_" + k in dict(schema_columns())})

        names = list(row)
        placeholders = ", ".join("?" for _ in names)
        cols = ", ".join(f'"{n}"' for n in names)
        verb = "INSERT OR REPLACE" if replace else "INSERT"
        self.conn.execute(
            f'{verb} INTO "{table}" ({cols}) VALUES ({placeholders})',
            [row[n] for n in names],
        )
        self.conn.commit()
        return table

    # -- lettura ----------------------------------------------------------- #
    def query(self, sql: str, params: Sequence[Any] = ()) -> list[sqlite3.Row]:
        return list(self.conn.execute(sql, params))

    def count(self, table: str = "runs") -> int:
        return int(self.conn.execute(f'SELECT COUNT(*) FROM "{table}"').fetchone()[0])

    def get(self, run_id: str) -> sqlite3.Row | None:
        table = "runs_dirty" if run_id.endswith("-dirty") else "runs"
        return self.conn.execute(
            f'SELECT * FROM "{table}" WHERE run_id = ?', (run_id,)
        ).fetchone()

    def pareto_candidates(self, limit: int = 1000) -> list[sqlite3.Row]:
        """Run fattibili, ordinate per spinta decrescente.

        Il fronte di Pareto vero lo calcola l'ottimizzatore: questa e' la query
        di sintesi che serve a guardare i risultati, e per costruzione **non
        vede** le run sporche.
        """
        return self.query(
            "SELECT * FROM runs WHERE feasible = 1 ORDER BY f_neg_thrust ASC LIMIT ?",
            (limit,),
        )


def _num(v: Any) -> float | None:
    return None if v is None else float(v)


# --------------------------------------------------------------------------- #
# Campi spaziali in Parquet (schema fisso, docs/architettura.md sezione 4.6)
# --------------------------------------------------------------------------- #
WALL_FIELD_COLUMNS: tuple[str, ...] = ("patch_id", "x", "y", "z", "value")


def write_wall_field(
    path: Path,
    patch_id: Iterable[str],
    x: Iterable[float],
    y: Iterable[float],
    z: Iterable[float],
    value: Iterable[float],
) -> Path:
    """Scrive un campo di parete con lo schema canonico (patch_id, x, y, z, value).

    E' l'unico formato che attraversa il confine CFD -> FEM. Lo schema e' fisso
    perche' i due solutori non parlano la stessa lingua e l'unico modo di
    verificare l'accoppiamento e' che il file in mezzo sia sempre lo stesso.
    """
    pa, pq = _pyarrow()
    table = pa.table({
        "patch_id": pa.array(list(patch_id), type=pa.string()),
        "x": pa.array(list(x), type=pa.float64()),
        "y": pa.array(list(y), type=pa.float64()),
        "z": pa.array(list(z), type=pa.float64()),
        "value": pa.array(list(value), type=pa.float64()),
    })
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, str(path), compression="zstd")
    return path


def read_wall_field(path: Path) -> dict[str, list[Any]]:
    pa, pq = _pyarrow()
    table = pq.read_table(str(path))
    if tuple(table.column_names) != WALL_FIELD_COLUMNS:
        raise ValueError(
            f"Schema Parquet inatteso in {path}: {table.column_names}, "
            f"attese {WALL_FIELD_COLUMNS}"
        )
    return {c: table.column(c).to_pylist() for c in table.column_names}


def _pyarrow():
    try:
        import pyarrow as pa
        import pyarrow.parquet as pq
    except ImportError as exc:  # pragma: no cover
        raise ImportError(
            "pyarrow non installato. Serve solo per i campi di parete di L1: "
            'installa con  pip install "zefiro[l1]"'
        ) from exc
    return pa, pq
