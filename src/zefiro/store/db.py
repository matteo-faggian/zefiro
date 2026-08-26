"""Indice SQLite delle run + campi in Parquet. STUB (fase 2).

Schema: una colonna scalare per ogni voce di x, di L0Result e di Objectives,
con `run_id` chiave primaria. Motivo: le query previste sono del tipo
    SELECT run_id FROM runs WHERE plug_trunc BETWEEN 0.2 AND 0.4 AND margin_yield > 2
cioe' SQL su scalari indicizzabili. Una colonna JSON renderebbe questa query
lenta e non indicizzabile.
"""
from __future__ import annotations

from pathlib import Path


def init_db(path: Path) -> None:
    raise NotImplementedError("Fase 2 della ROADMAP.")
