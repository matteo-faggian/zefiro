"""Criterio di chiusura della fase 2: stesso input -> stesso output, byte per byte.

Non e' un test di comodo. E' l'affermazione centrale del progetto
("deterministico e riproducibile") resa FALSIFICABILE: se questo test passa,
l'affermazione e' verificata sulla catena geometria; se fallisce, sappiamo
esattamente dove si rompe.
"""
from __future__ import annotations

import pytest

from zefiro.geometry import build_and_export
from zefiro.geometry.parameters import derive
from zefiro.store import make_run_id


def test_due_export_della_stessa_geometria_sono_byte_identici(
    operating_point, design_vector, tmp_path
):
    params, _ = derive(design_vector, operating_point)
    a = build_and_export(params, tmp_path / "a", "run-fisso")
    b = build_and_export(params, tmp_path / "b", "run-fisso")

    assert a.sha256_step == b.sha256_step, "STEP non riproducibile"
    assert a.sha256_stl == b.sha256_stl, "STL non riproducibile"
    assert a.volume == b.volume
    assert a.n_triangles == b.n_triangles


def test_lheader_step_e_normalizzato_e_porta_il_run_id(
    operating_point, design_vector, tmp_path
):
    """OCCT scrive l'ora di creazione nell'header, che rende due export della
    stessa geometria byte-diversi. La si sostituisce con un sentinella fisso e
    col run_id, cosi' il file si autoidentifica e resta riproducibile.
    L'ora vera non va persa: sta in `created_utc` nel database."""
    params, _ = derive(design_vector, operating_point)
    art = build_and_export(params, tmp_path, "run-xyz")
    header = art.step_path.read_text(errors="surrogateescape").split("ENDSEC;")[0]
    assert "zefiro run-xyz" in header
    assert "1970-01-01T00:00:00" in header
    # la normalizzazione tocca SOLO l'header: il file resta un STEP valido
    text = art.step_path.read_text(errors="surrogateescape")
    assert text.startswith("ISO-10303-21;")
    assert text.rstrip().endswith("END-ISO-10303-21;")
    assert "DATA;" in text and "APPLICATION_PROTOCOL_DEFINITION" in text


def test_il_run_id_non_dipende_dallistante_di_esecuzione(operating_point, design_vector):
    """Se dipendesse dal tempo, rilanciare la stessa run creerebbe un duplicato
    invece di riconoscerla."""
    a = make_run_id(design_vector, operating_point, code_rev="v1")
    b = make_run_id(design_vector, operating_point, code_rev="v1")
    assert a == b


def test_l0_e_deterministico(operating_point, design_vector):
    """Cantera risolve l'equilibrio numericamente: verifichiamo che due
    chiamate identiche diano bit identici, non solo valori vicini."""
    _, first = derive(design_vector, operating_point)
    _, second = derive(design_vector, operating_point)
    assert first.T_ad == second.T_ad
    assert first.c_star == second.c_star
    assert first.thrust == second.thrust
    assert first.X_eq == second.X_eq


@pytest.mark.parametrize("n", [8, 16])
def test_il_doe_e_riproducibile_end_to_end(operating_point, n):
    """Stesso seme -> stessi vettori di progetto -> stessi risultati L0."""
    from zefiro.doe import latin_hypercube
    from zefiro.geometry.parameters import (
        DESIGN_BOUNDS, INTEGER_PARAMETERS, bounds_per_impianto,
    )
    from zefiro.schemas import DesignVector

    # La scatola generica contiene p_c fino a 9 bar, che e' cio' che
    # l'architettura sa fare; questo impianto ne regge molte meno. Campionare
    # nella scatola generica genera punti in cui il combustibile non puo'
    # entrare in camera - ed e' giusto che evaluate_l0 li rifiuti, non che il
    # bound li nasconda.
    bounds = bounds_per_impianto(operating_point)

    def run(seed):
        out = []
        for v in latin_hypercube(n, bounds, seed, INTEGER_PARAMETERS):
            _, l0 = derive(DesignVector(values=v, bounds=bounds), operating_point)
            out.append(l0.thrust)
        return out

    assert run(3) == run(3)
