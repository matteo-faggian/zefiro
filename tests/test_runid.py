"""L'identita' di una run deve essere deterministica e sensibile agli input."""
from __future__ import annotations

from zefiro.geometry.parameters import DESIGN_BOUNDS, default_design_vector
from zefiro.schemas import DesignVector
from zefiro.store import env_fingerprint, make_run_id


def test_deterministico(operating_point, design_vector):
    a = make_run_id(design_vector, operating_point, code_rev="abc123")
    b = make_run_id(design_vector, operating_point, code_rev="abc123")
    assert a == b and len(a) == 32


def test_sensibile_a_una_perturbazione_minima(operating_point, design_vector):
    """Un bit di differenza sul vettore di progetto deve cambiare l'id."""
    vals = dict(design_vector.values)
    vals["p_c"] = vals["p_c"] * (1.0 + 1e-15)
    other = DesignVector(values=vals, bounds=DESIGN_BOUNDS)
    assert make_run_id(design_vector, operating_point, "r") != make_run_id(other, operating_point, "r")


def test_sensibile_alla_revisione_del_codice(operating_point, design_vector):
    assert make_run_id(design_vector, operating_point, "r1") != \
           make_run_id(design_vector, operating_point, "r2")


def test_il_working_tree_sporco_e_marcato(operating_point, design_vector):
    """Le run sporche servono a esplorare, non a concludere: devono essere
    riconoscibili senza consultare il database."""
    rid = make_run_id(design_vector, operating_point, code_rev="v0.1.0-dirty")
    assert rid.endswith("-dirty")


def test_env_fingerprint_contiene_le_versioni_che_contano():
    fp = env_fingerprint()
    for k in ("python", "platform", "cantera", "code_rev"):
        assert k in fp and fp[k]
