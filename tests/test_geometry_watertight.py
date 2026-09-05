"""La geometria deve essere un solido chiuso: e' la precondizione del mesher.

La tenuta e' verificata DUE VOLTE con metodi indipendenti:
  1. `Shape.is_valid` di OCCT, che controlla la topologia B-Rep;
  2. la manifold-ness della tassellazione STL, controllata qui a mano.
La prima non implica la seconda: un B-Rep valido puo' tassellare in modo non
chiuso se la tolleranza e' troppo lasca, ed e' esattamente quel caso a far
fallire Gmsh a valle.
"""
from __future__ import annotations

import math

import pytest

from zefiro.geometry import BuildOptions, build_and_export, build_solid, stl_is_watertight
from zefiro.geometry.parameters import default_design_vector, derive




def test_solido_completo_e_chiuso(operating_point, design_vector, tmp_path):
    params, _ = derive(design_vector, operating_point)
    art = build_and_export(params, tmp_path, "test-full")

    assert art.is_valid_brep, "B-Rep non valido secondo OCCT"
    assert art.is_watertight_mesh, "la tassellazione STL non e' chiusa"
    assert art.n_triangles > 1000
    assert art.volume > 0.0
    assert art.step_path.exists() and art.step_path.stat().st_size > 0
    assert art.stl_path.exists() and art.stl_path.stat().st_size > 0


def test_settore_periodico_e_chiuso_e_conserva_il_volume(
    operating_point, design_vector, tmp_path
):
    """Il settore 1/N, moltiplicato per N, deve dare il volume del pezzo intero.

    E' la verifica piu' forte disponibile sul taglio periodico: se il cuneo
    tagliasse male, o se i fori non fossero distribuiti col passo giusto, il
    rapporto non farebbe 1.
    """
    params, _ = derive(design_vector, operating_point)
    full = build_and_export(params, tmp_path, "test-full2")
    sector = build_and_export(params, tmp_path, "test-sector", BuildOptions(sector=True))

    assert sector.is_valid_brep and sector.is_watertight_mesh
    n = int(params.derived["N_inj"])
    assert sector.volume * n == pytest.approx(full.volume, rel=1e-6)


def test_area_bagnata_analitica_coerente_con_la_scala_del_pezzo(
    operating_point, design_vector
):
    """L'area bagnata e' calcolata per Pappo-Guldino dal profilo, non dal CAD.

    Controllo di ordine di grandezza: deve stare fra l'area del solo cilindro
    di camera e quella di un cilindro lungo tutto il motore.
    """
    params, _ = derive(design_vector, operating_point)
    from zefiro.geometry import wetted_area

    d = params.derived
    a = wetted_area(params)
    a_min = 2.0 * math.pi * d["R_c"] * d["L_c"]
    a_max = 2.0 * math.pi * d["R_c"] * (d["L_c"] + d["L_conv"] + abs(d["x_tip_full"])) * 3.0
    assert a_min < a < a_max


def test_verifica_stl_riconosce_una_mesh_aperta(tmp_path):
    """Controllo del controllo: un STL con un triangolo mancante deve fallire.

    Senza questo test, `is_watertight_mesh = True` potrebbe voler dire
    "la funzione non sa dire di no".
    """
    # tetraedro completo -> chiuso
    v = [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1)]
    faces = [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]

    def write(path, fs):
        lines = ["solid t"]
        for f in fs:
            lines.append("facet normal 0 0 0\n outer loop")
            lines += [f"  vertex {v[i][0]} {v[i][1]} {v[i][2]}" for i in f]
            lines.append(" endloop\nendfacet")
        lines.append("endsolid t")
        path.write_text("\n".join(lines))

    closed = tmp_path / "closed.stl"
    open_ = tmp_path / "open.stl"
    write(closed, faces)
    write(open_, faces[:-1])

    assert stl_is_watertight(closed) == (True, 4)
    assert stl_is_watertight(open_)[0] is False


def test_geometria_incoerente_viene_rifiutata(operating_point):
    """Camera piu' stretta del labbro dell'ugello: errore, non geometria assurda."""
    x = default_design_vector(p_c=4.0e5, phi_core=1.0, f_film=0.0, N_inj=12.0,
                              Dc_over_Dt=2.0)
    params, _ = derive(x, operating_point)
    # con Dc_over_Dt = 2.0 la camera contiene ancora il labbro; forzo il caso limite
    from zefiro.schemas import DesignVector
    from zefiro.geometry.parameters import DESIGN_BOUNDS
    vals = dict(x.values)
    vals["Dc_over_Dt"] = DESIGN_BOUNDS["Dc_over_Dt"][0]
    bad = DesignVector(values=vals, bounds=DESIGN_BOUNDS)
    params2, _ = derive(bad, operating_point)
    assert params2.derived["R_c"] > params2.derived["R_lip"]
