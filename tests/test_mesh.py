"""La mesh del settore periodico.

I test qui non chiedono "gira?" ma "e' la mesh giusta?". La differenza conta:
una mesh con le condizioni al contorno sulle facce sbagliate, o con le facce
periodiche non conformi, gira benissimo e da' un risultato falso.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from zefiro.geometry.parameters import default_design_vector, derive
from zefiro.geometry.profile import fluid_polygon, revolved_volume
from zefiro.l1.mesh import (
    MAX_SKEWNESS,
    REQUIRED_BOUNDARIES,
    MeshError,
    MeshOptions,
    generate_mesh,
)
from zefiro.schemas import CANONICAL_BOUNDARIES

GROSSA = MeshOptions(injector_size=5.0e-4, wall_size=1.5e-3,
                     farfield_size=2.0e-2, refine_distance=2.0e-3)


@pytest.fixture(scope="module")
def punto(request):
    from zefiro.schemas import FuelSpec, OperatingPoint
    op = OperatingPoint(
        p_amb=101325.0, p_air_supply=6.0e5, p_fuel_supply=8.0e5,
        fuel=FuelSpec(composition={"C3H8": 1.0}, phase_at_injection="gas",
                      thermo_source="gri30.yaml"),
        T_air_in=293.0, T_fuel_in=283.0, mdot_air_max=0.0718,
        cd_injector_ox=0.75, cd_injector_fuel=0.75,
    )
    # N_inj = 12 e d_ox_ratio = 0.05 e non il ginocchio del fronte: il
    # ginocchio (18 fori, d_ox_ratio 0.072) NON e' fabbricabile, e lo si e'
    # scoperto proprio qui. Vedi MIN_INJECTOR_LAND_FRACTION.
    x = default_design_vector(
        p_c=5.21e5, phi_core=0.84, f_film=0.02, N_inj=12.0, Dc_over_Dt=2.44,
        Lc_over_Dc=1.58, plug_trunc=0.785, t_wall=0.0009, d_ox_ratio=0.05,
        fuel_vel_ratio=1.06, conv_half_angle=0.678, theta_swirl=0.47,
    )
    return derive(x, op)


@pytest.fixture(scope="module")
def mesh(punto, tmp_path_factory):
    params, _ = punto
    return generate_mesh(params, tmp_path_factory.mktemp("mesh"), "test", GROSSA)


@pytest.mark.slow
def test_frontiere_obbligatorie_e_nessun_nome_inventato(mesh):
    assert REQUIRED_BOUNDARIES <= set(mesh.boundary_names)
    assert set(mesh.boundary_names) <= set(CANONICAL_BOUNDARIES)


@pytest.mark.slow
def test_le_facce_periodiche_corrispondono(mesh):
    """Il controllo decisivo: `setPeriodic` accetta senza protestare una
    trasformazione sbagliata, e il difetto si vede solo dentro OpenFOAM come
    squilibrio di massa. Qui e' verificato nodo per nodo."""
    assert mesh.is_periodic
    assert mesh.periodic_mismatch < 1.0e-5


@pytest.mark.slow
def test_qualita_entro_i_criteri_di_openfoam(mesh):
    assert mesh.max_skewness < MAX_SKEWNESS
    # la non-ortogonalita' puo' superare 70 gradi, ma deve restare LOCALE:
    # se fosse diffusa non sarebbe piu' una questione di correttori
    assert mesh.frac_non_orthogonal < 0.005, (
        f"{mesh.frac_non_orthogonal:.2%} delle facce oltre 70 gradi: "
        "il problema non e' piu' locale"
    )


@pytest.mark.slow
def test_il_volume_della_mesh_e_quello_analitico(mesh, punto):
    """Verifica INDIPENDENTE dal mesher: il volume dei tetraedri deve dare il
    volume del settore, che si calcola in forma chiusa dal poligono meridiano.
    Se il dominio fosse costruito male - un pezzo mancante, un booleano andato
    storto - questo test lo vedrebbe e nessun altro."""
    import gmsh

    params, _ = punto
    pts, _, _ = fluid_polygon(params.derived, params.plug_contour_x,
                              params.plug_contour_r)
    sector = 2.0 * math.pi / int(round(params.derived["N_inj"]))
    atteso = revolved_volume(pts) * sector / (2.0 * math.pi)

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.merge(str(mesh.msh_path))
        tipi, _, nodi_el = gmsh.model.mesh.getElements(3)
        conn = np.array(nodi_el[list(tipi).index(4)], dtype=np.int64).reshape(-1, 4)
        tag, coord, _ = gmsh.model.mesh.getNodes()
        xyz = np.array(coord, dtype=float).reshape(-1, 3)
        idx = np.zeros(int(np.max(tag)) + 1, dtype=np.int64)
        idx[np.array(tag, dtype=np.int64)] = np.arange(len(tag))
        c = xyz[idx[conn]]
    finally:
        gmsh.finalize()

    v = np.abs(np.einsum("ij,ij->i", c[:, 1] - c[:, 0],
                         np.cross(c[:, 2] - c[:, 0], c[:, 3] - c[:, 0]))) / 6.0
    ottenuto = float(v.sum())
    # la discrepanza residua e' la corda contro l'arco: la mesh e' fatta di
    # facce piatte e il settore e' curvo, quindi la mesh sta SOTTO
    assert ottenuto == pytest.approx(atteso, rel=0.02), (ottenuto, atteso)
    assert ottenuto <= atteso * 1.001, "la mesh non puo' contenere piu' del dominio"


@pytest.mark.slow
def test_gli_iniettori_ci_sono_e_hanno_l_area_giusta(mesh, punto):
    """I fori sono la ragione per cui questa mesh esiste: se l'imprint fosse
    andato storto, le patch esisterebbero comunque ma con l'area sbagliata, e
    la portata imposta darebbe la velocita' sbagliata."""
    import gmsh

    params, _ = punto
    d = params.derived

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.merge(str(mesh.msh_path))
        aree = {}
        for dim, tag in gmsh.model.getPhysicalGroups(2):
            nome = gmsh.model.getPhysicalName(dim, tag)
            tot = 0.0
            for sup in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
                tipi, _, nodi_el = gmsh.model.mesh.getElements(2, sup)
                if 2 not in tipi:
                    continue
                conn = np.array(nodi_el[list(tipi).index(2)],
                                dtype=np.int64).reshape(-1, 3)
                tag_n, coord, _ = gmsh.model.mesh.getNodes()
                xyz = np.array(coord, dtype=float).reshape(-1, 3)
                idx = np.zeros(int(np.max(tag_n)) + 1, dtype=np.int64)
                idx[np.array(tag_n, dtype=np.int64)] = np.arange(len(tag_n))
                t = xyz[idx[conn]]
                tot += float(np.linalg.norm(
                    np.cross(t[:, 1] - t[:, 0], t[:, 2] - t[:, 0]), axis=1).sum() / 2.0)
            aree[nome] = tot
    finally:
        gmsh.finalize()

    # un foro per settore: l'area della patch e' l'area di UN foro
    atteso_ox = math.pi / 4.0 * d["d_ox"] ** 2
    atteso_f = math.pi / 4.0 * d["d_fuel"] ** 2
    # 3 % e non 10 %: con MIN_SEGMENTS_PER_HOLE = 16 l'errore del poligono
    # inscritto vale 1.3 %, e il resto e' la curvatura della faccia anulare.
    # A portata imposta, un errore sull'area e' un errore sulla VELOCITA' di
    # iniezione, quindi sul rapporto delle quantita' di moto, quindi sulla
    # miscelazione: la grandezza per cui questa mesh esiste.
    assert aree["inlet_air"] == pytest.approx(atteso_ox, rel=0.03), (
        aree["inlet_air"], atteso_ox)
    assert aree["inlet_fuel_core"] == pytest.approx(atteso_f, rel=0.03), (
        aree["inlet_fuel_core"], atteso_f)


@pytest.mark.slow
def test_fori_troppo_grandi_per_il_settore_sollevano(punto):
    """Con pochi iniettori molto larghi i fori non ci stanno nel settore. Deve
    dirlo, non produrre una mesh con i fori tagliati dal piano periodico."""
    params, _ = punto
    d = dict(params.derived)
    d["d_ox"] = 20.0 * d["d_ox"]
    import dataclasses
    rotto = dataclasses.replace(params, derived=d)
    with pytest.raises(MeshError, match="non ci stanno"):
        generate_mesh(rotto, __import__("pathlib").Path("/tmp/zmesh"), "rotto", GROSSA)
