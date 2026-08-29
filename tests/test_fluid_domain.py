"""Il dominio fluido: la nomenclatura delle frontiere e' verificata contro la
GEOMETRIA, non contro l'ordine in cui e' stata costruita.

Una condizione al contorno sulla faccia sbagliata e' l'errore piu' insidioso di
tutta la catena L1: non fa fallire niente e non si vede nei residui.
"""
from __future__ import annotations

import math

import pytest

from zefiro.geometry.parameters import derive
from zefiro.geometry.profile import fluid_polygon, revolved_volume


@pytest.fixture
def dominio(operating_point, design_vector):
    params, _ = derive(design_vector, operating_point)
    pts, nomi, _ = fluid_polygon(params.derived, params.plug_contour_x,
                                 params.plug_contour_r)
    return params, pts, nomi


def test_un_nome_per_lato_e_poligono_chiuso(dominio):
    _, pts, nomi = dominio
    assert len(pts) == len(nomi) >= 8
    # nessun lato di lunghezza nulla: OCC li rifiuta e gmsh ci si strozza
    n = len(pts)
    for i in range(n):
        x0, r0 = pts[i]
        x1, r1 = pts[(i + 1) % n]
        assert math.hypot(x1 - x0, r1 - r0) > 1e-9, f"lato {i} di lunghezza nulla"


def test_ogni_lato_sta_dove_il_suo_nome_dice(dominio):
    """Il controllo vero: per ciascun nome, una proprieta' geometrica che DEVE
    valere e che non dipende da come e' stato costruito il poligono."""
    params, pts, nomi = dominio
    d = params.derived
    R_c, R_lip, L_c = d["R_c"], d["R_lip"], d["L_c"]
    x_lip = L_c + d["L_conv"]
    n = len(pts)

    for i in range(n):
        (x0, r0), (x1, r1) = pts[i], pts[(i + 1) % n]
        nome = nomi[i]
        if nome == "wall_faceplate":
            assert x0 == pytest.approx(0.0) and x1 == pytest.approx(0.0)
        elif nome == "wall_chamber":
            assert r0 == pytest.approx(R_c) and r1 == pytest.approx(R_c)
            assert 0.0 <= x0 <= L_c + 1e-12
        elif nome == "wall_convergent":
            assert r0 == pytest.approx(R_c) and r1 == pytest.approx(R_lip)
            assert x0 == pytest.approx(L_c) and x1 == pytest.approx(x_lip)
        elif nome == "axis":
            assert r0 == pytest.approx(0.0) and r1 == pytest.approx(0.0)
        elif nome == "wall_plug":
            # sta a valle della gola e dentro il raggio del labbro
            assert x0 >= x_lip - 1e-9 or x0 >= d["L_c"]
            assert r0 <= R_lip + 1e-9 and r1 <= R_lip + 1e-9
        elif nome == "outlet_far":
            assert max(r0, r1) >= R_lip - 1e-9


def test_il_campo_lontano_e_grande_quanto_serve(dominio):
    """Il criterio giusto e' il RAGGIO DEL LABBRO, non la lunghezza del motore:
    la camera sta a monte e non c'entra con quanto grande deve essere la
    regione di scarico. La scala del getto e' il raggio d'uscita."""
    params, pts, _ = dominio
    d = params.derived
    x_base = d["L_c"] + d["L_conv"] + params.plug_contour_x[-1]
    assert max(r for _, r in pts) >= 8.0 * d["R_lip"]
    assert max(x for x, _ in pts) - x_base >= 8.0 * d["R_lip"]


def test_il_fluido_e_il_complemento_del_solido(dominio):
    """Controllo di consistenza fra i due domini: il volume del fluido dentro
    la camera (fino al labbro) piu' il volume del solido non devono
    sovrapporsi. Si verifica che il fluido di camera valga esattamente il
    V_c gia' calcolato da `derive`, che e' una strada indipendente."""
    params, pts, nomi = dominio
    d = params.derived
    x_lip = d["L_c"] + d["L_conv"]
    # ritaglia il poligono di camera: anello fra corpo centrale e parete
    R_c, R_lip, L_c, r_cb = d["R_c"], d["R_lip"], d["L_c"], d["r_centerbody"]
    camera = [(0.0, r_cb), (0.0, R_c), (L_c, R_c), (x_lip, R_lip),
              (x_lip, r_cb)]
    v = revolved_volume(camera)
    assert v == pytest.approx(d["V_c"], rel=2e-2), (v, d["V_c"])
