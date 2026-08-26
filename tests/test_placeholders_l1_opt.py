"""Gli stub di L1 e opt devono FALLIRE, non restituire numeri finti.

Un risultato inventato a quel livello sarebbe indistinguibile da uno vero, e
finirebbe nel database delle run.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from zefiro.l1 import cfd, fem, mapping, mesh
from zefiro.opt import driver, surrogate


@pytest.mark.parametrize("fn,args", [
    (mesh.generate_mesh, (None, Path("."), 30.0, 5, 1e-3)),
    (cfd.run_reacting_foam, (None, None, None, Path("."), 4, 1000)),
    (fem.run_calculix, (None, None, {}, Path("."), True, 5.0)),
    (mapping.map_wall_loads, (Path("a"), Path("b"), Path("c"))),
    (surrogate.fit_gp, (None, None)),
    (driver.run_nsga2, (None, 40, 20, 0)),
])
def test_gli_stub_sollevano_not_implemented(fn, args):
    with pytest.raises(NotImplementedError):
        fn(*args)


def test_objectives_l0_gia_utilizzabile(operating_point, design_vector):
    """Il blocco obiettivi a L0 e' invece gia' operativo: serve a scartare
    candidati prima di spendere una run L1."""
    from zefiro.geometry.parameters import derive
    from zefiro.opt.objectives import objectives_l0

    _, l0 = derive(design_vector, operating_point)
    o = objectives_l0("test", l0)
    assert o.fidelity == "L0"
    assert o.f["neg_thrust"] == -l0.thrust
    assert o.f["neg_isp"] == -l0.Isp_s
