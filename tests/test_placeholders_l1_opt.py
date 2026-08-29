"""Gli stub di L1 e del surrogato devono FALLIRE, non restituire numeri finti.

Un risultato inventato a quel livello sarebbe indistinguibile da uno vero, e
finirebbe nel database delle run.

`driver.run_nsga2` e `mesh.generate_mesh` non sono piu' in questa lista: sono
implementati, e i loro test stanno in test_driver.py e test_mesh.py.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from zefiro.l1 import cfd, fem, mapping
from zefiro.opt import surrogate


@pytest.mark.parametrize("fn,args", [
    (cfd.run_reacting_foam, (None, None, None, Path("."), 4, 1000)),
    (fem.run_calculix, (None, None, {}, Path("."), True, 5.0)),
    (mapping.map_wall_loads, (Path("a"), Path("b"), Path("c"))),
    (surrogate.fit_gp, (None, None)),
])
def test_gli_stub_sollevano_not_implemented(fn, args):
    with pytest.raises(NotImplementedError):
        fn(*args)
