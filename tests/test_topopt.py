"""Ottimizzazione topologica: il FEM, i filtri, le derivate, l'ottimizzatore.

Ordine dei test = ordine della fiducia. Un errore nel FEM non produce un
risultato palesemente sbagliato: produce una struttura ottimizzata alla
perfezione per un problema che non esiste. Per questo si verifica dal basso:
prima che l'elemento finito sia esatto, poi che le derivate siano esatte, poi
che l'ottimizzatore trovi la soluzione nota.
"""
from __future__ import annotations

import numpy as np
import pytest

from zefiro.topopt import filters as flt
from zefiro.topopt.fem2d import (
    AXISYMMETRIC,
    PLANE_STRAIN,
    PLANE_STRESS,
    Grid,
    assemble,
    constitutive,
    solve,
)
from zefiro.topopt.problems import mbb_beam
from zefiro.topopt.simp import (
    Material,
    Options,
    Problem,
    compliance_and_sensitivity,
    optimize,
)
from zefiro.topopt.stress import build_B_set, stress_objective, thermal_free_strain


# --------------------------------------------------------------------------- #
# 1. L'elemento finito
# --------------------------------------------------------------------------- #
@pytest.mark.parametrize("mode", [PLANE_STRESS, PLANE_STRAIN])
def test_patch_test(mode):
    """Un campo di deformazione COSTANTE deve essere riprodotto esattamente.

    E' il controllo standard di un elemento finito. Se non lo passa, l'elemento
    non converge e ogni risultato costruito sopra e' privo di valore.
    """
    g = Grid(nx=4, ny=3, dx=0.7, dy=1.3, mode=mode)
    D = constitutive(210e9, 0.3, mode)
    K = assemble(g, g.element_stiffness_set(D), np.ones(g.n_elem), g.element_dofs())

    a, b = 1.0e-4, -3.0e-5
    coords = np.array([[i * g.dx, j * g.dy]
                       for i in range(g.nx + 1) for j in range(g.ny + 1)])
    uex = np.empty(g.n_dof)
    uex[0::2] = a * coords[:, 0]
    uex[1::2] = b * coords[:, 1]
    border = [n for n, (x, y) in enumerate(coords)
              if min(abs(x), abs(x - g.nx * g.dx), abs(y), abs(y - g.ny * g.dy)) < 1e-12]
    fixed = np.array(sorted([2 * n for n in border] + [2 * n + 1 for n in border]))
    f = np.asarray(-K[:, fixed] @ uex[fixed]).ravel()
    u = uex.copy()
    free = np.setdiff1d(np.arange(g.n_dof), fixed)
    import scipy.sparse.linalg as sla
    u[free] = sla.spsolve(K[free, :][:, free].tocsc(), f[free])
    assert np.max(np.abs(u - uex)) / np.max(np.abs(uex)) < 1e-12


def test_cilindro_in_pressione_coincide_con_lame():
    """Verifica dell'elemento ASSIALSIMMETRICO contro la soluzione esatta.

    E' il caso di Zefiro: la camera e' un solido di rivoluzione, e la rigidezza
    circonferenziale (eps_theta = u/r) e' proprio il termine che distingue
    l'assialsimmetrico da un problema piano.
    """
    E, nu, p = 200e9, 0.3, 5.0e6
    ri, ro, H = 0.050, 0.070, 0.004
    nx, ny = 40, 4
    g = Grid(nx=nx, ny=ny, dx=(ro - ri) / nx, dy=H / ny, mode=AXISYMMETRIC, x0=ri)
    D = constitutive(E, nu, AXISYMMETRIC)
    K = assemble(g, g.element_stiffness_set(D), np.ones(g.n_elem), g.element_dofs())
    f = np.zeros(g.n_dof)
    for j in range(ny):
        for node, w in ((g.node(0, j), 0.5), (g.node(0, j + 1), 0.5)):
            f[2 * node] += p * 2 * np.pi * ri * g.dy * w
    fixed = np.array([2 * g.node(i, 0) + 1 for i in range(nx + 1)]
                     + [2 * g.node(i, ny) + 1 for i in range(nx + 1)])
    u = solve(K, f, fixed)
    k = ro / ri
    ur_exact = p * ri / E * ((1 + nu) * ((1 - 2 * nu) + k ** 2) / (k ** 2 - 1))
    assert u[2 * g.node(0, ny // 2)] == pytest.approx(ur_exact, rel=1e-3)


# --------------------------------------------------------------------------- #
# 2. I filtri
# --------------------------------------------------------------------------- #
def test_il_filtro_di_densita_elimina_la_scacchiera():
    """La scacchiera e' un artefatto numerico: il FEM la valuta rigida quando
    non lo e'. Senza filtro l'ottimizzazione la produce sistematicamente."""
    nx = ny = 20
    H = flt.density_filter_matrix(nx, ny, rmin=1.6)
    cb = (np.indices((nx, ny)).sum(axis=0) % 2).reshape(-1).astype(float)
    out = H @ cb
    assert out.var() < cb.var() / 100.0
    assert out.mean() == pytest.approx(cb.mean(), abs=1e-9)   # conserva la massa


def test_derivata_di_heaviside():
    x = np.linspace(0.0, 1.0, 17)
    d = flt.heaviside_derivative(x, 8.0)
    fd = (flt.heaviside(x + 1e-6, 8.0) - flt.heaviside(x - 1e-6, 8.0)) / 2e-6
    assert np.max(np.abs(d - fd)) < 1e-6


def test_il_filtro_di_stampabilita_elimina_le_isole_sospese():
    """Un blocco di materiale che non poggia su nulla non si stampa: in SLM
    servirebbe un supporto, e dentro un condotto non si toglie piu'."""
    nx = ny = 20
    x = np.zeros((nx, ny))
    x[:, 0] = 1.0
    x[8:12, 10:14] = 1.0                    # isola sospesa
    xi = flt.printable(x.reshape(-1), nx, ny).reshape(nx, ny)
    assert xi[8:12, 10:14].max() < 1e-3


def test_il_filtro_di_stampabilita_conserva_una_rampa_a_45_gradi():
    """45 gradi e' il limite: sotto quello il materiale si sostiene da solo e
    il filtro non deve toglierlo."""
    nx = ny = 20
    y = np.zeros((nx, ny))
    y[:, 0] = 1.0
    for j in range(1, ny):
        y[max(0, 10 - j):min(nx, 10 + j + 1), j] = 1.0
    yi = flt.printable(y.reshape(-1), nx, ny)
    assert yi.sum() / y.sum() > 0.95
    assert flt.max_overhang_violation(yi, nx, ny) == 0.0


def test_adjoint_del_filtro_di_stampabilita():
    """Il filtro e' ricorsivo verso l'alto: la sua trasposta jacobiana scorre
    verso il basso. Se e' sbagliata, l'ottimizzatore converge tranquillamente
    verso qualcosa che non e' un ottimo, senza dare segnali."""
    rng = np.random.default_rng(0)
    nx, ny = 12, 10
    x = rng.uniform(0.05, 1.0, nx * ny)
    w = rng.normal(size=nx * ny)
    g = flt.printable_chain(x, nx, ny, w)
    h = 1e-6
    fd = np.empty_like(x)
    for k in range(x.size):
        xp = x.copy(); xp[k] += h
        xm = x.copy(); xm[k] -= h
        fd[k] = (w @ flt.printable(xp, nx, ny) - w @ flt.printable(xm, nx, ny)) / (2 * h)
    assert np.max(np.abs(g - fd)) / np.max(np.abs(fd)) < 1e-6


# --------------------------------------------------------------------------- #
# 3. Le derivate dell'obiettivo
# --------------------------------------------------------------------------- #
def _cantilever(mode=PLANE_STRESS, alpha=0.0):
    g = Grid(nx=8, ny=6, dx=1.0, dy=1.0, mode=mode)
    mat = Material(E=1.0, nu=0.3, alpha=alpha)
    fixed = np.array([2 * g.node(0, j) for j in range(g.ny + 1)]
                     + [2 * g.node(0, j) + 1 for j in range(g.ny + 1)])
    f = np.zeros(g.n_dof)
    f[2 * g.node(g.nx, 0) + 1] = -1.0
    return g, mat, fixed, f


@pytest.mark.parametrize("dT", [0.0, 300.0])
def test_sensibilita_della_compliance(dT):
    """Con carico termico la compliance NON e' autoaggiunta: il carico dipende
    dal progetto e la derivata acquista il termine 2 u' df/dx. Dimenticarlo
    e' l'errore classico, e questo test lo intercetta."""
    g, mat, fixed, f = _cantilever(alpha=1e-5 if dT else 0.0)
    prob = Problem(grid=g, material=mat, fixed_dofs=fixed, force=f, delta_T=dT)
    opt = Options(penal=3.0, thermal_penal=1.0, printable=False)
    D = constitutive(mat.E, mat.nu, g.mode)
    ke, ed = g.element_stiffness_set(D), g.element_dofs()
    fth = g.element_thermal_set(D, mat.alpha, dT) if dT else None
    rng = np.random.default_rng(1)
    xi = rng.uniform(0.3, 0.9, g.n_elem)

    _, dc, _ = compliance_and_sensitivity(prob, xi, opt, ke, fth, ed)
    h = 1e-7
    fd = np.empty_like(xi)
    for k in range(g.n_elem):
        a = xi.copy(); a[k] += h
        b = xi.copy(); b[k] -= h
        fd[k] = (compliance_and_sensitivity(prob, a, opt, ke, fth, ed)[0]
                 - compliance_and_sensitivity(prob, b, opt, ke, fth, ed)[0]) / (2 * h)
    assert np.max(np.abs(dc - fd)) / np.max(np.abs(fd)) < 1e-5


@pytest.mark.parametrize("mode", [PLANE_STRESS, AXISYMMETRIC])
def test_sensibilita_della_tensione(mode):
    """Derivata per variabile aggiunta della norma-P di von Mises, con carico
    termico. E' il pezzo piu' delicato del modulo."""
    g = Grid(nx=7, ny=6, dx=1e-3, dy=1e-3, mode=mode,
             x0=0.02 if mode == AXISYMMETRIC else 0.0)
    mat = Material(E=195e9, nu=0.3, alpha=17e-6)
    fixed = np.array([2 * g.node(0, j) for j in range(g.ny + 1)]
                     + [2 * g.node(0, j) + 1 for j in range(g.ny + 1)]
                     + [2 * g.node(g.nx, j) + 1 for j in range(g.ny + 1)])
    f = np.zeros(g.n_dof)
    f[2 * g.node(g.nx, g.ny // 2)] = 900.0
    prob = Problem(grid=g, material=mat, fixed_dofs=fixed, force=f, delta_T=300.0)
    opt = Options(penal=3.0, thermal_penal=1.0, printable=False)
    D = constitutive(mat.E, mat.nu, mode)
    ke, ed, B = g.element_stiffness_set(D), g.element_dofs(), build_B_set(g)
    fth = g.element_thermal_set(D, mat.alpha, 300.0)
    eps0 = thermal_free_strain(g, mat.alpha, 300.0)
    rng = np.random.default_rng(3)
    xi = rng.uniform(0.35, 0.95, g.n_elem)

    _, dJ, _ = stress_objective(prob, xi, opt, ke, fth, ed, D, B, 200e6, eps0)
    h = 1e-8
    fd = np.empty_like(xi)
    for k in range(g.n_elem):
        a = xi.copy(); a[k] += h
        b = xi.copy(); b[k] -= h
        fd[k] = (stress_objective(prob, a, opt, ke, fth, ed, D, B, 200e6, eps0)[0]
                 - stress_objective(prob, b, opt, ke, fth, ed, D, B, 200e6, eps0)[0]) / (2 * h)
    assert np.max(np.abs(dJ - fd)) / np.max(np.abs(fd)) < 1e-5


def test_la_tensione_esclude_la_dilatazione_libera():
    """sigma = D (B u - eps_termica). Un corpo LIBERO di dilatarsi non e'
    sollecitato: se il codice usasse B u e basta, gli attribuirebbe una tensione
    che non prova."""
    g = Grid(nx=6, ny=6, dx=1e-3, dy=1e-3, mode=PLANE_STRESS)
    mat = Material(E=195e9, nu=0.3, alpha=17e-6)
    # vincoli minimi: solo i moti rigidi
    fixed = np.array([2 * g.node(0, 0), 2 * g.node(0, 0) + 1, 2 * g.node(0, g.ny) ])
    prob = Problem(grid=g, material=mat, fixed_dofs=fixed,
                   force=np.zeros(g.n_dof), delta_T=250.0)
    opt = Options(penal=1.0, thermal_penal=1.0, printable=False)
    D = constitutive(mat.E, mat.nu, PLANE_STRESS)
    ke, ed, B = g.element_stiffness_set(D), g.element_dofs(), build_B_set(g)
    fth = g.element_thermal_set(D, mat.alpha, 250.0)
    eps0 = thermal_free_strain(g, mat.alpha, 250.0)
    _, _, vm = stress_objective(prob, np.ones(g.n_elem), opt, ke, fth, ed,
                                D, B, 200e6, eps0)
    libero = D[0, 0] * mat.alpha * 250.0        # ordine di grandezza se sbagliassimo
    assert vm.max() < 0.02 * libero


# --------------------------------------------------------------------------- #
# 4. L'ottimizzatore
# --------------------------------------------------------------------------- #
@pytest.mark.slow
def test_mbb_riproduce_il_traliccio_noto():
    """La trave MBB e' il caso di riferimento: la soluzione e' un traliccio con
    correnti superiore e inferiore continui e diagonali fra i due. Se non esce
    quella, l'ottimizzatore ha un problema."""
    p = mbb_beam(nx=90, ny=30)
    r = optimize(p, Options(volume_fraction=0.5, rmin=2.0, max_iter=60,
                            printable=False, beta=0.0, beta_max=0.0))
    assert r.compliance < 0.35 * r.history[0]["objective"]
    assert r.xi.mean() == pytest.approx(0.5, abs=0.02)

    X = r.xi.reshape(p.grid.nx, p.grid.ny)
    # I correnti devono essere continui e pieni sulla parte CARICATA della luce.
    # Verso l'appoggio di destra il corrente superiore si esaurisce, ed e'
    # corretto: li' il momento flettente e' nullo e quel materiale non serve.
    utile = int(0.7 * X.shape[0])
    assert X[:utile, 0].mean() > 0.9, "manca il corrente inferiore"
    assert X[:utile, -1].mean() > 0.9, "manca il corrente superiore"
    # la parte centrale deve essere in gran parte vuota: e' un traliccio, non
    # una piastra
    assert X[:, X.shape[1] // 2].mean() < 0.5


@pytest.mark.slow
def test_la_discesa_e_monotona_sulla_compliance():
    p = mbb_beam(nx=60, ny=20)
    r = optimize(p, Options(volume_fraction=0.5, rmin=2.0, max_iter=40,
                            printable=False, beta=0.0, beta_max=0.0))
    c = [h["objective"] for h in r.history]
    peggioramenti = sum(1 for a, b in zip(c[:-1], c[1:]) if b > a * 1.02)
    assert peggioramenti <= 2, f"{peggioramenti} peggioramenti: OC non sta convergendo"


@pytest.mark.slow
def test_il_vincolo_di_volume_vale_sulla_densita_fisica():
    """I filtri cambiano il volume: imporre il vincolo prima di attraversarli
    consegnerebbe un pezzo con una frazione di materiale diversa da quella
    richiesta. E' un errore che ho gia' commesso una volta."""
    p = mbb_beam(nx=60, ny=20)
    for vf in (0.35, 0.55):
        r = optimize(p, Options(volume_fraction=vf, rmin=2.4, max_iter=35,
                                printable=True, beta=1.0, beta_max=2.0))
        assert r.xi.mean() == pytest.approx(vf, abs=0.03), f"vf={vf}"


@pytest.mark.slow
def test_il_risultato_e_stampabile():
    p = mbb_beam(nx=60, ny=20)
    r = optimize(p, Options(volume_fraction=0.45, rmin=2.4, max_iter=35,
                            printable=True, beta=1.0, beta_max=4.0))
    assert r.overhang_violation < 0.02


def test_lobiettivo_di_tensione_richiede_lammissibile():
    p = mbb_beam(nx=20, ny=10)
    with pytest.raises(ValueError, match="sigma_allow"):
        optimize(p, Options(objective="stress", max_iter=1))


def test_obiettivo_non_riconosciuto():
    p = mbb_beam(nx=20, ny=10)
    with pytest.raises(ValueError, match="obiettivo"):
        optimize(p, Options(objective="fantasia", max_iter=1))


def test_il_solutore_generale_rispetta_il_vincolo_su_un_caso_noto():
    """Lagrangiano aumentato su un problema a soluzione nota: il minimo di
    SOMMA (x-0.9)^4 con somma vincolata e' la distribuzione uniforme."""
    from zefiro.topopt.nlp import NlpOptions, solve_augmented_lagrangian
    n, target = 40, 12.0
    obj = lambda x: (float(np.sum((x - 0.9) ** 4)), 4 * (x - 0.9) ** 3)
    vol = lambda x: (float(x.sum()), np.ones_like(x))
    x, _ = solve_augmented_lagrangian(np.full(n, 0.5), obj, vol, target,
                                      NlpOptions(outer=14, inner_maxiter=60))
    assert x.sum() == pytest.approx(target, rel=2e-3)
    assert x.std() < 1e-6


@pytest.mark.slow
def test_la_tensione_usa_il_solutore_generale_e_non_diverge():
    """Forzare OC sulla tensione fa divergere: la sensibilita' cambia segno.
    Con 'auto' il modulo sceglie il Lagrangiano aumentato, che regge."""
    p = mbb_beam(nx=40, ny=16, E=195e9)
    p.material.nu = 0.3
    r = optimize(p, Options(objective="stress", sigma_allow=200e6,
                            volume_fraction=0.4, rmin=2.0, printable=False,
                            beta=0.0, beta_max=0.0, outer_per_beta=5,
                            inner_maxiter=20))
    assert np.isfinite(r.compliance)
    assert r.xi.mean() == pytest.approx(0.4, abs=0.03)
