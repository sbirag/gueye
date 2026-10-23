"""Validations essentielles (version rapide des scripts de validation/) : noyau, 1D/2D/3D, convergence, PML, temps."""
import warnings
import numpy as np
import pytest
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import gueye as gy
from gueye.kernel import TridiagInverse, SingularOperator

warnings.simplefilter("ignore")
rel = lambda a, b: np.abs(a - b).max() / np.abs(b).max()


@pytest.mark.parametrize("N", [4, 5, 12, 40])
@pytest.mark.parametrize("coef", [(-1., 2., -1.), (-1., 1.3, -1.), (-1., 4., -1.), (-2., 5., 0.5), (-1. + .5j, 2. + .3j, -1. - .2j)])
@pytest.mark.parametrize("ends", [None, "N", "R"])
def test_noyau_vs_inverse_dense(N, coef, ends):
    l, d, u = coef; kw = {}
    if ends == "N": kw = dict(u1=2 * u, lN=2 * l)
    if ends == "R": kw = dict(d1=d + 0.7 * l, u1=u + l, dN=d - 0.3 * u, lN=l + u)
    if ends == "N" and (l, d, u) == (-1., 2., -1.):                                   # Neumann pur : singulier, doit etre signale
        with pytest.raises(SingularOperator): TridiagInverse(N, l, d, u, **kw)
        return
    T = TridiagInverse(N, l, d, u, **kw)
    assert rel(T.inverse(), np.linalg.inv(T.dense())) < 1e-11


def test_noyau_racine_double_et_grand_N():
    T = TridiagInverse(400000, -1.0, 2.0 - 0.25**2, -1.0)
    b = np.zeros(400000); b[200000] = 1
    x = T.solve(b)
    assert abs(T.entries(200007, 200000) - x[200007]) < 1e-9 and np.isfinite(T.column(200000)).all()


def test_noyau_singularites():
    with pytest.raises(SingularOperator): TridiagInverse(20, -1., 2., -1., u1=-2., lN=-2.)     # Neumann pur
    with pytest.raises(SingularOperator): TridiagInverse(21, -1., 0., -1.)                      # d = 0, N impair


_BC = {"DD": ("D", "D"), "DN": ("D", "N"), "ND": ("N", "D"), "NN": ("N", "N"), "RR": (("R", 1.7), ("R", 1.7)),
       "RN": (("R", 1.7), "N"), "NR": ("N", ("R", 1.7)), "DR": ("D", ("R", 1.7)), "RD": (("R", 1.7), "D")}


@pytest.mark.parametrize("name", list(_BC))
@pytest.mark.parametrize("k", [0.7, 3.3, 0.9 + 0.3j])
def test_1d_formule_vs_bande_vs_dense(name, k):
    pb = gy.problem((2.0,), (60,), _BC[name], k=k)
    P = pb.backend; f = lambda x: np.cos(2 * x) + 0.3
    rhs = P.rhs(f(pb.axes[0]))
    ug = P.K.solve(rhs, "green"); ub = P.K.solve(rhs, "banded"); ud = np.linalg.solve(P.dense(), rhs)
    assert rel(ug, ud) < 1e-9 and rel(ub, ud) < 1e-9


def _manufactured_error(name, M):
    L = 2.0; alpha = 1.7
    uex = lambda x: np.cos(2.3 * x) + 0.4 * x**2 - 0.3 * x
    du = lambda x: -2.3 * np.sin(2.3 * x) + 0.8 * x - 0.3
    d2 = lambda x: -2.3**2 * np.cos(2.3 * x) + 0.8
    k = 1.3; c1 = 0.0
    lt, rt = _BC[name]
    def spec(b, side):
        kind = b if isinstance(b, str) else b[0]
        s = -1 if side == "L" else 1; xb = 0.0 if side == "L" else L
        if kind == "D": return ("D", uex(xb))
        if kind == "N": return ("N", s * du(xb))
        return ("R", alpha, s * du(xb) + alpha * uex(xb))
    P = gy.Problem1D(L, M, spec(lt, "L"), spec(rt, "R"), c0=-k**2)
    return np.abs(P.solve(-d2(P.x) - k**2 * uex(P.x)) - uex(P.x)).max()


@pytest.mark.parametrize("name", list(_BC))
def test_1d_ordre_2(name):
    e1, e2 = _manufactured_error(name, 80), _manufactured_error(name, 160)
    assert 1.9 < np.log2(e1 / e2) < 2.1


@pytest.mark.parametrize("combo", [("D", "N", "D", "N"), ("N", "D", ("R", 1.3), "D"), (("R", 1.3), "D", "N", "N"), ("N", "N", "N", "N")])
@pytest.mark.parametrize("k", [1.1, 1.1 + 0.4j])
def test_2d_vs_creux(combo, k):
    pb = gy.problem((2.0, 1.4), (30, 24), [combo[:2], combo[2:]], k=k)
    P = pb.backend; f = lambda x, y: np.sin(x) * y + 0.2
    F = P.rhs(pb._f(f)); u = P.solve_reduced(F)
    ref = spla.spsolve(P.dense().tocsc(), F.ravel().astype(complex)).reshape(F.shape)
    assert rel(u, ref) < 1e-10


@pytest.mark.parametrize("combo", [("D", "N", "D", "N", "D", "N"), ("N", "D", ("R", 1.3), "D", "N", "D"), (("R", 1.3),) * 6])
def test_3d_vs_creux(combo):
    pb = gy.problem((2.0, 1.5, 1.2), (9, 8, 7), [combo[:2], combo[2:4], combo[4:]], k=1.1 + 0.2j)
    P = pb.backend; F = P.rhs(pb._f(lambda x, y, z: np.sin(x) + y * z))
    ref = spla.spsolve(P.dense().tocsc(), F.ravel().astype(complex)).reshape(F.shape)
    assert rel(P.solve_reduced(F), ref) < 1e-10


def test_2d_ordre_2_cl_inhomogenes():
    Lx, Ly, al = 2.0, 1.4, 1.6
    uex = lambda x, y: np.cos(1.7 * x + 0.3) * np.sin(2.1 * y + 0.2) + 0.3 * x * y + 0.2 * y**2
    ux = lambda x, y: -1.7 * np.sin(1.7 * x + 0.3) * np.sin(2.1 * y + 0.2) + 0.3 * y
    uy = lambda x, y: 2.1 * np.cos(1.7 * x + 0.3) * np.cos(2.1 * y + 0.2) + 0.3 * x + 0.4 * y
    lap = lambda x, y: -(1.7**2 + 2.1**2) * np.cos(1.7 * x + 0.3) * np.sin(2.1 * y + 0.2) + 0.4
    k = 1.1
    def err(M):
        bc = dict(left=("N", lambda y: -ux(0 * y, y)), right=("R", al, lambda y: ux(Lx + 0 * y, y) + al * uex(Lx + 0 * y, y)),
                  bottom=("D", lambda x: uex(x, 0 * x)), top=("N", lambda x: uy(x, Ly + 0 * x)))
        P = gy.Problem2D(Lx, Ly, M, int(M * Ly / Lx), c0=-k**2, **bc)
        Y, X = np.meshgrid(P.y, P.x, indexing="ij")
        return np.abs(P.solve(-lap(X, Y) - k**2 * uex(X, Y)) - uex(X, Y)).max()
    assert 1.9 < np.log2(err(40) / err(80)) < 2.1


def test_newmark_modal_vs_lu_creux():
    W = gy.wave((2.0, 1.4), (30, 22), [("D", ("R", 1.5)), ("N", "D")], c=1.0, eta=0.05)
    P = W.backend.P; dt, n = 0.02, 60
    sig = lambda t: (1 - 2 * (3 * (t - 0.8))**2) * np.exp(-(3 * (t - 0.8))**2)
    src = W.point_source((0.9, 0.7), sig)
    t, U = W.march(dt, n, source=src)[:2]
    A = P.dense().tocsc(); a0 = 1 + 0.05 * 0.5 * dt; s = dt**2 / 4
    lu = spla.splu((a0 * sp.identity(A.shape[0]) + s * A).tocsc())
    sl = (slice(P.Py.i0, P.Py.i1 + 1), slice(P.Px.i0, P.Px.i1 + 1))
    fv = lambda tt: src(tt)[sl].ravel()
    u = np.zeros(A.shape[0]); v = u.copy(); a = fv(0) - A @ u
    for i in range(n):
        ut = u + dt * v + dt**2 * 0.25 * a; vt = v + dt * 0.5 * a
        a = lu.solve(fv((i + 1) * dt) - A @ ut - 0.05 * vt); u = ut + s * a; v = vt + 0.5 * dt * a
    assert rel(U[-1][sl].ravel(), u) < 1e-11


def test_pml_1d_elimination_exacte_et_proche_de_transparente():
    k, L, M, m = 5.0, 2.0, 80, 12
    P = gy.Problem1D(L, M, ("P", m, 1e-6, 2), ("P", m, 1e-6, 2), c0=-k**2)
    T, nl, nr = P.extended_dense()
    b = np.zeros(T.shape[0], complex); b[nl + 25] = 1 / P.h
    uref = np.linalg.solve(T, b)[nl:nl + P.N]
    assert rel(P.green(25)[P.i0:P.i1 + 1], uref) < 1e-11
    M2, L2, k2 = 400, 4.0, 0.2 * 400 / 4.0
    ga = gy.Problem1D(L2, M2, ("A",), ("A",), c0=-k2**2).green(150)
    gp = gy.Problem1D(L2, M2, ("P", 32, 1e-6), ("P", 32, 1e-6), c0=-k2**2).green(150)
    assert rel(gp, ga) < 1e-4


def test_pml_2d_axe_ferme_vs_transparente():
    k = 14.0
    ga = gy.Problem2D(3.0, 3.0, 60, 60, ("N", 0.), ("N", 0.), ("A",), ("A",), c0=-k**2).green(30, 20)
    gp = gy.Problem2D(3.0, 3.0, 60, 60, ("N", 0.), ("N", 0.), ("P", 16, 1e-8), ("P", 16, 1e-8), c0=-k**2).green(30, 20)
    assert rel(gp, ga) < 2e-2


# ---------------------------------------------------------------- 0.6.0 : solveur modal vectorise, bc_time, temporel ouvert
def test_modestack_2d_egal_lu():
    from gueye.problem2d import Problem2D
    rng = np.random.default_rng(3)
    P = Problem2D(1, 1, 30, 25, left=("N", 0.), right=("R", 2., 0.), top=("N", 0.))
    f = rng.standard_normal(P.shape)
    a, b = P.solve(f), P.solve(f, method="lu")
    assert np.abs(a - b).max() < 1e-12 * np.abs(b).max()


def test_modestack_3d_egal_lu():
    from gueye.problem3d import Problem3D
    rng = np.random.default_rng(4)
    P = Problem3D(1, 1, 1, 12, 9, 10, left=("N", 0.), top=("R", 2., 0.))
    f = rng.standard_normal(P.shape)
    a, b = P.solve(f), P.solve(f, method="lu")
    assert np.abs(a - b).max() < 1e-12 * np.abs(b).max()


def test_bc_time_onde_plane_converge():
    from gueye.problem2d import Wave2D
    errs = []
    for n in (40, 80):
        W = Wave2D(1, 1, n, 4, left=("D", 1.0), right=("D", 0.0), bottom=("N", 0.), top=("N", 0.),
                   c=1.0, bc_time={"left": lambda t: np.sin(6 * t)})
        r = W.march(0.4 / n, n, every=n); t, x = r[0][-1], W.x
        m = x < t - 0.15
        errs.append(np.abs(r[1][-1][2] - np.sin(6 * (t - x)))[m].max())
    assert errs[1] < 0.6 * errs[0] and errs[1] < 5e-3


def test_bc_time_face_inconnue():
    from gueye.problem2d import Wave2D
    with pytest.raises(ValueError):
        Wave2D(1, 1, 10, 10, bc_time={"gauche": lambda t: t})


def test_time_response_ferme_vs_newmark():
    import gueye as gy
    from gueye.recipes import time_response
    f0 = 3.0
    ric = lambda t: (1 - 2 * (np.pi * f0 * (t - 0.4))**2) * np.exp(-(np.pi * f0 * (t - 0.4))**2)
    L, cells, src, rec = (1.0, 0.8), (40, 32), (0.3, 0.4), [(0.7, 0.5)]
    t = np.linspace(0, 2.5, 126)
    pb = gy.problem(size=L, cells=cells, k=1.0, bc=[("D", "D"), ("D", "D")])
    U = time_response(pb, src, rec, ric, t, wmax=80.0, nw=2400, shift=0.3)
    W = gy.wave(size=L, cells=cells, bc=[("D", "D"), ("D", "D")], c=1.0)
    tt, _, R = W.march(0.0025, 1000, source=W.point_source(src, ric), receivers=rec)
    assert np.abs(U[0] - np.interp(t, tt, np.array(R)[:, 0])).max() < 1e-3


def test_time_response_ouvert_pml_vs_grand_domaine():
    import gueye as gy
    from gueye.recipes import time_response
    ric = lambda t: (1 - 2 * (np.pi * (t - 1.0))**2) * np.exp(-(np.pi * (t - 1.0))**2)
    t = np.linspace(0, 6, 61); P = ("P", 16, 1e-6)
    pb = gy.problem(size=(2.0, 2.0), cells=(40, 40), k=1.0, bc=[(P, P), (P, P)])
    U = time_response(pb, (1.0, 1.0), [(1.8, 1.2)], ric, t, wmax=14.0, nw=140, shift=0.5)
    big = gy.problem(size=(8.0, 8.0), cells=(160, 160), k=1.0, bc=[("D", "D"), ("D", "D")])
    Ub = time_response(big, (4.0, 4.0), [(4.8, 4.2)], ric, t, wmax=14.0, nw=140, shift=0.5)
    assert np.abs(U - Ub).max() < 5e-3 * np.abs(Ub).max() + 1e-4


def test_time_response_fenetre_trop_courte():
    import gueye as gy
    from gueye.recipes import time_response
    pb = gy.problem(size=(1.0, 1.0), cells=(10, 10), k=1.0, bc=[("D", "D"), ("D", "D")])
    with pytest.raises(ValueError):
        time_response(pb, (0.3, 0.3), [(0.6, 0.6)], lambda t: t * 0, np.array([0.0, 10.0]), wmax=10.0, nw=5)
