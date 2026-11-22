import warnings
import numpy as np
import pytest
import scipy.sparse.linalg as spla
import gueye as gy
from gueye import recipes

warnings.simplefilter("ignore")
rel = lambda a, b: np.abs(a - b).max() / np.abs(b).max()


def test_dimension_et_noeuds():
    p1 = gy.problem((2.0,), (40,), ("D", "N"), k=1.0)
    p2 = gy.problem((2.0, 1.0), (40, 20), [("D", "N"), ("N", "D")], k=1.0)
    p3 = gy.problem((2.0, 1.0, 1.0), (20, 10, 10), [("D", "D")] * 3, k=1.0)
    assert (p1.dim, p2.dim, p3.dim) == (1, 2, 3)
    assert p2.node((1.0, 0.5)) == (10, 20)            # (iy, ix)
    assert p3.node((1.0, 0.5, 0.2)) == (2, 5, 10)     # (iz, iy, ix)


def test_solve_contre_creux_2d_inhomogene():
    pb = gy.problem((2.0, 1.4), (30, 22), [(("D", lambda y: 1 + y), ("R", 1.3, lambda y: 0.2 * y)), ("N", ("D", 0.5))], k=1.2 + 0.2j)
    u = pb.solve(lambda x, y: np.sin(x) * y)
    P = pb.backend
    F = P.rhs(pb._f(lambda x, y: np.sin(x) * y))
    ref = spla.spsolve(P.dense().tocsc(), F.ravel().astype(complex)).reshape(F.shape)
    assert rel(u[P.Py.i0:P.Py.i1 + 1, P.Px.i0:P.Px.i1 + 1], ref) < 1e-11


def test_green_et_value_coherents():
    pb = gy.problem((2.0, 1.4), (30, 22), [("D", ("R", 1.0)), ("N", "D")], k=2.0)
    g = pb.green((0.8, 0.6)); r = pb.node((1.3, 0.9))
    assert abs(pb.value((1.3, 0.9), (0.8, 0.6)) - g[r]) < 1e-12 * abs(g).max()


def test_sweep_equals_pointwise():
    pb = gy.problem((2.0, 1.4), (30, 22), [("D", "N"), ("N", "D")])
    ks = np.array([1.0, 2.3, 3.1 + 0.2j])
    G = pb.sweep(ks, (0.7, 0.5), [(1.5, 1.0)])
    for a, k in enumerate(ks):
        s = gy.problem((2.0, 1.4), (30, 22), [("D", "N"), ("N", "D")], k=k)
        assert abs(G[a, 0] - s.value((1.5, 1.0), (0.7, 0.5))) < 1e-10 * abs(G).max()


@pytest.mark.parametrize("bc,cells,size", [([("D", "N"), ("N", ("R", 1.3))], (14, 11), (2.0, 1.5)),
                                           ([("D", "D"), ("N", "N"), ("D", "N")], (8, 7, 6), (2.0, 1.5, 1.2))])
def test_resonances_vs_valeurs_propres_creuses(bc, cells, size):
    pb = gy.problem(size, cells, bc)
    kr = pb.resonances(5)
    import scipy.sparse as sp
    P = pb.backend; A = P.dense().toarray() if hasattr(P.dense(), "toarray") else P.dense()
    ev = np.sort(np.linalg.eigvals(A).real)[:5]
    assert np.allclose(kr, np.sqrt(np.maximum(ev, 0)), atol=1e-6)


def test_resonance_exacte_levee_par_l_api():
    k1 = gy.problem((1.0,), (40,), ("D", "D")).resonances(1)[0]
    with pytest.raises(gy.SingularOperator):
        gy.problem((1.0,), (40,), ("D", "D"), k=k1)
    with pytest.raises(gy.SingularOperator):
        gy.problem((1.0, 1.0), (20, 20), [("N", "N"), ("N", "N")], k=0.0)          # Neumann pur a k = 0
    with pytest.raises(gy.SingularOperator):
        gy.problem((1.0,), (40,), ("D", "D")).with_k(k1)
    pb = gy.problem((1.0,), (40,), ("D", "D"))
    G = pb.sweep(np.array([1.0, k1, 2.0]), (0.37,), [(0.6,)])
    assert np.isnan(G[1, 0]) and np.isfinite(G[[0, 2], 0]).all()                     # NaN a la resonance, valeurs ailleurs


def test_points_de_bord_et_hors_domaine():
    p2 = gy.problem((1.0, 1.0), (20, 20), [("D", "D"), ("D", "D")], k=3.0)
    for pt in [(0.0, 0.5), (0.5, 0.0), (1.0, 0.5)]:
        with pytest.raises(ValueError): p2.green(pt)                                 # source sur noeud de Dirichlet
    # Même un point hors domaine de moins d'une demi-maille doit être rejeté.
    for pt in [(1.7, 0.5), (-0.2, 0.5), (-0.01, 0.5), (1.01, 0.5), (0.5, np.nan)]:
        with pytest.raises(ValueError): p2.node(pt)                                  # hors domaine physique
    W = gy.wave((1.0,), (20,), ("N", "N"))
    for pt in [(-0.01,), (1.01,)]:
        with pytest.raises(ValueError): W.node(pt)
    assert p2.value((0.0, 0.5), (0.5, 0.5)) == 0.0                                  # recepteur sur Dirichlet : Green homogene nulle
    p3 = gy.problem((1.0, 1.0, 1.0), (8, 8, 8), [("D", "D")] * 3, k=2.0)
    with pytest.raises(ValueError): p3.green((0.5, 0.5, 0.0))
    p1 = gy.problem((1.0,), (20,), ("D", "N"), k=1.0)
    with pytest.raises(ValueError): p1.green((0.0,))
    v = gy.problem((1.0, 1.0), (20, 20), [(("D", 2.0), ("D", 3.0)), ("N", "N")], k=1.0).at([(0.0, 0.5), (1.0, 0.5)], 0.0)
    assert np.allclose(v, [2.0, 3.0])                                                # `at` sur face de Dirichlet : valeur imposee


def test_resonances_robin_complexe_ou_negatif():
    import warnings as w
    with w.catch_warnings(record=True) as rec:
        w.simplefilter("always")
        pb = gy.problem((1.0,), (60,), ("D", ("R", 1j)))
        kr = pb.resonances(2)
    assert np.iscomplexobj(kr) and np.abs(kr.imag).max() > 0.1 and rec                # pas de fausses resonances reelles
    ev = np.linalg.eigvals(gy.Problem1D(1.0, 60, ("D", 0.), ("R", 1j, 0.), check=False).dense())
    ev = ev[np.argsort(ev.real)][:2]
    assert np.allclose(kr, np.sqrt(ev + 0j))
    pn = gy.problem((1.0,), (60,), (("R", -5.0), "D"))                              # alpha < 0 : valeur propre negative
    with w.catch_warnings(record=True):
        w.simplefilter("always")
        assert pn.resonances(1)[0].imag > 0.1


@pytest.mark.parametrize("size,cells", [((1.0,), (0,)), ((-1.0,), (10,)), ((1.0,), (1,)), ((0.0,), (10,)), ((1.0,), (10.5,)),
                                        ((1.0, 1.0), (10,)), ((np.inf,), (10,))])
def test_entrees_invalides(size, cells):
    with pytest.raises(ValueError):
        gy.problem(size, cells, ("D", "D"))


def test_plus_petite_grille():
    u = gy.problem((1.0,), (2,), ("D", "D"), k=1.0).solve(1.0)
    assert np.isfinite(u).all() and u[0] == 0 and u[-1] == 0


def test_recepteur_sur_bord_en_temps():
    W = gy.wave((1.0, 1.0), (20, 20), [(("D", 0.0), "N"), ("N", "D")], c=1.0)
    src = W.point_source((0.5, 0.5), lambda t: np.exp(-((t - 0.3) / 0.1)**2))
    t, U, R = W.march(0.02, 30, source=src, receivers=[(0.0, 0.5), (0.5, 0.0), (0.5, 1.0)])
    assert np.allclose(R[:, 0], 0.0) and np.allclose(R[:, 2], 0.0) and np.abs(R[:, 1]).max() > 0     # D -> 0, N -> signal
    with pytest.raises(ValueError):
        W.point_source((0.0, 0.5), lambda t: 1.0)


def _explicit_green_2d(P, src):
    A, Nye, Nxe, oy, ox = P.explicit_dense()
    b = np.zeros((Nye, Nxe), complex); b[src[0] - P.Py.i0 + oy, src[1] - P.Px.i0 + ox] = 1.0
    u = spla.spsolve(A, b.ravel()).reshape(Nye, Nxe)[oy:oy + P.Ny, ox:ox + P.Nx]
    return u, P.green(*src)[P.Py.i0:P.Py.i1 + 1, P.Px.i0:P.Px.i1 + 1] * (P.hx * P.hy)


@pytest.mark.parametrize("k", [9.0, 9.0 + 0.5j])
def test_pml_toutes_faces_2d_vs_assemblage_explicite(k):
    P = gy.Problem2D(2.0, 1.5, 40, 30, ("P", 6, 1e-6), ("P", 6, 1e-6), ("P", 6, 1e-6), ("P", 6, 1e-6), c0=-k**2)
    u, g = _explicit_green_2d(P, (14, 18))
    assert np.abs(u - g).max() < 1e-11 * np.abs(u).max()


def test_pml_toutes_faces_3d_vs_assemblage_explicite():
    P = gy.Problem3D(1.2, 1.0, 0.9, 12, 10, 9, *[("P", 4, 1e-4)] * 6, c0=-36.0)
    A, (Nz, Ny, Nx), (oz, oy, ox) = P.explicit_dense()
    s = (4, 5, 6); b = np.zeros((Nz, Ny, Nx), complex); b[s[0] - P.Pz.i0 + oz, s[1] - P.Py.i0 + oy, s[2] - P.Px.i0 + ox] = 1
    u = spla.spsolve(A, b.ravel()).reshape(Nz, Ny, Nx)[oz:oz + P.Nz, oy:oy + P.Ny, ox:ox + P.Nx]
    g = P.green(*s)[P.Pz.i0:P.Pz.i1 + 1, P.Py.i0:P.Py.i1 + 1, P.Px.i0:P.Px.i1 + 1] * (P.hx * P.hy * P.hz)
    assert np.abs(u - g).max() < 1e-11 * np.abs(u).max()


def test_domaine_ouvert_toutes_faces_vs_grand_domaine():
    # meme maille h = 0.05 : petit domaine (PML partout) contre domaine 3x plus grand (PML partout), region commune
    k = 10.0; pml = 16
    _, gs = recipes.open_field((3.0, 3.0), (60, 60), k, (1.5, 1.5), pml=pml, R0=1e-8, all_open=True)
    _, gl = recipes.open_field((9.0, 9.0), (180, 180), k, (4.5, 4.5), pml=pml, R0=1e-8, all_open=True)
    assert rel(gs, gl[60:121, 60:121]) < 2e-2


def test_cache_lu_tridiagonal_complexe_et_reel():
    from gueye.kernel import TridiagInverse
    rng = np.random.default_rng(1)
    for l, d, u in [(-1.0, 2.3, -1.0), (-1.0 + 0.2j, 2.0 + 0.3j, -1.0)]:
        T = TridiagInverse(30, l, d, u, u1=2 * u, lN=2 * l)
        b = rng.standard_normal(30) + 1j * rng.standard_normal(30)
        for bb in (b, b.real):
            assert np.abs(T.dense() @ T.solve(bb) - bb).max() < 1e-12
def test_pml_proche_transparente():
    pb, g = recipes.open_field((3.0, 3.0), (90, 90), 10.0, (1.5, 1.5), pml=16, R0=1e-8)
    _, ga = recipes.open_field((3.0, 3.0), (90, 90), 10.0, (1.5, 1.5), exact=True)
    assert rel(g, ga) < 1e-2


def test_obstacle_dirichlet():
    pb0 = gy.problem((3.0, 3.0), (60, 60), [("D", "D"), ("P", "P")], k=8.0)       # PML : champ nul sur l'obstacle
    t0, (a, b_) = recipes.dirichlet_obstacle(pb0, (0.6, 1.5), lambda x, y: np.hypot(x - 1.8, y - 1.5) < 0.4)
    assert np.abs(t0[a, b_]).max() < 1e-9 * np.abs(t0).max()
    pb = gy.problem((3.0, 3.0), (60, 60), [("D", "D"), (("R", 8j), ("R", 8j))], k=8.0)   # Robin : comparable au creux
    tot, (Sy, Sx) = recipes.dirichlet_obstacle(pb, (0.6, 1.5), lambda x, y: np.hypot(x - 1.8, y - 1.5) < 0.4)
    assert np.abs(tot[Sy, Sx]).max() < 1e-9 * np.abs(tot).max()
    # equivalent creux : noeuds de l'obstacle elimines
    P = pb.backend
    A = P.dense().tocsc(); n = A.shape[0]
    idx = (Sy - P.Py.i0) * P.Nx + (Sx - P.Px.i0)
    free = np.setdiff1d(np.arange(n), idx)
    b = np.zeros(n, complex); i, j = pb.node((0.6, 1.5)); b[(i - P.Py.i0) * P.Nx + (j - P.Px.i0)] = 1 / (P.hx * P.hy)
    u = np.zeros(n, complex); u[free] = spla.spsolve(A[free][:, free].tocsc(), b[free])
    ref = u.reshape(P.Ny, P.Nx)
    assert rel(tot[P.Py.i0:P.Py.i1 + 1, P.Px.i0:P.Px.i1 + 1], ref) < 1e-10


def test_wave_point_source_conserve_energie_ordre():
    W = gy.wave((2.0, 1.4), (40, 28), [("D", "N"), ("N", "D")], c=1.0)
    src = W.point_source((0.8, 0.6), lambda t: np.exp(-((t - 0.5) / 0.15)**2))
    t, U, R = W.march(0.02, 100, source=src, receivers=[(1.5, 1.0)])
    assert np.isfinite(U).all() and np.abs(U).max() > 0 and R.shape == (101, 1)


def test_grid_dispersion():
    npw, err = recipes.grid_dispersion(5.0, 0.05)
    assert 24 < npw < 26 and 0 < err < 4e-3


def test_obstacle_geometrie_xy():
    pb = gy.problem((4.0, 2.0), (80, 40), [("D", "D"), ("D", "D")], k=3.0)
    _, (Sy, Sx) = recipes.dirichlet_obstacle(pb, (0.8, 1.0), lambda x, y: np.hypot(x - 2.6, y - 0.7) < 0.3)
    assert abs(pb.axes[0][Sx].mean() - 2.6) < 0.05 and abs(pb.axes[1][Sy].mean() - 0.7) < 0.05


def test_block_equals_entry():
    pb = gy.problem((2.0, 1.4), (30, 22), [("N", "D"), (("R", 1.1), "N")], k=1.7)
    P = pb.backend
    Sy = np.array([3, 5, 5, 9]); Sx = np.array([4, 4, 12, 20]); Ty = np.array([2, 8]); Tx = np.array([7, 15])
    B = P.block(Sy, Sx, Ty, Tx)
    E = P.entry(Sy[:, None], Sx[:, None], Ty[None, :], Tx[None, :])
    assert np.abs(B - E).max() < 1e-12 * np.abs(E).max()


@pytest.mark.parametrize("size,cells,bc", [
    ((1.0,), (20,), ("D", "D")),
    ((1.0, 1.0), (8, 7), [("D", "D"), ("D", "D")]),
    ((1.0, 1.0, 1.0), (6, 5, 4), [("D", "D")] * 3),
])
def test_resonances_incluent_mu(size, cells, bc):
    pb = gy.problem(size, cells, bc, mu=2.0, check=False)
    got = pb.resonances(3)
    A = pb.backend.dense()
    A = A.toarray() if hasattr(A, "toarray") else A
    expected = np.sqrt(np.sort(np.linalg.eigvals(A).real)[:3])
    assert np.allclose(got, expected, rtol=1e-9, atol=1e-10)


def test_dirichlet_complexe_1d_solve_et_at():
    pb = gy.problem((1.0,), (12,), (("D", 1 + 2j), ("D", 3 - 1j)), check=False)
    u = pb.solve(0.0)
    assert np.iscomplexobj(u)
    assert np.allclose(u[[0, -1]], [1 + 2j, 3 - 1j])
    assert np.allclose(pb.at([(0.0,), (1.0,)], 0.0), [1 + 2j, 3 - 1j])


def test_dirichlet_complexe_2d_solve_et_at():
    pb = gy.problem((1.0, 1.0), (10, 8), [(("D", 1 + 1j), ("D", 0.0)), ("N", "N")], check=False)
    u = pb.solve(0.0)
    assert np.iscomplexobj(u)
    assert np.allclose(u[:, 0], 1 + 1j)
    assert np.allclose(pb.at([(0.0, 0.5)], 0.0), [1 + 1j])


def test_dirichlet_complexe_3d_solve_et_at():
    bc = [(("D", 1 - 2j), ("D", 0.0)), ("N", "N"), ("N", "N")]
    pb = gy.problem((1.0, 1.0, 1.0), (6, 5, 4), bc, check=False)
    u = pb.solve(0.0)
    assert np.iscomplexobj(u)
    assert np.allclose(u[:, :, 0], 1 - 2j)
    assert np.allclose(pb.at([(0.0, 0.5, 0.5)], 0.0), [1 - 2j])


def test_dirichlet_complexe_fonctionnelle_2d_3d():
    pb2 = gy.problem((1.0, 1.0), (8, 6), [(("D", lambda y: (1 + 1j) * (1 + y)), ("D", 0.0)), ("N", "N")], check=False)
    u2 = pb2.solve(0.0)
    assert np.allclose(u2[:, 0], (1 + 1j) * (1 + pb2.axes[1]))

    bc3 = [(("D", lambda y, z: (1 - 2j) * (1 + y + z)), ("D", 0.0)), ("N", "N"), ("N", "N")]
    pb3 = gy.problem((1.0, 1.0, 1.0), (5, 4, 3), bc3, check=False)
    u3 = pb3.solve(0.0)
    expected = (1 - 2j) * (1 + pb3.axes[1][None, :, None] + pb3.axes[2][:, None, None])
    assert np.allclose(u3[:, :, 0], expected[:, :, 0])


def test_resonance_mu_nonzero_signalee_par_le_controle_singularite():
    args = ((1.0,), (20,), ("D", "D"))
    k1 = gy.problem(*args, mu=2.0, check=False).resonances(1)[0]
    with pytest.raises(gy.SingularOperator):
        gy.problem(*args, k=k1, mu=2.0, check=True)
