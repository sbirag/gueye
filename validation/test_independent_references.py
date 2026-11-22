"""Independent sparse-reference checks for Green columns, frequency sweeps and obstacles.

Run with: python validation/test_independent_references.py
The sparse reference is assembled independently from the public solver result.
"""
import numpy as np
import scipy.sparse.linalg as spla

import gueye as gy
from gueye import recipes


def relative_error(a, b):
    den = max(float(np.max(np.abs(b))), 1e-300)
    return float(np.max(np.abs(a - b)) / den)


def interior_data(pb):
    """Return sparse matrix, interior slices, interior shape, and cell volume."""
    P = pb.backend
    if pb.dim == 1:
        return P.dense(), (slice(P.i0, P.i1 + 1),), (P.N,), pb.h[0]
    if pb.dim == 2:
        return P.dense().tocsc(), (
            slice(P.Py.i0, P.Py.i1 + 1),
            slice(P.Px.i0, P.Px.i1 + 1),
        ), (P.Ny, P.Nx), pb.h[0] * pb.h[1]
    return P.dense().tocsc(), (
        slice(P.Pz.i0, P.Pz.i1 + 1),
        slice(P.Py.i0, P.Py.i1 + 1),
        slice(P.Px.i0, P.Px.i1 + 1),
    ), (P.Nz, P.Ny, P.Nx), pb.h[0] * pb.h[1] * pb.h[2]


def interior_index(pb, node):
    P = pb.backend
    if pb.dim == 1:
        return node - P.i0
    if pb.dim == 2:
        iy, ix = node
        return (iy - P.Py.i0) * P.Nx + (ix - P.Px.i0)
    iz, iy, ix = node
    return ((iz - P.Pz.i0) * P.Ny + (iy - P.Py.i0)) * P.Nx + (ix - P.Px.i0)


def green_sparse_check(dim, size, cells, bc, source):
    pb = gy.problem(size, cells, bc, k=1.17 + 0.19j)
    A, slices, shape, volume = interior_data(pb)
    node = pb.node(source)
    rhs = np.zeros(A.shape[0], dtype=complex)
    rhs[interior_index(pb, node)] = 1.0 / volume
    reference = spla.spsolve(A, rhs).reshape(shape)
    computed = pb.green(source)[slices]
    return relative_error(computed, reference)


def sweep_sparse_check(dim, size, cells, bc, source, receiver):
    ks = np.array([0.73 + 0.11j, 1.19 + 0.17j, 1.63 + 0.23j])
    base = gy.problem(size, cells, bc, k=ks[0])
    observed = base.sweep(ks, source, [receiver])[:, 0]
    errors = []
    for k, value in zip(ks, observed):
        pb = gy.problem(size, cells, bc, k=k)
        A, _, shape, volume = interior_data(pb)
        rhs = np.zeros(A.shape[0], dtype=complex)
        rhs[interior_index(pb, pb.node(source))] = 1.0 / volume
        ref = spla.spsolve(A, rhs).reshape(shape)
        rec = pb.node(receiver)
        P = pb.backend
        if dim == 1:
            expected = ref[rec - P.i0]
        elif dim == 2:
            expected = ref[rec[0] - P.Py.i0, rec[1] - P.Px.i0]
        else:
            expected = ref[rec[0] - P.Pz.i0, rec[1] - P.Py.i0, rec[2] - P.Px.i0]
        errors.append(abs(value - expected) / max(abs(expected), 1e-300))
    return float(max(errors))


def resonance_sparse_check(dim, size, cells):
    """Compare API wavenumbers with eigenvalues of an independently assembled matrix."""
    mu = 0.37
    pb = gy.problem(size, cells, [("D", "D")] * dim, mu=mu, k=0.0)
    A, _, _, _ = interior_data(pb)
    if dim == 1:
        eigenvalues = np.linalg.eigvalsh(np.asarray(A))
        reference = np.sqrt(np.sort(eigenvalues)[:6])
    else:
        eigenvalues = spla.eigsh(A, k=6, which="SM", return_eigenvectors=False)
        reference = np.sqrt(np.sort(eigenvalues))
    return relative_error(pb.resonances(6), reference)


def obstacle_sparse_check(M, geometry):
    Lx, Ly = 2.0, 1.0
    pb = gy.problem(
        (Lx, Ly), (2 * M, M),
        [("D", "D"), (("R", 8j), ("R", 8j))], k=5.4 + 0.16j
    )
    X, Y = pb.mesh()
    if geometry == "disk":
        inside = lambda x, y: (x - 1.35) ** 2 + (y - 0.5) ** 2 < 0.19 ** 2
    else:
        inside = lambda x, y: (1.15 < x) & (x < 1.55) & (0.34 < y) & (y < 0.66)
    source = (0.35, 0.5)
    total, (Sy, Sx) = recipes.dirichlet_obstacle(pb, source, inside)
    P = pb.backend
    A = P.dense().tocsc()
    n = A.shape[0]
    obstacle = (Sy - P.Py.i0) * P.Nx + (Sx - P.Px.i0)
    free = np.setdiff1d(np.arange(n), obstacle)
    rhs = np.zeros(n, dtype=complex)
    iy, ix = pb.node(source)
    rhs[(iy - P.Py.i0) * P.Nx + (ix - P.Px.i0)] = 1.0 / (P.hx * P.hy)
    reference = np.zeros(n, dtype=complex)
    reference[free] = spla.spsolve(A[free][:, free].tocsc(), rhs[free])
    reference = reference.reshape(P.Ny, P.Nx)
    interior = total[P.Py.i0:P.Py.i1 + 1, P.Px.i0:P.Px.i1 + 1]
    err = relative_error(interior, reference)
    obstacle_max = float(np.max(np.abs(total[Sy, Sx])))
    return err, obstacle_max, int(len(Sy))


def main():
    print("Independent Green columns: GUEYE versus independently assembled sparse systems")
    cases = [
        (1,  (1.0,),       (40,),       [("D", ("R", 1.3))],             (0.4,)),
        (2,  (1.0, 0.8),   (24, 20),    [("D", ("R", 1.3)), ("N", "D")], (0.4, 0.4)),
        (3,  (1.0, 0.8, 0.6), (12, 10, 8), [("D", ("R", 1.3)), ("N", "D"), ("D", "N")], (0.4, 0.4, 0.3)),
    ]
    for dim, size, cells, bc, src in cases:
        err = green_sparse_check(dim, size, cells, bc, src)
        print(f"  {dim}D: relative max error = {err:.3e}")
        assert err < 2e-10

    print("\nFrequency sweeps: sweep() versus independent sparse solve at each complex frequency")
    sweep_cases = [
        (1, (1.0,), (40,), [("D", ("R", 1.3))], (0.4,), (0.65,)),
        (2, (1.0, 0.8), (24, 20), [("D", ("R", 1.3)), ("N", "D")], (0.4, 0.4), (0.7, 0.4)),
        (3, (1.0, 0.8, 0.6), (12, 10, 8), [("D", ("R", 1.3)), ("N", "D"), ("D", "N")], (0.4, 0.4, 0.3), (0.7, 0.4, 0.3)),
    ]
    for dim, size, cells, bc, src, rec in sweep_cases:
        err = sweep_sparse_check(dim, size, cells, bc, src, rec)
        print(f"  {dim}D: max relative error over 3 frequencies = {err:.3e}")
        assert err < 2e-9

    print("\nDiscrete resonances: API versus eigenvalues of independently assembled sparse matrices")
    for dim, size, cells in [(1, (1.0,), (40,)), (2, (1.0, 0.8), (16, 12)), (3, (1.0, 0.8, 0.6), (10, 8, 6))]:
        err = resonance_sparse_check(dim, size, cells)
        print(f"  {dim}D: max relative eigen-wavenumber error = {err:.3e}")
        assert err < 2e-9

    print("\nDirichlet obstacles: capacitance formulation versus independent sparse elimination")
    for M in (32, 64):
        for geometry in ("disk", "rectangle"):
            err, obstacle_max, count = obstacle_sparse_check(M, geometry)
            print(f"  Mx={2*M}, My={M}, {geometry}: error={err:.3e}, max|u_obstacle|={obstacle_max:.3e}, nodes={count}")
            assert err < 2e-9
            assert obstacle_max < 2e-9


if __name__ == "__main__":
    main()
