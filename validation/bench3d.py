import sys, warnings, time; pass
import numpy as np, scipy.sparse.linalg as spla, scipy.sparse as sp
from scipy.fft import dstn
from gueye import Problem3D, Wave3D, SingularOperator
warnings.simplefilter("ignore")
rel = lambda a, b: np.abs(a - b).max() / np.abs(b).max()
Lx=Ly=Lz=1.0
print("\n5) Echelle : Poisson DDDDDD (cube) -- formule modale vs DST (FFT) vs LU creux")
print(f"{'N^3':>6} | formule init+solve | DST | LU creux | acces 1 point | ecart formule/DST")
for N in (24, 32, 48, 64, 96):
    t0 = time.perf_counter(); P = Problem3D(1, 1, 1, N + 1, N + 1, N + 1, *[("D", 0.)] * 6); ti = time.perf_counter() - t0
    rng = np.random.default_rng(0); F = rng.standard_normal((N, N, N))
    t0 = time.perf_counter(); u = P.solve_reduced(F); ts = time.perf_counter() - t0
    h = 1.0 / (N + 1); j = np.arange(1, N + 1); lam1 = (2 - 2 * np.cos(j * np.pi / (N + 1))) / h**2
    t0 = time.perf_counter(); Fh = dstn(F, type=1); Uh = Fh / (lam1[:, None, None] + lam1[None, :, None] + lam1[None, None, :]); ud = dstn(Uh, type=1) / (2 * (N + 1))**3; td = time.perf_counter() - t0
    tl = np.nan
    if N <= 32:
        A = P.dense().tocsc(); t0 = time.perf_counter(); ul = spla.spsolve(A, F.ravel()); tl = time.perf_counter() - t0
    t0 = time.perf_counter(); P.at([(N // 2, N // 3, N // 4)], lambda x, y, z: 1 + 0 * x); ta = time.perf_counter() - t0
    print(f"{N:>6} | {ti:6.2f} + {ts:6.3f} s | {td:.3f} s | {tl:6.2f} s | {ta*1e3:7.0f} ms | {rel(u, ud):.1e}")
