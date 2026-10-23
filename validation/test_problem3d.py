import sys, warnings, time; pass
import numpy as np, scipy.sparse.linalg as spla, scipy.sparse as sp
from scipy.fft import dstn
from gueye import Problem3D, Wave3D, SingularOperator
warnings.simplefilter("ignore")
rel = lambda a, b: np.abs(a - b).max() / np.abs(b).max()
Lx, Ly, Lz, al, k = 2.0, 1.5, 1.2, 1.6, 1.1
uex = lambda x, y, z: np.cos(1.3 * x + .2) * np.sin(1.7 * y + .1) * np.cos(1.1 * z) + 0.2 * x * y * z + 0.1 * z**2
gx = lambda x, y, z: -1.3 * np.sin(1.3 * x + .2) * np.sin(1.7 * y + .1) * np.cos(1.1 * z) + 0.2 * y * z
gy = lambda x, y, z: 1.7 * np.cos(1.3 * x + .2) * np.cos(1.7 * y + .1) * np.cos(1.1 * z) + 0.2 * x * z
gz = lambda x, y, z: -1.1 * np.cos(1.3 * x + .2) * np.sin(1.7 * y + .1) * np.sin(1.1 * z) + 0.2 * x * y + 0.2 * z
lap = lambda x, y, z: -(1.3**2 + 1.7**2 + 1.1**2) * np.cos(1.3 * x + .2) * np.sin(1.7 * y + .1) * np.cos(1.1 * z) + 0.2
def face(kd, name):
    f = {"left": (lambda y, z: uex(0 * y, y, z), lambda y, z: -gx(0 * y, y, z)),
         "right": (lambda y, z: uex(Lx + 0 * y, y, z), lambda y, z: gx(Lx + 0 * y, y, z)),
         "front": (lambda x, z: uex(x, 0 * x, z), lambda x, z: -gy(x, 0 * x, z)),
         "back": (lambda x, z: uex(x, Ly + 0 * x, z), lambda x, z: gy(x, Ly + 0 * x, z)),
         "bottom": (lambda x, y: uex(x, y, 0 * x), lambda x, y: -gz(x, y, 0 * x)),
         "top": (lambda x, y: uex(x, y, Lz + 0 * x), lambda x, y: gz(x, y, Lz + 0 * x))}[name]
    u, dn = f
    return {"D": ("D", u), "N": ("N", dn), "R": ("R", al, lambda a, b: dn(a, b) + al * u(a, b))}[kd]
names = ("left", "right", "front", "back", "bottom", "top")
def make(combo, Mx, My, Mz, kk=k):
    P = Problem3D(Lx, Ly, Lz, Mx, My, Mz, c0=-kk**2, **{n: face(c, n) for n, c in zip(names, combo)})
    return P, (lambda x, y, z: -lap(x, y, z) - kk**2 * uex(x, y, z))
combos = ["DDDDDD", "NNNNNN", "RRRRRR", "DNDNDN", "NDRDND", "RDNNDR", "DRRNND", "NRDDRN"]
print("1) Formule modale = assemblage creux, 9 x 8 x 7 ; ordre des faces : gauche droite avant arriere bas haut")
print(f"{'CL':>7} |  k reel | k complexe | acces selectif | entree = colonne")
for c in combos:
    row = []
    for kk in (k, k + 0.4j):
        P, f = make(c, 9, 8, 7, kk)
        F = P.rhs(f); u = P.solve_reduced(F)
        ref = spla.spsolve(P.dense().tocsc(), F.ravel().astype(complex)).reshape(F.shape)
        row.append(rel(u, ref))
    pts = [(P.Pz.i0 + 2, P.Py.i0 + 3, P.Px.i0 + 1), (P.Pz.i1 - 1, P.Py.i1 - 2, P.Px.i1)]
    sel = P.at(pts, f); selr = np.array([ref[a - P.Pz.i0, b - P.Py.i0, cc - P.Px.i0] for a, b, cc in pts])
    j = (P.Pz.i0 + 3, P.Py.i0 + 2, P.Px.i0 + 4); col = P.green(*j)
    e = P.entry(P.Pz.i0 + 1, P.Py.i0 + 4, P.Px.i0 + 2, *j)
    print(f"{c:>7} | {row[0]:7.1e} | {row[1]:10.1e} | {rel(sel, selr):14.1e} | {abs(e - col[P.Pz.i0 + 1, P.Py.i0 + 4, P.Px.i0 + 2]) / abs(e):.1e}")
print("\n2) Convergence vers la solution exacte (CL inhomogenes) M = 12, 24, 48 (x), proportionnel ensuite")
for c in ("DDDDDD", "NNNNNN", "RRRRRR", "NDRDND", "RDNNDR"):
    es = []
    for M in (12, 24, 48):
        P, f = make(c, M, int(M * Ly / Lx), int(M * Lz / Lx))
        U = P.solve(f); Z, Y, X = np.meshgrid(P.z, P.y, P.x, indexing="ij")
        es.append(np.abs(U - uex(X, Y, Z)).max())
    print(f"  {c} : {es[0]:.2e} {es[1]:.2e} {es[2]:.2e} ordre {np.log2(es[1]/es[2]):.2f}")
print("\n3) Singularite NNNNNN a k=0 ->", end=" ")
try: Problem3D(1, 1, 1, 10, 10, 10, *[("N", 0.)] * 6, c0=0.0); print("PAS DETECTE")
except SingularOperator as e: print("detecte :", str(e)[:50])
print("\n4) Temps : Wave3D mixte vs LU creux pas a pas (meme schema)")
bc = dict(left=("D", 0.), right=("R", 1.5, 0.), front=("N", 0.), back=("D", 0.), bottom=("N", 0.), top=("D", 0.))
W = Wave3D(2.0, 1.5, 1.2, 16, 12, 10, c=1.0, eta=0.05, **bc); P = W.P
src = (5, 6, 7); h3 = P.hx * P.hy * P.hz
rick = lambda t: (1 - 2 * (3 * (t - 0.8))**2) * np.exp(-(3 * (t - 0.8))**2)
def fs(t):
    F = np.zeros(P.shape); F[src] = rick(t) / h3; return F
dt, n = 0.03, 60
t, U = W.march(dt, n, f=fs)
A = P.dense().tocsc(); a0 = 1 + 0.05 * 0.5 * dt; s = dt**2 / 4
lu = spla.splu((a0 * sp.identity(A.shape[0]) + s * A).tocsc())
sl = (slice(P.Pz.i0, P.Pz.i1 + 1), slice(P.Py.i0, P.Py.i1 + 1), slice(P.Px.i0, P.Px.i1 + 1))
fv = lambda t: fs(t)[sl].ravel()
u = np.zeros(A.shape[0]); v = u.copy(); a = fv(0) - A @ u
for i in range(n):
    ut = u + dt * v + dt**2 * 0.25 * a; vt = v + dt * 0.5 * a
    a = lu.solve(fv((i + 1) * dt) - A @ ut - 0.05 * vt); u = ut + s * a; v = vt + 0.5 * dt * a
print(f"   ecart champ final : {rel(U[-1][sl].ravel(), u):.1e} (|u|max = {np.abs(U[-1]).max():.2e})")
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
