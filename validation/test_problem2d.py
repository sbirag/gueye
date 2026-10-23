import sys, warnings, time, itertools; pass
import numpy as np, scipy.sparse.linalg as spla
from gueye import Problem2D, Wave2D, SingularOperator
warnings.simplefilter("ignore")
rel = lambda a, b: np.abs(a - b).max() / np.abs(b).max()
Lx, Ly, al = 2.0, 1.4, 1.6
kind = ["D", "N", "R"]
# solution manufacturee 2D
uex = lambda x, y: np.cos(1.7 * x + 0.3) * np.sin(2.1 * y + 0.2) + 0.3 * x * y + 0.2 * y**2
ux = lambda x, y: -1.7 * np.sin(1.7 * x + 0.3) * np.sin(2.1 * y + 0.2) + 0.3 * y
uy = lambda x, y: 2.1 * np.cos(1.7 * x + 0.3) * np.cos(2.1 * y + 0.2) + 0.3 * x + 0.4 * y
lap = lambda x, y: -(1.7**2 + 2.1**2) * np.cos(1.7 * x + 0.3) * np.sin(2.1 * y + 0.2) + 0.4
def side(kd, name):
    if name == "left":   dn, xb = (lambda y: -ux(0, y)), (lambda y: uex(0 * y, y)); rob = lambda y: -ux(0, y) + al * uex(0 * y, y)
    if name == "right":  dn = lambda y: ux(Lx, y); xb = lambda y: uex(Lx + 0 * y, y); rob = lambda y: ux(Lx, y) + al * uex(Lx + 0 * y, y)
    if name == "bottom": dn = lambda x: -uy(x, 0); xb = lambda x: uex(x, 0 * x); rob = lambda x: -uy(x, 0) + al * uex(x, 0 * x)
    if name == "top":    dn = lambda x: uy(x, Ly); xb = lambda x: uex(x, Ly + 0 * x); rob = lambda x: uy(x, Ly) + al * uex(x, Ly + 0 * x)
    return {"D": ("D", xb), "N": ("N", dn), "R": ("R", al, rob)}[kd]
def make(combo, Mx, My, k=1.1):
    s = {n: side(kd, n) for n, kd in zip(("left", "right", "bottom", "top"), combo)}
    P = Problem2D(Lx, Ly, Mx, My, c0=-k**2, **s)
    f = lambda x, y: -lap(x, y) - k**2 * uex(x, y)
    return P, f

combos = [("D",)*4, ("N",)*4, ("R",)*4, ("D", "N", "D", "N"), ("N", "D", "R", "D"), ("R", "D", "N", "N"), ("D", "R", "R", "N"),
          ("N", "R", "D", "D"), ("R", "N", "R", "D"), ("D", "D", "N", "R"), ("N", "N", "D", "R"), ("R", "R", "N", "N")]
print("1) Formule modale = assemblage creux (resolution directe), Nx = 30, Ny = 24 ; k = 1.1 puis 1.1+0.4i ; ordre = 'gauche droite bas haut'")
print(f"{'CL':>9} |  k reel | k complexe | acces selectif | G(entry)=colonne")
for combo in combos:
    row = []
    for k in (1.1, 1.1 + 0.4j):
        P, f = make(combo, 30, 24, k)
        F = P.rhs(f); u = P.solve_reduced(F)
        ref = spla.spsolve(P.dense().tocsc(), F.ravel().astype(complex)).reshape(F.shape)
        row.append(rel(u, ref))
    pts = [(P.Py.i0 + 3, P.Px.i0 + 5), (P.Py.i1 - 2, P.Px.i1 - 4)]
    sel = P.at(pts, f); selr = np.array([ref[a - P.Py.i0, b - P.Px.i0] for a, b in pts])
    jy, jx = P.Py.i0 + 7, P.Px.i0 + 9
    col = P.green(jy, jx); e = P.entry(np.arange(P.Py.i0, P.Py.i1 + 1)[:, None], np.arange(P.Px.i0, P.Px.i1 + 1)[None, :], jy, jx)
    print(f"{''.join(combo):>9} | {row[0]:7.1e} | {row[1]:10.1e} | {rel(sel, selr):14.1e} | {rel(e, col[P.Py.i0:P.Py.i1+1, P.Px.i0:P.Px.i1+1]):.1e}")

print("\n2) Convergence vers la solution exacte (CL inhomogenes), Mx = 2 My = 20, 40, 80")
for combo in combos[:9]:
    es = []
    for M in (20, 40, 80):
        P, f = make(combo, M, int(M * Ly / Lx), 1.1)
        U = P.solve(f); Y, X = np.meshgrid(P.y, P.x, indexing="ij")
        es.append(np.abs(U - uex(X, Y)).max())
    print(f"  {''.join(combo)} : {es[0]:.2e} {es[1]:.2e} {es[2]:.2e} ordre {np.log2(es[1]/es[2]):.2f}")

print("\n3) Singularites : NNNN a k=0 ->", end=" ")
try: Problem2D(1, 1, 20, 20, ("N", 0.), ("N", 0.), ("N", 0.), ("N", 0.), c0=0.0); print("PAS DETECTE")
except SingularOperator as e: print("detecte (", str(e)[:60], "...)")
P = Problem2D(Lx, Ly, 40, 28, ("D", 0.), ("R", al, 0.), ("N", 0.), ("D", 0.), check=False)
kr = P.resonances(4).real; print("   resonances du schema DR-ND :", np.round(kr, 4))
for kk in kr[:2]:
    Pk = Problem2D(Lx, Ly, 40, 28, ("D", 0.), ("R", al, 0.), ("N", 0.), ("D", 0.), c0=-kk**2, check=False)
    print(f"   k = {kk:.4f} : plus petit determinant relatif des modes = {min(K.rcond for K in Pk.modes()):.1e}")

print("\n4) Temps : Wave2D mixte vs resolution creuse pas a pas (Newmark identique), 3 CL differentes")
W = Wave2D(2.0, 1.4, 40, 28, ("D", 0.), ("R", 1.5, 0.), ("N", 0.), ("D", 0.), c=1.0, eta=0.05)
src = (14, 12)
rick = lambda t: (1 - 2 * (3 * (t - 0.8))**2) * np.exp(-(3 * (t - 0.8))**2)
def fsrc(t):
    F = np.zeros(W.P.shape); F[src] = rick(t) / (W.P.hx * W.P.hy); return F
dt = 0.02
t, U, R = W.march(dt, 300, f=fsrc, receivers=[(10, 30), (20, 8)])
import scipy.sparse as sp
P = W.P; A = P.dense().tocsc(); a0 = 1 + 0.05 * 0.5 * dt; s = dt**2 / 4
lu = spla.splu((a0 * sp.identity(A.shape[0]) + s * A).tocsc())
iy0, ix0, Ny, Nx = P.Py.i0, P.Px.i0, P.Ny, P.Nx
u = np.zeros(A.shape[0]); v = u.copy()
def fvec(t): return fsrc(t)[iy0:iy0 + Ny, ix0:ix0 + Nx].ravel()
a = fvec(0) - A @ u
recs = []
for n in range(300):
    ut = u + dt * v + dt**2 * 0.25 * a; vt = v + dt * 0.5 * a
    a = lu.solve(fvec((n + 1) * dt) - A @ ut - 0.05 * vt)
    u = ut + s * a; v = vt + 0.5 * dt * a
recs = u.reshape(Ny, Nx)
print(f"   ecart champ final (formule modale / LU creux) : {rel(U[-1][iy0:iy0+Ny, ix0:ix0+Nx], recs):.1e} ; |u|max = {np.abs(U[-1]).max():.2e}")

print("\n5) Echelle : Poisson 2D (DN/ND), champ complet et acces selectif")
for M in (128, 256, 512):
    P = Problem2D(1.0, 1.0, M, M, ("D", 0.), ("N", 0.), ("N", 0.), ("D", 0.), check=False)
    f = lambda x, y: np.sin(3 * x) * np.cos(2 * y)
    t0 = time.perf_counter(); P2 = Problem2D(1.0, 1.0, M, M, ("D", 0.), ("N", 0.), ("N", 0.), ("D", 0.)); tc = time.perf_counter() - t0
    t0 = time.perf_counter(); u = P2.solve(f); ts = time.perf_counter() - t0
    t0 = time.perf_counter(); v = P2.at([(M // 3, M // 2)], f); ta = time.perf_counter() - t0
    Fr = P2.rhs(f); lu_t = time.perf_counter(); ref = spla.spsolve(P2.dense().tocsc(), Fr.ravel()); tl = time.perf_counter() - lu_t
    print(f"   {M}x{M} : init {tc:.2f} s, solve {ts:.2f} s, 1 point {ta*1e3:.0f} ms ; sparse LU {tl:.2f} s ; ecart {rel(u[P2.Py.i0:P2.Py.i1+1, P2.Px.i0:P2.Px.i1+1].ravel(), ref):.1e}")
