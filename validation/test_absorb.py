import sys, warnings, time; pass
import numpy as np
from gueye import Problem1D, Problem2D, Problem3D, Wave1D
warnings.simplefilter("ignore")
rel = lambda a, b: np.abs(a - b).max() / np.abs(b).max()
print("1) Exactitude de l'elimination de la couche PML : probleme 1D avec PML vs matrice complete (noeuds PML explicites)")
L, M, m, R0, p, k = 2.0, 80, 12, 1e-6, 2, 5.0
h = L / M
def full_pml(L, M, m, R0, p, k):
    delta = m * h; smax = (p + 1) * np.log(1 / R0) / (2 * delta)
    n = M + 1 + 2 * (m + 1)
    xs = (np.arange(n) - (m + 1)) * h
    def sig(x):
        t = np.where(x < 0, np.maximum(-x - h / 2, 0), np.maximum(x - L - h / 2, 0)); return smax * (t / delta) ** p
    s = lambda x: 1 - 1j * sig(x) / k
    A = np.zeros((n, n), complex)
    for i in range(n):
        xi = xs[i]; sp_, sm_ = s(xi + h / 2), s(xi - h / 2)
        # half-node : profile defined from the interface (t = 0 at x = 0 / L), nodes at distance t
        A[i, i] = (1 / (s(xi) * h**2)) * (1 / sp_ + 1 / sm_) - k**2
        if i + 1 < n: A[i, i + 1] = -1 / (s(xi) * sp_ * h**2)
        if i - 1 >= 0: A[i, i - 1] = -1 / (s(xi) * sm_ * h**2)
    return A[1:-1, 1:-1], xs[1:-1]
A, xs = full_pml(L, M, m, R0, p, k)
P = Problem1D(L, M, ("P", m, R0, p), ("P", m, R0, p), c0=-k**2)
src = 25; b = np.zeros(A.shape[0], complex); b[m + src] = 1 / h
uref = np.linalg.solve(A, b)[m:m + M + 1]
print(f"   ecart (champ physique) = {rel(P.green(src), uref):.1e}")

print("\n2) Reflexion residuelle 1D vs condition transparente exacte ('A'), k h = 0.2 : erreur relative max sur le champ de Green")
M = 400; L = 4.0; k = 0.2 * M / L; h = L / M
Pa = Problem1D(L, M, ("A",), ("A",), c0=-k**2); ga = Pa.green(150)
print(f"   (periode = {2*np.pi/k:.2f}, longueur d'onde = {2*np.pi/k/h:.0f} mailles)")
print(f"   {'m':>3} | R0=1e-2 | 1e-4 | 1e-6 | 1e-8")
for m in (4, 8, 16, 32):
    row = [rel(Problem1D(L, M, ("P", m, R0), ("P", m, R0), c0=-k**2).green(150), ga) for R0 in (1e-2, 1e-4, 1e-6, 1e-8)]
    print(f"   {m:>3} | " + " | ".join(f"{e:.1e}" for e in row))
Pr = Problem1D(L, M, ("R", 1j * k, 0.0), ("R", 1j * k, 0.0), c0=-k**2)
print(f"   Robin i k (incidence normale) : {rel(Pr.green(150), ga):.1e}")
print(f"   Dirichlet (aucune absorption) : {rel(Problem1D(L, M, ('D', 0.), ('D', 0.), c0=-k**2).green(150), ga):.1e}")

print("\n3) 2D : mode transverse nu de la boite (x) -> onde sortant en y avec angle theta ; erreur vs transparente exacte")
Lx, Ly, Mx, My = 3.0, 3.0, 120, 120
k = 14.0
def run(bc_y, nu):
    P = Problem2D(Lx, Ly, Mx, My, left=("D", 0.), right=("D", 0.), bottom=bc_y, top=bc_y, c0=-k**2)
    xs = P.x; f = np.zeros(P.shape); f[My // 2, :] = np.sin(nu * np.pi * xs / Lx) / P.hy
    return P.solve(f)
print(f"   {'nu':>3} | angle | PML m=8 (R0=1e-6) | PML m=16 | Robin i k | Dirichlet")
rows = []
for nu in (1, 3, 5, 7, 9, 11, 12):
    kx = nu * np.pi / Lx
    if kx >= k: continue
    th = np.degrees(np.arcsin(kx / k))
    ua = run(("A",), nu)
    e = [rel(run(("P", 8, 1e-6), nu), ua), rel(run(("P", 16, 1e-6), nu), ua), rel(run(("R", 1j * np.sqrt(k**2 - kx**2), 0.0), nu), ua), rel(run(("D", 0.), nu), ua)]
    rows.append((th, e)); print(f"   {nu:>3} | {th:4.0f}° | {e[0]:.1e} | {e[1]:.1e} | {e[2]:.1e} | {e[3]:.1e}")
print("   (Robin ici avec alpha = i k_y exact pour chaque mode ; avec alpha = i k fixe l'erreur croit avec l'angle :", end=" ")
print(", ".join(f"{rel(run(('R', 1j * k, 0.0), nu), run(('A',), nu)):.1e}" for nu in (1, 5, 9, 12)), ")")

print("\n4) 2D, source ponctuelle (tous les angles), PML sur les faces haut/bas, parois D a gauche/droite")
P_a = Problem2D(Lx, Ly, Mx, My, left=("D", 0.), right=("D", 0.), bottom=("A",), top=("A",), c0=-k**2)
ga = P_a.green(60, 40)
for m in (8, 16, 32):
    for R0 in (1e-4, 1e-8, 1e-12):
        Pp = Problem2D(Lx, Ly, Mx, My, left=("D", 0.), right=("D", 0.), bottom=("P", m, R0), top=("P", m, R0), c0=-k**2)
        print(f"   m = {m:2d}, R0 = {R0:.0e} : {rel(Pp.green(60, 40), ga):.1e}")
print(f"   verification : transparente exacte (A) = resolution creuse sur grand domaine ? -> ", end="")
# grand domaine : bandes Dirichlet lointaines + amortissement fort (k complexe) pour reference independante
import scipy.sparse.linalg as spla
Pg = Problem2D(Lx, 3 * Ly, Mx, 3 * My, left=("D", 0.), right=("D", 0.), bottom=("D", 0.), top=("D", 0.), c0=-(k + 2.5j)**2)
Pa2 = Problem2D(Lx, Ly, Mx, My, left=("D", 0.), right=("D", 0.), bottom=("A",), top=("A",), c0=-(k + 2.5j)**2)
gg = Pg.green(My + 60, 40)[My:2 * My + 1]; g2 = Pa2.green(60, 40)
print(f"{rel(gg, g2):.1e} (reference avec amortissement k+2.5i, boite 3x plus haute)")

print("\n5) 3D : PML sur z (haut/bas) vs transparente exacte, parois N sur x et y")
P3 = lambda bz: Problem3D(2.0, 2.0, 2.0, 40, 40, 40, left=("N", 0.), right=("N", 0.), front=("D", 0.), back=("D", 0.), bottom=bz, top=bz, c0=-9.0**2)
g3 = P3(("A",)).green(20, 20, 14)
for m in (6, 12):
    print(f"   m = {m:2d}, R0 = 1e-6 : {rel(P3(('P', m, 1e-6)).green(20, 20, 14), g3):.1e}")

print("\n6) Temps 1D ouvert (synthese de frequence avec PML) vs pas de temps sur grand domaine")
Lr, Mr = 3.0, 300; h = Lr / Mr
pass
wp = 8.0; ricker = lambda t: (1 - 2 * (0.5 * wp * (t - 1.2))**2) * np.exp(-(0.5 * wp * (t - 1.2))**2)
tt = np.linspace(0, 8, 321)
src, rec = 100, np.array([150, 240])
Us = {}
for name, bc in (("A exacte", ("A",)), ("PML m=12", ("P", 12, 1e-6)), ("PML m=24", ("P", 24, 1e-6)), ("PML m=24, R0=1e-10", ("P", 24, 1e-10))):
    Us[name] = Wave1D(Lr, Mr, bc, bc, c=1.0).impulse_response(src, rec, tt, wmax=24.0, nw=600, source=ricker, shift=0.05)
Lb = 3.0 + 2 * 6.0; Mb = int(Lb / h); off = int(6.0 / h)
Wb = Wave1D(Lb, Mb, ("D", 0.), ("D", 0.), c=1.0)
dt = 0.004
fs = lambda t: (lambda x: np.where(np.arange(Mb + 1) == off + src, ricker(t) / h, 0.0))
t, _, R = Wb.march(dt, int(8 / dt), f=fs, receivers=off + rec)
Rb = np.array([np.interp(tt, t, R[:, i]) for i in range(2)])
for name, U in Us.items():
    print(f"   {name:>20} : ecart max / amplitude = {np.abs(U - Rb).max() / np.abs(Rb).max():.1e}")
print(f"   (reference : domaine 5x plus grand, Dirichlet lointain, Newmark dt = {dt} ; sans absorption l'ecart serait {np.abs(Wave1D(Lr, Mr, ('D', 0.), ('D', 0.)).impulse_response(src, rec, tt, wmax=24.0, nw=600, source=ricker, shift=0.05) - Rb).max() / np.abs(Rb).max():.1e})")
