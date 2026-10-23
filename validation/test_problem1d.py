import sys, warnings; pass
import numpy as np
from gueye.problem1d import Problem1D
from gueye.kernel import SingularOperator
warnings.simplefilter("ignore")
rel = lambda a, b: np.abs(a - b).max() / np.abs(b).max()
L = 2.0
alpha = 1.7
BCs = {"DD": (("D",), ("D",)), "DN": (("D",), ("N",)), "ND": (("N",), ("D",)), "NN": (("N",), ("N",)),
       "RR": (("R", alpha), ("R", alpha)), "RN": (("R", alpha), ("N",)), "NR": (("N",), ("R", alpha)),
       "DR": (("D",), ("R", alpha)), "RD": (("R", alpha), ("D",))}
# solution manufacturee (pas de symetrie particuliere) ; CL inhomogenes lues sur la solution exacte
uex  = lambda x: np.cos(2.3 * x) + 0.4 * x**2 - 0.3 * x
du   = lambda x: -2.3 * np.sin(2.3 * x) + 0.8 * x - 0.3
d2u  = lambda x: -2.3**2 * np.cos(2.3 * x) + 0.8
def make(name, M, k=None, c1=0.0, c0=None):
    lt, rt = BCs[name]
    gl = lambda b: {"D": uex(0), "N": -du(0), "R": -du(0) + alpha * uex(0)}[b[0]]     # normale sortante a gauche = -d/dx
    gr = lambda b: {"D": uex(L), "N": du(L), "R": du(L) + alpha * uex(L)}[b[0]]
    left = lt + (gl(lt),); right = rt + (gr(rt),)
    c0 = (-k**2 if k is not None else 0.0) if c0 is None else c0
    P = Problem1D(L, M, left, right, c2=1.0, c1=c1, c0=c0)
    f = lambda x: -d2u(x) + c1 * du(x) + c0 * uex(x)
    return P, f

print("1) Formule (green) = inversion dense = bande, 9 CL x k, M = 60 (erreur relative)")
print(f"{'CL':>3} | " + " | ".join(f"k={k}" .rjust(9) for k in (0.0, 0.7, 3.3, "0.9+0.3i")))
for name in BCs:
    row = []
    for k in (0.0, 0.7, 3.3, 0.9 + 0.3j):
        try:
            P, f = make(name, 60, k)
            ug = P.solve(f, "green"); ub = P.solve(f, "banded")
            fr = f(P.x); rhs = P.rhs(f)
            ud = np.linalg.solve(P.dense(), rhs)
            row.append(f"{max(rel(ug[P.i0:P.i1+1], ud), rel(ub[P.i0:P.i1+1], ud)):.1e}")
        except SingularOperator:
            row.append("singulier")
    print(f"{name:>3} | " + " | ".join(r.rjust(9) for r in row))

print("\n2) Ordre de convergence vers la solution exacte du probleme continu (h -> h/2), k = 1.3, c1 = 0.4")
print(f"{'CL':>3} | err(M=40) | err(80) | err(160) | ordre")
for name in BCs:
    es = []
    for M in (40, 80, 160):
        P, f = make(name, M, 1.3, c1=0.4)
        es.append(np.abs(P.solve(f) - uex(P.x)).max())
    print(f"{name:>3} | {es[0]:9.2e} | {es[1]:7.2e} | {es[2]:8.2e} | {np.log2(es[1]/es[2]):.2f}")

print("\n3) Resonances exactes du schema vs continu, M = 200 (k_n h petit : ecart relatif <= (kh)^2/24)")
for name in ("DD", "DN", "ND", "NN", "RR", "DR"):
    lt, rt = BCs[name]
    P = Problem1D(L, 200, lt + (0.0,), rt + (0.0,), c2=1.0, check=False)
    kr = P.resonances(4).real
    # continu
    from scipy.optimize import brentq
    def cont(n):
        if name in ("DD",): return n * np.pi / L
        if name == "NN": return (n - 1) * np.pi / L
        if name in ("DN", "ND"): return (n - 0.5) * np.pi / L
        if name == "RR":     # ((k^2 - a^2) sin kL = 2 a k cos kL ... (R: u_n + a u = 0 des deux cotes)
            F = lambda k: (k * k - alpha**2) * np.sin(k * L) - 2 * alpha * k * np.cos(k * L)
        else:                # DR : u(0)=0, u'(L)+a u(L)=0 -> k cos kL + a sin kL = 0
            F = lambda k: k * np.cos(k * L) + alpha * np.sin(k * L)
        ks = np.linspace(1e-6, 12, 20000); v = F(ks); roots = [brentq(F, ks[i], ks[i+1]) for i in range(len(ks)-1) if v[i] * v[i+1] < 0]
        return roots[n - 1]
    ref = np.array([cont(n) for n in range(1, 5)]) if name not in ("NN",) else np.array([cont(n) for n in range(1, 5)])
    kr = np.sort(kr)
    kk = kr[: 4]
    if name == "RR" and kk[0] < 1e-9: pass
    singular = []
    for kn in kk[:3]:
        Pk = Problem1D(L, 200, lt + (0.0,), rt + (0.0,), c2=1.0, c0=-kn**2, check=False); singular.append(f"{Pk.K.rcond:.0e}")
    print(f"  {name}: k_n (schema) = {np.round(kk, 5)} ; continu = {np.round(ref[:len(kk)], 5)} ; ecarts relatifs = "
          f"{np.round(np.abs(kk[1:]-ref[1:len(kk)]) / ref[1:len(kk)], 7) if name=='NN' else np.round(np.abs(kk-ref[:len(kk)])/ref[:len(kk)],7)} ; det. relatif a k_n : {singular}")

print("\n4) Acces selectif : valeur en 3 points d'un probleme M = 10^6 (DN, k = 5) sans resoudre le champ")
import time
M = 1_000_000
P, f = make("DN", M, 5.0)
t0 = time.perf_counter(); v = P.at([M // 3, M // 2, M - 1], f); t1 = time.perf_counter() - t0
t0 = time.perf_counter(); ub = P.solve(f); t2 = time.perf_counter() - t0
print(f"   formule : {t1*1e3:.0f} ms ; resolution bande (champ complet) : {t2*1e3:.0f} ms ; ecart {np.abs(v - ub[[M//3, M//2, M-1]]).max():.1e}")
