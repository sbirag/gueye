import sys; pass
import numpy as np, warnings
from gueye.kernel import TridiagInverse, SingularOperator
rng = np.random.default_rng(1)
rel = lambda A, B: np.abs(A - B).max() / np.abs(B).max()
print("Noyau tridiagonal : erreur relative de T^-1 (formule vs inversion dense)")
cases = []
for N in (4, 5, 12, 40):
    for (l, d, u) in [(-1., 2., -1.), (-1., 1.3, -1.), (-1., 2.3, -1.), (-1., 4., -1.), (-2., 5., 0.5), (-1.+.5j, 2.+.3j, -1.-.2j), (-1., 2.+1e-9, -1.), (-1., 1e-9, -1.)]:
        for ends in (None, "N", "R"):
            kw = {}
            if ends == "N": kw = dict(u1=2*u, lN=2*l)
            if ends == "R": kw = dict(d1=d+0.7*l, u1=u+l, dN=d-0.3*u, lN=l+u)
            if d == 1e-9 and N % 2: continue
            try:
                T = TridiagInverse(N, l, d, u, **kw)
            except SingularOperator:
                continue
            Bd = np.linalg.inv(T.dense())
            cases.append((N, l, d, u, ends, rel(T.inverse(), Bd)))
errs = np.array([c[-1] for c in cases])
print(f"  {len(cases)} cas (N, regimes sinh/sin/racine double/complexe/non symetrique, extremites libres/N/R) : erreur max {errs.max():.1e}")
w = max(cases, key=lambda c: c[-1]); print("  pire cas :", w)
# singularites
for name, kw, N in [("Neumann pur d=2", dict(u1=-2., lN=-2.), 20), ("d=0, N impair", {}, 21)]:
    try:
        TridiagInverse(N, -1., 2. if "Neumann" in name else 0., -1., **kw); print("  PAS DETECTE", name)
    except SingularOperator as e: print(f"  singularite detectee : {name}")
# grand N stable
N = 400000
T = TridiagInverse(N, -1., 2.-0.25**2, -1.)
b = np.zeros(N); b[N//2] = 1
from scipy.linalg import solve_banded
x = T.solve(b)
print(f"  N = {N}, k=0.25 : |G(N/2,N/2)|={T.entries(N//2,N//2):.6f}, ecart a la resolution bande {abs(T.entries(N//2+7, N//2)-x[N//2+7]):.1e}, fini={np.isfinite(T.column(N//2)).all()}")
