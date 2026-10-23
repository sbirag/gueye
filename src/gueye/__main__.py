"""python -m gueye : auto-test rapide (quelques secondes) de l'installation."""
import time, warnings
import numpy as np
import scipy.sparse.linalg as spla
import gueye as gy

warnings.simplefilter("ignore")
t0 = time.perf_counter()
rel = lambda a, b: np.abs(a - b).max() / np.abs(b).max()
print(f"gueye {gy.__version__} : auto-test")
ok = True
cases = {1: ((2.0,), (60,), [("D", ("R", 1.5))]),
         2: ((2.0, 1.4), (28, 20), [("N", "D"), (("R", 1.2), "N")]),
         3: ((2.0, 1.5, 1.2), (12, 10, 9), [("D", "N"), ("N", ("R", 1.2)), ("D", "D")])}
for d, (size, cells, bc) in cases.items():
    for k in (1.3, 1.3 + 0.4j):
        pb = gy.problem(size, cells, bc, k=k)
        f = lambda *x: np.cos(2 * x[0]) + 0.3
        u = pb.solve(f)
        P = pb.backend
        if d == 1:
            ref = pb.solve(f); e = rel(P.K.solve(P.rhs(f(pb.axes[0])), "green"), P.K.solve(P.rhs(f(pb.axes[0])), "banded"))
        else:
            F = P.rhs(pb._f(f)); ref = spla.spsolve(P.dense().tocsc(), F.ravel().astype(complex)).reshape(F.shape)
            sl = tuple(slice(p.i0, p.i1 + 1) for p in ((P.Pz, P.Py, P.Px) if d == 3 else (P.Py, P.Px)))
            e = rel(u[sl], ref)
        good = e < 1e-10; ok &= good
        print(f"  {d}D  k = {k!s:12} formule vs resolution directe : {e:.1e}  {'OK' if good else 'ECHEC'}")
pb = gy.problem((2.0,), (400,), ["D", "N"])
kr = pb.resonances(3); kc = (np.arange(3) + 0.5) * np.pi / 2.0
good = np.abs(kr / kc - 1).max() < 1e-4; ok &= good
print(f"  resonances 1D (D-N) vs continu : {np.abs(kr / kc - 1).max():.1e}  {'OK' if good else 'ECHEC'}")
pb = gy.problem((3.0, 3.0), (90, 90), [("D", "D"), ("P", "P")], k=10.0)
pa = gy.problem((3.0, 3.0), (90, 90), [("D", "D"), ("A", "A")], k=10.0)
e = rel(pb.green((1.5, 1.5)), pa.green((1.5, 1.5))); good = e < 5e-2; ok &= good
print(f"  PML (16 mailles) vs transparente exacte, 2D : {e:.1e}  {'OK' if good else 'ECHEC'}")
# garde-fous
k1 = gy.problem((1.0,), (40,), ("D", "D")).resonances(1)[0]
try: gy.problem((1.0,), (40,), ("D", "D"), k=k1); good = False
except gy.SingularOperator: good = True
ok &= good; print(f"  resonance exacte signalee (SingularOperator) : {'OK' if good else 'ECHEC'}")
try: gy.problem((1.0, 1.0), (20, 20), [("D", "D")] * 2, k=3.0).green((0.0, 0.5)); good = False
except ValueError: good = True
ok &= good; print(f"  source sur bord de Dirichlet rejetee : {'OK' if good else 'ECHEC'}")
from gueye import recipes
_, gs = recipes.open_field((3.0, 3.0), (60, 60), 10.0, (1.5, 1.5), pml=16, R0=1e-8, all_open=True)
_, gl = recipes.open_field((9.0, 9.0), (180, 180), 10.0, (4.5, 4.5), pml=16, R0=1e-8, all_open=True)
e = rel(gs, gl[60:121, 60:121]); good = e < 2e-2; ok &= good
print(f"  domaine ouvert (PML sur 4 faces) vs domaine 3x plus grand : {e:.1e}  {'OK' if good else 'ECHEC'}")
print(f"{'TOUT EST OK' if ok else 'ECHEC'} ({time.perf_counter() - t0:.1f} s)")
raise SystemExit(0 if ok else 1)
