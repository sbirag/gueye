import sys, warnings, time; pass
import numpy as np
from gueye.wave1d import Wave1D
warnings.simplefilter("ignore")
L, alpha, c, eta, c1, c0 = 2.0, 1.7, 1.3, 0.25, 0.3, 0.4
# solution manufacturee u = cos(2.1 t + 0.3) * phi(x) + x t^2 ... (CL dependantes du temps)
uex  = lambda x, t: np.cos(2.1 * t + 0.3) * (np.cos(2.3 * x) + 0.4 * x**2 - 0.3 * x) + 0.2 * t**2 * x
ux   = lambda x, t: np.cos(2.1 * t + 0.3) * (-2.3 * np.sin(2.3 * x) + 0.8 * x - 0.3) + 0.2 * t**2
uxx  = lambda x, t: np.cos(2.1 * t + 0.3) * (-2.3**2 * np.cos(2.3 * x) + 0.8)
ut   = lambda x, t: -2.1 * np.sin(2.1 * t + 0.3) * (np.cos(2.3 * x) + 0.4 * x**2 - 0.3 * x) + 0.4 * t * x
utt  = lambda x, t: -2.1**2 * np.cos(2.1 * t + 0.3) * (np.cos(2.3 * x) + 0.4 * x**2 - 0.3 * x) + 0.4 * x
fex  = lambda x, t: utt(x, t) + eta * ut(x, t) - c**2 * uxx(x, t) + c1 * ux(x, t) + c0 * uex(x, t)
BC = {"DD": ("D", "D"), "DN": ("D", "N"), "ND": ("N", "D"), "NN": ("N", "N"), "RR": ("R", "R"), "RN": ("R", "N"),
      "NR": ("N", "R"), "DR": ("D", "R"), "RD": ("R", "D")}
def spec(kind, side):
    # normale sortante : gauche -d/dx, droite +d/dx
    s = -1 if side == "L" else 1; xb = 0.0 if side == "L" else L
    if kind == "D": return ("D", lambda t: uex(xb, t))
    if kind == "N": return ("N", lambda t: s * ux(xb, t))
    return ("R", alpha, lambda t: s * ux(xb, t) + alpha * uex(xb, t))
print("Wave1D : solution manufacturee (CL dependantes du temps), T = 1.5, dt = h/2 -> ordre 2 attendu en (h, dt)")
print(f"{'CL':>3} | err(M=40) | err(80) | err(160) | ordre | banded vs green")
for name, (kl, kr) in BC.items():
    es = []
    for M in (40, 80, 160):
        W = Wave1D(L, M, spec(kl, "L"), spec(kr, "R"), c=c, eta=eta, c1=c1, c0=c0)
        h = L / M; dt = h / 2; n = int(round(1.5 / dt))
        x = W.x
        t, U = W.march(dt, n, f=lambda t: (lambda x: fex(x, t)), u0=lambda x: uex(x, 0.0), v0=lambda x: ut(x, 0.0))
        es.append(np.abs(U[-1] - uex(x, t[-1])).max() / np.abs(uex(x, t[-1])).max())
    W = Wave1D(L, 60, spec(kl, "L"), spec(kr, "R"), c=c, eta=eta, c1=c1, c0=c0)
    f_ = lambda t: (lambda x: fex(x, t))
    _, U1 = W.march(0.02, 40, f=f_, u0=lambda x: uex(x, 0.0), v0=lambda x: ut(x, 0.0), solver="banded")
    _, U2 = W.march(0.02, 40, f=f_, u0=lambda x: uex(x, 0.0), v0=lambda x: ut(x, 0.0), solver="green")
    print(f"{name:>3} | {es[0]:9.2e} | {es[1]:7.2e} | {es[2]:8.2e} | {np.log2(es[1]/es[2]):5.2f} | {np.abs(U1-U2).max()/np.abs(U1).max():.1e}")

print("\nResonateur mixte (D-R), onde pure, impulsion de Ricker : synthese de frequence (formule) vs pas de temps (Newmark)")
Lr, Mr = 3.0, 600
Wr = Wave1D(Lr, Mr, ("D", 0.0), ("R", 1.2, 0.0), c=1.0, eta=0.0)
src, rec = 120, np.array([300, 480])
wp, t0 = 6.0, 1.0
ricker = lambda t: (1 - 2 * (0.5 * wp * (t - t0))**2) * np.exp(-(0.5 * wp * (t - t0))**2)
tt = np.linspace(0, 12, 481)
t_ = time.perf_counter()
Us = Wr.impulse_response(src, rec, tt, wmax=18.0, nw=900, source=ricker)
tsyn = time.perf_counter() - t_
h = Lr / Mr
res = {}
for dt in (0.02, 0.01, 0.005):
    n = int(round(12 / dt))
    fsrc = lambda t: (lambda x: np.where(np.arange(Mr + 1) == src, ricker(t) / h, 0.0))
    t, U, R = Wr.march(dt, n, f=fsrc, receivers=rec)
    Ri = np.array([np.interp(tt, t, R[:, i]) for i in range(2)])
    res[dt] = np.abs(Ri - Us).max() / np.abs(Us).max()
    print(f"   Newmark dt = {dt} : ecart max a la synthese = {res[dt]:.2e}")
dts = sorted(res)[::-1]
print("   ordre observe :", [f"{np.log2(res[dts[i]]/res[dts[i+1]]):.2f}" for i in range(2)], f" ; synthese : {tsyn:.2f} s pour 901 frequences x 2 recepteurs")
