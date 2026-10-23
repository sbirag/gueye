"""
Recettes d'ingenieur (chacune : une fonction, un resultat exploitable).

  resonance_table(pb, n)                        k_n exacts du schema + ecart au continu
  frequency_response(pb, src, rec, ks)          |G(k)| et pics detectes, compares aux resonances
  open_field(size, cells, k, src, ...)          source ponctuelle en milieu ouvert (PML / transparente) sur l'axe ferme
  dirichlet_obstacle(pb, src, inside)           diffusion par un obstacle conducteur (2D), en capacite embarquee (Green)
  grid_dispersion(k, h)                         erreur de dispersion de grille et nombre de mailles par longueur d'onde
"""
import numpy as np
from .api import Solver


def grid_dispersion(k, h):
    """(mailles par longueur d'onde, erreur relative de nombre d'onde de grille) pour -u'' - k^2 u (centre, ordre 2)."""
    kh = k * h
    kd = 2 / h * np.arcsin(kh / 2) if kh < 2 else np.nan
    return 2 * np.pi / kh, (kd / k - 1) if np.isfinite(kd) else np.nan


def resonance_table(pb, n=8):
    """Resonances exactes du schema ; pour une boite a CL D/N seulement, compare aux k du continu."""
    kd = pb.resonances(n)
    kc = []
    for k_ in kd:
        kc.append(np.nan)
    simple = all(lo[0] in "DN" and hi[0] in "DN" for lo, hi in pb.bc)
    if simple:
        modes = []
        grids = []
        for (lo, hi), L in zip(pb.bc, pb.size):
            kind = lo[0] + hi[0]
            f = {"DD": lambda m: m * np.pi / L, "NN": lambda m: m * np.pi / L, "DN": lambda m: (m + 0.5) * np.pi / L, "ND": lambda m: (m + 0.5) * np.pi / L}[kind]
            m0 = 1 if kind == "DD" else 0
            grids.append(np.array([f(m) for m in range(m0, m0 + 40 if kind != "DD" else m0 + 40)]))
        s = np.zeros(1)
        for g in grids: s = (s[:, None] ** 2 + g[None, :] ** 2).ravel() if False else (s[:, None] + g[None, :] ** 2).ravel()
        kc = np.sqrt(np.sort(s))[:n]
    return np.array([kd, np.array(kc) if len(kc) == len(kd) else np.full(len(kd), np.nan)]).T


def frequency_response(pb, src, rec, ks):
    """Reponse en frequence |G(rec ; src ; k)| et pics (indices) ; pb est un Solver."""
    from scipy.signal import find_peaks
    G = pb.sweep(ks, src, rec)
    pk, _ = find_peaks(np.log(np.abs(G[:, 0]) + 1e-300), prominence=0.5)
    return G, np.asarray(ks)[pk]


def open_field(size, cells, k, src, side_bc="N", pml=16, R0=1e-8, exact=False, all_open=False):
    """Source ponctuelle en milieu ouvert. all_open=False : parois `side_bc` sur les axes transverses, PML (ou transparente
    exacte) sur le DERNIER axe ; all_open=True : PML sur TOUTES les faces (la transparente exacte, si exact=True, reste
    limitee au dernier axe). Retourne (Solver, champ de Green)."""
    dim = len(size)
    last = ("A",) if exact else ("P", pml, R0)
    other = ("P", pml, R0) if all_open else side_bc
    bc = [(other, other)] * (dim - 1) + [(last, last)]
    pb = Solver(size, cells, bc, k=k)
    return pb, pb.green(src)


def dirichlet_obstacle(pb, src, inside):
    """Diffusion d'une source ponctuelle par un obstacle conducteur (Dirichlet sur les noeuds ou inside(x, y) est vrai), 2D.
    Capacite embarquee : q = solution de G_SS q = -u_libre|_S, correction = G q (une resolution). Retourne (champ total, noeuds S)."""
    if pb.dim != 2: raise ValueError("dirichlet_obstacle : 2D")
    if not pb.is_unknown(pb.node(src)): raise ValueError("source sur un noeud de Dirichlet")
    P = pb.backend
    X, Y = pb.mesh()
    m = inside(X, Y)
    iy, ix = np.nonzero(m)
    keep = (iy >= P.Py.i0) & (iy <= P.Py.i1) & (ix >= P.Px.i0) & (ix <= P.Px.i1)
    Sy, Sx = iy[keep], ix[keep]
    free = pb.green(src)
    vol = P.hx * P.hy
    GSS = P.block(Sy, Sx, Sy, Sx) * vol                                                   # (T^-1)[S, S]
    q = np.linalg.solve(GSS, -free[Sy, Sx] * vol)
    F = np.zeros((P.Ny, P.Nx), complex); F[Sy - P.Py.i0, Sx - P.Px.i0] = q
    corr = P.expand(P.solve_reduced(F))
    return (free * vol + corr) / vol, (Sy, Sx)


def time_response(pb, src, rec, signal, t, wmax, nw, c=1.0, shift=None):
    """Reponse TEMPORELLE d'un forcage ponctuel  u_tt - c^2 Lap u = signal(t) delta(x - src)  par synthese de Fourier :
    un balayage de frequences complexes w - i*shift (amortissement artificiel exp(-shift t) corrige a la fin).
    Fonctionne avec TOUTES les CL de `pb` (y compris PML / transparente : le domaine ouvert en temps, sans schema d'avancement).
    pb : gy.problem(...) ; signal : fonction du temps ; rec : liste de points ; t : instants ; wmax, nw : bande et nombre de frequences.
    Convention exp(+i w t). Retourne un tableau (len(rec), len(t))."""
    t = np.asarray(t, float); rec = list(rec)
    dw = wmax / nw; w = dw * np.arange(nw + 1); T = 2 * np.pi / dw
    if t.max() >= T: raise ValueError(f"nw trop petit : la fenetre de synthese 2 pi/dw = {T:.3g} doit depasser t max = {t.max():.3g}")
    sh = 3.0 / t.max() if shift is None else shift
    tq = np.arange(0, T, min(np.pi / wmax / 2, 0.02 * t.max()))
    fq = signal(tq) * np.exp(-sh * tq)
    Fh = np.array([np.sum(fq * np.exp(-1j * ww * tq)) * (tq[1] - tq[0]) for ww in w])
    G = pb.sweep((w - 1j * sh) / c, src, rec)                      # (nw+1, nrec), reponse a e_src/h^d pour -Lap - k^2
    if np.isnan(G).any(): raise ValueError("balayage singulier : augmenter `shift`")
    U = np.zeros((len(rec), len(t)))
    for j, ww in enumerate(w):
        wgt = 1.0 if j == 0 else 2.0
        U += wgt * np.real(np.outer(G[j] * Fh[j] / c**2, np.exp(1j * ww * t)))
    return np.exp(sh * t) * U / T
