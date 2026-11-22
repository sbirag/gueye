"""
Problem2D : -cx u_xx - cy u_yy + c0 u = f  sur [0,Lx]x[0,Ly], conditions D/N/R independantes sur les 4 cotes.

Methode (generalisation de Gueye 2022) : l'axe x est diagonalise (Tx = V Lam V^-1, symetrisable -> eigh_tridiagonal, O(Nx^2)),
l'axe y est inverse en forme fermee, mode par mode, par le noyau 1D avec decalage lambda_nu :
        U = [ sum_nu  (Ty + lambda_nu)^-1 ( F V^-T )_nu ] V^T,
        G(iy,ix ; jy,jx) = sum_nu V[ix,nu] V^-1[nu,jx] g_nu(iy,jy).
Cotes : left/right (x = 0, Lx), bottom/top (y = 0, Ly) ; valeurs g = reel ou fonction de la coordonnee le long du cote.
Normale sortante, memes conventions qu'en 1D.  Helmholtz : c0 = -k^2.
"""
import numpy as np
import scipy.sparse as sp
from scipy.linalg import eigh_tridiagonal, eig
from .kernel import TridiagInverse, SingularOperator, ModeStack
from .problem1d import Problem1D


def _eigdecomp(T):
    """T = V diag(lam) V^-1 pour une matrice tridiagonale (symetrisable en general)."""
    N = T.shape[0]
    up, lo, dg = np.diag(T, 1), np.diag(T, -1), np.diag(T)
    if not np.iscomplexobj(T) and np.all(up * lo > 0):
        s = np.concatenate([[1.0], np.cumprod(np.sqrt(lo / up))])
        w, Q = eigh_tridiagonal(dg, np.sign(up) * np.sqrt(up * lo))
        return w, s[:, None] * Q, Q.T / s[None, :]
    w, V = eig(T)
    return w, V, np.linalg.inv(V)


def _vals(g, coord):
    """Valeurs de Dirichlet le long d'une face (reelles ou complexes)."""
    values = np.asarray(g(coord) if callable(g) else g)
    return np.broadcast_to(values, np.shape(coord)) * np.ones_like(coord, dtype=np.result_type(values, coord))


def P_sides(P):
    return list(P.sides.values())


class Problem2D:
    def __init__(self, Lx, Ly, Mx, My, left=("D", 0.0), right=("D", 0.0), bottom=("D", 0.0), top=("D", 0.0),
                 cx=1.0, cy=1.0, c0=0.0, check=True):
        for _bc in (left, right):
            if str(_bc[0]).upper() == "A":
                raise NotImplementedError("transparente exacte ('A') : seulement sur le dernier axe ; sur les autres axes utiliser 'P' (PML).")
        self.cx, self.cy = cx, cy
        self.Px = Problem1D(Lx, Mx, left, right, c2=cx, c0=0.0, check=False, k_pml=float(np.real(np.sqrt(-c0 / cx + 0j))))
        self.Py = Problem1D(Ly, My, bottom, top, c2=cy, c0=c0, check=False)
        self.x, self.y = self.Px.x, self.Py.x
        self.hx, self.hy = self.Px.h, self.Py.h
        self.Nx, self.Ny = self.Px.N, self.Py.N
        self.shape = (My + 1, Mx + 1)
        self.sides = dict(left=left, right=right, bottom=bottom, top=top)
        self._eig_x()
        self._cache = {}
        self.modes(0.0, 1.0, check=check)

    def _eig_x(self):
        Tx, nl, nr = self.Px.extended_dense()               # avec PML sur x : couches explicites (dependent de k)
        self.offx, self.Nxe = nl, Tx.shape[0]
        self.lam, self.V, self.Vi = _eigdecomp(Tx)

    def _pad(self, F):
        Fe = np.zeros((F.shape[0], self.Nxe), dtype=F.dtype); Fe[:, self.offx:self.offx + self.Nx] = F
        return Fe

    # ------------------------------------------------------------------ modes
    def modes(self, sigma=0.0, mu=1.0, check=True):
        """Noyaux 1D (en y) de  sigma I + mu (Ty + lam_nu I), nu = 0..Nx-1  (mis en cache)."""
        key = (sigma, mu)
        if key in self._cache:
            return self._cache[key]
        K = self.Py.K
        out = []
        for nu, lm in enumerate(self.lam):
            try:
                d1m, dNm = self.Py.end_rows(lm)
                out.append(TridiagInverse(K.N, mu * K.l, sigma + mu * (K.d + lm), mu * K.u,
                                          d1=sigma + mu * d1m, u1=mu * K.u1, lN=mu * K.lN,
                                          dN=sigma + mu * dNm, check=check))
            except SingularOperator as e:
                raise SingularOperator(f"mode x nu={nu} (lambda={lm:.6g}) singulier : {e}") from None
        self._cache[key] = out
        return out

    # ------------------------------------------------------------------ second membre
    def bc_rhs(self):
        """Contribution (reduite aux inconnues) des valeurs aux limites."""
        Px, Py = self.Px, self.Py
        yu, xu = self.y[Py.i0:Py.i1 + 1], self.x[Px.i0:Px.i1 + 1]
        S = self.sides
        gv = lambda g, c: _vals(g, c)
        F = np.outer(gv(Px.gL, yu), Px.bl) + np.outer(gv(Px.gR, yu), Px.br)
        F = F + np.outer(Py.bl, gv(Py.gL, xu)) + np.outer(Py.br, gv(Py.gR, xu))
        return F

    def rhs(self, f):
        if callable(f):
            Y, X = np.meshgrid(self.y, self.x, indexing="ij"); f = f(X, Y)
        f = np.broadcast_to(np.asarray(f), self.shape)
        return f[self.Py.i0:self.Py.i1 + 1, self.Px.i0:self.Px.i1 + 1] + self.bc_rhs()

    def expand(self, U):
        Px, Py = self.Px, self.Py
        full = np.zeros(self.shape, dtype=U.dtype)
        full[Py.i0:Py.i1 + 1, Px.i0:Px.i1 + 1] = U
        S = self.sides
        if Px.kl == "D": full[:, 0] = _vals(Px.gL, self.y)
        if Px.kr == "D": full[:, -1] = _vals(Px.gR, self.y)
        if Py.kl == "D": full[0, :] = _vals(Py.gL, self.x)
        if Py.kr == "D": full[-1, :] = _vals(Py.gR, self.x)
        return full

    def with_c0(self, c0, check=False):
        """Meme probleme avec un autre c0 (balayage en frequence) : la diagonalisation de x est reutilisee."""
        import copy
        p = copy.copy(self)
        p.Py = Problem1D(self.Py.L, self.Py.M, self.sides["bottom"], self.sides["top"], c2=self.Py.c2, c0=c0, check=False)
        if getattr(self.Px, "_pmlp", None):                  # PML sur x : depend de k -> nouvelle diagonalisation
            p.Px = Problem1D(self.Px.L, self.Px.M, self.sides["left"], self.sides["right"], c2=self.cx, c0=0.0, check=False,
                             k_pml=float(np.real(np.sqrt(-c0 / self.cx + 0j))))
            p._eig_x()
        p._cache = {}
        p.modes(0.0, 1.0, check=check)
        return p

    # ------------------------------------------------------------------ resolution
    def _stack(self, sigma, mu):
        key = ("stack", sigma, mu)
        if key not in self._cache:
            ks = self.modes(sigma, mu)
            try: self._cache[key] = ModeStack(ks)
            except ValueError: self._cache[key] = None
        return self._cache[key]

    def solve_reduced(self, F, sigma=0.0, mu=1.0, method="banded"):
        Ft = self._pad(F) @ self.Vi.T
        S = self._stack(sigma, mu) if method in ("banded", "stack") else None
        if S is not None:                                   # tous les modes en un seul balayage vectorise
            Ut = S.solve(Ft)
            if not (np.iscomplexobj(Ft) or not S.is_real): Ut = Ut.real
        else:
            Ut = np.empty(Ft.shape, dtype=complex if (np.iscomplexobj(Ft) or any(k.is_real is False for k in self.modes(sigma, mu))) else float)
            for nu, K in enumerate(self.modes(sigma, mu)):
                Ut[:, nu] = K.solve(Ft[:, nu], method)
        return (Ut @ self.V.T)[:, self.offx:self.offx + self.Nx]

    def solve(self, f, method="banded"):
        return self.expand(self.solve_reduced(self.rhs(f), method=method))

    def green(self, jy, jx):
        """Fonction de Green discrete (second membre e_(jy,jx)/(hx hy)), CL homogenes, champ complet."""
        Px, Py = self.Px, self.Py
        a, b = jy - Py.i0, jx - Px.i0
        Ut = np.empty((self.Ny, self.Nxe), complex)
        for nu, K in enumerate(self.modes()):
            Ut[:, nu] = K.column(a) * self.Vi[nu, b + self.offx]
        U = (Ut @ self.V.T)[:, self.offx:self.offx + self.Nx] / (self.hx * self.hy)
        full = np.zeros(self.shape, dtype=U.dtype); full[Py.i0:Py.i1 + 1, Px.i0:Px.i1 + 1] = U
        return full.real if np.abs(full.imag).max() <= 1e-12 * max(np.abs(full.real).max(), 1e-300) else full

    def entry(self, iy, ix, jy, jx):
        """G[(iy,ix),(jy,jx)] pour des noeuds inconnus (tableaux diffusables) : O(Nx) par element."""
        Px, Py = self.Px, self.Py
        iy, ix, jy, jx = (np.asarray(v) for v in np.broadcast_arrays(iy, ix, jy, jx))
        out = np.zeros(iy.shape, complex)
        for nu, K in enumerate(self.modes()):
            out += self.V[ix - Px.i0 + self.offx, nu] * self.Vi[nu, jx - Px.i0 + self.offx] * K.entries(iy - Py.i0, jy - Py.i0)
        out = out / (self.hx * self.hy)
        return out.real if np.abs(out.imag).max() <= 1e-12 * max(np.abs(out.real).max(), 1e-300) else out

    def block(self, Sy, Sx, Ty, Tx):
        """G[(Sy,Sx),(Ty,Tx)] : bloc de la fonction de Green (noeuds inconnus) en regroupant les lignes y identiques :
        O(Nx (nlignes x ncolonnes) + Nx |S||T|), bien plus rapide que `entry` pour de gros blocs."""
        Sy, Sx, Ty, Tx = (np.asarray(v) for v in (Sy, Sx, Ty, Tx))
        ry, ia = np.unique(Sy, return_inverse=True); cy, ib = np.unique(Ty, return_inverse=True)
        out = np.zeros((len(Sy), len(Ty)), complex)
        for nu, K in enumerate(self.modes()):
            g = K.entries((ry - self.Py.i0)[:, None], (cy - self.Py.i0)[None, :])[ia][:, ib]
            out += np.outer(self.V[Sx - self.Px.i0 + self.offx, nu], self.Vi[nu, Tx - self.Px.i0 + self.offx]) * g
        out = out / (self.hx * self.hy)
        return out.real if np.abs(out.imag).max() <= 1e-12 * max(np.abs(out.real).max(), 1e-300) else out

    def at(self, pts, f):
        """Valeurs u(iy, ix) en des noeuds donnes (liste de (iy, ix)) sans reconstruire le champ."""
        F = self.rhs(f); Ft = self._pad(F) @ self.Vi.T
        modes = self.modes()
        out = []
        for iy, ix in pts:
            acc = 0j
            for nu, K in enumerate(modes):
                acc += self.V[ix - self.Px.i0 + self.offx, nu] * (K.entries(iy - self.Py.i0, np.arange(self.Ny)) @ Ft[:, nu])
            out.append(acc)
        out = np.array(out)
        return out.real if np.abs(out.imag).max() <= 1e-12 * max(np.abs(out.real).max(), 1e-300) else out

    def dense(self):
        if any(str(b[0]).upper() in ('P', 'A') for b in self.sides.values()):
            raise NotImplementedError('dense() : PML / transparente dependent du mode (utiliser des Robin pour comparer au creux)')
        """Matrice 2D complete (verification ; invalide avec PML / transparente, dependantes du mode) : Ty (x) I + I (x) Tx, indice n = iy*Nx + ix."""
        return sp.kron(sp.csr_matrix(self.Py.dense()), sp.identity(self.Nx)) + sp.kron(sp.identity(self.Ny), sp.csr_matrix(self.Px.dense()))

    def explicit_dense(self):
        """Systeme EXPLICITE complet (noeuds de PML inclus sur les deux axes) pour verification : renvoie (A, Nye, Nxe, offy, offx)."""
        Ty, nly, nry = self.Py.extended_dense(); Tx, nlx, nrx = self.Px.extended_dense()
        A = sp.kron(sp.csr_matrix(Ty), sp.identity(Tx.shape[0])) + sp.kron(sp.identity(Ty.shape[0]), sp.csr_matrix(Tx))
        return A.tocsc(), Ty.shape[0], Tx.shape[0], nly, nlx

    def resonances(self, n=8):
        """k_n (c0 hors jeu) : racines de  lambda_x + mu_y ; valeurs propres du schema, tri croissant (k^2 = lam + mu)."""
        mu = np.linalg.eigvals(Problem1D(self.Py.L, self.Py.M, self.sides["bottom"], self.sides["top"], c2=self.Py.c2, check=False).dense())
        s = (self.lam[:, None] + mu[None, :]).ravel()
        s = np.sort(s.real)[:n]
        return np.sqrt(s.astype(complex))


class Wave2D:
    """u_tt + eta u_t - cx u_xx - cy u_yy + c0 u = f(x,y,t) ; Newmark beta = 1/4 ; CL statiques (g(coordonnee)) ; une seule factorisation modale."""

    def __init__(self, Lx, Ly, Mx, My, left=("D", 0.0), right=("D", 0.0), bottom=("D", 0.0), top=("D", 0.0),
                 c=1.0, eta=0.0, c0=0.0, bc_time=None):
        self.P = Problem2D(Lx, Ly, Mx, My, left, right, bottom, top, cx=c**2, cy=c**2, c0=c0, check=False)
        self.eta = eta
        self.tbc = TimeBC(self.P, bc_time, 2)
        self.x, self.y = self.P.x, self.P.y

    def march(self, dt, nsteps, f=None, u0=None, v0=None, every=1, receivers=None):
        if any(str(b[0]).upper() in ('P', 'A') for b in P_sides(self.P)):
            raise NotImplementedError('Newmark : pas de PML (frequentiel seulement)')
        P = self.P; beta, gamma = 0.25, 0.5
        a0 = 1 + self.eta * gamma * dt; s = beta * dt**2
        P.modes(a0, s)
        A = P.dense().tocsr()
        iy0, ix0 = P.Py.i0, P.Px.i0
        shp = (P.Ny, P.Nx)
        sel = lambda a: np.asarray(a(*np.meshgrid(self.x, self.y)[::1]) if callable(a) else a, float)[iy0:iy0 + P.Ny, ix0:ix0 + P.Nx]
        u = np.zeros(shp); v = np.zeros(shp)
        if u0 is not None: u = sel(u0).copy()
        if v0 is not None: v = sel(v0).copy()
        Fbc = P.bc_rhs()
        Y, X = np.meshgrid(self.y, self.x, indexing="ij")
        def force(t):
            self.tbc.set(t); Fb = P.bc_rhs() if self.tbc.items else Fbc
            if f is None: ft = 0.0
            else:
                ft = f(t)
                if callable(ft): ft = ft(X, Y)
            return np.broadcast_to(np.asarray(ft, float), P.shape)[iy0:iy0 + P.Ny, ix0:ix0 + P.Nx] + Fb
        Av = lambda w: (A @ w.ravel()).reshape(shp)
        a = force(0.0) - Av(u) - self.eta * v
        ts, out = [0.0], [P.expand(u)]
        rec = None if receivers is None else [[P.expand(u)[r] for r in receivers]]
        for n in range(nsteps):
            t1 = (n + 1) * dt
            ut = u + dt * v + dt**2 * (0.5 - beta) * a
            vt = v + dt * (1 - gamma) * a
            a = P.solve_reduced(force(t1) - Av(ut) - self.eta * vt, a0, s)
            u = ut + beta * dt**2 * a; v = vt + gamma * dt * a
            if rec is not None:
                uf = P.expand(u); rec.append([uf[tuple(r)] for r in receivers])
            if (n + 1) % every == 0:
                ts.append(t1); out.append(P.expand(u))
        r = (np.array(ts), np.array(out))
        return r + ((np.array(rec),) if rec is not None else ())


class TimeBC:
    """Conditions aux limites dependant du temps : g(x, t) = g0(x) * s(t) par face (bc_time={'left': s, ...}).
    s : fonction du temps ; la valeur spatiale g0 est celle donnee dans la CL. Applique g0*s(t) dans bc_rhs et expand."""
    _AXES = {"left": ("Px", "gL"), "right": ("Px", "gR"), "front": ("Py", "gL"), "back": ("Py", "gR"),
             "bottom": ("Py", "gL"), "top": ("Py", "gR")}

    def __init__(self, P, bc_time, dim):
        amap = dict(TimeBC._AXES)
        if dim == 3: amap.update(bottom=("Pz", "gL"), top=("Pz", "gR"))
        bad = [k for k in (bc_time or {}) if k not in amap]
        if bad: raise ValueError(f"faces inconnues dans bc_time : {bad}")
        self.items = []
        for side, s in (bc_time or {}).items():
            if not callable(s): raise ValueError("bc_time : une fonction du temps par face")
            ax, attr = amap[side]; obj = getattr(P, ax)
            self.items.append((obj, attr, getattr(obj, attr), s))

    def set(self, t):
        for obj, attr, g0, s in self.items:
            st = float(s(t))
            if callable(g0): setattr(obj, attr, (lambda *c, g0=g0, st=st: np.asarray(g0(*c), float) * st))
            else: setattr(obj, attr, float(g0) * st)
