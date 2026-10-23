"""
Problem3D : -cx u_xx - cy u_yy - cz u_zz + c0 u = f  sur [0,Lx]x[0,Ly]x[0,Lz], D/N/R independants sur les 6 faces.
Axes x et y diagonalises (V_x, V_y), axe z en forme fermee (noyau 1D decale de lambda_x + lambda_y) :
    U = [ sum_(nx,ny) (Tz + lx + ly)^-1 ( F x_x Vx^-1 x_y Vy^-1 ) ] x_x Vx x_y Vy ,    tableaux U[iz, iy, ix].
Faces : left/right (x=0,Lx), front/back (y=0,Ly), bottom/top (z=0,Lz) ; valeur g = reel ou g(a, b) avec les deux
coordonnees le long de la face : (y,z) pour left/right, (x,z) pour front/back, (x,y) pour bottom/top.
"""
import numpy as np
import scipy.sparse as sp
from .problem2d import TimeBC
from .kernel import TridiagInverse, SingularOperator, ModeStack
from .problem1d import Problem1D
from .problem2d import _eigdecomp


def _face(g, a, b):
    A, B = np.meshgrid(a, b, indexing="ij")            # forme (len(a), len(b))
    return (np.asarray(g(A, B), float) if callable(g) else float(g)) * np.ones(A.shape)


class Problem3D:
    def __init__(self, Lx, Ly, Lz, Mx, My, Mz, left=("D", 0.), right=("D", 0.), front=("D", 0.), back=("D", 0.),
                 bottom=("D", 0.), top=("D", 0.), cx=1.0, cy=1.0, cz=1.0, c0=0.0, check=True):
        for _bc in (left, right, front, back):
            if str(_bc[0]).upper() == "A":
                raise NotImplementedError("transparente exacte ('A') : seulement sur le dernier axe ; sur les autres axes utiliser 'P' (PML).")
        self.cx, self.cy, self.cz = cx, cy, cz
        self.Px = Problem1D(Lx, Mx, left, right, c2=cx, check=False, k_pml=float(np.real(np.sqrt(-c0 / cx + 0j))))
        self.Py = Problem1D(Ly, My, front, back, c2=cy, check=False, k_pml=float(np.real(np.sqrt(-c0 / cy + 0j))))
        self.Pz = Problem1D(Lz, Mz, bottom, top, c2=cz, c0=c0, check=False)
        self.x, self.y, self.z = self.Px.x, self.Py.x, self.Pz.x
        self.hx, self.hy, self.hz = self.Px.h, self.Py.h, self.Pz.h
        self.Nx, self.Ny, self.Nz = self.Px.N, self.Py.N, self.Pz.N
        self.shape = (Mz + 1, My + 1, Mx + 1)
        self.sides = dict(left=left, right=right, front=front, back=back, bottom=bottom, top=top)
        self._eig_xy()
        self._cache = {}
        self.modes(0.0, 1.0, check=check)

    def _eig_xy(self):
        Tx, self.oxl, _ = self.Px.extended_dense(); Ty, self.oyl, _ = self.Py.extended_dense()
        self.Nxe, self.Nye = Tx.shape[0], Ty.shape[0]
        self.lx, self.Vx, self.Vxi = _eigdecomp(Tx)
        self.ly, self.Vy, self.Vyi = _eigdecomp(Ty)

    def _pad(self, F):
        Fe = np.zeros((F.shape[0], self.Nye, self.Nxe), dtype=F.dtype)
        Fe[:, self.oyl:self.oyl + self.Ny, self.oxl:self.oxl + self.Nx] = F
        return Fe

    def modes(self, sigma=0.0, mu=1.0, check=True):
        key = (sigma, mu)
        if key in self._cache:
            return self._cache[key]
        K = self.Pz.K
        out = []
        for ix, a in enumerate(self.lx):
            row = []
            for iy, b in enumerate(self.ly):
                lm = a + b
                try:
                    d1m, dNm = self.Pz.end_rows(lm)
                    row.append(TridiagInverse(K.N, mu * K.l, sigma + mu * (K.d + lm), mu * K.u, d1=sigma + mu * d1m,
                                              u1=mu * K.u1, lN=mu * K.lN, dN=sigma + mu * dNm, check=check))
                except SingularOperator as e:
                    raise SingularOperator(f"mode (nu_x={ix}, nu_y={iy}), lambda={lm:.6g} singulier : {e}") from None
            out.append(row)
        self._cache[key] = out
        return out

    def with_c0(self, c0, check=False):
        import copy
        p = copy.copy(self)
        p.Pz = Problem1D(self.Pz.L, self.Pz.M, self.sides["bottom"], self.sides["top"], c2=self.Pz.c2, c0=c0, check=False)
        if getattr(self.Px, "_pmlp", None) or getattr(self.Py, "_pmlp", None):          # PML sur x ou y : depend de k
            S = self.sides
            p.Px = Problem1D(self.Px.L, self.Px.M, S["left"], S["right"], c2=self.cx, check=False, k_pml=float(np.real(np.sqrt(-c0 / self.cx + 0j))))
            p.Py = Problem1D(self.Py.L, self.Py.M, S["front"], S["back"], c2=self.cy, check=False, k_pml=float(np.real(np.sqrt(-c0 / self.cy + 0j))))
            p._eig_xy()
        p._cache = {}
        p.modes(0.0, 1.0, check=check)
        return p

    # ------------------------------------------------------------------ second membre
    def bc_rhs(self):
        Px, Py, Pz, S = self.Px, self.Py, self.Pz, self.sides
        xu, yu, zu = (P.x[P.i0:P.i1 + 1] for P in (Px, Py, Pz))
        G = dict(left=Px.gL, right=Px.gR, front=Py.gL, back=Py.gR, bottom=Pz.gL, top=Pz.gR)
        F = np.zeros((self.Nz, self.Ny, self.Nx), dtype=np.result_type(Pz.bl, float))
        for side, vec in (("left", Px.bl), ("right", Px.br)):          # g(y,z) -> [iz, iy]
            F += _face(G[side], yu, zu).T[:, :, None] * vec[None, None, :]
        for side, vec in (("front", Py.bl), ("back", Py.br)):          # g(x,z) -> [iz, ix]
            F += _face(G[side], xu, zu).T[:, None, :] * vec[None, :, None]
        for side, vec in (("bottom", Pz.bl), ("top", Pz.br)):          # g(x,y) -> [iy, ix]
            F += _face(G[side], xu, yu).T[None, :, :] * vec[:, None, None]
        return F

    def rhs(self, f):
        if callable(f):
            Z, Y, X = np.meshgrid(self.z, self.y, self.x, indexing="ij"); f = f(X, Y, Z)
        f = np.broadcast_to(np.asarray(f), self.shape)
        i = lambda P: slice(P.i0, P.i1 + 1)
        return f[i(self.Pz), i(self.Py), i(self.Px)] + self.bc_rhs()

    def expand(self, U):
        Px, Py, Pz, S = self.Px, self.Py, self.Pz, self.sides
        full = np.zeros(self.shape, dtype=U.dtype)
        full[Pz.i0:Pz.i1 + 1, Py.i0:Py.i1 + 1, Px.i0:Px.i1 + 1] = U
        if Pz.kl == "D": full[0] = _face(Pz.gL, self.x, self.y).T
        if Pz.kr == "D": full[-1] = _face(Pz.gR, self.x, self.y).T
        if Py.kl == "D": full[:, 0, :] = _face(Py.gL, self.x, self.z).T
        if Py.kr == "D": full[:, -1, :] = _face(Py.gR, self.x, self.z).T
        if Px.kl == "D": full[:, :, 0] = _face(Px.gL, self.y, self.z).T
        if Px.kr == "D": full[:, :, -1] = _face(Px.gR, self.y, self.z).T
        return full

    # ------------------------------------------------------------------ resolution
    def _fwd(self, F):
        return np.einsum("ya,zab->zyb", self.Vyi, self._pad(F) @ self.Vxi.T)

    def _bwd(self, Ut):
        U = np.einsum("ya,zab->zyb", self.Vy, Ut @ self.Vx.T)
        return U[:, self.oyl:self.oyl + self.Ny, self.oxl:self.oxl + self.Nx]

    def _stack(self, sigma, mu):
        key = ("stack", sigma, mu)
        if key not in self._cache:
            try: self._cache[key] = ModeStack([k for row in self.modes(sigma, mu) for k in row])
            except ValueError: self._cache[key] = None
        return self._cache[key]

    def solve_reduced(self, F, sigma=0.0, mu=1.0, method="banded"):
        Ft = self._fwd(F); modes = self.modes(sigma, mu)
        S = self._stack(sigma, mu) if method in ("banded", "stack") else None
        if S is not None:                                   # Nxe*Nye modes en un balayage vectorise
            Nz = Ft.shape[0]
            X = S.solve(np.transpose(Ft, (0, 2, 1)).reshape(Nz, -1))      # colonne = ix*Nye + iy
            Ut = np.transpose(X.reshape(Nz, self.Nxe, self.Nye), (0, 2, 1))
            if not (np.iscomplexobj(Ft) or not S.is_real): Ut = Ut.real
            return self._bwd(np.ascontiguousarray(Ut))
        cplx = np.iscomplexobj(Ft) or not modes[0][0].is_real
        Ut = np.empty(Ft.shape, complex if cplx else float)
        for ix in range(self.Nxe):
            for iy in range(self.Nye):
                Ut[:, iy, ix] = modes[ix][iy].solve(Ft[:, iy, ix], method)
        return self._bwd(Ut)

    def solve(self, f, method="banded"):
        return self.expand(self.solve_reduced(self.rhs(f), method=method))

    def green(self, jz, jy, jx):
        Px, Py, Pz = self.Px, self.Py, self.Pz
        a, b, c = jz - Pz.i0, jy - Py.i0, jx - Px.i0
        Ut = np.empty((self.Nz, self.Nye, self.Nxe), complex); modes = self.modes()
        for ix in range(self.Nxe):
            for iy in range(self.Nye):
                Ut[:, iy, ix] = modes[ix][iy].column(a) * self.Vxi[ix, c + self.oxl] * self.Vyi[iy, b + self.oyl]
        U = self._bwd(Ut) / (self.hx * self.hy * self.hz)
        full = np.zeros(self.shape, complex)
        full[Pz.i0:Pz.i1 + 1, Py.i0:Py.i1 + 1, Px.i0:Px.i1 + 1] = U
        return full.real if np.abs(full.imag).max() <= 1e-12 * max(np.abs(full.real).max(), 1e-300) else full

    def entry(self, iz, iy, ix, jz, jy, jx):
        """Un element G[(iz,iy,ix),(jz,jy,jx)] : somme sur les Nx*Ny modes, O(Nx Ny)."""
        Px, Py, Pz = self.Px, self.Py, self.Pz
        modes = self.modes(); acc = 0j
        ox, oy = self.oxl, self.oyl
        for a in range(self.Nxe):
            for b in range(self.Nye):
                acc += (self.Vx[ix - Px.i0 + ox, a] * self.Vxi[a, jx - Px.i0 + ox] * self.Vy[iy - Py.i0 + oy, b] * self.Vyi[b, jy - Py.i0 + oy]
                        * modes[a][b].entries(iz - Pz.i0, jz - Pz.i0))
        acc /= self.hx * self.hy * self.hz
        return acc.real if abs(acc.imag) <= 1e-12 * abs(acc) else acc

    def at(self, pts, f):
        """u aux noeuds (iz,iy,ix) sans reconstruire le champ : somme sur les modes de la ligne z correspondante."""
        Ft = self._fwd(self.rhs(f)); modes = self.modes(); az = np.arange(self.Nz); out = []
        for iz, iy, ix in pts:
            acc = 0j
            for a in range(self.Nxe):
                for b in range(self.Nye):
                    acc += (self.Vx[ix - self.Px.i0 + self.oxl, a] * self.Vy[iy - self.Py.i0 + self.oyl, b]
                            * (modes[a][b].entries(iz - self.Pz.i0, az) @ Ft[:, b, a]))
            out.append(acc)
        out = np.array(out)
        return out.real if np.abs(out.imag).max() <= 1e-12 * max(np.abs(out.real).max(), 1e-300) else out

    def explicit_dense(self):
        """Systeme EXPLICITE (noeuds de PML inclus sur les trois axes) pour verification : (A, (Nze, Nye, Nxe), (offz, offy, offx))."""
        Tz, nlz, _ = self.Pz.extended_dense(); Ty, nly, _ = self.Py.extended_dense(); Tx, nlx, _ = self.Px.extended_dense()
        Ix, Iy, Iz = (sp.identity(T.shape[0]) for T in (Tx, Ty, Tz))
        A = sp.kron(sp.kron(sp.csr_matrix(Tz), Iy), Ix) + sp.kron(sp.kron(Iz, sp.csr_matrix(Ty)), Ix) + sp.kron(sp.kron(Iz, Iy), sp.csr_matrix(Tx))
        return A.tocsc(), (Tz.shape[0], Ty.shape[0], Tx.shape[0]), (nlz, nly, nlx)

    def dense(self):
        if any(str(b[0]).upper() in ('P', 'A') for b in self.sides.values()):
            raise NotImplementedError('dense() : PML / transparente dependent du mode (utiliser des Robin pour comparer au creux)')
        Ix, Iy, Iz = (sp.identity(n) for n in (self.Nx, self.Ny, self.Nz))
        Tx, Ty, Tz = (sp.csr_matrix(P.dense()) for P in (self.Px, self.Py, self.Pz))
        return sp.kron(sp.kron(Tz, Iy), Ix) + sp.kron(sp.kron(Iz, Ty), Ix) + sp.kron(sp.kron(Iz, Iy), Tx)


class Wave3D:
    """u_tt + eta u_t - c^2 Lap u + c0 u = f(x,y,z,t), Newmark beta = 1/4, CL statiques."""

    def __init__(self, Lx, Ly, Lz, Mx, My, Mz, c=1.0, eta=0.0, c0=0.0, bc_time=None, **bc):
        self.P = Problem3D(Lx, Ly, Lz, Mx, My, Mz, cx=c**2, cy=c**2, cz=c**2, c0=c0, check=False, **bc)
        self.eta = eta
        self.tbc = TimeBC(self.P, bc_time, 3)

    def march(self, dt, nsteps, f=None, every=1, receivers=None):
        if any(str(b[0]).upper() in ('P', 'A') for b in self.P.sides.values()):
            raise NotImplementedError('Newmark : pas de PML (frequentiel seulement)')
        P = self.P; beta, gamma = 0.25, 0.5
        a0 = 1 + self.eta * gamma * dt; s = beta * dt**2
        P.modes(a0, s)
        A = P.dense().tocsr(); shp = (P.Nz, P.Ny, P.Nx)
        sl = (slice(P.Pz.i0, P.Pz.i1 + 1), slice(P.Py.i0, P.Py.i1 + 1), slice(P.Px.i0, P.Px.i1 + 1))
        Fbc = P.bc_rhs()
        Z, Y, X = np.meshgrid(P.z, P.y, P.x, indexing="ij")
        def force(t):
            self.tbc.set(t); Fb = P.bc_rhs() if self.tbc.items else Fbc
            ft = 0.0 if f is None else f(t)
            if callable(ft): ft = ft(X, Y, Z)
            return np.broadcast_to(np.asarray(ft, float), P.shape)[sl] + Fb
        Av = lambda w: (A @ w.ravel()).reshape(shp)
        u = np.zeros(shp); v = np.zeros(shp)
        a = force(0.0) - Av(u)
        ts, out = [0.0], [P.expand(u)]
        rec = None if receivers is None else [[0.0 for _ in receivers]]
        for n in range(nsteps):
            t1 = (n + 1) * dt
            ut = u + dt * v + dt**2 * 0.25 * a; vt = v + dt * 0.5 * a
            a = P.solve_reduced(force(t1) - Av(ut) - self.eta * vt, a0, s)
            u = ut + s * a; v = vt + 0.5 * dt * a
            if rec is not None:
                uf = P.expand(u); rec.append([uf[tuple(r)] for r in receivers])
            if (n + 1) % every == 0:
                ts.append(t1); out.append(P.expand(u))
        r = (np.array(ts), np.array(out))
        return r + ((np.array(rec),) if rec is not None else ())
