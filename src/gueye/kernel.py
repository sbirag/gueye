"""
Noyau 1D : inverse exacte d'une matrice tridiagonale 'presque Toeplitz'

        | d1  u1                |
        | l   d   u             |
    T = |     l   d   u         |      (interieur constant (l, d, u), lignes d'extremite libres (d1,u1) et (lN,dN))
        |         ...           |
        |               lN  dN  |

Toutes les conditions aux limites D/N/R (schema a point fantome) donnent une matrice de cette forme.
Formule (Usmani) :  (T^-1)_ij = (-1)^(i+j) u_i..u_{j-1} theta_{i-1} phi_{j+1} / theta_N   (i<=j),  idem avec l (i>j).
Les mineurs theta_i sont obtenus en forme fermee a partir des racines r1, r2 de  x^2 - d x + l u = 0 ;
l'evaluation est faite en variables mises a l'echelle (theta_i = r1^i t_i, |r1|>=|r2|) : aucun debordement,
meme pour N ~ 10^6, et valable pour les regimes sinh (evanescent), sin (propagatif) et racine double (Poisson).
"""
import warnings
import numpy as np
from scipy.linalg import solve_banded


class SingularOperator(Exception):
    """L'operateur discret est singulier (resonance exacte, Neumann pur a k=0, d=0 avec N impair...)."""


def _S(m, q):
    """(1 - q^m)/(1 - q) stable pour q -> 1 (limite : m)."""
    m = np.asarray(m, float)
    x = q - 1.0
    if x == 0:
        return m + 0j
    return np.expm1(m * np.log1p(x)) / x


def _tvec(N, d1, u1, l, d, r1, q):
    """t_i = theta_i / r1^i, i = 0..N-1  (mineurs principaux de tete, ligne 1 = (d1, u1), interieur (l, d, .))."""
    p = r1 * r1 * q                       # = l*u = r1*r2
    i = np.arange(N)
    th1 = d1
    c1 = d * th1 - u1 * l                 # theta_2
    c2 = -p * th1
    t = np.empty(N, complex)
    t[0] = 1.0
    t[1] = th1 / r1
    ii = i[2:]
    t[2:] = c1 * _S(ii - 1, q) / r1**2 + c2 * _S(ii - 2, q) / r1**3
    return t


class TridiagInverse:
    """Inverse exacte de T (voir module).  Parametres : N, interieur (l, d, u), extremites (d1, u1) et (lN, dN)."""

    def __init__(self, N, l, d, u, d1=None, u1=None, lN=None, dN=None, check=True):
        self.N = int(N)
        self.l, self.d, self.u = l, d, u
        self.d1 = d if d1 is None else d1
        self.u1 = u if u1 is None else u1
        self.lN = l if lN is None else lN
        self.dN = d if dN is None else dN
        self.is_real = not any(np.iscomplexobj(np.asarray(v)) and np.iscomplexobj(v) or isinstance(v, complex)
                               for v in (l, d, u, self.d1, self.u1, self.lN, self.dN))
        if self.N < 4:                                        # petites matrices : inversion directe
            self._dense = np.linalg.inv(self.dense())
            self.rcond = 1.0
            return
        self._dense = None
        if l == 0 or u == 0:
            raise ValueError("l et u doivent etre non nuls (sinon : substitution directe).")
        p = complex(l) * complex(u)
        s = np.sqrt(p)
        disc = np.sqrt((d - 2 * s) * (d + 2 * s) + 0j)
        ra, rb = (d + disc) / 2, (d - disc) / 2
        self.r1, r2 = (ra, rb) if abs(ra) >= abs(rb) else (rb, ra)
        self.q = r2 / self.r1
        r1, q, N = self.r1, self.q, self.N
        self.t = _tvec(N, self.d1, self.u1, l, d, r1, q)                     # theta
        self.tp = _tvec(N, self.dN, self.lN, u, d, r1, q)                    # phi (matrice retournee)
        a = self.dN * self.t[N - 1] / r1
        b = u * self.lN * self.t[N - 2] / r1**2
        self.tN = a - b
        tl = max(abs(self.t[N - 1]), abs(self.t[N - 2]))
        scale = (abs(self.dN) / abs(r1) + abs(u * self.lN) / abs(r1) ** 2) * tl
        self.rcond = abs(self.tN) / scale if scale > 0 else 0.0
        tol = 1e-10 if check is True else (float(check) if check else 0.0)      # check=True : seuil 1e-10 ; nombre : seuil perso ; False : aucun
        if check and self.rcond < tol:
            raise SingularOperator(f"operateur singulier ou quasi singulier (determinant relatif {self.rcond:.1e} < {tol:.0e}) : resonance, "
                                   f"Neumann pur a k=0, ou d=0 avec N impair. (check=False pour forcer, a vos risques.)")
        if check and self.rcond < 1e-6:
            warnings.warn(f"proche d'une resonance : determinant relatif {self.rcond:.1e} "
                          f"(perte d'environ {-np.log10(self.rcond):.0f} chiffres)")
        if self.tN == 0:
            self.tN = 1e-300          # (check=False) matrice exactement singuliere : resultat sans signification
        self.wu = -u / r1
        self.wl = -l / r1
        self.rho = self.u1 / u
        self.sig = self.lN / l

    # ------------------------------------------------------------------ matrice
    def dense(self):
        N = self.N
        T = np.diag(np.full(N, self.d, complex)) + np.diag(np.full(N - 1, self.u, complex), 1) \
            + np.diag(np.full(N - 1, self.l, complex), -1)
        T[0, 0], T[N - 1, N - 1] = self.d1, self.dN
        if N > 1:
            T[0, 1], T[N - 1, N - 2] = self.u1, self.lN
        return T.real if self.is_real else T

    def banded(self):
        N = self.N
        ab = np.zeros((3, N), complex)
        ab[1, :] = self.d; ab[0, 1:] = self.u; ab[2, :-1] = self.l
        ab[1, 0], ab[1, -1] = self.d1, self.dN
        ab[0, 1], ab[2, N - 2] = self.u1, self.lN
        return ab.real if self.is_real else ab

    # ------------------------------------------------------------------ elements de l'inverse
    def entries(self, I, J):
        """(T^-1)[I, J] (indices 0-based, tableaux de meme forme ou diffusables)."""
        I, J = np.broadcast_arrays(np.asarray(I), np.asarray(J))
        if self._dense is not None:
            return self._dense[I, J]
        i, j, N = I + 1, J + 1, self.N
        up = i <= j
        out = np.empty(I.shape, complex)
        if up.any():
            iu, ju = i[up], j[up]
            out[up] = (self.wu ** (ju - iu)) * np.where((iu == 1) & (ju > iu), self.rho, 1.0) / self.r1 \
                * self.t[iu - 1] * self.tp[N - ju] / self.tN
        lo = ~up
        if lo.any():
            il, jl = i[lo], j[lo]
            out[lo] = (self.wl ** (il - jl)) * np.where(il == N, self.sig, 1.0) / self.r1 \
                * self.t[jl - 1] * self.tp[N - il] / self.tN
        return out.real if self.is_real else out

    def column(self, j):
        return self.entries(np.arange(self.N), j)

    def row(self, i):
        return self.entries(i, np.arange(self.N))

    def inverse(self):
        I, J = np.meshgrid(np.arange(self.N), np.arange(self.N), indexing="ij")
        return self.entries(I, J)

    # ------------------------------------------------------------------ resolution
    def solve(self, b, method="banded"):
        """T x = b.  'banded' : O(N) (Thomas/LU bande) ; 'green' : x = T^-1 b par la formule (O(N^2))."""
        b = np.asarray(b)
        if method == "green" or self._dense is not None:
            return self.inverse() @ b if self._dense is None else self._dense @ b
        return self._solve_lu(b)

    def _solve_lu(self, b):
        """Resolution tridiagonale avec factorisation LU (gttrf) calculee UNE fois puis reutilisee."""
        if getattr(self, "_lu", None) is None:
            ab = self.banded()
            from scipy.linalg.lapack import get_lapack_funcs
            gttrf, self._gttrs = get_lapack_funcs(("gttrf", "gttrs"), (ab,))
            dl, d, du = ab[2, :-1].copy(), ab[1].copy(), ab[0, 1:].copy()
            dl, d, du, du2, ipiv, info = gttrf(dl, d, du)
            if info != 0: raise SingularOperator("factorisation tridiagonale singuliere")
            self._lu = (dl, d, du, du2, ipiv); self._lu_cplx = np.iscomplexobj(ab)
        dl, d, du, du2, ipiv = self._lu
        def one(v):
            x, info = self._gttrs(dl, d, du, du2, ipiv, v)
            return x
        if self._lu_cplx:
            return one(b.astype(complex))
        if np.iscomplexobj(b):
            return one(b.real.astype(float)) + 1j * one(b.imag.astype(float))
        return one(b.astype(float))

    def apply_at(self, idx, b):
        """Valeurs (T^-1 b)[idx] sans resoudre tout le systeme : O(len(idx) N), acces selectif."""
        idx = np.atleast_1d(idx)
        G = self.entries(idx[:, None], np.arange(self.N)[None, :])
        return G @ np.asarray(b)


class ModeStack:
    """Resolution VECTORISEE de M systemes tridiagonaux (un par mode, meme N) avec la formule fermee :
    u_i = (1/(r1 tN)) [ t_{i-1} U_i + tp_{N-i} wl V_{i-1} ],  U_i = tp_{N-i} f_i + wu U_{i+1},  V_i = wl V_{i-1} + t_{i-1} f_i.
    O(N) iterations de vecteurs de longueur M (stable : |wu|, |wl| <= 1 pour l = u ; memes variables mises a l'echelle que le noyau)."""

    def __init__(self, kernels):
        K0 = kernels[0]
        if any(k._dense is not None for k in kernels) or any(k.N != K0.N for k in kernels):
            raise ValueError("ModeStack : noyaux N >= 4 de meme taille requis")
        self.N, self.M = K0.N, len(kernels)
        self.t = np.stack([k.t for k in kernels], axis=1)       # (N, M)
        self.tp = np.stack([k.tp for k in kernels], axis=1)
        g = lambda name: np.array([getattr(k, name) for k in kernels])
        self.pref = 1.0 / (g("r1") * g("tN")); self.wu, self.wl, self.rho, self.sig = g("wu"), g("wl"), g("rho"), g("sig")
        self.is_real = all(k.is_real for k in kernels)

    def solve(self, F):
        """F : (N, M) -> X : (N, M), colonne m = second membre du mode m."""
        N = self.N
        F = np.asarray(F, dtype=complex)
        U = np.empty((N, self.M), complex); acc = np.zeros(self.M, complex)
        for i in range(N - 1, -1, -1):                          # i 0-based : noeud i+1 ; tp[N-(i+1)] = tp[N-1-i]
            acc = self.tp[N - 1 - i] * F[i] + self.wu * acc
            U[i] = acc
        # cas i = 1 (0-based 0) : les termes j > 1 portent rho
        U[0] = self.tp[N - 1] * F[0] + self.rho * (U[0] - self.tp[N - 1] * F[0])
        X = np.empty((N, self.M), complex); V = np.zeros(self.M, complex)    # V = V_{i-1}
        for i in range(N):
            low = self.tp[N - 1 - i] * self.wl * V
            if i == N - 1: low = low * self.sig
            X[i] = self.pref * (self.t[i] * U[i] + low)
            V = self.wl * V + self.t[i] * F[i]
        return X
