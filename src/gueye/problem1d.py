"""
Problem1D : -c2 u'' + c1 u' + c0 u = f  sur [0, L], M intervalles (h = L/M), 9 combinaisons de conditions aux limites.

Conditions (normale SORTANTE n) :   ('D', g)         u = g
                                    ('N', g)         du/dn = g
                                    ('R', alpha, g)  du/dn + alpha u = g
Dirichlet : noeud de bord elimine. Neumann / Robin : noeud de bord inconnu, point fantome elimine (ordre 2).
Helmholtz  u'' + k^2 u = -f   <=>   Problem1D(L, M, left, right, c2=1, c0=-k**2).
g peut etre un reel ou une fonction du temps (utilise par Wave1D).
"""
import numpy as np
from scipy.linalg import eigvals
from .kernel import TridiagInverse, SingularOperator


def _parse(bc):
    kind = bc[0].upper()
    if kind == "D" or kind == "N":
        return kind, 0.0, bc[1] if len(bc) > 1 else 0.0
    if kind == "R":                       # alpha peut etre complexe (ex. i k : absorption a incidence normale)
        return "R", bc[1], bc[2] if len(bc) > 2 else 0.0
    if kind in ("P", "A"):                # ('P', m, R0[, p]) : PML ; ('A',) : condition transparente exacte
        return kind, 0.0, 0.0
    raise ValueError("condition inconnue : " + str(bc))


def pml_tail(c2, c0, h, k, m, R0, p, shift=0.0):
    """Terme de Schur S d'une couche PML de m noeuds (Dirichlet a l'exterieur) : la ligne du noeud de bord devient d - S.
    s(t) = 1 - i sigma(t)/k, sigma = sigma_max (t/delta)^p, delta = m h, sigma_max = (p+1) ln(1/R0)/(2 delta) (convention exp(+i w t))."""
    delta = m * h
    smax = (p + 1) * np.log(1.0 / R0) / (2 * delta)
    s = lambda t: 1.0 - 1j * smax * (t / delta) ** p / k
    j = np.arange(1, m + 1)
    sj = s((j - 0.5) * h); sp_, sm_ = s(j * h), s((j - 1) * h)
    dj = c2 / (sj * h**2) * (1 / sp_ + 1 / sm_) + c0 + shift
    up = -c2 / (sj * sp_ * h**2)                       # T[j, j+1]
    lo = np.empty(m); lo = -c2 / (np.concatenate([sj[1:], [1.0]]) * sp_ * h**2)   # T[j+1, j]  (dernier inutilise)
    D = dj[-1]
    for i in range(m - 2, -1, -1):
        D = dj[i] - up[i] * lo[i] / D
    return (-c2 / h**2) * (-c2 / (sj[0] * h**2)) / D


def exact_tail(c2, c0, h, shift=0.0):
    """Terme de Schur du demi-espace discret infini (onde sortante exp(+i w t)) : condition transparente EXACTE du schema."""
    a = -c2 / h**2; d = 2 * c2 / h**2 + c0 + shift
    disc = np.sqrt(d * d - 4 * a * a + 0j)
    best = None
    for D in ((d + disc) / 2, (d - disc) / 2):
        rho = -a / D
        if abs(abs(rho) - 1) < 1e-9: ok = rho.imag < 0 or abs(disc) < 1e-9 * abs(d)
        else: ok = abs(rho) < 1
        if ok: best = D
    return a * a / best


class Problem1D:
    def __init__(self, L, M, left=("D", 0.0), right=("D", 0.0), c2=1.0, c1=0.0, c0=0.0, check=True, k_pml=None):
        self.L, self.M, self.h = L, int(M), L / M
        self.left, self.right = left, right
        self.c2, self.c1, self.c0 = c2, c1, c0
        self.x = np.linspace(0, L, self.M + 1)
        (kl, al, self.gL), (kr, ar, self.gR) = _parse(left), _parse(right)
        self.pml_k = None
        self.kl, self.kr = kl, kr
        h = self.h
        d = 2 * c2 / h**2 + c0
        u = -c2 / h**2 + c1 / (2 * h)
        l = -c2 / h**2 - c1 / (2 * h)
        self._tail = {}
        for side, bc, kind in (("L", left, kl), ("R", right, kr)):
            if kind == "P":
                if c1 != 0: raise ValueError("PML : c1 = 0 requis")
                kk = np.sqrt(-(c0 / c2) + 0j).real if k_pml is None else float(np.real(k_pml))
                m, R0 = bc[1], bc[2]; pp = bc[3] if len(bc) > 3 else 2
                if kk * m * (L / M) < 1.0: kk = 1.0 / (m * (L / M))     # k trop petit : PML quasi-statique
                self.pml_k = kk
                self._pmlp = getattr(self, '_pmlp', {}); self._pmlp[side] = (m, R0, pp, kk)
                self._tail[side] = (lambda sh, m=m, R0=R0, pp=pp, kk=kk: pml_tail(c2, c0, L / M, kk, m, R0, pp, sh))
            elif kind == "A":
                self._tail[side] = (lambda sh: exact_tail(c2, c0, L / M, sh))
        self.i0 = 0 if kl != "D" else 1                 # premier noeud inconnu
        self.i1 = self.M if kr != "D" else self.M - 1   # dernier noeud inconnu
        N = self.i1 - self.i0 + 1
        kw = {}
        self._base = {}
        bl = np.zeros(N, complex if np.iscomplexobj(d) else float); br = bl.copy()
        if kl == "D":
            bl[0] = -l
        elif kl in ("P", "A"):
            pass                                         # ligne d'extremite ordinaire + terme de Schur (plus bas)
        else:
            kw.update(d1=d - 2 * h * al * l, u1=u + l); bl[0] = -2 * h * l
        if kr == "D":
            br[-1] = -u
        elif kr in ("P", "A"):
            pass
        else:
            kw.update(dN=d - 2 * h * ar * u, lN=l + u); br[-1] = -2 * h * u
        self.bl, self.br = bl, br
        self._d, self._u, self._l, self._kw = d, u, l, dict(kw)
        if "L" in self._tail: kw["d1"] = d - self._tail["L"](0.0)
        if "R" in self._tail: kw["dN"] = d - self._tail["R"](0.0)
        self.K = TridiagInverse(N, l, d, u, check=check, **kw)       # leve SingularOperator si resonance
        self.N = N
        self.unknown = np.arange(self.i0, self.i1 + 1)

    def extended_dense(self):
        """Matrice EXPLICITE : inconnues physiques + noeuds des couches PML (mur de Dirichlet au-dela). Sert aux axes diagonalises.
        Renvoie (T, nl, nr) : nombres de noeuds de couche a gauche / droite (0 sans PML)."""
        pml = getattr(self, "_pmlp", {})
        if not pml:
            return self.dense(), 0, 0
        h, c2, c0, N = self.h, self.c2, self.c0, self.N
        T0 = self.K.dense().astype(complex)
        nl = pml["L"][0] if "L" in pml else 0
        nr = pml["R"][0] if "R" in pml else 0
        T = np.zeros((nl + N + nr,) * 2, complex); T[nl:nl + N, nl:nl + N] = T0
        for side, b, sgn in (("L", nl, -1), ("R", nl + N - 1, +1)):
            if side not in pml: continue
            m, R0, pp, kk = pml[side]
            delta = m * h; smax = (pp + 1) * np.log(1.0 / R0) / (2 * delta)
            st = lambda t: 1.0 - 1j * smax * (t / delta) ** pp / kk
            j = np.arange(1, m + 1)
            sj, sp_ = st((j - 0.5) * h), st(j * h)
            T[b, b] = self._d                                  # diagonale physique (sans Schur)
            T[b + sgn, b] = -c2 / (sj[0] * h**2)               # noeud 1 <- noeud de bord (s inner = 1)
            T[b, b + sgn] = -c2 / h**2
            for jj in range(1, m + 1):
                g = b + sgn * jj
                T[g, g] = c2 / (sj[jj - 1] * h**2) * (1 / sp_[jj - 1] + 1 / (st((jj - 1) * h))) + c0
                if jj < m:
                    T[g, g + sgn] = -c2 / (sj[jj - 1] * sp_[jj - 1] * h**2)
                    T[g + sgn, g] = -c2 / (sj[jj] * sp_[jj - 1] * h**2)
        return T, nl, nr

    def end_rows(self, shift=0.0):
        """(d1, dN) des lignes d'extremite quand on ajoute `shift` a la diagonale (modes transverses en 2D/3D)."""
        d1 = self._kw.get("d1", self._d) + shift
        dN = self._kw.get("dN", self._d) + shift
        if "L" in self._tail: d1 = self._d + shift - self._tail["L"](shift)
        if "R" in self._tail: dN = self._d + shift - self._tail["R"](shift)
        return d1, dN

    # ------------------------------------------------------------------ outils
    def with_c0(self, c0, check=True):
        """Meme probleme, autre terme d'ordre 0 (balayage en frequence)."""
        return Problem1D(self.L, self.M, self.left, self.right, self.c2, self.c1, c0, check=check)

    def _g(self, g, t):
        return g(t) if callable(g) else g

    def rhs(self, f, t=None):
        """Second membre reduit aux inconnues (f : tableau des M+1 noeuds, fonction de x, ou scalaire)."""
        if callable(f):
            f = f(self.x)
        f = np.broadcast_to(np.asarray(f), self.x.shape)
        return f[self.i0:self.i1 + 1] + self._g(self.gL, t) * self.bl + self._g(self.gR, t) * self.br

    def expand(self, ured, t=None):
        gL = self._g(self.gL, t) if self.kl == "D" else 0.0
        gR = self._g(self.gR, t) if self.kr == "D" else 0.0
        u = np.zeros(self.M + 1, dtype=np.result_type(ured, gL, gR, float))
        u[self.i0:self.i1 + 1] = ured
        if self.kl == "D": u[0] = gL
        if self.kr == "D": u[-1] = gR
        return u

    # ------------------------------------------------------------------ resolution
    def solve(self, f, method="banded", t=None):
        """Champ complet sur les M+1 noeuds. method : 'banded' (O(N)) ou 'green' (formule semi-analytique)."""
        return self.expand(self.K.solve(self.rhs(f, t), method), t)

    def at(self, idx, f, t=None):
        """Valeurs aux noeuds idx par la formule (sans resoudre le champ entier)."""
        idx = np.atleast_1d(idx)
        ured = self.K.apply_at(idx - self.i0, self.rhs(f, t))
        out = np.array(ured, dtype=np.result_type(ured, self.gL, self.gR, float))
        for n, i in enumerate(idx):
            if i < self.i0: out[n] = self._g(self.gL, t)
            if i > self.i1: out[n] = self._g(self.gR, t)
        return out

    def green(self, s):
        """Fonction de Green discrete g(x_i ; x_s) : L g = delta(x - x_s) (second membre e_s / h), CL homogenes."""
        if not (self.i0 <= s <= self.i1):
            raise ValueError("source sur un noeud de Dirichlet")
        g = np.zeros(self.M + 1, dtype=np.result_type(self.K.column(0), float))
        g[self.i0:self.i1 + 1] = self.K.column(s - self.i0) / self.h
        return g

    def dense(self):
        return self.K.dense()

    # ------------------------------------------------------------------ resonances (c1 = 0)
    def resonances(self, n=8):
        """Premiers k_n tels que -u'' = k^2 u (sans c0, c1) avec ces CL (valeurs propres exactes du schema)."""
        if self.c1 != 0:
            raise ValueError("resonances : c1 = 0 requis")
        P = Problem1D(self.L, self.M, self.left, self.right, c2=self.c2, c0=0.0, check=False)
        ev = np.sort_complex(eigvals(P.dense()))
        ev = ev[np.argsort(ev.real)]
        return np.sqrt(ev[:n] / self.c2 + 0j) if np.iscomplexobj(ev) else np.sqrt(ev[:n] / self.c2)
