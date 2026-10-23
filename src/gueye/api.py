"""
API unifiee : on decrit la boite, la grille, les conditions aux limites et le nombre d'onde ; la dimension est deduite.
Toutes les entrees/sorties sont en COORDONNEES PHYSIQUES (x, y, z) ; les tableaux de champ sont ranges (z, y, x).

    import gueye as gy
    pb = gy.problem(size=(3.0, 2.0), cells=(300, 200), k=9.0,
                    bc=[("D", ("R", 3.0)), ("N", "D")])          # axe x : (gauche, droite) ; axe y : (bas, haut)
    u  = pb.solve(lambda x, y: np.exp(-((x - 1)**2 + (y - 1)**2) / 0.01))
    g  = pb.green((0.9, 0.7))                                    # reponse a une source ponctuelle unitaire
    pb.resonances(6)                                             # k_n exacts du schema, tri croissant
    pb.sweep(np.linspace(1, 12, 500), src=(0.9, 0.7), rec=[(2.0, 1.2)])   # reponse en frequence (rapide)

Conditions : 'D' (Dirichlet), 'N' (Neumann), ('R', alpha) (Robin du/dn + alpha u = 0, alpha complexe possible),
('D', g) / ('N', g) / ('R', alpha, g) avec g reel ou fonction des coordonnees de la face, 'A' (transparente exacte),
'P' ou ('P', m, R0) (PML de m mailles, reflexion visee R0, tous axes ; 'A' : dernier axe).
Operateur : -laplacien u + mu^2 u - k^2 u = f   (k complexe = absorption ; mu = masse / ecrantage).
"""
import numpy as np
from .problem1d import Problem1D
from .problem2d import Problem2D, Wave2D
from .problem3d import Problem3D, Wave3D
from .wave1d import Wave1D
from .kernel import SingularOperator

_FACES = {1: [("left", "right")], 2: [("left", "right"), ("bottom", "top")],
          3: [("left", "right"), ("front", "back"), ("bottom", "top")]}


def _norm(b):
    if isinstance(b, str):
        k = b.upper()
        if k in ("D", "N"): return (k, 0.0)
        if k == "A": return ("A",)
        if k == "P": return ("P", 16, 1e-6)
        raise ValueError(f"condition inconnue : {b!r}")
    b = tuple(b); k = str(b[0]).upper()
    if k in ("D", "N"): return (k, b[1] if len(b) > 1 else 0.0)
    if k == "R": return ("R", b[1], b[2] if len(b) > 2 else 0.0)
    if k == "P":
        m = int(b[1]) if len(b) > 1 else 16; R0 = float(b[2]) if len(b) > 2 else 1e-6
        if m < 2 or not (0 < R0 < 1): raise ValueError(f"PML : m >= 2 mailles et 0 < R0 < 1 attendus, recu m={m}, R0={R0}")
        return ("P", m, R0) + tuple(b[3:])
    if k == "A": return ("A",)
    raise ValueError(f"condition inconnue : {b!r}")


def _axis_bcs(bc, dim):
    if dim == 1 and len(bc) == 2:               # 1D : bc = (gauche, droite) accepte directement
        bc = [tuple(bc)]
    if len(bc) != dim:
        raise ValueError(f"bc : {dim} couple(s) (bas, haut) attendus, un par axe (x, y, z), recu {len(bc)}")
    return [(_norm(a), _norm(b)) for a, b in bc]


def _validate(size, cells):
    size = np.atleast_1d(np.asarray(size, dtype=float)); c = np.atleast_1d(np.asarray(cells))
    if size.ndim != 1 or len(size) not in (1, 2, 3) or len(c) != len(size):
        raise ValueError("size et cells doivent avoir le meme nombre d'elements (1, 2 ou 3 : x, y, z)")
    if not np.all(np.isfinite(size)) or np.any(size <= 0):
        raise ValueError(f"size : longueurs strictement positives et finies attendues, recu {tuple(size)}")
    if not np.all(np.isfinite(c.astype(float))) or np.any(c.astype(float) != np.round(c.astype(float))):
        raise ValueError(f"cells : entiers attendus, recu {tuple(c)}")
    c = c.astype(int)
    if np.any(c < 2):
        raise ValueError(f"cells : au moins 2 intervalles par axe (recu {tuple(c)})")
    return tuple(float(v) for v in size), tuple(int(v) for v in c)


class Solver:
    """Probleme de Helmholtz / Poisson / ecrante sur une boite 1D-3D, en coordonnees physiques.

    check : True (defaut) leve SingularOperator si l'operateur est singulier ou quasi singulier (determinant relatif < 1e-10,
    c.-a-d. resonance a ~10 chiffres : le resultat n'aurait plus de sens) et avertit au-dela de 1e-6 ; un nombre = seuil perso ;
    False desactive tout controle (a vos risques)."""

    def __init__(self, size, cells, bc, k=0.0, mu=0.0, check=True):
        size, cells = _validate(size, cells)
        self.dim = len(size)
        self.size, self.cells, self.k, self.mu, self.check = size, cells, k, mu, check
        self.bc = _axis_bcs(bc, self.dim)
        self.c0 = mu**2 - k**2
        kw = {}
        for (lo, hi), (n_lo, n_hi) in zip(self.bc, _FACES[self.dim]):
            kw[n_lo], kw[n_hi] = lo, hi
        if self.dim == 1:
            self.backend = Problem1D(size[0], cells[0], kw["left"], kw["right"], c0=self.c0, check=check)
            self.axes = [self.backend.x]
        elif self.dim == 2:
            self.backend = Problem2D(size[0], size[1], cells[0], cells[1], c0=self.c0, check=check, **kw)
            self.axes = [self.backend.x, self.backend.y]
        else:
            self.backend = Problem3D(size[0], size[1], size[2], cells[0], cells[1], cells[2], c0=self.c0, check=check, **kw)
            self.axes = [self.backend.x, self.backend.y, self.backend.z]
        self.h = tuple(s_ / c for s_, c in zip(size, cells))
        self._kw = kw

    # ------------------------------------------------------------------ coordonnees <-> noeuds
    def _axis_problems(self):
        b = self.backend
        return [b] if self.dim == 1 else ([b.Px, b.Py] if self.dim == 2 else [b.Px, b.Py, b.Pz])

    def node(self, point):
        """Noeud le plus proche d'un point (x[, y[, z]]) -> indice dans l'ordre du champ ((i), (iy, ix) ou (iz, iy, ix)).
        Toute coordonnee strictement hors de la boite physique [0, L] leve ValueError.
        Le point est ensuite associe au noeud de grille le plus proche."""
        p = np.atleast_1d(np.asarray(point, dtype=float))
        if len(p) != self.dim: raise ValueError(f"point a {self.dim} coordonnee(s) attendu")
        for c, L in zip(p, self.size):
            if not np.isfinite(c) or c < 0.0 or c > L:
                raise ValueError(f"point {tuple(p)} hors du domaine physique [0, {self.size}]")
        idx = [int(np.clip(round(c / h), 0, n)) for c, h, n in zip(p, self.h, self.cells)]
        return tuple(idx[::-1]) if self.dim > 1 else idx[0]

    def is_unknown(self, idx):
        """Vrai si le noeud est une inconnue (faux : noeud de Dirichlet elimine, valeur imposee)."""
        ix = (idx,) if self.dim == 1 else tuple(idx)[::-1]          # ordre (x, y, z)
        return all(P.i0 <= i <= P.i1 for P, i in zip(self._axis_problems(), ix))

    def mesh(self):
        """Coordonnees des noeuds (X, Y, Z) en ordre (z, y, x) ; 1D : x."""
        if self.dim == 1: return self.axes[0]
        return tuple(np.meshgrid(*self.axes, indexing="ij")[i].transpose(*range(self.dim)[::-1]) for i in range(self.dim))

    def _f(self, f):
        if callable(f):
            M = self.mesh(); return f(*M) if self.dim > 1 else f(M)
        return f

    # ------------------------------------------------------------------ resolution
    def solve(self, f=0.0):
        """Champ complet pour le second membre f (fonction (x, y, z), tableau nodal, ou scalaire). CL inhomogenes prises en compte."""
        return self.backend.solve(self._f(f))

    def _source_node(self, point):
        i = self.node(point)
        if not self.is_unknown(i):
            raise ValueError(f"source {tuple(np.atleast_1d(point))} sur un noeud de Dirichlet (valeur imposee) : "
                             "placer la source a au moins une maille du bord")
        return i

    def green(self, point):
        """Fonction de Green discrete : reponse a une source ponctuelle unitaire (second membre delta / volume de maille)."""
        i = self._source_node(point)
        return self.backend.green(*i) if self.dim > 1 else self.backend.green(i)

    def at(self, points, f):
        """Valeurs en quelques points sans reconstruire le champ (acces selectif ; pas toujours plus rapide que `solve`).
        Les points sur une face de Dirichlet renvoient la valeur imposee."""
        pts = [self.node(p) for p in points]
        out = np.zeros(len(pts), complex); inner = [j for j, q in enumerate(pts) if self.is_unknown(q)]
        if inner:
            sub = [pts[j] for j in inner]
            v = self.backend.at(np.array(sub) if self.dim == 1 else sub, self._f(f))
            out[inner] = v
        outer = [j for j in range(len(pts)) if j not in inner]
        if outer:
            full = self._boundary_field()
            for j in outer: out[j] = full[pts[j]]
        return out.real if np.abs(out.imag).max() <= 1e-12 * max(np.abs(out.real).max(), 1e-300) else out

    def _boundary_field(self):
        P = self.backend
        if self.dim == 1: return P.expand(np.zeros(P.N))
        shp = (P.Ny, P.Nx) if self.dim == 2 else (P.Nz, P.Ny, P.Nx)
        return P.expand(np.zeros(shp))

    def value(self, rec, src):
        """G(rec ; src) : un element de l'inverse, O(Nx) en 2D, O(Nx Ny) en 3D. Recepteur sur une face de Dirichlet : 0
        (fonction de Green homogene) ; source sur un noeud de Dirichlet : ValueError."""
        s = self._source_node(src); r = self.node(rec)
        if not self.is_unknown(r): return 0.0
        P = self.backend
        if self.dim == 1: return P.K.entries(r - P.i0, s - P.i0) / P.h
        return P.entry(*r, *s)

    def with_k(self, k=None, mu=None, check=None):
        k = self.k if k is None else k; mu = self.mu if mu is None else mu
        chk = self.check if check is None else check
        new = object.__new__(Solver); new.__dict__.update(self.__dict__)
        new.k, new.mu, new.c0, new.check = k, mu, mu**2 - k**2, chk
        new.backend = self.backend.with_c0(new.c0, check=chk)
        return new

    def sweep(self, ks, src, rec):
        """Reponse en frequence : G(rec_j ; src) pour chaque k de `ks` -> tableau (len(ks), len(rec)).
        k singulier (resonance exacte, determinant relatif < seuil) : NaN a cette frequence ; les avertissements de
        quasi-singularite sont regroupes en un seul."""
        import warnings
        ks = np.atleast_1d(ks); rec = list(rec)
        out = np.full((len(ks), len(rec)), np.nan + 0j)
        nearly = 0
        for a, k in enumerate(ks):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                try:
                    s = self.with_k(k)
                    for b, r in enumerate(rec):
                        out[a, b] = s.value(r, src)
                except SingularOperator:
                    pass
                nearly += bool(w)
        if nearly:
            warnings.warn(f"sweep : {nearly} frequence(s) proche(s) d'une resonance (perte de chiffres) ; NaN aux frequences singulieres")
        return out

    def resonances(self, n=8):
        """Premiers k_n tels que -laplacien u = k^2 u avec ces CL (valeurs propres exactes du schema ; axes sans PML/transparente).
        Spectre reel positif : tableau reel. Sinon (Robin complexe ou negatif, operateur non auto-adjoint) : k COMPLEXES renvoyes
        (avec avertissement), tries par partie reelle de la valeur propre."""
        import warnings
        ev = np.zeros(1, complex)
        for (lo, hi), L, M in zip(self.bc, self.size, self.cells):
            if lo[0] in ("A", "P") or hi[0] in ("A", "P"): raise ValueError("resonances : pas de PML / transparente")
            lam = np.linalg.eigvals(Problem1D(L, M, lo, hi, c0=0.0, check=False).dense())
            ev = (ev[:, None] + lam[None, :]).ravel()
            ev = ev[np.argsort(ev.real)[:max(n * 4, 400)]]
        ev = ev[np.argsort(ev.real)][:n]
        scale = max(1.0, np.abs(ev).max())
        if np.all(np.abs(ev.imag) <= 1e-9 * scale) and np.all(ev.real >= -1e-9 * scale):
            return np.sqrt(np.maximum(ev.real, 0.0))
        warnings.warn("spectre non reel positif (Robin complexe ou alpha < 0 : operateur non auto-adjoint) : k complexes renvoyes")
        return np.sqrt(ev)

    def __repr__(self):
        return f"Solver({self.dim}D, size={self.size}, cells={self.cells}, k={self.k}, bc={self.bc})"


class WaveSolver:
    """u_tt + eta u_t - c^2 laplacien u = f(x, t) : pas de temps (Newmark beta = 1/4) en 1D, 2D, 3D."""

    def __init__(self, size, cells, bc, c=1.0, eta=0.0):
        size, cells = _validate(size, cells)
        self.dim = len(size); self.size, self.cells = size, cells
        bcs = _axis_bcs(bc, self.dim)
        kw = {}
        for (lo, hi), (n_lo, n_hi) in zip(bcs, _FACES[self.dim]): kw[n_lo], kw[n_hi] = lo, hi
        if self.dim == 1: self.backend = Wave1D(size[0], cells[0], kw["left"], kw["right"], c=c, eta=eta)
        elif self.dim == 2: self.backend = Wave2D(size[0], size[1], cells[0], cells[1], c=c, eta=eta, **kw)
        else: self.backend = Wave3D(size[0], size[1], size[2], cells[0], cells[1], cells[2], c=c, eta=eta, **kw)
        self.h = tuple(s / n for s, n in zip(size, cells)); self.c = c

    def node(self, point):
        p = np.atleast_1d(np.asarray(point, dtype=float))
        if len(p) != self.dim: raise ValueError(f"point a {self.dim} coordonnee(s) attendu")
        for c, L in zip(p, self.size):
            if not np.isfinite(c) or c < 0.0 or c > L:
                raise ValueError(f"point {tuple(p)} hors du domaine physique [0, {self.size}]")
        idx = [int(np.clip(round(c / h), 0, n)) for c, h, n in zip(p, self.h, self.cells)]
        return tuple(idx[::-1]) if self.dim > 1 else idx[0]

    def point_source(self, point, signal):
        """Terme source ponctuel f(t) = signal(t) delta(x - point) -> fonction du temps compatible avec march."""
        idx = self.node(point); vol = float(np.prod(self.h))
        ax = [self.backend.space] if self.dim == 1 else ([self.backend.P.Px, self.backend.P.Py] if self.dim == 2 else [self.backend.P.Px, self.backend.P.Py, self.backend.P.Pz])
        xi = (idx,) if self.dim == 1 else tuple(idx)[::-1]
        if not all(Pa.i0 <= i <= Pa.i1 for Pa, i in zip(ax, xi)):
            raise ValueError("source sur un noeud de Dirichlet (valeur imposee) : placer la source a au moins une maille du bord")
        P = self.backend.space if self.dim == 1 else self.backend.P
        shape = (self.cells[0] + 1,) if self.dim == 1 else tuple(c + 1 for c in self.cells[::-1])
        def f(t):
            F = np.zeros(shape); F[idx] = signal(t) / vol
            return F
        return f

    def march(self, dt, nsteps, source=None, receivers=None, every=1, **kw):
        recs = None if receivers is None else [self.node(r) for r in receivers]
        if self.dim == 1:
            return self.backend.march(dt, nsteps, f=source, receivers=None if recs is None else np.array(recs), every=every, **kw)
        return self.backend.march(dt, nsteps, f=source, receivers=recs, every=every, **kw)


def problem(size, cells, bc, k=0.0, mu=0.0, check=True):
    """Fabrique : -laplacien u + mu^2 u - k^2 u = f sur une boite ; dimension deduite de len(size). Voir le module pour bc."""
    return Solver(size, cells, bc, k, mu, check)


def wave(size, cells, bc, c=1.0, eta=0.0):
    return WaveSolver(size, cells, bc, c, eta)
