"""
Wave1D : u_tt + eta u_t - c^2 u_xx + c1 u_x + c0 u = f(x, t)   (onde, telegraphe, Klein-Gordon), 9 CL, g(t) possibles.

Temps : Newmark beta = 1/4, gamma = 1/2 (implicite, inconditionnellement stable, ordre 2). Avec les accelerations comme inconnues,
chaque pas resout  [(1 + eta gamma dt) I + beta dt^2 T] a_{n+1} = ...,  matrice tridiagonale du MEME type (une seule inversion
pour dt constant). Frequence : (T - w^2 + i w eta) u = F, convention exp(+i w t) ; un element par frequence par la formule.
"""
import numpy as np
from .problem1d import Problem1D
from .kernel import TridiagInverse


class Wave1D:
    def __init__(self, L, M, left=("D", 0.0), right=("D", 0.0), c=1.0, eta=0.0, c1=0.0, c0=0.0):
        self.space = Problem1D(L, M, left, right, c2=c**2, c1=c1, c0=c0, check=False)
        self.c, self.eta, self.L, self.M, self.x = c, eta, L, M, self.space.x

    # ------------------------------------------------------------------ domaine temporel
    def march(self, dt, nsteps, f=None, u0=None, v0=None, solver="banded", every=1, receivers=None):
        """Retourne (t, U) avec U[k] = champ complet au temps t[k] (tous les `every` pas), et receivers (noeuds) -> (nsteps+1, nrec)."""
        S = self.space
        if S.kl in ('P', 'A') or S.kr in ('P', 'A'):
            raise NotImplementedError("PML / condition transparente : domaine frequentiel seulement (utiliser impulse_response)")
        beta, gamma = 0.25, 0.5
        K = S.K
        # matrice dynamique : (1 + eta gamma dt) I + beta dt^2 T  (memes extremites modifiees, memes CL)
        a0 = 1 + self.eta * gamma * dt
        s = beta * dt**2
        Kd = TridiagInverse(S.N, s * K.l, a0 + s * K.d, s * K.u,
                            d1=a0 + s * K.d1, u1=s * K.u1, lN=s * K.lN, dN=a0 + s * K.dN)
        T = K.banded()
        def Tmul(v):
            out = T[1] * v
            out[:-1] += T[0, 1:] * v[1:]
            out[1:] += T[2, :-1] * v[:-1]
            return out
        fx = (lambda t: 0.0) if f is None else f
        def force(t):
            ft = fx(t) if callable(fx) else fx
            if callable(ft): ft = ft(self.x)
            return S.rhs(ft, t)
        N = S.N
        u = np.zeros(N); v = np.zeros(N)
        if u0 is not None: u = np.array(u0(self.x) if callable(u0) else u0, float)[S.i0:S.i1 + 1]
        if v0 is not None: v = np.array(v0(self.x) if callable(v0) else v0, float)[S.i0:S.i1 + 1]
        a = force(0.0) - Tmul(u) - self.eta * v
        ts = [0.0]; out = [S.expand(u, 0.0)]
        rec = None if receivers is None else [S.expand(u, 0.0)[receivers]]
        for n in range(nsteps):
            t1 = (n + 1) * dt
            ut = u + dt * v + dt**2 * (0.5 - beta) * a
            vt = v + dt * (1 - gamma) * a
            rhs = force(t1) - Tmul(ut) - self.eta * vt
            a = Kd.solve(rhs, solver)
            u = ut + beta * dt**2 * a
            v = vt + gamma * dt * a
            if rec is not None: rec.append(S.expand(u, t1)[receivers])
            if (n + 1) % every == 0:
                ts.append(t1); out.append(S.expand(u, t1))
        res = (np.array(ts), np.array(out))
        return res + ((np.array(rec),) if rec is not None else ())

    # ------------------------------------------------------------------ domaine frequentiel
    def helmholtz(self, w):
        """Probleme spatial a la pulsation w : (T - w^2 + i w eta) u = F."""
        return self.space.with_c0(self.space.c0 - w**2 + 1j * w * self.eta)

    def transfer(self, w, src, rec):
        """G(rec, src ; w) pour un forçage ponctuel unitaire (rhs e_src / h) : un element de l'inverse par frequence."""
        P = self.helmholtz(w)
        return P.K.entries(np.asarray(rec) - P.i0, src - P.i0) / P.h

    def impulse_response(self, src, rec, t, wmax, nw, shift=None, source=None):
        """Reponse temporelle par synthese de Fourier (frequences w_j = j dw, j = 0..nw), shift = pole complexe  w -> w - i*shift
        (amortissement artificiel exp(-shift t) corrige a la fin).  source(t) : fonction du temps (Ricker...), echantillonnee finement."""
        t = np.asarray(t); rec = np.atleast_1d(rec)
        dw = wmax / nw; w = dw * np.arange(nw + 1)
        T = 2 * np.pi / dw
        sh = (3.0 / (t.max())) if shift is None else shift
        tq = np.arange(0, T, min(np.pi / wmax / 2, 0.02 * t.max()))
        fq = source(tq) * np.exp(-sh * tq)
        Fh = np.array([np.sum(fq * np.exp(-1j * ww * tq)) * (tq[1] - tq[0]) for ww in w])
        U = np.zeros((len(rec), len(t)))
        for j, ww in enumerate(w):
            G = self.transfer(ww - 1j * sh, src, rec)
            Uh = G * Fh[j]
            wgt = 1.0 if j == 0 else 2.0
            U += wgt * np.real(np.outer(Uh, np.exp(1j * ww * t)))
        return np.exp(sh * t) * U / T
