"""Interface graphique interactive de gueye (navigateur, serveur local).

    python -m gueye.gui            # ouvre http://127.0.0.1:<port>/ dans le navigateur
    python -m gueye.gui --port 8765 --no-browser

Aucune dependance en plus de NumPy/SciPy : serveur HTTP de la bibliotheque standard + page HTML/JS. Le serveur n'ecoute que
127.0.0.1 et refuse les requetes dont l'en-tete Host n'est pas local. Les calculs sont ceux de gueye (memes noyaux, memes garde-fous).
"""
import argparse, base64, json, threading, time, warnings, webbrowser
from urllib.parse import urlsplit
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
import numpy as np

from . import api, recipes
from .kernel import SingularOperator
from ._gui_page import PAGE

MAX_CELLS = {1: 20000, 2: 400 * 400, 3: 64 ** 3}
_LOCK = threading.Lock()
_CACHE = {}


def _bc_from_json(b):
    t = str(b.get("type", "D")).upper()
    if t in ("D", "N"): return (t, float(b.get("g", 0.0)))
    if t == "R": return ("R", complex(float(b.get("re", 0.0)), float(b.get("im", 0.0))), float(b.get("g", 0.0)))
    if t == "A": return ("A",)
    if t == "P": return ("P", int(b.get("m", 16)), float(b.get("R0", 1e-6)))
    raise ValueError(f"condition inconnue : {t}")


def _parse(req):
    size = [float(v) for v in req["size"]]; cells = [int(v) for v in req["cells"]]
    dim = len(size)
    if len(cells) != dim or dim not in (1, 2, 3): raise ValueError("dimension invalide")
    if int(np.prod(cells)) > MAX_CELLS[dim]:
        raise ValueError(f"grille trop grande pour l'interface ({dim}D : au plus {MAX_CELLS[dim]} mailles)")
    bc = [(_bc_from_json(a), _bc_from_json(b)) for a, b in req["bc"]]
    return size, cells, bc


def _solver(req):
    """Solveur mis en cache par (geometrie, CL) ; k est applique par with_k (la diagonalisation est reutilisee)."""
    size, cells, bc = _parse(req)
    key = json.dumps([size, cells, req["bc"], req.get("mu", 0.0)], sort_keys=True)
    if key not in _CACHE:
        if len(_CACHE) > 6: _CACHE.clear()
        _CACHE[key] = api.Solver(size, cells, bc, k=0.0, mu=float(req.get("mu", 0.0)), check=False)
    return _CACHE[key].with_k(float(req.get("k", 0.0)), check=True)


def _pack(a):
    return base64.b64encode(np.ascontiguousarray(a, dtype=np.float32).tobytes()).decode()


def _field(sol, req):
    src = [float(v) for v in req["src"]]
    obst = req.get("obstacle")
    if obst and sol.dim == 2:
        cx, cy, r = float(obst["cx"]), float(obst["cy"]), float(obst["r"])
        tot, _ = recipes.dirichlet_obstacle(sol, src, lambda x, y: np.hypot(x - cx, y - cy) < r)
        return np.asarray(tot)
    return np.asarray(sol.green(src))


def handle(path, req):
    t0 = time.perf_counter()
    if path == "/api/solve":
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            sol = _solver(req)
            F = _field(sol, req)
        out = dict(dim=sol.dim, shape=list(F.shape), axes=[list(map(float, a[[0, -1]])) for a in sol.axes],
                   re=_pack(F.real), im=_pack(F.imag), warnings=sorted({str(x.message)[:200] for x in w}),
                   ms=round(1e3 * (time.perf_counter() - t0), 1), h=list(sol.h))
        return out
    if path == "/api/sweep":
        n = max(2, min(int(req.get("n", 200)), 600))
        sol = _solver(dict(req, k=float(req["kmin"])))
        ks = np.linspace(float(req["kmin"]), float(req["kmax"]), n)
        if sol.dim > 1 and any(b[0] in ("P", "A") for pair in sol.bc for b in pair): ks = ks[: min(n, 80)]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            G = sol.sweep(ks, [float(v) for v in req["src"]], [[float(v) for v in req["rec"]]])
        mag = np.abs(G[:, 0])
        res = None
        try:
            r = sol.resonances(12)
            res = [float(np.real(v)) for v in r if kmin_ok(v, ks)]
        except Exception:
            pass
        return dict(ks=[float(v) for v in ks], mag=[None if not np.isfinite(v) else float(v) for v in mag],
                    resonances=res, ms=round(1e3 * (time.perf_counter() - t0), 1))
    if path == "/api/resonances":
        sol = _solver(req)
        r = np.asarray(sol.resonances(int(req.get("n", 10))))
        return dict(re=[float(np.real(v)) for v in r], im=[float(np.imag(v)) for v in r])
    if path == "/api/wave":
        size, cells, bc = _parse(req)
        if len(size) == 1 or int(np.prod(cells)) > (200 * 200 if len(size) == 2 else 40 ** 3):
            raise ValueError("onde : 2D jusqu'a 200x200 ou 3D jusqu'a 40^3")
        W = api.wave(size, cells, bc, c=float(req.get("c", 1.0)), eta=float(req.get("eta", 0.0)))
        f0 = float(req.get("f0", 2.0)); t0r = 1.2 / f0
        ric = lambda t: (1 - 2 * (np.pi * f0 * (t - t0r)) ** 2) * np.exp(-(np.pi * f0 * (t - t0r)) ** 2)
        T = float(req.get("T", 3.0)); nf = max(5, min(int(req.get("frames", 60)), 120))
        h = min(s / c for s, c in zip(size, cells)); dt = 0.5 * h / float(req.get("c", 1.0))
        nsteps = int(np.ceil(T / dt)); every = max(1, nsteps // nf)
        src = W.point_source([float(v) for v in req["src"]], ric)
        res = W.march(dt, nsteps, source=src, every=every); tt, U = res[0], res[1]
        U = np.asarray(U)
        return dict(dim=len(size), shape=list(U.shape), t=[float(v) for v in tt], frames=_pack(U), ms=round(1e3 * (time.perf_counter() - t0), 1))
    raise ValueError("route inconnue")


def kmin_ok(v, ks):
    return abs(np.imag(v)) < 1e-9 and ks[0] <= np.real(v) <= ks[-1]


class _H(BaseHTTPRequestHandler):
    def _local_authority_ok(self, value):
        """Accept only a loopback authority addressed to this server's port."""
        if not value or any(ch in value for ch in "\r\n\\"):
            return False
        try:
            parsed = urlsplit("//" + value)
            if parsed.username is not None or parsed.password is not None or parsed.path or parsed.query or parsed.fragment:
                return False
            host = (parsed.hostname or "").lower()
            port = parsed.port
        except ValueError:
            return False
        server_port = int(self.server.server_address[1])
        if host not in ("127.0.0.1", "localhost"):
            return False
        return port == server_port or (port is None and server_port == 80)

    def _host_ok(self):
        return self._local_authority_ok(self.headers.get("Host", ""))

    def _origin_ok(self):
        # Non-browser clients often omit Origin; if present, it must be local
        # and target this exact server port to mitigate cross-origin requests.
        origin = self.headers.get("Origin")
        if origin is None:
            return True
        try:
            parsed = urlsplit(origin)
            if parsed.scheme != "http" or parsed.username or parsed.password:
                return False
            authority = parsed.netloc
        except ValueError:
            return False
        return self._local_authority_ok(authority)

    def _send(self, code, body, ctype="application/json"):
        b = body if isinstance(body, bytes) else body.encode()
        self.send_response(code); self.send_header("Content-Type", ctype + "; charset=utf-8")
        self.send_header("Content-Length", str(len(b))); self.send_header("Cache-Control", "no-store"); self.end_headers(); self.wfile.write(b)

    def do_GET(self):
        if not self._host_ok(): return self._send(403, "{}")
        if self.path in ("/", "/index.html"): return self._send(200, PAGE, "text/html")
        self._send(404, "{}")

    def do_POST(self):
        if not self._host_ok() or not self._origin_ok(): return self._send(403, "{}")
        try:
            raw_length = self.headers.get("Content-Length", "0")
            n = int(raw_length)
            if n < 0 or n > 1_000_000: raise ValueError("taille de requete invalide ou trop grande")
            if self.headers.get("Transfer-Encoding"):
                raise ValueError("Transfer-Encoding non pris en charge")
            req = json.loads(self.rfile.read(n) or b"{}")
            with _LOCK:
                out = handle(self.path, req)
            self._send(200, json.dumps(out))
        except SingularOperator as e:
            self._send(200, json.dumps(dict(error="Resonance : operateur singulier a cette frequence. " + str(e)[:300], kind="singular")))
        except Exception as e:                                           # message lisible cote interface
            self._send(200, json.dumps(dict(error=f"{type(e).__name__}: {str(e)[:400]}")))

    def log_message(self, *a): pass


def serve(port=0, open_browser=True):
    srv = ThreadingHTTPServer(("127.0.0.1", port), _H)
    url = f"http://127.0.0.1:{srv.server_address[1]}/"
    print(f"gueye GUI : {url}   (Ctrl+C pour arreter)")
    if open_browser: threading.Timer(0.4, lambda: webbrowser.open(url)).start()
    try: srv.serve_forever()
    except KeyboardInterrupt: pass
    finally: srv.server_close()
    return srv


def main(argv=None):
    p = argparse.ArgumentParser(prog="python -m gueye.gui", description="Interface graphique interactive de gueye")
    p.add_argument("--port", type=int, default=0); p.add_argument("--no-browser", action="store_true")
    a = p.parse_args(argv); serve(a.port, not a.no_browser)


if __name__ == "__main__":
    main()
