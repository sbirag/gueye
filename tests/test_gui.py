"""Interface graphique : routes du serveur local (sans navigateur)."""
import base64, json, threading, urllib.request
import numpy as np
import pytest
import gueye as gy
from gueye import gui

D = {"type": "D"}
BC = [[D, D], [D, D]]


def _field(r):
    return np.frombuffer(base64.b64decode(r["re"]), np.float32).reshape(r["shape"]), np.frombuffer(base64.b64decode(r["im"]), np.float32).reshape(r["shape"])


def test_solve_egal_api():
    r = gui.handle("/api/solve", dict(size=[3, 2], cells=[30, 20], bc=BC, k=5.0, src=[0.9, 0.7]))
    re, im = _field(r)
    g = gy.problem(size=(3, 2), cells=(30, 20), k=5.0, bc=[("D", "D"), ("D", "D")]).green((0.9, 0.7))
    assert np.abs(re - g.real).max() < 1e-6 * np.abs(g.real).max()


def test_resonance_signalee():
    k0 = gy.problem(size=(3, 2), cells=(30, 20), k=0.0, bc=[("D", "D"), ("D", "D")]).resonances(1)[0]
    with pytest.raises(gy.SingularOperator):
        gui.handle("/api/solve", dict(size=[3, 2], cells=[30, 20], bc=BC, k=float(k0), src=[0.9, 0.7]))


def test_source_sur_bord_dirichlet_refusee():
    with pytest.raises(ValueError):
        gui.handle("/api/solve", dict(size=[3, 2], cells=[30, 20], bc=BC, k=5.0, src=[0.0, 0.7]))


def test_grille_trop_grande_refusee():
    with pytest.raises(ValueError):
        gui.handle("/api/solve", dict(size=[3, 2], cells=[2000, 2000], bc=BC, k=5.0, src=[0.9, 0.7]))


def test_balayage_et_resonances():
    r = gui.handle("/api/sweep", dict(size=[3, 2], cells=[30, 20], bc=BC, kmin=1, kmax=6, n=40, src=[0.9, 0.7], rec=[2, 1.2]))
    assert len(r["ks"]) == 40 and r["resonances"] and abs(r["resonances"][0] - 1.88) < 0.1


def test_onde_refuse_pml():
    P = {"type": "P", "m": 8, "R0": 1e-4}
    with pytest.raises(NotImplementedError):
        gui.handle("/api/wave", dict(size=[3, 2], cells=[30, 20], bc=[[P, P], [P, P]], src=[1, 1]))


def test_serveur_http_et_host_local():
    srv = gui.ThreadingHTTPServer(("127.0.0.1", 0), gui._H)
    port = srv.server_address[1]
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        page = urllib.request.urlopen(f"http://127.0.0.1:{port}/").read().decode()
        assert "gueye" in page
        body = json.dumps(dict(size=[3, 2], cells=[20, 10], bc=BC, k=3.0, src=[0.9, 0.7])).encode()
        r = json.load(urllib.request.urlopen(urllib.request.Request(f"http://127.0.0.1:{port}/api/solve", body, method="POST")))
        assert "re" in r
        bad = urllib.request.Request(f"http://127.0.0.1:{port}/", headers={"Host": "evil.example"})
        with pytest.raises(urllib.error.HTTPError):
            urllib.request.urlopen(bad)
    finally:
        srv.shutdown(); srv.server_close()


def test_authorite_locale_verifie_port_et_origine():
    srv = gui.ThreadingHTTPServer(("127.0.0.1", 0), gui._H)
    port = srv.server_address[1]
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        import http.client
        conn = http.client.HTTPConnection("127.0.0.1", port, timeout=3)
        conn.request("GET", "/", headers={"Host": f"localhost:{port + 1}"})
        response = conn.getresponse()
        assert response.status == 403
        response.read(); conn.close()

        conn = http.client.HTTPConnection("127.0.0.1", port, timeout=3)
        conn.request("GET", "/", headers={"Host": f"attacker@localhost:{port}"})
        response = conn.getresponse()
        assert response.status == 403
        response.read(); conn.close()

        conn = http.client.HTTPConnection("127.0.0.1", port, timeout=3)
        conn.request("GET", "/", headers={"Host": f"127.0.0.1:{port}"})
        response = conn.getresponse()
        assert response.status == 200
        response.read(); conn.close()
    finally:
        srv.shutdown(); srv.server_close()


def test_post_origine_externe_refusee():
    srv = gui.ThreadingHTTPServer(("127.0.0.1", 0), gui._H)
    port = srv.server_address[1]
    threading.Thread(target=srv.serve_forever, daemon=True).start()
    try:
        body = json.dumps(dict(size=[3, 2], cells=[20, 10], bc=BC, k=3.0, src=[0.9, 0.7])).encode()
        req = urllib.request.Request(
            f"http://127.0.0.1:{port}/api/solve", body, method="POST",
            headers={"Host": f"127.0.0.1:{port}", "Origin": "https://evil.example"})
        with pytest.raises(urllib.error.HTTPError) as exc:
            urllib.request.urlopen(req)
        assert exc.value.code == 403
    finally:
        srv.shutdown(); srv.server_close()
