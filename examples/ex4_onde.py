"""Onde 2D : impulsion de Ricker dans une boite a CL mixtes, sismogrammes et champ."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import gueye as gy

W = gy.wave(size=(3.0, 2.0), cells=(120, 80), bc=[("D", ("R", 1.5)), ("N", "D")], c=1.0, eta=0.02)
wp = 8.0
ricker = lambda t: (1 - 2 * (0.5 * wp * (t - 1.0))**2) * np.exp(-(0.5 * wp * (t - 1.0))**2)
rec = [(1.0, 1.0), (2.2, 0.5), (2.6, 1.6)]
t, U, R = W.march(0.02, 400, source=W.point_source((0.9, 0.7), ricker), receivers=rec, every=1)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
tr = np.arange(R.shape[0]) * 0.02
for i in range(len(rec)): ax[0].plot(tr, R[:, i] / np.abs(R).max() + 1.6 * i, lw=1)
ax[0].set_xlabel("t"); ax[0].set_yticks([]); ax[0].set_title("signaux aux 3 recepteurs")
for a_, tt in zip(ax[1:], (2.0, 4.5)):
    z = U[int(tt / 0.02)]; v = np.abs(z).max() * 0.8
    a_.imshow(z, origin="lower", extent=[0, 3, 0, 2], cmap="RdBu_r", vmin=-v, vmax=v); a_.set_title(f"t = {tt}", fontsize=9)
plt.tight_layout(); plt.savefig("ex4_onde.png", dpi=120)
