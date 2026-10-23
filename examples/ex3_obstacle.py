"""Diffusion par un obstacle conducteur (Dirichlet) en guide ouvert : capacite embarquee, une seule resolution."""
import time, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import gueye as gy
from gueye import recipes

pb = gy.problem(size=(4.0, 3.0), cells=(240, 180), k=16.0, bc=[("D", "D"), ("P", "P")])
src = (2.0, 0.7)
t0 = time.perf_counter()
tot, (Sy, Sx) = recipes.dirichlet_obstacle(pb, src, lambda x, y: np.hypot(x - 2.0, y - 1.7) < 0.35)
print(f"{len(Sy)} noeuds d'obstacle, calcul total : {time.perf_counter() - t0:.2f} s")
free = pb.green(src)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
for a_, g, nm in zip(ax, (free, tot, tot - free), ("champ libre", "champ total (obstacle de Dirichlet)", "champ diffuse")):
    v = np.percentile(np.abs(g.real), 99); a_.imshow(g.real, origin="lower", extent=[0, 4, 0, 3], cmap="RdBu_r", vmin=-v, vmax=v); a_.set_title(nm, fontsize=9)
ax[1].plot(pb.axes[0][Sx], pb.axes[1][Sy], "k.", ms=1.2)
plt.tight_layout(); plt.savefig("ex3_obstacle.png", dpi=120)
