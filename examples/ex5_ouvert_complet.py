"""Domaine ouvert sur TOUTES les faces (PML sur x et y) : diffusion d'une source par un obstacle conducteur, sans paroi parasite."""
import time, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import gueye as gy
from gueye import recipes

pb = gy.problem(size=(4.0, 3.0), cells=(200, 150), k=16.0, bc=[(("P", 16, 1e-8),) * 2, (("P", 16, 1e-8),) * 2])
src = (1.2, 1.0)
t0 = time.perf_counter()
tot, (Sy, Sx) = recipes.dirichlet_obstacle(pb, src, lambda x, y: np.hypot(x - 2.4, y - 1.6) < 0.35)
print(f"{len(Sy)} noeuds d'obstacle, domaine ouvert 4 faces : {time.perf_counter() - t0:.1f} s")
free = pb.green(src)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.4))
for a_, g, nm in zip(ax, (free, tot, tot - free), ("champ libre (ouvert)", "champ total", "champ diffuse")):
    v = np.percentile(np.abs(g.real), 99); a_.imshow(g.real, origin="lower", extent=[0, 4, 0, 3], cmap="RdBu_r", vmin=-v, vmax=v); a_.set_title(nm, fontsize=9)
ax[1].plot(pb.axes[0][Sx], pb.axes[1][Sy], "k.", ms=1.2)
plt.tight_layout(); plt.savefig("ex5_ouvert_complet.png", dpi=120)
