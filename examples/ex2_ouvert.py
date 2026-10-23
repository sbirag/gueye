"""Rayonnement en milieu ouvert : PML sur l'axe y, parois de Neumann sur x ; comparaison avec la transparente exacte et une cavite fermee."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import gueye as gy
from gueye import recipes

size, cells, k, src = (3.0, 3.0), (150, 150), 14.0, (1.5, 1.5)
print("mailles par longueur d'onde, erreur de dispersion :", recipes.grid_dispersion(k, size[0] / cells[0]))
_, closed = recipes.open_field(size, cells, k, src, pml=2, R0=1e-1)       # quasiment pas d'absorption
_, pml = recipes.open_field(size, cells, k, src, pml=16, R0=1e-8)
_, exact = recipes.open_field(size, cells, k, src, exact=True)
fig, ax = plt.subplots(1, 4, figsize=(18, 4.4))
for a_, g, nm in zip(ax[:3], (closed, pml, exact), ("couche de 2 mailles (peu absorbante)", "PML 16 mailles, R0 = 1e-8", "transparente exacte")):
    v = np.percentile(np.abs(g.real), 99); a_.imshow(g.real, origin="lower", extent=[0, 3, 0, 3], cmap="RdBu_r", vmin=-v, vmax=v); a_.set_title(nm, fontsize=9)
ax[3].plot(np.linspace(0, 3, cells[1] + 1), np.abs(pml - exact)[:, 100] / np.abs(exact).max(), label="|PML - exacte|, x = 2")
ax[3].set_yscale("log"); ax[3].set_xlabel("y"); ax[3].legend(); ax[3].set_title("ecart residuel (reflexion + discretisation)", fontsize=9)
plt.tight_layout(); plt.savefig("ex2_ouvert.png", dpi=120)
print("ecart max PML/exacte :", np.abs(pml - exact).max() / np.abs(exact).max())
