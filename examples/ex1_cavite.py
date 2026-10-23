"""Cavite 2D a conditions mixtes : resonances exactes du schema et reponse en frequence (balayage de 600 k en ~1 s)."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
import gueye as gy
from gueye import recipes

pb = gy.problem(size=(3.0, 2.0), cells=(120, 80), bc=[("D", ("R", 2.0)), ("N", "D")])
kres = pb.resonances(40)
ks = np.linspace(0.5, 7.0, 600)
G, pics = recipes.frequency_response(pb, src=(0.8, 0.6), rec=[(2.1, 1.3)], ks=ks)
fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
ax[0].semilogy(ks, np.abs(G[:, 0]), lw=0.9)
for k in kres[kres < 7]: ax[0].axvline(k, color="r", lw=0.3)
ax[0].set_xlabel("k"); ax[0].set_ylabel("|G|"); ax[0].set_title(f"Reponse en frequence ; {len(pics)} pics, {np.sum((kres>0.5)&(kres<7))} resonances exactes (rouge)")
for a_, kk, nm in zip(ax[1:], (kres[2] * 1.0005, 4.5), ("juste au-dessus de la resonance n=3", "hors resonance, k = 4.5")):
    s = pb.with_k(kk); g = s.green((0.8, 0.6)).real; v = np.percentile(np.abs(g), 99)
    a_.imshow(g, origin="lower", extent=[0, 3, 0, 2], cmap="RdBu_r", vmin=-v, vmax=v); a_.set_title(f"{nm}\nk = {kk:.3f}", fontsize=9)
plt.tight_layout(); plt.savefig("ex1_cavite.png", dpi=120)
print("resonances :", np.round(kres[:6], 4))
