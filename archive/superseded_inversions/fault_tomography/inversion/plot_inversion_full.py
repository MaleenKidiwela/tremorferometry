#!/usr/bin/env python
"""Full inversion figure panel from inversion_4d.npz (era-split v1). Annual delta-beta/beta maps (deep cells),
posterior sigma_m map, ETS composite map, detectability bound. Honest: shows the un-thresholded model with a
sigma_m mask overlay (thresholded maps are empty at v1 sensitivity), per Merlin."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "fault_tomography/inversion/res_catalog"
d = np.load(f"{OUT}/inversion_4d.npz", allow_pickle=True)
wins = d["wins"]; MODEL = d["MODEL"]; SIGM = d["SIGM"]; ets = d["ets_map"]; clat = d["clat"]; clon = d["clon"]
deep = d["deep"]; idx = d["idx"]; bk = d["bound_keys"]; bv = d["bound_vals"]
lo, la = clon[deep], clat[deep]
fin = [int(np.isfinite(MODEL[i][deep]).sum()) for i in range(len(wins))]
print("finite deep cells per year:", fin)

# panel 1: annual maps grid (2011-2025, 15 panels) + ETS composite + sigma map
yrs = list(range(2011, 2026))
fig, ax = plt.subplots(3, 6, figsize=(19, 9))
vlim = 0.5
for k, y in enumerate(yrs):
    wi = list(wins).index(y); a = ax.flat[k]; m = MODEL[wi]
    sc = a.scatter(lo, la, c=m[deep], cmap="RdBu_r", vmin=-vlim, vmax=vlim, s=11, marker="s", edgecolors="none")
    a.set_title(f"{y}  (idx {idx[wi]:+.2f}%)", fontsize=9); a.set_xticks([]); a.set_yticks([])
# ETS composite
a = ax.flat[15]; sc2 = a.scatter(lo, la, c=ets[deep], cmap="RdBu_r", vmin=-vlim, vmax=vlim, s=13, marker="s")
a.set_title("ETS composite\n(≈0)", fontsize=9); a.set_xticks([]); a.set_yticks([])
# sigma_m map (recent)
a = ax.flat[16]; smrec = np.nanmedian(np.dstack([SIGM[list(wins).index(y)] for y in [2023, 2024, 2025]]), 2)
sc3 = a.scatter(lo, la, c=smrec[deep], cmap="viridis_r", vmin=0, vmax=2, s=13, marker="s")
a.set_title("posterior σ_m\n(2023-25, model %)", fontsize=9); a.set_xticks([]); a.set_yticks([]); plt.colorbar(sc3, ax=a, fraction=.05)
# bound text
a = ax.flat[17]; a.axis("off")
txt = "DETECTABILITY BOUND\n(smallest coherent deep patch\nrecoverable above noise)\n\n"
for k, key in enumerate(bk):
    T, n95, Amin, Amd = bv[k]; txt += f"{key}: {Amin:.2f}% model / {Amd*1000:.1f} milli-% data\n"
txt += f"\nindex max {np.nanmax(np.abs(idx)):.2f}% @ {float(d['idx_pctile']):.0f}th pctile\nETS composite ≈ 0\n(v1 raw tensor — noise-field\nuncorrected; see mirror test)"
a.text(0, .5, txt, fontsize=9, va="center", family="monospace")
fig.suptitle("4-D deep-interface δβ/β inversion (era-split v1, raw tensor) — annual maps + ETS composite + σ_m + bound", fontsize=13, y=1.0)
cax = fig.add_axes([1.005, 0.4, 0.01, 0.4]); plt.colorbar(sc, cax=cax, label="δβ/β (%, model)")
fig.tight_layout(); fig.savefig(f"{OUT}/inversion_full_panel.png", dpi=120, bbox_inches="tight")
print("wrote inversion_full_panel.png")
