#!/usr/bin/env python
"""Figures for the 4-D dv/v inversion (reads inversion_4d.npz). Anomaly movie + ETS composite + index-vs-year."""
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
OUT = "fault_tomography/inversion/res_catalog"
d = np.load(f"{OUT}/inversion_4d.npz", allow_pickle=True)
wins = d["wins"]; MODEL = d["MODEL"]; SIGM = d["SIGM"]; clat = d["clat"]; clon = d["clon"]; deep = d["deep"]
ets = d["ets_map"]; ets_sig = d["ets_sig"]; idx = d["idx"]; NULL = float(d["null"]); fmed = 0.049
lo, la = clon[deep], clat[deep]
def masked(m, s):
    mm = m.copy().astype(float); mm[(2*s > 0.5) | ~np.isfinite(s)] = np.nan; return mm

# (1) annual anomaly movie
fig, ax = plt.subplots(4, 5, figsize=(16, 12))
for k, y in enumerate(wins):
    a = ax.flat[k]; mm = masked(MODEL[k], SIGM[k])
    sc = a.scatter(lo, la, c=mm[deep], cmap="RdBu_r", vmin=-0.3, vmax=0.3, s=12, marker="s")
    a.set_title(f"{y}", fontsize=10); a.set_xticks([]); a.set_yticks([])
for k in range(len(wins), 20):
    ax.flat[k].axis("off")
fig.suptitle("Deep interface δβ/β ANOMALY by year (v1, raw tensor; masked where 2σ>0.5%). "
             "Anomalies are WITHIN the scrambled-year null → velocity-stable.", fontsize=13, y=0.995)
cax = fig.add_axes([0.93, 0.3, 0.012, 0.4]); plt.colorbar(sc, cax=cax, label="δβ/β anomaly (%)")
fig.tight_layout(rect=[0, 0, 0.92, 0.99]); fig.savefig(f"{OUT}/inversion_annual.png", dpi=120)

# (2) ETS composite + (3) index vs year
fig2, ax2 = plt.subplots(1, 2, figsize=(13, 5))
em = masked(ets, ets_sig)
sc = ax2[0].scatter(lo, la, c=em[deep], cmap="RdBu_r", vmin=-0.3, vmax=0.3, s=22, marker="s")
ax2[0].set_title("ETS-phase COMPOSITE difference (ETS − inter-ETS)\ndeep interface δβ/β, masked"); ax2[0].set_xlabel("lon"); ax2[0].set_ylabel("lat")
plt.colorbar(sc, ax=ax2[0], label="δβ/β (%)")
ax2[1].axhspan(-NULL, NULL, color="grey", alpha=0.25, label=f"scrambled-year null (±{NULL:.2f}%)")
ax2[1].plot(wins, idx, "o-", color="#1a1a2e", label="deep network-mean anomaly")
ax2[1].axhline(0, color="k", lw=0.6)
ax2[1].set_xlabel("year"); ax2[1].set_ylabel("deep δβ/β anomaly (%, model)")
ax2[1].set_title("Deep index vs year\n(within null → no resolvable secular change)"); ax2[1].legend(fontsize=8); ax2[1].grid(alpha=.3)
fig2.tight_layout(); fig2.savefig(f"{OUT}/inversion_ets_index.png", dpi=130)
print("wrote inversion_annual.png, inversion_ets_index.png")
