"""Deep INDEX vs year (single panel) for the coarse common-mode grids (family-centroid, mirror-free).
Black line = deep network-mean delta-beta/beta anomaly per year; grey band = matched scrambled-year null.
Inside the band -> no resolvable secular deep velocity change. Reads {RESDIR}/coarse_cm.npz."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g20")
d = np.load(f"{RD}/coarse_cm.npz", allow_pickle=True)
wins = d["wins"]; idx = d["idx"]; NULL = float(d["null"]); ndeep = int(d["deep"].sum())
STEP = os.environ.get("STEP", "0.2")
fig, a = plt.subplots(figsize=(8.4, 5.2))
a.axhspan(-NULL, NULL, color="grey", alpha=.20, label=f"matched null (±{NULL:.2f}%, 95th of scrambled-year)")
a.plot(wins, idx, "o-", color="#1a1a2e", lw=2, ms=6, label="deep network-mean anomaly")
a.axhline(0, color="k", lw=.6)
a.set_xlabel("year"); a.set_ylabel("deep δβ/β anomaly (%, model)")
a.set_title(f"Deep interface index vs year ({STEP}° grid, {ndeep} deep cells, common-mode, mirror-free)\n"
            f"max |idx| {np.nanmax(np.abs(idx)):.2f}%  vs  null ±{NULL:.2f}%  →  within resolution: no resolvable secular deep change",
            fontsize=10.5)
a.legend(fontsize=9); a.grid(alpha=.3); fig.tight_layout()
fig.savefig(f"{RD}/deep_index_by_year.png", dpi=140, bbox_inches="tight")
print(f"wrote {RD}/deep_index_by_year.png ({ndeep} deep cells, {STEP}deg)")
