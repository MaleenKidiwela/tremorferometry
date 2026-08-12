import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RD = "fault_tomography/inversion/res_catalog_g50"
d = np.load(f"{RD}/coarse_cm.npz", allow_pickle=True)
wins = d["wins"]; MODEL = d["MODEL"]; idx = d["idx"]; null = float(d["null"]); clat = d["clat"]; clon = d["clon"]; deep = d["deep"]; cm = d["cm_net"]
# reconstruct an example station's families + common mode for the illustration panel
pm = pd.read_parquet(f"{RD}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["t"] = pd.PeriodIndex(pm.ym, freq="M").to_timestamp(); pm["yr"] = pd.PeriodIndex(pm.ym, freq="M").year
ex = pm.tag.value_counts().index[3]; sub = pm[pm.tag == ex]
piv = sub.pivot_table(index="t", columns="cell", values="dvv_month", aggfunc="mean")
cmode = piv.mean(1)

fig = plt.figure(figsize=(15, 9))
# A: example station families + common mode
a = fig.add_axes([0.05, 0.55, 0.42, 0.4])
for c in piv.columns: a.plot(piv.index, piv[c], color="#bbb", lw=.5)
a.plot(cmode.index, cmode.values, color="#c0392b", lw=2, label="COMMON MODE (mean across families)")
a.plot([], [], color="#bbb", lw=.5, label="individual families (they differ → residuals)")
a.set_title(f"station {ex}: families (grey) + extracted common mode (red)"); a.set_ylabel("dv/v (%)"); a.legend(fontsize=8); a.grid(alpha=.3)
# B: residuals for that station (family - common mode)
b = fig.add_axes([0.55, 0.55, 0.42, 0.4])
res = piv.sub(cmode, axis=0)
for c in res.columns: b.plot(res.index, res[c], lw=.6)
b.axhline(0, color="k", lw=.5); b.set_title("residuals = family − common mode (the deep, family-specific signal)"); b.set_ylabel("residual dv/v (%)"); b.grid(alpha=.3)
# C: deep index vs year + null
c = fig.add_axes([0.05, 0.08, 0.42, 0.38])
c.axhspan(-null, null, color="grey", alpha=.2, label=f"noise floor (±{null:.2f}%, 95th of scrambles)")
c.plot(wins[:-1], idx[:-1], "o-", color="#c0392b", lw=1.8, label="deep index (residual, coarse grid)")
c.axhline(0, color="k", lw=.5); c.set_title(f"deep INDEX vs year — max {np.nanmax(np.abs(idx)):.2f}% (borderline; swings 44th–100th pctile w/ method/λ)")
c.set_xlabel("year"); c.set_ylabel("deep δβ/β (%, model)"); c.legend(fontsize=8); c.grid(alpha=.3)
# D: annual deep maps (coarse grid)
yrs = [2015, 2016, 2020, 2021, 2024, 2025]
lo, la = clon[deep], clat[deep]
for k, y in enumerate(yrs):
    ax = fig.add_axes([0.55 + (k % 3)*0.145, 0.28 - (k//3)*0.20, 0.13, 0.17])
    wi = list(wins).index(y); sc = ax.scatter(lo, la, c=MODEL[wi][deep], cmap="RdBu_r", vmin=-0.6, vmax=0.6, s=40, marker="s")
    ax.set_title(f"{y}", fontsize=9); ax.set_xticks([]); ax.set_yticks([])
fig.text(0.77, 0.48, "deep δβ/β maps (0.5° coarse grid, 32 deep cells)", ha="center", fontsize=10)
cax = fig.add_axes([0.985, 0.1, 0.008, 0.3]); plt.colorbar(sc, cax=cax, label="δβ/β (%)")
fig.suptitle("Coarse-grid common-mode-subtracted DEEP inversion (mirror-free): method + index + maps", fontsize=13)
fig.savefig(f"{RD}/coarse_cm_panel.png", dpi=120, bbox_inches="tight")
print("wrote coarse_cm_panel.png")
