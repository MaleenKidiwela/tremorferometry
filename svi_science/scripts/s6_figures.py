"""S6 figures: size distributions, size vs time since front arrival, size over the years per group."""
import json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = "/home/jovyan/tremorferometry"; OUT = f"{ROOT}/svi_science/results/s6"; FIG = f"{ROOT}/svi_science/figures"
S = pd.read_parquet(f"{OUT}/sizes_annotated.parquet")
ALL = pd.read_parquet(f"{OUT}/sizes.parquet")
NE = ALL[ALL.kind == "nonevent"]
A = S[(S.P_same >= 0.75) & S.keep_xfam8 & S.dM.notna()]
H = A[A.complete.fillna(False)]
C1, C2, C3, GR = "#2a6f97", "#c8553d", "#7a7a7a", "#bbbbbb"
rng = np.random.default_rng(0)

# ---- 1. size distributions
fig, ax = plt.subplots(1, 3, figsize=(14, 4.2))
b = np.arange(-1.5, 1.6, 0.05)
ax[0].hist(A.dM, b, color=GR, label=f"all P≥0.75 once-counted (n={len(A)})")
ax[0].hist(H.dM, b, color=C1, label=f"above completeness (n={len(H)})")
ax[0].hist(NE.dM.dropna(), b, histtype="step", color=C2, lw=1.5, label=f"non-events P<0.05 (n={NE.dM.notna().sum()})")
ax[0].set_xlabel("ΔM = log10 amplitude / family median"); ax[0].set_ylabel("detections"); ax[0].legend(fontsize=8)
ax[0].set_title("Relative size of each LFE")
b2 = np.arange(0.8, 3.2, 0.05)
ax[1].hist(A.M_approx.dropna(), b2, color=GR, label="all"); ax[1].hist(H.M_approx.dropna(), b2, color=C1, label="above completeness")
ax[1].set_xlabel("approx. magnitude (Bostock family median + ΔM)"); ax[1].set_ylabel("detections"); ax[1].legend(fontsize=8)
ax[1].set_title("Approximate absolute size")
F = pd.read_csv(f"{OUT}/family_sizes.csv")
ax[2].scatter(F.dep, F.mag_med, s=14, color=C3, label="Bostock mag_med")
ax[2].scatter(F.dep, F.M_f_rel + F.mag_med.median(), s=14, color=C1, label="ours: distance-corrected (shifted)")
ax[2].set_xlabel("family depth (km)"); ax[2].set_ylabel("family size (magnitude units)"); ax[2].legend(fontsize=8)
ax[2].set_title("Family size vs depth")
fig.tight_layout(); fig.savefig(f"{FIG}/s6_size_distributions.png", dpi=140); plt.close(fig)

# ---- 2. size vs time since front arrival
L = pd.read_csv(f"{OUT}/t1_lag_rows.csv")
edges = np.array([-3, -1, 0, 1, 2, 3, 5, 7, 10, 15, 25])
fig, ax = plt.subplots(1, 1, figsize=(7.5, 4.5))
for lab, d, col in [("ΔM, above completeness", L[L.complete == True], C1), ("ΔM, all rows", L, GR), ("noise level ΔN (all rows)", L, C2)]:
    y = "dN" if "noise" in lab else "dM"
    mids, med, lo, hi = [], [], [], []
    for a, bb in zip(edges[:-1], edges[1:]):
        v = d[(d.lag_d >= a) & (d.lag_d < bb)][y].dropna().values
        if len(v) < 15: continue
        bs = [np.median(rng.choice(v, len(v))) for _ in range(500)]
        mids.append((a + bb) / 2); med.append(np.median(v)); lo.append(np.percentile(bs, 2.5)); hi.append(np.percentile(bs, 97.5))
    ax.errorbar(mids, med, yerr=[np.array(med) - lo, np.array(hi) - med], marker="o", color=col, label=lab, capsize=2)
ax.axvspan(-1, 2, color="#f2e8cf", zorder=0, label="'at the front' (−1…2 d)")
ax.axhline(0, color="k", lw=0.5)
ax.set_xlabel("days since the tremor front reached the family (10 SVI ETS)"); ax.set_ylabel("median log10 amplitude (family-relative)")
ax.legend(fontsize=8); ax.set_title("LFE size vs position in the episode")
fig.tight_layout(); fig.savefig(f"{FIG}/s6_size_vs_front.png", dpi=140); plt.close(fig)

# ---- 3. over the years per group
G = sorted(H.group.dropna().unique())
fig, axs = plt.subplots(3, 4, figsize=(15, 8.5), sharex=True, sharey=True)
for ax, g in zip(axs.flat, G):
    d = H[H.group == g]; dA = A[A.group == g]
    for dd, col, lab, y in [(d, C1, "ΔM above completeness", "dM"), (dA, C2, "noise ΔN", "dN")]:
        yr = dd.groupby(dd.time_utc.dt.year)[y].agg(["median", "count"]); yr = yr[yr["count"] >= 10]
        ax.plot(yr.index, yr["median"], "o-", color=col, ms=3, label=lab)
    ax.axhline(0, color="k", lw=0.5); ax.set_title(f"group {int(g)} (n={len(d)})", fontsize=9)
    for c in ["2015-01-12", "2017-02-23", "2017-03-03", "2017-08-10", "2019-05-30"]:
        ax.axvline(pd.Timestamp(c).year + pd.Timestamp(c).dayofyear / 365, color=GR, lw=0.5)
axs.flat[0].legend(fontsize=7)
for ax in axs.flat[len(G):]: ax.axis("off")
fig.supylabel("yearly median log10 amplitude (family-relative)"); fig.supxlabel("year (grey lines: station changes)")
fig.suptitle("LFE size over 2014–2026 per family group (≥ 10 events per year)")
fig.tight_layout(); fig.savefig(f"{FIG}/s6_size_by_year.png", dpi=130); plt.close(fig)
print("figures written")
