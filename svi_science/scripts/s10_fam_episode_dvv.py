"""Per-episode dv/v across 2003-2026 for one family: one stack per (station, epoch, episode) from the pilot
(tables/episodes.parquet; stretch vs the station-epoch reference, 2-4 s and 3-5 s after picked S). dv/v = -eps.
No error bars; colour = CC. User request 2026-10-11."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
PIL = "/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot"; FAM = int(sys.argv[1]) if len(sys.argv) > 1 else 10
STAS = ["CN.PGC", "PB.B011", "PB.B009", "PB.B010", "PB.B926", "CN.VGZ", "CN.SNB"]
d = pd.read_parquet(f"{PIL}/tables/episodes.parquet"); d = d[(d.fam == FAM) & d.netsta.isin(STAS)].copy()
d["date"] = pd.to_datetime(d.date_mid); d["dvv24"] = -d.eps24 * 100; d["dvv35"] = -d.eps35 * 100
d.loc[d.rail24.astype(bool), "dvv24"] = np.nan; d.loc[d.rail35.astype(bool), "dvv35"] = np.nan
fig, axs = plt.subplots(len(STAS) + 1, 1, figsize=(13, 2.1 * (len(STAS) + 1)), sharex=True)
for ax, ns in zip(axs, STAS):
    h = d[d.netsta == ns]
    sc = ax.scatter(h.date, h.dvv24, c=h.cc24, cmap="viridis", vmin=0.4, vmax=1, s=45, marker="o", edgecolor="k", lw=0.4)
    ax.scatter(h.date, h.dvv35, c=h.cc35, cmap="viridis", vmin=0.4, vmax=1, s=30, marker="^", edgecolor="grey", lw=0.3)
    for e, he in h.groupby("epoch"):
        ax.axvline(he.date.min() - pd.Timedelta(days=60), color="lime", lw=1)
    ax.axhline(0, color="grey", lw=0.6); ax.set_ylim(-0.6, 0.6); ax.set_ylabel(f"{ns.split('.')[1]}\ndv/v (%)", fontsize=9)
    ax.text(0.005, 0.85, f"{h.dvv24.notna().sum()} episodes; epochs {', '.join(sorted(h.epoch.unique()))}", transform=ax.transAxes, fontsize=8)
# pooled over stations (2-4 s, CC^2-weighted), per episode cluster
p = d.dropna(subset=["dvv24"]).groupby("cluster").apply(lambda g: pd.Series(dict(date=g.date.median(),
        dvv=np.average(g.dvv24, weights=g.cc24.clip(0) ** 2), cc=g.cc24.mean(), nst=g.netsta.nunique())), include_groups=False)
ax = axs[-1]; sc2 = ax.scatter(p.date, p.dvv, c=p.cc, cmap="viridis", vmin=0.4, vmax=1, s=20 + 15 * p.nst, edgecolor="k", lw=0.5)
ax.axhline(0, color="grey", lw=0.6); ax.set_ylim(-0.4, 0.4); ax.set_ylabel("all stations\ndv/v (%)", fontsize=9)
ax.text(0.005, 0.85, "CC²-weighted mean over stations (2–4 s); marker size = number of stations", transform=ax.transAxes, fontsize=8)
fig.colorbar(sc, ax=axs, shrink=0.4, pad=0.01, label="CC with station-epoch reference")
axs[0].set_title(f"Family {FAM}: dv/v per slow-slip episode (one stack per episode per station); ● 2–4 s, ▲ 3–5 s after picked S; green = sensor-epoch start; each epoch its own zero", fontsize=10)
out = f"{PIL}/figures/s10_fam{FAM}_episode_dvv.png"; fig.savefig(out, dpi=115, bbox_inches="tight"); print(out)
p.to_csv(f"{PIL}/tables/fam{FAM}_episode_dvv_pooled.csv")
print(p.round(3).to_string())
