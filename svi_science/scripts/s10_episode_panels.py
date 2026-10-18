"""One subplot per episode: daily dv/v (raw stacks, 2-4 s, -eps) vs day within the episode, per station.
Episodes = runs of active days with gaps <= 20 d. User request 2026-10-11."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
PIL = "/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot"
FAM = int(sys.argv[1]) if len(sys.argv) > 1 else 10
STA = ["CN.PGC", "PB.B011", "PB.B009", "PB.B010", "PB.B926"]
COL = dict(zip(STA, ["k", "tab:blue", "tab:orange", "tab:green", "tab:purple"]))
MK = dict(zip(STA, ["s", "o", "^", "v", "D"]))
d = pd.read_parquet(f"{PIL}/tables/daily.parquet")
d = d[(d.fam == FAM) & d.netsta.isin(STA) & ~d.rail24.astype(bool)].copy()
d["date"] = pd.to_datetime(d.date); d["dvv"] = -d.eps24 * 100; d["sig"] = d.sig24 * 100
days = np.sort(d.date.unique()); ep = np.r_[0, np.cumsum(np.diff(days).astype("timedelta64[D]").astype(int) > 20)]
d["ep"] = d.date.map(dict(zip(days, ep)))
n = d.ep.nunique(); nc = 6; nr = int(np.ceil(n / nc))
fig, axs = plt.subplots(nr, nc, figsize=(3.3 * nc, 2.7 * nr + 0.8), sharey=True, squeeze=False, layout="constrained")
for k, (e, g) in enumerate(d.groupby("ep")):
    ax = axs.flat[k]; t0 = g.date.min()
    for j, ns in enumerate(STA):
        h = g[g.netsta == ns]
        if len(h) == 0: continue
        x = (h.date - t0).dt.days + (j - 2) * 0.08
        sc = ax.scatter(x, h.dvv, c=h.cc24, cmap="viridis", vmin=0.4, vmax=1.0, s=22, marker=MK[ns], edgecolor="k", linewidth=0.3)
    m = np.average(g.dvv, weights=g.cc24.clip(0, None) ** 2)
    ax.axhline(0, color="grey", lw=0.5); ax.axhline(m, color="r", lw=1)
    ax.set_title(f"{t0:%Y-%m-%d}  ({g.date.nunique()} d)  mean {m:+.2f} %", fontsize=8)
    ax.set_ylim(-1.1, 1.1); ax.tick_params(labelsize=7)
for ax in axs.flat[n:]: ax.axis("off")
for ax in axs[:, 0]: ax.set_ylabel("dv/v (%)", fontsize=8)
for ax in axs[-1, :]: ax.set_xlabel("days since episode start", fontsize=8)
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([],[],marker=MK[ns],ls="",color="grey",markeredgecolor="k",label=ns.split(".")[1]) for ns in STA]+[Line2D([],[],color="r",label="episode mean (CC²-weighted)")], loc="outside lower center", ncol=6, fontsize=9)
cb = fig.colorbar(sc, ax=axs, shrink=0.5, pad=0.01, aspect=30); cb.set_label("CC with reference (2–4 s)")
fig.suptitle(f"Family {FAM}: daily dv/v (2–4 s after S, raw day stacks), one panel per episode — colour = CC with reference; red = CC²-weighted episode mean", fontsize=11)
out = f"{PIL}/figures/s10_episode_panels_fam{FAM}.png"; fig.savefig(out, dpi=110); print(out)
