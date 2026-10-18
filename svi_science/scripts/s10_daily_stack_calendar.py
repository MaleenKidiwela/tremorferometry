"""Calendar-time version of the family daily-stack plot (user request 2026-10-11).
x = calendar date; y = time after S. Each active day is drawn as a strip WIDTH_D days wide (centred on its date,
so it is visible at a 23-year scale; days within one episode overlap, later days on top). Bottom: repeats per day."""
import sys
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt, matplotlib.dates as mdates
PIL = "/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot"
FAM = int(sys.argv[1]) if len(sys.argv) > 1 else 10
STAS = sys.argv[2].split(",") if len(sys.argv) > 2 else ["CN.PGC", "PB.B011"]
WIDTH_D = 25
rows = pd.read_parquet(f"{PIL}/tables/rows.parquet")
rows = rows[(rows.fam == FAM) & (rows.series == "data") & (rows.P >= 0.5)]
rs = pd.read_parquet(f"{PIL}/tables/row_station.parquet")
fig, axes = plt.subplots(2 * len(STAS), 1, figsize=(16, 6.2 * len(STAS)),
                         gridspec_kw=dict(height_ratios=[4, 1] * len(STAS)), sharex=True)
for j, ns in enumerate(STAS):
    z = np.load(f"{PIL}/win/{ns}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}; tau = z["tau"]
    g = rows.merge(rs[rs.netsta == ns][["row_id", "epoch"]], on="row_id")
    g = g[g.row_id.isin(pos)].reset_index(drop=True); g["date"] = pd.to_datetime(g.date)
    W = z["W"][[pos[r] for r in g.row_id], 0, :].astype(np.float64); W[~np.isfinite(W)] = 0
    ax, axn = axes[2 * j], axes[2 * j + 1]
    days, nrep = [], []
    for d, gg in g.groupby("date"):
        w = gg.P.to_numpy()
        if w.sum() < 3:
            continue
        s = (W[gg.index.to_numpy()] * w[:, None]).sum(0) / w.sum(); s = s / (np.abs(s).max() + 1e-12)
        x0 = mdates.date2num(d - pd.Timedelta(days=WIDTH_D / 2)); x1 = mdates.date2num(d + pd.Timedelta(days=WIDTH_D / 2))
        ax.pcolormesh([x0, x1], np.r_[tau - 0.0125, tau[-1] + 0.0125], s[:, None], cmap="seismic", vmin=-1, vmax=1,
                      shading="flat")
        days.append(d); nrep.append(len(gg))
    ax.axhline(0, color="k", lw=0.6, ls="--"); ax.axhline(2, color="k", lw=0.6); ax.axhline(4, color="k", lw=0.6)
    ax.set_ylim(10, -14); ax.set_ylabel("time after S (s)")
    for e, d0 in g.groupby("epoch").date.min().items():
        ax.axvline(mdates.date2num(d0), color="lime", lw=1.5); ax.text(mdates.date2num(d0), -13.3, f" epoch {e}", color="green", fontsize=9)
    ax.set_title(f"family {FAM} at {ns} (Z, 2–8 Hz): {len(days)} active days on calendar time "
                 f"(each strip drawn {WIDTH_D} days wide so it is visible; empty = no repeats)", fontsize=10)
    axn.bar(days, nrep, width=6, color="k"); axn.set_ylabel("repeats/day"); axn.set_yscale("log")
axes[-1].xaxis.set_major_locator(mdates.YearLocator(2)); axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
axes[-1].set_xlim(pd.Timestamp("2002-10-01"), pd.Timestamp("2026-12-31"))
fig.tight_layout(); out = f"{PIL}/figures/s10_daily_stacks_fam{FAM}_calendar.png"; fig.savefig(out, dpi=120); print(out)
