"""Daily stack plot for the best repeating family (user request 2026-10-11).

Same band as the original single-station repeating-family dv/v (densify_gnw_gpu.py: 2-8 Hz, 40 Hz).
Pilot windows are already 2-8 Hz, 40 Hz, despiked, RMS-normalised, tau = time relative to predicted S.
For each active day (data rows, P >= 0.5, counted once): P-weighted mean of the day's windows (Z).
No realignment, no dv/v fitting: the raw daily stacks as they are.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PIL = "/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot"
FAM = int(sys.argv[1]) if len(sys.argv) > 1 else 10
STAS = sys.argv[2].split(",") if len(sys.argv) > 2 else ["PB.B011", "CN.PGC"]
XL0, XL1 = -1, 5
OUT = f"{PIL}/figures/s10_daily_stacks_fam{FAM}_m1to5.png"

SER = pd.read_csv(f"{PIL}/tables/series.csv")
rows = pd.read_parquet(f"{PIL}/tables/rows.parquet")
rows = rows[(rows.fam == FAM) & (rows.series == "data") & (rows.P >= 0.5)]
rs = pd.read_parquet(f"{PIL}/tables/row_station.parquet")
ets = pd.read_csv("/home/jovyan/tremorferometry/svi_science/results/backward_scope/truth_episodes.csv",
                  parse_dates=["start", "end"])

fig, axes = plt.subplots(2, len(STAS), figsize=(7.5 * len(STAS), 13),
                         gridspec_kw=dict(height_ratios=[1, 4]), squeeze=False)
for j, ns in enumerate(STAS):
    z = np.load(f"{PIL}/win/{ns}.npz")
    pos = {r: i for i, r in enumerate(z["row_id"])}
    tau = z["tau"]
    g = rows.merge(rs[rs.netsta == ns][["row_id", "epoch"]], on="row_id")
    g = g[g.row_id.isin(pos)].reset_index(drop=True)
    W = z["W"][[pos[r] for r in g.row_id], 0, :].astype(np.float64)  # Z component
    W[~np.isfinite(W)] = 0.0
    g["date"] = pd.to_datetime(g.date)
    days, stacks, nrep = [], [], []
    for d, gg in g.groupby("date"):
        idx = gg.index.to_numpy()
        w = gg.P.to_numpy()
        if w.sum() < 3:          # same active-day rule as the pilot (sum P >= 3)
            continue
        s = (W[idx] * w[:, None]).sum(0) / w.sum()
        days.append(d); stacks.append(s); nrep.append(len(idx))
    S = np.array(stacks)
    grand = (W * g.P.to_numpy()[:, None]).sum(0) / g.P.sum()
    tp = float(SER[(SER.fam == FAM) & (SER.netsta == ns)].tS_pick.median())  # picked S used by the stretch
    ax = axes[0, j]
    ax.plot(tau, grand / np.abs(grand).max(), "k", lw=0.8)
    ax.axvspan(tp + 2, tp + 4, color="tab:orange", alpha=0.2, label=f"stretched window: 2–4 s after picked S (S pick at {tp:+.2f} s)")
    ax.axvline(tp, color="tab:orange", lw=1.2)
    ax.axvline(0, color="r", lw=0.8, ls="--")
    ax.set_xlim(XL0, XL1); ax.set_title(f"family {FAM} at {ns} (Z, 2–8 Hz): grand stack of {len(g)} repeats")
    ax.legend(loc="upper left", fontsize=8)
    ax = axes[1, j]
    norm = S / (np.abs(S).max(1, keepdims=True) + 1e-12)
    ax.imshow(norm, aspect="auto", cmap="seismic", vmin=-1, vmax=1, interpolation="nearest",
              extent=[tau[0], tau[-1], len(days) - 0.5, -0.5])
    ax.axvline(0, color="k", lw=0.5, ls="--"); ax.axvline(tp, color="orange", lw=1.0); ax.axvline(tp + 2, color="k", lw=0.8); ax.axvline(tp + 4, color="k", lw=0.8)
    yt = [i for i in range(len(days)) if i == 0 or days[i].year != days[i - 1].year]
    ax.set_yticks(yt); ax.set_yticklabels([f"{days[i]:%Y-%m-%d}" for i in yt], fontsize=7)
    ax.set_xlim(XL0, XL1)
    ax.set_xlabel("time relative to predicted S (s); orange line = picked S")
    ax.set_title(f"{len(days)} daily stacks (one row per active day, ΣP ≥ 3); repeats/day median {int(np.median(nrep))}",
                 fontsize=9)
    ep = g.groupby("epoch").date.min()
    for e, d0 in ep.items():
        k = np.searchsorted(np.array(days, dtype="datetime64[ns]"), np.datetime64(d0))
        ax.axhline(k - 0.5, color="lime", lw=1.2); ax.text(-13.8, k + 1, f"epoch {e}", color="lime", fontsize=8)
fig.suptitle(f"Family {FAM}: daily P-weighted stacks, 2003–2026 (rows ordered by date; each row normalised)", y=0.995)
fig.tight_layout()
fig.savefig(OUT, dpi=130)
print(OUT)
