"""Grand stacks Z/N/E (and envelopes) for one family at given stations: is it P then S? (user question 2026-10-11)"""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.signal import hilbert
PIL = "/home/jovyan/tremorferometry/svi_science/results/s10_dvv/pilot"
FAM = int(sys.argv[1]) if len(sys.argv) > 1 else 10
STAS = ["PB.B011", "CN.PGC"]
rows = pd.read_parquet(f"{PIL}/tables/rows.parquet"); rows = rows[(rows.fam == FAM) & (rows.series == "data") & (rows.P >= 0.5)]
off = pd.read_csv(f"{PIL}/tables/offsets_sp.csv"); ser = pd.read_csv(f"{PIL}/tables/series.csv")
fig, axs = plt.subplots(2, 1, figsize=(13, 8), sharex=True)
for ax, ns in zip(axs, STAS):
    z = np.load(f"{PIL}/win/{ns}.npz"); pos = {r: i for i, r in enumerate(z["row_id"])}; tau = z["tau"]
    g = rows[rows.row_id.isin(pos)]
    W = z["W"][[pos[r] for r in g.row_id]].astype(np.float64); W[~np.isfinite(W)] = 0; w = g.P.to_numpy()
    sp = float(off[(off.fam == FAM) & (off.netsta == ns)].sp.iloc[0]); tp = float(ser[(ser.fam == FAM) & (ser.netsta == ns)].tS_pick.median())
    for c, (lab, col, dy) in enumerate([("Z", "k", 2.4), ("N/1", "tab:blue", 1.2), ("E/2", "tab:red", 0.0)]):
        s = (W[:, c] * w[:, None]).sum(0) / w.sum()
        a = np.abs(s).max() + 1e-12
        env = np.abs(hilbert(s))
        ax.plot(tau, s / a * 0.5 + dy, color=col, lw=0.8)
        ax.fill_between(tau, dy, env / a * 0.5 + dy, color=col, alpha=0.15)
        ax.text(-13.8, dy + 0.2, f"{lab}  (peak |amp| rel. to Z: {a / (np.abs((W[:, 0] * w[:, None]).sum(0) / w.sum()).max()):.2f})", fontsize=8, color=col)
    ax.axvline(0, color="grey", ls="--", lw=0.8); ax.text(0.05, 3.05, "predicted S", fontsize=8, color="grey")
    ax.axvline(-sp, color="green", ls="--", lw=1); ax.text(-sp + 0.05, 3.05, f"predicted P (S−P {sp:.2f} s)", fontsize=8, color="green")
    ax.axvline(tp, color="orange", lw=1); ax.axvspan(tp + 2, tp + 4, color="orange", alpha=0.15)
    ax.set_title(f"family {FAM} at {ns}: grand stacks per component (2–8 Hz), {len(g)} repeats; each trace normalised; shaded = envelope; orange = picked S and stretched window", fontsize=9)
    ax.set_yticks([]); ax.set_ylim(-0.7, 3.3)
axs[-1].set_xlabel("time relative to predicted S (s)"); axs[-1].set_xlim(-14, 10)
fig.tight_layout(); out = f"{PIL}/figures/s10_grand_3c_fam{FAM}.png"; fig.savefig(out, dpi=120); print(out)
