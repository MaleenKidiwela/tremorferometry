#!/usr/bin/env python
"""D3 stage 5 figure: what the network sum and the template length actually buy."""
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

R = Path("lfe_catalog_repair/results")
S = pd.read_csv(R / "d3_stage5_summary.csv")
P = pd.read_csv(R / "pgconly" / "d3_stage5_summary.csv") if (R / "pgconly" / "d3_stage5_summary.csv").exists() else None
ST = pd.read_csv(R / "d3_stage5_stacks.csv")
fig, ax = plt.subplots(1, 4, figsize=(19, 4.4))

# A sigma vs template length, measured vs predicted
Ls = sorted(S.L.unique())
sig = [float(S[S.L == L].sigma_net.iloc[0]) for L in Ls]
pred1 = [1 / np.sqrt(2 * 7 * L) for L in Ls]
ax[0].plot(Ls, pred1, "k--o", label="1 channel, theory 1/sqrt(2BL)")
ax[0].plot(Ls, [p / np.sqrt(21) for p in pred1], "--s", color="#1b7f4b",
           label="21 channels, ideal sqrt(N)")
ax[0].plot(Ls, sig, "-o", color="#c0392b", lw=2, label="21 channels, MEASURED")
if P is not None:
    ax[0].plot([6], [float(P[P.L == 6].sigma_net.iloc[0])], "*", ms=14, color="#1b6ca8",
               label="PGC only (3 ch), measured")
ax[0].axhline(0.078 / (0.6745 * 9), color="0.4", ls=":",
              label="sigma needed for k=9 on a median LFE")
ax[0].axhline(0.078 / (0.6745 * 6), color="0.7", ls=":", label="... for k=6")
ax[0].set_xlabel("template length L (s)")
ax[0].set_ylabel(r"network cc noise $\sigma$")
ax[0].set_yscale("log")
ax[0].legend(fontsize=7)
ax[0].grid(alpha=.3)
ax[0].set_title("A  the network sum and template length both work,\n"
                "and both fall short of k=9", fontsize=10, loc="left")

# B detections vs k
for L, c in zip(Ls, ["#4a90c2", "#1b7f4b", "#c0392b"]):
    g = S[S.L == L].sort_values("k")
    ax[1].semilogy(g.k, np.maximum(g.n_det, 0.5), "-o", color=c, label=f"L = {L:.0f} s")
ax[1].set_xlabel("threshold k (median + k*MAD of the network cc)")
ax[1].set_ylabel("network detections in 34 days")
ax[1].legend(fontsize=8)
ax[1].grid(alpha=.3)
ax[1].set_title("B  the detector goes silent between k=6 and k=8", fontsize=10, loc="left")

# C precision vs the shifted control
g = S[S.L == 6].sort_values("k")
ax[2].plot(g.k, g.precision, "-o", color="#c0392b", label="precision")
ax[2].plot(g.k, g.precision_shift, "-^", color="0.5", label="+600 s shifted control")
ax[2].set_xlabel("threshold k")
ax[2].set_ylabel("precision vs held-out refined Lin picks")
ax[2].legend(fontsize=8)
ax[2].grid(alpha=.3)
ax[2].set_title("C  L = 6 s: precision never leaves the control", fontsize=10, loc="left")

# D per-station stack SNR vs null
ST2 = ST.sort_values("snr_h", ascending=False)
xx = np.arange(len(ST2))
ax[3].bar(xx - 0.2, ST2.snr_h, 0.4, color="#1b7f4b", label="real stack (388 picks)")
ax[3].bar(xx + 0.2, ST2.null_snr_h, 0.4, color="0.6", label="random-time null")
ax[3].set_xticks(xx)
ax[3].set_xticklabels(ST2.station, rotation=45, fontsize=8)
ax[3].set_ylabel("stack SNR on |H|")
ax[3].legend(fontsize=8)
ax[3].grid(alpha=.3, axis="y")
ax[3].set_title("D  only PGC and OCP clear their own null", fontsize=10, loc="left")

fig.suptitle("D3 stage 5 — network matched filter on the Lin/Bostock ground truth "
             "(2010-08-10..09-13, 9 stations / 21 channels)", fontsize=12)
fig.tight_layout(rect=[0, 0, 1, 0.94])
out = "lfe_catalog_repair/figures/d3_stage5_network.png"
fig.savefig(out, dpi=125, bbox_inches="tight")
print("->", out)
