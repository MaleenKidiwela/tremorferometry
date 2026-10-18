#!/usr/bin/env python
"""P1a baseline figure: 2011 / 2013 own recall against the 2012 S15 ablation reference."""
import json, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
RES = "/home/jovyan/tremorferometry/lfe_picker/results"
FIG = "/home/jovyan/tremorferometry/lfe_picker/figures"
REF = {"T_core": 0.1463, "T_ext": 0.1841}
C = {"T_core": "#2a78d6", "T_ext": "#6da7ec"}
R = pd.concat([pd.read_csv(f"{RES}/o2_recall_{y}.csv") for y in (2011, 2013)])
R = R[R.scope == "data_window"]
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
yrs = [2011, 2013]; x = np.arange(len(yrs)); w = 0.36
for i, t in enumerate(("T_core", "T_ext")):
    v = [float(R[(R.year == y) & (R.tier == t)].recall.iloc[0]) for y in yrs]
    lo = [float(R[(R.year == y) & (R.tier == t)].recall_lo.iloc[0]) for y in yrs]
    hi = [float(R[(R.year == y) & (R.tier == t)].recall_hi.iloc[0]) for y in yrs]
    ax[0].bar(x + (i - 0.5) * w, v, w, color=C[t], label=t,
              yerr=[np.array(v) - lo, np.array(hi) - np.array(v)], capsize=3, ecolor="#52514e")
    ax[0].axhline(REF[t], color=C[t], ls="--", lw=1)
    ax[0].axhspan(0.7 * REF[t], 1.3 * REF[t], color=C[t], alpha=0.10)
ax[0].set_xticks(x); ax[0].set_xticklabels([str(y) for y in yrs])
ax[0].set_ylabel("own recall of Bostock (+-2 s, one-to-one)")
ax[0].set_title("P1a: frozen template detector\ndashed = 2012 S15 reference, band = 0.7-1.3x accept")
ax[0].legend(fontsize=8, frameon=False)
F = pd.concat([pd.read_csv(f"{RES}/o2_recall_family_{y}.csv") for y in (2011, 2013)])
for y, c in ((2011, "#2a78d6"), (2013, "#eb6834")):
    g = F[(F.year == y) & (F.tier == "T_core") & (F.n_events >= 20)].recall.to_numpy()
    ax[1].hist(g, bins=np.arange(0, 0.75, 0.04), histtype="step", lw=1.8, color=c,
               label=f"{y}  (n={len(g)} families, median {np.median(g):.3f})")
ax[1].set_xlabel("per-family own recall at T_core"); ax[1].set_ylabel("families")
ax[1].set_title("per-family recall, families with >= 20 events")
ax[1].legend(fontsize=8, frameon=False)
fig.tight_layout(); fig.savefig(f"{FIG}/o2_p1a_baseline.png", dpi=150)
print("ok")
