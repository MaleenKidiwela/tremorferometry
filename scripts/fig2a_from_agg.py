#!/usr/bin/env python
"""Fig-2A-style full-record figure from data/mf_certified_agg.npz (certified-only, vetoes applied):
time x along-strike density of ALL matched-filter repeats at cc>=0.8, + daily-rate strip."""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

z = np.load("data/mf_certified_agg.npz", allow_pickle=True)
H = np.asarray(z["H_cc80"], float)              # (days, sbins)
sb = z["sbins"]; T0 = pd.Timestamp(str(z["t0"]))
days = pd.date_range(T0, periods=H.shape[0], freq="D")
tot = H.sum()
# trim empty along-strike edges
occ = H.sum(0) > 0
i0, i1 = np.argmax(occ), len(occ) - np.argmax(occ[::-1])
fig = plt.figure(figsize=(20, 8.5))
axh = fig.add_axes([0.05, 0.88, 0.92, 0.08])
axa = fig.add_axes([0.05, 0.07, 0.92, 0.78])
tb = mdates.date2num(days)
tedges = np.concatenate([tb - 0.5, [tb[-1] + 0.5]])
pc = axa.pcolormesh(tedges, sb[i0:i1+1], np.log10(H[:, i0:i1].T + 1), cmap="inferno", shading="flat")
cb = fig.colorbar(pc, ax=axa, pad=0.005); cb.set_label("log10 (repeats / day / 10 km)")
axa.set_ylabel("along-strike distance (km, + = north)"); axa.set_xlabel("time"); axa.xaxis_date()
rate = H.sum(1)
axh.fill_between(tb, rate, step="mid", color="#555")
axh.set_xlim(axa.get_xlim()); axh.set_xticks([]); axh.set_ylabel("repeats/day", fontsize=8)
axh.set_yscale("log")
axh.set_title(f"FULL repeating-LFE record (matched-filter, cc≥0.8): {tot:,.0f} detections, certified "
              f"families only, vetoes applied — daily × 10-km bins", fontsize=11)
fig.savefig("figures/lfe_fronts_full_mf.png", dpi=140, bbox_inches="tight")
print(f"wrote figures/lfe_fronts_full_mf.png | total certified cc>=0.8 detections {tot:,.0f}")
