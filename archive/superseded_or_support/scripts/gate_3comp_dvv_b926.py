#!/usr/bin/env python
"""Gate B926 dv/v to physically-certified families per component and plot the 3-component result:
  Z  = 110 fwd-vs-rev certified families (ratio>1.5)
  H2 = 23 families passing the purity gate (caus>1.5 & fwd/rev>1.5)
  H1 = FAILS the physical-reality gate (anti-causal, 0 families) -> reported, not plotted
Certified dv/v = per-date median over the certified families (robust; ringing averages to ~0)."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

zc = set(pd.read_csv("data/b926_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)
h2c = set(pd.read_csv("data/b926_h2_pass_families.csv").fam)


def certified_median(csv, fams, cc_min=0.6):
    df = pd.read_csv(csv)
    df = df[df.patch.isin(fams) & (df.cc_max >= cc_min)]
    df["date"] = pd.to_datetime(df.date)
    med = df.groupby("date").dvv.median()
    n = df.groupby("date").patch.nunique()
    med = med[n >= 3]                      # >=3 families contributing that day
    return med.rolling(15, center=True, min_periods=5).median()   # light smooth for display


zmed = certified_median("data/daily_dvv_B926p90f40_Z_2to4.csv", zc)
h2med = certified_median("data/daily_dvv_B926p90f40_H2_2to4.csv", h2c)
print(f"Z  certified dv/v: {len(zmed)} days, std {zmed.std()*100:.2f}% ({len(zc)} families)")
print(f"H2 certified dv/v: {len(h2med)} days, std {h2med.std()*100:.2f}% ({len(h2c)} families)")
# Z-H2 agreement on overlapping dates (reported check, NOT a gate)
j = pd.concat([zmed.rename("Z"), h2med.rename("H2")], axis=1).dropna()
if len(j) > 30:
    print(f"Z vs H2 dv/v correlation (overlap {len(j)} days): r = {j.Z.corr(j.H2):.2f}")

fig, ax = plt.subplots(figsize=(13, 4.5))
ax.plot(zmed.index, zmed.values*100, lw=1.4, color="#1a1a2e", label=f"Z  (110 certified families)")
ax.plot(h2med.index, h2med.values*100, lw=1.2, color="#c1440e", alpha=0.85, label=f"H2 (23 certified families)")
ax.axhline(0, color="0.6", lw=0.6)
ax.set_ylabel("dv/v  (%)"); ax.set_title("B926 certified coda dv/v — Z + H2 (2–4 s, origin-anchored)\nH1 fails physical-reality gate (anti-causal coda) — not shown")
ax.legend(loc="upper right", frameon=False); ax.grid(alpha=0.2)
ax.set_ylim(-1, 1)
fig.tight_layout(); fig.savefig("figures/b926_certified_3comp_gated_dvv.png", dpi=130)
print("-> figures/b926_certified_3comp_gated_dvv.png")
