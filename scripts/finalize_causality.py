#!/usr/bin/env python
"""Z-only causality certification + gated dv/v (the retired-reverse fleet method).
Certify reliable families from the FORWARD stack alone: causality = RMS(coda 2-4s)/RMS(mirror -2..0s) > 1.5
(verified equivalent to fwd-vs-reversed at 99%, corr 1.00). No reversed stacks needed.
Usage: finalize_causality.py <STA> <TAG>   needs long_window_daily_<TAG>_Z.npz + daily_dvv_<TAG>_Z_2to4.csv
Writes data/<sta>_causality_cert.csv, data/<sta>_3comp_summary.json, figures/<sta>_certified_dvv.png"""
import sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

STA, TAG = sys.argv[1], sys.argv[2]; sta = STA.lower()
d = np.load(f"data/long_window_daily_{TAG}_Z.npz", allow_pickle=True)
t, p, nd, S = d["t"], d["patches"], d["n_det"].astype(float), d["stacks"]
coda = (t >= 2) & (t <= 4); mir = (t >= -2) & (t < 0)
rms = lambda x: float(np.sqrt(np.mean(x**2)))

rows = []
for fam in pd.unique(p):
    m = p == fam; w = nd[m]
    if w.sum() <= 0: continue
    g = (S[m] * w[:, None]).sum(0) / w.sum()
    fm = rms(g[mir])
    rows.append((fam, rms(g[coda]) / fm if fm else np.nan))
C = pd.DataFrame(rows, columns=["fam", "causality"]).dropna()
C["reliable"] = C.causality > 1.5
C.to_csv(f"data/{sta}_causality_cert.csv", index=False)
cert = set(C[C.reliable].fam)
print(f"[{STA}] causality certification: {len(cert)}/{len(C)} reliable (caus>1.5)")

df = pd.read_csv(f"data/daily_dvv_{TAG}_Z_2to4.csv")
df = df[df.patch.isin(cert) & (df.cc_max >= 0.6)]; df["date"] = pd.to_datetime(df.date)
med = df.groupby("date").dvv.median(); n = df.groupby("date").patch.nunique()
zmed = med[n >= 3].rolling(15, center=True, min_periods=5).median()

fig, ax = plt.subplots(figsize=(13, 4.3))
ax.plot(zmed.index, zmed.values*100, lw=1.4, color="#1a1a2e", label=f"Z ({len(cert)} causality-certified families)")
ax.axhline(0, color="0.6", lw=0.6); ax.set_ylim(-1, 1); ax.grid(alpha=0.2)
ax.set_ylabel("dv/v (%)"); ax.legend(loc="upper right", frameon=False)
ax.set_title(f"{STA} certified coda dv/v (2-4 s, origin-anchored) — Z-only, causality cert")
fig.tight_layout(); fig.savefig(f"figures/{sta}_certified_dvv.png", dpi=130)
summ = dict(station=STA, tag=TAG, method="causality", z_certified=len(cert),
            z_std_pct=round(float(zmed.std()*100), 3))
json.dump(summ, open(f"data/{sta}_3comp_summary.json", "w"), indent=2)
print(f"[{STA}] Z dv/v std {summ['z_std_pct']}% -> figures/{sta}_certified_dvv.png ; SUMMARY {json.dumps(summ)}")
