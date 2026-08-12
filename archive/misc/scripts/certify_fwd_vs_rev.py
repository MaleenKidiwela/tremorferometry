#!/usr/bin/env python
"""Reusable per-station Z certification: fwd-vs-reversed grand-stack coda ratio.
Real LFE coda >> matched-filter ringing (which the time-flipped template reproduces symmetrically).
Per family: ratio = RMS(fwd Z grand coda 2-4s) / RMS(rev Z grand coda 2-4s). Keep ratio>1.5.
Grand stack = n_det-weighted mean of daily stacks. Usage:
  python certify_fwd_vs_rev.py <STA_TAG>  e.g. B926p90f40   (reads long_window_daily_<tag>_Z.npz + <tag>rev_Z.npz)
Writes data/<sta_lower>_fwd_vs_rev_coda.csv (cols fam,fwd,rev,ratio)."""
import sys, re, numpy as np, pandas as pd

tag = sys.argv[1]                      # e.g. B926p90f40
sta = re.match(r"[A-Za-z]+\d+", tag).group(0).lower()
CODA = (2.0, 4.0)


def grands(path):
    d = np.load(path, allow_pickle=True)
    t, p, nd, S = d["t"], d["patches"], d["n_det"].astype(float), d["stacks"]
    coda = (t >= CODA[0]) & (t <= CODA[1])
    out = {}
    for fam in pd.unique(p):
        m = p == fam; w = nd[m]
        if w.sum() <= 0: continue
        g = (S[m] * w[:, None]).sum(0) / w.sum()
        out[fam] = float(np.sqrt(np.mean(g[coda] ** 2)))
    return out


fwd = grands(f"data/long_window_daily_{tag}_Z.npz")
rev = grands(f"data/long_window_daily_{tag}rev_Z.npz")
rows = []
for fam, fc in fwd.items():
    rc = rev.get(fam)
    if rc is None or rc <= 0: continue
    rows.append(dict(fam=fam, fwd=fc, rev=rc, ratio=fc / rc))
R = pd.DataFrame(rows).sort_values("ratio", ascending=False)
out = f"data/{sta}_fwd_vs_rev_coda.csv"
R.to_csv(out, index=False)
n_cert = int((R.ratio > 1.5).sum())
print(f"[{sta}] {len(R)} families with fwd+rev stacks | CERTIFIED (ratio>1.5): {n_cert} "
      f"| ratio median {R.ratio.median():.2f} -> {out}")
