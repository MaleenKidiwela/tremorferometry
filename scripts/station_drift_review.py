#!/usr/bin/env python
"""Station review for the family-differential drift problem (Merlin split ruling 2026-09-02, step 2). NO rerun.
Evidence only — vetoes are the user's decision.

(a) Per-station yearly table for named stations (default B017, B023): station-mean raw eps, rows/yr, mean cc,
    per-family W2-W1 spread (raw and station-common-removed dw).
(b) Every era_table.csv boundary at a NON-whitelisted station: family-differential step test =
    per-family change of dw across the boundary (mean[b, b+365d) - mean[b-365d, b)), summarised as
    (i) support-vs-nonsupport (q16 pattern>0.3) difference of means, (ii) across-family SD of the change,
    each compared with the same statistics at CONTROL dates b-2y and b+2y (where data allow).
Writes data/station_drift_review_boundaries.csv; prints tables.
"""
import os, sys
import numpy as np, pandas as pd

RD = "fault_tomography/inversion/res_catalog"
STAS = sys.argv[1:] or ["B017p90f40", "B023p90f40"]
gz = np.load(f"{RD}/G.npz", allow_pickle=True)
keep = gz["captured"].astype(float) > 0
ptag, pcell = gz["pair_tag"][keep], gz["pair_cell"][keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}
zt = np.load(f"{RD}/transmission_T.npz", allow_pickle=True)
pat16 = np.asarray(zt["patterns"], float)[16]
cellpos = {c: i for i, c in enumerate(gz["cell"])}
insup_row = pat16[np.array([cellpos[c] for c in pcell])] > 0.3

fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "fam", "cell", "date", "dvv", "cc", "era", "src"])
fd["grow"] = [row_of.get((t, c), -1) for t, c in zip(fd.tag, fd.cell)]
fd = fd[fd.grow >= 0].copy(); fd["date"] = pd.to_datetime(fd.date); fd["yr"] = fd.date.dt.year
fd["dw"] = fd.dvv - fd.groupby(["tag", "date"]).dvv.transform("mean")
fd["insup"] = insup_row[fd.grow.values]

print("=== (a) station yearly tables ===")
for sta in STAS:
    x = fd[fd.tag == sta]
    if not len(x):
        print(sta, "no rows"); continue
    y = x.groupby("yr").agg(rows=("dvv", "size"), eps_mean=("dvv", "mean"), cc_mean=("cc", "mean"),
                            nfam=("fam", "nunique"))
    print(f"\n{sta} (src {x.src.iloc[0]}):"); print(y.round(3).to_string())
    w1 = x[(x.date >= "2011-01-01") & (x.date <= "2016-12-31")].groupby("fam").agg(n1=("dvv", "size"), r1=("dvv", "mean"), d1=("dw", "mean"))
    w2 = x[x.date >= "2022-01-01"].groupby("fam").agg(n2=("dvv", "size"), r2=("dvv", "mean"), d2=("dw", "mean"))
    b = w1.join(w2, how="inner"); b = b[(b.n1 >= 60) & (b.n2 >= 60)]
    if len(b):
        print(f"  bridge families {len(b)}: raw W2-W1 min/med/max {(b.r2-b.r1).min():+.2f}/{(b.r2-b.r1).median():+.2f}/"
              f"{(b.r2-b.r1).max():+.2f} %; station-common-removed D min/med/max {(b.d2-b.d1).min():+.2f}/"
              f"{(b.d2-b.d1).median():+.2f}/{(b.d2-b.d1).max():+.2f} %; |D|>0.3: {((b.d2-b.d1).abs()>0.3).sum()}")

print("\n=== (b) era boundaries at NON-whitelisted stations: family-differential step test ===")
et = pd.read_csv(f"{RD}/era_table.csv")
wl = {x.strip().upper() for x in open("data/era_split_whitelist.txt") if x.strip()}
rows = []
def step_stats(x, b):
    pre = x[(x.date >= b - pd.Timedelta(days=365)) & (x.date < b)]
    post = x[(x.date >= b) & (x.date < b + pd.Timedelta(days=365))]
    a = pre.groupby("fam").agg(n0=("dw", "size"), m0=("dw", "mean"), s=("insup", "first"))
    c = post.groupby("fam").agg(n1=("dw", "size"), m1=("dw", "mean"))
    j = a.join(c, how="inner"); j = j[(j.n0 >= 30) & (j.n1 >= 30)]
    if len(j) < 4:
        return dict(nfam=len(j), sd=np.nan, sup_minus_non=np.nan, nsup=int(j.s.sum()) if len(j) else 0)
    ch = j.m1 - j.m0
    smn = (ch[j.s].mean() - ch[~j.s].mean()) if (j.s.sum() >= 2 and (~j.s).sum() >= 2) else np.nan
    return dict(nfam=len(j), sd=float(ch.std()), sup_minus_non=float(smn) if smn == smn else np.nan, nsup=int(j.s.sum()))
for _, r in et.iterrows():
    if str(r.tag).upper().replace("P90F40", "") in wl or not isinstance(r.boundaries, str) or r.boundaries == "nan":
        continue
    x = fd[fd.tag == r.tag]
    if not len(x):
        continue
    for bs in r.boundaries.split(";"):
        if not bs: continue
        b = pd.Timestamp(bs)
        s0 = step_stats(x, b)
        ctl = [step_stats(x, b + pd.Timedelta(days=d)) for d in (-730, 730)]
        ctl_sd = np.nanmean([c["sd"] for c in ctl]) if any(np.isfinite(c["sd"]) for c in ctl) else np.nan
        ctl_smn = np.nanmean([c["sup_minus_non"] for c in ctl]) if any(np.isfinite(c["sup_minus_non"]) for c in ctl) else np.nan
        rows.append(dict(tag=r.tag, boundary=str(b.date()), nfam=s0["nfam"], nsup=s0["nsup"], fam_sd_step=s0["sd"],
                         ctl_sd=ctl_sd, sd_ratio=s0["sd"] / ctl_sd if ctl_sd and np.isfinite(ctl_sd) else np.nan,
                         sup_minus_non=s0["sup_minus_non"], ctl_sup_minus_non=ctl_smn))
tab = pd.DataFrame(rows)
tab.to_csv("data/station_drift_review_boundaries.csv", index=False)
print(tab.sort_values("sd_ratio", ascending=False).round(3).to_string(index=False))
print("\nwrote data/station_drift_review_boundaries.csv")
