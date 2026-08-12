#!/usr/bin/env python
"""Merlin's two gates before any resolution verdict.
GATE A (confound): on the 7 stations with a MIRROR product, is the 0.34% post-site family-specific residual
  mostly the removable mirror artifact (collapses toward sigma_meas ~0.13% after fwd - beta*mirror) or a real
  irreducible floor (stays >=0.3%)?  raw vs corrected, SAME families/months.
GATE B (statics): split each pair's post-site residual time series into a STATIC per-pair offset (cancels in
  epoch-difference maps) + a TIME-VARYING part; report variance shares + lag-1 month autocorr (-> annual
  stacking gain sqrt(N_indep)).  Uses the existing tensor (res_catalog/pair_months.parquet).
"""
import os, glob
import numpy as np, pandas as pd

CC = 0.70
OUT = "fault_tomography/inversion/res_catalog"


def certfams(stem):
    f = f"data/{stem}_causality_cert.csv"
    if os.path.exists(f):
        c = pd.read_csv(f); return set(c[c.reliable].fam)
    f = f"data/{stem}_fwd_vs_rev_coda.csv"
    if os.path.exists(f):
        c = pd.read_csv(f); return set(c[c.ratio > 1.5].fam)
    return set()


def deseason(dates, v):
    v = np.asarray(v, float); ok = ~np.isnan(v)
    if ok.sum() < 90: return v
    t = (pd.to_datetime(dates) - pd.to_datetime(dates[0])).days.values.astype(float); yr = 365.25
    X = np.column_stack([np.ones_like(t), t, np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr),
                         np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
    b, *_ = np.linalg.lstsq(X[ok], v[ok], rcond=None)
    return v - X @ b


print("=== GATE A: raw vs mirror-corrected post-site residual (7 mirror stations) ===")
rowsA = []
for mf in sorted(glob.glob("data/daily_dvv_*_MIRROR_2to4.csv")):
    tag = os.path.basename(mf).split("daily_dvv_")[1].split("_MIRROR")[0]
    raw_f = f"data/daily_dvv_{tag}_Z_2to4.csv"
    if not os.path.exists(raw_f):
        continue
    stem = tag.lower().replace("p90f40", "")
    fams = certfams(stem)
    R = pd.read_csv(raw_f); Mr = pd.read_csv(mf)
    R = R[(R.cc_max > CC) & (R.patch.isin(fams))]; Mr = Mr[(Mr.cc_max > CC) & (Mr.patch.isin(fams))]
    if not len(R) or not len(Mr):
        continue
    rp = R.pivot_table(index="date", columns="patch", values="dvv", aggfunc="first") * 100
    mp = Mr.pivot_table(index="date", columns="patch", values="dvv", aggfunc="first") * 100
    days = sorted(set(rp.index) & set(mp.index)); rp = rp.reindex(days); mp = mp.reindex(days)
    common = [c for c in rp.columns if c in mp.columns]
    rp = rp[common]; mp = mp[common]
    # station-level beta from deseasoned medians (FINAL_PIPELINE mirror-clean)
    cc_c = deseason(days, rp.median(axis=1).values); cc_m = deseason(days, mp.median(axis=1).values)
    jj = np.isfinite(cc_c) & np.isfinite(cc_m)
    beta = float(np.polyfit(cc_m[jj], cc_c[jj], 1)[0]) if jj.sum() > 30 else 0.0
    corr = rp - beta * mp
    # post-site residual per month = across-family std after removing the station-month mean
    def postsite_std(piv):
        idx = pd.to_datetime(piv.index).to_period("M")
        s = []
        for m, g in piv.groupby(idx):
            mo = g.mean(axis=0)                      # family monthly mean (percent)
            mo = mo.dropna()
            if len(mo) >= 5:
                s.append((mo - mo.mean()).std())     # remove station-month mean (the site term)
        return np.nanmedian(s) if s else np.nan
    rawstd, corstd = postsite_std(rp), postsite_std(corr)
    rowsA.append((tag, beta, len(common), rawstd, corstd))
    print(f"  {tag:14s} beta={beta:5.2f}  fams={len(common):3d}  post-site resid RAW {rawstd:.3f}% -> CORRECTED {corstd:.3f}%  "
          f"({100*(1-corstd/rawstd):+3.0f}%)")
A = pd.DataFrame(rowsA, columns=["tag", "beta", "nfam", "raw_resid", "corr_resid"])
print(f"  MEDIAN: raw {A.raw_resid.median():.3f}%  corrected {A.corr_resid.median():.3f}%  (sigma_meas ~0.13%)")
print(f"  -> {'OUTCOME A: artifact-dominated, mirror correction MANDATORY, map reopens' if A.corr_resid.median()<0.20 else 'OUTCOME B: real family-specific floor, index verdict stands'}")

print("\n=== GATE B: static vs time-varying decomposition of post-site residual (full tensor) ===")
pm = pd.read_parquet(f"{OUT}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["site_mean"] = pm.groupby(["tag", "ym"]).dvv_month.transform("mean")
pm["resid"] = pm.dvv_month - pm.site_mean
pm["pair"] = pm.tag + "|" + pm.cell
stat, tvar, ac1 = [], [], []
for p, g in pm.groupby("pair"):
    if len(g) < 6:
        continue
    g = g.sort_values("ym"); r = g.resid.values
    s = r.mean(); tv = r - s
    stat.append(abs(s)); tvar.append(tv.std())
    if len(tv) > 8 and tv.std() > 0:
        ac1.append(np.corrcoef(tv[:-1], tv[1:])[0, 1])
stat, tvar, ac1 = np.array(stat), np.array(tvar), np.array(ac1)
print(f"  pairs with >=6 months: {len(stat)}")
print(f"  |static offset|  median {np.median(stat):.3f}%   time-varying std median {np.median(tvar):.3f}%   (sigma_meas ~0.13%)")
vs = np.median(stat)**2; vt = np.median(tvar)**2
print(f"  variance share: static {100*vs/(vs+vt):.0f}%  time-varying {100*vt/(vs+vt):.0f}%")
print(f"  lag-1 monthly autocorr of time-varying part: median {np.nanmedian(ac1):.2f}  -> N_indep/yr ~ {12*(1-np.nanmedian(ac1))/(1+np.nanmedian(ac1)):.1f} (annual gain ~x{np.sqrt(12*(1-np.nanmedian(ac1))/(1+np.nanmedian(ac1))):.1f})")
print(f"  => transient/difference maps see sigma ~ time-varying std {np.median(tvar):.3f}% (statics cancel)")
