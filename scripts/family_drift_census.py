#!/usr/bin/env python
"""Family-differential slow-drift census (Merlin split-network ruling 2026-09-02, prescription step 1). NO rerun.

For every family-era in res_catalog/family_days.parquet that maps to a G row:
  dw  = dvv - station daily family-mean  (the data-level analogue of the per-station-day site term s_j(t))
  D   = mean(dw, W2=2022-26) - mean(dw, W1=2011-16)   where >= 60 rows in BOTH windows  ("bridge" family-eras)
  R365 = range (max-min) of the 365-d rolling mean of dw over the family-era span (monthly means, 12-month window
         with >= 60 rows in the window), for ALL family-eras
Flag: |D| > 0.3 % or R365 > 0.4 %  (Merlin's bars; a -2 % deep field predicts |d eps| ~0.08 % at median capture).

Writes NEW files only: data/family_drift_census.csv (one row per family-era) and prints the station table.
Gate values can be changed downstream (0.2 / 0.3 / 0.5 %) from the same CSV.
"""
import os, sys
import numpy as np, pandas as pd

RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog")
OUTCSV = os.environ.get("DRIFT_CSV", "data/family_drift_census.csv")
D_BAR, R_BAR = 0.3, 0.4
W1 = (pd.Timestamp("2011-01-01"), pd.Timestamp("2016-12-31"))
W2 = (pd.Timestamp("2022-01-01"), pd.Timestamp("2026-12-31"))

gz = np.load(f"{RD}/G.npz", allow_pickle=True)
keep = gz["captured"].astype(float) > 0
ptag, pcell = gz["pair_tag"][keep], gz["pair_cell"][keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}
zt = np.load(f"{RD}/transmission_T.npz", allow_pickle=True)
pat16 = np.asarray(zt["patterns"], float)[16]
cellpos = {c: i for i, c in enumerate(gz["cell"])}
pc = np.array([cellpos[c] for c in pcell])
insup_row = pat16[pc] > 0.3

fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "fam", "cell", "date", "dvv", "era", "src"])
fd["grow"] = [row_of.get((t, c), -1) for t, c in zip(fd.tag, fd.cell)]
fd = fd[fd.grow >= 0].copy()
fd["date"] = pd.to_datetime(fd.date)
fd["dw"] = fd.dvv - fd.groupby(["tag", "date"]).dvv.transform("mean")
fd["fe"] = fd.tag + "|" + fd.fam + "|e" + fd.era.astype(str)
print(f"rows {len(fd):,} | family-eras {fd.fe.nunique()} | stations {fd.tag.nunique()}", flush=True)

# --- two-window level change D (bridge family-eras) ---
m1 = (fd.date >= W1[0]) & (fd.date <= W1[1]); m2 = (fd.date >= W2[0]) & (fd.date <= W2[1])
g1 = fd[m1].groupby("fe").agg(n1=("dw", "size"), mean1=("dw", "mean"))
g2 = fd[m2].groupby("fe").agg(n2=("dw", "size"), mean2=("dw", "mean"))
meta = fd.groupby("fe").agg(tag=("tag", "first"), fam=("fam", "first"), era=("era", "first"), src=("src", "first"),
                            grow=("grow", "first"), n=("dw", "size"), t0=("date", "min"), t1=("date", "max"))
meta["insup16"] = insup_row[meta.grow.values]
tab = meta.join(g1, how="left").join(g2, how="left")
bridge = (tab.n1 >= 60) & (tab.n2 >= 60)
tab["D"] = np.where(bridge, tab.mean2 - tab.mean1, np.nan)

# --- 365-d rolling-mean range (monthly means, 12-month window, >= 60 rows in window) ---
fd["ym"] = fd.date.dt.to_period("M")
mo = fd.groupby(["fe", "ym"]).dw.agg(["sum", "size"]).reset_index()
mo["ym"] = mo.ym.dt.to_timestamp()
R365 = {}
for fe, g in mo.groupby("fe", sort=False):
    g = g.set_index("ym").sort_index()
    full = pd.date_range(g.index.min(), g.index.max(), freq="MS")
    s = g["sum"].reindex(full, fill_value=0.0); c = g["size"].reindex(full, fill_value=0)
    rs = s.rolling(12, min_periods=1).sum(); rc = c.rolling(12, min_periods=1).sum()
    rm = (rs / rc.replace(0, np.nan))[rc >= 60]
    R365[fe] = (float(rm.max() - rm.min()) if len(rm) >= 2 else np.nan, int(len(rm)))
tab["R365"] = tab.index.map(lambda k: R365.get(k, (np.nan, 0))[0])
tab["n_R365"] = tab.index.map(lambda k: R365.get(k, (np.nan, 0))[1])
tab["flag_D"] = tab.D.abs() > D_BAR
tab["flag_R"] = tab.R365 > R_BAR
tab["flagged"] = tab.flag_D | tab.flag_R
tab.index.name = "fe"
tab.round(4).to_csv(OUTCSV)

b = tab[bridge]
print(f"\nBRIDGE family-eras (>=60 rows in both W1 and W2): {len(b)} at {b.tag.nunique()} stations")
print(f"  median |D| {b.D.abs().median():.3f}%  (boreholes {b[b.tag.str.contains('p90f40')].D.abs().median():.3f}, "
      f"broadband {b[~b.tag.str.contains('p90f40')].D.abs().median():.3f})  p75 {b.D.abs().quantile(.75):.3f}  "
      f"p90 {b.D.abs().quantile(.9):.3f}  | |D|>0.3%: {(b.D.abs() > 0.3).sum()}/{len(b)}")
print(f"  by src: " + ", ".join(f"{s} med|D| {v:.3f}" for s, v in b.groupby('src').D.apply(lambda v: v.abs().median()).items()))
print(f"\nALL family-eras: {len(tab)} | R365 available {tab.R365.notna().sum()} | median R365 {tab.R365.median():.3f}% "
      f"| R365>0.4%: {(tab.R365 > 0.4).sum()}")
print(f"FLAGGED (|D|>{D_BAR} or R365>{R_BAR}): {tab.flagged.sum()}/{len(tab)} = {tab.flagged.mean():.2f} "
      f"| among bridge: {b.flagged.sum()}/{len(b)} = {b.flagged.mean():.2f} "
      f"| among q16-support: {tab[tab.insup16].flagged.sum()}/{tab.insup16.sum()} = {tab[tab.insup16].flagged.mean():.2f}")
for gate in (0.2, 0.3, 0.5):
    fl = (tab.D.abs() > gate) | (tab.R365 > gate * 4 / 3)
    print(f"  gate |D|>{gate} or R365>{gate*4/3:.2f}: flagged {fl.sum()} ({fl.mean():.2f}); rows dropped "
          f"{tab.n[fl].sum() / tab.n.sum():.2f}")
st = tab.groupby("tag").agg(src=("src", "first"), n_fe=("D", "size"), n_bridge=("D", lambda v: v.notna().sum()),
                            medD=("D", lambda v: np.nanmedian(np.abs(v)) if v.notna().any() else np.nan),
                            medR=("R365", "median"), frac_flag=("flagged", "mean"),
                            sup16=("insup16", "sum")).sort_values("medD", ascending=False)
print("\nSTATION TABLE (sorted by median |D| among bridge family-eras; top 40):")
print(st.head(40).round(3).to_string())
print(f"\nwrote {OUTCSV}")
