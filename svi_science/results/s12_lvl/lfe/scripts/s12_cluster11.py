#!/usr/bin/env python
"""s12_cluster11 -- PLAN §7 follow-up (coordinator, ~10:50 am PT): are the 11 WA-NW front families that carry the -0.33 %/yr
recharge slope a spatial cluster? Map vs the other 39; scrambled control = 1000 random 11-family sets from the 50 with
between-ETS data (compactness = RMS distance to the set centroid, km; slope = same §7 fit without bootstrap).
    S11_REGION=WA_NW python s12_cluster11.py
"""
import os, sys, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s12_recharge2 as R2
from s12_epoch import C, OUT, TB, log
rng = np.random.default_rng(1150)
F = pd.read_parquet(f"{TB}/b_days.parquet"); F["date"] = pd.to_datetime(F.date)
F["site"] = [C.SITE.get(n, n) for n in F.netsta]
ph = 2 * np.pi * F.date.dt.dayofyear.to_numpy() / 365.25
F["cs"], F["sn"], F["cs2"], F["sn2"] = np.cos(ph), np.sin(ph), np.cos(2 * ph), np.sin(2 * ph); F["ty"] = F.prev / 365.25
d = R2.prep(F[F.betw & (F.prev < 2000)], "dvv_2-5")
fam = pd.read_csv(f"{C.MTB}/families.csv").set_index("fam")
fams = np.array(sorted(d.fam.unique())); fams = fams[np.isin(fams, fam.index)]
f11 = np.array(sorted(d[d.has_front].fam.unique()))
lat0 = fam.loc[fams, "lat"].mean()
xy = lambda fs: np.c_[(fam.loc[fs, "lon"] - fam.loc[fams, "lon"].mean()) * 111.19 * np.cos(np.radians(lat0)), (fam.loc[fs, "lat"] - lat0) * 111.19]
def compact(fs):
    p = xy(fs); return float(np.sqrt(((p - p.mean(0)) ** 2).sum(1).mean()))
def slope(fs):
    return R2.fit(d[d.fam.isin(fs)], "dvv_2-5", boot=False)["b"]
c0, s0 = compact(f11), slope(f11)
cs, ss = [], []
for k in range(1000):
    r = rng.choice(fams, len(f11), replace=False); cs.append(compact(r)); ss.append(slope(r))
cs, ss = np.array(cs), np.array(ss)
PGC = (48.6495, -123.4513)
def dist(fs, la, lo):
    return np.hypot((fam.loc[fs, "lon"] - lo) * 111.19 * np.cos(np.radians(la)), (fam.loc[fs, "lat"] - la) * 111.19)
oth = np.setdiff1d(fams, f11)
out = dict(n11=len(f11), n_other=len(oth), compact11_km=c0, compact_scr_median=float(np.median(cs)),
           p_compact=float((np.sum(cs <= c0) + 1) / 1001), slope11=s0, slope_scr_sd=float(ss.std()),
           p_slope=float((np.sum(ss <= s0) + 1) / 1001),
           p_joint=float((np.sum((cs <= c0) & (ss <= s0)) + 1) / 1001),
           lat11_med=float(fam.loc[f11, "lat"].median()), lon11_med=float(fam.loc[f11, "lon"].median()),
           depth11_med=float(fam.loc[f11, "depth_km"].median()), depth11_range=[float(fam.loc[f11, "depth_km"].min()), float(fam.loc[f11, "depth_km"].max())],
           lat_other_med=float(fam.loc[oth, "lat"].median()), lon_other_med=float(fam.loc[oth, "lon"].median()),
           depth_other_med=float(fam.loc[oth, "depth_km"].median()), depth_other_q90=float(fam.loc[oth, "depth_km"].quantile(.9)),
           dPGC11_med_km=float(dist(f11, *PGC).median()), dPGC_other_med_km=float(dist(oth, *PGC).median()),
           n11_lat_lt_48_4=int((fam.loc[f11, "lat"] < 48.4).sum()), n_other_lat_lt_48_4=int((fam.loc[oth, "lat"] < 48.4).sum()),
           families11=[int(x) for x in f11])
json.dump(out, open(f"{TB}/r2_cluster11.json", "w"), indent=1)
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].scatter(fam.loc[oth, "lon"], fam.loc[oth, "lat"], c=fam.loc[oth, "depth_km"], cmap="viridis", vmin=25, vmax=45, marker="o", s=30, label="other 39")
sc = ax[0].scatter(fam.loc[f11, "lon"], fam.loc[f11, "lat"], c=fam.loc[f11, "depth_km"], cmap="viridis", vmin=25, vmax=45, marker="*", s=180, edgecolor="r", label="11 with fronts (−0.33 %/yr)")
ax[0].plot(PGC[1], PGC[0], "k^", ms=9, label="PGC"); plt.colorbar(sc, ax=ax[0], label="depth (km)")
ax[0].set_xlabel("lon"); ax[0].set_ylabel("lat"); ax[0].legend(fontsize=7); ax[0].set_title("WA-NW families with between-ETS data")
ax[1].scatter(cs, ss, s=6, c="0.6", label="random 11-family sets"); ax[1].scatter([c0], [s0], c="r", s=80, marker="*", label="the 11")
ax[1].set_xlabel("RMS distance to centroid (km)"); ax[1].set_ylabel("slope (%/yr)"); ax[1].legend(fontsize=7)
ax[1].set_title(f"compactness p {out['p_compact']:.3f}, slope p {out['p_slope']:.3f}, joint p {out['p_joint']:.3f}")
fig.tight_layout(); fig.savefig(f"{OUT}/figures/fig_cluster11.png", dpi=120)
log(f"[cluster11] {out}")
