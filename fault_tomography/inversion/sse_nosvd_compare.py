#!/usr/bin/env python
"""Mechanism test: does the station-specific ETS-locked apparent dv/v shift survive WITHOUT the SVD-Wiener filter?
Per station with a data/daily_dvv_<TAG>_Z_2to4_nosvd.csv: per-family episode contrast (ON[+15,+60]-PRE[-90,-30],
monthly-deseasoned, admitted episodes within 150 km) with the production (SVD-Wiener) file vs the no-SVD file, same
families/rows (cc>0.7 in both). Shrinkage-release predicts the contrast SHRINKS toward zero (or vanishes) without SVD;
a stack-composition/noise-field effect predicts it PERSISTS. Also the whole-record corr with local slip-in-window.
"""
import numpy as np, pandas as pd, os, glob, sys
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; gi = pd.Series(np.arange(G), index=days)
clat, clon, dep = zc["cell_lat"], zc["cell_lon"], zc["iface_depth_km"]; inc = zc["slip_thrust_inc30"].astype(float) * 1e3
km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); XY = km(clat, clon)
ep = pd.read_csv(f"{SD}/episodes.csv"); ep["key"] = ep.seg + "_" + ep.peak; okep = ep[ep.ok & ep.seg.isin(["N", "C", "S"])]; zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True)
pairs = pd.read_csv(f"{RD}/pairs.csv"); sll = pairs.groupby("tag")[["sta_lat", "sta_lon"]].first()
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "fam"]).drop_duplicates(); famset = {t: set(g.fam) for t, g in fd.groupby("tag")}
win = np.asarray((days >= "2011-01-01") & (days <= "2026-07-20"))
def contr(x, ions):
    piv = x.pivot_table(index="date", columns="fam", values="v").reindex(days); piv = piv - piv.groupby(piv.index.month).transform("mean"); out = {}
    for fam in piv.columns:
        y = piv[fam].values.astype(float); vals = []
        for i in ions:
            a = y[max(i + 15, 0):i + 61]; b = y[max(i - 90, 0):i - 29]
            if np.isfinite(a).mean() > .6 and np.isfinite(b).mean() > .6: vals.append(np.nanmean(a) - np.nanmean(b))
        if len(vals) >= 4: out[fam] = np.mean(vals)
    return pd.Series(out, dtype=float), piv
rows = []
for f in sorted(glob.glob("data/daily_dvv_*_Z_2to4_nosvd.csv")):
    tag = os.path.basename(f).split("daily_dvv_")[1].split("_Z_2to4")[0]
    if tag not in sll.index: continue
    sxy = km(np.array([sll.sta_lat[tag]]), np.array([sll.sta_lon[tag]]))[0]
    ions = [int(gi[pd.Timestamp(r.fp_onset)]) for _, r in okep.iterrows() if np.sqrt(((XY[zf[f"fp_{r.key}"]] - sxy) ** 2).sum(1)).min() <= 150]
    near = np.isfinite(dep) & (dep >= 30) & (dep <= 45) & (np.sqrt(((XY - sxy) ** 2).sum(1)) <= 60); local = np.nanmean(inc[:, near], 1) if near.sum() >= 3 else None
    prod = pd.read_csv(f"data/daily_dvv_{tag}_Z_2to4.csv", usecols=["patch", "date", "dvv", "cc_max"]); nos = pd.read_csv(f, usecols=["patch", "date", "dvv", "cc_max"])
    for x in (prod, nos): x.rename(columns={"patch": "fam"}, inplace=True); x["date"] = pd.to_datetime(x.date); x["v"] = -x.dvv * 100
    keep = prod.merge(nos, on=["fam", "date"], suffixes=("_p", "_n")); keep = keep[(keep.cc_max_p > 0.7) & (keep.cc_max_n > 0.7) & keep.fam.isin(famset[tag])]
    cp, pp = contr(keep.rename(columns={"v_p": "v"})[["fam", "date", "v"]], ions); cn, pn = contr(keep.rename(columns={"v_n": "v"})[["fam", "date", "v"]], ions)
    j = pd.concat([cp.rename("svd"), cn.rename("nosvd")], axis=1).dropna()
    yp = pp.mean(axis=1).values; yn = pn.mean(axis=1).values
    def d30(y): o = np.full_like(y, np.nan); o[30:] = y[30:] - y[:-30]; return o
    def corr(a, b): m = np.isfinite(a) & np.isfinite(b) & win; return np.corrcoef(a[m], b[m])[0, 1] if m.sum() > 500 else np.nan
    rows.append(dict(tag=tag, n_fam=len(j), n_ep=len(ions), contrast_svd=j.svd.mean(), contrast_nosvd=j.nosvd.mean(), fpos_svd=(j.svd > 0).mean(), fpos_nosvd=(j.nosvd > 0).mean(),
                     corr_fam=j.corr().iloc[0, 1], cc_med_svd=float(keep.cc_max_p.median()), cc_med_nosvd=float(keep.cc_max_n.median()),
                     r_local_svd=corr(d30(yp), local) if local is not None else np.nan, r_local_nosvd=corr(d30(yn), local) if local is not None else np.nan,
                     sd_daily_svd=float(np.nanstd(np.diff(yp))), sd_daily_nosvd=float(np.nanstd(np.diff(yn)))))
T = pd.DataFrame(rows); T.to_csv(f"{SD}/nosvd_compare.csv", index=False)
print(T.round(3).to_string(index=False))
print(f"\nratio of ETS contrasts nosvd/svd (mean over stations of |c_nosvd|/|c_svd|): {np.mean(T.contrast_nosvd.abs() / T.contrast_svd.abs()):.2f}; sign agreement {100*np.mean(np.sign(T.contrast_nosvd) == np.sign(T.contrast_svd)):.0f}%")
