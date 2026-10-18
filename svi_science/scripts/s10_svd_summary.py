#!/usr/bin/env python
"""s10_svd_summary -- raw vs SVD tables and pass verdicts for the S10 SVD-Wiener step (rules: svd/SVD.md §3-5)."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S

TAG = sys.argv[1] if len(sys.argv) > 1 else ""
T = f"{S.OUT}/svd/tables"
SER = pd.read_csv(f"{T}/series_svd{TAG}.csv")
INJ = pd.read_csv(f"{T}/injection_svd{TAG}.csv")
CT = pd.read_csv(f"{T}/control_svd{TAG}.csv")
D = pd.read_parquet(f"{T}/daily_svd{TAG}.parquet")
KEY = ["fam", "netsta", "epoch"]

# ---- injection verdicts (§4)
st = INJ[INJ.kind == "step"]
rows = []
for (f, n, e, p), g in st.groupby(KEY + ["path"]):
    r = dict(fam=f, netsta=n, epoch=e, path=p, n_dates=g.boundary.nunique(),
             rec_min=g.recovery.min(), rec_max=g.recovery.max(), rec_med=g.recovery.median(),
             leak_fixed_maxabs=g.leak_fixed.abs().max(), leak_contrast_maxabs=g.leak_contrast.abs().max(),
             rec_fixed_med=g.recovery_fixed.median())
    lc = g.leak_contrast.dropna().abs()
    r["inj_pass"] = bool(g.recovery.between(0.8, 1.2).all() and (g.leak_fixed.abs() < 0.2).all() and (lc < 0.2).all())
    for kind in ("episode", "day"):
        h = INJ[(INJ.kind == kind) & (INJ.fam == f) & (INJ.netsta == n) & (INJ.epoch == e) & (INJ.path == p)]
        r[f"transient_{kind}_rec24"] = float(h[h.win == 24].recovery.mean()) if len(h) else np.nan
        r[f"transient_{kind}_rec35"] = float(h[h.win == 35].recovery.mean()) if len(h) else np.nan
    rows.append(r)
IV = pd.DataFrame(rows)
IV.to_csv(f"{T}/injection_summary{TAG}.csv", index=False)

# ---- scrambled control (§5)
crow = []
for (f, n, e), g in CT.groupby(KEY):
    for v in ("raw", "own", "ins"):
        if f"{v}_eps24" not in g:
            continue
        r = dict(fam=f, netsta=n, epoch=e, variant=v, n_ctrl_days=len(g))
        for k in (24, 35):
            z = (g[f"{v}_eps{k}"] / g[f"{v}_sig{k}"]).to_numpy(); z = z[np.isfinite(z)]
            r[f"n{k}"] = len(z); r[f"frac_absz_lt2_{k}"] = float(np.mean(np.abs(z) < 2)) if len(z) else np.nan
            r[f"mean_z{k}"] = float(np.mean(z)) if len(z) else np.nan
            r[f"cc{k}_med"] = float(g[f"{v}_cc{k}"].median()) if g[f"{v}_cc{k}"].notna().any() else np.nan
            r[f"rail{k}"] = float(g[f"{v}_eps{k}"].isna().mean())
            r[f"pass{k}"] = None if len(z) < 3 else bool(r[f"frac_absz_lt2_{k}"] >= 0.9 and abs(r[f"mean_z{k}"]) <= 0.5)
        crow.append(r)
CV = pd.DataFrame(crow)
CV.to_csv(f"{T}/control_summary{TAG}.csv", index=False)

# ---- comparison table (§3)
C = SER.copy()
for k in (24, 35):
    C[f"sig{k}_ratio"] = C[f"svd_sig{k}_med"] / C[f"raw_sig{k}_med"]
    C[f"sigJ{k}_ratio"] = C[f"svd_sigJ{k}_med"] / C[f"raw_sigJ{k}_med"]
for p in ("raw", "svd"):
    m = IV[IV.path == p].set_index(KEY)
    C = C.merge(m[["rec_med", "rec_min", "rec_max", "leak_fixed_maxabs", "leak_contrast_maxabs", "inj_pass",
                   "transient_episode_rec24", "transient_day_rec24"]].add_prefix(f"{p}_"),
                left_on=KEY, right_index=True, how="left")
    # post-hoc diagnostic (not in SVD.md): error in units of a true change = sigma / step recovery
    C[f"{p}_sig24_eff"] = C[f"{p}_sig24_med"] / C[f"{p}_rec_med"]
C["svd_shrinks"] = (C.sig24_ratio < 0.8) & (C.svd_half_ratio24 <= 1.5) & (C.svd_inj_pass == True)
C.to_csv(f"{T}/compare_raw_vs_svd{TAG}.csv", index=False)

# ---- pooled day-level numbers
res = {}
for k in (24, 35):
    for p in ("raw", "svd"):
        res[f"{p}_sig{k}_median_days"] = float(np.nanmedian(D[f"{p}_sig{k}"]))
        res[f"{p}_sigJ{k}_median_days"] = float(np.nanmedian(D[f"{p}_sigJ{k}"]))
        res[f"{p}_cc{k}_median_days"] = float(np.nanmedian(D[f"{p}_cc{k}"]))
        res[f"{p}_rail{k}_frac_days"] = float(D[f"{p}_eps{k}"].isna().mean())
        res[f"{p}_sd_eps{k}_days"] = float(np.nanstd(D[f"{p}_eps{k}"]))
    both = D[f"raw_sig{k}"].notna() & D[f"svd_sig{k}"].notna()
    res[f"sig{k}_ratio_median_days"] = float(np.median(D.loc[both, f"svd_sig{k}"] / D.loc[both, f"raw_sig{k}"]))
for p in ("raw", "svd"):
    res[f"{p}_lapse_pass"] = int((C[f"{p}_lapse_pass"] == True).sum()); res["n_series"] = len(C)
    res[f"{p}_lapse_slope_median"] = float(C[f"{p}_lapse_slope"].median())
    res[f"{p}_lapse_sign_median"] = float(C[f"{p}_lapse_sign"].median())
    res[f"{p}_inj_pass"] = int((C[f"{p}_inj_pass"] == True).sum())
    res[f"{p}_rec_median"] = float(IV[IV.path == p].rec_med.median())
    res[f"{p}_leak_fixed_maxabs_median"] = float(IV[IV.path == p].leak_fixed_maxabs.median())
    res[f"{p}_half_ratio24_median"] = float(C[f"{p}_half_ratio24"].median())
    res[f"{p}_transient_episode_rec24_median"] = float(C[f"{p}_transient_episode_rec24"].median())
    res[f"{p}_transient_day_rec24_median"] = float(C[f"{p}_transient_day_rec24"].median())
    res[f"{p}_sig24_eff_median"] = float(C[f"{p}_sig24_eff"].median())
for v in ("raw", "own", "ins"):
    g = CV[CV.variant == v]
    res[f"ctrl_{v}_pass24"] = f"{int((g.pass24 == True).sum())}/{int(g.pass24.notna().sum())}"
    res[f"ctrl_{v}_cc24_median"] = float(g.cc24_med.median())
res["rank_gt_nmf_median"] = float(C.rank_gt_nmf.median()); res["sum_w_median"] = float(C.sum_w.median())
res["svd_shrinks_series"] = int(C.svd_shrinks.sum())
json.dump(res, open(f"{T}/summary_svd{TAG}.json", "w"), indent=1)
pd.set_option("display.width", 250)
print(json.dumps(res, indent=1))
print(C[KEY + ["n_days", "rank_gt_nmf", "sum_w", "raw_sig24_med", "svd_sig24_med", "sig24_ratio", "raw_cc24_med", "svd_cc24_med",
               "raw_rail24", "svd_rail24", "raw_lapse_slope", "svd_lapse_slope", "raw_lapse_sign", "svd_lapse_sign",
               "raw_half_ratio24", "svd_half_ratio24", "raw_rec_med", "svd_rec_med", "svd_leak_fixed_maxabs",
               "svd_transient_day_rec24", "raw_sig24_eff", "svd_sig24_eff"]].round(3).to_string())
print(CV.round(3).to_string())
