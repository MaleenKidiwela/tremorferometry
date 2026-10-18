#!/usr/bin/env python
"""False share (scrambled control / data) by confidence cut-off, quiet vs tremor days, for Shelly and Ducellier (2026-10-02).
Quiet day = no PNSN tremor window within R km of the region's family centroid that day (PNSN catalogue starts 2010-01-06)."""
import numpy as np, pandas as pd
ROOT = "/home/jovyan/tremorferometry"; R = f"{ROOT}/lfe_south/results"
T = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
T["date"] = pd.to_datetime(T.time, utc=True).dt.tz_localize(None).dt.normalize()
fam = pd.read_csv(f"{R}/families_south.csv")
def hav(la1, lo1, la2, lo2):
    p = np.pi / 180; a = np.sin((la2 - la1) * p / 2) ** 2 + np.cos(la1 * p) * np.cos(la2 * p) * np.sin((lo2 - lo1) * p / 2) ** 2
    return 12742 * np.arcsin(np.sqrt(a))
def tremor_days(src, rkm):
    f = fam[fam.source.str.contains(src, case=False)]
    d = np.min(np.vstack([hav(T.lat.values, T.lon.values, la, lo) for la, lo in zip(f.lat, f.lon)]), 0)
    return T[d <= rkm].groupby("date").size()
def table(D, pcol, trem, label, t0="2010-01-06", t1="2026-05-27"):
    D = D[(D.date >= t0) & (D.date <= t1)].copy()
    D["tremor"] = D.date.map(trem).fillna(0)
    out = []
    for cls, g in [("no tremor", D[D.tremor == 0]), ("tremor day", D[D.tremor > 0])]:
        nd = g[g.series == "data"].date.nunique(); nc = g[g.series == "c1"].date.nunique()
        for th in [0.5, 0.6, 0.7, 0.75, 0.8, 0.9]:
            a = g[(g.series == "data") & (g[pcol] >= th)]; c = g[(g.series == "c1") & (g[pcol] >= th)]
            cn = len(c) * (nd / nc if nc else np.nan)   # per-day normalised to the data's days
            out.append(dict(region=label, days=cls, n_days=nd, cutoff=th, data=len(a), control=len(c),
                            false_share=round(cn / max(len(a), 1), 3), families=a.fam.nunique()))
    return pd.DataFrame(out)
res = []
# Shelly: KMR-fixed v2 model (P_same_v2) on all v3b + forward candidates
S = pd.read_parquet(f"{R}/sh_kmrfix/conf_v2_rows.parquet", columns=["series", "date", "fam", "P_same_v2", "P_same_FULL"])
S["date"] = pd.to_datetime(S.date)
for rkm in (30, 60):
    res.append(table(S, "P_same_v2", tremor_days("shelly", rkm), f"Shelly (PNSN ≤{rkm} km)"))
# Ducellier: final model P_same
D = pd.read_parquet(f"{R}/duc_fin/candidates_conf.parquet", columns=["series", "date", "fam", "P_same", "role"])
D["date"] = pd.to_datetime(D.date)
for rkm in (30, 60):
    res.append(table(D, "P_same", tremor_days("ducellier|plourde", rkm), f"Ducellier (PNSN ≤{rkm} km)"))
X = pd.concat(res); X.to_csv(f"{R}/cutoff_check/false_share_by_cutoff.csv", index=False)
pd.set_option("display.width", 200); print(X.to_string(index=False))
