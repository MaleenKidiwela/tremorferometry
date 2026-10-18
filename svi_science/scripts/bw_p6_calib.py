#!/usr/bin/env python
"""bw_p6_calib -- P6 label-free calibration check (rules: results/backward/p6/P6.md, "Calibration check").
Control (series c1) fraction with P > 0.5 per detector x set x year, all days and between days (> 3 d from every Bostock
window), 2000-draw day bootstrap 95 % CI; "clearly above 3 %" = point > 4.5 % and lower bound > 3 %.  P* for flagged sets.
Also the KELB-gain diagnostic (era vs era-noamp on 2003-07-31..2004-06-03) and the unrefit frozen-S15 P on CN4/B1 rows.
Out: p6/tables/calib_check.csv, calib_decision.json, kelb_diag.csv.  No truth read.
"""
import os, sys, json, glob
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_p4_common as P4

BW = "/home/jovyan/tremorferometry/svi_science/results/backward"
P6 = f"{BW}/p6"
PROD = [f"p5f_{y}" for y in range(2003, 2014)] + ["p5e_2003_kelb", "p5e_2004", "p5e_2005", "p5e_2006"]
REF, CLEAR = 0.03, 0.045


def boot_frac(v, days, nb=2000, seed=1):
    g = pd.DataFrame(dict(v=v, d=days)).groupby("d").v.agg(["sum", "count"])
    s, n = g["sum"].to_numpy(float), g["count"].to_numpy(float)
    rng = np.random.default_rng(seed); out = []
    for _ in range(nb):
        i = rng.integers(0, len(s), len(s)); out.append(s[i].sum() / max(n[i].sum(), 1))
    return np.percentile(out, [2.5, 97.5]) if len(s) else (np.nan, np.nan)


def main():
    A = pd.concat([pd.read_parquet(f"{P6}/P/{r}_P.parquet") for r in PROD], ignore_index=True)
    A["detector"] = np.where(A.run.str.startswith("p5e_"), "era", "frozen"); A["year"] = A.date.str[:4]
    near = dict(zip(sorted(A.date.unique()), P4.bostock_near(sorted(A.date.unique()), pad=3)))
    A["between"] = ~A.date.map(near).astype(bool)
    fw = pd.read_parquet("/home/jovyan/tremorferometry/lfe_extension/results/refine/conf/P_by_candidate.parquet", columns=["series", "P_same"])
    ref075 = float((fw[fw.series == "c1"].P_same >= 0.75).mean()); ref05 = float((fw[fw.series == "c1"].P_same > 0.5).mean())
    rows = []
    C = A[A.series == "c1"]
    for (det, ns), g0 in C.groupby(["detector", "net_set"]):
        for yr, g1 in list(g0.groupby("year")) + [("all", g0)]:
            for cls, g in (("all_days", g1), ("between", g1[g1.between])):
                if not len(g): continue
                v = (g.P > 0.5).to_numpy(float); lo, hi = boot_frac(v, g.date.to_numpy())
                r = dict(detector=det, net_set=ns, year=yr, days=cls, n_ctl=len(g), n_days=g.date.nunique(), frac_P_gt05=v.mean(),
                         lo95=lo, hi95=hi, frac_P_ge075=float((g.P >= 0.75).mean()), clearly_above=bool(v.mean() > CLEAR and lo > REF))
                if "P_same_frozenS15" in g and g.P_same_frozenS15.notna().any():
                    r["frac_P_gt05_frozenS15_model"] = float((g.P_same_frozenS15 > 0.5).mean())
                rows.append(r)
    T = pd.DataFrame(rows); T.to_csv(f"{P6}/tables/calib_check.csv", index=False)
    dec = dict(rule="clearly above = point > 4.5 % and day-bootstrap 95 % lower bound > 3 % (all days, all years of the set)",
               forward_prod_ctl_frac_P_ge075=ref075, forward_prod_ctl_frac_P_gt05=ref05, sets={})
    grid = np.round(np.arange(0.5, 1.0, 0.01), 2)
    for (det, ns), g in T[(T.year == "all") & (T.days == "all_days")].groupby(["detector", "net_set"]):
        r = g.iloc[0]; d = dict(frac_P_gt05=r.frac_P_gt05, lo95=r.lo95, hi95=r.hi95, clearly_above=bool(r.clearly_above))
        if r.clearly_above:
            pc = C[(C.detector == det) & (C.net_set == ns)].P.to_numpy()
            ok = [p for p in grid if (pc >= p).mean() <= ref075 + 1e-12]
            d["P_star"] = float(ok[0]) if ok else None
            d["outcome"] = ("clearly above: P is the per-set G1-recipe refit (never the frozen S15 model); rates as ΣP excess only; "
                            "rows flagged ctl_heavy_set; headline count also at P*") if ns in ("CN4", "B1") else \
                           "clearly above: flagged (ctl_heavy_set); rates as ΣP excess only; headline count also at P*"
        else:
            d["outcome"] = "not clearly above: P as is"
        dec["sets"][f"{det}:{ns}"] = d
    json.dump(dec, open(f"{P6}/tables/calib_decision.json", "w"), indent=1, default=float)
    # KELB diagnostic
    E = A[(A.detector == "era") & (A.date >= "2003-07-31") & (A.date < "2004-06-04")]
    kd = []
    for s, g in E.groupby("series"):
        kd.append(dict(series=s, n=len(g), sumP_full=g.P_era_full.sum(), sumP_noamp=g.P_era_noamp.sum(),
                       ratio_noamp_over_full=g.P_era_noamp.sum() / max(g.P_era_full.sum(), 1e-9),
                       n_ge075_full=int((g.P_era_full >= 0.75).sum()), n_ge075_noamp=int((g.P_era_noamp >= 0.75).sum())))
    E2 = A[(A.detector == "era") & (A.date < "2003-07-31")]
    for s, g in E2.groupby("series"):
        kd.append(dict(series=s + " (before 2003-07-31, KELB gain low; noamp used)", n=len(g), sumP_full=g.P_era_full.sum(),
                       sumP_noamp=g.P_era_noamp.sum(), ratio_noamp_over_full=g.P_era_noamp.sum() / max(g.P_era_full.sum(), 1e-9),
                       n_ge075_full=int((g.P_era_full >= 0.75).sum()), n_ge075_noamp=int((g.P_era_noamp >= 0.75).sum())))
    pd.DataFrame(kd).to_csv(f"{P6}/tables/kelb_diag.csv", index=False)
    pd.set_option("display.width", 250)
    print(T[(T.year == "all")].round(4).to_string(index=False)); print(json.dumps(dec, indent=1, default=float)); print(pd.DataFrame(kd).round(3).to_string())


if __name__ == "__main__":
    main()
