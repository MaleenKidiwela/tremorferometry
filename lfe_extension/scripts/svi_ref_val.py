#!/usr/bin/env python
"""SVI_REF inter-ETS validation (answer-key-free; Merlin).  Bostock's catalogue covers ETS windows only, so P between
episodes cannot be checked against it.  Three checks that need no answer key:

  1 withheld S-energy   On quiet days (no PNSN tremor in the family box that day, >= 3 days from any tremor day): mean
                        S-window energy excess of the 2 LOWEST-WEIGHT sites (conf_ze: in no model) per bin of the headline
                        P, data vs the scrambled control.  Real repeats should carry S energy at sites that hardly
                        triggered them; chance coincidences should not.  Pass = data conf_ze at P >= 0.5 exceeds the
                        control's (all control rows) with a day-bootstrap 95 % CI above 0, and rises with P (Spearman > 0).
  2 control sum of P    per month: the scrambled control's sum of P (expected about 0 and flat).  Pass = in every network
                        epoch (E0-E3) the control sum of P per day is <= 5 % of the data's, and the monthly control series
                        has no trend (|Spearman vs time| < 0.3).
  3 station-set         the S13 re-run of E0 months (LRIV and SNB removed while they were working) scored by the same model:
    invariance          sum of P(same) x his/sumP(S13) against the production sum of P(same) x his/sumP(S15), pooled
                        ratio and median monthly ratio both within 0.8-1.2.
If any check fails, inter-ETS rates are reported only as a control-subtracted bracket (svi_ref_rates.py reads the verdict).
Outputs OUT/val/: val_energy_bins.csv, val_ctl_monthly.csv, val_invariance.csv, val_summary.json.

  python scripts/svi_ref_val.py
"""
from __future__ import annotations
import os, sys, json, warnings
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R
import c_common as C

OV = f"{R.OUT}/val"
EDGES = [0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0001]


def boot_mean_diff(a_val, a_day, b_val, b_day, nb=300, seed=5):
    rs = np.random.default_rng(seed)
    da = pd.Series(a_val).groupby(a_day).agg(["sum", "count"]); db = pd.Series(b_val).groupby(b_day).agg(["sum", "count"])
    days = np.array(sorted(set(da.index) | set(db.index)))
    da = da.reindex(days).fillna(0); db = db.reindex(days).fillna(0)
    out = []
    for _ in range(nb):
        ii = rs.integers(0, len(days), len(days))
        sa, na = da["sum"].values[ii].sum(), da["count"].values[ii].sum(); sb, nb_ = db["sum"].values[ii].sum(), db["count"].values[ii].sum()
        if na and nb_: out.append(sa / na - sb / nb_)
    return (float(np.mean(out)), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))) if out else (np.nan,) * 3


def main():
    warnings.filterwarnings("ignore"); R.ensure(OV)
    from scipy.stats import spearmanr
    P = pd.read_parquet(f"{R.OUT}/conf/P_by_candidate.parquet")
    P = P[P.run == json.load(open(f"{R.OUT}/conf/runs.json"))["prod"]]
    summ = {}
    # ---- 1 withheld S-energy on quiet days
    Q = P[P.quiet.fillna(False).astype(bool) & P.conf_ze.notna()].copy()
    Q["bin"] = pd.cut(Q.P_same, EDGES, right=False)
    rows = []
    for b, g in Q.groupby("bin", observed=False):
        d = g[g.series == "data"]; c = g[g.series == "c1"]
        rows.append(dict(P_bin=str(b), n_data=len(d), n_ctl=len(c), conf_ze_data=float(d.conf_ze.mean()) if len(d) else np.nan,
                         conf_ze_ctl=float(c.conf_ze.mean()) if len(c) else np.nan, n_days_data=d.date.nunique()))
    Eb = pd.DataFrame(rows); Eb.to_csv(f"{OV}/val_energy_bins.csv", index=False)
    hi = Q[(Q.series == "data") & (Q.P_same >= 0.5)]; ctl = Q[Q.series == "c1"]
    d = boot_mean_diff(hi.conf_ze.values, hi.date.values, ctl.conf_ze.values, ctl.date.values) if len(hi) >= 30 else (np.nan,) * 3
    ok = Eb[Eb.n_data >= 30]
    rho = float(spearmanr(np.arange(len(ok)), ok.conf_ze_data).correlation) if len(ok) >= 3 else np.nan
    summ["energy"] = dict(n_quiet_data=int((Q.series == "data").sum()), n_quiet_ctl=len(ctl), n_data_P_ge_05=len(hi),
                          diff_hiP_minus_ctl=d, spearman_bins=rho,
                          **{"pass": bool(np.isfinite(d[1]) and d[1] > 0 and np.isfinite(rho) and rho > 0)})
    # ---- 2 scrambled-control sum of P per month
    g = P.groupby(["date", "series"]).P_same.sum().unstack(fill_value=0.0).reset_index()
    for c_ in ("data", "c1"):
        if c_ not in g: g[c_] = 0.0
    g["month"] = g.date.str[:7]
    Mo = g.groupby("month").agg(days=("date", "nunique"), sumP_data=("data", "sum"), sumP_ctl=("c1", "sum")).reset_index()
    Mo["data_per_day"] = Mo.sumP_data / Mo.days; Mo["ctl_per_day"] = Mo.sumP_ctl / Mo.days
    Mo["ctl_over_data"] = Mo.sumP_ctl / Mo.sumP_data.replace(0, np.nan)
    Mo["epoch"] = Mo.month.map(C.epoch_of_month)
    Mo.to_csv(f"{OV}/val_ctl_monthly.csv", index=False)
    ep = Mo.groupby("epoch").agg(data=("sumP_data", "sum"), ctl=("sumP_ctl", "sum"), days=("days", "sum"))
    ep["ctl_over_data"] = ep.ctl / ep.data.replace(0, np.nan)
    rho2 = float(spearmanr(np.arange(len(Mo)), Mo.ctl_per_day).correlation) if len(Mo) >= 6 else np.nan
    summ["control_sumP"] = dict(by_epoch=ep.round(4).reset_index().to_dict("records"), spearman_time=rho2,
                                **{"pass": bool((ep.ctl_over_data <= 0.05).all() and np.isfinite(rho2) and abs(rho2) < 0.3)})
    # ---- 3 station-set invariance (S13 re-run of E0 months)
    cfg = json.load(open(f"{R.OUT}/conf/runs.json"))
    A = pd.read_parquet(f"{R.OUT}/conf/candidates_conf.parquet",
                        columns=["run", "kind", "sset_run", "series", "date", "year", "row_role", "P_same", "m_same"])
    fit = A[(A.row_role == "fit") & (A.series == "data") & A.m_same.notna()]
    kset = fit.groupby("sset_run").apply(lambda x: x.m_same.sum() / max(x.P_same.sum(), 1e-9)).to_dict()
    inv = dict(available=False)
    if cfg.get("prod_S13"):
        S = A[(A.run == cfg["prod_S13"]) & (A.series == "data")]
        Pm = P[(P.series == "data") & P.date.isin(set(S.date))]
        a1 = S.groupby(S.date.str[:7]).P_same.sum() * kset.get("S13", np.nan)
        a0 = Pm.groupby(Pm.date.str[:7]).P_same.sum() * kset.get("S15", np.nan)
        T = pd.DataFrame(dict(S13_scaled=a1, prod_scaled=a0)).dropna()
        T["ratio"] = T.S13_scaled / T.prod_scaled
        T.to_csv(f"{OV}/val_invariance.csv")
        pooled = float(T.S13_scaled.sum() / T.prod_scaled.sum()) if len(T) else np.nan
        inv = dict(available=True, k_S15=kset.get("S15"), k_S13=kset.get("S13"), n_months=len(T), pooled_ratio=pooled,
                   median_monthly_ratio=float(T.ratio.median()) if len(T) else np.nan)
        inv["pass"] = bool(np.isfinite(pooled) and 0.8 <= pooled <= 1.2 and 0.8 <= inv["median_monthly_ratio"] <= 1.2)
    else:
        inv["pass"] = False
    summ["invariance"] = inv
    summ["all_pass"] = bool(summ["energy"]["pass"] and summ["control_sumP"]["pass"] and inv["pass"])
    summ["inter_ets_mode"] = "rate" if summ["all_pass"] else "control_subtracted_bracket"
    json.dump(summ, open(f"{OV}/val_summary.json", "w"), indent=1, default=float)
    R.log(f"[val] {json.dumps(summ, default=float)[:2000]}")


if __name__ == "__main__":
    main()
