#!/usr/bin/env python
"""E6 steps 1-3 -- the FA / recall decision table for one e4_template_test.py --budget --sites run (CPU, seconds).
Per threshold T: detections per template per day (data, station-rolled, site-rolled) on the QUIET window (2012-08-15..08-29, doy
228-242) and the POST-SSE window (09-26..10-05, doy 270-279); ETS-day (>= 100 Bostock events) rates, rolled fractions and own
precision; own-family recall (+-2 s) of the Bostock TEST-half events overall and by magnitude bin (like e5r_score.py, from the
detections parquet); the pre-registered E6 bars (audit prescription 2) evaluated on BOTH windows and the operating threshold =
smallest T in {7.5, 8, 8.5, 9} meeting all of them.  Bars are REPORTED, never adjusted.
  FA bars   quiet <= 1 detection / template / day;  site-rolled / data <= 1.3 (the literal wording; the symmetric version
            max(site/data, data/site) <= 1.3 is reported alongside as 'ratio_sym');  ETS-day site-rolled fraction <= 5 %.
  recall    own >= 0.25 overall, >= 0.20 at M1.6-1.8 (deciding); M1.8-2.0 >= 0.35 and M>=2.0 >= 0.60 on record.
Writes results/e6_fa_table_<tag>.csv and results/e6_recall_<tag>.csv; prints both."""
from __future__ import annotations
import argparse, json
import numpy as np, pandas as pd
RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
QUIET = set(range(228, 243)); POST = set(range(270, 280))
MAGBINS = [(-9, 1.6, "M<1.6"), (1.6, 1.8, "M1.6-1.8"), (1.8, 2.0, "M1.8-2.0"), (2.0, 9, "M>=2.0")]
T_OPS = [7.5, 8.0, 8.5, 9.0]


def min_absdiff(a, b_sorted):
    a = np.asarray(a, float)
    if len(a) == 0:
        return np.zeros(0)
    if len(b_sorted) == 0:
        return np.full(len(a), np.inf)
    j = np.searchsorted(b_sorted, a)
    lo = b_sorted[np.clip(j - 1, 0, len(b_sorted) - 1)]; hi = b_sorted[np.clip(j, 0, len(b_sorted) - 1)]
    return np.minimum(np.abs(a - lo), np.abs(a - hi))


def own_recall(det, test, T):
    hit = np.zeros(len(test), bool)
    for fk, g in test.groupby("bkey"):
        td = np.sort(det[(det.fam == float(fk)) & (det.net_mad >= T)].t.to_numpy())
        hit[test.index.get_indexer(g.index)] = min_absdiff(g.t_origin.to_numpy(), td) <= 2.0
    return hit


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--L", type=float, default=8.0)
    a = ap.parse_args()
    gt = pd.read_parquet(f"{RES}/e5r_gt_bostock2012.parquet").reset_index(drop=True)
    gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    per_day = gt.groupby("doy").size(); ETS = set(per_day[per_day >= 100].index)
    test = gt[(gt.half == "test") & (~gt.guard)].copy().reset_index(drop=True)
    bud = pd.read_csv(f"{RES}/e5_budget_{a.tag}.csv"); bud = bud[bud.L == a.L]
    det = pd.read_parquet(f"{RES}/e5r_detections_{a.tag}.parquet"); det = det[det.L == a.L]
    has_site = "n_det_roll_site" in bud.columns
    rows = []
    for T in sorted(bud["T"].unique()):
        b = bud[bud["T"] == T]
        r = dict(tag=a.tag, L=a.L, T=T, n_family_days=int(len(b)))
        for name, days in (("quiet", QUIET), ("post", POST), ("ets", ETS)):
            x = b[b.doy.isin(days)]
            r[f"{name}_ndays"] = int(x.doy.nunique())
            r[f"{name}_data"] = float(x.n_det.mean()); r[f"{name}_data_p90"] = float(x.n_det.quantile(.9))
            r[f"{name}_roll_sta"] = float(x.n_det_roll.mean())
            r[f"{name}_roll_site"] = float(x.n_det_roll_site.mean()) if has_site else np.nan
            r[f"{name}_site_over_data"] = r[f"{name}_roll_site"] / max(r[f"{name}_data"], 1e-9)
            r[f"{name}_data_over_site"] = r[f"{name}_data"] / max(r[f"{name}_roll_site"], 1e-9)
            r[f"{name}_ratio_sym"] = max(r[f"{name}_site_over_data"], r[f"{name}_data_over_site"])
            r[f"{name}_sta_over_data"] = r[f"{name}_roll_sta"] / max(r[f"{name}_data"], 1e-9)
            r[f"{name}_frac_roll_site"] = float(x.n_det_roll_site.sum() / max(x.n_det.sum(), 1)) if has_site else np.nan
            r[f"{name}_frac_roll_sta"] = float(x.n_det_roll.sum() / max(x.n_det.sum(), 1))
            r[f"{name}_precision_own"] = float(x.n_match_own.sum() / max(x.n_det.sum(), 1))
            r[f"{name}_precision_any"] = float(x.n_match_any.sum() / max(x.n_det.sum(), 1))
        # pre-registered bars on BOTH windows (quiet and post-SSE) + ETS site-rolled fraction
        r["bar_quiet_le1"] = bool(r["quiet_data"] <= 1.0); r["bar_post_le1"] = bool(r["post_data"] <= 1.0)
        r["bar_quiet_site_over_data_le1.3"] = bool(r["quiet_site_over_data"] <= 1.3); r["bar_post_site_over_data_le1.3"] = bool(r["post_site_over_data"] <= 1.3)
        r["bar_quiet_sym_le1.3"] = bool(r["quiet_ratio_sym"] <= 1.3); r["bar_post_sym_le1.3"] = bool(r["post_ratio_sym"] <= 1.3)
        r["bar_ets_site_frac_le5pct"] = bool(r["ets_frac_roll_site"] <= 0.05)
        r["fa_bars_met_literal"] = bool(r["bar_quiet_le1"] and r["bar_post_le1"] and r["bar_quiet_site_over_data_le1.3"] and r["bar_post_site_over_data_le1.3"] and r["bar_ets_site_frac_le5pct"])
        r["fa_bars_met_sym"] = bool(r["bar_quiet_le1"] and r["bar_post_le1"] and r["bar_quiet_sym_le1.3"] and r["bar_post_sym_le1.3"] and r["bar_ets_site_frac_le5pct"])
        if T in T_OPS or T in (7.0, 10.0):
            hit = own_recall(det, test, T)
            r["recall_own_all"] = float(hit.mean())
            for lo, hi, nm in MAGBINS:
                k = (test.mag >= lo) & (test.mag < hi); r[f"recall_own_{nm}"] = float(hit[k].mean())
            r["n_catalog_ge_T"] = int((det.net_mad >= T).sum())
            r["bar_recall_all_ge0.25"] = bool(r["recall_own_all"] >= 0.25); r["bar_recall_M16_18_ge0.20"] = bool(r["recall_own_M1.6-1.8"] >= 0.20)
            r["bar_recall_M18_20_ge0.35"] = bool(r["recall_own_M1.8-2.0"] >= 0.35); r["bar_recall_Mge2_ge0.60"] = bool(r["recall_own_M>=2.0"] >= 0.60)
            r["merge_kept_sensitivity_ge0.20"] = bool(r["recall_own_all"] >= 0.20)
        rows.append(r)
    D = pd.DataFrame(rows)
    ops_lit = [t for t in T_OPS if bool(D[D["T"] == t].fa_bars_met_literal.iloc[0])]
    ops_sym = [t for t in T_OPS if bool(D[D["T"] == t].fa_bars_met_sym.iloc[0])]
    D["operating_T_literal"] = ops_lit[0] if ops_lit else np.nan; D["operating_T_sym"] = ops_sym[0] if ops_sym else np.nan
    D.to_csv(f"{RES}/e6_fa_table_{a.tag}.csv", index=False)
    pd.set_option("display.width", 300)
    cols = ["T", "quiet_data", "quiet_roll_sta", "quiet_roll_site", "quiet_site_over_data", "post_data", "post_roll_sta", "post_roll_site", "post_site_over_data",
            "ets_data", "ets_frac_roll_sta", "ets_frac_roll_site", "ets_precision_own", "fa_bars_met_literal", "fa_bars_met_sym"]
    print(f"\n=== {a.tag} L{a.L:g}: detections / template / day (data | station-rolled | site-rolled); ETS = {len(ETS)} days with >= 100 Bostock events")
    print(D[cols].round(3).to_string(index=False))
    rc = [c for c in D.columns if c.startswith("recall_own") or c.startswith("bar_recall") or c == "merge_kept_sensitivity_ge0.20"]
    print(D.loc[D.recall_own_all.notna(), ["T"] + rc].round(3).to_string(index=False))
    print(f"operating T (literal site/data bar): {ops_lit[0] if ops_lit else 'NONE below 9.5'};  (symmetric ratio bar): {ops_sym[0] if ops_sym else 'NONE'}")
    print(f"[write] results/e6_fa_table_{a.tag}.csv")


if __name__ == "__main__":
    main()
