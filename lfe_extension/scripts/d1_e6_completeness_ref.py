#!/usr/bin/env python
"""E6 step 5 -- per-family COMPLETENESS table for the 2012 catalog at the adopted operating point: own-family recall of the Bostock
TEST-half events (+-2 s) per family with a Wilson 95 % interval, the family's template coverage, build-event count, median magnitude
and hypocentre; plus the family's quiet-day detection rate (the floor its completeness-weighted rate must subtract).
Writes results/e6_family_completeness_2012.csv (from results/e6_score_<tag>_recall_family.csv and e5_budget_<tag>.csv)."""
from __future__ import annotations
import argparse
import numpy as np, pandas as pd
RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"


def wilson(k, n, z=1.96):
    k = np.asarray(k, float); n = np.asarray(n, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d
        h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return np.where(n > 0, c - h, np.nan), np.where(n > 0, c + h, np.nan)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True); ap.add_argument("--T", type=float, required=True)
    a = ap.parse_args()
    rf = pd.read_csv(f"{RES}/e6_score_{a.tag}_recall_family.csv"); rf = rf[(rf.L == 8.0) & (rf["T"] == a.T)].copy()
    bud = pd.read_csv(f"{RES}/e5_budget_{a.tag}.csv"); b = bud[(bud.L == 8.0) & (bud["T"] == a.T)]
    quiet = b[b.doy.isin(range(228, 243))].groupby("fam").n_det.mean(); post = b[b.doy.isin(range(270, 280))].groupby("fam").n_det.mean()
    qsite = b[b.doy.isin(range(228, 243))].groupby("fam").n_det_roll_site.mean() if "n_det_roll_site" in b.columns else None
    F = pd.read_csv(f"{RES}/e5r_families_2012.csv").set_index("bkey")
    k = np.rint(rf.recall_own_2s.fillna(0) * rf.n).astype(int); lo, hi = wilson(k, rf.n)
    out = pd.DataFrame(dict(bkey=rf.bkey.astype(int), sb_family_id=rf.bkey.map(F.sb_fam_id), lat=rf.bkey.map(F.lat), lon=rf.bkey.map(F.lon), dep_km=rf.bkey.map(F.dep),
                            n_test_events=rf.n.astype(int), n_own_hits=k, own_recall=rf.recall_own_2s, own_recall_ci_lo=lo, own_recall_ci_hi=hi,
                            n_stations_template=rf.nsta.astype(int), n_build_events=rf.n_build.astype(int), mag_median_test=rf.mag_med,
                            quiet_fa_per_day=rf.bkey.astype(float).map(quiet), quiet_site_rolled_per_day=rf.bkey.astype(float).map(qsite) if qsite is not None else np.nan,
                            post_sse_per_day=rf.bkey.astype(float).map(post), operating_tag=a.tag, operating_T=a.T, template_len_s=8))
    out = out.sort_values("bkey").reset_index(drop=True)
    out.to_csv(f"{RES}/e6_family_completeness_2012.csv", index=False)
    t = out[out.n_stations_template > 0]
    print(f"[completeness] {len(out)} families ({len(t)} templated) at {a.tag} T {a.T:g}: own recall p10/50/90 {t.own_recall.quantile([.1, .5, .9]).round(3).tolist()}, "
          f"event-weighted {np.average(t.own_recall, weights=t.n_test_events):.3f}; families >= 0.5: {int((t.own_recall >= 0.5).sum())}, <= 0.05: {int((t.own_recall <= 0.05).sum())}; "
          f"median CI half-width {((t.own_recall_ci_hi - t.own_recall_ci_lo) / 2).median():.3f}; quiet FA / day p50/90 {t.quiet_fa_per_day.quantile([.5, .9]).round(2).tolist()}")
    print(f"[write] results/e6_family_completeness_2012.csv")


if __name__ == "__main__":
    main()
