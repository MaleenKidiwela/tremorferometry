#!/usr/bin/env python
"""WA pair-stage statistics.  Prefix `l_`.

There is no catalogue truth in Washington, so the two north statistics that need truth (true-pair
fraction, families seen) are replaced by the two that do not:
  * the TREMOR-vs-QUIET pair-count ratio at 8 x MAD -- detections per station-hour on PNSN tremor
    hours against the 40 hour-matched, earthquake-screened quiet hours of WA_V2_INPUTS sec.2.1;
  * the nmad distribution and the top-tier counts the north reports.
Also merges the per-hour parquets into one detection file, applying k_group's own block-quality
rule (a block whose 8 x MAD bar is not > 3 sd is a data-quality failure, not a detection set).
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import l_common as L                                                        # noqa: E402


def merge(d):
    A = [pd.read_parquet(f) for f in sorted(glob.glob(f"{d}/*.parquet"))
         if not os.path.basename(f).startswith("_")]
    A = [x for x in A if len(x)]
    R = pd.concat(A, ignore_index=True) if A else pd.DataFrame(
        columns=["i", "j", "cc", "nmad", "t_i", "t_j", "date", "hour"])
    S = pd.read_csv(f"{d}/_stats.csv")
    bad = {(r.date, int(r.hour)) for r in S.itertuples() if not (r.nmad_in_sd > 3)}
    if bad and len(R):
        R = R[~pd.Series(list(zip(R.date, R.hour)), index=R.index).isin(bad)]
    R = R[np.isfinite(R.nmad)].reset_index(drop=True)
    return R, S, bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tremor", required=True)
    ap.add_argument("--quiet", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    RT, ST, badT = merge(a.tremor)
    RT.to_parquet(f"{a.out}_all.parquet", index=False)
    okT = ST[~ST.set_index(["date", "hour"]).index.isin(badT)]
    rep = dict(n_blocks_tremor=int(len(ST)), n_blocks_tremor_clean=int(len(okT)),
               n_blocks_tremor_rejected=int(len(badT)),
               n_det_tremor=int(len(RT)),
               nmad_in_sd_median=float(ST.nmad_in_sd.median()),
               sta_per_block_median=float(ST.n_sta.median()),
               det_per_hour_tremor=float(len(RT) / max(len(okT), 1)))
    if len(RT):
        rep["nmad"] = {k: float(np.quantile(RT.nmad, q)) for k, q in
                       (("p50", .5), ("p90", .9), ("p99", .99), ("max", 1.0))}
        for n in (100, 500, 2000, 10000, 50000):
            sub = RT.nlargest(min(n, len(RT)), "nmad")
            rep[f"top{n}_nmad_min"] = float(sub.nmad.min())
            rep[f"top{n}_hours"] = int(sub.groupby(["date", "hour"]).ngroups)
            rep[f"top{n}_days"] = int(sub.date.nunique())
    if a.quiet and os.path.exists(f"{a.quiet}/_stats.csv"):
        RQ, SQ, badQ = merge(a.quiet)
        RQ.to_parquet(f"{a.out}_quiet.parquet", index=False)
        okQ = SQ[~SQ.set_index(["date", "hour"]).index.isin(badQ)]
        dq = len(RQ) / max(len(okQ), 1)
        rep.update(n_blocks_quiet=int(len(SQ)), n_blocks_quiet_clean=int(len(okQ)),
                   n_det_quiet=int(len(RQ)), det_per_hour_quiet=float(dq),
                   tremor_quiet_ratio=float(rep["det_per_hour_tremor"] / dq) if dq > 0 else np.nan,
                   sta_per_block_median_quiet=float(SQ.n_sta.median()),
                   nmad_in_sd_median_quiet=float(SQ.nmad_in_sd.median()))
        # the same ratio restricted to the strongest tier, where the north's purity lives
        for n in (2000,):
            thr = float(RT.nlargest(min(n, len(RT)), "nmad").nmad.min()) if len(RT) else np.nan
            if np.isfinite(thr):
                rep[f"tremor_quiet_ratio_above_nmad{n}"] = float(
                    ((RT.nmad >= thr).sum() / max(len(okT), 1))
                    / max((RQ.nmad >= thr).sum() / max(len(okQ), 1), 1e-9))
    # every pair is inside ONE hour block by construction, so the pair-lag/tremor-hour
    # coincidence is 1.0 at this stage; the informative version is the grown-detection one in
    # l_verify.py (fraction of matched-filter detections inside a PNSN tremor-active hour).
    rep["pair_lag_in_tremor_hour"] = 1.0
    rep["pair_lag_in_tremor_hour_note"] = "1.0 by construction: an autocorrelation block is one hour"
    json.dump(rep, open(f"{a.out}_pairstats.json", "w"), indent=1, default=float)
    L.log(json.dumps(rep, indent=1, default=float))


if __name__ == "__main__":
    main()
