#!/usr/bin/env python
"""ADDENDUM (reporting only): agent Q's tier A/B/C screen applied to agent O's ARRAY-based 2010
run, so the array-vs-member comparison is direct.

`results/q_screenval/wa2010_{real,null}_rescreen.csv` were computed on agent L's MEMBER-WISE 2010
chain (`results/l_wa2010*/verify_it0.csv`), not on the array chain of `results/o_oly2010/`.  The
tier definitions are agent Q's, unchanged (SCREEN_VALIDATION_NORTH.md sec.5.1):

  tier A = located on >= 5 stations AND |dz| <= 19 km of Slab2 AND >= 3 distinct tremor days
           AND quiet-day FA <= 0.25/template-day with the tremor-day rate >= 3x it
  tier B = tier A + the shuffled-station control
  tier C = tier B + tremor co-timing above the base-rate-matched tremor-hour set

Agent L's bed never had a quiet-day re-fire, so his tier A was computed WITHOUT clause (a); agent
O's bed does have one (and a high one: 6.0 detections/template-day on the real chain against 0.33
on the null).  Both variants are therefore reported here, so the effect of clause (a) is visible.

This is a comparison table, not a verdict: the verdict comes from the literature-grounded
dedupe + location + calibrated-confirmation path in `sf_accept.py` / `sf_confirm.py`.
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sf_common as SF                                                      # noqa: E402


def tiers(V, dz_max, min_days, fa_max, with_a):
    V = V.copy()
    V["t_c1"] = V.c_c1.astype(bool)
    V["t_c4"] = np.array([bool(np.isfinite(x) and abs(x) <= dz_max) for x in V.dz]) & V.t_c1
    V["t_d"] = V.n_days >= min_days
    if with_a:
        V["t_a"] = ((V.fa_quiet_per_tpl_day <= fa_max)
                    & (V.rate_tremor_per_tpl_day >= 3 * V.fa_quiet_per_tpl_day)
                    & (V.n_det > 0))
    else:
        V["t_a"] = True
    V["tierA"] = V.t_c1 & V.t_c4 & V.t_d & V.t_a
    V["tierB"] = V.tierA & V.c_c2.astype(bool)
    V["tierC"] = V.tierB & V.c_b.astype(bool)
    return V


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dz-max", type=float, default=19.0)
    ap.add_argument("--min-days", type=int, default=3)
    ap.add_argument("--fa-max", type=float, default=0.25)
    ap.add_argument("--out", default=f"{SF.RES}/sf_final/tiers_array2010.json")
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    R = dict(dz_max=a.dz_max, min_days=a.min_days, fa_max=a.fa_max)
    for tag, root in (("real", f"{SF.RES}/o_oly2010/2010"),
                      ("null", f"{SF.RES}/o_oly2010/2010_null")):
        V0 = pd.read_csv(f"{root}/screen_final.csv")
        out = {}
        for with_a, key in ((True, "with_a"), (False, "without_a")):
            V = tiers(V0, a.dz_max, a.min_days, a.fa_max, with_a)
            if with_a:
                V.to_csv(f"{SF.RES}/sf_final/tiers_array2010_{tag}.csv", index=False)
            out[key] = {c: int(V[c].sum()) for c in
                        ["t_a", "t_c1", "t_c4", "t_d", "tierA", "tierB", "tierC"]}
        out["n"] = int(len(V0))
        out["as_coded_candidate"] = int(V0.candidate.astype(bool).sum()) \
            if "candidate" in V0.columns else None
        out["network_fa_per_tpl_day"] = float(V0.fa_quiet_per_tpl_day.mean())
        out["network_tremor_per_tpl_day"] = float(V0.rate_tremor_per_tpl_day.mean())
        R[tag] = out
    json.dump(R, open(a.out, "w"), indent=1, default=float)
    print(json.dumps(R, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
