#!/usr/bin/env python
"""STEP 1c' STAGE 5 -- T_pair by the frozen FA rule, with a MANDATORY leave-one-hour-out column.

The rule is unchanged (pooled <= 1 false pair per quiet hour AND zero quiet components of size
>= 5, 0.25 grid from the 5.0 engineering floor).  What is new is the LOO column Merlin requires:
T recomputed with each quiet hour removed in turn, so a single transient can never set T unseen
again (in 1c one hour set T = 12.50 alone, LOO 7.50).

ACCEPT: T binds by the FA rule between 6.5 and 8.5 AND the LOO range is <= 1.0.
Otherwise the queue stops and names the dominating hour.
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
import w_common as W                                                        # noqa: E402
import w_cluster as C                                                       # noqa: E402
import y_step1diag as Y                                                     # noqa: E402

T_LO, T_HI = 6.5, 8.5
LOO_RANGE_MAX = 1.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quiet-pairs", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-quiet-hours", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    P = pd.read_parquet(a.quiet_pairs)
    nq = a.n_quiet_hours or json.load(
        open(a.quiet_pairs.replace(".parquet", "_meta.json")))["n_hours"]
    T, Tt, binds = C.threshold(P, nq, out=a.out)

    hours = sorted({(str(d), int(h)) for d, h in zip(P.date, P.hour)}) if len(P) else []
    QH = pd.read_csv(W.QUIET_HOURS_CSV) if (W.QUIET_HOURS_CSV
                                            and os.path.exists(W.QUIET_HOURS_CSV)) else None
    allh = ([(str(r.date), int(r.hour)) for r in QH.itertuples()] if QH is not None
            else hours)
    rows = []
    for (d, h) in allh:
        Q = P[~((P.date == d) & (P.hour == h))] if len(P) else P
        t_loo, b_loo = Y.fast_threshold(Q, max(nq - 1, 1))
        rows.append(dict(date=d, hour=int(h),
                         n_pairs_this_hour=int(((P.date == d) & (P.hour == h)).sum())
                         if len(P) else 0,
                         T_without_this_hour=float(t_loo), binds_without=b_loo))
    L = pd.DataFrame(rows)
    L.to_csv(a.out + "_loo.csv", index=False)
    lo = float(L.T_without_this_hour.min()) if len(L) else float("nan")
    hi = float(L.T_without_this_hour.max()) if len(L) else float("nan")
    rng = hi - lo if np.isfinite(lo) else float("nan")
    dom = None
    if len(L):
        k = int(np.argmin(L.T_without_this_hour.to_numpy()))
        if L.T_without_this_hour.iloc[k] < T - 1e-9:
            dom = f"{L.date.iloc[k]} {int(L.hour.iloc[k]):02d}h " \
                  f"(T without it {L.T_without_this_hour.iloc[k]:.2f} vs {T:.2f})"
    okT = (T_LO - 1e-9) <= T <= (T_HI + 1e-9) and binds == "fa_rule"
    okL = np.isfinite(rng) and rng <= LOO_RANGE_MAX + 1e-9
    j = json.load(open(a.out + ".json"))
    j.update(loo_T_min=lo, loo_T_max=hi, loo_range=rng, dominating_hour=dom,
             n_loo_hours=int(len(L)),
             accept=dict(T_range=[T_LO, T_HI], loo_range_le=LOO_RANGE_MAX),
             verdict=("ACCEPT" if (okT and okL) else "BELOW_ACCEPT"),
             reason=("" if (okT and okL) else
                     (f"T = {T:.2f} ({binds}) outside [{T_LO}, {T_HI}]; " if not okT else "")
                     + (f"LOO range {rng:.2f} > {LOO_RANGE_MAX}"
                        + (f"; dominated by {dom}" if dom else "") if not okL else "")))
    json.dump(j, open(a.out + ".json", "w"), indent=1, default=str)
    print(json.dumps(j, indent=1, default=str))


if __name__ == "__main__":
    main()
