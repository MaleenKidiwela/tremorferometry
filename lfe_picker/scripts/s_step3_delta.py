#!/usr/bin/env python
"""P3D STEP 3 analysis -- fit delta from the tail quantiles and apply the PRE-REGISTERED Step-6 rule.

Method (P3 RESULTS 5.5): the null network statistic is standardised by its own spread sigma; the
same-family distribution is the null shifted by delta.  At a null quantile q the same-family rate
above that threshold is r_q, so delta_q = Phi^-1(1-q) - Phi^-1(1-r_q).  delta is the (weighted)
combination over the four tail quantiles actually populated.

PRE-REGISTERED RULE: Step 6 (the template-free pairwise discovery search) runs only if
predicted delta_full = 0.38 sigma x (delta_dense / delta_old) >= 1.5 sigma, i.e. ratio >= 4.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402

QS = [("q_0.001", 1e-3), ("q_0.0001", 1e-4), ("q_1e-5", 1e-5), ("q_1e-6", 1e-6)]
DELTA_OLD_P3 = 0.38          # P3 RESULTS 5.5, 7 sites, 2010+2012 Bostock candidates


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", default="cal_null")
    ap.add_argument("--true", default="cal_true")
    ap.add_argument("--out", default="step3")
    ap.add_argument("--dense-key", default="dense")
    ap.add_argument("--old-key", default="old7")
    a = ap.parse_args()
    null = json.load(open(f"{S.P3D}/pairlink3_{a.null}.json"))
    true = json.load(open(f"{S.P3D}/pairlink3_{a.true}.json"))
    n_same = int(true["same_tot"])
    n_far = int(true["far_tot"])
    rows = []
    for c, sites in true["configs"].items():
        for qn, q in QS:
            hit = int(true["same_hit"][c][qn])
            far = int(true["far_hit"][c][qn])
            r = hit / max(n_same, 1)
            if hit < 3 or r <= 0 or r >= 1:
                d = np.nan
            else:
                d = float(norm.isf(q) - norm.isf(r))
            rows.append(dict(config=c, sites="+".join(sites), q=q, thr=null["thresholds"][c][qn],
                             n_same_hit=hit, n_same_pairs=n_same, rate_same=r,
                             n_far_hit=far, n_far_pairs=n_far, rate_far=far / max(n_far, 1),
                             enrichment=(r / (far / max(n_far, 1))) if far else np.nan,
                             delta_sigma=d))
    D = pd.DataFrame(rows)
    D.to_csv(f"{S.P3D}/{a.out}_delta.csv", index=False)
    fit = D.dropna(subset=["delta_sigma"]).groupby("config").agg(
        delta=("delta_sigma", "median"), n_q=("delta_sigma", "size"),
        sites=("sites", "first")).reset_index()
    fit.to_csv(f"{S.P3D}/{a.out}_delta_fit.csv", index=False)
    f = fit.set_index("config")
    d_old = float(f.delta.get(a.old_key, np.nan))
    d_den = float(f.delta.get(a.dense_key, np.nan))
    ratio = d_den / d_old if np.isfinite(d_old) and d_old > 0 else np.nan
    pred = DELTA_OLD_P3 * ratio
    verdict = "STEP6_LICENSED" if np.isfinite(pred) and pred >= 1.5 else "STEP6_NOT_LICENSED"
    out = dict(variant=true.get("variant", "v1"), episodes=true["episodes"],
               delta_old7=d_old, delta_dense=d_den, ratio=ratio,
               delta_old_p3_reference=DELTA_OLD_P3, predicted_delta_full_sigma=pred,
               rule="Step 6 only if predicted delta_full >= 1.5 sigma (ratio >= 4)",
               verdict=verdict, n_events=true["n_events"], n_pairs=true["n_pairs"],
               n_same_pairs=n_same,
               per_config={r.config: dict(sites=r.sites, delta=float(r.delta), n_q=int(r.n_q))
                           for r in fit.itertuples()})
    json.dump(out, open(f"{S.P3D}/{a.out}_verdict.json", "w"), indent=1, default=str)
    print(D.to_string(index=False))
    print()
    print(fit.sort_values("delta", ascending=False).to_string(index=False))
    print(f"\n[step3] delta_old7 {d_old:.3f} sigma; delta_dense {d_den:.3f} sigma; "
          f"ratio {ratio:.2f}; predicted delta_full = 0.38 x {ratio:.2f} = {pred:.2f} sigma "
          f"-> {verdict}")


if __name__ == "__main__":
    main()
