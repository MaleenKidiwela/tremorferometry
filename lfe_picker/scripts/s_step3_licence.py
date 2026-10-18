#!/usr/bin/env python
"""P3D Step-3 LICENCE statistic with the protocol fix (P3D agent, 2026-09-19 addendum).

The licence rule as originally written — Step 6 proceeds if
delta_full = 0.38 sigma x (delta_dense / delta_old) >= 1.5 sigma — is evaluated on ONE null
realisation.  This session measured that on a 1,384-candidate episode the denominator swings +-49 %
between null realisations and the ratio crosses the bar of 4 in 2 of 5 draws with no change in the
data.  The fix, applied here:

  (i)  MINIMUM POPULATED-QUANTILE COUNT: a configuration's delta is admissible in a realisation only
       if at least `--min-quantiles` of the four tail quantiles have >= `--min-hits` same-family
       pairs above threshold (default 3 quantiles, >= 3 pairs).  Otherwise that realisation
       contributes nothing for that configuration and is reported as unusable.
  (ii) NULL-REALISATION AVERAGE: delta is the mean over admissible realisations, and the licence
       statistic is mean(delta_dense) / mean(delta_old), with the per-realisation spread reported.
       The rule FIRES only if the bar is met by the mean AND by the 2.5th percentile of the
       bootstrap over realisations, so a single lucky draw cannot license the arm.
  (iii) MINIMUM ADMISSIBLE REALISATIONS: at least `--min-runs` (default 3) realisations must be
       admissible before the rule may be adjudicated at all.  Without this, (ii) degenerates -- a
       bootstrap over one realisation returns that realisation, which is exactly the failure mode
       the fix exists to prevent.

  python scripts/s_step3_licence.py --runs step3_v2dip_s18,step3_v2dip_s18_sh1200,... --label V2_S18
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402

DELTA_OLD_P3 = 0.38
BAR_SIGMA = 1.5


def delta_of(df, cfg, min_q, min_hits):
    g = df[(df.config == cfg) & df.delta_sigma.notna() & (df.n_same_hit >= min_hits)]
    if len(g) < min_q:
        return np.nan, len(g)
    return float(g.delta_sigma.median()), len(g)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", required=True)
    ap.add_argument("--label", default="run")
    ap.add_argument("--dense-key", default="dense")
    ap.add_argument("--old-key", default="old7")
    ap.add_argument("--min-quantiles", type=int, default=3)
    ap.add_argument("--min-hits", type=int, default=3)
    ap.add_argument("--min-runs", type=int, default=3)
    a = ap.parse_args()
    rows = []
    for t in a.runs.split(","):
        p = f"{S.P3D}/{t}_delta.csv"
        if not os.path.exists(p):
            print(f"  [licence] {t}: missing")
            continue
        D = pd.read_csv(p)
        d_old, nq_old = delta_of(D, a.old_key, a.min_quantiles, a.min_hits)
        d_den, nq_den = delta_of(D, a.dense_key, a.min_quantiles, a.min_hits)
        rows.append(dict(run=t, delta_old=d_old, nq_old=nq_old, delta_dense=d_den,
                         nq_dense=nq_den,
                         ratio=(d_den / d_old if np.isfinite(d_old) and d_old > 0 else np.nan),
                         admissible=bool(np.isfinite(d_old) and np.isfinite(d_den))))
    R = pd.DataFrame(rows)
    R.to_csv(f"{S.P3D}/step3_licence_{a.label}.csv", index=False)
    ok = R[R.admissible]
    out = dict(label=a.label, n_runs=len(R), n_admissible=int(len(ok)),
               min_quantiles=a.min_quantiles, min_hits=a.min_hits,
               per_run=R.to_dict("records"))
    if len(ok) < a.min_runs:
        out["verdict"] = "UNADJUDICABLE"
        out["reason"] = (f"only {len(ok)} of {len(R)} null realisations have >= {a.min_quantiles} "
                         f"populated quantiles (>= {a.min_hits} same-family pairs) for BOTH "
                         f"configurations; >= {a.min_runs} are required before the ratio may be "
                         f"adjudicated. The denominator's tail is not estimable on this candidate set.")
        out["delta_old_admissible"] = [float(x) for x in ok.delta_old]
        out["delta_dense_all"] = [float(x) for x in R.delta_dense.dropna()]
    else:
        mo, md = float(ok.delta_old.mean()), float(ok.delta_dense.mean())
        ratio = md / mo if mo > 0 else np.nan
        rng = np.random.default_rng(7)
        boot = [float(ok.delta_dense.iloc[k].mean() / max(ok.delta_old.iloc[k].mean(), 1e-9))
                for k in (rng.integers(0, len(ok), len(ok)) for _ in range(4000))]
        lo = float(np.percentile(boot, 2.5))
        out.update(delta_old_mean=mo, delta_dense_mean=md,
                   delta_old_range=[float(ok.delta_old.min()), float(ok.delta_old.max())],
                   delta_dense_range=[float(ok.delta_dense.min()), float(ok.delta_dense.max())],
                   ratio_mean=ratio, ratio_boot_lb_2p5=lo,
                   ratio_per_run=[float(x) for x in ok.ratio],
                   predicted_delta_full_sigma=DELTA_OLD_P3 * ratio,
                   predicted_delta_full_lb=DELTA_OLD_P3 * lo,
                   rule=(f"licence requires 0.38 x ratio >= {BAR_SIGMA} sigma for BOTH the mean and "
                         "the 2.5 % bootstrap lower bound over null realisations"),
                   verdict=("STEP6_LICENSED"
                            if (DELTA_OLD_P3 * ratio >= BAR_SIGMA
                                and DELTA_OLD_P3 * lo >= BAR_SIGMA)
                            else "STEP6_NOT_LICENSED"))
    json.dump(out, open(f"{S.P3D}/step3_licence_{a.label}.json", "w"), indent=1, default=str)
    print(R.to_string(index=False))
    print()
    for k in ("n_admissible", "delta_old_admissible", "delta_dense_all", "delta_old_mean", "delta_dense_mean", "ratio_mean",
              "ratio_boot_lb_2p5", "predicted_delta_full_sigma", "predicted_delta_full_lb",
              "verdict", "reason"):
        if k in out:
            v = out[k]
            print(f"  {k:28s} {v:.3f}" if isinstance(v, float) else f"  {k:28s} {v}")
    print(f"  {'bar':28s} 0.38 x ratio >= {BAR_SIGMA} sigma (ratio >= {BAR_SIGMA/DELTA_OLD_P3:.1f}), "
          f"mean AND bootstrap LB, over >= {a.min_runs} admissible null realisations")


if __name__ == "__main__":
    main()
