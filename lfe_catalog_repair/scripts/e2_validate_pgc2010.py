#!/usr/bin/env python
"""E2 step 2 — single-station validation of the LFEDetect picks at PGC, 2010-05-01..2010-11-01.

Ground truth
  * Bostock et al. (2015) origin times + family hypocentres (Slow Earthquake DB).  The DB record
    does not state its time convention, so the offset to a true origin time is calibrated
    INDEPENDENTLY of our picks from the Lin catalogue (which publishes a grid-searched origin time)
    on co-detected events -- see e2_groundtruth.py.
    Predicted PGC S arrival, two variants:
      (bostock_stack) t + C + ttP_model(PGC) + (S-P) read off the SeisBench Bostock LFE stacks
                      for that family at PGC   <- primary; no S-velocity assumption
      (bostock_1d)    t + C + ttS_model(PGC)                        <- pure Crosson-1D
  * Lin et al. (2024) located events: OT + ttS_model(PGC).  Lin depths are unreliable (a third of
    the catalogue is pinned to his grid bottom), so a depth-fixed variant at 33 km is also scored.

Chance control.  At ~10^3 picks/day a +-3 s window covers a non-negligible fraction of the day, so
every recall/precision number is reported next to a null obtained by shifting the ground-truth
times by +-600 s (and by a random uniform shift), which is the recall/precision an unrelated pick
series of the same density would score.

Outputs `results/e2_validation_pgc2010.csv` and `results/e2_pgc2010_daily.csv`.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
PGC = dict(lat=48.6498, lon=-123.4521)
TOLS = (3.0, 5.0)
SWEEP = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)
NULL_SHIFTS = (-600.0, 600.0, 1800.0)
W0, W1 = "2010-05-01", "2010-11-01"


def epoch(s):
    return pd.to_datetime(s).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9


def match_frac(a, b, tol):
    """fraction of `a` with a `b` within +-tol, and fraction of `b` with an `a` within +-tol."""
    if len(a) == 0 or len(b) == 0:
        return 0.0, 0.0
    a = np.sort(a); b = np.sort(b)
    i = np.searchsorted(b, a)
    d = np.full(len(a), np.inf)
    for off in (-1, 0):
        j = np.clip(i + off, 0, len(b) - 1)
        d = np.minimum(d, np.abs(b[j] - a))
    fa = float((d <= tol).mean())
    i = np.searchsorted(a, b)
    d = np.full(len(b), np.inf)
    for off in (-1, 0):
        j = np.clip(i + off, 0, len(a) - 1)
        d = np.minimum(d, np.abs(a[j] - b))
    fb = float((d <= tol).mean())
    return fa, fb


def load_picks():
    fr = []
    for m in ("cascadia", "jcms"):
        p = f"{RES}/e2_picks_pgc2010_{m}.parquet"
        if os.path.exists(p):
            fr.append(pd.read_parquet(p))
    df = pd.concat(fr, ignore_index=True)
    df["t"] = epoch(df.time)
    return df


def gt_tables():
    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    lin = pd.read_parquet(f"{RES}/e2_gt_lin.parquet")
    lin = lin[lin.window == "pgc2010"].copy()
    cal = pd.read_csv(f"{RES}/e2_stack_calib.csv")
    dt = pd.read_csv(f"{RES}/e2_bostock_lin_dt.csv")
    C = float(np.median(dt.dt))
    # (S-P) at PGC per Bostock family, from the SeisBench stacks
    sp = (cal[cal.station_code == "PGC"].set_index("source_family_id").sp_obs_s)
    bo["sp_pgc"] = bo.fam_id.map(sp)
    bo["t_cat"] = epoch(bo.time)
    bo["pred_stack"] = bo.t_cat + C + bo.ttP_PGC + bo.sp_pgc
    bo["pred_1d"] = bo.t_cat + C + bo.ttS_PGC
    lin["t_cat"] = epoch(lin.OT)
    lin["pred_1d"] = lin.t_cat + lin.ttS_PGC
    # depth-fixed variant: recompute ttS at 33 km using the cached table
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import e2_groundtruth as gtm
    lin["pred_z33"] = lin.t_cat + gtm.tt_vec(lin.epi_PGC.values, np.full(len(lin), 33.0), "S")
    return bo, lin, C


def main():
    pk = load_picks()
    bo, lin, C = gt_tables()
    print(f"[calib] Bostock time convention offset C = {C:+.3f} s (from Lin OT on co-detections)")
    print(f"[gt] Bostock {len(bo)} events ({bo.pred_stack.notna().sum()} with a PGC stack S-P), "
          f"Lin {len(lin)} events in {W0}..{W1}")

    linq = lin[(lin.residual < 0.5) & (lin.N >= 4)]          # Lin's own "high quality" cut
    gts = {
        "bostock_stack": bo.pred_stack.dropna().to_numpy(),
        "bostock_1d": bo.pred_1d.dropna().to_numpy(),
        "bostock_near40": bo[bo.epi_PGC < 40].pred_stack.dropna().to_numpy(),
        "lin_1d": lin.pred_1d.dropna().to_numpy(),
        "lin_z33": lin.pred_z33.dropna().to_numpy(),
        "lin_hq": linq.pred_1d.dropna().to_numpy(),
        "lin_near40": lin[lin.epi_PGC < 40].pred_1d.dropna().to_numpy(),
    }
    gts["either"] = np.sort(np.concatenate([gts["bostock_stack"], gts["lin_1d"]]))
    print("[gt sizes] " + ", ".join(f"{k}={len(v)}" for k, v in gts.items()))

    # ---- pick-time offset per ground truth
    # BUGFIX 2026-09-10 (resumed run): the previous estimator was "median residual of the NEAREST S
    # pick within 10 s".  At 150-900 picks/day almost every ground-truth event has a pick within
    # 10 s BY CHANCE (16,026 of 18,114 Bostock events did), so that median is pulled to ~0 whatever
    # the true offset is -- it reported C = +0.09 s for the Bostock catalogue whose true offset to
    # our picks is +12.4 s.  Replaced by the PEAK of the residual histogram over +-40 s, built from
    # ALL picks in the window (0.25 s bins, background = median bin count > 5 s from the peak),
    # which is the standard cross-correlation estimator and is chance-immune.
    off, offdiag = {}, []
    S01 = np.sort(pk[(pk.model == "cascadia") & (pk.phase == "S") & (pk.thr_first <= 0.3)
                     & (pk.prob >= 0.3)].t.to_numpy())
    import e2_resid_hist as RH
    for k, g in gts.items():
        d = RH.pairs_within(np.sort(g), S01, 40.0)
        e = np.arange(-40.0, 40.001, 0.25)
        h, _ = np.histogram(d, bins=e)
        c = 0.5 * (e[1:] + e[:-1])
        i = int(np.argmax(h))
        bg = max(float(np.median(h[np.abs(c - c[i]) > 5.0])), 1e-9)
        off[k] = float(c[i])
        offdiag.append(dict(gt=k, peak_lag_s=c[i], peak_over_bg=h[i] / bg,
                            peak_counts=int(h[i]), bg_per_bin=bg, n_gt=len(g)))
        print(f"[offset] {k}: peak of (pick - predicted) at {off[k]:+.2f} s, "
              f"{h[i]/bg:.1f}x background ({h[i]} counts in a 0.25 s bin, bg {bg:.1f})")
    pd.DataFrame(offdiag).to_csv(f"{RES}/e2_pgc2010_offsets.csv", index=False)

    rows = []
    rng = np.random.default_rng(0)
    for model in sorted(pk.model.unique()):
        for phase in ("S", "P"):
            for thr in SWEEP:
                sub = pk[(pk.model == model) & (pk.phase == phase)
                         & (pk.thr_first <= thr) & (pk.prob >= thr)]
                t = np.sort(sub.t.to_numpy())
                for gk, g in gts.items():
                    gg = np.sort(g + off[gk])
                    for tol in TOLS:
                        rec, prec = match_frac(gg, t, tol)
                        nulls_r, nulls_p = [], []
                        for sh in NULL_SHIFTS:
                            r, p = match_frac(gg + sh, t, tol)
                            nulls_r.append(r); nulls_p.append(p)
                        r, p = match_frac(gg + rng.uniform(3000, 40000), t, tol)
                        nulls_r.append(r); nulls_p.append(p)
                        rows.append(dict(
                            model=model, phase=phase, threshold=thr, gt=gk, tol_s=tol,
                            n_gt=len(gg), n_picks=len(t), picks_per_day=len(t) / 184.0,
                            recall=rec, precision=prec,
                            recall_null=float(np.mean(nulls_r)),
                            precision_null=float(np.mean(nulls_p)),
                            recall_excess=rec - float(np.mean(nulls_r)),
                            offset_s=off[gk]))
    val = pd.DataFrame(rows)
    val.to_csv(f"{RES}/e2_validation_pgc2010.csv", index=False)
    print(f"[write] {RES}/e2_validation_pgc2010.csv  {len(val)} rows")

    show = val[(val.phase == "S") & (val.tol_s == 3.0)
               & (val["gt"].isin(["bostock_stack", "lin_1d"]))]  # BUGFIX 2026-09-10: column "gt" shadows DataFrame.gt
    print(show.pivot_table(index=["gt", "threshold"], columns="model",
                           values=["recall", "recall_null", "precision", "picks_per_day"]
                           ).round(3).to_string())

    # ---- daily counts vs PNSN tremor within 40 km of PGC
    tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time).astype("datetime64[ns]")
    d = np.hypot((tr.lat - PGC["lat"]) * 111.195,
                 (tr.lon - PGC["lon"]) * 111.195 * np.cos(np.radians(PGC["lat"])))
    tr = tr[(d < 40) & (tr.time >= W0) & (tr.time < W1)]
    trd = tr.groupby(tr.time.dt.date).size().rename("tremor")
    daily = []
    for model in sorted(pk.model.unique()):
        for thr in (0.1, 0.3, 0.5):
            s = pk[(pk.model == model) & (pk.phase == "S") & (pk.thr_first <= thr)
                   & (pk.prob >= thr)]
            c = s.groupby(s.day).size().rename("picks")
            df = pd.concat([c, trd], axis=1).fillna(0.0)
            df.index.name = "date"
            df = df.reset_index()
            df["model"] = model
            df["threshold"] = thr
            daily.append(df)
    daily = pd.concat(daily, ignore_index=True)
    daily.to_csv(f"{RES}/e2_pgc2010_daily.csv", index=False)
    from scipy.stats import spearmanr
    print("\n[daily] ETS / quiet ratio of daily S-pick counts (ETS = 30-d tremor > p90, "
          "quiet = < p25, D1 convention) and Spearman vs daily tremor")
    for (model, thr), g in daily.groupby(["model", "threshold"]):
        g = g.sort_values("date")
        roll = g.tremor.rolling(30, min_periods=10, center=True).sum()
        hi, lo = np.nanpercentile(roll, 90), np.nanpercentile(roll, 25)
        ets = g.picks[roll >= hi].mean()
        qui = g.picks[roll <= lo].mean()   # BUGFIX 2026-09-10: strict < gave an EMPTY quiet set
        #   (>25 % of days have zero tremor within 40 km -> p25 == min == 0, `roll < 0` never true)
        rho, pv = spearmanr(g.picks, g.tremor)
        print(f"  {model:9s} thr {thr}: ETS {ets:7.1f}/d  quiet {qui:7.1f}/d  "
              f"ratio {ets/max(qui,1e-9):5.2f}  spearman {rho:+.3f} (p={pv:.1e})")


if __name__ == "__main__":
    main()
