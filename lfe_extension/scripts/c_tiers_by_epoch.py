#!/usr/bin/env python
"""Track C — re-derive the TWO TIERS PER NETWORK EPOCH, by the frozen false-alarm-only rule.

WHY.  T_core and T_ext are false-alarm thresholds, and the false-alarm rate at a fixed MAD is a
property of the NETWORK: a smaller network is noisier, so the same MAD buys a worse catalog (Track D
measured ~3x the quiet-day rate at the frozen T after SNB was lost in 2021-03).  Holding 8.25 / 8.75
fixed across the record would therefore mean a catalog whose purity silently changes with the station
list — exactly what the FA-only rule exists to prevent.

THE RULE (unchanged from Track A step 1; only the data it is measured on changes):
    T_core = smallest T on a 0.25-MAD grid, at or above the 8.0 recording floor, with
             quiet-day FA <= 0.1 per template-day
    T_ext  = smallest T with quiet-day FA <= 1.0 per template-day AND post-SSE FA <= 1.5
Measured on THAT epoch's own days, corrected by nothing (the raw rate, as in Track A):
    quiet family-day   observable, no PNSN tremor in the family box that day, and >= 3 days from any
                       tremor day  (the `quiet_flag` already stored in family_days_*.parquet)
    post-SSE family-day  observable, 1..10 days after the last tremor day of a tremor episode
                       (an episode ends when >= 3 tremor-free days follow) - the production analogue
                       of Track A's 2012 doy 270-279 window
    FA(T)              detections with net_mad >= T AND hf_pass, on those family-days, divided by the
                       number of those family-days.  The HF veto is applied because the frozen point
                       was derived with it; cross-family assignment is NOT, because the rate is per
                       template per day.
net_mad, the templates and the HF cut are NOT touched by any of this.

  python scripts/c_tiers_by_epoch.py            # derive and write results/trackC/tiers_by_epoch.csv
  python scripts/c_tiers_by_epoch.py --report   # derive and print, write nothing
"""
from __future__ import annotations
import os, sys, glob, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C

POST_LO, POST_HI = 1, 10          # days after an episode ends
EPISODE_GAP = 3                   # tremor-free days that end an episode


def post_sse_dates(trem_days, all_days):
    """dates that are POST_LO..POST_HI days after the last tremor day of an episode."""
    td = sorted(pd.to_datetime(sorted(trem_days)))
    if not td:
        return set()
    ends = [td[i] for i in range(len(td))
            if i == len(td) - 1 or (td[i + 1] - td[i]).days >= EPISODE_GAP]
    out = set()
    for e in ends:
        for k in range(POST_LO, POST_HI + 1):
            out.add((e + pd.Timedelta(days=k)).strftime("%Y-%m-%d"))
    return out - {t.strftime("%Y-%m-%d") for t in td}


def load(out_dir):
    FDs, DETs = [], []
    for f in sorted(glob.glob(f"{out_dir}/family_days_*.parquet")):
        ym = os.path.basename(f).rsplit("_", 1)[1][:-8]
        d = pd.read_parquet(f); d["ym"] = ym
        FDs.append(d)
        p = f"{out_dir}/detections_{ym}.parquet"
        if os.path.exists(p):
            e = pd.read_parquet(p); e["ym"] = ym
            DETs.append(e[["fam", "date", "net_mad", "hf_pass", "ym"]])
    if not FDs:
        raise SystemExit("[tiers] no processed months")
    return pd.concat(FDs, ignore_index=True), (pd.concat(DETs, ignore_index=True)
                                               if DETs else pd.DataFrame())


def fa_curve(FDsel, DET, grid):
    """FA per template-day on the given family-day set, at each threshold."""
    n = len(FDsel)
    if n == 0 or not len(DET):
        return np.full(len(grid), np.nan), 0
    key = set(zip(FDsel.fam.astype(int), FDsel.date))
    d = DET[[(f, dt) in key for f, dt in zip(DET.fam.astype(int), DET.date)]]
    d = d[d.hf_pass]
    m = d.net_mad.to_numpy()
    return np.array([float((m >= T).sum()) / n for T in grid]), n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default="")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--tmin", type=float, default=8.0)
    ap.add_argument("--tmax", type=float, default=15.0)
    ap.add_argument("--step", type=float, default=0.25)
    ap.add_argument("--min-quiet", type=int, default=200,
                    help="family-days needed before an epoch's tiers are derived rather than anchored")
    ap.add_argument("--min-post", type=int, default=100)
    ap.add_argument("--rederive", action="store_true",
                    help="recompute epochs that are already marked `derived` in the existing file "
                         "(default: a derived row is FROZEN, so the tiers never drift as months arrive)")
    a = ap.parse_args()
    OUT = a.out_dir or C.TRC
    fz = C.frozen()
    T_core0, T_ext0 = float(fz["T_core"]), float(fz["T_ext"])
    FD, DET = load(OUT)
    FD = FD[FD.observable].copy()
    FD["epoch"] = FD.ym.map(C.epoch_of_month)
    trem = set(FD[FD.n_tremor_day > 0].date.unique())
    post = post_sse_dates(trem, set(FD.date.unique()))
    grid = np.round(np.arange(a.tmin, a.tmax + 1e-9, a.step), 4)
    rows = []
    for name, m0, m1, why in C.NETWORK_EPOCHS:
        E = FD[FD.epoch == name]
        if not len(E):
            rows.append(dict(epoch=name, months=f"{m0}..{m1}", note=why, status="no data",
                             T_core=np.nan, T_ext=np.nan, n_quiet_family_days=0,
                             n_post_family_days=0))
            continue
        D = DET[DET.ym.map(C.epoch_of_month) == name] if len(DET) else DET
        q, nq = fa_curve(E[E.quiet_flag], D, grid)
        p_, npo = fa_curve(E[E.date.isin(post)], D, grid)
        tc = next((float(T) for T, v in zip(grid, q) if np.isfinite(v) and v <= 0.1), np.nan)
        te = next((float(T) for T, v, w in zip(grid, q, p_)
                   if np.isfinite(v) and v <= 1.0 and (not np.isfinite(w) or w <= 1.5)), np.nan)
        status = "derived"
        # E0 ANCHOR (coordinator): if E0's own quiet days reproduce the frozen 8.25 / 8.75 within
        # 0.25 MAD, the frozen pair STAYS - the operating point the whole method was frozen on is
        # not moved by a re-measurement that agrees with it.
        anchored = ""
        if name == "E0" and np.isfinite(tc) and np.isfinite(te) and nq >= a.min_quiet \
                and abs(tc - T_core0) <= 0.25 and abs(te - T_ext0) <= 0.25:
            anchored = f" (derived {tc:.2f}/{te:.2f}, within 0.25 -> frozen anchor kept)"
            tc, te = T_core0, T_ext0
        if nq < a.min_quiet or not np.isfinite(tc) or not np.isfinite(te):
            status = f"insufficient (quiet {nq} < {a.min_quiet})" if nq < a.min_quiet else "grid exhausted"
            tc, te = (T_core0, T_ext0)
        elif npo < a.min_post:
            status = f"derived, post-SSE thin ({npo} family-days): T_ext from the quiet bar alone"
        status += anchored
        gi = lambda T: int(np.argmin(np.abs(grid - T))) if np.isfinite(T) else 0
        rows.append(dict(epoch=name, months=f"{m0}..{m1}", note=why, status=status,
                         T_core=tc, T_ext=te,
                         n_quiet_family_days=nq, n_post_family_days=npo,
                         quiet_fa_at_T_core=float(q[gi(tc)]) if nq else np.nan,
                         quiet_fa_at_T_ext=float(q[gi(te)]) if nq else np.nan,
                         post_fa_at_T_ext=float(p_[gi(te)]) if npo else np.nan,
                         quiet_fa_at_frozen_T_core=float(q[gi(T_core0)]) if nq else np.nan,
                         quiet_fa_at_frozen_T_ext=float(q[gi(T_ext0)]) if nq else np.nan,
                         months_used=",".join(sorted(E.ym.unique()))))
    R = pd.DataFrame(rows)
    pd.set_option("display.width", 250)
    print(R[["epoch", "months", "status", "T_core", "T_ext", "n_quiet_family_days",
             "n_post_family_days", "quiet_fa_at_T_core", "quiet_fa_at_T_ext", "post_fa_at_T_ext",
             "quiet_fa_at_frozen_T_core", "quiet_fa_at_frozen_T_ext"]].round(4).to_string(index=False))
    e0 = R[R.epoch == "E0"]
    if len(e0) and e0.status.iloc[0] == "derived":
        dc = abs(float(e0.T_core.iloc[0]) - T_core0); de = abs(float(e0.T_ext.iloc[0]) - T_ext0)
        print(f"\n[anchor] E0 re-derived on its OWN quiet days: T_core {e0.T_core.iloc[0]:.2f} "
              f"(frozen {T_core0:.2f}, |d| {dc:.2f}), T_ext {e0.T_ext.iloc[0]:.2f} "
              f"(frozen {T_ext0:.2f}, |d| {de:.2f}) -> "
              f"{'REPRODUCES the frozen anchor within 0.25' if max(dc, de) <= 0.25 else 'DOES NOT reproduce it within 0.25'}")
    if not a.report:
        prev_p = f"{OUT}/tiers_by_epoch.csv"
        if os.path.exists(prev_p) and not a.rederive:
            prev = pd.read_csv(prev_p)
            keep = prev[prev.status.astype(str).str.startswith("derived")]
            if len(keep):
                R = pd.concat([R[~R.epoch.isin(set(keep.epoch))], keep], ignore_index=True)
                R = R.sort_values("epoch")
                print(f"[tiers] {len(keep)} epoch(s) already derived are FROZEN "
                      f"({', '.join(keep.epoch)}); pass --rederive to recompute them")
        R.to_csv(prev_p, index=False)
        print(f"\n[tiers] -> {OUT}/tiers_by_epoch.csv")


if __name__ == "__main__":
    main()
