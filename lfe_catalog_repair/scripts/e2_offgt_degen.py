#!/usr/bin/env python
"""E2 step 4b (added 2026-09-10, resumed run) — the degenerate/zero-amplitude-window QC hole for
the LFEDetect picks at the off-ground-truth stations.

D3 found that 58 % of WISH's and 27 % of YACT's raw matched-filter re-detections sat in
zero-amplitude windows (NOTES, D3 stage 3).  In THESE day files the mechanism is identifiable: the
raw miniSEED contains no exact-zero samples at all (checked on WISH 2017 days: 0.000 zero
fraction), so a zero-amplitude window can only come from a GAP that the reader zero-fills.  The
measurement is therefore exact and needs no waveform pass: obspy `headonly=True` gives the covered
intervals per component per day, and a pick is classified by how much of its 60 s LFEDetect input
window is actually covered on ALL THREE components.

    gap_window      any second of the 60 s window missing on >= 1 component
    half_gap        more than half of the window missing on >= 1 component
    empty_window    the whole window missing on >= 1 component (a literally zero-amplitude input)

The same fractions are computed over every 60 s window of every processed day (the pick-independent
baseline), so enrichment = pick fraction / window fraction is visible.  A --amp-days sample also
re-reads that many days in full and confirms the header-based classification against the actual
sample values (max|x| == 0).
Writes results/e2_offgt_degen.csv (station x model x threshold) and e2_offgt_degen_days.csv.
"""
from __future__ import annotations

import argparse
import glob
import os

import numpy as np
import obspy
import pandas as pd

RES = "lfe_catalog_repair/results"
TRACES = "lfe_catalog_repair/traces"
STAS = {"WISH": "UW.WISH", "YACT": "UW.YACT", "B943": "PB.B943"}
THRS = (0.1, 0.3, 0.5)
WIN = 60.0


def coverage_mask(path):
    """per-second coverage of a day file: True where ALL components have data (86400 booleans)."""
    try:
        st = obspy.read(path, headonly=True)
    except Exception:
        return None
    day0 = obspy.UTCDateTime(st[0].stats.starttime.date).timestamp
    per = {}
    for tr in st:
        ch = tr.stats.channel[-1]
        m = per.setdefault(ch, np.zeros(86400, bool))
        a = int(np.floor(tr.stats.starttime.timestamp - day0))
        b = int(np.ceil(tr.stats.endtime.timestamp - day0))
        m[max(a, 0):min(b, 86400)] = True
    if len(per) < 3:
        return np.zeros(86400, bool) if not per else np.logical_and.reduce(
            list(per.values()) + [np.zeros(86400, bool)])
    return np.logical_and.reduce(list(per.values()))


def classify(cov, sec0, sec1):
    s0 = max(int(np.floor(sec0)), 0); s1 = min(int(np.ceil(sec1)), 86400)
    if s1 <= s0:
        return 1.0
    return 1.0 - float(cov[s0:s1].mean())      # missing fraction


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", default="WISH,YACT,B943")
    ap.add_argument("--amp-days", type=int, default=8,
                    help="days to re-read in full to confirm the header classification")
    args = ap.parse_args()
    fr = [pd.read_parquet(p) for p in sorted(glob.glob(f"{RES}/e2_picks_offgt_*.parquet"))]
    pk = pd.concat(fr, ignore_index=True)
    pk["day"] = pd.to_datetime(pk.day)
    pk["sec"] = (pd.to_datetime(pk.time).astype("datetime64[ns]").astype("int64") / 1e9) % 86400
    rows, drows = [], []
    rng = np.random.default_rng(0)
    for sta in args.stations.split(","):
        sub = pk[pk.sta == sta]
        if len(sub) == 0:
            print(f"[skip] {sta}: no picks on disk")
            continue
        days = np.sort(sub.day.unique())
        amp_days = set(pd.Timestamp(d) for d in
                       rng.choice(days, size=min(args.amp_days, len(days)), replace=False))
        acc = {}
        for d in days:
            dd = pd.Timestamp(d)
            p = f"{TRACES}/{STAS[sta]}/{dd.year}/{dd.dayofyear:03d}.mseed"
            if not os.path.exists(p):
                continue
            cov = coverage_mask(p)
            if cov is None:
                continue
            nw = int(86400 // WIN)
            wmiss = np.array([classify(cov, k * WIN, (k + 1) * WIN) for k in range(nw)])
            dayp = sub[sub.day == dd]
            pmiss = np.array([classify(cov, s - WIN / 2, s + WIN / 2) for s in dayp.sec.to_numpy()])
            drows.append(dict(station=sta, day=dd.date(), n_picks=len(dayp),
                              day_coverage=float(cov.mean()),
                              win_gap=float((wmiss > 0).mean()),
                              win_half_gap=float((wmiss > 0.5).mean()),
                              win_empty=float((wmiss >= 1.0).mean()),
                              pick_gap=float((pmiss > 0).mean()) if len(pmiss) else np.nan,
                              pick_half_gap=float((pmiss > 0.5).mean()) if len(pmiss) else np.nan,
                              pick_empty=float((pmiss >= 1.0).mean()) if len(pmiss) else np.nan))
            for (mo, ph), g in dayp.groupby(["model", "phase"]):
                if ph != "S":
                    continue
                i = dayp.index.get_indexer(g.index)
                acc.setdefault(mo, []).append(
                    pd.DataFrame(dict(prob=g.prob.values, thr_first=g.thr_first.values,
                                      miss=pmiss[i])))
            if dd in amp_days:
                st = obspy.read(p)
                st.merge(method=1, fill_value=0)
                z = min(float((np.abs(tr.data) == 0).mean()) for tr in st)
                drows[-1]["amp_zero_frac_min_comp"] = z
        dsum = pd.DataFrame([r for r in drows if r["station"] == sta])
        for mo, lst in sorted(acc.items()):
            arr = pd.concat(lst, ignore_index=True)
            for thr in THRS:
                q = arr[(arr.thr_first <= thr) & (arr.prob >= thr)]
                if not len(q):
                    continue
                rows.append(dict(station=sta, model=mo, threshold=thr, n_days=len(dsum),
                                 n_picks=len(q),
                                 pick_gap=float((q.miss > 0).mean()),
                                 pick_half_gap=float((q.miss > 0.5).mean()),
                                 pick_empty=float((q.miss >= 1.0).mean()),
                                 win_gap=dsum.win_gap.mean(),
                                 win_half_gap=dsum.win_half_gap.mean(),
                                 win_empty=dsum.win_empty.mean(),
                                 enrich_gap=float((q.miss > 0).mean())
                                 / max(dsum.win_gap.mean(), 1e-9),
                                 day_coverage_mean=dsum.day_coverage.mean(),
                                 days_below_99pct=int((dsum.day_coverage < 0.99).sum())))
        print(f"[{sta}] {len(dsum)} days, mean coverage {dsum.day_coverage.mean():.4f}, "
              f"{int((dsum.day_coverage < 0.99).sum())} days below 99 %", flush=True)
    pd.DataFrame(drows).to_csv(f"{RES}/e2_offgt_degen_days.csv", index=False)
    out = pd.DataFrame(rows)
    out.to_csv(f"{RES}/e2_offgt_degen.csv", index=False)
    print(out.round(5).to_string(index=False))
    print(f"[write] {RES}/e2_offgt_degen.csv, {RES}/e2_offgt_degen_days.csv")


if __name__ == "__main__":
    main()
