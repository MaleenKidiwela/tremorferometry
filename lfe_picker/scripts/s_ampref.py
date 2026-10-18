#!/usr/bin/env python
"""P3D amplitude reference for the new stations (the `rel_amp` denominator of o_run_month).

The frozen table `e5r_amp_ref_2012.csv` is the median sliding 3-C window energy over the 2012 quiet
days; the new stations have no 2012.  The same quantity is computed here over the P3D quiet months,
with the frozen window lengths and 1-s stride.  It enters NOTHING but the diagnostic `rel_amp`
column -- the detection statistic never sees it -- and is written to results/p3d/amp_ref_p3d.csv.
"""
from __future__ import annotations
import os, sys, argparse
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s_common as S                                                       # noqa: E402

LENS = (8.0, 12.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--months", default=",".join(S.QUIET_MONTHS))
    ap.add_argument("--max-days", type=int, default=6)
    a = ap.parse_args()
    ST = S.station_table()
    stas = [s for s in ST.index if ST.role[s] in ("new",) or s == "LRIV"]
    days = set()
    for m in a.months.split(","):
        for t in pd.date_range(f"{m}-01", pd.Timestamp(f"{m}-01") + pd.offsets.MonthEnd(0)):
            days.add((int(t.year), int(t.dayofyear)))
    idx = S.cache_index(years=sorted({y for y, _ in days}), stas=stas)
    rows = []
    for s in stas:
        got = [(y, d) for (ss, y, d) in sorted(idx) if ss == s and (y, d) in days][:a.max_days]
        if not got:
            print(f"  [ampref] {s}: no quiet-month cache")
            continue
        for L in LENS:
            n = int(round(L * S.FS))
            step = int(S.FS)
            vals = []
            for (y, d) in got:
                v = S.load_day_zne(s, y, d, idx)
                if v is None:
                    continue
                x = np.asarray(v[1], np.float64)
                if x.shape[0] != 3:
                    continue
                e = (x ** 2).sum(0)
                cs = np.concatenate([[0.0], np.cumsum(e)])
                st = np.arange(0, len(e) - n, step)
                w = cs[st + n] - cs[st]
                w = w[w > 0]
                if len(w):
                    vals.append(np.median(w))
            if vals:
                rows.append(dict(sta=s, L=L, e3_ref=float(np.median(vals)), n_days=len(got)))
                print(f"  [ampref] {s} L={L}: {rows[-1]['e3_ref']:.4g} ({len(got)} quiet days)")
    pd.DataFrame(rows).to_csv(f"{S.P3D}/amp_ref_p3d.csv", index=False)
    print(f"[ampref] -> {S.P3D}/amp_ref_p3d.csv")


if __name__ == "__main__":
    main()
