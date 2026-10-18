#!/usr/bin/env python
"""D3 stage 3 addendum — the STORED catalogue as a baseline, on exactly the pilot window.

`d3_validate.py` compares the new MAD rules against `cc_stack >= 0.80`, which isolates the
threshold rule from the 3-C stacking. The third leg of the comparison is what is actually ON DISK:
the densify detections at a fixed **Z-only** cc >= 0.80 with a 6 s gap (`data/mf_<sta>p90f40_<yr>.csv`,
read-only). This script runs the same ETS/quiet and tremor-correlation statistics on those rows,
restricted to the pilot window, so the three rules are measured on the same days and cells.

Note the stored files cover only the 300 coverage-selected families per station, not all of them.
"""
from __future__ import annotations

import argparse
import glob
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from d3_validate import (cell_of, load_tremor, cell_tremor_daily, ets_quiet_mask)   # noqa: E402
from scipy.stats import spearmanr                                                    # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sta", required=True, help="station tag used in the mf filenames, e.g. wish")
    ap.add_argument("--label", required=True, help="row label, e.g. WISH")
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--glob", default="")
    ap.add_argument("--out", default="lfe_catalog_repair/results/d3_storedrule_baseline.csv")
    a = ap.parse_args()
    pat = a.glob or f"data/mf_{a.sta}p90f40_*.csv"
    files = sorted(glob.glob(pat))
    d0, d1 = pd.Timestamp(a.start), pd.Timestamp(a.end)
    keep = [f for f in files if d0.year <= int(f.rsplit("_", 1)[1][:4]) <= d1.year]
    print(f"[stored {a.label}] {len(keep)} year files of {len(files)}: "
          f"{[Path(f).name for f in keep]}", flush=True)
    if not keep:
        print("  no stored files -> skip")
        return
    parts = []
    for f in keep:
        for ch in pd.read_csv(f, usecols=["template", "time", "cc"], chunksize=4_000_000,
                              parse_dates=["time"]):
            ch = ch[(ch.time >= d0) & (ch.time < d1)]
            if len(ch):
                parts.append(ch)
    det = pd.concat(parts, ignore_index=True)
    del parts
    det["day"] = det.time.dt.floor("D")
    print(f"  {len(det):,} stored detections in {a.start}..{a.end}, "
          f"{det.template.nunique()} templates", flush=True)
    all_days = pd.DatetimeIndex(sorted(det.day.unique()))
    tremor = load_tremor()
    cells, rows = {}, []
    day_index = {d: i for i, d in enumerate(all_days)}
    for fam, g in det.groupby("template", sort=False):
        lat, lon = cell_of(str(fam))
        ck = (round(lat, 3), round(lon, 3))
        if ck not in cells:
            cells[ck] = ets_quiet_mask(cell_tremor_daily(tremor, lat, lon), all_days)
        ets_m, qui_m, roll = cells[ck]
        cnt = np.zeros(len(all_days))
        for d, n in g.day.value_counts().items():
            j = day_index.get(d)
            if j is not None:
                cnt[j] = n
        re_, rq = cnt[ets_m].mean(), cnt[qui_m].mean()
        rho, pv = (np.nan, np.nan)
        if np.std(cnt) > 0 and np.std(roll) > 0:
            rho, pv = spearmanr(cnt, roll)
        rows.append(dict(template=str(fam), lat=lat, lon=lon, n_det=len(g),
                         rate_ets=re_, rate_quiet=rq,
                         ratio=(re_ / rq) if rq > 0 else np.nan, rho=rho, p=pv,
                         cc_med=float(g.cc.median())))
    T = pd.DataFrame(rows)
    summ = dict(station=a.label, rule="stored_Zonly_cc0.80",
                n_det=int(T.n_det.sum()), n_templates=len(T),
                det_per_tmpl_day=float(T.n_det.sum() / max(len(T) * len(all_days), 1)),
                ratio_median=float(T.ratio.median()),
                frac_ratio_ge2=float((T.ratio >= 2).mean()),
                frac_ratio_ge5=float((T.ratio >= 5).mean()),
                rho_median=float(T.rho.median()),
                frac_rho_pos=float((T.rho > 0).mean()),
                frac_rho_sig=float((T.p < 0.05).mean()),
                cc_median=float(T.cc_med.median()))
    out = Path(a.out)
    S = pd.DataFrame([summ])
    if out.exists():
        prev = pd.read_csv(out)
        prev = prev[prev.station != a.label]
        S = pd.concat([prev, S], ignore_index=True)
    S.to_csv(out, index=False)
    T.to_csv(str(out).replace(".csv", f"_{a.label}_templates.csv"), index=False)
    print(pd.DataFrame([summ]).T.to_string())


if __name__ == "__main__":
    main()
