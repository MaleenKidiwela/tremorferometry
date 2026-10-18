#!/usr/bin/env python
"""D3 stage 3b — network coincidence between two stations for families in the SAME source cell.

For every source cell carried by families at BOTH stations:
  1. estimate the inter-station lag from the cross-correlation of the two stations' daily-binned
     detection-time series during ETS (0.5 s bins, +-60 s scan; computed as the pair-difference
     histogram, which is the same thing but runnable) -- the same construction D1 used;
  2. count detections at station A that have a same-cell detection at station B within +-3 s of
     (t_A + lag);
  3. repeat with the B series shifted by +600 s = the chance-level control.
Reported per rule (old fixed cc>=0.80 vs new MAD k=8/9/10) and separately for ETS and quiet days.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

sys_path_hint = "run from the repo root"
PNSN = "catalogs/pnsn_tremor_cascadia_full.csv"
KS = (8, 9, 10)
SHIFT_CONTROL_S = 600.0
COINC_S = 3.0
LAG_SCAN_S = 60.0
BIN_S = 0.5


def cell_of(fam):
    return fam.split("__")[0]


def cell_latlon(cell):
    a = cell.split("_")
    return float(a[0]), float(a[1])


def ets_days(tremor, lat, lon, all_days, radius_km=40.0):
    d = np.hypot((tremor.lat.values - lat) * 111.2,
                 (tremor.lon.values - lon) * 111.2 * np.cos(np.radians(lat)))
    s = tremor.loc[d <= radius_km].groupby("day").size()
    c = s.reindex(all_days, fill_value=0).astype(float)
    roll = c.rolling(30, min_periods=1).sum()
    p90, p25 = np.percentile(roll.values, 90), np.percentile(roll.values, 25)
    return set(all_days[roll.values > p90]), set(all_days[roll.values < p25])


def xcorr_lag(ta, tb, max_lag=LAG_SCAN_S, binsec=BIN_S):
    """Lag (s) maximising the cross-correlation of the two binned detection-time series.

    Computed as the histogram of pair time differences within +-max_lag, which IS the
    cross-correlation of the delta-binned series but costs O(n log n + n_pairs) instead of
    O(n_bins * n_lags) -- a 2-year span at 1 s bins is 6e7 bins, so the direct product form is
    not runnable.
    """
    if len(ta) < 5 or len(tb) < 5:
        return 0.0, 0.0
    a = np.sort(np.asarray(ta, float))
    b = np.sort(np.asarray(tb, float))
    i0 = np.searchsorted(b, a - max_lag)
    i1 = np.searchsorted(b, a + max_lag)
    n = int((i1 - i0).sum())
    if n == 0:
        return 0.0, 0.0
    if n > 40_000_000:                      # subsample A if the pair count explodes
        keep = np.random.default_rng(0).choice(len(a), max(1000, int(len(a) * 4e7 / n)),
                                              replace=False)
        a = np.sort(a[keep])
        i0 = np.searchsorted(b, a - max_lag)
        i1 = np.searchsorted(b, a + max_lag)
    diffs = np.concatenate([b[x:y] - t for t, x, y in zip(a, i0, i1) if y > x])
    if diffs.size == 0:
        return 0.0, 0.0
    bins = np.arange(-max_lag, max_lag + binsec, binsec)
    h, _ = np.histogram(diffs, bins=bins)
    if h.max() < 3:
        return 0.0, 0.0
    lag = float(bins[int(h.argmax())] + binsec / 2)
    # "xcorr" quality = peak height over the mean of the rest of the scan (1 = flat/chance)
    base = float((h.sum() - h.max()) / max(len(h) - 1, 1))
    return lag, float(h.max() / base) if base > 0 else 0.0


def coincidence(ta, tb, lag, tol=COINC_S):
    if len(ta) == 0 or len(tb) == 0:
        return 0, 0.0
    b = np.sort(tb)
    x = np.sort(ta) + lag
    i = np.searchsorted(b, x)
    d = np.full(len(x), np.inf)
    ok = i < len(b)
    d[ok] = np.abs(b[i[ok]] - x[ok])
    ok2 = i > 0
    d[ok2] = np.minimum(d[ok2], np.abs(b[i[ok2] - 1] - x[ok2]))
    n = int((d <= tol).sum())
    return n, n / len(x)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--a", required=True, help="station A code")
    ap.add_argument("--b", required=True, help="station B code")
    ap.add_argument("--resdir", default="lfe_catalog_repair/results")
    ap.add_argument("--out", default="")
    ap.add_argument("--min-det", type=int, default=20)
    a = ap.parse_args()
    R = Path(a.resdir)
    A = pd.read_parquet(R / f"d3_detections_{a.a}.parquet")
    B = pd.read_parquet(R / f"d3_detections_{a.b}.parquet")
    for D in (A, B):
        D["cell"] = D["template"].astype(str).map(cell_of)
        D["day"] = pd.to_datetime(D["t"], unit="s").dt.floor("D")
    tremor = pd.read_csv(PNSN, usecols=["time", "lat", "lon"], parse_dates=["time"])
    tremor["day"] = tremor.time.dt.floor("D")
    all_days = pd.DatetimeIndex(sorted(set(A.day.unique()) | set(B.day.unique())))
    cells = sorted(set(A.cell) & set(B.cell))
    Ag = {k: v for k, v in A.groupby("cell", sort=False)}
    Bg = {k: v for k, v in B.groupby("cell", sort=False)}
    print(f"[{a.a} x {a.b}] {len(cells)} shared source cells; "
          f"{A.template.nunique()} vs {B.template.nunique()} templates", flush=True)
    rows = []
    for cell in cells:
        lat, lon = cell_latlon(cell)
        ets, qui = ets_days(tremor, lat, lon, all_days)
        sa_all, sb_all = Ag.get(cell), Bg.get(cell)
        if sa_all is None or sb_all is None:
            continue
        for rule in ["old"] + [f"new{k}" for k in KS]:
            sa = sa_all[sa_all[rule]] if rule in sa_all else sa_all
            sb = sb_all[sb_all[rule]] if rule in sb_all else sb_all
            ta_e = sa.loc[sa.day.isin(ets), "t"].values
            tb_e = sb.loc[sb.day.isin(ets), "t"].values
            ta_q = sa.loc[sa.day.isin(qui), "t"].values
            tb_q = sb.loc[sb.day.isin(qui), "t"].values
            if len(ta_e) < a.min_det or len(tb_e) < a.min_det:
                continue
            lag, cc = xcorr_lag(ta_e, tb_e)
            n_e, f_e = coincidence(ta_e, tb_e, lag)
            n_q, f_q = coincidence(ta_q, tb_q, lag)
            n_s, f_s = coincidence(ta_e, tb_e + SHIFT_CONTROL_S, lag)
            nd_e = max(len(set(sa.day) & ets), 1)
            nd_q = max(len(set(sa.day) & qui), 1)
            rows.append(dict(cell=cell, lat=lat, lon=lon, rule=rule, lag_s=lag, xcorr=cc,
                             n_a_ets=len(ta_e), n_b_ets=len(tb_e),
                             n_a_quiet=len(ta_q), n_b_quiet=len(tb_q),
                             coinc_ets=n_e, coinc_quiet=n_q, coinc_shift=n_s,
                             frac_ets=f_e, frac_quiet=f_q, frac_shift=f_s,
                             rate_ets=n_e / nd_e, rate_quiet=n_q / nd_q, rate_shift=n_s / nd_e))
    out = pd.DataFrame(rows)
    f = Path(a.out) if a.out else R / f"d3_coincidence_{a.a}_{a.b}.csv"
    out.to_csv(f, index=False)
    if len(out):
        for rule, g in out.groupby("rule"):
            print(f"  {rule:6s} cells={len(g):4d} lag med {g.lag_s.median():+6.1f} s  "
                  f"frac_ets {g.frac_ets.mean():.4f}  frac_quiet {g.frac_quiet.mean():.4f}  "
                  f"frac_shift {g.frac_shift.mean():.4f}  "
                  f"ETS/shift {g.frac_ets.mean()/max(g.frac_shift.mean(),1e-9):.1f}  "
                  f"ETS/quiet {g.frac_ets.mean()/max(g.frac_quiet.mean(),1e-9):.2f}")
    print(f"  -> {f}")


if __name__ == "__main__":
    main()
