#!/usr/bin/env python
"""VLFE scan pilot 2 (PLAN_SCAN.md S7.1): reference-station QC and the GNW quality check.
- per station-month: Z median of daily median 400-s band RMS; month excluded if > 3 x the station's 2011 median of daily
  medians; finite fraction and glitch-mask seconds per month;
- amplitude vs CN.PGC on 2011 events M >= 6.8 at D >= 1,500 km: peak |Z| (0.02-0.05 Hz, bp files) in the surface-wave group
  window [t + D/4.5 - 300, t + D/2.8 + 900] s, both >= 95 % finite there; ratio ref/PGC; station used only if the median
  ratio is in [0.5, 2] (no event -> used, 'unchecked').
Writes tables/scan2_refqc_month.csv, tables/scan2_refqc_events.csv, tables/scan2_refqc.csv (per station: use, reason)."""
import os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S
import scan_refs as R

BPR = f'{S.SD}/bp2_ref'
BPN = f'{S.SD}/bp_{S.YEAR}g'


def bp_window(netsta, a, b, pad=3600):
    """Band-passed (scan recipe, scan_common.bandpass_runs) Z between year-seconds a and b, computed from the v3 day files
    WITHOUT the S4 glitch mask (the 200-sigma mask removes large teleseismic surface waves, which this check needs)."""
    lo, hi = max(a - pad, 0), min(b + pad, S.NDAY * 86400)
    x = np.full(hi - lo, np.nan)
    for d in range(lo // 86400, (hi - 1) // 86400 + 1):
        day = S.Y0 + pd.Timedelta(days=d)
        f = f'{S.V3}/{netsta}/{day.year}/{day.dayofyear:03d}.npz'
        if os.path.isfile(f):
            z = np.load(f)['d'][2].astype(np.float64)
            s0 = d * 86400
            i0, i1 = max(lo, s0), min(hi, s0 + 86400)
            x[i0 - lo:i1 - lo] = z[i0 - s0:i1 - s0]
    return S.bandpass_runs(x)[a - lo:b - lo]


def main():
    T = R.load_table()
    T = T[T.bp_exists]
    stas = list(T.netsta)
    y0 = S.Y0.value / 1e9
    E = pd.read_csv(S.tab('scan2_eq_events.csv'))
    E = E[(E.mag >= 6.8) & (E.dist_km > 1500) & (E.t >= y0) & (E.t < y0 + S.NDAY * 86400)].sort_values('mag', ascending=False)
    E = E.loc[~(E.t.round(-2).duplicated())].sort_values('t')
    pgc = np.load(f'{BPN}/CN.PGC.npy', mmap_mode='r')
    mrows, erows, srows = [], [], []
    months = pd.date_range(f'{S.YEAR}-01-01', periods=S.NDAY, freq='D').month.values
    for s in stas:
        X = np.load(f'{BPR}/{s}.npy', mmap_mode='r')
        DS = np.load(f'{BPR}/{s}_daystats.npz')
        med = DS['med400'][:, 2] * 1e9; gl = DS['glitch_s'][:, 2]
        fin = np.isfinite(np.asarray(X[2])).reshape(S.NDAY, 86400).mean(axis=1)
        smed = np.nanmedian(med) if np.isfinite(med).any() else np.nan
        for m in range(1, 13):
            k = months == m
            mm = np.nanmedian(med[k]) if np.isfinite(med[k]).any() else np.nan
            mrows.append(dict(sta=s, month=m, z_med400_nm=mm, sta_2011_median_nm=smed, finite_frac=float(fin[k].mean()),
                              glitch_s=float(np.nansum(gl[k])), exclude=bool(np.isfinite(mm) and mm > 3 * smed),
                              live=bool(np.isfinite(mm))))
        rat = []
        for e in E.itertuples():
            a = int(e.t - y0 + e.dist_km / 4.5 - 300); b = int(e.t - y0 + e.dist_km / 2.8 + 900)
            if a < 0 or b > S.NDAY * 86400:
                continue
            xr = bp_window(s, a, b); xp = bp_window('CN.PGC', a, b)
            if np.isfinite(xr).mean() < 0.95 or np.isfinite(xp).mean() < 0.95:
                continue
            r = np.nanmax(np.abs(xr)) / np.nanmax(np.abs(xp))
            rat.append(r)
            erows.append(dict(sta=s, event_time=e.time, mag=e.mag, dist_km=e.dist_km, peak_ref_nm=np.nanmax(np.abs(xr)) * 1e9,
                              peak_pgc_nm=np.nanmax(np.abs(xp)) * 1e9, ratio=r))
        medr = float(np.median(rat)) if rat else np.nan
        use = (not rat) or (0.5 <= medr <= 2)
        srows.append(dict(sta=s, n_events=len(rat), median_ratio_vs_PGC=medr, ratio_min=min(rat) if rat else np.nan,
                          ratio_max=max(rat) if rat else np.nan, use=bool(use),
                          reason=('unchecked (no event)' if not rat else ('ratio in [0.5, 2]' if use else 'ratio outside [0.5, 2]')),
                          finite_frac_2011=float(fin.mean()), days_with_data=int((fin > 0).sum())))
    M = pd.DataFrame(mrows); M.to_csv(S.tab('scan2_refqc_month.csv'), index=False)
    pd.DataFrame(erows).to_csv(S.tab('scan2_refqc_events.csv'), index=False)
    Q = pd.DataFrame(srows); Q.to_csv(S.tab('scan2_refqc.csv'), index=False)
    print(Q.round(3).to_string())
    print(M[M.sta == 'UW.GNW'].round(2).to_string())


def usable_month(Q=None, M=None):
    """dict sta -> bool[12] (month usable) for QC-passed stations."""
    Q = pd.read_csv(S.tab('scan2_refqc.csv')) if Q is None else Q
    M = pd.read_csv(S.tab('scan2_refqc_month.csv')) if M is None else M
    out = {}
    for r in Q.itertuples():
        mm = M[M.sta == r.sta].sort_values('month')
        out[r.sta] = (~mm.exclude.values.astype(bool)) & bool(r.use)
    return out


if __name__ == '__main__':
    main()
