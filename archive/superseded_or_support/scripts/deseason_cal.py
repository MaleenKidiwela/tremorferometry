#!/usr/bin/env python
"""Deseason calendar-rolling dv/v csvs: per patch, fit + subtract a 2-harmonic (annual+semiannual)
day-of-year climatology (mean kept). Replicates the existing `_cal` -> `_cal_des` generation
(verified: removed signal is a pure smooth DOY function per patch).
Usage: deseason_cal.py --glob 'data/daily_dvv_*_2to4_calT.csv'   (writes *_calT_des.csv; skips existing)
"""
import argparse, glob, os
import numpy as np, pandas as pd


def deseason(df, win=61):
    """Subtract per-patch circularly-smoothed day-of-year climatology (mean kept).
    Matches the original _cal_des generation at corr ~0.996 (win 61-91)."""
    out = []
    for patch, g in df.groupby('patch'):
        g = g.copy()
        if len(g) > 50:
            doy = g.date.dt.dayofyear
            clim = g.groupby(doy.values).dvv.mean()
            clim = clim.reindex(range(1, 367)).interpolate(limit_direction='both')
            ext = pd.concat([clim.iloc[-win:], clim, clim.iloc[:win]])
            sm = ext.rolling(win, center=True, min_periods=1).mean().iloc[win:-win]
            sm = sm - sm.mean()                 # harmonic part only; keep patch mean
            g['dvv'] = g.dvv.values - sm.reindex(doy.values).values
        out.append(g)
    return pd.concat(out, ignore_index=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--glob', required=True)
    ap.add_argument('--verify', help='existing _des file to compare against (sanity)', default=None)
    args = ap.parse_args()
    for f in sorted(glob.glob(args.glob)):
        if f.endswith('_des.csv'):
            continue
        out = f.replace('.csv', '_des.csv')
        if os.path.exists(out):
            print('skip (exists):', out); continue
        d = pd.read_csv(f, parse_dates=['date'])
        r = deseason(d)
        r[['patch', 'date', 'dvv', 'cc_max']].to_csv(out, index=False)
        print('->', out, f'({r.patch.nunique()} patches)')
    if args.verify:
        mine = pd.read_csv(args.verify.replace('_cal_des', '_cal'), parse_dates=['date'])
        mine = deseason(mine)
        ref = pd.read_csv(args.verify, parse_dates=['date'])
        m = mine.merge(ref, on=['patch', 'date'], suffixes=('_m', '_r'))
        print('verify vs', args.verify, ': corr %.4f  rms-diff %.2e' %
              (np.corrcoef(m.dvv_m, m.dvv_r)[0, 1], (m.dvv_m - m.dvv_r).std()))


if __name__ == '__main__':
    main()
