#!/usr/bin/env python
"""Per-family detection-time fingerprint -> physical-identity label (anthropogenic vs natural).
From mf_*_all.csv detection times only (free, on disk). Computes per (station,family):
  weekday/weekend ratio, work-hours fraction (09-17 local, baseline 0.375), day/night ratio,
  midday-weekday concentration (blast signature). Classifies the contaminant axis:
    CULTURAL  : weekday-skewed (wkratio>1.15) AND work-hours-heavy (wh>0.45)
    BLAST     : strong weekday MIDDAY spike (wkratio>1.3 AND midday_wd>0.30)
    NATURAL   : week-flat (0.9-1.1) AND night-favored (wh<0.36)  [detectability signature of real EQ/LFE]
    WIND/THERMAL/OTHER : diurnal but week-flat, not night-favored
(Local time approximated as UTC-8 for PST; Cascadia.)
"""
import glob, re, sys, warnings
warnings.filterwarnings('ignore')
import numpy as np, pandas as pd

OFF = 8  # UTC -> PST hours
rows = []
for mf in sorted(glob.glob('data/mf_*_all.csv')):
    sta = re.match(r'data/mf_(.+)_all', mf).group(1).upper()
    if sta.endswith('REV') or 'EXP' in sta: continue
    try:
        m = pd.read_csv(mf, usecols=['template','time'])
    except Exception:
        continue
    t = pd.to_datetime(m.time, errors='coerce'); m = m[t.notna()]; t = t[t.notna()]
    lt = t - pd.Timedelta(hours=OFF)
    m = m.assign(hr=lt.dt.hour.values, wd=(lt.dt.dayofweek.values < 5))
    print(f'[{sta}] {len(m):,} det, {m.template.nunique()} fams', flush=True)
    for fam, g in m.groupby('template'):
        n = len(g)
        if n < 200: continue
        wk = g.wd.mean()                                  # fraction on weekdays (baseline 5/7=0.714)
        wkratio = (wk/0.714) / ((1-wk)/0.286 + 1e-9)      # weekday-rate / weekend-rate
        wh = ((g.hr >= 9) & (g.hr < 17)).mean()           # work-hours frac (baseline 8/24=0.333)
        dn = ((g.hr >= 12) & (g.hr <= 21)).sum() / max(((g.hr >= 0) & (g.hr <= 9)).sum(), 1)
        midday_wd = (g.wd & (g.hr >= 9) & (g.hr < 16)).mean()
        rows.append(dict(station=sta, fam=fam, n=n, wkratio=wkratio, workhr=wh, dn=dn, midday_wd=midday_wd))
A = pd.DataFrame(rows)

def classify(r):
    if r.wkratio > 1.3 and r.midday_wd > 0.30: return 'BLAST'
    if r.wkratio > 1.15 and r.workhr > 0.42:   return 'CULTURAL'
    if 0.9 <= r.wkratio <= 1.12 and r.dn < 1.1: return 'NATURAL-like'
    if r.dn > 1.4:                              return 'DIURNAL-other'  # wind/thermal/hydro
    return 'FLAT-mixed'
A['identity'] = A.apply(classify, axis=1)
A.to_csv('data/family_fingerprint.csv', index=False)
print('\n=== identity counts ===')
print(A.identity.value_counts().to_string())
print('\n=== per-station identity mix ===')
print(pd.crosstab(A.station, A.identity).to_string())
print(f'\n-> data/family_fingerprint.csv ({len(A)} families)')
print('FINGERPRINT_DONE')
