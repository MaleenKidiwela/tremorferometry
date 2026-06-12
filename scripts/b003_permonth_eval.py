#!/usr/bin/env python
"""B003 expansion: which densifying way is correct? Judge by RESOLVED CELL-MONTHS (per-month coverage),
not the static pooled count. Three B003 variants, each combined with the other 34 stations' existing
coverage; count 0.1-deg cells reaching >=3 stations WITH DATA in each month. Split ETS vs quiet months.
Anchor-independent (pure coverage bookkeeping). Run after the pilot's expanded dv/v exists.
"""
import pandas as pd, numpy as np, glob, re

GRID = 0.10
def cellof(p):
    la, lo = p.split('__')[0].split('_'); return (round(round(float(la)/GRID)*GRID,3), round(round(float(lo)/GRID)*GRID,3))

def cov_from(csv, sta):
    """-> set of (cell, 'YYYY-MM') that station images (has a dv/v point) per month."""
    d = pd.read_csv(csv, usecols=['patch','date']); d['ym'] = pd.to_datetime(d.date).dt.to_period('M').astype(str)
    d['cell'] = d.patch.map(cellof)
    return set(map(tuple, d[['cell','ym']].drop_duplicates().itertuples(index=False, name=None)))

# other 34 stations' per-month coverage (their existing 2-4s products)
others = {}
for f in sorted(glob.glob('data/daily_dvv_*_2to4_cal.csv')):
    s = re.match(r'.*daily_dvv_(.+)_2to4_cal', f).group(1)
    if s in ('B003','B933'): continue
    for cm in cov_from(f, s):
        others.setdefault(cm, set()).add(s)

# B003 variants
def b003_variant(csv, contfilter=None):
    d = pd.read_csv(csv, parse_dates=['date'])
    if contfilter is not None:
        d = d[d.patch.isin(contfilter)]
    d['ym'] = d.date.dt.to_period('M').astype(str); d['cell'] = d.patch.map(cellof)
    return set(map(tuple, d[['cell','ym']].drop_duplicates().itertuples(index=False, name=None)))

# continuity classification of the EXPANDED B003 family set (from the expanded dv/v itself)
def continuous_set(csv):
    d = pd.read_csv(csv, parse_dates=['date']); ad = pd.Index(d.date.unique()); keep=[]
    for p,g in d.groupby('patch'):
        dt=g.date.sort_values(); own=ad[(ad>=dt.iloc[0])&(ad<=dt.iloc[-1])]
        if len(dt)/max(len(own),1) >= 0.5 and g.cc_max.median() >= 0.6: keep.append(p)
    return set(keep)

BASE = 'data/daily_dvv_B003_2to4_calT.csv'        # current 81 (t0-anchor ok; coverage anchor-independent)
EXP  = 'data/daily_dvv_B003_2to4_calT_exp.csv'    # expanded (81 + 346 genuine)
cont = continuous_set(EXP)
variants = {
    'current (81)':        b003_variant(BASE),
    '+continuous (S1+S2)': b003_variant(EXP, contfilter=cont),
    '+all genuine (S1)':   b003_variant(EXP),
}
# ETS months for B003's latitude (~40N, southern): derive from tremor catalog rate near B003
cat = pd.read_csv('catalogs/pnsn_tremor_cascadia_full.csv', usecols=['time','lat','lon'])
cat['t']=pd.to_datetime(cat.time,errors='coerce'); cat=cat.dropna(subset=['t'])
near = cat[(cat.lat.between(39.5,41.0))]
mrate = near.groupby(near.t.dt.to_period('M').astype(str)).size()
ets_months = set(mrate[mrate >= mrate.quantile(0.80)].index)

print(f"{'variant':22s} {'B003 fam':>8s} {'cell-mo @>=3sta':>15s} {'ETS-mo':>7s} {'quiet-mo':>8s}")
n_b003_fam = {'current (81)':81, '+continuous (S1+S2)':len(cont), '+all genuine (S1)':pd.read_csv(EXP).patch.nunique()}
for name, b3 in variants.items():
    cov = {k:set(v) for k,v in others.items()}
    for cm in b3: cov.setdefault(cm, set()).add('B003')
    resolved = [cm for cm,ss in cov.items() if len(ss)>=3]
    ets = sum(1 for (c,m) in resolved if m in ets_months)
    quiet = len(resolved)-ets
    print(f"{name:22s} {n_b003_fam[name]:>8d} {len(resolved):>15d} {ets:>7d} {quiet:>8d}")
print("\n(resolved cell-MONTHS = #(0.1-deg cell, month) reaching >=3 stations WITH DATA that month;")
print(" ETS months = top-20% tremor-rate months at 39.5-41N. This is the per-month metric the advisor flagged.)")
