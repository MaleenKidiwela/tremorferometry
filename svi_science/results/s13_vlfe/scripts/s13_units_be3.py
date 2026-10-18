#!/usr/bin/env python
"""S13 Amendment 7 (2026-10-15 PT): between-ETS stacks at candidate times (Merlin's scan fallback). New files only.

usage: s13_units_be3.py days   -> tables/be3_blocks.csv (drawn 7-day blocks), tables/be3_candidates.csv (deduped between-ETS
                                  candidates P >= 0.1, all days, with class / drawn flag), tables/be3_pnsn_epochs.csv
                                  (PNSN epochs in the footprint, between ETS, with LFE-free flags), tables/be3_prep_days.csv
                                  (drawn days +-1 and PNSN days +-1), tables/be3_env_days.csv (PNSN days, for envelopes)
       s13_units_be3.py rows   -> data/rows_be3/be3_{hi,c030,c010,pnsn}.parquet (needs tables/be3_pnsn_origins.csv from
                                  s13_be3_env.py)
Rules: PLAN.md Amendment 7."""
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
from s13_units_v2 import nms_dedupe, slab2_depth

RES = f'{C.ROOT}/svi_science/results'
ROWS = f'{C.DATA}/rows_be3'
T0, T1 = '2003-01-01', '2026-09-15'
CLASSES = [('C075', 0.75, 1.01), ('C050', 0.50, 0.75), ('C030', 0.30, 0.50), ('C010', 0.10, 0.30)]
ERAS = [('e1', 2003, 2008), ('e2', 2009, 2013), ('e3', 2014, 2026)]
SEED = 20261020
FRAC = 0.20
COLS = ['unit', 'to', 'P', 'lat', 'lon', 'dep', 'fam', 'date', 'full']
MIN_COL = 30


def majors():
    b = pd.read_csv(f'{RES}/b_ets_record/episodes_b.csv')
    b = b[(b.region == 'SVI') & b.major.astype(bool)][['episode', 'start', 'end']]
    p = pd.read_csv(f'{RES}/p3_bostock/episodes_bostock_final.csv')
    p = p[(p.start < '2010-01-01') & p.major.astype(bool)][['episode', 'start', 'end']]
    return pd.concat([p, b], ignore_index=True).sort_values('start').reset_index(drop=True)


def excluded_dates():
    out = set()
    for e in majors().itertuples():
        for d in pd.date_range(pd.Timestamp(e.start) - pd.Timedelta(days=2), pd.Timestamp(e.end) + pd.Timedelta(days=2)):
            out.add(d.strftime('%Y-%m-%d'))
    return out


def days():
    ex = excluded_dates()
    cal = pd.date_range(T0, T1)
    btw = np.array([d.strftime('%Y-%m-%d') not in ex for d in cal])
    blk = np.arange(len(cal)) // 7
    B = pd.DataFrame(dict(date=cal.strftime('%Y-%m-%d'), block=blk, between=btw))
    elig = B.groupby('block').between.any()
    elig = elig[elig].index.values
    rng = np.random.default_rng(SEED)
    drawn = np.sort(rng.choice(elig, int(round(FRAC * len(elig))), replace=False))
    B['drawn_block'] = B.block.isin(drawn)
    B['drawn_day'] = B.drawn_block & B.between
    B.to_csv(f'{C.TB}/be3_blocks.csv', index=False)
    dd = set(B.date[B.drawn_day])
    print('between days', int(btw.sum()), 'eligible blocks', len(elig), 'drawn', len(drawn), 'drawn between days', len(dd))
    # candidates
    R = C.svi_rows()
    R = R[(R.series == 'data') & (R.P >= 0.1) & (R.date >= T0) & (R.date <= T1) & ~R.date.isin(ex)].copy()
    sh = C.origin_shift_svi(sorted(R.fam.unique()))
    R['to'] = R.t + R.fam.map(sh).fillna(0.0)
    n0 = len(R)
    keep = nms_dedupe(R.to.values, R.P.values)
    R = R[keep].copy()
    F = C.svi_families()
    R = R.join(F[['lat', 'lon', 'depth_km']], on='fam').rename(columns={'depth_km': 'dep'}).dropna(subset=['lat'])
    R['cls'] = ''
    for nm, a, b in CLASSES:
        R.loc[(R.P >= a) & (R.P < b), 'cls'] = nm
    R['drawn'] = R.date.isin(dd)
    R[['date', 'fam', 'P', 't', 'to', 'lat', 'lon', 'dep', 'cls', 'drawn', 'src']].to_csv(f'{C.TB}/be3_candidates.csv', index=False)
    print('candidates P>=0.1 between ETS', n0, 'after +-2 s dedupe', int(keep.sum()), '; per class (all / drawn):',
          {c: (int((R.cls == c).sum()), int(((R.cls == c) & R.drawn).sum())) for c, *_ in CLASSES})
    # PNSN epochs
    from matplotlib.path import Path
    T = pd.read_csv(f'{C.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv')
    poly = pd.read_csv(f'{C.TB}/scan_footprint_polygon.csv')
    T = T[Path(poly[['lon', 'lat']].values).contains_points(T[['lon', 'lat']].values)].copy()
    T['date'] = pd.to_datetime(T.time).dt.strftime('%Y-%m-%d')
    T = T[(T.date >= T0) & (T.date <= T1) & ~T.date.isin(ex)].copy()
    T['t'] = pd.to_datetime(T.time).astype('datetime64[ns]').astype('int64').values / 1e9
    Rall = C.svi_rows(); Rall = Rall[Rall.series == 'data']
    for thr, nm in ((0.75, 'noLFE'), (0.1, 'noCand')):
        ts = np.sort(Rall[Rall.P >= thr].t.values)
        n = np.searchsorted(ts, T.t.values + 300, 'right') - np.searchsorted(ts, T.t.values - 300, 'left')
        T[nm] = n == 0
    T['dep'] = slab2_depth(T.lat.values, T.lon.values)
    T.to_csv(f'{C.TB}/be3_pnsn_epochs.csv', index=False)
    pdays = sorted(set(T.date[T.noLFE]))
    pd.Series(pdays, name='date').to_csv(f'{C.TB}/be3_env_days.csv', index=False)
    alld = set()
    for d in sorted(dd | set(pdays)):
        for k in (-1, 0, 1):
            alld.add((pd.Timestamp(d) + pd.Timedelta(days=k)).strftime('%Y-%m-%d'))
    pd.Series(sorted(alld), name='date').to_csv(f'{C.TB}/be3_prep_days.csv', index=False)
    print('PNSN epochs between ETS in footprint', len(T), 'noLFE', int(T.noLFE.sum()), 'noCand', int(T.noCand.sum()),
          'PNSN days', len(pdays), '; prep days', len(alld))


def units_for(X, name):
    """main class unit (full) + era / tremor-day column units (thrust only)."""
    parts = [X.assign(unit=name, full=True)]
    yr = X.date.str[:4].astype(int)
    for tag, a, b in ERAS:
        g = X[(yr >= a) & (yr <= b)]
        if len(g) >= MIN_COL:
            parts.append(g.assign(unit=f'{name}_{tag}', full=False))
    tdays = set(pd.read_csv(f'{C.TB}/be3_pnsn_epochs.csv').date)
    m = yr >= 2010
    for tag, sel in (('trem', X.date.isin(tdays)), ('notrem', ~X.date.isin(tdays))):
        g = X[m & sel]
        if len(g) >= MIN_COL:
            parts.append(g.assign(unit=f'{name}_{tag}', full=False))
    return parts


def rows():
    os.makedirs(ROWS, exist_ok=True)
    R = pd.read_csv(f'{C.TB}/be3_candidates.csv')
    R = R[R.drawn]
    files = {'be3_hi': ['C075', 'C050'], 'be3_c030': ['C030'], 'be3_c010': ['C010']}
    for fn, cls in files.items():
        parts = []
        for c in cls:
            parts += units_for(R[R.cls == c], c)
        X = pd.concat(parts, ignore_index=True)[COLS]
        X.to_parquet(f'{ROWS}/{fn}.parquet', index=False)
        print(fn, X.groupby('unit', sort=False).size().to_dict())
    O = pd.read_csv(f'{C.TB}/be3_pnsn_origins.csv')
    parts = []
    for nm, flag in (('PNSN_noLFE', 'noLFE'), ('PNSN_noCand', 'noCand')):
        g = O[O[flag].astype(bool)]
        g = g[nms_dedupe(g.to.values, g.env_peak.values, win=60.0)]
        parts.append(pd.DataFrame(dict(unit=nm, to=g.to.values, P=1.0, lat=g.lat.values, lon=g.lon.values, dep=g.dep.values,
                                       fam=-1, date=g.date.values, full=True)))
    BE = pd.read_parquet(f'{C.DATA}/rows_v2/svi_be_v2.parquet')
    parts.append(BE[BE.unit == 'BE_all'].assign(unit='BE_broad', full=True)[COLS])
    X = pd.concat(parts, ignore_index=True)[COLS]
    X.to_parquet(f'{ROWS}/be3_pnsn.parquet', index=False)
    print('be3_pnsn', X.groupby('unit', sort=False).size().to_dict())


if __name__ == '__main__':
    {'days': days, 'rows': rows}[sys.argv[1]]()
