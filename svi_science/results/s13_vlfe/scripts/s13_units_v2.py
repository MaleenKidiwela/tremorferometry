#!/usr/bin/env python
"""S13 v2 row files (PLAN.md Amendment 5, 2026-10-14 PT; Merlin scale-up review steps 2-5). New files only:
data/rows_v2/*.parquet; data/rows/ and s13_units.py are not touched (the fill/refit queues read them).

usage: s13_units_v2.py wa | svi_pbin | svi_be
  wa        WA-NW rows: origin to = catalogue time + 0.5 s (no wa_shift; Merlin F1), cross-family dedupe (+-2 s, keep the
            highest P_own_pnet; Merlin F2), family depth < 20 km or missing -> Slab2 at the family lat/lon (flag; F3).
            Writes rows_v2/wa_<episode>.parquet, rows_v2/wa_between.parquet, tables/wa2_episodes.csv,
            tables/wa2_families.csv, tables/wa2_dedupe.csv, tables/wa2_prep_days.csv (LFE days +-6 for the +-5-day controls).
  svi_pbin  Merlin step 4: for SVI 2004-07, 2005-09, 2012-08, 2013-09 the pooled rows split by P in [0.75,0.90),
            [0.90,0.97), [0.97,1]; 2004-07 also top-3,000 by P and random-3,000. rows_v2/svi_pbin_<episode>.parquet.
  svi_be    Merlin step 5: SVI BE_broad variants (all with real P_same; minus 2020-09-13/14; minus 2015-10-08/10; minus the
            top-5 days; tau_next >= 30 d; 2014-2019; 2020-2026). rows_v2/svi_be_v2.parquet, tables/svi_be_burst_days.csv."""
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_units as U1          # read-only import (along(), CENTRE, MINN_FAM)

RES = f'{C.ROOT}/svi_science/results'
ROWS1 = f'{C.DATA}/rows'
ROWS = f'{C.DATA}/rows_v2'
PRE = 0.5                        # lfe_extension_wa/scripts/wa_common.py l.34 (template starts PRE s before S)
DEDUPE_S = 2.0
DEP_MIN = 20.0
COLS = ['unit', 'to', 'P', 'lat', 'lon', 'dep', 'fam', 'date', 'full']


def slab2_depth(lat, lon):
    from scipy.interpolate import RegularGridInterpolator
    S = pd.read_csv(f'{C.ROOT}/data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz', header=None, names=['lon', 'lat', 'z'])
    lo = np.sort(S.lon.unique()); la = np.sort(S.lat.unique())
    Z = S.pivot(index='lat', columns='lon', values='z').reindex(index=la, columns=lo).values
    f = RegularGridInterpolator((la, lo), Z, bounds_error=False, fill_value=np.nan)
    return -f(np.c_[np.asarray(lat), np.mod(np.asarray(lon), 360.0)])


def wa_families_v2():
    F = pd.read_csv(f'{C.ROOT}/lfe_extension_wa/results/wa_prod/setup/families.csv')
    F = F.rename(columns={'fid': 'fam', 'lat_use': 'lat', 'lon_use': 'lon'}).set_index('fam')
    F['depth_cat'] = F.depth_km
    F['slab2_km'] = slab2_depth(F.lat.values, F.lon.values)
    bad = F.depth_km.isna() | (F.depth_km < DEP_MIN)
    F['dep_slab2'] = bad
    F.loc[bad, 'depth_km'] = F.loc[bad, 'slab2_km']
    return F


def nms_dedupe(t, P, win=DEDUPE_S):
    """Keep the highest-P row, drop every other row within +-win s of a kept row (greedy, descending P; ties -> earlier)."""
    order = np.lexsort((t, -P))
    ts = np.sort(t); kept_sorted = np.zeros(len(t), bool)
    rank = np.empty(len(t), np.int64); rank[np.argsort(t, kind='stable')] = np.arange(len(t))
    keep = np.zeros(len(t), bool)
    for i in order:
        lo = np.searchsorted(ts, t[i] - win, 'left'); hi = np.searchsorted(ts, t[i] + win, 'right')
        if not kept_sorted[lo:hi].any():
            keep[i] = True; kept_sorted[rank[i]] = True
    return keep


def wa_dupC(det_ids):
    """Variant-C cross-family duplicate flag from post_v2/dedup/flags (column only)."""
    import glob
    want = pd.Index(det_ids)
    out = []
    for f in sorted(glob.glob(f'{C.ROOT}/lfe_extension_wa/results/wa_prod/post_v2/dedup/flags/flags_*.parquet')):
        x = pd.read_parquet(f, columns=['det_id', 'dup_C'])
        out.append(x[x.det_id.isin(want)])
    X = pd.concat(out).drop_duplicates('det_id').set_index('det_id').dup_C
    return X.reindex(want).values


def episode_units(R, E, reg, ROWS_OUT, prefix):
    ep_rows, days = [], set()
    for e in E.itertuples():
        a = (pd.Timestamp(e.start) - pd.Timedelta(days=2)).strftime('%Y-%m-%d')
        b = (pd.Timestamp(e.end) + pd.Timedelta(days=2)).strftime('%Y-%m-%d')
        Uu = R[(R.date >= a) & (R.date <= b)]
        if len(Uu) < 30:
            ep_rows.append(dict(episode=e.episode, start=e.start, end=e.end, n=len(Uu), built=False)); continue
        parts = [Uu.assign(unit='pooled', full=True)]
        dd = pd.to_datetime(Uu.date).dt.day
        parts += [Uu[dd % 2 == 1].assign(unit='pooled_odd', full=False), Uu[dd % 2 == 0].assign(unit='pooled_even', full=False)]
        for z in (1, 2, 3):
            if (Uu.zone == z).sum() >= U1.MINN_FAM:
                parts.append(Uu[Uu.zone == z].assign(unit=f'zone{z}', full=False))
        for f, g in Uu.groupby('fam'):
            if len(g) >= U1.MINN_FAM:
                parts.append(g.assign(unit=f'fam{f}', full=False))
        X = pd.concat(parts, ignore_index=True)[COLS]
        X.to_parquet(f'{ROWS_OUT}/{prefix}_{e.episode}.parquet', index=False)
        for d in np.floor(Uu.to.values / 86400).astype(int):
            days.add(int(d))
        ep_rows.append(dict(episode=e.episode, start=e.start, end=e.end, n=len(Uu), n_units=X.unit.nunique(), built=True,
                            w0=a, w1=b, n_dep_slab2=int(Uu.dep_slab2.sum()), n_dupC=int(Uu.dupC.fillna(False).sum())))
    return pd.DataFrame(ep_rows), days


def wa():
    os.makedirs(ROWS, exist_ok=True)
    c = pd.read_parquet(f'{C.ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet',
                        columns=['det_id', 'time', 'date', 'family', 'region', 'P_own_pnet', 'selected_v2'])
    c = c[c.selected_v2 & (c.region == 'NW')]
    F = wa_families_v2()
    R = pd.DataFrame(dict(det_id=c.det_id.values, fam=c.family.astype(int).values, P=c.P_own_pnet.values,
                          t=c.time.astype('datetime64[ns]').astype('int64').values / 1e9, date=c.date.astype(str).values))
    n0 = len(R)
    R = R[R.fam.isin(F.index)].copy()
    R['to'] = R.t + PRE
    R['dupC'] = wa_dupC(R.det_id.values)
    keep = nms_dedupe(R.to.values, R.P.values)
    D = R.assign(kept=keep)
    # dedupe log per year
    D.assign(year=D.date.str[:4]).groupby('year').agg(n=('kept', 'size'), kept=('kept', 'sum'),
                                                       dupC=('dupC', lambda x: int(x.fillna(False).sum()))).to_csv(f'{C.TB}/wa2_dedupe.csv')
    R = R[keep].copy()
    R = R.join(F[['lat', 'lon', 'depth_km', 'dep_slab2']], on='fam').rename(columns={'depth_km': 'dep'})
    R = R.dropna(subset=['lat', 'dep'])
    fx = pd.Series(U1.along(F.lat.values, F.lon.values, 'wa'), index=F.index)
    q = np.quantile(fx.loc[sorted(R.fam.unique())], [1 / 3, 2 / 3])
    R['zone'] = 1 + (R.fam.map(fx) > q[0]).astype(int) + (R.fam.map(fx) > q[1]).astype(int)
    E = pd.read_csv(f'{C.ROOT}/lfe_extension_wa/results/episode_rates/episodes_nw_final.csv')
    E = E[E.major.astype(bool)][['episode', 'start', 'end']].sort_values('start').reset_index(drop=True)
    ep, days = episode_units(R, E, 'wa', ROWS, 'wa')
    ep.to_csv(f'{C.TB}/wa2_episodes.csv', index=False)
    F.loc[sorted(R.fam.unique()), ['family_uid', 'lat', 'lon', 'depth_cat', 'slab2_km', 'depth_km', 'dep_slab2']].assign(
        along_km=fx, zone=R.groupby('fam').zone.first(), n_rows=R.groupby('fam').size()).to_csv(f'{C.TB}/wa2_families.csv')
    # between-ETS sets: same origin convention, same dedupe within each set
    parts = []
    for name, f in (('WA_BE', f'{RES}/a_wa_between/between_events.csv'),
                    ('WA_BE_broad', f'{RES}/a_wa_between/between_broad_events.parquet')):
        S = pd.read_csv(f) if f.endswith('.csv') else pd.read_parquet(f)
        t = pd.to_datetime(S.time_utc).astype('datetime64[ns]').astype('int64').values / 1e9 + PRE
        Pp = S.P_own_WAVE_PNET.values if 'P_own_WAVE_PNET' in S else np.ones(len(S))
        X = pd.DataFrame(dict(unit=name, fam=S.fid.astype(int).values, P=Pp, to=t, full=True,
                              date=pd.to_datetime(S.time_utc).dt.strftime('%Y-%m-%d').values))
        k = nms_dedupe(X.to.values, X.P.values)
        print(name, 'rows', len(X), 'kept after +-2 s dedupe', int(k.sum()))
        X = X[k].join(F[['lat', 'lon', 'depth_km']], on='fam').rename(columns={'depth_km': 'dep'}).dropna(subset=['lat', 'dep'])
        parts.append(X)
        for d in np.floor(X.to.values / 86400).astype(int):
            days.add(int(d))
    X = pd.concat(parts, ignore_index=True)[COLS]
    X.to_parquet(f'{ROWS}/wa_between.parquet', index=False)
    alld = sorted({d + k for d in days for k in range(-6, 7)})
    pd.Series([pd.Timestamp(d * 86400, unit='s').strftime('%Y-%m-%d') for d in alld], name='date').to_csv(
        f'{C.TB}/wa2_prep_days.csv', index=False)
    print('wa2: catalogue rows', n0, 'with family', int(D.shape[0]), 'kept', int(keep.sum()),
          f'({keep.mean():.3f}); dupC among input {int(D.dupC.fillna(False).sum())}; families depth->slab2',
          int(F.loc[sorted(R.fam.unique()), 'dep_slab2'].sum()), 'of', R.fam.nunique(),
          '; episodes built', int(ep.built.sum()), 'of', len(ep), '; prep days', len(alld))


def svi_pbin():
    os.makedirs(ROWS, exist_ok=True)
    rng = np.random.default_rng(C.SEED + 4)
    bins = [(0.75, 0.90, 'pbin_lo'), (0.90, 0.97, 'pbin_mid'), (0.97, 1.0001, 'pbin_hi')]
    for ep in ('2004-07', '2005-09', '2012-08', '2013-09'):
        R = pd.read_parquet(f'{ROWS1}/svi_{ep}.parquet')
        Pq = R[R.unit == 'pooled'].copy()
        parts = [Pq[(Pq.P >= a) & (Pq.P < b)].assign(unit=nm, full=False) for a, b, nm in bins]
        if ep == '2004-07':
            parts.append(Pq.nlargest(3000, 'P').assign(unit='top3000', full=False))
            parts.append(Pq.iloc[rng.choice(len(Pq), 3000, replace=False)].assign(unit='rand3000', full=False))
        X = pd.concat(parts, ignore_index=True)[COLS]
        X.to_parquet(f'{ROWS}/svi_pbin_{ep}.parquet', index=False)
        print(ep, X.groupby('unit').agg(n=('P', 'size'), Pmed=('P', 'median')).to_dict('index'))


def svi_be():
    os.makedirs(ROWS, exist_ok=True)
    S = pd.read_parquet(f'{RES}/s7/between_broad_events.parquet')
    S = S[S.series.astype(str) == 'data'] if 'series' in S else S
    F = C.svi_families()
    sh = C.origin_shift_svi(sorted(S.fam.astype(int).unique()))
    t = pd.to_datetime(S.time_utc).astype('datetime64[ns]').astype('int64').values / 1e9
    X = pd.DataFrame(dict(fam=S.fam.astype(int).values, P=S.P_same.values, date=S.date.astype(str).values,
                          to=t + pd.Series(S.fam.astype(int).values).map(sh).fillna(0.0).values,
                          tau_next=S.tau_next.values))
    X = X.join(F[['lat', 'lon', 'depth_km']], on='fam').rename(columns={'depth_km': 'dep'}).dropna(subset=['lat'])
    # check vs the first-pass rows (same events, same origin rule)
    R1 = pd.read_parquet(f'{ROWS1}/svi_between.parquet')
    R1 = R1[R1.unit == 'SVI_BE_broad']
    print('first-pass BE_broad rows', len(R1), 'v2', len(X), 'max |to diff| (sorted)',
          float(np.abs(np.sort(R1.to.values) - np.sort(X.to.values)).max()) if len(R1) == len(X) else 'n differs')
    cnt = X.groupby('date').size().sort_values(ascending=False)
    cnt.rename('n').to_csv(f'{C.TB}/svi_be_burst_days.csv')
    top5 = list(cnt.index[:5])
    yr = X.date.str[:4].astype(int)
    sets = {'BE_all': X,
            'BE_minus_2020-09-13_14': X[~X.date.isin(['2020-09-13', '2020-09-14'])],
            'BE_minus_2015-10-08_10': X[~X.date.isin(['2015-10-08', '2015-10-09', '2015-10-10'])],
            'BE_minus_top5days': X[~X.date.isin(top5)],
            'BE_tau_next_ge30': X[X.tau_next >= 30],
            'BE_2014_2019': X[yr <= 2019],
            'BE_2020_2026': X[yr >= 2020]}
    parts = [v.assign(unit=k, full=True) for k, v in sets.items()]
    Y = pd.concat(parts, ignore_index=True)[COLS]
    Y.to_parquet(f'{ROWS}/svi_be_v2.parquet', index=False)
    print('top5 days', top5, cnt.iloc[:5].tolist(), '; units', Y.groupby('unit').size().to_dict())


if __name__ == '__main__':
    {'wa': wa, 'svi_pbin': svi_pbin, 'svi_be': svi_be}[sys.argv[1]]()
