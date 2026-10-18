#!/usr/bin/env python
"""S13 scale-up: build the rows files (one per episode / burst set) for s13_stack2.py (PLAN.md Amendment 1).

usage: s13_units.py {svi,wa}
Writes data/rows/<region>_<episode>.parquet, data/rows/<region>_between.parquet, tables/<region>_episodes.csv
(episode, start, end, n_pooled, days list for prep: tables/<region>_prep_days.csv)."""
import os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_green as G

RES = f'{C.ROOT}/svi_science/results'
ROWS = f'{C.DATA}/rows'
MINN_FAM = 20
CENTRE = {'svi': (48.40, -123.60), 'wa': (47.75, -123.10)}


def along(lat, lon, reg):
    n, e = G.ned_offset(CENTRE[reg][0], CENTRE[reg][1], np.asarray(lat), np.asarray(lon))
    s = np.radians(315.0)
    return (n * np.cos(s) + e * np.sin(s)) / 1e3          # km toward N45W


def wa_families():
    F = pd.read_csv(f'{C.ROOT}/lfe_extension_wa/results/wa_prod/setup/families.csv')
    F = F.rename(columns={'fid': 'fam', 'lat_use': 'lat', 'lon_use': 'lon'})
    F['dep_missing'] = F.depth_km.isna()
    F['depth_km'] = F.depth_km.fillna(35.0)
    return F.set_index('fam')


def wa_shift(F):
    T = pd.read_csv(f'{C.ROOT}/lfe_extension_wa/results/wa_prod/setup/tpl_library.csv')
    co = pd.read_csv(f'{C.TB}/station_coords_all.csv', index_col=0)
    out = {}
    for f, g in T.groupby('fid'):
        if f not in F.index:
            continue
        v = []
        for r in g.itertuples():
            if r.netsta in co.index and np.isfinite(r.s_off):
                n, e = G.ned_offset(F.loc[f, 'lat'], F.loc[f, 'lon'], co.loc[r.netsta, 'lat'], co.loc[r.netsta, 'lon'])
                v.append(r.s_off - np.sqrt(n ** 2 + e ** 2 + (F.loc[f, 'depth_km'] * 1e3) ** 2) / 1e3 / C.VS_KMS)
        out[f] = float(np.median(v)) if v else 0.0
    return pd.Series(out)


def load(reg):
    if reg == 'svi':
        R = C.svi_rows()
        R = R[(R.series == 'data') & R.counted_once & (R.P >= 0.75)].copy()
        F = C.svi_families()
        sh = C.origin_shift_svi(sorted(R.fam.unique()))
        R['to'] = R.t + R.fam.map(sh).fillna(0.0)
        b = pd.read_csv(f'{RES}/b_ets_record/episodes_b.csv')
        b = b[(b.region == 'SVI') & b.major.astype(bool)][['episode', 'start', 'end']]
        p = pd.read_csv(f'{RES}/p3_bostock/episodes_bostock_final.csv')
        p = p[(p.start < '2010-01-01') & p.major.astype(bool)][['episode', 'start', 'end']]
        E = pd.concat([p, b], ignore_index=True)
    else:
        c = pd.read_parquet(f'{C.ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet',
                            columns=['time', 'date', 'family', 'region', 'P_own_pnet', 'selected_v2'])
        c = c[c.selected_v2 & (c.region == 'NW')]
        F = wa_families()
        R = pd.DataFrame(dict(fam=c.family.astype(int).values, P=c.P_own_pnet.values,
                              t=c.time.astype('datetime64[ns]').astype('int64').values / 1e9,
                              date=c.date.astype(str).values))
        R = R[R.fam.isin(F.index)]
        sh = wa_shift(F)
        R['to'] = R.t + R.fam.map(sh).fillna(0.0)
        E = pd.read_csv(f'{C.ROOT}/lfe_extension_wa/results/episode_rates/episodes_nw_final.csv')
        E = E[E.major.astype(bool)][['episode', 'start', 'end']]
    R = R.join(F[['lat', 'lon', 'depth_km']], on='fam').rename(columns={'depth_km': 'dep'})
    R = R.dropna(subset=['lat'])
    return R, F, E.sort_values('start').reset_index(drop=True)


def main(reg):
    os.makedirs(ROWS, exist_ok=True)
    R, F, E = load(reg)
    fx = pd.Series(along(F.lat.values, F.lon.values, reg), index=F.index)
    q = np.quantile(fx.loc[sorted(R.fam.unique())], [1 / 3, 2 / 3])
    R['zone'] = 1 + (R.fam.map(fx) > q[0]).astype(int) + (R.fam.map(fx) > q[1]).astype(int)
    ep_rows, days = [], set()
    for e in E.itertuples():
        a = (pd.Timestamp(e.start) - pd.Timedelta(days=2)).strftime('%Y-%m-%d')
        b = (pd.Timestamp(e.end) + pd.Timedelta(days=2)).strftime('%Y-%m-%d')
        U = R[(R.date >= a) & (R.date <= b)]
        if len(U) < 30:
            ep_rows.append(dict(episode=e.episode, start=e.start, end=e.end, n=len(U), built=False)); continue
        parts = [U.assign(unit='pooled', full=True)]
        dd = pd.to_datetime(U.date).dt.day
        parts += [U[dd % 2 == 1].assign(unit='pooled_odd', full=False), U[dd % 2 == 0].assign(unit='pooled_even', full=False)]
        for z in (1, 2, 3):
            if (U.zone == z).sum() >= MINN_FAM:
                parts.append(U[U.zone == z].assign(unit=f'zone{z}', full=False))
        for f, g in U.groupby('fam'):
            if len(g) >= MINN_FAM:
                parts.append(g.assign(unit=f'fam{f}', full=False))
        X = pd.concat(parts, ignore_index=True)[['unit', 'to', 'P', 'lat', 'lon', 'dep', 'fam', 'date', 'full']]
        X.to_parquet(f'{ROWS}/{reg}_{e.episode}.parquet', index=False)
        for d in pd.date_range(a, b):
            for k in (-1, 0, 1):
                days.add((d + pd.Timedelta(days=k)).strftime('%Y-%m-%d'))
        ep_rows.append(dict(episode=e.episode, start=e.start, end=e.end, n=len(U), n_units=X.unit.nunique(), built=True,
                            w0=a, w1=b))
    pd.DataFrame(ep_rows).to_csv(f'{C.TB}/{reg}_episodes.csv', index=False)
    # between-ETS sets (q5 / s7 / a_wa_between inputs)
    if reg == 'svi':
        sets = {'SVI_BE': pd.read_csv(f'{RES}/s7/between_events.csv'),
                'SVI_BE_broad': pd.read_parquet(f'{RES}/s7/between_broad_events.parquet')}
        famcol = 'fam'
    else:
        sets = {'WA_BE': pd.read_csv(f'{RES}/a_wa_between/between_events.csv'),
                'WA_BE_broad': pd.read_parquet(f'{RES}/a_wa_between/between_broad_events.parquet')}
        famcol = 'fid'
    parts = []
    for name, S in sets.items():
        t = pd.to_datetime(S.time_utc).astype('datetime64[ns]').astype('int64').values / 1e9
        fam = S[famcol].astype(int).values
        U = pd.DataFrame(dict(unit=name, fam=fam, P=1.0, date=pd.to_datetime(S.time_utc).dt.strftime('%Y-%m-%d').values,
                              to=t + pd.Series(fam).map(R.groupby('fam').apply(lambda g: (g.to - g.t).median())).fillna(0).values,
                              full=True))
        U = U.join(F[['lat', 'lon', 'depth_km']], on='fam').rename(columns={'depth_km': 'dep'}).dropna(subset=['lat'])
        parts.append(U)
        for d in U.date.unique():
            for k in (-1, 0, 1):
                days.add((pd.Timestamp(d) + pd.Timedelta(days=k)).strftime('%Y-%m-%d'))
    X = pd.concat(parts, ignore_index=True)[['unit', 'to', 'P', 'lat', 'lon', 'dep', 'fam', 'date', 'full']]
    X.to_parquet(f'{ROWS}/{reg}_between.parquet', index=False)
    pd.Series(sorted(days), name='date').to_csv(f'{C.TB}/{reg}_prep_days.csv', index=False)
    print(reg, 'episodes built', sum(r['built'] for r in ep_rows), 'of', len(ep_rows), 'prep days', len(days),
          'between', X.groupby('unit').size().to_dict())


if __name__ == '__main__':
    main(sys.argv[1])
