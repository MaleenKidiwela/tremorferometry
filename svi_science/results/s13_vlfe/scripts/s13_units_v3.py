#!/usr/bin/env python
"""S13 Amendment 6 step 3 (PLAN.md, 2026-10-14 PT): WA rows v3 -> data/rows_v3/ (new files; rows/ and rows_v2/ untouched).
to_f = catalogue time + 0.5 + offset_f (tables/wa_family_offsets.csv), greedy +-2 s dedupe on to_f over all families
(keep the highest P_own_pnet; same code as v2), Green's-function location from tables/wa3_location_flags.csv (duplicate
cluster -> best-located member; depth-suspect singletons -> Slab2 depth), zones = terciles of the along-strike position of
the GF locations of the kept families. Same majors, +-2 d windows, >= 30 LFE build rule and units as v2. Between: WA_BE,
WA_BE_broad with the same offsets and dedupe within each set.
Also: n_nb (pooled rows, self excluded) per episode for the step-3 acceptance, tables/wa3_dedupe.csv, wa3_episodes.csv,
wa3_families.csv, wa3_prep_days.csv, wa3_prep_missing.csv (days not yet in data/disp1sps_v2 for the WA stations)."""
import os, sys, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_units as U1
from s13_units_v2 import nms_dedupe, episode_units, wa_dupC, COLS, PRE
from s13_stack2_v2 import nb_counts

RES = f'{C.ROOT}/svi_science/results'
ROWS = f'{C.DATA}/rows_v3'
DISP2 = f'{C.DATA}/disp1sps_v2'


def main():
    os.makedirs(ROWS, exist_ok=True)
    O = pd.read_csv(f'{C.TB}/wa_family_offsets.csv', index_col=0)
    L = pd.read_csv(f'{C.TB}/wa3_location_flags.csv', index_col=0)
    G = L[['gf_lat', 'gf_lon', 'gf_dep', 'gf_source']].rename(columns={'gf_lat': 'lat', 'gf_lon': 'lon', 'gf_dep': 'depth_km'})
    G['dep_slab2'] = G.gf_source.str.contains('slab2')
    c = pd.read_parquet(f'{C.ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet',
                        columns=['det_id', 'time', 'date', 'family', 'region', 'P_own_pnet', 'selected_v2'])
    c = c[c.selected_v2 & (c.region == 'NW')]
    R = pd.DataFrame(dict(det_id=c.det_id.values, fam=c.family.astype(int).values, P=c.P_own_pnet.values,
                          t=c.time.astype('datetime64[ns]').astype('int64').values / 1e9, date=c.date.astype(str).values))
    n0 = len(R)
    R = R[R.fam.isin(O.index)].copy()
    R['offset'] = R.fam.map(O.offset_f)
    R['to'] = R.t + PRE + R.offset
    R['date'] = pd.to_datetime(R.to, unit='s').dt.strftime('%Y-%m-%d')   # date of the corrected origin
    R['dupC'] = wa_dupC(R.det_id.values)
    keep = nms_dedupe(R.to.values, R.P.values)
    D = R.assign(kept=keep)
    D.assign(year=D.date.str[:4]).groupby('year').agg(n=('kept', 'size'), kept=('kept', 'sum')).to_csv(f'{C.TB}/wa3_dedupe.csv')
    R = R[keep].copy()
    R = R.join(G[['lat', 'lon', 'depth_km', 'dep_slab2']], on='fam').rename(columns={'depth_km': 'dep'})
    R = R.dropna(subset=['lat', 'dep'])
    fx = pd.Series(U1.along(G.lat.values, G.lon.values, 'wa'), index=G.index)
    q = np.quantile(fx.loc[sorted(R.fam.unique())], [1 / 3, 2 / 3])
    R['zone'] = 1 + (R.fam.map(fx) > q[0]).astype(int) + (R.fam.map(fx) > q[1]).astype(int)
    E = pd.read_csv(f'{C.ROOT}/lfe_extension_wa/results/episode_rates/episodes_nw_final.csv')
    E = E[E.major.astype(bool)][['episode', 'start', 'end']].sort_values('start').reset_index(drop=True)
    ep, days = episode_units(R, E, 'wa', ROWS, 'wa')
    # step-3 acceptance numbers: n_nb per built episode (pooled rows, self excluded), same function as the stacker
    nb = []
    for e in ep[ep.built].episode:
        X = pd.read_parquet(f'{ROWS}/wa_{e}.parquet')
        to = X[X.unit == 'pooled'].to.values
        nb.append(dict(episode=e, **nb_counts(to, to)))
    NB = pd.DataFrame(nb)
    ep = ep.merge(NB, on='episode', how='left')
    ep.to_csv(f'{C.TB}/wa3_episodes.csv', index=False)
    fam = G.loc[sorted(R.fam.unique())].assign(along_km=fx, zone=R.groupby('fam').zone.first(), n_rows=R.groupby('fam').size(),
                                               offset_f=O.offset_f, flag=O.flag)
    fam.to_csv(f'{C.TB}/wa3_families.csv')
    # between-ETS sets
    parts = []
    for name, f in (('WA_BE', f'{RES}/a_wa_between/between_events.csv'),
                    ('WA_BE_broad', f'{RES}/a_wa_between/between_broad_events.parquet')):
        S = pd.read_csv(f) if f.endswith('.csv') else pd.read_parquet(f)
        S = S[S.fid.isin(O.index)]
        t = pd.to_datetime(S.time_utc).astype('datetime64[ns]').astype('int64').values / 1e9 + PRE + \
            S.fid.map(O.offset_f).values
        Pp = S.P_own_WAVE_PNET.values if 'P_own_WAVE_PNET' in S else np.ones(len(S))
        X = pd.DataFrame(dict(unit=name, fam=S.fid.astype(int).values, P=Pp, to=t, full=True,
                              date=pd.to_datetime(t, unit='s').strftime('%Y-%m-%d')))
        k = nms_dedupe(X.to.values, X.P.values)
        print(name, 'rows', len(X), 'kept after +-2 s dedupe', int(k.sum()))
        X = X[k].join(G[['lat', 'lon', 'depth_km']], on='fam').rename(columns={'depth_km': 'dep'}).dropna(subset=['lat', 'dep'])
        parts.append(X)
        for d in np.floor(X.to.values / 86400).astype(int):
            days.add(int(d))
    X = pd.concat(parts, ignore_index=True)[COLS]
    X.to_parquet(f'{ROWS}/wa_between.parquet', index=False)
    alld = sorted({d + k for d in days for k in range(-6, 7)})
    pd.Series([pd.Timestamp(d * 86400, unit='s').strftime('%Y-%m-%d') for d in alld], name='date').to_csv(
        f'{C.TB}/wa3_prep_days.csv', index=False)
    # days not yet prepped (any WA station that had a file on a neighbouring prepped day is not checked: we list days
    # with no file for any WA station at all, plus per-station counts)
    sta = [s for s in open(f'{C.OUT}/scripts/s13_fill_queue.sh').read().split('WA_STA=')[1].split('\n')[0].split(',')]
    have = set()
    for s in sta:
        for d in alld:
            day = pd.Timestamp(d * 86400, unit='s')
            if os.path.isfile(f'{DISP2}/{s}/{day.year}/{day.dayofyear:03d}.npz'):
                have.add(d)
    v2days = set(pd.to_datetime(pd.read_csv(f'{C.TB}/wa2_prep_days.csv').date).astype('datetime64[ns]').astype('int64') // 86400_000_000_000)
    new = [d for d in alld if d not in v2days]
    pd.Series([pd.Timestamp(d * 86400, unit='s').strftime('%Y-%m-%d') for d in new], name='date').to_csv(
        f'{C.TB}/wa3_prep_missing.csv', index=False)
    built = ep[ep.built]
    v = dict(catalogue_rows=n0, rows_with_offset=int(len(D)), kept=int(keep.sum()), removed_frac=float(1 - keep.mean()),
             episodes_built=int(len(built)), episodes_ge300=int((built.n >= 300).sum()),
             nb2_mean=float(built.n_nb2.mean()), nb10_mean=float(built.n_nb10.mean()), nb10_median=float(built.n_nb10.median()),
             nb_kern_mean=float(built.nb_kern.mean()), nb10_max=float(built.n_nb10.max()),
             prep_days=len(alld), prep_days_not_in_v2_list=len(new), days_with_any_file=len(have))
    v['accept_nb2'] = bool(v['nb2_mean'] == 0); v['accept_nb10'] = bool(v['nb10_mean'] < 0.15)
    v['accept_removed'] = bool(0.25 <= v['removed_frac'] <= 0.35)
    v['step3_accept'] = bool(v['accept_nb2'] and v['accept_nb10'] and v['accept_removed'])
    json.dump(v, open(f'{C.TB}/wa3_rows_verdict.json', 'w'), indent=1)
    print(json.dumps(v, indent=1))


if __name__ == '__main__':
    main()
