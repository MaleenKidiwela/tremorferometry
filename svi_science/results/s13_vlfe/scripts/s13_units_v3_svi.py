#!/usr/bin/env python
"""S13 Amendment 6 step 5 (PLAN.md, 2026-10-14 PT): SVI net_mad terciles within P >= 0.90 for 2004-07, 2005-09,
2012-08, 2013-09. Pooled rows of data/rows/svi_<ep>.parquet (unchanged) joined to the catalogue net_mad (same row
construction as s13_common.svi_rows + s13_units.load('svi')). Writes data/rows_v3/svi_nmt_<ep>.parquet (units nm_lo,
nm_mid, nm_hi) and tables/nmt_svi_composition.csv."""
import os, sys
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
from s13_units_v2 import COLS

RES = f'{C.ROOT}/svi_science/results'


def rows_with_netmad():
    fc = ['run', 'date', 'series', 'fam', 'P_same', 'keep_xfam8', 'time_utc', 'net_mad']
    A = pd.read_parquet(f'{RES}/s0/analysis_table.parquet', columns=fc)
    A = A[(A.run.astype(str) == 'prod') & A.keep_xfam8].copy()
    A = A.assign(series=A.series.astype(str), P=A.P_same, counted_once=True,
                 t=A.time_utc.astype('datetime64[ns]').astype('int64') / 1e9)[['series', 'fam', 'P', 't', 'counted_once', 'net_mad']]
    B = pd.read_parquet(f'{C.ROOT}/data/catalogue/svi_catalogue_v1_4_backward.parquet',
                        columns=['series', 'fam', 'P', 't', 'primary', 'keep_xfam8', 'counted_once', 'net_mad'])
    B = B[B.primary & (B.keep_xfam8 == 1)].copy().assign(series=lambda x: x.series.astype(str))
    R = pd.concat([B[['series', 'fam', 'P', 't', 'counted_once', 'net_mad']], A], ignore_index=True)
    R = R[(R.series == 'data') & R.counted_once & (R.P >= 0.75)].copy()
    R['fam'] = R.fam.astype(int)
    sh = C.origin_shift_svi(sorted(R.fam.unique()))
    R['to'] = R.t + R.fam.map(sh).fillna(0.0)
    return R


def main():
    os.makedirs(f'{C.DATA}/rows_v3', exist_ok=True)
    R = rows_with_netmad()
    R['key'] = list(zip(R.fam, np.round(R.to, 2)))
    nm = R.drop_duplicates('key').set_index('key').net_mad
    comp = []
    for ep in ('2004-07', '2005-09', '2012-08', '2013-09'):
        X = pd.read_parquet(f'{C.DATA}/rows/svi_{ep}.parquet')
        Pq = X[X.unit == 'pooled'].copy()
        Pq['net_mad'] = [nm.get(k, np.nan) for k in zip(Pq.fam, np.round(Pq.to, 2))]
        matched = float(Pq.net_mad.notna().mean())
        H = Pq[(Pq.P >= 0.90) & Pq.net_mad.notna()].copy()
        q = np.quantile(H.net_mad, [1 / 3, 2 / 3])
        H['tert'] = np.where(H.net_mad <= q[0], 'nm_lo', np.where(H.net_mad <= q[1], 'nm_mid', 'nm_hi'))
        parts = [H[H.tert == u].assign(unit=u, full=False) for u in ('nm_lo', 'nm_mid', 'nm_hi')]
        Y = pd.concat(parts, ignore_index=True)[COLS]
        Y.to_parquet(f'{C.DATA}/rows_v3/svi_nmt_{ep}.parquet', index=False)
        rho = spearmanr(Pq.P, Pq.net_mad, nan_policy='omit').statistic
        rho9 = spearmanr(H.P, H.net_mad).statistic
        for u in ('nm_lo', 'nm_mid', 'nm_hi'):
            g = H[H.tert == u]
            comp.append(dict(episode=ep, unit=u, n=len(g), net_mad_med=g.net_mad.median(), net_mad_min=g.net_mad.min(),
                             net_mad_max=g.net_mad.max(), P_med=g.P.median(), P_mean=g.P.mean(), matched_frac=matched,
                             n_pooled=len(Pq), n_P90=len(H), rho_P_netmad_all=rho, rho_P_netmad_P90=rho9))
    Cm = pd.DataFrame(comp); Cm.to_csv(f'{C.TB}/nmt_svi_composition.csv', index=False)
    print(Cm.to_string())


if __name__ == '__main__':
    main()
