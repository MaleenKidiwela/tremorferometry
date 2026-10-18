"""S15 Q4 summary: recovery of GF sources per node (2012, 2024) vs GF C0/C1 references and the empirical C0/C1 sources at the
same times; real VLFE counts near each node."""
import json, sys
import numpy as np, pandas as pd
import s15_common as C
sys.path.insert(0, f'{C.CAT}/scripts')
import cat_common as K

OUT = f'{C.S15}/data/q4'
RUN = {2012: 'prod2012', 2024: 'prodfy2024'}


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def main():
    Q1 = pd.read_csv(f'{C.OUT_T}/q1_events.csv')
    rows = []
    for y, run in RUN.items():
        L = pd.read_parquet(f'{OUT}/inj_list_{y}.parquet'); R = pd.read_parquet(f'{OUT}/inj_{y}.parquet')
        J = L.merge(R, on='inj_id', how='left')
        n_err = int(J['err'].notna().sum()) if 'err' in J else 0
        J = K.rescore_inj(J)
        J['rec'] = ((J.mf_best_z >= 5.5) & (J.mf_second_z >= 5.5) & J.battery_e2.fillna(False).astype(bool)).astype(int)
        J['mf_only'] = ((J.mf_best_z >= 5.5) & (J.mf_second_z >= 5.5)).astype(int)
        J.to_parquet(f'{OUT}/inj_scored_{y}.parquet')
        M = pd.read_csv(f'{OUT}/sources_{y}.csv').set_index('source')
        # empirical C0 at the same times (base_inj_id) and empirical C1 at all times, Mw 3.6/3.8
        E = pd.read_parquet(f'{C.CAT}/data/{run}/inj_list.parquet').merge(pd.read_parquet(f'{C.CAT}/data/{run}/inj_v2.parquet'), on='inj_id')
        E = K.rescore_inj(E)
        E['rec'] = ((E.mf_best_z >= 5.5) & (E.mf_second_z >= 5.5) & E.battery_e2.fillna(False).astype(bool)).astype(int)
        for src in list(M.index) + ['C0_emp', 'C1_emp']:
            for mw in (3.6, 3.8):
                for per in ('ETS', 'between', 'all'):
                    if src.endswith('_emp'):
                        s = E[(E.source == src[:2]) & (E.Mw == mw) & E.period.isin(['ETS', 'between'])]    # whole-year empirical source
                    else:
                        s = J[(J.source == src) & (J.Mw == mw)]
                    if per != 'all':
                        s = s[s.period == per]
                    k, n = int(s.rec.sum()), len(s)
                    lo, hi = wilson(k, n)
                    lat = M.loc[src, 'lat'] if src in M.index else np.nan; lon = M.loc[src, 'lon'] if src in M.index else np.nan
                    rows.append(dict(year=y, source=src, Mw=mw, period=per, n=n, recovered=k, recovery=k / n if n else np.nan, ci_lo=lo, ci_hi=hi,
                                     mf_only=float(s.mf_only.mean()) if 'mf_only' in s and n else np.nan, lat=lat, lon=lon,
                                     node=M.loc[src, 'node'] if src in M.index else np.nan))
        print(y, 'errors', n_err)
    T = pd.DataFrame(rows)
    # ratio to the GF C0/C1 reference mean, same year/Mw/period
    ref = T[T.source.isin(['G0', 'G1'])].groupby(['year', 'Mw', 'period']).recovery.mean().rename('ref_G01')
    T = T.merge(ref, on=['year', 'Mw', 'period'], how='left'); T['ratio_to_G01'] = T.recovery / T.ref_G01
    # real catalogue VLFEs within 16 km of each node (VLFE location) and with concurrent-LFE anchor within 16 km
    with_lfe = Q1.n_lfe_10m.values >= 1
    nn = []
    for src in T.source.unique():
        if src.endswith('_emp'):
            continue
        r = T[T.source == src].iloc[0]
        x, y_ = C.xy(r.lat, r.lon)
        d = np.hypot(Q1.vx - x, Q1.vy - y_); da = np.hypot(Q1.anchor_x - x, Q1.anchor_y - y_)
        nn.append(dict(source=src, n_vlfe_loc_16km=int((d <= 16).sum()), n_vlfe_anchor_16km=int(((da <= 16) & with_lfe).sum())))
    T = T.merge(pd.DataFrame(nn), on='source', how='left')
    T.to_csv(f'{C.OUT_T}/q4_recovery.csv', index=False)
    pd.set_option('display.width', 250)
    print(T[T.period == 'all'][['year', 'source', 'node', 'lat', 'lon', 'Mw', 'n', 'recovery', 'ci_lo', 'ci_hi', 'mf_only', 'ratio_to_G01',
                                'n_vlfe_loc_16km', 'n_vlfe_anchor_16km']].round(3).to_string())
    print(T[(T.period == 'ETS') & (T.Mw == 3.8)][['year', 'source', 'n', 'recovery', 'ratio_to_G01']].round(3).to_string())


if __name__ == '__main__':
    main()
