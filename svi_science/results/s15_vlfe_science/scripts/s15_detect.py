"""S15 Q1 detection column (PLAN §3.9; column only per §0): recovery of C0/C1 injections vs days since local LFE onset."""
import json
import numpy as np, pandas as pd
import s15_common as C
import cat_common as K

SRC = {'C0': (48.141757352279875, -123.52741640459564), 'C1': (48.29464879935246, -123.61541420640553)}


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def main():
    F, _ = C.families(); L = C.lfes(); E = C.episodes()
    byk = {k: g.t.values for k, g in L.groupby('key')}
    near = {}
    for s, (la, lo) in SRC.items():
        x, y = C.xy(la, lo)
        near[s] = F.key.values[np.hypot(F.x - x, F.y - y) <= 20]
    rows = []
    for r in C.RUNS:
        J = pd.read_parquet(f'{C.CAT}/data/{r}/inj_list.parquet').merge(pd.read_parquet(f'{C.CAT}/data/{r}/inj_v2.parquet'), on='inj_id')
        J = K.rescore_inj(J)
        J = J[J.source.isin(['C0', 'C1']) & (J.period == 'ETS')].copy()
        J['rec'] = ((J.mf_best_z >= 5.5) & (J.mf_second_z >= 5.5) & J.battery_e2.fillna(False).astype(bool)).astype(int)
        J['run'] = r
        rows.append(J[['run', 'source', 'Mw', 'epoch', 'rec']])
    J = pd.concat(rows, ignore_index=True)
    ep = np.full(len(J), None, object); onset = np.full(len(J), np.nan)
    for _, e in E.iterrows():
        m = (J.epoch >= e.t0) & (J.epoch < e.t1)
        ep[m.values] = e.episode
        for s in SRC:
            ms = m.values & (J.source.values == s)
            if not ms.any():
                continue
            f = [byk[k][(byk[k] >= e.t0) & (byk[k] < e.t1)].min() for k in near[s] if k in byk and np.any((byk[k] >= e.t0) & (byk[k] < e.t1))]
            if f:
                onset[ms] = min(f)
    J['episode'] = ep; J['onset'] = onset
    J['days_since_local_onset'] = (J.epoch - J.onset) / C.DAY
    J['bin'] = pd.cut(J.days_since_local_onset, [-1e9, 0, 2, 5, 1e9], right=False, labels=['<0', '0-2', '2-5', '>=5'])
    J.loc[J.onset.isna(), 'bin'] = np.nan
    J.to_csv(f'{C.OUT_T}/detect_injections.csv', index=False)
    out = []
    for lab, sub in [('Mw3.6', J[J.Mw == 3.6]), ('Mw3.8', J[J.Mw == 3.8]), ('Mw3.6+3.8', J[J.Mw.isin([3.6, 3.8])]),
                     ('Mw3.8_big9', J[(J.Mw == 3.8) & J.episode.isin(C.BIG9)])]:
        for b in ['<0', '0-2', '2-5', '>=5', 'no_local_LFE']:
            s = sub[sub.bin.astype(str) == b] if b != 'no_local_LFE' else sub[sub.onset.isna() & sub.episode.notna()]
            lo, hi = wilson(s.rec.sum(), len(s))
            out.append(dict(set=lab, bin=b, n=len(s), recovery=s.rec.mean() if len(s) else np.nan, ci_lo=lo, ci_hi=hi))
        a = sub[sub.bin.astype(str) == '0-2'].rec; b_ = sub[sub.bin.astype(str).isin(['2-5', '>=5'])].rec
        out.append(dict(set=lab, bin='rel_diff_0-2_vs_>=2', n=len(a) + len(b_),
                        recovery=(a.mean() - b_.mean()) / b_.mean() if len(b_) else np.nan, ci_lo=np.nan, ci_hi=np.nan))
    T = pd.DataFrame(out)
    T.to_csv(f'{C.OUT_T}/detect_recovery_vs_local_day.csv', index=False)
    # per-episode recovery (episodes with VLFE runs; Mw 3.8)
    P = J[J.Mw == 3.8].groupby(['episode', 'bin'], observed=True).rec.agg(['size', 'mean']).reset_index()
    P.to_csv(f'{C.OUT_T}/detect_recovery_per_episode.csv', index=False)
    print(T.round(3).to_string())


if __name__ == '__main__':
    main()
