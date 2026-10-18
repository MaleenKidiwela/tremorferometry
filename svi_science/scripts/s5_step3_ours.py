"""S5 step 3 — tidal sensitivity in our 2014-2026 SVI catalogue (only run after step 2 PASS).
Rules: results/s5/PLAN_S5.md section 5.  Outputs results/s5/step3_*.
"""
import json, sys, numpy as np, pandas as pd
from pathlib import Path
import s5_lib as L
from s5_step2_bostock import DT_EDGES, DT_NAMES, hist_specs

ROOT = Path('/home/jovyan/tremorferometry')
R = ROOT / 'svi_science/results'
OUT = R / 's5'
COMPS = ['tau_Pa', 'sigma_n_Pa', 'tau_updip_Pa', 'coulomb_0p4', 'body_tau_Pa', 'ocean_tau_Pa']


def codes_dt(x):
    c = np.digitize(x, DT_EDGES[1:-1]); c[~np.isfinite(x)] = -1
    return c


def front_arrivals(fam_as):
    """PNSN tremor arrival per family per migrating SVI episode (s3b rule: 10th pct per 5-km along-strike bin, >=20)."""
    T = pd.read_parquet(R / 's3/tremor_svi.parquet')
    K = pd.read_csv(R / 's3/episodes_kinematics.csv'); K = K[(K.region == 'SVI') & (K.migrates == True)]
    rows = []
    for e in K.episode:
        tr = T[T.ep == e].assign(b=lambda d: np.floor(d.as_km / 5) * 5)
        g = tr.groupby('b').time.apply(lambda s: s.sort_values().iloc[int(0.1 * (len(s) - 1))] if len(s) >= 20 else pd.NaT).dropna()
        for f, a in fam_as.items():
            b = np.floor(a / 5) * 5
            if b in g.index:
                rows.append(dict(fam=f, episode=e, t_front=g[b]))
    return pd.DataFrame(rows), K.episode.tolist()


def main():
    s2 = json.load(open(OUT / 'step2_result.json'))
    assert s2['STEP2_PASS'], 'step 2 did not pass: step 3 must not run'
    F = pd.read_csv(R / 's0/families.csv')
    q1, q2 = F.downdip_km.quantile([1 / 3, 2 / 3])
    F['d3'] = np.where(F.downdip_km >= q2, 2, np.where(F.downdip_km <= q1, 0, 1))
    A = pd.read_parquet(R / 's0/analysis_table.parquet', columns=['series', 'fam', 'P_same', 'keep_xfam8', 'time_utc', 'episode', 'major',
                                                                   'flag_ep_2017_03', 'flag_ep_2025_10'])
    A['episode'] = A.episode.astype(str).replace('nan', '')
    S = L.Stress(pd.read_parquet(R / 's5_setup/tidal_stress_svi_hourly.parquet', columns=['fid', 'time'] + COMPS), COMPS)
    fam_as = F.set_index('fid').along_strike_km.to_dict()
    FR, mig = front_arrivals(fam_as)
    FR.to_csv(OUT / 'step3_front_arrivals.csv', index=False)
    names_all = {}

    def subset(series, pmin, once):
        m = (A.series == series) & (A.P_same >= pmin)
        if once:
            m &= A.keep_xfam8
        d = A[m].rename(columns={'fam': 'fid', 'time_utc': 't'})[['fid', 't', 'episode', 'major', 'flag_ep_2017_03', 'flag_ep_2025_10', 'P_same']]
        return d.reset_index(drop=True)

    out_tabs = []; summary = {}
    sets = [('data_075_once', 'data', 0.75, True), ('data_050_once', 'data', 0.5, True), ('data_085_once', 'data', 0.85, True),
            ('data_075_raw', 'data', 0.75, False), ('c1_075_once', 'c1', 0.75, True), ('c1_050_once', 'c1', 0.5, True),
            ('c1_all_raw', 'c1', -1.0, False)]
    want = sys.argv[1:] or [x[0] for x in sets]
    for name, series, pmin, once in sets:
        if name not in want:
            continue
        ev = subset(series, pmin, once)
        es = L.EventSet(ev, S, COMPS); ev = es.ev
        n = len(ev)
        yr = ev.t.dt.year.values
        half = np.where(yr <= 2019, 0, 1)
        d3 = ev.fid.map(F.set_index('fid').d3).values.astype(int)
        grp = ev.fid.map(F.set_index('fid').group).values.astype(int) - 1
        # own onset (10th pct) in major episodes
        mj = ev.major.values & (ev.episode.values != '')
        o = ev[mj].groupby(['fid', 'episode']).t.quantile(0.10).rename('t0').reset_index()
        ev = ev.merge(o, on=['fid', 'episode'], how='left')
        dt_own = np.where(mj, (ev.t - ev.t0).dt.total_seconds().values / 86400, np.nan)
        ev = ev.merge(FR.rename(columns={'fam': 'fid'}), on=['fid', 'episode'], how='left')
        dt_front = (ev.t - ev.t_front).dt.total_seconds().values / 86400
        in_mig = ev.episode.isin(mig).values
        dt_own_mig = np.where(in_mig, dt_own, np.nan)
        noflag = (~(ev.flag_ep_2017_03 | ev.flag_ep_2025_10)).values.astype(int)
        cf = codes_dt(dt_front); co = codes_dt(dt_own); com = codes_dt(dt_own_mig)
        labels = {'all': (np.zeros(n, int), 1), 'half': (half, 2), 'depth3': (d3, 3), 'group': (grp, 12),
                  'noflag': (noflag, 2), 'major': (mj.astype(int), 2), 'depth_half': (d3 * 2 + half, 6),
                  'dtbin': (cf, 7), 'dtown_major': (co, 7), 'dtown_mig': (com, 7),
                  'dtfront_half': (np.where(cf >= 0, cf * 2 + half, -1), 14),
                  'dtfront_depth': (np.where(cf >= 0, cf * 3 + d3, -1), 21)}
        nm = {'all': ['all'], 'half': ['2014-2019', '2020-2026'], 'depth3': ['shallow', 'mid', 'deep'],
              'group': [f'G{i+1}' for i in range(12)], 'noflag': ['flagged ep', 'not flagged'], 'major': ['outside major', 'major'],
              'depth_half': [f'{a}|{b}' for a in ('shallow', 'mid', 'deep') for b in ('2014-19', '2020-26')],
              'dtbin': DT_NAMES, 'dtown_major': DT_NAMES, 'dtown_mig': DT_NAMES,
              'dtfront_half': [f'{a}|{b}' for a in DT_NAMES for b in ('2014-19', '2020-26')],
              'dtfront_depth': [f'{a}|{b}' for a in DT_NAMES for b in ('shallow', 'mid', 'deep')]}
        import pickle
        ck = OUT / f'step3_surrogates_{name}.pkl'
        if ck.exists():
            obs, surr, hobs, hsur = pickle.load(open(ck, 'rb'))
        else:
            obs, surr, hobs, hsur = L.run(es, labels, nsurr=1000, seed=5, hist_specs=hist_specs())
            pickle.dump((obs, surr, hobs, hsur), open(ck, 'wb'))
        T = pd.concat([L.summarize(obs, surr, ln, nm[ln], const=['m2', 'o1', 'k1', 's2']) for ln in labels], ignore_index=True)
        T.insert(0, 'set', name); out_tabs.append(T)
        q0 = es.quantities(np.zeros(n, np.int64))
        th = L.astro_th(((es.day * L.NS_D) + (np.random.default_rng(11).random(n) * L.NS_D).astype(np.int64)).astype('datetime64[ns]'))[:, 0] + es.arg['m2']
        sm = dict(n=n, random_time_M2_rayleigh=L.rayleigh_p(np.cos(th), np.sin(th)),
                  observed_M2_rayleigh=L.rayleigh_p(q0['cos_m2'], q0['sin_m2']))
        if name.startswith('data_075'):
            fe_front = np.where(np.isfinite(dt_front), ev.fid.values * 1000 + pd.factorize(ev.episode)[0], -1)
            fe_own = np.where(np.isfinite(dt_own), ev.fid.values * 1000 + pd.factorize(ev.episode)[0], -1)
            fe_own_mig = np.where(np.isfinite(dt_own_mig), fe_own, -1)
            sm['rise_front'] = L.rise_perm(es, np.nan_to_num(dt_front, nan=-99), fe_front, q0['cos_m2'])
            sm['rise_own_major'] = L.rise_perm(es, np.nan_to_num(dt_own, nan=-99), fe_own, q0['cos_m2'])
            sm['rise_own_mig'] = L.rise_perm(es, np.nan_to_num(dt_own_mig, nan=-99), fe_own_mig, q0['cos_m2'])
            for h in (0, 1):
                fh = np.where(half == h, fe_front, -1)
                sm[f'rise_front_half{h}'] = L.rise_perm(es, np.nan_to_num(dt_front, nan=-99), fh, q0['cos_m2'])
            sm['diff_half'] = L.diff_test(obs, surr, 'half', 1, 0)
            sm['diff_deep_vs_shallow'] = L.diff_test(obs, surr, 'depth3', 2, 0)
            for d in (0, 2):
                sm[f'diff_half_depth{d}'] = L.diff_test(obs, surr, 'depth_half', d * 2 + 1, d * 2)
            # late-minus-first difference per half from surrogate-corrected bins, for power
            if name == 'data_075_once':
                np.savez_compressed(OUT / 'step3_hists.npz', **{f'obs_{k}': v for k, v in hobs.items()}, **{f'sur_{k}': v for k, v in hsur.items()})
        summary[name] = sm
        T.to_csv(OUT / f'step3_table_{name}.csv', index=False)
        json.dump(sm, open(OUT / f'step3_result_{name}.json', 'w'), indent=1, default=float)
        print(name, n, json.dumps({k: v for k, v in sm.items()}, default=float)[:1500], flush=True)
    names = [x[0] for x in sets]
    if all((OUT / f'step3_result_{n}.json').exists() for n in names):
        pd.concat([pd.read_csv(OUT / f'step3_table_{n}.csv') for n in names], ignore_index=True).to_csv(OUT / 'step3_modulation_table.csv', index=False)
        json.dump({n: json.load(open(OUT / f'step3_result_{n}.json')) for n in names}, open(OUT / 'step3_result.json', 'w'), indent=1, default=float)
        print('merged step3 outputs')


if __name__ == '__main__':
    main()
