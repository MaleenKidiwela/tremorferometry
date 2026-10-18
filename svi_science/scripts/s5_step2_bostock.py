"""S5 step 2 — POSITIVE CONTROL: Royer et al. (2015) on Bostock's 2003-2013 LFEs with our tide model.
Rules frozen in results/s5/PLAN_S5.md section 4.  Outputs results/s5/step2_*.
"""
import json, numpy as np, pandas as pd
from pathlib import Path
import s5_lib as L

ROOT = Path('/home/jovyan/tremorferometry')
OUT = ROOT / 'svi_science/results/s5'
COMPS = ['tau_Pa', 'sigma_n_Pa', 'tau_updip_Pa', 'coulomb_0p4', 'body_tau_Pa', 'ocean_tau_Pa']
DT_EDGES = [-np.inf, 1, 2, 3, 4, 6, 10, np.inf]
DT_NAMES = ['<1d', '1-2d', '2-3d', '3-4d', '4-6d', '6-10d', '>=10d']
FRAME = dict(lat0=48.48373913043479, lon0=-123.74852173913044, az=51.15105595138788)


def downdip(lat, lon):
    x = (lon - FRAME['lon0']) * 111.195 * np.cos(np.radians(FRAME['lat0'])); y = (lat - FRAME['lat0']) * 111.195
    a = np.radians(FRAME['az'])
    return x * np.sin(a) + y * np.cos(a)


def depth_thirds():
    F = pd.read_csv(ROOT / 'svi_science/results/s0/families.csv')
    q1, q2 = F.downdip_km.quantile([1 / 3, 2 / 3])
    return q1, q2


def dt_codes(ev, onset_col):
    x = (ev.t - ev[onset_col]).dt.total_seconds().values / 86400
    c = np.digitize(x, DT_EDGES[1:-1])
    c[~np.isfinite(x)] = -1
    return c


def hist_specs():
    ph = lambda q: np.arctan2(q['sin_m2'], q['cos_m2'])
    return {'tauz_all': (lambda q: q['tau_Pa_z'], np.linspace(-2.5, 2.5, 11), 'all'),
            'tauz_dt': (lambda q: q['tau_Pa_z'], np.linspace(-2.5, 2.5, 11), 'dtbin'),
            'phM2_all': (ph, np.linspace(-np.pi, np.pi, 13), 'all'),
            'phM2_dt': (ph, np.linspace(-np.pi, np.pi, 13), 'dtbin')}


def load_bostock():
    c = pd.read_csv(ROOT / 'lfe_catalog_repair/results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv')
    e = pd.read_csv(ROOT / 'lfe_catalog_repair/results/e5r_families_2012.csv')
    c = c.merge(e[['bkey', 'lat', 'lon', 'dep']], on=['lat', 'lon', 'dep'], how='left')
    assert c.bkey.notna().all()
    t = pd.to_datetime(dict(year=c.year, month=c.month, day=c.day, hour=c.hour, minute=c['min'])) + pd.to_timedelta(c.sec, unit='s')
    ev = pd.DataFrame({'fid': c.bkey.astype(int), 't': t.astype('datetime64[ns]'), 'lat': c.lat, 'lon': c.lon})
    W = pd.read_csv(OUT.parent / 'p3_bostock/episodes_bostock_final.csv')
    ev['ep'] = -1
    for i, w in W.iterrows():
        m = (ev.t >= pd.Timestamp(w.start)) & (ev.t < pd.Timestamp(w.end) + pd.Timedelta(days=1))
        ev.loc[m, 'ep'] = i
    return ev, W


def onsets(ev):
    g = ev[ev.ep >= 0].groupby(['fid', 'ep']).t
    o = pd.DataFrame({'t0_p10': g.quantile(0.10), 't0_first': g.min(), 'n_fe': g.size()}).reset_index()
    return ev.merge(o, on=['fid', 'ep'], how='left')


def make_labels(ev, ep_names, depth_cut):
    q1, q2 = depth_cut
    dd = downdip(ev.lat.values, ev.lon.values)
    d3 = np.where(dd >= q2, 2, np.where(dd <= q1, 0, 1))
    prod = set(pd.read_csv(ROOT / 'svi_science/results/s0/families.csv').fid)
    dt = dt_codes(ev, 't0_p10'); dtf = dt_codes(ev, 't0_first')
    L_ = {'all': (np.zeros(len(ev), int), 1),
          'fam115': (ev.fid.isin(prod).values.astype(int), 2),
          'episode': (ev.ep.values.astype(int), len(ep_names)),
          'dtbin': (dt, 7), 'dtbin_first': (dtf, 7),
          'depth3': (d3, 3), 'dt_depth': (np.where(dt >= 0, dt * 3 + d3, -1), 21),
          'dt_fam115': (np.where(ev.fid.isin(prod).values, dt, -1), 7)}
    return L_, dd


def main():
    acc = json.load(open(OUT / 'step1_acceptance.json'))
    assert acc['pass'], 'step 1 acceptance failed'
    ev, W = load_bostock()
    ev = onsets(ev)
    print('events', len(ev), 'families', ev.fid.nunique(), 'outside windows', int((ev.ep < 0).sum()), flush=True)
    S = L.Stress(pd.read_parquet(OUT / 'tidal_stress_bostock_hourly.parquet', columns=['fid', 'time'] + COMPS), COMPS)
    es = L.EventSet(ev, S, COMPS)
    ev = es.ev
    labels, dd = make_labels(ev, W.episode.tolist(), depth_thirds())
    # --- random continuous times on the same days: flatness of M2 phase
    rng = np.random.default_rng(11)
    tr = (es.day * L.NS_D + (rng.random(len(ev)) * L.NS_D).astype(np.int64)).astype('datetime64[ns]')
    th = L.astro_th(tr)[:, 0] + es.arg['m2']
    p_rand, R_rand = L.rayleigh_p(np.cos(th), np.sin(th))
    # --- main surrogate run
    import pickle
    ck = OUT / 'step2_surrogates.pkl'
    if ck.exists():
        obs, surr, hobs, hsur = pickle.load(open(ck, 'rb'))
    else:
        obs, surr, hobs, hsur = L.run(es, labels, nsurr=1000, seed=5, hist_specs=hist_specs())
        pickle.dump((obs, surr, hobs, hsur), open(ck, 'wb'))
    tabs = []
    names = {'all': ['all'], 'fam115': ['15 non-production', '115 production'], 'episode': W.episode.tolist(),
             'dtbin': DT_NAMES, 'dtbin_first': DT_NAMES, 'depth3': ['shallow', 'mid', 'deep'],
             'dt_depth': [f'{a}|{b}' for a in DT_NAMES for b in ('shallow', 'mid', 'deep')], 'dt_fam115': DT_NAMES}
    for ln in labels:
        tabs.append(L.summarize(obs, surr, ln, names[ln], const=['m2', 'o1', 'k1', 's2']))
    T = pd.concat(tabs, ignore_index=True)
    T.to_csv(OUT / 'step2_modulation_table.csv', index=False)
    # phase statistics for sigma_n (column)
    o = obs['all']; s = surr['all']; n = o['n'][0]
    sn = {c: dict(C=o[f'cos_sn_{c}'][0] / n, Cex=o[f'cos_sn_{c}'][0] / n - s[f'cos_sn_{c}'][:, 0].mean() / n,
                  Sex=o[f'sin_sn_{c}'][0] / n - s[f'sin_sn_{c}'][:, 0].mean() / n,
                  p_two=float((1 + np.sum(np.abs(s[f'cos_sn_{c}'][:, 0] / n - s[f'cos_sn_{c}'][:, 0].mean() / n) >=
                                          abs(o[f'cos_sn_{c}'][0] / n - s[f'cos_sn_{c}'][:, 0].mean() / n))) / 1001))
          for c in ('m2', 'o1')}
    # --- P2 rise (whole-day permutation within family-episode)
    fe = np.where(ev.ep.values >= 0, ev.fid.values * 100 + ev.ep.values, -1)
    q0 = es.quantities(np.zeros(len(ev), np.int64))
    dtd = (ev.t - ev.t0_p10).dt.total_seconds().values / 86400
    dtf = (ev.t - ev.t0_first).dt.total_seconds().values / 86400
    rise = {'p10_onset': L.rise_perm(es, dtd, fe, q0['cos_m2']),
            'first_event_onset': L.rise_perm(es, dtf, fe, q0['cos_m2']),
            'p10_onset_tau_z': L.rise_perm(es, dtd, fe, np.nan_to_num(q0['tau_Pa_z']))}
    np.savez_compressed(OUT / 'step2_hists.npz', **{f'obs_{k}': v for k, v in hobs.items()}, **{f'sur_{k}': v for k, v in hsur.items()})
    a = T[T.label == 'all'].iloc[0]
    dtt = T[T.label == 'dtbin'].set_index('name')
    P1 = bool(a.p_C_m2 < 0.01 and a.Cex_m2 > 0 and a.p_tau < 0.01 and a.dmean_tau_Pa > 0)
    P2 = bool(rise['p10_onset']['p_one'] < 0.01 and rise['p10_onset']['diff'] > 0)
    res = dict(n_events=int(len(ev)), n_families=int(ev.fid.nunique()), random_time_M2_rayleigh_p=p_rand, random_time_M2_R=R_rand,
               surrogate_mean_R_M2=float(np.hypot(s['cos_m2'][:, 0].mean() / n, s['sin_m2'][:, 0].mean() / n)),
               headline=a.to_dict(), sigma_n_phase=sn, rise=rise,
               theta0_by_dt=dtt[['n', 'Cex_m2', 'a_m2_pct', 'theta0_m2_deg', 'p_C_m2', 'dmean_tau_Pa', 'asig_tau_kPa']].to_dict('index'),
               P1_sign_pass=P1, P2_rise_pass=P2, STEP2_PASS=bool(P1 and P2))
    json.dump(res, open(OUT / 'step2_result.json', 'w'), indent=1, default=float)
    print(json.dumps({k: res[k] for k in ('n_events', 'random_time_M2_rayleigh_p', 'surrogate_mean_R_M2', 'rise', 'P1_sign_pass', 'P2_rise_pass', 'STEP2_PASS')}, indent=1, default=float))
    cols = ['label', 'name', 'n', 'Cex_m2', 'Csurr_sd_m2', 'a_m2_pct', 'theta0_m2_deg', 'p_C_m2', 'Cex_o1', 'a_o1_pct', 'theta0_o1_deg', 'p_C_o1',
            'dmean_tau_Pa', 'p_tau', 'asig_tau_kPa', 'dmean_sigma_n_Pa', 'p_sigma_n', 'dmean_tau_updip_Pa', 'p_tau_updip']
    with pd.option_context('display.width', 250, 'display.max_columns', 30):
        print(T[cols].round(4).to_string())


if __name__ == '__main__':
    main()
