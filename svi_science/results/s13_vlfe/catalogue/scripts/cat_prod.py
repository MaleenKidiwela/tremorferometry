#!/usr/bin/env python
"""Production 2010-2018 (PLAN_CAT Amendment S) helpers.
  stations Y   -> appends/updates tables/prod_stations.csv: 2011-bank carriers with >= 180 v3 day files in Y
  summary Y    -> tables/prodY_episodes.csv, prodY_ledger.csv, prodY_summary.json (needs the cat_pilot/cat_r outputs)
Run with PROD_YEAR=Y PILOT_SFX=_v2 for summary."""
import os, sys, json, glob
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K

CARRIERS = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.YOUB', 'CN.LZB', 'UW.LRIV', 'UW.DOSE', 'UW.GNW', 'UW.SP2', 'TA.A04D',
            'TA.B05D', 'TA.D03D', 'TA.D04D', 'UW.WISH', 'US.NLWA', 'CC.OBSR', 'UW.ELW', 'UW.LEBA', 'UW.RATT', 'CN.PFB',
            'UW.OTR', 'UW.HDW']                         # Amendment X: OTR (broadband 2020-10+), HDW (2023+) for relay banks


def run_stations(y):
    f = f'{K.TB}/prod_stations.csv'
    T = pd.read_csv(f) if os.path.isfile(f) else pd.DataFrame(columns=['year', 'sta', 'v3_days'])
    T = T[T.year != y]
    rows = []
    for s in CARRIERS:
        n = len(glob.glob(f'{K.S13}/data/disp1sps_v3/{s}/{y}/*.npz'))
        if n >= 180:
            rows.append(dict(year=y, sta=s, v3_days=n))
    T = pd.concat([T, pd.DataFrame(rows)], ignore_index=True).sort_values(['year', 'sta'])
    T.to_csv(f, index=False)
    print(y, [r['sta'] for r in rows])


def run_missing(y):
    """carriers with >= 180 broadband days on disk (station-era table) but < 180 v3 day files -> print for prep."""
    E = pd.read_csv(f'{K.TB}/era_station_year.csv')
    e = E[(E.year == y) & (E.days_bb >= 180) & E.sta.isin(CARRIERS)]
    miss = [s for s in e.sta if len(glob.glob(f'{K.S13}/data/disp1sps_v3/{s}/{y}/*.npz')) < 180]
    print(','.join(miss))


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def run_summary(y):
    import cat_pilot as P
    SFX = P.SFX
    D = pd.read_parquet(f'{P.PD}/rows_scored{SFX}.parquet')
    R = json.load(open(f'{K.TB}/{P.NAME}{y}_R_verdict{SFX}.json'))
    wins = P.CFG.get('wins', [])
    t = pd.to_datetime(D.epoch, unit='s')
    D['episode'] = ''
    for a, b, ep in wins:
        D.loc[(t >= pd.Timestamp(a)) & (t < pd.Timestamp(b)), 'episode'] = ep

    def dd(x):
        tt = np.sort(x.i.values)
        return int(1 + np.sum(np.diff(tt) > 30)) if len(tt) else 0
    J = pd.read_parquet(f'{P.PD}/inj_list.parquet').merge(pd.read_parquet(f'{P.PD}/inj{SFX}.parquet'), on='inj_id', how='left')
    J = K.rescore_inj(J)                                     # Amendment Y14: f'
    J['rec'] = (J.mf_best_z >= P.LEVEL) & (J.mf_second_z >= P.LEVEL) & J.battery_e2.fillna(False).astype(bool)
    tj = pd.to_datetime(J.epoch, unit='s')
    rows = []
    for a, b, ep in wins:
        x = D[D.tier & (D.episode == ep)]
        real = x[x.source == 'real'].sort_values('epoch')
        ev = real.assign(g=(real.i.diff().fillna(1e9) > 30).cumsum()).sort_values('zmad').groupby('g').tail(1) if len(real) else real
        m0 = ev.m0_cols.abs() if 'm0_cols' in ev else pd.Series(dtype=float)
        cc = [dd(x[(x.source == s) & (x.draw == q)]) for s in ('A1', 'A2') for q in range(P.NDRAW)]
        o = dict(year=y, episode=ep, major=ep in P.CFG.get('majors', []), start=a, end=b, tier_events=len(ev), ctrl_mean=float(np.mean(cc)), ctrl_max=int(np.max(cc)),
                 sum_M0_Nm=float(m0.sum()) if len(m0) else 0.0,
                 Mw_of_sum=float((np.log10(m0.sum()) - 9.1) / 1.5) if len(m0) and m0.sum() > 0 else np.nan,
                 with_LFEs=int((ev.lfe_status == 'with LFEs').sum()) if len(ev) else 0)
        for mw in (3.4, 3.6, 3.8, 4.0):
            g = J[(J.Mw == mw) & (tj >= pd.Timestamp(a)) & (tj < pd.Timestamp(b))]
            k, n = int(g.rec.sum()), len(g); lo, hi = wilson(k, n)
            o[f'compl_{mw}'] = k / n if n else np.nan; o[f'compl_{mw}_ci'] = f'{lo:.2f}-{hi:.2f}' if n else ''; o[f'n_inj_{mw}'] = n
        rows.append(o)
    E = pd.DataFrame(rows); E.to_csv(f'{K.TB}/{P.NAME}{y}_episodes{SFX}.csv', index=False)
    # events table
    ev = D[D.tier & (D.source == 'real')].sort_values('epoch')
    ev = ev.assign(g=(ev.i.diff().fillna(1e9) > 30).cumsum()).sort_values('zmad').groupby('g').tail(1).sort_values('epoch')
    ev = ev.assign(time_pt=pd.to_datetime(ev.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z'),
                   Mw_scan=(np.log10(ev.m0_cols.abs()) - 9.1) / 1.5)
    ev[['time_pt', 'period', 'episode', 'template', 'agree', 'zmad', 'z_best', 'Mw_scan', 'ratio_far', 'vr_thrust', 'lat_best',
        'lon_best', 'nearest_family', 'lfe_status', 'n_lfe_10m', 'i', 'epoch']].to_csv(f'{K.TB}/{P.NAME}{y}_events{SFX}.csv', index=False)
    # known-answer ledger
    known = []
    V = pd.read_csv(f'{K.TB}/e3_events.csv')
    for r in V.itertuples():
        known.append(('2011-12 ' + r.role, float(r.epoch), r.time_pt))
    for py in (2013, 2018):
        for sfx in ('_v2', ''):
            f = f'{K.TB}/pilot{py}_ets_events{sfx}.csv'
            if os.path.isfile(f):
                for r in pd.read_csv(f).itertuples():
                    known.append((f'pilot{py}{sfx} ETS tier', float(r.epoch), r.time_pt))
                break
    L = []
    allr = D[(D.source == 'real')]
    for kind, ep_, tpt in known:
        if pd.Timestamp(ep_, unit='s').year != y:
            continue
        x = allr[(allr.epoch - ep_).abs() <= 30]
        L.append(dict(year=y, known=kind, time_pt=tpt, found_tier=bool(x.tier.any()), best_mf_z=float(x.zmad.max()) if len(x) else np.nan,
                      battery_e2=bool(x.battery_e2.any()) if len(x) else False, agree=bool((x.agree != '').any()) if len(x) else False))
    LG = pd.DataFrame(L); LG.to_csv(f'{K.TB}/{P.NAME}{y}_ledger{SFX}.csv', index=False)
    S = dict(year=y, stations=P.STAS, n_episodes=len(wins), tier_events_episodes=int(E.tier_events.sum()) if len(E) else 0,
             between_raw=R['between']['raw'], between_tierA=R['between']['tierA'], controls=R['controls'],
             completeness=R['completeness'], T1=R['T1'], T2=R['T2'], T3=R['T3'],
             ledger_found=f"{int(LG.found_tier.sum()) if len(LG) else 0}/{len(LG)}",
             flag_relay=bool(isinstance(R['completeness']['at_3p8'].get('ETS', 0), float) and min(R['completeness']['at_3p8'].values()) < 0.5))
    json.dump(S, open(f'{K.TB}/{P.NAME}{y}_summary{SFX}.json', 'w'), indent=1, default=str)
    pd.set_option('display.width', 250)
    print(json.dumps(S, indent=1, default=str)); print(E.round(3).to_string()); print(LG.to_string())


if __name__ == '__main__':
    m, y = sys.argv[1], int(sys.argv[2])
    {'stations': run_stations, 'summary': run_summary, 'missing': run_missing}[m](y)
