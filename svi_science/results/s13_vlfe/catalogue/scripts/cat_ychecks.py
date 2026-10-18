#!/usr/bin/env python
"""Amendment Y3/Y4: pooled scrambled-control tier rate (2011 bank, final runs 2010-2018) with Poisson interval, A1/A2
split, dispersion test; Q2/Tier A' on every control tier event. -> tables/y3_controls.csv, y3_verdict.json,
y4_control_tierA.csv, y4_verdict.json"""
import os, sys, json, subprocess
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K
from scipy.stats import chi2, poisson

RUNS = {2010: 'prod', 2011: 'prod', 2012: 'prod', 2013: 'prodq', 2014: 'prod', 2015: 'prod', 2016: 'prodq', 2017: 'prodq',
        2018: 'prodq'}
NEAR = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV', 'CN.LZB']
FAR = ['IU.COR', 'US.NEW', 'US.HAWA', 'US.MSO', 'US.WVOR']


def between_years(y, nm):
    env = dict(os.environ, PROD_YEAR=str(y), PILOT_SFX='_v2')
    if nm.endswith('q'):
        env['PROD_QC'] = '1'
    out = subprocess.run([sys.executable, '-c', 'import sys; sys.path.insert(0, "%s"); import cat_pilot as P; print(P.period_years()["between"])'
                          % os.path.dirname(os.path.abspath(__file__))], env=env, capture_output=True, text=True)
    return float(out.stdout.strip().splitlines()[-1])


def events(D, noise=False):
    if noise:                                                # Amendment Y7: noise-template (noise-pair) tier events
        x = D[D.noise_pair_tier & (D.period == 'between')].sort_values('epoch').assign(source='N', draw=lambda d: d.template.str[1:].astype(int))
    else:
        x = D[D.tier & D.source.isin(['A1', 'A2']) & (D.period == 'between')].sort_values('epoch')
    rows = []
    for (s, d), g in x.groupby(['source', 'draw']):
        g = g.assign(grp=(g.i.diff().fillna(1e9) > 30).cumsum())
        rows.append(g.sort_values('zmad').groupby('grp').tail(1))
    return pd.concat(rows) if rows else x.iloc[:0]


def amp(y, s, ep, refbp):
    SD = f'{K.S13}/data/scan' if y == 2011 else f'{K.S13}/data/scan_{y}'
    f = f'{SD}/bp_{y}g/{s}.npy' if s in NEAR else f'{refbp or SD + "/bp2_ref"}/{s}.npy'
    if not os.path.isfile(f):
        return np.nan
    X = np.load(f, mmap_mode='r'); ts = int(round(ep - pd.Timestamp(f'{y}-01-01').value / 1e9))
    x = np.asarray(X[2, ts:ts + 180], float) * 1e9
    return float(np.nanmax(np.abs(x))) if np.isfinite(x).sum() > 150 else np.nan


def main():
    import scan_refs as SR
    cc = SR.all_coords()
    per, allev = [], []
    for y, nm in RUNS.items():
        D = pd.read_parquet(f'{K.DATA}/{nm}{y}/rows_scored_v2.parquet')
        nn = int(D[D.kind == 'noise'].template.nunique())
        T = between_years(y, nm)
        evs = events(D)
        evn = events(D, noise=True)
        ev = pd.concat([evs, evn], ignore_index=True)
        ev['ctrl'] = np.where(ev.source == 'N', 'noise_template', 'scrambled')
        per.append(dict(year=y, run=nm, between_yr=T, n_noise_templates=nn, n_noise_events=len(evn),
                        noise_rate_per_tpl_yr=len(evn) / (nn * T) if nn else np.nan,
                        n_events=len(evs), n_A1=int((evs.source == 'A1').sum()), n_A2=int((evs.source == 'A2').sum()),
                        rate=len(evs) / (16 * T)))
        refbp = (f'{K.S13}/data/scan/bp2_ref_v2' if y == 2011 else f'{K.S13}/data/scan_2012/bp2_ref_v2') if y in (2011, 2012) else None
        st = pd.read_csv(f'{K.TB}/{nm}{y}_R_stormdays_v2.csv'); sd = dict(zip(st.day, st.storm_day))
        ev = ev.assign(year=y, run=nm, day=pd.to_datetime(ev.epoch, unit='s').dt.date.astype(str))
        ev['cluster_day'] = ev.groupby(['source', 'draw', 'day']).i.transform('count') >= 2 if len(ev) else []
        ev['storm_day'] = ev.day.map(sd).fillna(False)
        fo = []
        for r in ev.itertuples():
            live = [s for s in NEAR if s in cc.index and np.isfinite(amp(y, s, r.epoch, refbp))]
            n2 = sorted(live, key=lambda s: K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], r.lat_best, r.lon_best))[:2]
            nm_ = max([amp(y, s, r.epoch, refbp) for s in n2] or [np.nan]); fm = np.nanmax([amp(y, s, r.epoch, refbp) for s in FAR])
            fo.append(fm / nm_ if nm_ and np.isfinite(nm_) else np.nan)
        ev['far_over_near2'] = fo
        ev['q2'] = (ev.vr_thrust >= 0.12) & (ev.far_over_near2 < 0.5) & np.where(ev.lat_best > 48.1, ev.ratio_GNW <= 0.551, True)
        ev['tierA'] = ev.q2 & ~(ev.storm_day & ev.cluster_day)
        allev.append(ev)
    P = pd.DataFrame(per)
    tot, T16 = P.n_events.sum(), (16 * P.between_yr).sum()
    rate = tot / T16
    lo, hi = poisson.ppf(0.025, tot) / T16 if tot else 0.0, poisson.ppf(0.975, tot) / T16
    lo = (chi2.ppf(0.025, 2 * tot) / 2) / T16 if tot else 0.0; hi = (chi2.ppf(0.975, 2 * tot + 2) / 2) / T16
    E = rate * 16 * P.between_yr
    X2 = float(((P.n_events - E) ** 2 / E).sum()); pdisp = float(chi2.sf(X2, len(P) - 1))
    P['expected'] = E
    P.to_csv(f'{K.TB}/y3_controls.csv', index=False)
    V3 = dict(pooled_rate_per_yr=rate, ci95=[lo, hi], n_events=int(tot), catalogue_years=float(T16),
              rate_A1=float(P.n_A1.sum() / T16), rate_A2=float(P.n_A2.sum() / T16), chi2=X2, df=len(P) - 1, p_dispersion=pdisp,
              verdict='per-year differences consistent with noise; 2016-17 = low-completeness years' if pdisp > 0.05 else '2016-17 flagged elevated noise')
    json.dump(V3, open(f'{K.TB}/y3_verdict.json', 'w'), indent=1)
    A = pd.concat(allev, ignore_index=True)
    A['time_pt'] = pd.to_datetime(A.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
    A[['year', 'run', 'ctrl', 'time_pt', 'source', 'draw', 'template', 'zmad', 'z_best', 'vr_thrust', 'ratio_far', 'ratio_GNW', 'lat_best',
       'far_over_near2', 'storm_day', 'cluster_day', 'q2', 'tierA']].to_csv(f'{K.TB}/y4_control_tierA.csv', index=False)
    S_ = A[A.ctrl == 'scrambled']; N_ = A[A.ctrl == 'noise_template']
    py = S_.groupby('year').tierA.sum().reindex(P.year, fill_value=0)
    ntA = int(S_.tierA.sum()); rA = ntA / T16
    TN = float((P.n_noise_templates * P.between_yr).sum()); nN = int(P.n_noise_events.sum()); nNA = int(N_.tierA.sum())
    ci = lambda k, T: [(chi2.ppf(0.025, 2 * k) / 2) / T if k else 0.0, (chi2.ppf(0.975, 2 * k + 2) / 2) / T]
    P['noise_tierA'] = N_.groupby('year').tierA.sum().reindex(P.year, fill_value=0).values
    P['noise_tierA_rate_per_tpl_yr'] = P.noise_tierA / (P.n_noise_templates * P.between_yr)
    P['scr_tierA'] = py.values
    P.to_csv(f'{K.TB}/y3_controls.csv', index=False)
    VN = dict(primary='noise-template control (Amendment Y7)', noise_template_years=TN,
              noise_raw_tier_events=nN, noise_raw_tier_rate_per_tpl_yr=nN / TN, noise_raw_ci95=ci(nN, TN),
              noise_tierA_events=nNA, noise_tierA_rate_per_tpl_yr=nNA / TN, noise_tierA_ci95=ci(nNA, TN),
              per_year_noise_raw={int(r.year): int(r.n_noise_events) for r in P.itertuples()},
              per_year_noise_tierA={int(r.year): int(r.noise_tierA) for r in P.itertuples()})
    json.dump(VN, open(f'{K.TB}/y7_noise_template_verdict.json', 'w'), indent=1)
    print(json.dumps(VN, indent=1))
    V4 = dict(control_tierA_events=ntA, control_tierA_rate_per_yr=rA,
              ci95=[(chi2.ppf(0.025, 2 * ntA) / 2) / T16 if ntA else 0.0, (chi2.ppf(0.975, 2 * ntA + 2) / 2) / T16],
              per_year={int(k): int(v) for k, v in py.items()}, accept=bool(rA <= 0.3),
              q2_pass_fraction=float(A.q2.mean()) if len(A) else np.nan)
    json.dump(V4, open(f'{K.TB}/y4_verdict.json', 'w'), indent=1)
    pd.set_option('display.width', 250)
    print(P.round(3).to_string(index=False)); print(json.dumps(V3, indent=1)); print(json.dumps(V4, indent=1))


if __name__ == '__main__':
    main()
