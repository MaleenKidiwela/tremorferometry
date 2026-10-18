#!/usr/bin/env python
"""Amendment R re-summary of a pilot (T1', T2', T3', controls, completeness with Wilson CIs and MF/battery split, storm /
cluster / Q2 columns, Tier A'). usage: PILOT=2013 python cat_r.py  -> tables/pilotY_R_*.csv|json, figures/pilotY_R_*.png"""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
import cat_common as K
import cat_pilot as P

Y = P.PILOT
SFX = P.SFX
PD = P.PD
NEAR = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV', 'CN.LZB']
REG = ['UW.GNW', 'UW.DOSE', 'UW.SP2', 'TA.D03D', 'TA.B05D', 'UW.WISH', 'US.NLWA']
FAR = ['IU.COR', 'US.NEW', 'US.HAWA', 'US.MSO', 'US.WVOR']
SD = f'{K.S13}/data/scan_{Y}' if Y != 2011 else f'{K.S13}/data/scan'
Y0 = pd.Timestamp(f'{Y}-01-01').value / 1e9


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


_A = {}


def arr(s):
    if s not in _A:
        f = f'{SD}/bp_{Y}g/{s}.npy' if s in NEAR else f'{SD}/bp2_ref/{s}.npy'
        _A[s] = np.load(f, mmap_mode='r') if os.path.isfile(f) else None
    return _A[s]


def amplitudes(epoch):
    ts = int(round(epoch - Y0)); out = {}
    for s in NEAR + REG + FAR:
        X = arr(s)
        if X is None or ts + 180 > X.shape[1]:
            continue
        x = np.asarray(X[2, ts:ts + 180], float) * 1e9
        if np.isfinite(x).sum() > 150:
            out[s] = float(np.nanmax(np.abs(x)))
    return out


def t3_event(ev, cc):
    from scipy.stats import spearmanr
    A = amplitudes(ev.epoch)
    R = {s: float(K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], ev.lat_best, ev.lon_best)) for s in A if s in cc.index}
    nm = [s for s in A if s in R and s not in FAR]
    near = [s for s in nm if R[s] <= 40]
    if len(near) < 2:
        near = sorted(nm, key=R.get)[:2]
    mid = [s for s in nm if 40 < R[s] <= 150 and s not in near]
    far = [s for s in FAR if s in A]
    r1 = np.median([A[s] for s in near]) / np.median([A[s] for s in mid]) if near and mid else np.nan
    r2 = np.median([A[s] for s in far]) / max(A[s] for s in near) if far and near else np.nan
    allS = [s for s in near + mid + far if s in R]
    rho = spearmanr([R[s] for s in allS], [A[s] for s in allS]).correlation if len(allS) >= 4 else np.nan
    ok = bool(r1 >= 1.5 and r2 <= 0.5 and rho <= -0.4)
    return dict(near=','.join(near), n_mid=len(mid), n_far=len(far), near_mid_ratio=r1, far_over_nearmax=r2, spearman=rho, pass_=ok,
                amps=json.dumps({s: round(A[s], 1) for s in allS}), dists=json.dumps({s: round(R[s]) for s in allS}))


def storm_table():
    from scipy.signal import butter, sosfiltfilt
    sos = butter(2, (0.1, 0.2), btype='band', fs=1.0, output='sos')
    rows = []
    for k in range(P.NDAY):
        d = P.T0 + pd.Timedelta(days=k)
        f = f'{K.S13}/data/disp1sps_v3/CN.PGC/{d.year}/{d.dayofyear:03d}.npz'
        v = np.nan
        if os.path.isfile(f):
            x = np.load(f)['d'][2].astype(float)
            if np.isfinite(x).mean() > 0.8:
                x = np.where(np.isfinite(x), x, 0.0); x = x - x.mean()
                y = sosfiltfilt(sos, x) * 1e9
                v = float(np.median([np.sqrt(np.mean(y[i:i + 600] ** 2)) for i in range(0, 86400, 600)]))
        rows.append(dict(day=str(d.date()), pgc_5_10s_rms_nm=v))
    S = pd.DataFrame(rows)
    S['pct'] = S.pgc_5_10s_rms_nm.rank(pct=True)
    S['storm_day'] = S.pct >= 0.85
    return S


def main():
    import hg_r3 as X
    import scan_refs as SR
    D = pd.read_parquet(f'{PD}/rows_scored{SFX}.parquet')
    V = json.load(open(f'{K.TB}/pilot{Y}_verdict{SFX}.json'))
    C = V['counts']
    PY = P.period_years()
    out = dict(pilot=Y, statistic=V['statistic'])
    # ---- events (dedupe +-30 s, keep highest MF z)
    def events(per):
        x = D[D.tier & (D.period == per) & (D.source == 'real')].sort_values('epoch')
        if x.empty:
            return x
        x = x.assign(g=(x.i.diff().fillna(1e9) > 30).cumsum()).sort_values('zmad').groupby('g').tail(1).sort_values('epoch')
        return x
    ets, btw = events('ETS'), events('between')
    # T1'
    rate_b = C['ctrl_mean_between'] / PY['between']
    lam = max(rate_b * PY['ETS'], 1e-6)
    from scipy.stats import poisson
    out['T1'] = dict(n_ETS=len(ets), lam=lam, p=float(poisson.sf(len(ets) - 1, lam)) if len(ets) else 1.0,
                     ets_window_ctrl_mean=C['ctrl_mean_ETS'])
    out['T1']['pass'] = bool(len(ets) >= 5 and out['T1']['p'] < 0.01)
    # T2'
    R = __import__('s13_common').svi_rows()
    R = R[(R.series == 'data') & R.counted_once & (R.P >= 0.75)]
    F = K.families(); F = F[F.region == 'S_Vancouver_Is'].copy()
    F['family_id'] = pd.to_numeric(F.family_id, errors='coerce'); F = F[F.family_id.notna()]
    F['family_id'] = F.family_id.astype(int); F = F.set_index('family_id')        # svi_rows.fam is int
    LF = X.lfe_lists()
    n_any = X.cnt(LF['svi075'], ets.epoch.values, 600) + X.cnt(LF['wa_nw'], ets.epoch.values, 600) if len(ets) else np.array([])
    t2a = float((n_any >= 1).mean()) if len(ets) else np.nan
    rng = np.random.default_rng(P.SEED + 11)
    obs, null = 0, np.zeros(1000)
    for e in ets.itertuples():
        fam = F.index[K.dkm(F.lat.values, F.lon.values, e.lat_best, e.lon_best) <= 20]
        tt = np.sort(R[R.fam.isin(fam)].t.values)
        obs += int(X.cnt(tt, [e.epoch], 600)[0])
        day0 = np.floor(e.epoch / 86400) * 86400
        null += X.cnt(tt, day0 + rng.uniform(0, 86400, 1000), 600)
    p2b = float(np.mean(null >= obs)) if len(ets) else 1.0
    out['T2'] = dict(T2a=t2a, T2b_obs=int(obs), T2b_null_mean=float(null.mean()), T2b_null_p99=float(np.percentile(null, 99)),
                     T2b_p=p2b)
    out['T2']['pass'] = bool(np.isfinite(t2a) and t2a >= 0.7 and obs > np.percentile(null, 99))
    # T3'
    cc = SR.all_coords()
    top = ets.sort_values('zmad').tail(3)
    t3 = []
    for e in top.itertuples():
        r = t3_event(e, cc); r['time_pt'] = str(pd.Timestamp(e.epoch, unit='s').tz_localize('UTC').tz_convert('America/Los_Angeles'))
        r['zmad'] = e.zmad; t3.append(r)
    T3 = pd.DataFrame(t3); T3.to_csv(f'{K.TB}/pilot{Y}_R_T3{SFX}.csv', index=False)
    out['T3'] = dict(n_pass=int(T3.pass_.sum()) if len(T3) else 0, n=len(T3)); out['T3']['pass'] = out['T3']['n_pass'] >= 2
    # controls
    out['controls'] = dict(ctrl_rate_between_per_yr=C['ctrl_rate_between_per_yr'], noise_pair_rate_between_per_yr=C['noise_pair_rate_between_per_yr'])
    out['controls']['pass'] = bool(C['ctrl_rate_between_per_yr'] <= 1 and C['noise_pair_rate_between_per_yr'] <= 1)
    # completeness
    J = pd.read_parquet(f'{PD}/inj_list.parquet').merge(pd.read_parquet(f'{PD}/inj{SFX}.parquet'), on='inj_id', how='left')
    J['mf2'] = (J.mf_best_z >= P.LEVEL) & (J.mf_second_z >= P.LEVEL)
    J['bat'] = J.battery_e2.fillna(False).astype(bool)
    J['rec'] = J.mf2 & J.bat
    rows = []
    for (per, mw), g in J.groupby(['period', 'Mw']):
        k, n = int(g.rec.sum()), len(g); lo, hi = wilson(k, n)
        km = int(g.mf2.sum()); kb = int((g.mf2 & g.bat).sum())
        rows.append(dict(period=per, Mw=mw, n=n, recovered=k / n, ci_lo=lo, ci_hi=hi, mf_only=km / n,
                         battery_given_mf=kb / km if km else np.nan, battery_alone=float(g.bat.mean())))
    CP = pd.DataFrame(rows); CP.to_csv(f'{K.TB}/pilot{Y}_R_completeness{SFX}.csv', index=False)
    mwc = {}
    for per, g in CP.groupby('period'):
        ok = g[g.recovered >= 0.8]
        mwc[per] = float(ok.Mw.min()) if len(ok) else '>4.0'
    out['completeness'] = dict(Mw_c=mwc, at_3p8={r.period: round(r.recovered, 3) for r in CP[CP.Mw == 3.8].itertuples()})
    # between-ETS columns, Tier A'
    S = storm_table(); S.to_csv(f'{K.TB}/pilot{Y}_R_stormdays{SFX}.csv', index=False)
    sd = dict(zip(S.day, S.storm_day)); sp = dict(zip(S.day, S.pct))
    bt = btw.copy()
    if len(bt):
        bt['day'] = pd.to_datetime(bt.epoch, unit='s').dt.date.astype(str)
        bt['storm_day'] = bt.day.map(sd).fillna(False); bt['storm_pct'] = bt.day.map(sp)
        bt['cluster_day'] = bt.groupby('day').i.transform('count') >= 2
        q2 = []
        for e in bt.itertuples():
            A = amplitudes(e.epoch)
            nm = max([v for s, v in A.items() if s in NEAR] or [np.nan]); fm = max([v for s, v in A.items() if s in FAR] or [np.nan])
            c1 = e.vr_thrust >= 0.12; c2 = fm < 0.5 * nm; c3 = (e.ratio_GNW <= 0.551) if e.lat_best > 48.1 else True
            q2.append(dict(far_over_near=fm / nm, q2_real_looking=bool(c1 and c2 and c3)))
        bt = pd.concat([bt.reset_index(drop=True), pd.DataFrame(q2)], axis=1)
        bt['tierA'] = bt.q2_real_looking & ~(bt.storm_day & bt.cluster_day)
        bt['time_pt'] = pd.to_datetime(bt.epoch, unit='s').dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.strftime('%Y-%m-%d %H:%M:%S %Z')
        bt[['time_pt', 'template', 'agree', 'zmad', 'z_best', 'vr_thrust', 'ratio_far', 'ratio_GNW', 'far_over_near', 'lat_best',
            'lon_best', 'lfe_status', 'n_lfe_10m', 'storm_pct', 'storm_day', 'cluster_day', 'q2_real_looking', 'tierA']].to_csv(
            f'{K.TB}/pilot{Y}_R_between{SFX}.csv', index=False)
    out['between'] = dict(raw=len(bt), tierA=int(bt.tierA.sum()) if len(bt) else 0)
    out['GO_items'] = dict(T1=out['T1']['pass'], T2=out['T2']['pass'], T3=out['T3']['pass'], controls=out['controls']['pass'])
    out['GO'] = bool(all(out['GO_items'].values()))
    json.dump(out, open(f'{K.TB}/pilot{Y}_R_verdict{SFX}.json', 'w'), indent=1, default=str)
    pd.set_option('display.width', 250)
    print(json.dumps(out, indent=1, default=str)); print(T3.drop(columns=['amps', 'dists']).round(2).to_string())
    print(CP.round(3).to_string())
    if len(bt):
        print(bt[['time_pt', 'vr_thrust', 'ratio_GNW', 'far_over_near', 'storm_pct', 'cluster_day', 'q2_real_looking', 'tierA', 'lfe_status']].round(2).to_string())
    figs(CP)


def figs(CP):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7, 4.5))
    for per, c in (('ETS', 'tab:red'), ('between', 'k')):
        g = CP[CP.period == per]
        ax.plot(g.Mw, g.recovered, 'o-', color=c, label=f'{per}: tier')
        ax.fill_between(g.Mw, g.ci_lo, g.ci_hi, color=c, alpha=0.15, lw=0)
        ax.plot(g.Mw, g.mf_only, 's--', color=c, alpha=0.6, label=f'{per}: MF two templates')
        ax.plot(g.Mw, g.battery_given_mf, '^:', color=c, alpha=0.6, label=f'{per}: battery′ given MF')
    ax.axhline(0.8, color='k', lw=0.5); ax.set_ylim(0, 1.02); ax.set_xlabel('injected Mw'); ax.set_ylabel('fraction')
    ax.set_title(f'Pilot {Y} completeness (Amendment R; band = Wilson 95 %)', fontsize=9); ax.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{K.FIG}/pilot{Y}_R_completeness{SFX}.png', dpi=130)


if __name__ == '__main__':
    main()
