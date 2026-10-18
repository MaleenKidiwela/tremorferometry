"""S15 Q3 — slip partition per episode (PLAN Amendment 4): LFE self-moment, continuous (S13), discrete (catalogue),
GNSS; raw and corrected side by side; size test on the matched-completeness subset."""
import json, glob
import numpy as np, pandas as pd
from scipy.stats import spearmanr
import s15_common as C

EPS = ['2010-08b1', '2011-08', '2012-08', '2013-09', '2018-06', '2019-08', '2019-11', '2020-11', '2023-04', '2024-08', '2025-09b1',
       '2015-12', '2021-01', '2021-10', '2022-03']
MATCHED = ['2011-08', '2012-08', '2013-09', '2018-06', '2019-08', '2020-11', '2024-08']
FLAGGED = ['2023-04', '2025-09b1', '2015-12']
M0 = lambda mw: 10 ** (1.5 * np.asarray(mw, float) + 9.1)
C0C1_MID = ((48.141757 + 48.294649) / 2, (-123.527416 + -123.615414) / 2)


def compl_table():
    rows = []
    for r in C.RUNS:
        e = pd.read_csv(f'{C.TB}/{r}_episodes_v2.csv'); e['run'] = r; rows.append(e)
    E = pd.concat(rows)
    E['span'] = (pd.to_datetime(E.end) - pd.to_datetime(E.start)).dt.days + 1
    g = {}
    for ep, d in E.groupby('episode'):
        w = d.span.values.astype(float)
        g[ep] = {c: float(np.average(d[c].values, weights=w)) for c in ('compl_3.4', 'compl_3.6', 'compl_3.8')}
    return g


def main():
    Q1 = pd.read_csv(f'{C.OUT_T}/q1_events.csv')
    F, _ = C.families(); L = C.lfes(); E = C.episodes().set_index('episode')
    fam = pd.read_csv(f'{C.RES}/s0/families.csv').set_index('fid')
    B = pd.read_csv(f'{C.ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv').set_index('bkey')
    mag = {f'S:{f}': fam.mag_med.get(f, np.nan) for f in fam.index}
    for f in B.index:
        if f'S:{f}' not in mag or not np.isfinite(mag[f'S:{f}']):
            mag[f'S:{f}'] = B.mag_med.get(f, np.nan)
    L['mag'] = L.key.map(mag)
    S13 = pd.read_csv(f'{C.RES}/s13_vlfe/tables/scale_svi_episodes.csv').set_index('episode')
    G = pd.read_csv(f'{C.RES}/s1/episode_geodetic.csv').set_index('episode')
    SL = pd.read_parquet(f'{C.ROOT}/lfe_extension/results/slip_compare/slip_svi_family_daily.parquet', columns=['fid', 'date', 'slip_mm'])
    SL['date'] = pd.to_datetime(SL.date)
    reg = SL.groupby('date').slip_mm.mean()
    comp = compl_table()
    relay = pd.read_csv(f'{C.TB}/prod2018_relay_events.csv').epoch.values
    is_relay = np.array([np.any(np.abs(relay - t) <= 30) for t in Q1.epoch.values])
    # n_sta_day per LFE (backward + forward) for the size-test control
    import s13_common as S
    A = pd.read_parquet(f'{C.RES}/s0/analysis_table.parquet', columns=['run', 'series', 'P_same', 'keep_xfam8', 'time_utc', 'n_sta_day'])
    A = A[(A.run.astype(str) == 'prod') & A.keep_xfam8 & (A.series.astype(str) == 'data') & (A.P_same >= 0.75)]
    At = A.time_utc.astype('datetime64[ns]').astype('int64').values / 1e9
    Bk = pd.read_parquet(f'{C.ROOT}/data/catalogue/svi_catalogue_v1_4_backward.parquet', columns=['series', 't', 'P', 'primary', 'keep_xfam8', 'counted_once', 'n_sta_day'])
    Bk = Bk[Bk.primary & (Bk.keep_xfam8 == 1) & (Bk.series.astype(str) == 'data') & Bk.counted_once & (Bk.P >= 0.75)]
    Bt = Bk.t.values.astype(float)                                  # epoch seconds (trap: not a datetime)
    nsd_t = np.r_[Bt, At]; nsd_v = np.r_[Bk.n_sta_day.values, A.n_sta_day.values].astype(float)
    with_lfe = Q1.n_lfe_10m.values >= 1
    ax = np.where(with_lfe, Q1.anchor_x, Q1.vx).astype(float); ay = np.where(with_lfe, Q1.anchor_y, Q1.vy).astype(float)
    rows = []
    for e in EPS:
        t0, t1 = E.loc[e, 't0'], E.loc[e, 't1']
        m = (Q1.episode == e).values
        if m.sum():
            cx, cy = np.median(ax[m]), np.median(ay[m]); cen = 'VLFE anchors'
        else:
            cx, cy = C.xy(*C0C1_MID); cx, cy = float(cx), float(cy); cen = 'C0/C1 midpoint (no VLFEs)'
        zone = set(F.key.values[np.hypot(F.x - cx, F.y - cy) <= 20])
        Le = L[(L.t >= t0) & (L.t < t1)]
        svi = Le.key.str.startswith('S:')
        zs = Le.key.isin(zone) & svi
        o = dict(episode=e, centroid=cen, n_vlfe=int(m.sum()), n_vlfe_with_lfe=int((m & with_lfe).sum()),
                 n_lfe_svi=int(svi.sum()), n_lfe_zone_svi=int(zs.sum()), n_lfe_zone_all=int(Le.key.isin(zone).sum()),
                 lfe_self_m0_svi=float(M0(Le.mag[svi]).sum()), lfe_self_m0_zone=float(M0(Le.mag[zs]).sum()),
                 lfe_self_m0_per_lfe=float(M0(Le.mag[svi]).mean()) if svi.any() else np.nan)
        if e in S13.index:
            r = S13.loc[e]
            o.update(s13_n=int(r.n), thrust_m0=float(r.thrust_m0), thrust_z=float(r.thrust_z), cont_m0=float(r.m0_total),
                     cont_m0_zone=float(r.thrust_m0) * o['n_lfe_zone_svi'])
        mw = Q1.Mw_scan.values[m]; rl = is_relay[m]
        c = comp.get(e, {})
        o.update(compl_3_4=c.get('compl_3.4'), compl_3_6=c.get('compl_3.6'), compl_3_8=c.get('compl_3.8'),
                 disc_m0_raw=float(M0(mw).sum()), disc_m0_no_relay=float(M0(mw[~rl]).sum()), n_relay_self=int(rl.sum()))
        o['disc_m0_div_c34'] = o['disc_m0_no_relay'] / c['compl_3.4'] if c.get('compl_3.4') else np.nan
        b1 = (mw >= 3.5) & (mw < 3.7) & ~rl; b2 = (mw >= 3.7) & ~rl
        o['disc_m0_ge35_bincorr'] = (float(M0(mw[b1]).sum()) / c['compl_3.6'] + float(M0(mw[b2]).sum()) / c['compl_3.8']) if c.get('compl_3.6') else np.nan
        o['disc_frac_m0_ge35'] = float(M0(mw[b1 | b2]).sum() / max(M0(mw[~rl]).sum(), 1))
        if e in G.index:
            g = G.loc[e]
            o.update(gnss_D=float(g.D), gnss_SN=float(g.SN), gnss_tierB_M0=float(g.tierB_M0_Nm))
        d0 = pd.Timestamp(E.loc[e, 'start']) - pd.Timedelta(days=2); d1 = pd.Timestamp(E.loc[e, 'end']) + pd.Timedelta(days=2)
        rr = reg[(reg.index >= d0) & (reg.index <= d1)]
        o['reg_slip_mm'] = float(rr.iloc[-1] - rr.iloc[0]) if len(rr) > 1 else np.nan
        o['reg_M0'] = o['reg_slip_mm'] * 1e-3 * 30e9 * (100e3 * 60e3) if np.isfinite(o['reg_slip_mm']) else np.nan   # nominal 100x60 km SVI area, mu 30 GPa (scale only)
        k = (nsd_t >= t0) & (nsd_t < t1)
        o['n_sta_day_median'] = float(np.median(nsd_v[k])) if k.any() else np.nan
        rows.append(o)
    P = pd.DataFrame(rows)
    P['vlfe_per_100_zone_lfe_raw'] = 100 * (P.n_vlfe - P.n_relay_self) / P.n_lfe_zone_svi.replace(0, np.nan)
    P['vlfe_per_100_zone_lfe_corr'] = P.vlfe_per_100_zone_lfe_raw / P.compl_3_4
    P['disc_over_cont_raw'] = P.disc_m0_raw / P.cont_m0
    P['disc_over_cont_corr'] = P.disc_m0_div_c34 / P.cont_m0
    P['disc_over_cont_zone_corr'] = P.disc_m0_div_c34 / P.cont_m0_zone
    P['cont_over_lfe_self'] = P.cont_m0 / P.lfe_self_m0_svi
    P['disc_over_lfe_self_corr'] = P.disc_m0_div_c34 / P.lfe_self_m0_svi
    P['lfe_self_over_tierB'] = P.lfe_self_m0_svi / P.gnss_tierB_M0
    P['cont_over_tierB'] = P.cont_m0 / P.gnss_tierB_M0
    P['disc_corr_over_tierB'] = P.disc_m0_div_c34 / P.gnss_tierB_M0
    P['seismic_slow_over_tierB'] = (P.lfe_self_m0_svi + P.cont_m0 + P.disc_m0_div_c34) / P.gnss_tierB_M0
    P['subset'] = np.where(P.episode.isin(MATCHED), 'matched', np.where(P.episode.isin(FLAGGED), 'flagged', 'other'))
    P.to_csv(f'{C.OUT_T}/q3_partition.csv', index=False)
    # size test
    M = P[P.episode.isin(MATCHED)]
    res = {}
    for ycol in ('vlfe_per_100_zone_lfe_corr', 'vlfe_per_100_zone_lfe_raw', 'disc_over_cont_corr'):
        o = {}
        for xcol in ('reg_slip_mm', 'gnss_D', 'compl_3_4', 'n_sta_day_median'):
            d = M[[ycol, xcol]].dropna()
            if len(d) >= 4:
                s = spearmanr(d[ycol], d[xcol]); o[xcol] = dict(n=len(d), rho=float(s.statistic), p=float(s.pvalue))
        y = M[ycol].dropna(); o['max_over_min'] = float(y.max() / y.min()) if len(y) and y.min() > 0 else np.nan
        res[ycol] = o
    head = res['vlfe_per_100_zone_lfe_corr']
    comp_ok = all(abs(head.get(k, {}).get('rho', 0)) < 0.6 for k in ('compl_3_4', 'n_sta_day_median'))
    size = head.get('reg_slip_mm', {}).get('rho', np.nan)
    if not comp_ok:
        v = 'uninterpretable (productivity tracks completeness or station count, |rho| >= 0.6)'
    elif abs(size) >= 0.79 and head['max_over_min'] >= 3:
        v = 'size effect (rho %.2f, factor %.1f)' % (size, head['max_over_min'])
    else:
        v = 'flat / not claimable (rho %.2f, factor %.1f)' % (size, head['max_over_min'])
    res['verdict'] = v
    json.dump(res, open(f'{C.OUT_T}/q3_size_test.json', 'w'), indent=1, default=float)
    pd.set_option('display.width', 320)
    cols = ['episode', 'n_vlfe', 'n_lfe_svi', 'n_lfe_zone_svi', 'lfe_self_m0_svi', 'cont_m0', 'disc_m0_raw', 'disc_m0_div_c34', 'disc_m0_ge35_bincorr',
            'compl_3_4', 'gnss_D', 'gnss_tierB_M0', 'reg_slip_mm', 'n_sta_day_median', 'vlfe_per_100_zone_lfe_raw', 'vlfe_per_100_zone_lfe_corr',
            'disc_over_cont_raw', 'disc_over_cont_corr', 'cont_over_lfe_self', 'seismic_slow_over_tierB']
    print(P[cols].to_string(float_format=lambda x: f'{x:.3g}'))
    print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
