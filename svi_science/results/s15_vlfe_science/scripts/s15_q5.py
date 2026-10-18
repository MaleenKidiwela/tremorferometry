"""S15 Q5 — VLFE productivity, RTRs and tidal gating vs migration rate (PLAN Amendment 4).
Local passage pace = Theil–Sen slope of family onset (10 % quantile of the family's episode LFE times, >= 5 LFEs) vs
along-strike, over families within ±10 km (>= 4 families), folded to the episode front direction (+ = forward).
Ranking variable 'pace_f' = max(folded pace, 0) d/km: small = fast (or simultaneous / locally reversed), large = slow.
Daily front speed = Theil–Sen (km/d, folded) of the daily median along-strike of families whose onset falls that day,
over ±1.5 d."""
import json
import numpy as np, pandas as pd
from scipy.stats import theilslopes
import s15_common as C
import s15_q1 as Q
import s15_q2 as Q2

rng = np.random.default_rng(5015)
ND = 1000


def ts(x, y):
    if len(x) < 2 or np.ptp(x) == 0:
        return np.nan
    return theilslopes(y, x)[0]


def episode_kinematics(L, F, E, D, eps):
    al = F.set_index('key').along
    out, grids, daily = {}, {}, {}
    for e in eps:
        r = E.loc[e]
        Le = L[(L.t >= r.t0) & (L.t < r.t1)]
        g = Le.groupby('key').t.agg(['size', lambda x: np.quantile(x, 0.1)])
        g.columns = ['n', 'onset']; g = g[g.n >= 5]
        g['along'] = g.index.map(al); g = g.dropna()
        sgn = 1.0 if D.get(e) == 'NW-ward' else -1.0
        # episode main speed (km/d) along vs onset, and pace
        v_main = ts(g.onset.values / C.DAY, g.along.values)
        grid = np.arange(-80, 81, 1.0); pace = np.full(len(grid), np.nan); nf = np.zeros(len(grid), int)
        for j, a in enumerate(grid):
            m = np.abs(g.along.values - a) <= 10
            nf[j] = m.sum()
            if m.sum() >= 4:
                pace[j] = ts(g.along.values[m], g.onset.values[m] / C.DAY) * sgn
        grids[e] = (grid, pace, nf)
        # daily front speed
        day = np.floor((g.onset.values - r.t0) / C.DAY).astype(int)
        dd = pd.DataFrame(dict(day=day, along=g.along.values)).groupby('day').along.median()
        dsp = {}
        for d0 in range(int(np.floor((r.t1 - r.t0) / C.DAY)) + 1):
            m = (dd.index >= d0 - 1.5) & (dd.index <= d0 + 1.5)
            dsp[d0] = ts(dd.index.values[m].astype(float), dd.values[m]) * sgn if m.sum() >= 3 else np.nan
        daily[e] = dsp
        out[e] = dict(episode=e, direction=D.get(e), n_fam_onset=len(g), v_main_kmd=v_main * sgn if np.isfinite(v_main) else np.nan,
                      pace_grid_median=float(np.nanmedian(pace)))
    return pd.DataFrame(out.values()), grids, daily


def pace_at(grids, e, a):
    grid, pace, _ = grids[e]
    a = np.asarray(a, float)
    j = np.clip(np.round(a - grid[0]).astype(int), 0, len(grid) - 1)
    p = pace[j]; p[~np.isfinite(a)] = np.nan
    return p


def daily_at(daily, E, e, t):
    d = np.floor((np.asarray(t, float) - E.loc[e, 't0']) / C.DAY).astype(int)
    return np.array([daily[e].get(int(x), np.nan) for x in d])


def tercile(x, edges):
    x = np.maximum(np.asarray(x, float), 0)
    out = np.where(x <= edges[0], 0, np.where(x <= edges[1], 1, 2)).astype(float)
    out[~np.isfinite(np.asarray(x, float))] = np.nan
    return out                                                           # 0 fast, 1 medium, 2 slow


def main():
    V, F, L, D_, pools, pools_loc, finfo = Q.build()
    for c in ('epoch', 'anchor_along', 'n_lfe_10m', 'Mw_scan'):
        V[c] = V[c].astype(float)
    E = C.episodes().set_index('episode')
    D = D_.set_index('episode').direction.to_dict()
    A = pd.read_parquet(f'{C.OUT_T}/q1b_pseudo_anchor_lfe.parquet')         # per-LFE pseudo-anchor Δt, x_ex (along), eligible
    eps = sorted(V.episode[V.n_lfe_10m >= 1].unique())
    K, grids, daily = episode_kinematics(L, F, E, D, eps)
    s3 = pd.read_csv(f'{C.RES}/s3/episodes_kinematics.csv'); s3 = s3[s3.region == 'SVI'].set_index('episode')
    K['s3_lfe_speed'] = K.episode.map(s3.lfe_speed)
    K.to_csv(f'{C.OUT_T}/q5_episode_speeds.csv', index=False)
    with_lfe = V.n_lfe_10m.values >= 1
    V['pace_local'] = np.nan; V['daily_speed'] = np.nan
    for e in eps:
        m = (V.episode == e).values & with_lfe
        V.loc[m, 'pace_local'] = pace_at(grids, e, V.anchor_along.values[m])
        V.loc[m, 'daily_speed'] = daily_at(daily, E, e, V.epoch.values[m])
    # pool-LFE distribution (all pool LFEs of the analysed episodes, at their own family along) -> tercile edges
    al = F.set_index('key').along
    Lp = []
    for e in eps:
        r = E.loc[e]; Le = L[(L.t >= r.t0) & (L.t < r.t1)].copy()
        Le['episode'] = e; Le['pace_local'] = pace_at(grids, e, Le.key.map(al).values); Le['daily_speed'] = daily_at(daily, E, e, Le.t.values)
        Lp.append(Le)
    Lp = pd.concat(Lp, ignore_index=True)
    edges = np.nanpercentile(np.maximum(Lp.pace_local.values, 0), [33.333, 66.667])
    dedges = np.nanpercentile(Lp.daily_speed.values, [33.333, 66.667])
    V['pace_terc'] = tercile(V.pace_local, edges)
    V['daily_terc'] = np.where(V.daily_speed <= dedges[0], 0, np.where(V.daily_speed <= dedges[1], 1, 2)).astype(float)
    V.loc[V.daily_speed.isna(), 'daily_terc'] = np.nan
    # matched re-anchored draws (±24 h, same concurrent families, eligible pseudo-VLFEs)
    Ag = {e: g for e, g in A.groupby('episode')}
    n = len(V)
    Pp = np.full((n, ND), np.nan); Pt = np.full((n, ND), np.nan); Px = np.full((n, ND), np.nan)
    Pw = np.full((n, ND), np.nan)   # whole-episode variant: any eligible pool LFE of the episode
    for i in np.where(with_lfe)[0]:
        e = V.episode.values[i]; g = Ag[e]
        ks = set(pools[i][1])
        p = g[g.key.isin(ks) & g.eligible & (np.abs(g.t - V.epoch.values[i]) <= C.DAY)]
        if len(p) == 0:
            continue
        j = rng.integers(0, len(p), ND)
        Pt[i] = p.t.values[j]; Px[i] = p.x_ex.values[j]
        Pp[i] = pace_at(grids, e, Px[i])
        q = g[g.eligible]
        j2 = rng.integers(0, len(q), ND)
        Pw[i] = pace_at(grids, e, q.x_ex.values[j2])
    ok = with_lfe & np.isfinite(Pt[:, 0]) & np.isfinite(V.pace_local.values)
    res = dict(pace_edges_d_per_km=edges.tolist(), speed_edges_kmd=[float(1 / x) if x > 0 else np.inf for x in edges],
               daily_speed_edges_kmd=dedges.tolist(), n_used=int(ok.sum()), n_with_lfe=int(with_lfe.sum()),
               frac_vlfe_pace_le0=float(np.mean(V.pace_local.values[ok] <= 0)), frac_lfe_pace_le0=float(np.nanmean(Lp.pace_local.values <= 0)))

    def contrast(tobs, tdraw, rows):
        fo = np.array([np.mean(tobs[rows] == k) for k in range(3)])
        fd = np.array([[np.nanmean(tdraw[rows, d] == k) for k in range(3)] for d in range(tdraw.shape[1])])   # (ND,3)
        fm = fd.mean(0)
        ratio = fo / fm
        Robs = ratio[2] / ratio[0]
        Rd = (fd[:, 2] / fm[2]) / (fd[:, 0] / fm[0])
        return dict(n=int(rows.sum()), frac_obs=fo.tolist(), frac_ctrl=fm.tolist(), per_lfe_ratio=ratio.tolist(),
                    slow_over_fast=float(Robs), ctrl_p01=float(np.nanpercentile(Rd, 1)), ctrl_p99=float(np.nanpercentile(Rd, 99)),
                    pct=float(100 * np.nanmean(Rd < Robs)))
    Tm = np.full_like(Pp, np.nan); Tm[ok] = tercile(Pp[ok].ravel(), edges).reshape(Pp[ok].shape)
    Tw = np.full_like(Pw, np.nan); Tw[ok] = tercile(Pw[ok].ravel(), edges).reshape(Pw[ok].shape)
    to = V.pace_terc.values
    res['main_matched_24h'] = contrast(to, Tm, ok)
    res['variant_whole_episode'] = contrast(to, Tw, ok)
    per = []
    for e in C.BIG9:
        m = ok & (V.episode == e).values
        if m.sum() >= 4:
            c = contrast(to, Tm, m); per.append(dict(episode=e, direction=D.get(e), **{k: c[k] for k in ('n', 'slow_over_fast', 'pct')},
                                                      frac_obs=';'.join(f'{x:.2f}' for x in c['frac_obs']), frac_ctrl=';'.join(f'{x:.2f}' for x in c['frac_ctrl'])))
    PE = pd.DataFrame(per); PE.to_csv(f'{C.OUT_T}/q5_per_episode.csv', index=False)
    for d in ('NW-ward', 'SE-ward'):
        m = ok & (V.direction == d).values
        res[f'dir_{d}'] = contrast(to, Tm, m)
    s = np.sign(np.log(res['main_matched_24h']['slow_over_fast']))
    n_same = int((np.sign(np.log(PE.slow_over_fast.replace(0, np.nan))) == s).sum())
    dir_same = all(np.sign(np.log(res[f'dir_{d}']['slow_over_fast'])) == s for d in ('NW-ward', 'SE-ward'))
    R = res['main_matched_24h']['slow_over_fast']
    big = (R >= 1.5 and res['main_matched_24h']['pct'] > 99) or (R <= 1 / 1.5 and res['main_matched_24h']['pct'] < 1)
    res['verdict'] = dict(slow_over_fast=R, pct=res['main_matched_24h']['pct'], big9_same_sign=n_same, n_big9=len(PE),
                          both_directions=bool(dir_same),
                          outcome=('more VLFEs per LFE where slow' if R > 1 else 'more VLFEs per LFE where fast') if (big and n_same >= 6 and dir_same) else 'flat (rule not met)')
    # Mw by tercile (column)
    res['median_Mw_by_tercile'] = [float(np.nanmedian(V.Mw_scan.values[ok & (to == k)])) for k in range(3)]
    # RTR vs ±24 h matched draws, and fast-forward (top daily-speed tercile)
    Rt = pd.read_csv(f'{C.RES}/q8_rtr/tables/rtr_list_frozen_d2.csv')
    Rt = Rt[(Rt.cat == 'bostock') & Rt.ep.isin(['2010-08', '2011-08', '2012-08', '2013-09'])]
    rs = pd.to_datetime(Rt.start).values.astype('datetime64[ns]').astype('int64') / 1e9
    re_ = pd.to_datetime(Rt.end).values.astype('datetime64[ns]').astype('int64') / 1e9
    xlo = np.minimum(Rt.x_start, Rt.x_end).values - 20; xhi = np.maximum(Rt.x_start, Rt.x_end).values + 20

    def inside(t, x):
        t = np.asarray(t, float)[..., None]; x = np.asarray(x, float)[..., None]
        return np.any((t >= rs) & (t <= re_) & (x >= xlo) & (x <= xhi), -1)
    m = ok & (V.year.values <= 2013)
    io = inside(V.epoch.values[m], V.anchor_along.values[m]); idr = inside(Pt[m], Px[m]).mean(0)
    res['rtr_vs_matched_24h'] = dict(n=int(m.sum()), frac_obs=float(io.mean()), ctrl_mean=float(idr.mean()),
                                     ctrl_p99=float(np.percentile(idr, 99)), pct=float(100 * np.mean(idr < io.mean())))
    m2 = ok & np.isfinite(V.daily_speed.values)
    dsd = np.full_like(Pt, np.nan)
    for i in np.where(m2)[0]:
        dsd[i] = daily_at(daily, E, V.episode.values[i], Pt[i])
    ff_o = V.daily_speed.values[m2] > dedges[1]; ff_d = np.nanmean(dsd[m2] > dedges[1], 0)
    res['fast_forward_vs_matched_24h'] = dict(n=int(m2.sum()), frac_obs=float(ff_o.mean()), ctrl_mean=float(ff_d.mean()),
                                              ctrl_p01=float(np.percentile(ff_d, 1)), ctrl_p99=float(np.percentile(ff_d, 99)),
                                              pct=float(100 * np.mean(ff_d < ff_o.mean())))
    to_d = V.daily_terc.values; Td = np.full_like(dsd, np.nan)
    Td[m2] = np.where(dsd[m2] <= dedges[0], 0, np.where(dsd[m2] <= dedges[1], 1, 2))
    Td[~np.isfinite(dsd)] = np.nan
    cd = contrast(np.where(np.isfinite(to_d), 2 - to_d, np.nan), np.where(np.isfinite(Td), 2 - Td, np.nan), m2)   # reorder so index 2 = slow daily speed
    res['daily_speed_tercile_matched_24h'] = cd
    # tidal a by local-pace tercile (declustered, matched LFE draws from q2_phases)
    Z = np.load(f'{C.OUT_T}/q2_phases.npz'); okc = Z['okc']
    first = np.r_[True, np.diff(V.epoch.values) > 3600]
    tid = []
    for k, nm in enumerate(['fast', 'medium', 'slow']):
        for lab, msk in (('declustered_1h', first), ('all', np.ones(n, bool))):
            mm = okc & ok & (to == k) & msk
            th = Z['vlfe'][mm]; a = 2 * abs(np.mean(np.exp(1j * th)))
            ac = 2 * np.abs(np.nanmean(np.exp(1j * Z['lfe'][mm]), 0))
            tid.append(dict(tercile=nm, set=lab, n=int(mm.sum()), a_vlfe=a, theta0=float(np.degrees(np.angle(np.mean(np.exp(1j * th))))),
                            a_lfe_median=float(np.median(ac)), a_lfe_p99=float(np.percentile(ac, 99))))
    TI = pd.DataFrame(tid); TI.to_csv(f'{C.OUT_T}/q5_tide_by_tercile.csv', index=False)
    # Poisson: daily-speed tercile x tide interaction (hourly table from q2c)
    import statsmodels.api as sm
    H = pd.read_parquet(f'{C.OUT_T}/q2c_hourly.parquet')
    H['dspeed'] = np.nan
    for e in H.episode.unique():
        if e in daily:
            m = (H.episode == e).values
            H.loc[m, 'dspeed'] = daily_at(daily, E, e, H.t_mid.values[m])
    H = H[np.isfinite(H.dspeed)].copy()
    H['dterc'] = np.where(H.dspeed <= dedges[0], 'fast0', np.where(H.dspeed <= dedges[1], 'med1', 'slow2'))
    H['dterc'] = np.where(H.dspeed <= dedges[0], 'slowest', np.where(H.dspeed <= dedges[1], 'medium', 'fastest'))
    pr = []
    for ycol, rate in (('n_vlfe', True), ('n_zone', False)):
        X = pd.get_dummies(H.episode, prefix='ep').astype(float)
        X = pd.concat([X, pd.get_dummies(H.onset_bin, prefix='on').astype(float).drop(columns=['on_0-1'], errors='ignore')], axis=1)
        dt = pd.get_dummies(H.dterc, prefix='sp').astype(float).drop(columns=['sp_medium'])
        X = pd.concat([X, dt], axis=1)
        for k in ('slowest', 'medium', 'fastest'):
            ind = (H.dterc == k).astype(float).values
            X[f'cos_{k}'] = np.cos(H.th_m2) * ind; X[f'sin_{k}'] = np.sin(H.th_m2) * ind
        if rate:
            X['lograte'] = np.log1p(H.n_zone)
        r = sm.GLM(H[ycol].values, X.values, family=sm.families.Poisson()).fit()
        p = dict(zip(X.columns, r.params)); cov = pd.DataFrame(r.cov_params(), index=X.columns, columns=X.columns)
        for k in ('slowest', 'medium', 'fastest'):
            c, s_ = p[f'cos_{k}'], p[f'sin_{k}']
            pr.append(dict(model=ycol, daily_speed_tercile=k, hours=int((H.dterc == k).sum()), events=int(H[ycol][H.dterc == k].sum()),
                           A=float(np.hypot(c, s_)), theta0=float(np.degrees(np.arctan2(s_, c))),
                           A_se=float(np.sqrt((c * c * cov.loc[f'cos_{k}', f'cos_{k}'] + s_ * s_ * cov.loc[f'sin_{k}', f'sin_{k}']
                                               + 2 * c * s_ * cov.loc[f'cos_{k}', f'sin_{k}'])) / max(np.hypot(c, s_), 1e-9))))
    PR = pd.DataFrame(pr); PR.to_csv(f'{C.OUT_T}/q5_poisson_speed_tide.csv', index=False)
    # detection vs speed: C0/C1 injections
    J = pd.read_csv(f'{C.OUT_T}/detect_injections.csv')
    J = J[J.episode.isin(eps)].copy()
    SRC = {'C0': (48.141757352279875, -123.52741640459564), 'C1': (48.29464879935246, -123.61541420640553)}
    coef = finfo['along_fit_coef']
    J['along'] = [coef[0] + coef[1] * C.xy(*SRC[s])[0] + coef[2] * C.xy(*SRC[s])[1] for s in J.source]
    J['pace_local'] = np.nan; J['daily_speed'] = np.nan
    for e in J.episode.unique():
        m = (J.episode == e).values
        J.loc[m, 'pace_local'] = pace_at(grids, e, J.along.values[m]); J.loc[m, 'daily_speed'] = daily_at(daily, E, e, J.epoch.values[m])
    J['pace_terc'] = tercile(J.pace_local, edges)
    J['daily_terc'] = np.where(J.daily_speed <= dedges[0], 0, np.where(J.daily_speed <= dedges[1], 1, 2)).astype(float)
    J.loc[J.daily_speed.isna(), 'daily_terc'] = np.nan
    det = []
    for col, names in (('pace_terc', ['fast', 'medium', 'slow']), ('daily_terc', ['slow daily', 'medium daily', 'fast daily'])):
        for mw in (3.6, 3.8):
            s_ = J[J.Mw == mw]
            rec = [s_.rec[s_[col] == k].mean() for k in range(3)]; nn = [int((s_[col] == k).sum()) for k in range(3)]
            rel = (rec[2] - rec[0]) / rec[0] if rec[0] else np.nan
            det.append(dict(by=col, Mw=mw, **{f'rec_{names[k]}': rec[k] for k in range(3)}, **{f'n_{names[k]}': nn[k] for k in range(3)},
                            rel_diff_outer=rel, flag_30pct=bool(abs(rel) > 0.3)))
    DT = pd.DataFrame(det); DT.to_csv(f'{C.OUT_T}/q5_detection_by_speed.csv', index=False)
    V[['time_pt', 'epoch', 'episode', 'direction', 'n_lfe_10m', 'anchor_along', 'pace_local', 'pace_terc', 'daily_speed', 'daily_terc', 'Mw_scan']].to_csv(
        f'{C.OUT_T}/q5_events.csv', index=False)
    np.savez_compressed(f'{C.OUT_T}/q5_draws.npz', pace_draw=Pp.astype(np.float32), ok=ok, lfe_pace=Lp.pace_local.values.astype(np.float32),
                        lfe_daily=Lp.daily_speed.values.astype(np.float32), edges=edges, dedges=dedges)
    json.dump(res, open(f'{C.OUT_T}/q5_summary.json', 'w'), indent=1, default=float)
    pd.set_option('display.width', 250)
    print(K.round(3).to_string()); print(json.dumps(res, indent=1, default=float)); print(PE.round(3).to_string())
    print(TI.round(3).to_string()); print(PR.round(3).to_string()); print(DT.round(3).to_string())


if __name__ == '__main__':
    main()
