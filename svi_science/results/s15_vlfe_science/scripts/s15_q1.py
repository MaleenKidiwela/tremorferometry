"""S15 Q1 — VLFE timing against the local LFE sequence (day / minutes / hours scales). Rules: ../PLAN.md §2–3."""
import json
import numpy as np, pandas as pd
import s15_common as C

NDRAW = 1000
rng = np.random.default_rng(15)
DAY = C.DAY


def build():
    V = C.vlfes(); F, finfo = C.families(); L = C.lfes(); E = C.episodes()
    D = C.directions(L, F, E[E.episode.isin(V.episode.unique())])
    Fi = F.set_index('key')
    byk = {k: g.t.values for k, g in L.groupby('key')}
    tl = L.t.values; kl = L.key.values
    rows, pools, pools_loc = [], [], []
    for _, v in V.iterrows():
        t, a, b = v.epoch, v.t0_ep, v.t1_ep
        o = dict(epoch=t)
        # concurrent families (+-600 s)
        lo, hi = np.searchsorted(tl, t - 600), np.searchsorted(tl, t + 600, side='right')
        ks, cn = np.unique(kl[lo:hi], return_counts=True)
        o['n_conc_fam'] = len(ks); o['n_conc_lfe'] = int(cn.sum())
        if len(ks):
            first, p5, ins, keys = [], [], [], []
            for k in ks:
                x = byk[k]; x = x[(x >= a) & (x < b)]
                first.append(x.min()); p5.append(np.percentile(x, 5)); ins.append(x); keys.append(np.full(len(x), k, object))
            first, p5, w = np.array(first), np.array(p5), cn.astype(float)
            Fw = np.sum(w * first) / w.sum()
            o.update(onset_w=Fw, dt_d=(t - Fw) / DAY, dt_min_d=(t - first.max()) / DAY,
                     dt_median_d=np.median(t - first) / DAY, dt_earliest_d=(t - first.min()) / DAY,
                     dt_p5_d=(t - np.sum(w * p5) / w.sum()) / DAY,
                     first_on_window_day1=bool(np.any(first < a + DAY)),
                     anchor_x=np.sum(w * Fi.loc[ks].x.values) / w.sum(), anchor_y=np.sum(w * Fi.loc[ks].y.values) / w.sum(),
                     anchor_along=np.sum(w * Fi.loc[ks].along.values) / w.sum())
            pt, pk = np.concatenate(ins), np.concatenate(keys)
            o['n_pool'] = len(pt); o['u_rank'] = float(np.mean(pt < t))
            pools.append((pt, pk, Fw))
        else:
            pools.append(None)
        # VLFE-location variant: families within 20 km (30 km fallback) with >= 1 LFE in the episode
        d = np.hypot(F.x.values - v.vx, F.y.values - v.vy)
        pl = None
        for R in (20.0, 30.0):
            near = F.key.values[d <= R]
            ins, keys, first, w = [], [], [], []
            for k in near:
                if k not in byk:
                    continue
                x = byk[k]; x = x[(x >= a) & (x < b)]
                if len(x):
                    ins.append(x); keys.append(np.full(len(x), k, object)); first.append(x.min()); w.append(len(x))
            if ins:
                first, w = np.array(first), np.array(w, float)
                Fw2 = np.sum(w * first) / w.sum()
                pt, pk = np.concatenate(ins), np.concatenate(keys)
                o.update(loc_R_km=R, loc_n_fam=len(first), loc_n_pool=len(pt), loc_onset_w=Fw2,
                         loc_dt_d=(t - Fw2) / DAY, loc_u_rank=float(np.mean(pt < t)))
                pl = (pt, pk, Fw2)
                break
        pools_loc.append(pl)
        rows.append(o)
    P = pd.DataFrame(rows)
    V = pd.concat([V.reset_index(drop=True), P.drop(columns='epoch')], axis=1)
    V['direction'] = V.episode.map(D.set_index('episode').direction)
    V['vlfe_along'] = np.c_[np.ones(len(V)), V.vx, V.vy] @ np.array(finfo['along_fit_coef'])
    return V, F, L, D, pools, pools_loc, finfo


def draws(pools):
    """NDRAW draw indices per event -> times and keys (None rows stay NaN)."""
    n = len(pools)
    T = np.full((n, NDRAW), np.nan); K = np.full((n, NDRAW), None, object)
    for i, p in enumerate(pools):
        if p is None:
            continue
        j = rng.integers(0, len(p[0]), NDRAW)
        T[i] = p[0][j]; K[i] = p[1][j]
    return T, K


def pct(obs, ctrl):
    return float(100.0 * np.mean(ctrl < obs) + 50.0 * np.mean(ctrl == obs))


def day_scale(V, obs_col, onset_col, T, mask, label):
    """median Δt observed vs 1,000 control medians on events in mask; per-episode and per-direction tables."""
    m = mask & V[obs_col].notna().values
    obs = V.loc[m, obs_col].values
    ctrl = (T[m] - V.loc[m, onset_col].values[:, None]) / DAY          # (n, NDRAW)
    cm = np.median(ctrl, 0)
    out = dict(label=label, n=int(m.sum()), obs_median_d=float(np.median(obs)), ctrl_median_of_medians_d=float(np.median(cm)),
               ctrl_p01=float(np.percentile(cm, 1)), ctrl_p99=float(np.percentile(cm, 99)), obs_pct=pct(np.median(obs), cm),
               obs_mean_d=float(obs.mean()), ctrl_mean_d=float(ctrl.mean()))
    out['diff_d'] = out['obs_median_d'] - out['ctrl_median_of_medians_d']
    ep = []
    Vm = V[m].reset_index(drop=True)
    for e, g in Vm.groupby('episode'):
        idx = g.index.values
        c = np.median(ctrl[idx], 0)
        ep.append(dict(label=label, episode=e, direction=g.direction.iloc[0], n=len(idx), obs_median_d=float(np.median(obs[idx])),
                       ctrl_median_d=float(np.median(c)), ctrl_p01=float(np.percentile(c, 1)), ctrl_p99=float(np.percentile(c, 99)),
                       obs_pct=pct(np.median(obs[idx]), c), mean_u=float(g.get('u_rank', pd.Series(np.nan)).mean())
                       if 'u_rank' in g else np.nan))
    ep = pd.DataFrame(ep); ep['diff_d'] = ep.obs_median_d - ep.ctrl_median_d
    dr = []
    for d_, g in Vm.groupby('direction'):
        idx = g.index.values
        c = np.median(ctrl[idx], 0)
        dr.append(dict(label=label, direction=d_, n=len(idx), episodes=';'.join(sorted(g.episode.unique())),
                       obs_median_d=float(np.median(obs[idx])), ctrl_median_d=float(np.median(c)),
                       ctrl_p01=float(np.percentile(c, 1)), ctrl_p99=float(np.percentile(c, 99)), obs_pct=pct(np.median(obs[idx]), c)))
    dr = pd.DataFrame(dr); dr['diff_d'] = dr.obs_median_d - dr.ctrl_median_d
    return out, ep, dr, ctrl


def verdict(out, ep, dr):
    big = ep[ep.episode.isin(C.BIG9)]
    s = np.sign(out['diff_d'])
    n_same = int((np.sign(big.diff_d) == s).sum())
    dir_same = {r.direction: bool(np.sign(r.diff_d) == s) for r in dr.itertuples()}
    late = out['obs_pct'] > 99 and s > 0 and n_same >= 6 and all(dir_same.get(k, False) for k in ('SE-ward', 'NW-ward'))
    early = out['obs_pct'] < 1 and s < 0 and n_same >= 6 and all(dir_same.get(k, False) for k in ('SE-ward', 'NW-ward'))
    return dict(verdict='PASS-late' if late else ('PASS-early' if early else 'NULL (rule)'), pooled_sign=int(s),
                big9_present=int(len(big)), big9_same_sign=n_same, direction_same_sign=dir_same, obs_pct=out['obs_pct'])


def minutes(V, L, T, mask):
    tl = L.t.values
    edges = np.array([-600, -120, 0, 120, 600.0])
    fine = np.arange(-600, 601, 20.0)

    def counts(times, excl_self=False):
        c = np.zeros((len(times), 4)); h = np.zeros(len(fine) - 1)
        for i, t in enumerate(times):
            lo, hi = np.searchsorted(tl, t - 600), np.searchsorted(tl, t + 600, side='right')
            d = tl[lo:hi] - t
            if excl_self:
                d = d[np.abs(d) > 10]
            c[i] = np.histogram(d, edges)[0]
            h += np.histogram(d, fine)[0]
        return c, h

    t_obs = V.epoch.values[mask]
    c_obs, h_obs = counts(t_obs)
    # control: 200 of the 1,000 draws (LFE-clustering shape), drawn LFE and anything within 10 s excluded
    Tm = T[mask]
    cc, hh = [], np.zeros(len(fine) - 1)
    for k in range(200):
        c, h = counts(Tm[:, k], excl_self=True); cc.append(c.sum(0)); hh += h
    cc = np.array(cc); hh /= 200
    from scipy.stats import binomtest
    S = c_obs.sum(0)
    res = dict(n_events=int(mask.sum()), counts=S.tolist(), per_min=(S / np.array([8, 2, 2, 8])).tolist(),
               inner_after_over_before=float(S[2] / max(S[1], 1)), outer_after_over_before=float(S[3] / max(S[0], 1)),
               inner_binom_p=float(binomtest(int(S[2]), int(S[1] + S[2]), 0.5).pvalue) if S[1] + S[2] else np.nan,
               outer_binom_p=float(binomtest(int(S[3]), int(S[0] + S[3]), 0.5).pvalue) if S[0] + S[3] else np.nan,
               ctrl_counts_mean=cc.mean(0).tolist(),
               ctrl_inner_ratio_mean=float(np.mean(cc[:, 2] / np.maximum(cc[:, 1], 1))),
               ctrl_inner_ratio_p01_p99=[float(np.percentile(cc[:, 2] / np.maximum(cc[:, 1], 1), q)) for q in (1, 99)],
               ctrl_outer_ratio_mean=float(np.mean(cc[:, 3] / np.maximum(cc[:, 0], 1))),
               ctrl_outer_ratio_p01_p99=[float(np.percentile(cc[:, 3] / np.maximum(cc[:, 0], 1), q)) for q in (1, 99)],
               ctrl_peak_over_outer_mean=float(np.mean((cc[:, 1] + cc[:, 2]) / 4 / np.maximum((cc[:, 0] + cc[:, 3]) / 16, 1e-9))),
               obs_peak_over_outer=float((S[1] + S[2]) / 4 / max((S[0] + S[3]) / 16, 1e-9)),
               frac_events_more_after_inner=float(np.mean(c_obs[:, 2] > c_obs[:, 1])),
               frac_events_more_before_inner=float(np.mean(c_obs[:, 2] < c_obs[:, 1])))
    return res, c_obs, h_obs, hh, fine


def rtr(V, F, T, K, Tl, Kl):
    R = pd.read_csv(f'{C.RES}/q8_rtr/tables/rtr_list_frozen_d2.csv')
    R = R[(R.cat == 'bostock') & R.ep.isin(['2010-08', '2011-08', '2012-08', '2013-09'])]
    assert len(R) == 30, len(R)
    rs = pd.to_datetime(R.start).values.astype('datetime64[ns]').astype('int64') / 1e9
    re_ = pd.to_datetime(R.end).values.astype('datetime64[ns]').astype('int64') / 1e9
    xlo = np.minimum(R.x_start, R.x_end).values - 20; xhi = np.maximum(R.x_start, R.x_end).values + 20
    al = F.set_index('key').along

    def inside(t, x):
        t = np.asarray(t, float)[..., None]; x = np.asarray(x, float)[..., None]
        return np.any((t >= rs) & (t <= re_) & (x >= xlo) & (x <= xhi), -1)

    m = V.year.between(2010, 2013).values
    x_obs = np.where(V.anchor_along.notna(), V.anchor_along, V.vlfe_along)
    obs = inside(V.epoch.values[m], x_obs[m])
    obs_loc = inside(V.epoch.values[m], V.vlfe_along.values[m])
    # control: concurrent pool if present, location pool otherwise; x = drawn family's along
    Tc = np.where(np.isfinite(T), T, Tl); Kc = np.where(np.isfinite(T), K, Kl)
    Tc, Kc = Tc[m], Kc[m]
    xc = np.vectorize(lambda k: al.get(k, np.nan) if k is not None else np.nan)(Kc)
    ok = np.isfinite(Tc[:, 0])
    ins = inside(Tc[ok], xc[ok])                                 # (n_ok, NDRAW)
    f_ctrl = ins.mean(0)
    # time-only control variant: drawn time with the VLFE's own x
    ins_t = inside(Tc[ok], np.repeat(x_obs[m][ok][:, None], NDRAW, 1))
    res = dict(n_vlfe=int(m.sum()), n_with_pool=int(ok.sum()), frac_obs=float(obs[ok].mean()), n_inside=int(obs[ok].sum()),
               frac_obs_vlfe_location_x=float(obs_loc[ok].mean()),
               ctrl_mean=float(f_ctrl.mean()), ctrl_p01=float(np.percentile(f_ctrl, 1)), ctrl_p99=float(np.percentile(f_ctrl, 99)),
               obs_pct=pct(obs[ok].mean(), f_ctrl),
               ctrl_timeonly_mean=float(ins_t.mean()), ctrl_timeonly_p01_p99=[float(np.percentile(ins_t.mean(0), q)) for q in (1, 99)],
               obs_pct_timeonly=pct(obs[ok].mean(), ins_t.mean(0)))
    per_ep = pd.DataFrame(dict(episode=V.episode.values[m][ok], inside=obs[ok], ctrl=ins.mean(1)))
    res['per_episode'] = per_ep.groupby('episode').agg(n=('inside', 'size'), frac_obs=('inside', 'mean'),
                                                       frac_ctrl=('ctrl', 'mean')).reset_index().to_dict('records')
    return res, obs, m


def s3_returns(V, F, T, K, Tl, Kl):
    B = pd.read_csv(f'{C.RES}/s3/bursts.csv')
    B = B[B.region == 'SVI']
    days = np.array(B.return_day_PT.astype(str).tolist(), dtype='U10'); bc = B.bin_center_km.values.astype(float)
    al = F.set_index('key').along

    def inside(t, x):
        d = pd.to_datetime(np.asarray(t, float).ravel(), unit='s', utc=True).tz_convert('America/Los_Angeles').strftime('%Y-%m-%d')
        d = np.array(list(d), dtype='U10').reshape(np.shape(t))
        return np.any((d[..., None] == days) & (np.abs(np.asarray(x)[..., None] - bc) <= 10), -1)

    x_obs = np.where(V.anchor_along.notna(), V.anchor_along, V.vlfe_along)
    m = V.episode.isin(B.episode.unique()).values
    obs = inside(V.epoch.values[m], x_obs[m])
    Tc = np.where(np.isfinite(T), T, Tl)[m]; Kc = np.where(np.isfinite(T), K, Kl)[m]
    ok = np.isfinite(Tc[:, 0])
    xc = np.vectorize(lambda k: al.get(k, np.nan) if k is not None else np.nan)(Kc[ok][:, :200])
    ins = inside(Tc[ok][:, :200], xc)
    f = ins.mean(0)
    res = dict(n_vlfe=int(m.sum()), n_with_pool=int(ok.sum()), episodes=sorted(V.episode[m].unique().tolist()),
               frac_obs=float(obs[ok].mean()), n_inside=int(obs[ok].sum()), ctrl_mean=float(f.mean()),
               ctrl_p01_p99=[float(np.percentile(f, q)) for q in (1, 99)], obs_pct=pct(obs[ok].mean(), f), n_draws=200)
    per = pd.DataFrame(dict(episode=V.episode.values[m][ok], inside=obs[ok], ctrl=ins.mean(1)))
    res['per_episode'] = per.groupby('episode').agg(n=('inside', 'size'), frac_obs=('inside', 'mean'), frac_ctrl=('ctrl', 'mean')).reset_index().to_dict('records')
    return res


def main():
    V, F, L, D, pools, pools_loc, finfo = build()
    T, K = draws(pools); Tl, Kl = draws(pools_loc)
    with_lfe = V.n_lfe_10m.values >= 1
    assert with_lfe.sum() == 431
    S, EP, DR = [], [], []
    res = {}
    for lab, col, on, TT, mk in [('A_concurrent_431', 'dt_d', 'onset_w', T, with_lfe),
                                 ('B_location_all', 'loc_dt_d', 'loc_onset_w', Tl, np.ones(len(V), bool)),
                                 ('B_location_431', 'loc_dt_d', 'loc_onset_w', Tl, with_lfe),
                                 ('C_location_109', 'loc_dt_d', 'loc_onset_w', Tl, ~with_lfe),
                                 ('A_concurrent_431_no2016_17', 'dt_d', 'onset_w', T, with_lfe & ~V.flag_2016_17.values)]:
        out, ep, dr, ctrl = day_scale(V, col, on, TT, mk, lab)
        out.update(verdict(out, ep, dr)); S.append(out); EP.append(ep); DR.append(dr)
        if lab == 'A_concurrent_431':
            np.save(f'{C.OUT_T}/q1_ctrl_dt_concurrent.npy', ctrl.astype(np.float32))
        if lab == 'C_location_109':
            np.save(f'{C.OUT_T}/q1_ctrl_dt_loc109.npy', ctrl.astype(np.float32))
    # robust-onset column variant (p5 onset; control with the same p5 onset)
    m = with_lfe
    on_p5 = V.epoch.values - V.dt_p5_d.values * DAY
    ctrl = (T[m] - on_p5[m][:, None]) / DAY
    cm = np.median(ctrl, 0)
    res['robust_p5_onset'] = dict(obs_median_d=float(np.nanmedian(V.dt_p5_d[m])), ctrl_median_d=float(np.median(cm)),
                                  obs_pct=pct(np.nanmedian(V.dt_p5_d[m]), cm))
    res['u_rank'] = dict(mean_u_431=float(V.u_rank[m].mean()), median_u_431=float(V.u_rank[m].median()),
                         mean_loc_u_109=float(V.loc_u_rank[~m].mean()), n_loc_109=int(V.loc_u_rank[~m].notna().sum()),
                         frac_u_gt_half=float((V.u_rank[m] > 0.5).mean()))
    pd.DataFrame(S).to_csv(f'{C.OUT_T}/q1_day_summary.csv', index=False)
    pd.concat(EP).to_csv(f'{C.OUT_T}/q1_day_per_episode.csv', index=False)
    pd.concat(DR).to_csv(f'{C.OUT_T}/q1_day_per_direction.csv', index=False)
    D.to_csv(f'{C.OUT_T}/q1_episode_directions.csv', index=False)
    mres, c_obs, h_obs, h_ctrl, fine = minutes(V, L, T, with_lfe)
    mres_all, _, h_all, _, _ = minutes(V, L, np.where(np.isfinite(T), T, Tl), np.isfinite(np.where(np.isfinite(T), T, Tl)[:, 0]))
    res['minutes_431'] = mres; res['minutes_all_with_pool'] = mres_all
    np.savez(f'{C.OUT_T}/q1_minutes_hist.npz', fine=fine, h_obs=h_obs, h_ctrl=h_ctrl, c_obs=c_obs)
    V[['m10_2', 'm2_0', 'p0_2', 'p2_10']] = np.nan
    V.loc[with_lfe, ['m10_2', 'm2_0', 'p0_2', 'p2_10']] = c_obs
    rres, rin, rm = rtr(V, F, T, K, Tl, Kl)
    res['rtr_2010_13'] = rres
    V['in_rtr'] = np.nan; V.loc[rm, 'in_rtr'] = rin.astype(float)
    res['s3_returns'] = s3_returns(V, F, T, K, Tl, Kl)
    res['families_info'] = finfo
    V['time_pt'] = C.to_pt(V.epoch.values)
    V.to_csv(f'{C.OUT_T}/q1_events.csv', index=False)
    np.save(f'{C.OUT_T}/q1_draw_times_concurrent.npy', T.astype(np.float64))
    np.save(f'{C.OUT_T}/q1_draw_times_location.npy', Tl.astype(np.float64))
    json.dump(res, open(f'{C.OUT_T}/q1_other_scales.json', 'w'), indent=1, default=float)
    pd.set_option('display.width', 250)
    print(pd.DataFrame(S).round(3).to_string())
    print(pd.concat(EP)[lambda x: x.episode.isin(C.BIG9)].round(3).to_string())
    print(pd.concat(DR).round(3).to_string())
    print(json.dumps({k: v for k, v in res.items() if k != 'families_info'}, indent=1, default=float)[:6000])


if __name__ == '__main__':
    main()
