"""S15 Q1, PLAN Amendment 1 — re-anchored control: pseudo-VLFEs at random LFE times of the same families/episode, with
their own concurrent set, weights, first-LFE times and Δt recomputed exactly as for a real VLFE."""
import json
import numpy as np, pandas as pd
import s15_common as C
import s15_q1 as Q

DAY = C.DAY
rng = np.random.default_rng(1515)


def per_lfe_anchor(L, F, E, episodes):
    """for every pool LFE inside the listed episodes: Δt if it were the anchor (self excluded / included),
    eligibility (>= 1 other LFE within +-600 s), anchor along-strike."""
    al = F.set_index('key').along
    out = []
    for e in episodes:
        r = E[E.episode == e].iloc[0]
        m = (L.t >= r.t0) & (L.t < r.t1)
        Le = L[m].reset_index(drop=True)
        t = Le.t.values; k = Le.key.values
        first = Le.groupby('key').t.min()
        fk = first.reindex(k).values
        alk = al.reindex(k).values
        lo = np.searchsorted(t, t - 600); hi = np.searchsorted(t, t + 600, side='right')
        dt_ex = np.full(len(t), np.nan); dt_in = np.full(len(t), np.nan); x_ex = np.full(len(t), np.nan)
        nfam = np.zeros(len(t), int)
        for j in range(len(t)):
            idx = np.r_[lo[j]:j, j + 1:hi[j]]
            if len(idx):
                w = np.ones(len(idx))
                dt_ex[j] = (t[j] - np.sum(fk[idx]) / len(idx)) / DAY        # count-weighted over neighbours = mean of first_f per LFE
                x_ex[j] = np.nanmean(alk[idx])
                nfam[j] = len(np.unique(k[idx]))
            idx2 = np.arange(lo[j], hi[j])
            dt_in[j] = (t[j] - np.mean(fk[idx2])) / DAY
        out.append(pd.DataFrame(dict(episode=e, t=t, key=k, dt_ex=dt_ex, dt_in=dt_in, x_ex=x_ex, nfam_ex=nfam,
                                     eligible=np.isfinite(dt_ex))))
    return pd.concat(out, ignore_index=True)


def main():
    V, F, L, D, pools, pools_loc, finfo = Q.build()
    E = C.episodes()
    with_lfe = V.n_lfe_10m.values >= 1
    eps = sorted(V.episode[with_lfe].unique())
    A = per_lfe_anchor(L, F, E, eps)
    # sanity: the same formula applied to the real VLFEs reproduces dt_d (count-weighted mean of first_f over LFEs)
    res = {}
    n = len(V)
    Tex = np.full((n, Q.NDRAW), np.nan); Tin = np.full((n, Q.NDRAW), np.nan); Xex = np.full((n, Q.NDRAW), np.nan)
    Tdraw = np.full((n, Q.NDRAW), np.nan)
    Ag = {e: g for e, g in A.groupby('episode')}
    for i in np.where(with_lfe)[0]:
        v = V.iloc[i]
        g = Ag[v.episode]
        ks = set(pools[i][1])
        p = g[g.key.isin(ks)]
        pe = p[p.eligible]
        if len(pe) == 0:
            continue
        j = rng.integers(0, len(pe), Q.NDRAW)
        Tex[i] = pe.dt_ex.values[j]; Xex[i] = pe.x_ex.values[j]; Tdraw[i] = pe.t.values[j]
        j2 = rng.integers(0, len(p), Q.NDRAW)
        Tin[i] = p.dt_in.values[j2]
    okR = with_lfe & np.isfinite(Tex[:, 0])
    res['n_events_no_eligible_pseudo'] = int((with_lfe & ~okR).sum())
    S, EP, DR = [], [], []
    for lab, M in [('R_reanchored_excl_self', Tex), ('R_reanchored_incl_self', Tin)]:
        # day_scale expects draw TIMES and an onset column; pass Δt directly via onset 0
        V2 = V.copy(); V2['zero'] = 0.0
        out, ep, dr, ctrl = Q.day_scale(V2, 'dt_d', 'zero', M * DAY, okR, lab)
        out.update(Q.verdict(out, ep, dr)); S.append(out); EP.append(ep); DR.append(dr)
        if lab == 'R_reanchored_excl_self':
            np.save(f'{C.OUT_T}/q1_ctrl_dt_reanchored.npy', ctrl.astype(np.float32))
    # single-concurrent-family events under the fixed-weight control (no intersection effect)
    T = np.load(f'{C.OUT_T}/q1_draw_times_concurrent.npy')
    one = with_lfe & (V.n_conc_fam.values == 1)
    out, ep, dr, _ = Q.day_scale(V, 'dt_d', 'onset_w', T, one, 'A_fixed_single_family_only')
    out.update(Q.verdict(out, ep, dr)); S.append(out); EP.append(ep); DR.append(dr)
    out, ep, dr, _ = Q.day_scale(V.assign(zero=0.0), 'dt_d', 'zero', Tex * DAY, one & okR, 'R_reanchored_single_family_only')
    out.update(Q.verdict(out, ep, dr)); S.append(out); EP.append(ep); DR.append(dr)
    # u-rank against the re-anchored pseudo-VLFE times of the same pool
    u = np.full(n, np.nan)
    for i in np.where(okR)[0]:
        u[i] = np.mean(Tdraw[i] < V.epoch.values[i]) + 0.5 * np.mean(Tdraw[i] == V.epoch.values[i])
    V['u_rank_reanchored'] = u
    res['u_rank_reanchored'] = dict(mean=float(np.nanmean(u)), median=float(np.nanmedian(u)), frac_gt_half=float(np.nanmean(u[with_lfe] > 0.5)))
    res['pseudo_vlfe_mean_conc_fam'] = float(A[A.eligible].nfam_ex.mean())
    res['vlfe_mean_conc_fam'] = float(V.n_conc_fam[with_lfe].mean())
    res['n_single_family_events'] = int(one.sum())
    # RTR against re-anchored pseudo-VLFEs (time + pseudo anchor x)
    R = pd.read_csv(f'{C.RES}/q8_rtr/tables/rtr_list_frozen_d2.csv')
    R = R[(R.cat == 'bostock') & R.ep.isin(['2010-08', '2011-08', '2012-08', '2013-09'])]
    rs = pd.to_datetime(R.start).values.astype('datetime64[ns]').astype('int64') / 1e9
    re_ = pd.to_datetime(R.end).values.astype('datetime64[ns]').astype('int64') / 1e9
    xlo = np.minimum(R.x_start, R.x_end).values - 20; xhi = np.maximum(R.x_start, R.x_end).values + 20

    def inside(t, x):
        t = np.asarray(t, float)[..., None]; x = np.asarray(x, float)[..., None]
        return np.any((t >= rs) & (t <= re_) & (x >= xlo) & (x <= xhi), -1)
    m = V.year.between(2010, 2013).values & okR
    obs = inside(V.epoch.values[m], V.anchor_along.values[m])
    ins = inside(Tdraw[m], Xex[m]); f = ins.mean(0)
    res['rtr_2010_13_reanchored'] = dict(n=int(m.sum()), frac_obs=float(obs.mean()), n_inside=int(obs.sum()), ctrl_mean=float(f.mean()),
                                         ctrl_p01=float(np.percentile(f, 1)), ctrl_p99=float(np.percentile(f, 99)), obs_pct=Q.pct(obs.mean(), f),
                                         per_episode=pd.DataFrame(dict(ep=V.episode.values[m], o=obs, c=ins.mean(1))).groupby('ep')
                                         .agg(n=('o', 'size'), frac_obs=('o', 'mean'), frac_ctrl=('c', 'mean')).reset_index().to_dict('records'))
    pd.DataFrame(S).to_csv(f'{C.OUT_T}/q1b_day_summary.csv', index=False)
    pd.concat(EP).to_csv(f'{C.OUT_T}/q1b_day_per_episode.csv', index=False)
    pd.concat(DR).to_csv(f'{C.OUT_T}/q1b_day_per_direction.csv', index=False)
    A.to_parquet(f'{C.OUT_T}/q1b_pseudo_anchor_lfe.parquet')
    np.save(f'{C.OUT_T}/q1b_draw_times_reanchored.npy', Tdraw)
    V[['epoch', 'u_rank_reanchored']].to_csv(f'{C.OUT_T}/q1b_u_rank.csv', index=False)
    json.dump(res, open(f'{C.OUT_T}/q1b_other.json', 'w'), indent=1, default=float)
    pd.set_option('display.width', 250)
    print(pd.DataFrame(S).round(3).to_string())
    x = pd.concat(EP); print(x[x.episode.isin(C.BIG9)].round(3).to_string())
    print(pd.concat(DR).round(3).to_string())
    print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
