"""S15 PLAN Amendment 3 (Merlin audit): Q2 reporting columns, O1 cross-check, and the Poisson-regression gate
(hourly VLFE counts vs zone-LFE rate + M2 tide + days since local onset, episode offsets; same model for LFEs)."""
import os, json
os.environ.setdefault('OMP_NUM_THREADS', '4'); os.environ.setdefault('OPENBLAS_NUM_THREADS', '4')
import numpy as np, pandas as pd
import statsmodels.api as sm
from scipy.special import i0, i1
from scipy.optimize import brentq
import s15_common as C
import s15_q1 as Q
import s15_q2 as Q2

rng = np.random.default_rng(4015)
Q2.CONST = ['m2', 'k1', 'o1']
NB = 500


def A1inv(R):
    if R <= 0:
        return 0.0
    return brentq(lambda k: i1(k) / i0(k) - R, 1e-6, 50)


def circ_a(th):
    th = th[np.isfinite(th)]; z = np.mean(np.exp(1j * th)); return 2 * abs(z), float(np.degrees(np.angle(z))), len(th)


def columns(V, pools, pools_loc, Sb, Ss, F):
    """κ, peak/trough, LFE baseline without the VLFE's own ±1 h, a by Δt / burst count / Mw, O1."""
    E = pd.read_csv(f'{C.OUT_T}/q2_events.csv')
    assert np.allclose(E.epoch.values, V.epoch.values)
    th = np.radians(E.th_tau_updip_Pa_m2_deg.values)
    Z = np.load(f'{C.OUT_T}/q2_phases.npz'); okc = Z['okc']
    out = {}
    t = V.epoch.values
    first1h = np.r_[True, np.diff(t) > 3600] & okc
    for lab, m in (('all', okc), ('declustered_1h', first1h)):
        a, t0, n = circ_a(th[m]); k = A1inv(a / 2)
        out[f'kappa_{lab}'] = dict(n=n, a=a, theta0=t0, kappa=k, peak_over_trough=float(np.exp(2 * k)))
    # matched draws again, (a) as before, (b) excluding |dt| <= 1 h; plus O1 for VLFE / LFE / random
    Fi = F.set_index('key')
    with_lfe = V.n_lfe_10m.values >= 1
    n = len(V); ND = 1000
    Ct = np.full((n, ND), np.nan); Ce = np.full((n, ND), np.nan); Cf = np.zeros((n, ND), int); Cfe = np.zeros((n, ND), int)
    for i in range(n):
        p = pools[i] if with_lfe[i] else pools_loc[i]
        if p is None:
            continue
        d = np.abs(p[0] - t[i])
        for w in (1.0, 2.0, None):
            m = d <= w * C.DAY if w else np.ones(len(d), bool)
            if m.any():
                break
        j = rng.choice(np.where(m)[0], ND); Ct[i] = p[0][j]; Cf[i] = Fi.svi_fid.reindex(p[1][j]).values.astype(int)
        me = m & (d > 3600)
        if me.any():
            j = rng.choice(np.where(me)[0], ND); Ce[i] = p[0][j]; Cfe[i] = Fi.svi_fid.reindex(p[1][j]).values.astype(int)
    sf = E.stress_fid.values
    okx = okc & np.isfinite(Ce[:, 0])
    Uf = np.repeat(sf[:, None], ND, 1)
    qv = Q2.quantities(t, sf, Sb, Ss)
    qc = Q2.quantities(np.where(okc[:, None], Ct, np.nan), np.where(okc[:, None], Cf, Uf), Sb, Ss)
    qe = Q2.quantities(np.where(okx[:, None], Ce, np.nan), np.where(okx[:, None], Cfe, Uf), Sb, Ss)
    qu = Q2.quantities(t[:, None] + rng.uniform(-C.DAY, C.DAY, (n, ND)), Uf, Sb, Ss)
    for k_ in ('m2', 'o1'):
        key = f'th_tau_updip_Pa_{k_}'
        for lab, m in (('all', okc), ('declustered_1h', first1h)):
            av = circ_a(qv[key][m])
            ac = 2 * np.abs(np.nanmean(np.exp(1j * qc[key][m]), 0))
            ae = 2 * np.abs(np.nanmean(np.exp(1j * qe[key][m & okx]), 0))
            au = 2 * np.abs(np.nanmean(np.exp(1j * qu[key][m]), 0))
            out[f'{k_}_{lab}'] = dict(n=av[2], a_vlfe=av[0], theta0_vlfe=av[1], a_lfe_median=float(np.median(ac)),
                                      a_lfe_p01_p99=[float(np.percentile(ac, 1)), float(np.percentile(ac, 99))],
                                      a_lfe_excl_own_1h_median=float(np.median(ae)), a_lfe_excl_own_1h_p99=float(np.percentile(ae, 99)),
                                      a_rand_mean=float(au.mean()), a_rand_p99=float(np.percentile(au, 99)),
                                      pct_vs_lfe=float(100 * np.mean(ac < av[0])), pct_vs_rand=float(100 * np.mean(au < av[0])))
    # a by Δt, burst intensity, Mw
    Q1 = pd.read_csv(f'{C.OUT_T}/q1_events.csv')
    dt = Q1.dt_d.values; nb = V.n_lfe_10m.values; mw = V.Mw_scan.values
    rows = []
    for lab, m in [('dt<1', dt < 1), ('dt 1-3', (dt >= 1) & (dt < 3)), ('dt>=3', dt >= 3)]:
        rows.append(('days since local onset', lab) + circ_a(th[m & okc]))
    q = np.nanpercentile(nb[with_lfe], [33.33, 66.67])
    for lab, m in [('count tercile 1', with_lfe & (nb <= q[0])), ('count tercile 2', with_lfe & (nb > q[0]) & (nb <= q[1])),
                   ('count tercile 3', with_lfe & (nb > q[1])), ('no LFEs', ~with_lfe)]:
        rows.append(('±10-min LFE count', lab) + circ_a(th[m & okc]))
    med = np.nanmedian(mw)
    for lab, m in [(f'Mw_scan < {med:.2f}', mw < med), (f'Mw_scan >= {med:.2f}', mw >= med)]:
        rows.append(('size', lab) + circ_a(th[m & okc]))
    B = pd.DataFrame(rows, columns=['by', 'bin', 'a', 'theta0_deg', 'n'])
    B.to_csv(f'{C.OUT_T}/q2c_a_by_bins.csv', index=False)
    out['count_terciles_edges'] = q.tolist()
    return out, B


def hourly(V, F, L, Sb, Ss):
    E = C.episodes().set_index('episode')
    Q1 = pd.read_csv(f'{C.OUT_T}/q1_events.csv')
    with_lfe = V.n_lfe_10m.values >= 1
    ax = np.where(with_lfe, Q1.anchor_x, V.vx).astype(float); ay = np.where(with_lfe, Q1.anchor_y, V.vy).astype(float)
    cnt = V.episode.value_counts()
    eps = sorted(cnt[cnt >= 5].index)
    S115 = F[F.src == 'svi']
    tl = L.t.values; kl = L.key.values
    rows, info = [], []
    for e in eps:
        m = (V.episode == e).values
        cx, cy = np.median(ax[m]), np.median(ay[m])
        zone = set(F.key.values[np.hypot(F.x - cx, F.y - cy) <= 20])
        sfid = int(S115.fid.values[np.argmin(np.hypot(S115.x - cx, S115.y - cy))])
        t0, t1 = E.loc[e, 't0'], E.loc[e, 't1']
        h0 = np.floor(t0 / 3600) * 3600; edges = np.arange(h0, t1 + 3600, 3600)
        nv = np.histogram(V.epoch.values[m], edges)[0]
        inw = (tl >= t0) & (tl < t1)
        zm = inw & np.isin(kl, list(zone))
        nz = np.histogram(tl[zm], edges)[0]; na = np.histogram(tl[inw], edges)[0]
        onset = tl[zm].min() if zm.any() else np.nan
        mid = edges[:-1] + 1800
        ds = (mid - onset) / C.DAY
        q = Q2.quantities(mid, np.full(len(mid), sfid), Sb, Ss)
        flag = bool(V.flag_2016_17.values[m].any())
        rows.append(pd.DataFrame(dict(episode=e, t_mid=mid, n_vlfe=nv, n_zone=nz, n_all=na, days_since_onset=ds,
                                      th_m2=q['th_tau_updip_Pa_m2'], th_o1=q['th_tau_updip_Pa_o1'], flag_2016_17=flag)))
        info.append(dict(episode=e, n_vlfe=int(m.sum()), n_vlfe_in_window=int(nv.sum()), centroid_lat=C.LAT0 + cy / C.KM,
                         centroid_lon=C.LON0 + cx / (C.KM * np.cos(np.radians(C.LAT0))), n_zone_fam=len(zone),
                         n_zone_lfe=int(nz.sum()), stress_fid=sfid, onset_pt=C.to_pt([onset])[0] if np.isfinite(onset) else None, hours=len(mid)))
    H = pd.concat(rows, ignore_index=True)
    b = pd.cut(H.days_since_onset.fillna(-1), [-1e9, 0, 1, 3, 6, 1e9], right=False, labels=['pre', '0-1', '1-3', '3-6', '>=6'])
    H['onset_bin'] = b.astype(str)
    H['day_block'] = H.episode + '_' + (np.floor((H.t_mid - H.groupby('episode').t_mid.transform('min')) / C.DAY)).astype(int).astype(str)
    return H, pd.DataFrame(info)


def design(H, rate_col=None, const='m2'):
    X = pd.get_dummies(H.episode, prefix='ep', drop_first=False).astype(float)
    X = pd.concat([X, pd.get_dummies(H.onset_bin, prefix='on', drop_first=False).astype(float).drop(columns=['on_0-1'], errors='ignore')], axis=1)
    X[f'cos'] = np.cos(H[f'th_{const}']); X['sin'] = np.sin(H[f'th_{const}'])
    if rate_col:
        X['lograte'] = np.log1p(H[rate_col])
    return X


def fit(H, ycol, rate_col=None, const='m2'):
    X = design(H, rate_col, const)
    r = sm.GLM(H[ycol].values, X.values, family=sm.families.Poisson()).fit(maxiter=100)
    p = dict(zip(X.columns, r.params))
    A = float(np.hypot(p['cos'], p['sin'])); th0 = float(np.degrees(np.arctan2(p['sin'], p['cos'])))
    return dict(A=A, theta0=th0, beta=p.get('lograte', np.nan), c=p['cos'], s=p['sin'])


def gate(H, rate_col='n_zone', const='m2', label='headline', nb=NB):
    fv = fit(H, 'n_vlfe', rate_col, const); fl = fit(H, rate_col if rate_col else 'n_zone', None, const)
    blocks = H.day_block.unique(); idx = {b: np.where(H.day_block.values == b)[0] for b in blocks}
    BV, BL, BB = [], [], []
    for _ in range(nb):
        pick = rng.choice(blocks, len(blocks))
        Hb = H.iloc[np.concatenate([idx[b] for b in pick])]
        try:
            a = fit(Hb, 'n_vlfe', rate_col, const); b_ = fit(Hb, rate_col if rate_col else 'n_zone', None, const)
        except Exception:
            continue
        BV.append(a['A']); BL.append(b_['A']); BB.append(a['beta'])
    BV, BL, BB = map(np.array, (BV, BL, BB))
    ratio = BV / BL; diff = BV - (BB if rate_col else 0) * BL
    res = dict(label=label, const=const, rate_term=rate_col, hours=len(H), n_vlfe=int(H.n_vlfe.sum()),
               A_vlfe=fv['A'], theta0_vlfe=fv['theta0'], beta=fv['beta'], A_lfe=fl['A'], theta0_lfe=fl['theta0'],
               ratio=fv['A'] / fl['A'], ratio_ci=[float(np.percentile(ratio, 2.5)), float(np.percentile(ratio, 97.5))],
               A_vlfe_ci=[float(np.percentile(BV, 2.5)), float(np.percentile(BV, 97.5))],
               A_lfe_ci=[float(np.percentile(BL, 2.5)), float(np.percentile(BL, 97.5))],
               beta_times_A_lfe=(fv['beta'] if rate_col else np.nan) * fl['A'],
               diff_ci=[float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))], n_boot=len(BV))
    acc = res['ratio'] >= 2 and res['ratio_ci'][0] > 1
    fail = (res['ratio_ci'][0] <= 1) and (res['diff_ci'][0] <= 0 <= res['diff_ci'][1]) if rate_col else res['ratio_ci'][0] <= 1
    res['outcome'] = 'ACCEPT' if acc else ('FAILURE signature' if fail else 'intermediate')
    return res


def main():
    V, F, L, D, pools, pools_loc, finfo = Q.build()
    for c in ('epoch', 'vx', 'vy', 'Mw_scan', 'n_lfe_10m'):
        V[c] = V[c].astype(float)
    Sb, Ss = Q2.load_stress()
    cols, B = columns(V, pools, pools_loc, Sb, Ss, F)
    H, info = hourly(V, F, L, Sb, Ss)
    H.to_parquet(f'{C.OUT_T}/q2c_hourly.parquet'); info.to_csv(f'{C.OUT_T}/q2c_episodes.csv', index=False)
    G = [gate(H, 'n_zone', 'm2', 'headline: zone LFE rate term'),
         gate(H, None, 'm2', 'no LFE-rate term'),
         gate(H.assign(n_zone=H.n_all), 'n_zone', 'm2', 'all pool LFEs as rate term (LFE model on all pool LFEs)'),
         gate(H[H.n_zone >= 1], 'n_zone', 'm2', 'hours with >= 1 zone LFE'),
         gate(H[~H.flag_2016_17], 'n_zone', 'm2', 'without 2016-17 episodes'),
         gate(H, 'n_zone', 'o1', 'O1, zone LFE rate term')]
    G = pd.DataFrame(G); G.to_csv(f'{C.OUT_T}/q2c_gate.csv', index=False)
    json.dump(cols, open(f'{C.OUT_T}/q2c_columns.json', 'w'), indent=1, default=float)
    pd.set_option('display.width', 300)
    print(info.to_string()); print(G.round(3).to_string()); print(B.round(3).to_string()); print(json.dumps(cols, indent=1, default=float))


if __name__ == '__main__':
    main()
