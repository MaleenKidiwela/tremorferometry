#!/usr/bin/env python
"""S13 Amendment 5 summary (2026-10-14 PT): WA corrected run (wa2_*), SVI step 4 (pbin_*), step 5 (svibe2), step 6
(fixed n_nb for SVI from the row files), Q1c PGC-out and the single-channel Z table (refits of the stored guarded SVI
stacks). Reads only; writes tables/v2_*.csv, tables/v2_verdict.json, figures/s13_v2_*.png.

usage: s13_v2_summary.py [wa] [svi]   (default both)"""
import os, sys, glob, json, pickle
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr, chi2 as CH
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
from s13_stack import FITM, LAGS, fit
from s13_stack2_v2 import nb_counts, vetoes

RES = f'{C.ROOT}/svi_science/results'
RNG = np.random.default_rng(C.SEED)
GREY = '#9a9a9a'
KS = [k for k in range(-5, 6) if k != 0]
VETO_LOG = []


def perm_spearman(x, y, n=5000):
    m = np.isfinite(x) & np.isfinite(y)
    x, y = np.asarray(x)[m], np.asarray(y)[m]
    if len(x) < 5:
        return dict(n=int(len(x)), rho=np.nan, p=np.nan)
    r = spearmanr(x, y).statistic
    rr = np.array([spearmanr(x, RNG.permutation(y)).statistic for _ in range(n)])
    return dict(n=int(len(x)), rho=float(r), p=float((1 + (np.abs(rr) >= abs(r)).sum()) / (n + 1)))


def wmean(m, s):
    m, s = np.asarray(m, float), np.asarray(s, float)
    k = np.isfinite(m) & np.isfinite(s) & (s > 0)
    w = 1 / s[k] ** 2
    mu = (m[k] * w).sum() / w.sum()
    chi = ((m[k] - mu) ** 2 * w).sum(); dof = max(k.sum() - 1, 1)
    return dict(n=int(k.sum()), wmean=float(mu), sd=float(1 / np.sqrt(w.sum())), chi2_dof=float(chi / dof),
                p_const=float(CH.sf(chi, dof)))


# ------------------------------------------------------------------------------------------------------------- WA
def wa():
    E = pd.read_csv(f'{C.TB}/wa2_episodes.csv')
    U = []
    for e in E[E.built].episode:
        f = f'{C.TB}/wa2_{e}_units.csv'
        if os.path.isfile(f):
            U.append(pd.read_csv(f).assign(episode=e))
    U = pd.concat(U, ignore_index=True)
    U.to_csv(f'{C.TB}/v2_wa_units_all.csv', index=False)
    for k in KS:
        U[f'zk{k:+d}'] = (U[f'thrust_m0_k{k:+d}'] - U.thrust_ctrl_mean) / U.thrust_ctrl_sd
    Pq = U[U.unit == 'pooled'].set_index('episode')
    v = {'n_built': int(E.built.sum()), 'n_stacked': int(len(Pq))}
    # W1 lag
    lags = Pq.thrust_lag
    vc = lags.value_counts()
    v['W1_lag_counts'] = {int(k): int(x) for k, x in vc.sort_index().items()}
    v['W1_lag_mode'] = int(vc.idxmax()); v['W1_frac_at_pm8'] = float((lags.abs() == 8).mean())
    v['W1_pass'] = bool(abs(v['W1_lag_mode']) <= 4 and v['W1_frac_at_pm8'] <= 0.10)
    for g in ('c6_', 'wide_'):
        lv = Pq[f'thrust_{g}lag'].value_counts()
        v[f'col_{g}lag_mode'] = int(lv.idxmax()); v[f'col_{g}lag_median'] = float(Pq[f'thrust_{g}lag'].median())
    # W2 duplicates
    v['W2_nb2_mean'] = float(Pq.n_nb2.mean()); v['W2_pass'] = bool(v['W2_nb2_mean'] < 0.05)
    v['nb_pooled_mean'] = {c: float(Pq[c].mean()) for c in ('n_nb2', 'n_nb10', 'n_nb20', 'n_nb60', 'nb_kern')}
    # W3 calibration
    zk = Pq[[f'zk{k:+d}' for k in KS]].values.ravel(); zk = zk[np.isfinite(zk)]
    sdk = float(np.std(zk)); v['W3_zk_sd'] = sdk; v['W3_zk_mean'] = float(np.mean(zk)); v['W3_n'] = int(len(zk))
    v['W3_zk_absge3'] = int((np.abs(zk) >= 3).sum())
    zka = U[[f'zk{k:+d}' for k in KS]].values.ravel(); zka = zka[np.isfinite(zka)]
    v['W3_allunits_zk_sd'] = float(np.std(zka))
    v['W3_pass'] = bool(0.85 <= sdk <= 1.15)
    scale = 1.0 if v['W3_pass'] else sdk
    v['sigma_scale_applied'] = scale
    U['sd_s'] = U.thrust_ctrl_sd * scale; U['z_s'] = (U.thrust_m0 - U.thrust_ctrl_mean) / U.sd_s
    Pq = U[U.unit == 'pooled'].set_index('episode')
    # W4 zones
    big = Pq.index[Pq.n >= 300]
    Z = U[U.unit.str.startswith('zone')]
    zt = Z.pivot_table(index='episode', columns='unit', values='thrust_m0')
    zz = Z.pivot_table(index='episode', columns='unit', values='z_s')
    v['W4_n_big'] = int(len(big))
    for z in ('zone1', 'zone2', 'zone3'):
        x = zt.reindex(big)[z] if z in zt else pd.Series(dtype=float)
        v[f'W4_{z}_pos'] = f'{int((x > 0).sum())}/{int(x.notna().sum())}'
        Zu = Z[Z.unit == z]
        v[f'{z}_wmean'] = wmean(Zu.thrust_m0, Zu.sd_s)
        v[f'{z}_nz3'] = int((Zu.z_s >= 3).sum()); v[f'{z}_nzm3'] = int((Zu.z_s <= -3).sum())
        v[f'{z}_lag_mode'] = int(Zu.thrust_lag.value_counts().idxmax()) if len(Zu) else None
    ok = []
    for z in ('zone2', 'zone3'):
        x = zt.reindex(big)[z].dropna() if z in zt else pd.Series(dtype=float)
        ok.append(len(x) > 0 and (x > 0).sum() > len(x) / 2)
    v['W4_pass'] = bool(all(ok))
    # W5 halves
    zo = U[U.unit == 'pooled_odd'].set_index('episode').z_s; ze = U[U.unit == 'pooled_even'].set_index('episode').z_s
    j = zo.index.intersection(ze.index)
    r = spearmanr(zo[j], ze[j]).statistic
    v['W5_rho_odd_even'] = float(r); v['W5_n'] = int(len(j)); v['W5_pass'] = bool(r > 0.4)
    det = Pq[Pq.z_s >= 3].index
    v['halves_both_z2_in_detected'] = f'{int(((zo.reindex(det) >= 2) & (ze.reindex(det) >= 2)).sum())}/{len(det)}'
    # A1-A6
    v['A1_detected'] = f'{len(det)}/{len(Pq)}'; v['A1_list'] = list(det)
    v['A1_neg3'] = list(Pq.index[Pq.z_s <= -3]); v['A1_positive'] = f'{int((Pq.thrust_m0 > 0).sum())}/{len(Pq)}'
    v['A1_big_detected'] = f'{int((Pq.loc[big].z_s >= 3).sum())}/{len(big)}'
    v['A2_all'] = wmean(Pq.thrust_m0, Pq.sd_s)
    v['A2_detected'] = wmean(Pq.loc[det].thrust_m0, Pq.loc[det].sd_s) if len(det) > 1 else None
    v['A2_median_detected'] = float(Pq.loc[det].thrust_m0.median()) if len(det) else None
    v['A2_median_all'] = float(Pq.thrust_m0.median())
    ep = E.set_index('episode').join(Pq[['thrust_m0', 'z_s', 'sd_s', 'thrust_lag', 'thrust_vr', 'thrust_c6_z', 'thrust_c6_m0',
                                         'thrust_wide_lag', 'thrust_wide_m0', 'thrust_wide_z', 'n_nb2', 'n_nb10', 'nb_kern',
                                         'n_sta', 'guard_dropped', 'max_share', 'max_share_ch', 'loso_m0_min',
                                         'loso_m0_max', 'thrust_zk_rms']], rsuffix='_u')
    ep['z_odd'] = zo; ep['z_even'] = ze
    ep['m0_total'] = ep.thrust_m0 * ep.n
    sz = pd.read_csv(f'{RES}/s11_depth/round6/wa_nw/tables/r6_sizes.csv').set_index('episode')
    ep = ep.join(sz[[c for c in ('G_WA', 'T', 'Du', 'L') if c in sz.columns]])
    dd = pd.read_csv(f'{RES}/s11_depth/round7/wa_nw/tables/r7_hint_dD.csv').set_index('episode')
    ep = ep.join(dd[['dD']].rename(columns={'dD': 'dvv_dD'}))
    for z in ('zone1', 'zone2', 'zone3'):
        if z in zt:
            ep[f'm0_{z}'] = zt[z]; ep[f'z_{z}'] = zz[z]
    ep.to_csv(f'{C.TB}/v2_wa_episodes.csv')
    b = ep[ep.thrust_m0.notna()]
    for nm, col in (('tremor', 'T'), ('duration', 'Du'), ('lfecount', 'L'), ('gnss_wa', 'G_WA')):
        if col in b:
            v[f'A3_total_vs_{nm}'] = perm_spearman(b.m0_total.values, b[col].values)
            v[f'A3_perlfe_vs_{nm}'] = perm_spearman(b.thrust_m0.values, b[col].values)
    Fm = U[U.unit.str.startswith('fam') & U.thrust_m0.notna() & (U.sd_s > 0)].copy()
    Fm['w'] = 1 / Fm.sd_s ** 2
    fam = Fm.groupby('unit').apply(lambda g: pd.Series(dict(
        n_ep=len(g), n_lfe=g.n.sum(), m0_wmean=(g.thrust_m0 * g.w).sum() / g.w.sum(), m0_sd=1 / np.sqrt(g.w.sum()),
        chi2_dof=((g.thrust_m0 - (g.thrust_m0 * g.w).sum() / g.w.sum()) ** 2 * g.w).sum() / max(len(g) - 1, 1),
        lat=g.src_lat.mean(), lon=g.src_lon.mean(), dep=g.src_dep.mean())))
    fam['z_wmean'] = fam.m0_wmean / fam.m0_sd
    fam.to_csv(f'{C.TB}/v2_wa_families.csv')
    v['A4_families'] = int(len(fam)); v['A4_n_z3'] = int((fam.z_wmean >= 3).sum()); v['A4_n_zm3'] = int((fam.z_wmean <= -3).sum())
    fw = 1 / fam.m0_sd ** 2; fmu = (fam.m0_wmean * fw).sum() / fw.sum()
    v['A4_between_family_chi2_dof'] = float(((fam.m0_wmean - fmu) ** 2 * fw).sum() / max(len(fam) - 1, 1))
    v['A4_within_family_chi2_dof_median'] = float(fam.chi2_dof.median())
    v['A6_total_vs_dvv'] = perm_spearman(b.m0_total.values, b.dvv_dD.values)
    v['A6_perlfe_vs_dvv'] = perm_spearman(b.thrust_m0.values, b.dvv_dD.values)
    # A5 between
    fb = f'{C.TB}/wa2_between_units.csv'
    if os.path.isfile(fb):
        B = pd.read_csv(fb)
        B['sd_s'] = B.thrust_ctrl_sd * scale; B['z_s'] = (B.thrust_m0 - B.thrust_ctrl_mean) / B.sd_s
        for k in KS:
            B[f'zk{k:+d}'] = (B[f'thrust_m0_k{k:+d}'] - B.thrust_ctrl_mean) / B.thrust_ctrl_sd
        B.to_csv(f'{C.TB}/v2_wa_between.csv', index=False)
        v['A5'] = {r.unit: dict(n=int(r.n), m0=float(r.thrust_m0), z=float(r.thrust_z), z_scaled=float(r.z_s),
                                lag=int(r.thrust_lag), vr=float(r.thrust_vr), m0_Pw=float(r.thrust_m0_Pw),
                                zk=[round(float(B.loc[i, f'zk{k:+d}']), 2) for k in KS],
                                ratio_to_ets=float(r.thrust_m0 / v['A2_all']['wmean']),
                                ratio_sd=float(np.hypot(r.sd_s / v['A2_all']['wmean'],
                                                        r.thrust_m0 * v['A2_all']['sd'] / v['A2_all']['wmean'] ** 2)))
                    for i, r in zip(B.index, B.itertuples())}
    # zone-1 grand stack wide scan (W1 failure branch, also a column): weighted mean of zone-1 wide-lag m0
    Z1 = Z[Z.unit == 'zone1']
    v['zone1_wide_lag_counts'] = {int(k): int(x) for k, x in Z1.thrust_wide_lag.value_counts().sort_index().items()}
    wa_fig(ep, Z)
    return v


def wa_fig(ep, Z):
    b = ep[ep.thrust_m0.notna()].copy()
    b['mid'] = pd.to_datetime(b.start) + (pd.to_datetime(b.end) - pd.to_datetime(b.start)) / 2
    fig, ax = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    ax[0].plot(b.mid, b.thrust_m0, color=GREY, lw=0.8)
    sc = ax[0].scatter(b.mid, b.thrust_m0, c=b.z_s.clip(-2, 8), cmap='viridis', vmin=-2, vmax=8,
                       s=12 + 60 * np.sqrt(b.n / b.n.max()), zorder=3)
    ax[0].axhline(0, color=GREY, lw=0.6); ax[0].set_ylabel('VLF moment per LFE (N m)')
    ax[0].set_title('S13 WA-NW corrected (origin = time + 0.5 s, deduped ±2 s, Slab2 depths): pooled, per major episode '
                    '(size ∝ √n)', fontsize=9)
    for i, z in enumerate(('zone1', 'zone2', 'zone3')):
        if f'm0_{z}' in b:
            ax[1].scatter(b.mid + pd.Timedelta(days=20 * (i - 1)), b[f'm0_{z}'], c=b[f'z_{z}'].clip(-2, 8), cmap='viridis',
                          vmin=-2, vmax=8, s=18, marker='osD'[i], label=z, zorder=3)
    ax[1].axhline(0, color=GREY, lw=0.6); ax[1].set_ylabel('per LFE, by zone (N m)')
    ax[1].legend(fontsize=7, title='o zone1 SE, s zone2, D zone3 NW', title_fontsize=7)
    ax[2].scatter(b.mid, b.thrust_lag, c=b.z_s.clip(-2, 8), cmap='viridis', vmin=-2, vmax=8, s=20, label='±8 s grid')
    ax[2].scatter(b.mid, b.thrust_wide_lag, facecolors='none', edgecolors='k', s=20, label='±40 s grid (column)')
    ax[2].axhspan(-4, 4, color='#dddddd', zorder=0); ax[2].set_ylabel('best lag (s)'); ax[2].legend(fontsize=7)
    fig.colorbar(sc, ax=ax, label='z vs scrambled control (σ-scaled if W3 failed)')
    fig.savefig(f'{C.FIG}/s13_v2_wa_episodes.png', dpi=130); plt.close(fig)


# ------------------------------------------------------------------------------------------------------------- SVI
def svi_pick(ep):
    for pre in ('scalef', 'refit', 'scale'):
        f = f'{C.TB}/{pre}_svi_{ep}_units.csv'
        if os.path.isfile(f):
            return pre
    return None


def refit_unit(st, S, drop_sta=(), dates=None):
    """Refit a stored unit with guard2 weights from its stacomp rows S; drop_sta weights set to 0.
    Returns m0, lag, analytic sigma, per-station Z single-channel (m0, sigma_white) at the unit lag."""
    sig = {s: g.set_index('comp').reindex(list('NEZ')).sigma_ctrl_nm.values * 1e-9 for s, g in S.groupby('sta', sort=False)}
    nus = {s: g.set_index('comp').reindex(list('NEZ')).n_used.values for s, g in S.groupby('sta', sort=False)}
    sts = [s for s in st['stk'] if s in sig]
    sig = {s: sig[s] for s in sts}; nus = {s: nus[s] for s in sts}
    wts, drop, _ = C.guard2(sig, nus)
    for s in drop_sta:
        if s in wts:
            wts[s] = np.zeros(3)
    if dates is not None:      # instrument vetoes (Amendment 5 addendum): zero a vetoed component if the unit overlaps it
        d0, d1 = (int(pd.Timestamp(x).value // 86400e9) for x in dates)
        for s in sts:
            for a0, a1, comps in vetoes(s):
                if a0 <= d1 and d0 <= a1:
                    wts[s] = wts[s].copy(); wts[s][comps] = 0.0
                    VETO_LOG.append(dict(sta=s, comps=str(comps), d0=dates[0], d1=dates[1]))
    if sum(w.sum() for w in wts.values()) == 0:
        return np.nan, np.nan, np.nan, {}
    stk = {s: st['stk'][s] for s in sts}; T = {s: st['tpl'][s] for s in sts}
    m0, L, _ = fit(stk, T, wts)
    num = den = 0.0
    for s in sts:
        tt = np.roll(T[s], L, axis=1)[:, FITM]
        num += (wts[s][:, None] ** 2 * tt ** 2 * sig[s][:, None] ** 2).sum()
        den += (wts[s][:, None] * tt ** 2).sum()
    single = {}
    Lu = st['lag']
    for s in sts:
        if drop[s][2]:
            continue
        tz = np.roll(T[s], Lu, axis=1)[2, FITM]; dz = st['stk'][s][2, FITM]
        single[s] = ((tz * dz).sum() / (tz * tz).sum(), sig[s][2] / np.sqrt((tz * tz).sum()))
    return m0, L, np.sqrt(num) / den, single


def svi():
    v = {}
    E = pd.read_csv(f'{C.TB}/svi_episodes.csv'); E = E[E.built]
    rows, srows, nbrows = [], [], []
    for ep in E.episode:
        pre = svi_pick(ep)
        if pre is None:
            continue
        Ut = pd.read_csv(f'{C.TB}/{pre}_svi_{ep}_units.csv').set_index('unit')
        S = pd.read_csv(f'{C.TB}/{pre}_svi_{ep}_stacomp.csv')
        st = pickle.load(open(f'{C.DATA}/stacks/{pre}_svi_{ep}.pkl', 'rb'))['units']
        if 'pooled' not in st:
            continue
        Sp = S[S.unit == 'pooled']
        dts = (r0.date0, r0.date1) if 'date0' in (r0 := Ut.loc['pooled']) else None
        m_all, L_all, sa_all, single = refit_unit(st['pooled'], Sp, dates=dts)
        m_out, L_out, sa_out, _ = refit_unit(st['pooled'], Sp, drop_sta=('CN.PGC',), dates=dts)
        v2tab = 'guard_rule' in Ut.columns and (Ut.guard_rule == C.GUARD_RULE).all()
        r = Ut.loc['pooled']
        ratio = float(r.thrust_ctrl_sd / sa_all) if v2tab else np.nan
        rows.append(dict(episode=ep, table=pre, guard_v2_table=v2tab, n=int(r.n), m0_stored=r.thrust_m0,
                         z_stored=r.thrust_z, ctrl_sd=r.thrust_ctrl_sd, m0_refit_all=m_all, lag_all=L_all,
                         sig_white_all=sa_all, m0_pgc_out=m_out, lag_pgc_out=L_out, sig_white_pgc_out=sa_out,
                         ctrl_over_white=ratio, pgc_share=float(r.get('max_share', np.nan)) if str(r.get('max_share_ch', '')).startswith('CN.PGC') else np.nan,
                         max_share_ch=r.get('max_share_ch', '')))
        if ep < '2014':
            for s, (m, sw) in single.items():
                srows.append(dict(episode=ep, sta=s, m0_Z=m, sig_white=sw))
        # fixed neighbour counts from the row file (pooled rows only, self excluded)
        R = pd.read_parquet(f'{C.DATA}/rows/svi_{ep}.parquet')
        nbrows.append(dict(episode=ep, **nb_counts(R[R.unit == 'pooled'].to.values, R[R.unit == 'pooled'].to.values),
                           n_nb10_old=float(r.n_nb10) if 'n_nb10' in r else np.nan))
    D = pd.DataFrame(rows)
    rmed = float(D.ctrl_over_white.median())
    D['ratio_used'] = D.ctrl_over_white.fillna(rmed)
    D['sd_all'] = D.sig_white_all * D.ratio_used; D['sd_pgc_out'] = D.sig_white_pgc_out * D.ratio_used
    D['z_pgc_out'] = D.m0_pgc_out / D.sd_pgc_out
    D.to_csv(f'{C.TB}/v2_svi_pgc_out.csv', index=False)
    v['ctrl_over_white_median'] = rmed
    v['refit_reproduces_stored_v2'] = float((D[D.guard_v2_table].m0_refit_all / D[D.guard_v2_table].m0_stored - 1).abs().max()) \
        if D.guard_v2_table.any() else None
    for lab, sel in (('2008_2026', D.episode >= '2008'), ('all', D.episode > '')):
        x = D[sel]
        v[f'pgc_all_{lab}'] = wmean(x.m0_refit_all, x.sd_all)
        v[f'pgc_out_{lab}'] = wmean(x.m0_pgc_out, x.sd_pgc_out)
        v[f'pgc_out_change_{lab}'] = v[f'pgc_out_{lab}']['wmean'] / v[f'pgc_all_{lab}']['wmean'] - 1
    v['veto_refit_drops'] = VETO_LOG
    v['pgc_out_positive'] = f'{int((D.m0_pgc_out > 0).sum())}/{len(D)}'
    v['pgc_out_z3'] = int((D.z_pgc_out >= 3).sum())
    Sz = pd.DataFrame(srows)
    if len(Sz):
        Sz = Sz.merge(D[['episode', 'ratio_used']], on='episode')
        Sz['sd_Z'] = Sz.sig_white * Sz.ratio_used
        Sz.to_csv(f'{C.TB}/v2_svi_single_Z_2003_2013.csv', index=False)
    NB = pd.DataFrame(nbrows); NB.to_csv(f'{C.TB}/v2_svi_nb_fixed.csv', index=False)
    v['nb_fixed_median'] = {c: float(NB[c].median()) for c in ('n_nb2', 'n_nb10', 'n_nb20', 'n_nb60', 'nb_kern')}
    v['nb_fixed_by_ep'] = NB.set_index('episode')[['n_nb10', 'nb_kern']].round(3).to_dict('index')
    # step 4 P-bins
    pb = []
    for ep in ('2004-07', '2005-09', '2012-08', '2013-09'):
        f = f'{C.TB}/pbin_svi_{ep}_units.csv'
        if os.path.isfile(f):
            pb.append(pd.read_csv(f).assign(episode=ep))
    if pb:
        PB = pd.concat(pb, ignore_index=True)
        PB.to_csv(f'{C.TB}/v2_svi_pbin.csv', index=False)
        t = PB.set_index(['episode', 'unit'])
        v['pbin'] = {f'{e}/{u}': dict(n=int(r.n), m0=float(r.thrust_m0), sd=float(r.thrust_ctrl_sd), z=float(r.thrust_z),
                                      lag=int(r.thrust_lag), nb10=float(r.n_nb10)) for (e, u), r in t.iterrows()}
        g = lambda e, u: t.loc[(e, u), 'thrust_m0'] if (e, u) in t.index else np.nan
        gs = lambda e, u: t.loc[(e, u), 'thrust_ctrl_sd'] if (e, u) in t.index else np.nan
        e4 = '2004-07'
        b12 = [g('2012-08', u) for u in ('pbin_lo', 'pbin_mid', 'pbin_hi')]
        spread12 = max(b12) / min(b12) if min(b12) > 0 else np.inf
        a = (g(e4, 'pbin_hi') >= 3.7e12 and g(e4, 'top3000') >= 3.7e12 and g(e4, 'pbin_lo') < 2e12 and spread12 < 1.5)
        bb = all(1e12 <= g(e4, u) <= 3e12 for u in ('pbin_lo', 'pbin_mid', 'pbin_hi', 'top3000', 'rand3000')) and \
            abs(g(e4, 'top3000') - g(e4, 'rand3000')) < 2 * np.hypot(gs(e4, 'top3000'), gs(e4, 'rand3000'))
        hl = lambda e: g(e, 'pbin_hi') / g(e, 'pbin_lo') if g(e, 'pbin_lo') > 0 else np.inf
        c = (max(hl('2004-07'), hl('2005-09')) >= 1.5) and (max(hl('2012-08'), hl('2013-09')) >= 1.5)
        v['step4_spread_2012_08'] = float(spread12)
        v['step4_hi_over_lo'] = {e: float(hl(e)) for e in ('2004-07', '2005-09', '2012-08', '2013-09')}
        v['step4_outcome'] = 'a_dilution' if a else ('b_era_difference' if bb else ('c_scales_with_strength' if c else 'mixed'))
        v['step4_flags'] = dict(a=bool(a), b=bool(bb), c=bool(c))
        pbin_fig(PB)
    # step 5 between variants
    fb = f'{C.TB}/svibe2_units.csv'
    if os.path.isfile(fb):
        B = pd.read_csv(fb)
        B.to_csv(f'{C.TB}/v2_svi_between.csv', index=False)
        ets = v['pgc_all_2008_2026']
        v['step5'] = {r.unit: dict(n=int(r.n), m0=float(r.thrust_m0), z=float(r.thrust_z), lag=int(r.thrust_lag),
                                   m0_Pw=float(r.thrust_m0_Pw), sd=float(r.thrust_ctrl_sd),
                                   ratio_to_ets=float(r.thrust_m0 / ets['wmean']),
                                   ratio_sd=float(np.hypot(r.thrust_ctrl_sd / ets['wmean'],
                                                           r.thrust_m0 * ets['sd'] / ets['wmean'] ** 2)))
                      for r in B.itertuples()}
        s5 = v['step5']
        v['step5_pass'] = bool(all(x['m0'] > 0 for x in s5.values()) and s5['BE_minus_top5days']['z'] >= 2
                               and s5['BE_2014_2019']['m0'] > 0 and s5['BE_2020_2026']['m0'] > 0)
    return v


def pbin_fig(PB):
    fig, ax = plt.subplots(figsize=(7, 4.2))
    xs = {'pbin_lo': 0, 'pbin_mid': 1, 'pbin_hi': 2, 'top3000': 3.2, 'rand3000': 3.8}
    for i, (e, g) in enumerate(PB.groupby('episode')):
        g = g[g.unit.isin(xs)]
        x = g.unit.map(xs) + 0.06 * (i - 1.5)
        bins = g[g.unit.str.startswith('pbin')].sort_values('unit', key=lambda s: s.map(xs))
        ax.plot(bins.unit.map(xs) + 0.06 * (i - 1.5), bins.thrust_m0, color=GREY, lw=0.8)
        sc = ax.scatter(x, g.thrust_m0, c=g.thrust_z.clip(-2, 8), cmap='viridis', vmin=-2, vmax=8, s=30,
                        marker='os^D'[i], label=e, zorder=3)
    ax.axhline(0, color=GREY, lw=0.6)
    ax.set_xticks(list(xs.values())); ax.set_xticklabels(['P .75–.90', 'P .90–.97', 'P .97–1', 'top 3k', 'random 3k'])
    ax.set_ylabel('VLF moment per LFE (N m)'); ax.legend(fontsize=7)
    fig.colorbar(sc, ax=ax, label='z vs scrambled control')
    ax.set_title('S13 SVI step 4: moment per LFE by detection probability', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{C.FIG}/s13_v2_svi_pbin.png', dpi=130); plt.close(fig)


def main():
    which = sys.argv[1:] or ['wa', 'svi']
    f = f'{C.TB}/v2_verdict.json'
    V = json.load(open(f)) if os.path.isfile(f) else {}
    if 'svi' in which:
        V['svi'] = svi()
    if 'wa' in which:
        V['wa'] = wa()
    json.dump(V, open(f, 'w'), indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print(json.dumps(V, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))


if __name__ == '__main__':
    main()
