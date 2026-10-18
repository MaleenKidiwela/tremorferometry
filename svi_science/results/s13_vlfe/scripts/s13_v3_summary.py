#!/usr/bin/env python
"""S13 Amendment 6 summary (PLAN.md Amendment 6 + addendum 6a, 2026-10-14 PT). Reads only. WA v3 (wa3_*: per-family
origin offsets, re-deduped, cluster GF locations), the per-family lag-prediction check (catalogue-only offsets vs the v2
stacks, which never used them), and SVI step 5 (nmt_svi_*: net_mad terciles within P >= 0.9).
Writes tables/v3_verdict.json, tables/v3_wa_units_all.csv, v3_wa_episodes.csv, v3_wa_families.csv, v3_wa_between.csv,
v3_family_lag_check.csv, v3_wa_{pooled,zone1}_grand_lagscan.csv, v3_svi_nmt.csv; figures s13_v3_wa_episodes.png,
s13_v3_family_lags.png, s13_v3_svi_nmt.png.

usage: s13_v3_summary.py"""
import os, sys, json, pickle, glob
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
from s13_stack import FITM
from s13_v2_summary import perm_spearman, wmean, KS, GREY

RES = f'{C.ROOT}/svi_science/results'
WIDE = np.arange(-40, 41)
VIR = dict(cmap='viridis', vmin=-2, vmax=8)


def grand_scan(units, lags=WIDE):
    """Joint fit over episodes: units = list of (n, store dict). m(L) = sum_e n_e sum_s w T_L stk / sum_e n_e sum_s w T_L^2."""
    sd = np.zeros(len(lags)); ss = np.zeros(len(lags)); num = 0.0
    for n, st in units:
        for s, stk in st['stk'].items():
            w = st['wts'][s]
            if not np.any(w > 0):
                continue
            x = np.nan_to_num(stk[:, FITM]) * (w[:, None] > 0)
            num += n * (w[:, None] * x ** 2).sum()
            for i, L in enumerate(lags):
                tt = np.roll(st['tpl'][s], L, axis=1)[:, FITM]
                sd[i] += n * (w[:, None] * tt * x).sum()
                ss[i] += n * (w[:, None] * tt * tt).sum()
    m = np.where(ss > 0, sd / np.where(ss > 0, ss, 1), np.nan)
    vr = np.where(ss > 0, sd ** 2 / np.where(ss > 0, ss, 1), np.nan) / num if num > 0 else m * np.nan
    return m, vr


def load_stores(prefix, eps):
    out = {}
    for e in eps:
        f = f'{C.DATA}/stacks/{prefix}_{e}.pkl'
        if os.path.isfile(f):
            out[e] = pickle.load(open(f, 'rb'))['units']
    return out


# ------------------------------------------------------------------------------------------------------------- WA
def wa(v):
    E = pd.read_csv(f'{C.TB}/wa3_episodes.csv')
    U = []
    for e in E[E.built].episode:
        f = f'{C.TB}/wa3_{e}_units.csv'
        if os.path.isfile(f):
            U.append(pd.read_csv(f).assign(episode=e))
    U = pd.concat(U, ignore_index=True)
    U.to_csv(f'{C.TB}/v3_wa_units_all.csv', index=False)
    for k in KS:
        U[f'zk{k:+d}'] = (U[f'thrust_m0_k{k:+d}'] - U.thrust_ctrl_mean) / U.thrust_ctrl_sd
    Pq = U[U.unit == 'pooled'].set_index('episode')
    v['n_built'] = int(E.built.sum()); v['n_stacked'] = int(len(Pq))
    # W1 (Amendment 6: mode in +1..+4 s, <= 10 % at +-8 s) -- the headline gate
    lags = Pq.thrust_lag
    vc = lags.value_counts()
    v['W1_lag_counts'] = {int(k): int(x) for k, x in vc.sort_index().items()}
    v['W1_lag_mode'] = int(vc.idxmax()); v['W1_frac_at_pm8'] = float((lags.abs() == 8).mean())
    v['W1_lag_median'] = float(lags.median())
    v['W1_pass'] = bool(1 <= v['W1_lag_mode'] <= 4 and v['W1_frac_at_pm8'] <= 0.10)
    ofs = json.load(open(f'{C.TB}/wa3_offsets_verdict.json'))
    sh = ofs['cvm_constant_shift_s']
    v['W1_cvm_column'] = dict(shift_s=sh, lag_mode_cvm=float(v['W1_lag_mode'] - sh),
                              note='under the CVM v1.7 anchor every origin moves by shift_s, so every lag by -shift_s')
    v['col_wide_lag_mode'] = int(Pq.thrust_wide_lag.value_counts().idxmax())
    v['col_wide_lag_median'] = float(Pq.thrust_wide_lag.median())
    # W2
    v['W2_nb2_mean'] = float(Pq.n_nb2.mean()); v['W2_pass'] = bool(v['W2_nb2_mean'] < 0.05)
    v['nb_pooled_mean'] = {c: float(Pq[c].mean()) for c in ('n_nb2', 'n_nb10', 'n_nb20', 'n_nb60', 'nb_kern')}
    # W3
    zk = Pq[[f'zk{k:+d}' for k in KS]].values.ravel(); zk = zk[np.isfinite(zk)]
    sdk = float(np.std(zk)); v['W3_zk_sd'] = sdk; v['W3_zk_mean'] = float(np.mean(zk)); v['W3_n'] = int(len(zk))
    v['W3_zk_absge3'] = int((np.abs(zk) >= 3).sum())
    v['W3_pass'] = bool(0.85 <= sdk <= 1.15)
    scale = 1.0 if v['W3_pass'] else sdk
    v['sigma_scale_applied'] = scale
    U['sd_s'] = U.thrust_ctrl_sd * scale; U['z_s'] = (U.thrust_m0 - U.thrust_ctrl_mean) / U.sd_s
    U['sd_f2'] = U.thrust_fix2_ctrl_sd * scale; U['z_f2'] = (U.thrust_fix2_m0 - U.thrust_fix2_ctrl_mean) / U.sd_f2
    Pq = U[U.unit == 'pooled'].set_index('episode')
    big = Pq.index[Pq.n >= 300]
    v['n_big'] = int(len(big))
    Z = U[U.unit.str.startswith('zone')]
    out_zone = {}
    for lab, mcol, zcol, scol in (('grid', 'thrust_m0', 'z_s', 'sd_s'), ('fix2', 'thrust_fix2_m0', 'z_f2', 'sd_f2')):
        zt = Z.pivot_table(index='episode', columns='unit', values=mcol)
        ok = []
        for z in ('zone1', 'zone2', 'zone3'):
            x = zt.reindex(big)[z].dropna() if z in zt else pd.Series(dtype=float)
            Zu = Z[Z.unit == z]
            out_zone[f'{lab}_{z}'] = dict(pos=f'{int((x > 0).sum())}/{int(len(x))}', wmean=wmean(Zu[mcol], Zu[scol]),
                                          nz3=int((Zu[zcol] >= 3).sum()), nzm3=int((Zu[zcol] <= -3).sum()))
            if z != 'zone1':
                ok.append(len(x) > 0 and (x > 0).sum() > len(x) / 2)
        out_zone[f'{lab}_W4_pass'] = bool(all(ok))
    v['zones'] = out_zone
    v['W4_pass'] = out_zone['grid_W4_pass']
    # W5
    zo = U[U.unit == 'pooled_odd'].set_index('episode').z_s; ze = U[U.unit == 'pooled_even'].set_index('episode').z_s
    j = zo.index.intersection(ze.index)
    v['W5_rho_odd_even'] = float(spearmanr(zo[j], ze[j]).statistic); v['W5_n'] = int(len(j))
    v['W5_pass'] = bool(v['W5_rho_odd_even'] > 0.4)
    zo2 = U[U.unit == 'pooled_odd'].set_index('episode').z_f2; ze2 = U[U.unit == 'pooled_even'].set_index('episode').z_f2
    v['fix2_rho_odd_even'] = float(spearmanr(zo2[j], ze2[j]).statistic)
    # A1 / A2 on the headline grid and on the fixed +2 s column
    for lab, mcol, zcol, scol in (('grid', 'thrust_m0', 'z_s', 'sd_s'), ('fix2', 'thrust_fix2_m0', 'z_f2', 'sd_f2')):
        det = Pq.index[Pq[zcol] >= 3]
        v[f'{lab}_A1_detected'] = f'{len(det)}/{len(Pq)}'
        v[f'{lab}_A1_big_detected'] = f'{int((Pq.loc[big, zcol] >= 3).sum())}/{len(big)}'
        v[f'{lab}_A1_neg3'] = list(Pq.index[Pq[zcol] <= -3]); v[f'{lab}_A1_positive'] = f'{int((Pq[mcol] > 0).sum())}/{len(Pq)}'
        v[f'{lab}_A2_all'] = wmean(Pq[mcol], Pq[scol])
        v[f'{lab}_A2_big'] = wmean(Pq.loc[big, mcol], Pq.loc[big, scol])
        v[f'{lab}_A2_median_all'] = float(Pq[mcol].median())
    svi = json.load(open(f'{C.TB}/v2_verdict.json'))['svi']['pgc_all_2008_2026']
    v['svi_ref_2008_2026'] = svi
    for lab in ('grid', 'fix2'):
        a = v[f'{lab}_A2_all']
        v[f'{lab}_ratio_to_svi'] = float(a['wmean'] / svi['wmean'])
        v[f'{lab}_ratio_to_svi_sd'] = float(np.hypot(a['sd'] / svi['wmean'], a['wmean'] * svi['sd'] / svi['wmean'] ** 2))
    v['size_pass'] = bool(v['grid_A2_all']['wmean'] >= 3.5e12)
    # episodes table + A3 / A6
    ep = E.set_index('episode').join(Pq[['thrust_m0', 'z_s', 'sd_s', 'thrust_lag', 'thrust_vr', 'thrust_fix2_m0', 'z_f2',
                                         'sd_f2', 'thrust_wide_lag', 'thrust_wide_m0', 'thrust_wide_z', 'n_sta',
                                         'guard_dropped', 'max_share', 'max_share_ch', 'loso_m0_min', 'loso_m0_max',
                                         'thrust_zk_rms']], rsuffix='_u')
    ep['z_odd'] = zo; ep['z_even'] = ze
    ep['m0_total'] = ep.thrust_m0 * ep.n; ep['m0_total_f2'] = ep.thrust_fix2_m0 * ep.n
    for z in ('zone1', 'zone2', 'zone3'):
        Zu = Z[Z.unit == z].set_index('episode')
        ep[f'm0_{z}'] = Zu.thrust_m0; ep[f'z_{z}'] = Zu.z_s; ep[f'm0f2_{z}'] = Zu.thrust_fix2_m0; ep[f'zf2_{z}'] = Zu.z_f2
    sz = pd.read_csv(f'{RES}/s11_depth/round6/wa_nw/tables/r6_sizes.csv').set_index('episode')
    ep = ep.join(sz[[c for c in ('G_WA', 'T', 'Du', 'L') if c in sz.columns]])
    dd = pd.read_csv(f'{RES}/s11_depth/round7/wa_nw/tables/r7_hint_dD.csv').set_index('episode')
    ep = ep.join(dd[['dD']].rename(columns={'dD': 'dvv_dD'}))
    ep.to_csv(f'{C.TB}/v3_wa_episodes.csv')
    b = ep[ep.thrust_m0.notna()]
    for lab, tot, per in (('grid', 'm0_total', 'thrust_m0'), ('fix2', 'm0_total_f2', 'thrust_fix2_m0')):
        for nm, col in (('tremor', 'T'), ('duration', 'Du'), ('lfecount', 'L'), ('gnss_wa', 'G_WA'), ('dvv', 'dvv_dD')):
            if col in b:
                v[f'{lab}_A3_total_vs_{nm}'] = perm_spearman(b[tot].values, b[col].values)
                v[f'{lab}_A3_perlfe_vs_{nm}'] = perm_spearman(b[per].values, b[col].values)
    # A4 families (headline grid and fix2)
    Fm = U[U.unit.str.startswith('fam') & U.thrust_m0.notna() & (U.sd_s > 0)].copy()
    rows = []
    for u, g in Fm.groupby('unit'):
        w = 1 / g.sd_s ** 2; w2 = 1 / g.sd_f2 ** 2
        mu = (g.thrust_m0 * w).sum() / w.sum(); mu2 = (g.thrust_fix2_m0 * w2).sum() / w2.sum()
        rows.append(dict(unit=u, n_ep=len(g), n_lfe=int(g.n.sum()), m0_wmean=mu, m0_sd=1 / np.sqrt(w.sum()),
                         m0f2_wmean=mu2, m0f2_sd=1 / np.sqrt(w2.sum()), lag_mode=int(g.thrust_lag.value_counts().idxmax()),
                         lat=g.src_lat.mean(), lon=g.src_lon.mean(), dep=g.src_dep.mean()))
    fam = pd.DataFrame(rows).set_index('unit')
    fam['z_wmean'] = fam.m0_wmean / fam.m0_sd; fam['zf2_wmean'] = fam.m0f2_wmean / fam.m0f2_sd
    fam.to_csv(f'{C.TB}/v3_wa_families.csv')
    v['A4_families'] = int(len(fam)); v['A4_n_z3'] = int((fam.z_wmean >= 3).sum()); v['A4_n_zm3'] = int((fam.z_wmean <= -3).sum())
    v['A4_fix2_n_z3'] = int((fam.zf2_wmean >= 3).sum()); v['A4_fix2_n_zm3'] = int((fam.zf2_wmean <= -3).sum())
    # A5 between
    fb = f'{C.TB}/wa3_between_units.csv'
    if os.path.isfile(fb):
        B = pd.read_csv(fb)
        B['sd_s'] = B.thrust_ctrl_sd * scale; B['z_s'] = (B.thrust_m0 - B.thrust_ctrl_mean) / B.sd_s
        B['z_f2'] = (B.thrust_fix2_m0 - B.thrust_fix2_ctrl_mean) / (B.thrust_fix2_ctrl_sd * scale)
        B.to_csv(f'{C.TB}/v3_wa_between.csv', index=False)
        a2 = v['grid_A2_all']
        v['A5'] = {r.unit: dict(n=int(r.n), m0=float(r.thrust_m0), z=float(r.z_s), lag=int(r.thrust_lag),
                                m0_fix2=float(r.thrust_fix2_m0), z_fix2=float(r.z_f2), m0_Pw=float(r.thrust_m0_Pw),
                                ratio_to_wa_ets=float(r.thrust_m0 / a2['wmean']))
                   for r in B.itertuples()}
    # grand scans (v3, joint fit over episodes) for pooled and zone1
    S3 = load_stores('wa3', list(Pq.index))
    for un in ('pooled', 'zone1'):
        m, vr = grand_scan([(int(Pq.loc[e, 'n']) if un == 'pooled' else int(U[(U.episode == e) & (U.unit == un)].n.iloc[0]), st[un])
                            for e, st in S3.items() if un in st])
        pd.DataFrame(dict(lag_s=WIDE, m0_per_lfe=m, vr=vr)).to_csv(f'{C.TB}/v3_wa_{un}_grand_lagscan.csv', index=False)
        k = int(np.nanargmax(m)); v[f'grand_{un}_peak'] = dict(lag=int(WIDE[k]), m0=float(m[k]),
                                                              m0_at_2=float(m[WIDE == 2][0]),
                                                              neg_lobe_lag=int(WIDE[int(np.nanargmin(m))]),
                                                              neg_lobe_m0=float(np.nanmin(m)))
    wa_fig(ep, Z)
    fam_check(v)
    return v


def fam_check(v):
    """Per-family lag prediction: catalogue-only offset_f + 2 s vs the family's v2 grand-stack lag (origin time + 0.5;
    the v2 stacks never saw the offsets), and the v3 family grand-stack lag (expected near +2)."""
    O = pd.read_csv(f'{C.TB}/wa_family_offsets.csv', index_col=0)
    L = pd.read_csv(f'{C.TB}/wa3_location_flags.csv', index_col=0)
    rows = []
    for tag, tb in (('v2', 'v2_wa_units_all.csv'), ('v3', 'v3_wa_units_all.csv')):
        U = pd.read_csv(f'{C.TB}/{tb}')
        Fm = U[U.unit.str.startswith('fam')]
        stores = load_stores('wa2' if tag == 'v2' else 'wa3', sorted(Fm.episode.unique()))
        for u, g in Fm.groupby('unit'):
            items = [(int(r.n), stores[r.episode][u]) for r in g.itertuples() if r.episode in stores and u in stores[r.episode]]
            if not items:
                continue
            m, vr = grand_scan(items)
            sdw = np.sqrt(np.sum((g.n.values * g.thrust_wide_ctrl_sd.values) ** 2)) / g.n.sum()
            kp = int(np.nanargmax(m)); ke = int(np.nanargmax(vr))
            rows.append(dict(fam=int(u[3:]), stack=tag, n_ep=len(items), n_lfe=int(g.n.sum()), lag_pos_peak=int(WIDE[kp]),
                             m0_pos_peak=float(m[kp]), z_approx=float(m[kp] / sdw) if sdw > 0 else np.nan,
                             lag_best_energy=int(WIDE[ke]), m0_best_energy=float(m[ke]), m0_at_pred=np.nan))
    T = pd.DataFrame(rows)
    T = T.pivot_table(index='fam', columns='stack', values=['n_ep', 'n_lfe', 'lag_pos_peak', 'm0_pos_peak', 'z_approx',
                                                             'lag_best_energy'])
    T.columns = [f'{a}_{b}' for a, b in T.columns]
    T = T.join(O[['offset_f', 'flag', 'n_clear', 'delta_f', 'suspect']]).join(L[['gf_source', 'cluster', 'sus_codetect']])
    T['pred_v2_lag'] = T.offset_f + 2.0
    T['resid_v2'] = T.lag_pos_peak_v2 - T.pred_v2_lag
    T.to_csv(f'{C.TB}/v3_family_lag_check.csv')
    s = T[(T.n_ep_v2 >= 3) & (T.z_approx_v2 >= 3)]
    v['famcheck_n'] = int(len(s))
    v['famcheck_resid_median'] = float(s.resid_v2.median()) if len(s) else None
    v['famcheck_resid_mad'] = float((s.resid_v2 - s.resid_v2.median()).abs().median()) if len(s) else None
    v['famcheck_within_2s'] = f'{int((s.resid_v2.abs() <= 2).sum())}/{len(s)}'
    v['famcheck_spearman_pred_obs'] = float(spearmanr(s.pred_v2_lag, s.lag_pos_peak_v2).statistic) if len(s) > 4 else None
    v['famcheck_by_flag'] = {f: dict(n=int(len(g)), resid_med=float(g.resid_v2.median()), within2=int((g.resid_v2.abs() <= 2).sum()))
                             for f, g in s.groupby('flag')}
    s3 = T[(T.n_ep_v3 >= 3) & (T.z_approx_v3 >= 3)]
    v['famcheck_v3_n'] = int(len(s3)); v['famcheck_v3_lag_median'] = float(s3.lag_pos_peak_v3.median()) if len(s3) else None
    v['famcheck_v3_within_pm2_of_2'] = f'{int(((s3.lag_pos_peak_v3 - 2).abs() <= 2).sum())}/{len(s3)}'
    fam_fig(T)


def wa_fig(ep, Z):
    b = ep[ep.thrust_m0.notna()].copy()
    b['mid'] = pd.to_datetime(b.start) + (pd.to_datetime(b.end) - pd.to_datetime(b.start)) / 2
    fig, ax = plt.subplots(4, 1, figsize=(10, 12), sharex=False)
    sz = 12 + 60 * np.sqrt(b.n / b.n.max())
    ax[0].plot(b.mid, b.thrust_m0, color=GREY, lw=0.8)
    sc = ax[0].scatter(b.mid, b.thrust_m0, c=b.z_s.clip(-2, 8), s=sz, zorder=3, **VIR)
    ax[0].axhline(0, color=GREY, lw=0.6); ax[0].set_ylabel('per LFE, ±8 s grid (N m)')
    ax[0].set_title('S13 WA-NW v3 (per-family origin offsets, deduped ±2 s, cluster GF locations): pooled per major '
                    '(size ∝ √n)', fontsize=9)
    ax[1].plot(b.mid, b.thrust_fix2_m0, color=GREY, lw=0.8)
    ax[1].scatter(b.mid, b.thrust_fix2_m0, c=b.z_f2.clip(-2, 8), s=sz, zorder=3, **VIR)
    ax[1].axhline(0, color=GREY, lw=0.6); ax[1].set_ylabel('per LFE, fixed lag +2 s (N m)')
    for i, z in enumerate(('zone1', 'zone2', 'zone3')):
        ax[2].scatter(b.mid + pd.Timedelta(days=20 * (i - 1)), b[f'm0_{z}'], c=b[f'z_{z}'].clip(-2, 8), s=18,
                      marker='osD'[i], label=z, zorder=3, **VIR)
    ax[2].axhline(0, color=GREY, lw=0.6); ax[2].set_ylabel('per LFE by zone, ±8 s grid (N m)')
    ax[2].legend(fontsize=7, title='o zone1 SE, s zone2, D zone3 NW', title_fontsize=7)
    ax[3].scatter(b.mid, b.thrust_lag, c=b.z_s.clip(-2, 8), s=20, label='±8 s grid', **VIR)
    ax[3].scatter(b.mid, b.thrust_wide_lag, facecolors='none', edgecolors='k', s=20, label='±40 s grid (column)')
    ax[3].axhspan(1, 4, color='#dddddd', zorder=0, label='W1 window +1…+4 s'); ax[3].set_ylabel('best lag (s)')
    ax[3].set_ylim(-42, 42); ax[3].legend(fontsize=7)
    fig.colorbar(sc, ax=ax, label='z vs scrambled control (σ-scaled if W3 failed)')
    fig.savefig(f'{C.FIG}/s13_v3_wa_episodes.png', dpi=130); plt.close(fig)
    # grand scans v2 vs v3
    fig, ax = plt.subplots(figsize=(8, 4))
    for f, c, ls in (('v2_wa_pooled', 'C1', ':'), ('v2_wa_zone1', 'C0', ':'), ('v3_wa_pooled', 'C1', '-'), ('v3_wa_zone1', 'C0', '-')):
        p = f'{C.TB}/{f}_grand_lagscan.csv'
        if os.path.isfile(p):
            G = pd.read_csv(p); ax.plot(G.lag_s, G.m0_per_lfe, color=c, ls=ls, label=f.replace('_wa_', ' '))
    ax.axhline(0, color=GREY, lw=0.6); ax.axvspan(1, 4, color='#eeeeee', zorder=0)
    ax.set_xlabel('lag (s): v2 relative to time + 0.5; v3 relative to time + 0.5 + offset_f')
    ax.set_ylabel('M0 per LFE at that lag (N m)'); ax.legend(fontsize=7)
    ax.set_title('S13 WA grand stacks (all episodes, joint fit): v2 (dotted) vs v3 (solid)', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{C.FIG}/s13_v3_wa_grandscan.png', dpi=130); plt.close(fig)


def fam_fig(T):
    s = T[T.n_ep_v2 >= 3].copy()
    fig, ax = plt.subplots(figsize=(6, 5.5))
    sc = ax.scatter(s.pred_v2_lag, s.lag_pos_peak_v2, c=s.z_approx_v2.clip(-2, 8), s=12 + 3 * np.sqrt(s.n_lfe_v2),
                    edgecolors=np.where(s.flag == 'network', 'k', 'r'), **VIR)
    lim = [-25, 10]
    ax.plot(lim, lim, color=GREY, lw=0.8); ax.set_xlim(lim); ax.set_ylim(-42, 42)
    for r in s[s.z_approx_v2 >= 5].itertuples():
        ax.annotate(str(r.Index), (r.pred_v2_lag, r.lag_pos_peak_v2), fontsize=6)
    ax.set_xlabel('predicted lag in the v2 stacks = catalogue offset_f + 2 s'); ax.set_ylabel('observed v2 family grand-stack lag (s)')
    ax.set_title('S13 WA: catalogue-only offsets vs observed VLF lag per family\n(black edge = network offset, red = clear-only/default)',
                 fontsize=8)
    fig.colorbar(sc, ax=ax, label='approx. z of the family grand stack')
    fig.tight_layout(); fig.savefig(f'{C.FIG}/s13_v3_family_lags.png', dpi=130); plt.close(fig)


# ------------------------------------------------------------------------------------------------------------- SVI
def svi(v):
    rows = []
    for ep in ('2004-07', '2005-09', '2012-08', '2013-09'):
        f = f'{C.TB}/nmt_svi_{ep}_units.csv'
        if os.path.isfile(f):
            rows.append(pd.read_csv(f).assign(episode=ep))
    if not rows:
        return v
    N = pd.concat(rows, ignore_index=True)
    Cm = pd.read_csv(f'{C.TB}/nmt_svi_composition.csv')
    N = N.merge(Cm[['episode', 'unit', 'net_mad_med', 'P_med']], on=['episode', 'unit'], how='left')
    N.to_csv(f'{C.TB}/v3_svi_nmt.csv', index=False)
    t = N.set_index(['episode', 'unit'])
    res, n15, n_eq = {}, 0, 0
    for ep in N.episode.unique():
        lo, hi = t.loc[(ep, 'nm_lo'), 'thrust_m0'], t.loc[(ep, 'nm_hi'), 'thrust_m0']
        r = (np.inf if hi > 0 else np.nan) if lo <= 0 else hi / lo
        res[ep] = dict(lo=float(lo), mid=float(t.loc[(ep, 'nm_mid'), 'thrust_m0']), hi=float(hi), ratio=float(r),
                       z_lo=float(t.loc[(ep, 'nm_lo'), 'thrust_z']), z_mid=float(t.loc[(ep, 'nm_mid'), 'thrust_z']),
                       z_hi=float(t.loc[(ep, 'nm_hi'), 'thrust_z']),
                       sd_lo=float(t.loc[(ep, 'nm_lo'), 'thrust_ctrl_sd']), sd_hi=float(t.loc[(ep, 'nm_hi'), 'thrust_ctrl_sd']),
                       hi_minus_lo_sigma=float((hi - lo) / np.hypot(t.loc[(ep, 'nm_lo'), 'thrust_ctrl_sd'],
                                                                    t.loc[(ep, 'nm_hi'), 'thrust_ctrl_sd'])))
        n15 += bool(r >= 1.5); n_eq += bool(0.67 <= r < 1.5)
    v['nmt'] = res
    v['nmt_n_ratio_ge_1p5'] = n15; v['nmt_n_ratio_0p67_1p5'] = n_eq
    v['nmt_outcome'] = 'scales_with_detection_amplitude' if n15 >= 3 else ('does_not_scale_at_fixed_P' if n_eq >= 3 else 'mixed')
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    xs = {'nm_lo': 0, 'nm_mid': 1, 'nm_hi': 2}
    for i, (e, g) in enumerate(N.groupby('episode')):
        g = g.sort_values('unit', key=lambda s: s.map(xs))
        x = g.unit.map(xs) + 0.06 * (i - 1.5)
        ax.plot(x, g.thrust_m0, color=GREY, lw=0.8)
        sc = ax.scatter(x, g.thrust_m0, c=g.thrust_z.clip(-2, 8), s=30, marker='os^D'[i], label=e, zorder=3, **VIR)
    ax.axhline(0, color=GREY, lw=0.6)
    ax.set_xticks([0, 1, 2]); ax.set_xticklabels(['net_mad low third', 'middle', 'high third'])
    ax.set_ylabel('VLF moment per LFE (N m)'); ax.legend(fontsize=7)
    fig.colorbar(sc, ax=ax, label='z vs scrambled control')
    ax.set_title('S13 SVI step 5 (Amendment 6): moment per LFE by detection amplitude, P ≥ 0.9 only', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{C.FIG}/s13_v3_svi_nmt.png', dpi=130); plt.close(fig)
    return v


def main():
    V = {'offsets': json.load(open(f'{C.TB}/wa3_offsets_verdict.json')),
         'rows': json.load(open(f'{C.TB}/wa3_rows_verdict.json'))}
    V['offsets'].pop('clusters', None)
    caps = []
    if not V['offsets']['step1_accept']:
        caps.append('step1: leave-one-out network-arrival SD %.2f s (rule < 3 s)' % V['offsets']['net_primary']['loo_dev_sd'])
    if not V['rows']['step3_accept']:
        caps.append('step3: removed %.1f %% (rule 25-35 %%), n_nb10 %.3f (rule < 0.15)' % (
            100 * V['rows']['removed_frac'], V['rows']['nb10_mean']))
    V['wa'] = wa({})
    V['svi'] = svi({})
    w = V['wa']
    V['headline_caps'] = caps
    if V['offsets'].get('step1_fail_stop'):
        V['wa_headline'] = 'stopped at step 1'
    elif w['W1_pass'] and not caps:
        V['wa_headline'] = 'headline'
    else:
        V['wa_headline'] = 'present, size provisional'
    json.dump(V, open(f'{C.TB}/v3_verdict.json', 'w'), indent=1,
              default=lambda o: o.item() if hasattr(o, 'item') else str(o))
    print(json.dumps(V, indent=1, default=lambda o: o.item() if hasattr(o, 'item') else str(o)))


if __name__ == '__main__':
    main()
