#!/usr/bin/env python
"""Step 2-3 analysis: H&G thresholds, reproduced catalogue, known-answer check K1-K4 (PLAN.md §1), H&G-style noise rate.
Inputs data/mf/real_T*.npz, noise_*.npz. Outputs tables/thresholds.csv, tables/catalogue_{prim,snr}.csv,
tables/period_counts.csv, tables/known_answer.csv, figures/fig6_repro.png, figures/fig7_repro.png,
figures/separation_real_vs_noise.png.  Counting uses UTC dates; detections at a template's own origin (+-90 s) are not
'new' and are excluded."""
import os, sys, glob, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hg_common as H

VARS = ['prim', 'snr']


def load():
    R = {K: dict(np.load(f'{H.MF}/real_T{K}.npz', allow_pickle=True)) for K in range(1, 6)}
    NZ = [dict(np.load(f, allow_pickle=True)) for f in sorted(glob.glob(f'{H.MF}/noise_[0-9][0-9][0-9].npz'))]
    return R, NZ


def noise_peaks(z, name):
    """positive peaks (idx, val) of noise run z under config `name` (e.g. 'T1_prim')."""
    for q in range(int(z['ncfg'])):
        if name in list(z[f'c{q}_names']):
            return z[f'c{q}_pi'], z[f'c{q}_pv'], int(z[f'c{q}_n'])
    raise KeyError(name)


def thresholds(R, NZ):
    rows = []
    for K in range(1, 6):
        for v in VARS:
            third_a, third_b = [], []
            n = len(R[K]['prim']) if v == 'prim' else len(R[K]['snr_chans'])
            key = 'real' if v == 'prim' else 'snr'
            med, mad = float(R[K][f'{key}_med']), float(R[K][f'{key}_mad'])
            for z in NZ:
                pi, pv, n = noise_peaks(z, f'T{K}_{v}')
                o = np.argsort(pv)[::-1]
                third_a.append(pv[o[2]] if o.size > 2 else np.nan)
                ns = np.abs(pi - int(z['i0'])) > H.L
                pvn = np.sort(pv[ns])[::-1]
                third_b.append(pvn[2] if pvn.size > 2 else np.nan)
            pool = []
            for z in NZ:
                pi, pv, _ = noise_peaks(z, f'T{K}_{v}')
                pool.append(pv[np.abs(pi - int(z['i0'])) > H.L])
            pool = np.sort(np.concatenate(pool))[::-1]
            kth = int(round(0.79 * H.NDAY * len(NZ)))                 # 0.79 per template per day over the full run
            thr_rate = float(pool[kth - 1]) if 0 < kth <= pool.size else np.nan
            rows.append({'template': K, 'variant': v, 'n_chan': n, 'n_noise': len(NZ),
                         'thr_6MAD': med + 6 * mad, 'thr_rate079': thr_rate,
                         'thrA_3rd_incl_self': np.nanmean(third_a), 'thrA_sd': np.nanstd(third_a),
                         'thrB_3rd_excl_self': np.nanmean(third_b), 'median_S': med, 'MAD_S': mad,
                         'six_MAD': 6 * mad, 'med_plus_6MAD': med + 6 * mad,
                         'thrA_in_MAD': (np.nanmean(third_a) - med) / mad if mad > 0 else np.nan,
                         'thrA_mean_cc': np.nanmean(third_a) / n})
    return pd.DataFrame(rows)


def catalogue(R, TH, v='prim', name=None, thr_col='thrA_3rd_incl_self'):
    """Detections above threshold for sum `name` (default real/snr) of each template."""
    out = []
    TT = H.templates()
    for K in range(1, 6):
        key = name or ('real' if v == 'prim' else 'snr')
        thr = float(TH[(TH.template == K) & (TH.variant == v)][thr_col].iloc[0])
        pi, pv, pnl = R[K][f'{key}_pi'], R[K][f'{key}_pv'], R[K][f'{key}_pnl']
        sel = pv > thr
        if key in ('real', 'snr'):
            sel &= np.abs(pi - int(R[K]['i0'])) > H.L
        out.append(pd.DataFrame({'template': K, 'i': pi[sel], 'S': pv[sel], 'n_live': pnl[sel], 'thr': thr}))
    C = pd.concat(out, ignore_index=True)
    C['time_utc'] = H.itime(C.i.values)
    C['date'] = C.time_utc.dt.floor('D')
    return C


def period_counts(C, label):
    rows = []
    for K in range(1, 6):
        c = C[C.template == K]
        for p in H.PER:
            a, b = H.per_idx(p)
            n = int(((c.i >= a) & (c.i < b)).sum())
            rows.append({'set': label, 'template': K, 'period': p, 'n': n, 'days': (b - a) / 86400,
                         'rate_per_day': n / ((b - a) / 86400)})
    return pd.DataFrame(rows)


def noise_rate(NZ, TH, v='prim', thk='thrA_3rd_incl_self'):
    """H&G-style false rate: noise-template detections above the template threshold (self excluded), per day,
    averaged over the noise templates; per template configuration and period."""
    rows = []
    for K in range(1, 6):
        thr = float(TH[(TH.template == K) & (TH.variant == v)][thk].iloc[0])
        for p in ['all18', 'ETS2011', 'interETS', 'ETS2012', 'win2011']:
            a, b = H.per_idx(p)
            cnt = []
            for z in NZ:
                pi, pv, n = noise_peaks(z, f'T{K}_{v}')
                s = (pv > thr) & (np.abs(pi - int(z['i0'])) > H.L) & (pi >= a) & (pi < b)
                cnt.append(s.sum())
            rows.append({'template': K, 'variant': v, 'thr_kind': thk, 'period': p,
                         'noise_rate_per_day': np.mean(cnt) / ((b - a) / 86400),
                         'noise_rate_sd_over_templates': np.std(cnt) / ((b - a) / 86400)})
    return pd.DataFrame(rows)


def tremor_daily(lat, lon, rkm):
    T = pd.read_csv(f'{H.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv', usecols=['time', 'lat', 'lon'])
    T['time'] = pd.to_datetime(T.time)
    T = T[(T.time >= H.T0) & (T.time < H.T1)]
    d = 6371 * np.hypot(np.radians(T.lat - lat), np.radians(T.lon - lon) * np.cos(np.radians(lat)))
    s = T[d <= rkm].time.dt.floor('D').value_counts().sort_index()
    return s.reindex(pd.date_range(H.T0, H.T1 - pd.Timedelta(days=1)), fill_value=0)


THK = ['thrA_3rd_incl_self', 'thr_6MAD', 'thr_rate079']
RUNS = [('prim', t) for t in THK] + [('snr', 'thrA_3rd_incl_self')]


def known_answer(P, NR, v, thk):
    F7 = pd.read_csv(f'{H.TB}/paper_fig7_period_fractions.csv')
    tot = P.groupby('period').agg(n=('n', 'sum'), days=('days', 'first'))
    tot['rate'] = tot.n / tot.days
    r_all, r_e11, r_int, r_e12 = (tot.loc[p, 'rate'] for p in ['all18', 'ETS2011', 'interETS', 'ETS2012'])
    ets_ratio = ((tot.loc['ETS2011', 'n'] + tot.loc['ETS2012', 'n']) / (tot.loc['ETS2011', 'days'] + tot.loc['ETS2012', 'days'])
                 / r_int if r_int > 0 else np.nan)
    nz = NR[(NR.period == 'all18')].noise_rate_per_day
    rel = {}
    for K in range(1, 6):
        p = P[P.template == K].set_index('period')
        mean = p.loc['all18', 'rate_per_day']
        rel[K] = tuple(p.loc[q, 'rate_per_day'] / mean if mean > 0 else np.nan for q in ['ETS2011', 'interETS', 'ETS2012'])
    rr = {K: rel[K][0] / rel[K][1] if rel[K][1] > 0 else np.nan for K in rel}
    k4 = max(rr, key=lambda K: -np.inf if not np.isfinite(rr[K]) else rr[K])
    base = {'variant': v, 'thr_kind': thk}
    ka = [
        {**base, 'item': 'K1 pooled rate 1 Jul 2011-1 Dec 2012 (/day, all 5)', 'paper': 8.50, 'ours': r_all,
         'pass': bool(8.5 / 2 <= r_all <= 8.5 * 2)},
        {**base, 'item': 'rate 2011 ETS (/day)', 'paper': 10.69, 'ours': r_e11, 'pass': None},
        {**base, 'item': 'rate inter-ETS 6 Sep 2011-30 Aug 2012 (/day)', 'paper': 8.24, 'ours': r_int, 'pass': None},
        {**base, 'item': 'rate 2012 ETS (/day)', 'paper': 10.26, 'ours': r_e12, 'pass': None},
        {**base, 'item': 'K2 ETS / inter-ETS rate ratio', 'paper': 1.27, 'ours': ets_ratio,
         'pass': bool(1.0 <= ets_ratio <= 1.8)},
        {**base, 'item': 'K3 noise-template rate per template per day', 'paper': '0.6 (text) / 0.79 (Fig 6)',
         'ours': nz.mean(), 'pass': bool(0.3 <= nz.mean() <= 1.6)},
        {**base, 'item': 'noise-template rate summed over 5 (/day)', 'paper': 3.94, 'ours': nz.sum(), 'pass': None},
        {**base, 'item': 'K4 template with largest 2011-ETS/inter-ETS ratio', 'paper': 'T1', 'ours': f'T{k4}',
         'pass': bool(k4 == 1)},
        {**base, 'item': 'count 1 Jul-1 Oct 2011 (all 5)', 'paper': '1394 (text run) / ~776 (Fig 6)',
         'ours': int(tot.loc['win2011', 'n']), 'pass': None},
    ]
    for K in range(1, 6):
        f = F7[F7.template == K].set_index('period').rate_rel_to_mean
        ka.append({**base, 'item': f'T{K} rel. rate ETS2011 / inter / ETS2012',
                   'paper': f"{f['ETS2011 (07-23..09-06)']:.2f} / {f['interETS_a (09-06-2011..08-30-2012)']:.2f} / "
                            f"{f['ETS2012 (08-30..10-11)']:.2f}",
                   'ours': ' / '.join(f'{x:.2f}' for x in rel[K]), 'pass': None})
    return ka


def main():
    R, NZ = load()
    TH = thresholds(R, NZ)
    TH.to_csv(f'{H.TB}/thresholds.csv', index=False)
    print(TH.round(3).to_string(index=False))
    PCs, NRs, KAs, CATS = [], [], [], {}
    for v, thk in RUNS:
        C = catalogue(R, TH, v, thr_col=thk)
        C.to_csv(f'{H.TB}/catalogue_{v}_{thk}.csv', index=False)
        if (v, thk) == ('prim', 'thrA_3rd_incl_self'):
            C.to_csv(f'{H.TB}/catalogue_prim.csv', index=False)
        CATS[(v, thk)] = C
        P = period_counts(C, f'{v}|{thk}'); P['variant'] = v; P['thr_kind'] = thk
        NR = noise_rate(NZ, TH, v, thk)
        PCs.append(P); NRs.append(NR)
        KAs += known_answer(P, NR, v, thk)
    PC, NR, KA = pd.concat(PCs), pd.concat(NRs), pd.DataFrame(KAs)
    PC.to_csv(f'{H.TB}/period_counts.csv', index=False); NR.to_csv(f'{H.TB}/noise_rate.csv', index=False)
    KA.to_csv(f'{H.TB}/known_answer.csv', index=False)
    pd.set_option('display.width', 250)
    print(KA.to_string(index=False))
    print(PC.pivot_table(index=['thr_kind', 'variant', 'period'], columns='template', values='n'))
    # best alternative = |log(rate/8.5)| smallest among prim thresholds (PLAN Amendment 2)
    k1 = KA[KA.item.str.startswith('K1') & (KA.variant == 'prim')].copy()
    k1['err'] = np.abs(np.log(np.maximum(k1.ours.astype(float), 1e-6) / 8.5))
    best = k1.sort_values('err').thr_kind.iloc[0]
    json.dump({'best_thr_kind_K1': best, 'primary': 'thrA_3rd_incl_self'}, open(f'{H.TB}/best_threshold.json', 'w'))
    print('best threshold kind for K1:', best)
    for thk in dict.fromkeys(['thrA_3rd_incl_self', best]):
        figs(CATS[('prim', thk)], NR[(NR.variant == 'prim') & (NR.thr_kind == thk)], TH, R, NZ, thk)


def figs(C, NR, TH, R, NZ, thk='thrA_3rd_incl_self'):
    TT = H.templates()
    days = pd.date_range(H.T0, H.T1 - pd.Timedelta(days=1))
    lat0, lon0 = TT.lat.mean(), TT.lon.mean()
    tr100 = tremor_daily(lat0, lon0, 100)
    daily = C.groupby('date').size().reindex(days, fill_value=0)
    bgn = NR[NR.period == 'all18'].noise_rate_per_day.sum()
    F6 = pd.read_csv(f'{H.TB}/paper_fig6_daily_digitised.csv', parse_dates=['date'])
    fig, ax = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)
    for a, (lab, y) in zip(ax, [('ours (reproduced)', daily), ('paper Fig. 6 (digitised)', None)]):
        a2 = a.twinx()
        a.bar(tr100.index, tr100.values, width=1, color='#3b4cc0', alpha=0.6, label='PNSN tremor ≤100 km')
        if y is not None:
            a2.plot(y.index, y.values, color='red', lw=0.8, label='VLFE detections/day (all 5)')
            a2.axhline(bgn, color='c', lw=1.5, label=f'noise templates (sum of 5) {bgn:.1f}/day')
        else:
            a2.plot(F6.date, F6.vlfe_per_day, color='red', lw=0.6, label='paper VLFE/day')
            a2.axhline(3.94, color='c', lw=1.5, label='paper BGN 3.9/day')
        for p in ['ETS2011', 'ETS2012']:
            a.axvspan(*pd.to_datetime(H.PER[p]), color='0.85', zorder=0)
        a.set_ylabel('tremor / day'); a2.set_ylabel('VLFE / day', color='red')
        a.set_title(lab, loc='left', fontsize=10)
        a2.legend(loc='upper right', fontsize=8)
    ax[1].set_xlim(pd.Timestamp('2011-07-01'), pd.Timestamp('2012-12-01'))
    fig.suptitle(f'H&G 2019 Fig. 6 reproduction: daily detections, five 2011 templates (all 9 stations; threshold {thk})')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/fig6_repro_{thk}.png', dpi=130); plt.close(fig)
    # Fig 7 analogue
    F7 = pd.read_csv(f'{H.TB}/paper_fig7_cumulative.csv', parse_dates=['date'])
    fig, ax = plt.subplots(5, 1, figsize=(11, 12), sharex=True)
    for K in range(1, 6):
        r = TT[TT.template == K].iloc[0]
        tr = tremor_daily(r.lat, r.lon, 100)
        a = ax[K - 1]; a2 = a.twinx()
        a.bar(tr.index, tr.values, width=1, color='#3b4cc0', alpha=0.6)
        c = C[(C.template == K) & (C.time_utc >= '2011-07-01') & (C.time_utc < '2012-12-01')].sort_values('i')
        a2.plot(c.time_utc, np.arange(1, len(c) + 1) / max(len(c), 1), color='orangered', lw=1.5, label=f'ours (n={len(c)})')
        f = F7[F7.template == K]
        a2.plot(f.date, f.cumfrac_mono, color='k', lw=1, ls='--', label='paper Fig. 7 (digitised)')
        a2.set_ylim(0, 1.02); a2.legend(loc='upper left', fontsize=8)
        a.set_ylabel(f'T{K} tremor/day')
        for p in ['ETS2011', 'ETS2012']:
            a.axvspan(*pd.to_datetime(H.PER[p]), color='0.88', zorder=0)
    ax[-1].set_xlim(pd.Timestamp('2011-07-01'), pd.Timestamp('2012-12-01'))
    fig.suptitle(f'H&G 2019 Fig. 7 reproduction ({thk}): cumulative detections per template (red) vs paper (dashed)')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/fig7_repro_{thk}.png', dpi=120); plt.close(fig)
    # separation: count above threshold vs threshold, real vs noise templates (mean CC units)
    fig, ax = plt.subplots(1, 5, figsize=(16, 3.6), sharey=True)
    a0, b0 = H.per_idx('all18')
    for K in range(1, 6):
        n = len(R[K]['prim'])
        x = np.linspace(0, 0.8, 161)
        pv, pi = R[K]['real_pv'], R[K]['real_pi']
        m = (pi >= a0) & (pi < b0) & (np.abs(pi - int(R[K]['i0'])) > H.L)
        ax[K - 1].semilogy(x, [(pv[m] / n > t).sum() / ((b0 - a0) / 86400) for t in x], 'r', label='real template')
        nz = []
        for z in NZ:
            p_i, p_v, _ = noise_peaks(z, f'T{K}_prim')
            mm = (p_i >= a0) & (p_i < b0) & (np.abs(p_i - int(z['i0'])) > H.L)
            nz.append([(p_v[mm] / n > t).sum() / ((b0 - a0) / 86400) for t in x])
        ax[K - 1].semilogy(x, np.mean(nz, axis=0), 'c', label='noise templates (mean)')
        for tk, ls in [('thrA_3rd_incl_self', ':'), ('thr_6MAD', '--'), ('thr_rate079', '-.')]:
            thr = float(TH[(TH.template == K) & (TH.variant == 'prim')][tk].iloc[0]) / n
            ax[K - 1].axvline(thr, color='k', ls=ls, lw=0.8, label=tk)
        ax[K - 1].set_title(f'T{K} ({n} ch)'); ax[K - 1].set_xlabel('summed CC / n_chan')
    ax[0].set_ylabel('detections per day'); ax[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{H.FIG}/separation_real_vs_noise_{thk}.png', dpi=120); plt.close(fig)


if __name__ == '__main__':
    main()
