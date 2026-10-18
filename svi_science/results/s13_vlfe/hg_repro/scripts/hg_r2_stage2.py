#!/usr/bin/env python
"""Amendment R2.4: at the chosen (calibrated) setting — H&G-style noise templates, our scrambled controls, far/near flag,
survival, LFE coincidence, ETS vs inter-ETS, and the inter-ETS excess above the noise templates.
Inputs: tables/r2_choice.json, tables/r2_thresholds.csv (Table-S1 thresholds, target_set S1, for the chosen variant; Amendment R2a); real + control peaks and
noise-template peaks for the chosen spec (R1 files data/mf/real_T*.npz, noise_*.npz when the spec is V0 — same code path;
else data/r2/real_<tag>_T*.npz, noise_<tag>_*.npz).
Outputs tables/r2_excess.csv, r2_controls.csv, r2_lfe.csv, r2_noise_counts.csv; figures/r2_controls_daily.png,
r2_excess_hist.png, r2_lfe.png."""
import os, sys, json, glob
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hg_common as H
import hg_r2 as R
import hg_r2_calib as RC
import hg_analyze as A
import hg_controls as C
import hg_lfe as LF

NDRAW = 8
PERS = ['ETS2011', 'interETS', 'interETS_b', 'ETS2012', 'all18']


def load(ch):
    spec = ch['chosen_spec']
    base, mask = spec.split('+')[0], spec.endswith('+mask')
    if base == RC.SPEC['V0'] and ch['chosen_variant'] in ('V0', 'V8'):
        real = {K: dict(np.load(f'{H.MF}/real_T{K}.npz', allow_pickle=True)) for K in range(1, 6)}
        nz = []
        for f in sorted(glob.glob(f'{H.MF}/noise_[0-9][0-9][0-9].npz')):
            z = dict(np.load(f, allow_pickle=True))
            d = {'i0': int(z['i0']), **{K: A.noise_peaks(z, f'T{K}_prim')[:2] for K in range(1, 6)}}
            for q in range(int(z['ncfg'])):                 # R1 noise peaks were stored above median + 3 MAD
                for nm in z[f'c{q}_names']:
                    if str(nm).endswith('_prim'):
                        d[f'floor{str(nm)[1]}'] = float(z[f'c{q}_med'] + 3 * z[f'c{q}_mad'])
            nz.append(d)
        src = 'R1 files (identical code path for V0)'
    else:
        tg = R.tag(base)
        real = {K: dict(np.load(f'{R.R2}/real_{tg}_T{K}.npz', allow_pickle=True)) for K in range(1, 6)}
        nz = []
        for f in sorted(glob.glob(f'{R.R2}/noise_{tg}_[0-9][0-9][0-9].npz')):
            z = dict(np.load(f, allow_pickle=True))
            d = {'i0': int(z['i0'])}
            for q in range(int(z['ncfg'])):
                for K in z[f'c{q}_templates']:
                    d[int(K)] = (z[f'c{q}_pi'], z[f'c{q}_pv'])
            nz.append(d)
        src = f'data/r2/*_{tg}_*'
    return real, nz, mask, src


def main():
    ch = json.load(open(f'{H.TB}/r2_choice.json'))
    TH = pd.read_csv(f'{H.TB}/r2_thresholds.csv')
    th = TH[(TH.variant == ch['chosen_variant']) & (TH.target_set == 'S1')].set_index('template').thr
    real, NZ, mask, src = load(ch)
    print('source', src, 'noise templates', len(NZ), 'thresholds', th.round(3).to_dict())
    fl = [(K, d[f'floor{K}']) for d in NZ for K in range(1, 6) if f'floor{K}' in d and d[f'floor{K}'] > th[K]]
    print(f'noise templates whose stored-peak floor exceeds the threshold: {len(fl)} of {5 * len(NZ)} (must be 0)')
    keepm = (lambda i: ~RC.eqmask(i)) if mask else (lambda i: np.ones(len(i), bool))
    # ---- catalogues: real, A1_r, A2_r
    rows = []
    for K in range(1, 6):
        z = real[K]; i0 = int(z['i0'])
        for s in ['real'] + [f'A1_{r}' for r in range(NDRAW)] + [f'A2_{r}' for r in range(NDRAW)]:
            pi, pv = z[f'{s}_pi'], z[f'{s}_pv']
            sel = pv > th[K]
            if s == 'real':
                sel &= np.abs(pi - i0) > H.L
            sel &= keepm(pi)
            rows.append(pd.DataFrame({'template': K, 'set': s, 'i': pi[sel], 'S': pv[sel]}))
    AC = pd.concat(rows, ignore_index=True)
    fl = []
    for K in range(1, 6):
        m = AC.template == K
        chans = real[K]['chans'] if 'chans' in real[K] else real[K]['prim']
        f = C.farnear(AC[m].i.values, C.zstas_of(chans)); f.index = AC[m].index; fl.append(f)
    AC = AC.join(pd.concat(fl))
    AC['time_utc'] = H.itime(AC.i.values)
    AC.to_parquet(f'{R.R2}/r2_controls_catalogue.parquet')
    # check: real 2011-ETS counts equal the C1 targets
    a, b = H.per_idx('ETS2011')
    print('2011 ETS real counts', AC[(AC.set == 'real') & (AC.i >= a) & (AC.i < b)].groupby('template').size().to_dict())
    # ---- noise templates
    nrows = []
    for j, d in enumerate(NZ):
        for K in range(1, 6):
            pi, pv = d[K]
            s = (pv > th[K]) & (np.abs(pi - d['i0']) > H.L)
            s &= keepm(pi)
            for p in PERS:
                a, b = H.per_idx(p)
                nrows.append({'j': j, 'template': K, 'period': p, 'n': int((s & (pi >= a) & (pi < b)).sum())})
    NC = pd.DataFrame(nrows); NC.to_csv(f'{H.TB}/r2_noise_counts.csv', index=False)
    # ---- excess above noise, survival vs scrambled controls
    ex, cr = [], []
    for p in PERS:
        a, b = H.per_idx(p); nd = (b - a) / 86400
        d = AC[(AC.i >= a) & (AC.i < b)]
        for K in list(range(1, 6)) + ['all']:
            dd = d if K == 'all' else d[d.template == K]
            nr = int((dd.set == 'real').sum()); nru = int(((dd.set == 'real') & ~dd.flag_far).sum())
            nn = NC[(NC.period == p) & ((NC.template == K) if K != 'all' else True)].groupby('j').n.sum().values
            e = (nr - nn) / nr if nr else np.full(nn.shape, np.nan)
            a1 = np.array([((dd.set == f'A1_{r}')).sum() for r in range(NDRAW)])
            a2 = np.array([((dd.set == f'A2_{r}')).sum() for r in range(NDRAW)])
            a1u = np.array([((dd.set == f'A1_{r}') & ~dd.flag_far).sum() for r in range(NDRAW)])
            a2u = np.array([((dd.set == f'A2_{r}') & ~dd.flag_far).sum() for r in range(NDRAW)])
            ntpl = 5 if K == 'all' else 1
            ex.append({'period': p, 'template': K, 'days': nd, 'N_real': nr,
                       'real_per_tpl_day': nr / nd / ntpl, 'noise_per_tpl_day_mean': nn.mean() / nd / ntpl,
                       'noise_N_mean': nn.mean(), 'noise_N_p16': np.percentile(nn, 16), 'noise_N_p84': np.percentile(nn, 84),
                       'excess_frac_mean': (nr - nn.mean()) / nr if nr else np.nan,
                       'excess_frac_p16': np.percentile(e, 16), 'excess_frac_p84': np.percentile(e, 84),
                       'sqrtN_real_frac': np.sqrt(nr) / nr if nr else np.nan,
                       'frac_noise_templates_ge_real': float((nn >= nr).mean())})
            surv = max(0.0, nru - max(a1u.mean(), a2u.mean()))
            cr.append({'period': p, 'template': K, 'N_real': nr, 'N_real_unflagged': nru,
                       'A1_mean': a1.mean(), 'A1_range': f'{a1.min()}-{a1.max()}', 'A2_mean': a2.mean(),
                       'A2_range': f'{a2.min()}-{a2.max()}', 'surv_frac_scrambled': surv / nr if nr else np.nan,
                       'surv_frac_noise_templates': max(0.0, nru - nn.mean()) / nr if nr else np.nan,
                       'real_flag_frac': 1 - nru / nr if nr else np.nan})
    EX, CR = pd.DataFrame(ex), pd.DataFrame(cr)
    EX.to_csv(f'{H.TB}/r2_excess.csv', index=False); CR.to_csv(f'{H.TB}/r2_controls.csv', index=False)
    pd.set_option('display.width', 250)
    print(EX.round(3).to_string(index=False)); print(CR.round(3).to_string(index=False))
    # ---- LFE coincidence (+-60 s): real, scrambled, noise-template detections, random times
    L = LF.lfe_times()
    rng = np.random.default_rng(H.SEED + 31)
    lr = []
    for p in ['ETS2011', 'interETS', 'ETS2012']:
        a, b = H.per_idx(p)
        d = AC[(AC.i >= a) & (AC.i < b)]
        ni = []
        for j, dz in enumerate(NZ):
            for K in range(1, 6):
                pi, pv = dz[K]
                s = (pv > th[K]) & (np.abs(pi - dz['i0']) > H.L) & (pi >= a) & (pi < b)
                s &= keepm(pi)
                ni.append(pi[s])
        ni = np.concatenate(ni)
        for cat in ['both', 'SVI', 'WA_NW']:
            Lc = L[cat]
            sets = {'real': d[d.set == 'real'].i.values, 'polarity-scrambled': d[d.set.str.startswith('A2')].i.values,
                    'station-day-shift': d[d.set.str.startswith('A1')].i.values, 'noise templates': ni,
                    'random times': rng.integers(a, b, 20000)}
            for k, t in sets.items():
                c = LF.counts(Lc, t)
                lr.append({'period': p, 'lfe_cat': cat, 'set': k, 'n': len(t),
                           'frac_with_lfe': (c > 0).mean() if len(c) else np.nan, 'mean_lfe': c.mean() if len(c) else np.nan})
    LR = pd.DataFrame(lr)
    rnd = LR[LR.set == 'random times'].set_index(['period', 'lfe_cat']).mean_lfe
    LR['ratio_vs_random'] = [r.mean_lfe / rnd[(r.period, r.lfe_cat)] if rnd[(r.period, r.lfe_cat)] > 0 else np.nan
                             for r in LR.itertuples()]
    LR.to_csv(f'{H.TB}/r2_lfe.csv', index=False)
    print(LR[LR.lfe_cat == 'both'].round(4).to_string(index=False))
    figs(AC, NC, EX, LR, ch)


def figs(AC, NC, EX, LR, ch):
    days = pd.date_range('2011-07-01', '2012-12-01', inclusive='left')
    a18, b18 = H.per_idx('all18')
    TT = H.templates(); tr = A.tremor_daily(TT.lat.mean(), TT.lon.mean(), 100)
    def daily(sets, n):
        x = AC[AC.set.isin(sets) & (AC.i >= a18) & (AC.i < b18)]
        return x.groupby(x.time_utc.dt.floor('D')).size().reindex(days, fill_value=0) / n
    nzr = EX[(EX.template == 'all') & (EX.period == 'all18')].noise_per_tpl_day_mean.iloc[0] * 5
    fig, ax = plt.subplots(2, 1, figsize=(12, 7), sharex=True)
    ax[0].plot(days, daily(['real'], 1), 'r', lw=0.8, label='real templates (calibrated H&G-style)')
    ax[0].plot(days, daily([f'A2_{r}' for r in range(NDRAW)], NDRAW), color='purple', lw=0.8, label='polarity-scrambled (mean of 8)')
    ax[0].plot(days, daily([f'A1_{r}' for r in range(NDRAW)], NDRAW), color='green', lw=0.8, label='station-day-shift (mean of 8)')
    ax[0].axhline(nzr, color='c', lw=1.5, label=f'H&G-style noise templates (mean of 100) {nzr:.1f}/day')
    ax[0].axhline(3.94, color='c', lw=1, ls='--', label='paper noise background 3.9/day')
    ax[0].set_ylabel('detections/day (5 templates)'); ax[0].legend(fontsize=7, ncol=2)
    ax[1].bar(tr.index, tr.values, width=1, color='#3b4cc0'); ax[1].set_ylabel('PNSN tremor/day ≤100 km')
    for g in ax:
        for p in ['ETS2011', 'ETS2012']:
            g.axvspan(*pd.to_datetime(H.PER[p]), color='0.88', zorder=0)
    ax[1].set_xlim(days[0], days[-1])
    fig.suptitle(f'R2 chosen setting {ch["chosen_variant"]} ({ch["chosen_spec"]}): real vs noise templates vs scrambled controls')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2_controls_daily.png', dpi=130); plt.close(fig)
    # raw distribution of noise-template counts vs the real count, per period (pooled 5 templates)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6))
    for g, p in zip(ax, ['ETS2011', 'interETS', 'ETS2012']):
        nn = NC[NC.period == p].groupby('j').n.sum().values
        nr = EX[(EX.template == 'all') & (EX.period == p)].N_real.iloc[0]
        g.hist(nn, bins=30, color='c', label='100 noise templates (each: 5 template configs summed)')
        g.axvline(nr, color='r', lw=2, label=f'real templates {nr}')
        g.set_title(p); g.set_xlabel('detections in period')
    ax[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2_excess_hist.png', dpi=130); plt.close(fig)
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.6), sharey=True)
    for g, p in zip(ax, ['ETS2011', 'interETS', 'ETS2012']):
        x = LR[(LR.period == p) & (LR.lfe_cat == 'both')]
        g.bar(np.arange(len(x)), x.frac_with_lfe, color=['r', 'purple', 'green', 'c', '0.6'])
        g.set_xticks(np.arange(len(x))); g.set_xticklabels(x.set, rotation=30, fontsize=7); g.set_title(p)
    ax[0].set_ylabel('fraction with ≥1 of our LFEs within ±60 s')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2_lfe.png', dpi=130); plt.close(fig)


if __name__ == '__main__':
    main()
