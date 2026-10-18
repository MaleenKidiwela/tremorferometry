#!/usr/bin/env python
"""Amendment R2.2 held-out test. Requires tables/r2_choice.json (made by hg_r2_calib.py before any held-out count).
Counts per period with the calibrated thresholds (tables/r2_thresholds.csv), every variant and target set, pooled and per
template, vs the paper (tables/r2_targets.csv). Pass rule (pre-set): chosen variant, C1 targets, pooled inter-ETS and 2012
ETS counts both within +-25 % (Table S1 targets, Amendment R2a). Outputs tables/r2_heldout.csv, tables/r2_heldout_pass.json, figures/r2_heldout_counts.png,
figures/r2_fig6_chosen.png, figures/r2_fig7_chosen.png."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hg_common as H
import hg_r2_calib as RC
import hg_analyze as A

PERS = ['pre2011', 'ETS2011', 'interETS', 'interETS_b', 'ETS2012', 'post2012']
HELD = ['interETS', 'interETS_b', 'ETS2012', 'post2012']


def main():
    ch = json.load(open(f'{H.TB}/r2_choice.json'))
    TG = pd.read_csv(f'{H.TB}/r2_targets_S1.csv'); TG['template'] = TG.template.astype(str)
    TH = pd.read_csv(f'{H.TB}/r2_thresholds.csv')
    extra = ch.get('V7_file')
    rows = []
    cats = {}
    for v in TH.variant.unique():
        P = RC.load_variant(v, extra if v == 'V7' else None)
        for tset in TH.target_set.unique():
            th = TH[(TH.variant == v) & (TH.target_set == tset)].set_index('template').thr
            for K in range(1, 6):
                pi, pv = P[K][0], P[K][1]
                s = pv > th[K]
                if v == ch['chosen_variant'] and tset == 'S1':
                    cats[K] = pd.DataFrame({'template': K, 'i': pi[s], 'S': pv[s]})
                for p in PERS:
                    a, b = H.per_idx(p)
                    n = int((s & (pi >= a) & (pi < b)).sum())
                    tg = TG[(TG.target_set == ('S1' if tset == 'C0_commonMAD' else tset)) & (TG.template == str(K))
                            & (TG.period == p)].n_target
                    rows.append({'variant': v, 'target_set': tset, 'template': K, 'period': p, 'n_ours': n,
                                 'n_paper_template': float(tg.iloc[0]) if len(tg) else np.nan})
    D = pd.DataFrame(rows)
    pool = D.groupby(['variant', 'target_set', 'period']).n_ours.sum().reset_index()
    pp = TG[TG.target_set == 'pooled_S1'].set_index('period').n_target
    pool['n_paper'] = pool.period.map(pp)
    pool['ratio'] = pool.n_ours / pool.n_paper
    pool['within25'] = (pool.ratio >= 0.75) & (pool.ratio <= 1.25)
    D['ratio_template'] = D.n_ours / D.n_paper_template
    D.to_csv(f'{H.TB}/r2_heldout_per_template.csv', index=False)
    pool.to_csv(f'{H.TB}/r2_heldout.csv', index=False)
    pd.set_option('display.width', 250)
    print(pool.pivot_table(index=['target_set', 'variant'], columns='period', values='ratio').round(2).to_string())
    c = pool[(pool.variant == ch['chosen_variant']) & (pool.target_set == 'S1')].set_index('period')
    d = D[(D.variant == ch['chosen_variant']) & (D.target_set == 'S1')]
    ntpl = {p: int(((d.period == p) & (d.ratio_template.between(0.75, 1.25))).sum()) for p in HELD}
    res = {'chosen_variant': ch['chosen_variant'], 'spec': ch['chosen_spec'],
           **{f'{p}_ours': int(c.loc[p, 'n_ours']) for p in PERS}, **{f'{p}_paper': int(c.loc[p, 'n_paper']) for p in PERS},
           'pass_interETS': bool(c.loc['interETS', 'within25']), 'pass_ETS2012': bool(c.loc['ETS2012', 'within25']),
           'templates_within25': ntpl}
    res['PASS'] = res['pass_interETS'] and res['pass_ETS2012']
    json.dump(res, open(f'{H.TB}/r2_heldout_pass.json', 'w'), indent=1)
    print(res)
    print(d.pivot(index='template', columns='period', values='n_ours'))
    figs(pool, D, ch, cats)


def figs(pool, D, ch, cats):
    # 1. pooled counts per period: chosen (bars) vs paper (black), other variants (grey dots)
    fig, ax = plt.subplots(1, 1, figsize=(10, 4.2))
    x = np.arange(len(PERS))
    c1 = pool[pool.target_set == 'S1']
    for v, g in c1.groupby('variant'):
        g = g.set_index('period').reindex(PERS)
        if v == ch['chosen_variant']:
            ax.bar(x + 0.2, g.n_ours, 0.4, color='r', label=f'ours, chosen {v} ({ch["chosen_spec"]})')
        else:
            ax.plot(x + 0.2, g.n_ours, 'o', color='0.6', ms=4)
            ax.text(x[-1] + 0.45, g.n_ours.iloc[-1], v, fontsize=6, color='0.4')
    ax.plot([], [], 'o', color='0.6', label='other variants (thresholds fitted to Table S1)')
    g = c1[c1.variant == ch['chosen_variant']].set_index('period').reindex(PERS)
    ax.bar(x - 0.2, g.n_paper, 0.4, color='k', label='paper (Table S1)')
    ax.axvspan(-0.5, 1.5, color='#dde8ff', zorder=0, label='calibration period')
    ax.set_xticks(x); ax.set_xticklabels(['pre-ETS\n1-23 Jul 11', '2011 ETS', 'inter-ETS\n6 Sep 11-30 Aug 12',
                                          'inter-ETS b\n1 Oct-30 Aug', '2012 ETS', 'post\n11 Oct-1 Dec 12'])
    ax.set_yscale('log'); ax.set_ylabel('detections (five templates)'); ax.legend(fontsize=7, loc='upper right')
    ax.set_title('R2: thresholds fitted to the 2011-ETS count per template; everything right of the blue band is held out',
                 fontsize=9)
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2_heldout_counts.png', dpi=130); plt.close(fig)
    # 2. Fig 6 / Fig 7 style for the chosen setting
    C = pd.concat(cats.values(), ignore_index=True)
    C['time_utc'] = H.itime(C.i.values); C['date'] = C.time_utc.dt.floor('D')
    C.to_csv(f'{H.TB}/r2_catalogue_chosen.csv', index=False)
    TT = H.templates()
    tr = A.tremor_daily(TT.lat.mean(), TT.lon.mean(), 100)
    days = pd.date_range('2011-07-01', '2012-12-01', inclusive='left')
    daily = C.groupby('date').size().reindex(days, fill_value=0)
    S1 = pd.read_csv(f'{H.TB}/hg2019_tableS1.csv', parse_dates=['t'])
    s1d = S1.groupby(S1.t.dt.floor('D')).size().reindex(days, fill_value=0)
    fig, ax = plt.subplots(2, 1, figsize=(12, 6.5), sharex=True)
    for a, lab in zip(ax, ['ours (R2a chosen setting, calibrated on 1 Jul-6 Sep 2011 only)', 'paper Table S1']):
        a2 = a.twinx()
        a.bar(tr.index, tr.values, width=1, color='#3b4cc0', alpha=0.5)
        if a is ax[0]:
            a2.plot(daily.index, daily.values, 'r', lw=0.8)
        else:
            a2.plot(s1d.index, s1d.values, 'r', lw=0.8)
        a2.set_ylim(0, 30); a2.set_ylabel('VLFE / day', color='r'); a.set_ylabel('tremor / day ≤100 km')
        a.axvspan(pd.Timestamp('2011-07-01'), pd.Timestamp('2011-09-06'), color='#dde8ff', zorder=0)
        a.set_title(lab, loc='left', fontsize=10)
    ax[1].set_xlim(days[0], days[-1])
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2_fig6_chosen.png', dpi=130); plt.close(fig)
    F7 = pd.read_csv(f'{H.TB}/paper_fig7_cumulative.csv', parse_dates=['date'])
    fig, ax = plt.subplots(5, 1, figsize=(11, 11), sharex=True)
    for K in range(1, 6):
        c = C[(C.template == K) & (C.time_utc >= '2011-07-01') & (C.time_utc < '2012-12-01')].sort_values('i')
        ax[K - 1].plot(c.time_utc, np.arange(1, len(c) + 1) / max(len(c), 1), color='orangered', lw=1.5, label=f'ours (n={len(c)})')
        f = S1[S1.template == K].sort_values('t')
        ax[K - 1].plot(f.t, np.arange(1, len(f) + 1) / max(len(f), 1), 'k--', lw=1, label=f'paper Table S1 (n={len(f)})')
        ax[K - 1].axvspan(pd.Timestamp('2011-07-01'), pd.Timestamp('2011-09-06'), color='#dde8ff', zorder=0)
        ax[K - 1].set_ylabel(f'T{K} cumulative'); ax[K - 1].legend(fontsize=7, loc='upper left')
    ax[-1].set_xlim(days[0], days[-1])
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2_fig7_chosen.png', dpi=120); plt.close(fig)


if __name__ == '__main__':
    main()
