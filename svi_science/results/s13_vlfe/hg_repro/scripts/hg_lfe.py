#!/usr/bin/env python
"""Coordinator request 2026-10-15 (after CONTROLS/STACKINV): (1) do H&G-style detections carry LFEs? Our LFE catalogues
within +-60 s of each detection time vs random times in the same period and vs the scrambled-control detections.
LFEs: SVI rows (s13_common.svi_rows, data series, counted once, P >= 0.5) + WA-NW (catalogue_v2, selected_v2, region NW;
detection window times, no origin correction - a few s, small vs 60 s). 1 Jun 2011 - 31 Dec 2012.
(2) selection-bias check on the stacks: waveform CC (0-90 s after the origin, all channels concatenated) between each
stack and the window it was selected with (real stacks vs the template; noise-template stacks vs their own noise
window; random stacks vs the template).
Outputs tables/lfe_coincidence_{thk}.csv, tables/stack_selfsimilarity_{thk}.csv, figures/lfe_coincidence.png."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/scripts')
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hg_common as H

W = 60
NDRAW = 8


def lfe_times():
    import s13_common as C
    R = C.svi_rows()
    R = R[(R.series == 'data') & R.counted_once & (R.P >= 0.5)]
    t_svi = R.t.values
    c = pd.read_parquet(f'{H.ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet',
                        columns=['time', 'region', 'selected_v2'])
    c = c[c.selected_v2 & (c.region == 'NW')]
    t_wa = c.time.astype('datetime64[ns]').astype('int64').values / 1e9
    t0 = (H.T0 - pd.Timestamp('1970-01-01')).total_seconds()
    out = {}
    for k, t in (('SVI', t_svi), ('WA_NW', t_wa)):
        i = np.sort(np.round(t - t0).astype(np.int64))
        out[k] = i[(i >= 0) & (i < H.tidx('2013-01-01'))]
    out['both'] = np.sort(np.concatenate([out['SVI'], out['WA_NW']]))
    return out


def counts(L, t):
    t = np.asarray(t, np.int64)
    return np.searchsorted(L, t + W, side='right') - np.searchsorted(L, t - W, side='left')


def lfe_table(thk, LF, rng):
    AC = pd.read_parquet(f'{H.MF}/controls_catalogue_{thk}.parquet')
    rows = []
    for p in ['ETS2011', 'interETS', 'ETS2012', 'all18']:
        a, b = H.per_idx(p)
        for K in list(range(1, 6)) + ['all']:
            d = AC[(AC.i >= a) & (AC.i < b)] if K == 'all' else AC[(AC.template == K) & (AC.i >= a) & (AC.i < b)]
            for cat, L in LF.items():
                r = {'thr_kind': thk, 'period': p, 'template': K, 'lfe_cat': cat}
                re_ = d[d.set == 'real'].i.values
                rnd = rng.integers(a, b, 20000)
                cr, cn = counts(L, re_), counts(L, rnd)
                c2 = counts(L, d[d.set.str.startswith('A2')].i.values)
                c1 = counts(L, d[d.set.str.startswith('A1')].i.values)
                r.update({'n_det': len(re_), 'frac_det_with_lfe': (cr > 0).mean() if len(cr) else np.nan,
                          'frac_random_with_lfe': (cn > 0).mean(),
                          'frac_A2_with_lfe': (c2 > 0).mean() if len(c2) else np.nan,
                          'frac_A1_with_lfe': (c1 > 0).mean() if len(c1) else np.nan,
                          'mean_lfe_det': cr.mean() if len(cr) else np.nan, 'mean_lfe_random': cn.mean(),
                          'mean_lfe_A2': c2.mean() if len(c2) else np.nan})
                r['ratio_det_vs_random'] = r['mean_lfe_det'] / r['mean_lfe_random'] if r['mean_lfe_random'] > 0 else np.nan
                rows.append(r)
    return pd.DataFrame(rows)


def selfsim(thk):
    z = np.load(f'{H.MF}/stacks_{thk}.npz')
    rows = []
    for K in range(1, 6):
        if f'T{K}_template' not in z.files:
            continue
        tpl = z[f'T{K}_template'][:, 205:295]
        def cc(a, b):
            ok = np.all(np.isfinite(a), 1) & np.all(np.isfinite(b), 1)
            x, y = a[ok].ravel(), b[ok].ravel()
            x = x - x.mean(); y = y - y.mean()
            return float(x @ y / np.sqrt((x @ x) * (y @ y)))
        for key in z.files:
            if not key.startswith(f'T{K}_') or key.endswith('_single') or key == f'T{K}_template':
                continue
            st = z[key][:, 205:295]
            name = key.split('_', 1)[1]
            if name.startswith('noise'):
                ref = z[f'T{K}_{name.replace("_stack", "_single")}'][:, 205:295]
                rows.append({'thr_kind': thk, 'template': K, 'set': name, 'cc_with_selecting_window': cc(st, ref),
                             'cc_with_template': cc(st, tpl)})
            else:
                rows.append({'thr_kind': thk, 'template': K, 'set': name, 'cc_with_selecting_window':
                             cc(st, tpl) if not name.startswith('random') else np.nan, 'cc_with_template': cc(st, tpl)})
    return pd.DataFrame(rows)


def main():
    LF = lfe_times()
    print({k: len(v) for k, v in LF.items()})
    rng = np.random.default_rng(H.SEED + 21)
    outs, sims = [], []
    for thk in ['thrA_3rd_incl_self', 'thr_6MAD', 'thr_rate079']:
        T = lfe_table(thk, LF, rng); T.to_csv(f'{H.TB}/lfe_coincidence_{thk}.csv', index=False); outs.append(T)
        if os.path.isfile(f'{H.MF}/stacks_{thk}.npz'):
            S = selfsim(thk); S.to_csv(f'{H.TB}/stack_selfsimilarity_{thk}.csv', index=False); sims.append(S)
    pd.set_option('display.width', 250)
    T = pd.concat(outs)
    print(T[(T.template == 'all') & (T.lfe_cat == 'both')].round(3).to_string(index=False))
    for S in sims:
        print(S.round(2).to_string(index=False))
    # figure: fraction of detections with >=1 LFE within 60 s, real vs A2 vs random, per period (thr_rate079, both cats)
    t = T[(T.thr_kind == 'thr_rate079') & (T.lfe_cat == 'both') & (T.template != 'all')]
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    for g, p in zip(ax, ['ETS2011', 'interETS', 'ETS2012']):
        x = t[t.period == p]
        k = np.arange(len(x))
        g.bar(k - 0.27, x.frac_det_with_lfe, 0.27, color='r', label='H&G-style detections')
        g.bar(k, x.frac_A2_with_lfe, 0.27, color='purple', label='polarity-scrambled detections')
        g.bar(k + 0.27, x.frac_random_with_lfe, 0.27, color='0.6', label='random times')
        g.set_xticks(k); g.set_xticklabels([f'T{v}' for v in x.template]); g.set_title(p)
    ax[0].set_ylabel('fraction with ≥1 of our LFEs within ±60 s'); ax[0].legend(fontsize=7)
    fig.suptitle('Do H&G-style VLFE detections (threshold thr_rate079) carry LFEs? SVI + WA-NW catalogues')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/lfe_coincidence.png', dpi=120); plt.close(fig)


if __name__ == '__main__':
    main()
