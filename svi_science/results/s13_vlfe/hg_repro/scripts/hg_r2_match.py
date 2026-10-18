#!/usr/bin/env python
"""Amendment R2a.4: detection-by-detection match of our calibrated catalogue (chosen setting, Table-S1 thresholds) vs
Table S1, per template and period, at +-30 s, raw +-5 s and +-5 s after removing one constant offset (median ours - theirs
over the +-30 s matches in the calibration period 1 Jul - 6 Sep 2011). Matched vs unmatched Table-S1 detections compared
on our summed CC at their time (largest of our declustered peaks within +-45 s, in MAD units above the median), Pacific
hour of day, period, month; random times give the baseline CC.
Outputs tables/r2_match.csv, tables/r2_match_S1_detail.csv, tables/r2_match_offset.json, figures/r2_match.png."""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hg_common as H
import hg_r2_calib as RC

PERS = ['pre2011', 'ETS2011', 'interETS', 'interETS_b', 'ETS2012', 'post2012', 'all18']
PCOL = {'pre2011': '0.5', 'ETS2011': 'r', 'interETS': 'b', 'ETS2012': 'm', 'post2012': 'g'}


def nearest(x, ref):
    """signed distance from each x to the nearest ref (ref sorted)."""
    if ref.size == 0:
        return np.full(len(x), np.inf)
    j = np.clip(np.searchsorted(ref, x), 1, ref.size - 1) if ref.size > 1 else np.zeros(len(x), int)
    if ref.size == 1:
        return x - ref[0]
    d1, d0 = x - ref[j], x - ref[j - 1]
    return np.where(np.abs(d1) < np.abs(d0), d1, d0)


def per_name(i):
    out = np.full(len(i), 'other', object)
    for p in ['pre2011', 'ETS2011', 'interETS', 'ETS2012', 'post2012']:
        a, b = H.per_idx(p)
        out[(i >= a) & (i < b)] = p
    return out


def main():
    ch = json.load(open(f'{H.TB}/r2_choice.json'))
    C = pd.read_csv(f'{H.TB}/r2_catalogue_chosen.csv')
    S = pd.read_csv(f'{H.TB}/hg2019_tableS1.csv')
    P = RC.load_variant(ch['chosen_variant'], ch.get('V7_file') if ch['chosen_variant'] == 'V7' else None)
    a0, a1 = H.tidx('2011-07-01'), RC.CAL_END
    # offset from the calibration period
    offs = []
    for K in range(1, 6):
        ours = np.sort(C[C.template == K].i.values)
        th = S[(S.template == K) & (S.i >= a0) & (S.i < a1)].i.values
        d = -nearest(th, ours)                      # ours - theirs
        offs.append(d[np.abs(d) <= 30])
    offs = np.concatenate(offs)
    off0 = float(np.median(offs)) if offs.size else 0.0
    json.dump({'offset_ours_minus_theirs_s': off0, 'n_matches_used': int(offs.size),
               'p16': float(np.percentile(offs, 16)), 'p84': float(np.percentile(offs, 84))},
              open(f'{H.TB}/r2_match_offset.json', 'w'), indent=1)
    print('offset ours - theirs (s): median', off0, 'n', offs.size, 'p16/p84', np.percentile(offs, [16, 84]))
    rows, det = [], []
    for K in range(1, 6):
        ours = np.sort(C[C.template == K].i.values)
        allS = S[S.template == K].copy()
        th_all = allS.i.values
        d_all = -nearest(th_all, ours)
        pi, pv, med, mad = P[K]
        o = np.argsort(pi); pi, pv = pi[o], pv[o]
        lo = np.searchsorted(pi, th_all - 45); hi = np.searchsorted(pi, th_all + 45, side='right')
        ccv = np.array([pv[a:b].max() if b > a else med for a, b in zip(lo, hi)])
        allS['cc_mad_ours'] = (ccv - med) / mad
        allS['d_nearest_ours_s'] = d_all
        allS['matched30'] = np.abs(d_all) <= 30
        allS['period'] = per_name(th_all)
        det.append(allS)
        for p in PERS:
            a, b = H.per_idx(p)
            th = th_all[(th_all >= a) & (th_all < b)]
            ou = ours[(ours >= a) & (ours < b)]
            dt = -nearest(th, ou) if ou.size else np.full(th.size, np.inf)
            do = nearest(ou, np.sort(th)) if th.size else np.full(ou.size, np.inf)
            r = {'template': K, 'period': p, 'n_theirs': th.size, 'n_ours': ou.size}
            for lab, w, sh in [('30', 30, 0.0), ('5raw', 5, 0.0), ('5off', 5, off0)]:
                r[f'recall_{lab}'] = float((np.abs(dt - sh) <= w).mean()) if th.size else np.nan
                r[f'precision_{lab}'] = float((np.abs(do - sh) <= w).mean()) if ou.size else np.nan
            rows.append(r)
    M = pd.DataFrame(rows)
    # pooled rows
    for p in PERS:
        m = M[M.period == p]
        r = {'template': 'all', 'period': p, 'n_theirs': m.n_theirs.sum(), 'n_ours': m.n_ours.sum()}
        for lab in ['30', '5raw', '5off']:
            r[f'recall_{lab}'] = (m[f'recall_{lab}'] * m.n_theirs).sum() / max(m.n_theirs.sum(), 1)
            r[f'precision_{lab}'] = (m[f'precision_{lab}'] * m.n_ours).sum() / max(m.n_ours.sum(), 1)
        rows.append(r)
    M = pd.DataFrame(rows)
    M.to_csv(f'{H.TB}/r2_match.csv', index=False)
    D = pd.concat(det, ignore_index=True)
    D['t'] = H.itime(D.i.values)
    D['hour_PT'] = D.t.dt.tz_localize('UTC').dt.tz_convert('America/Los_Angeles').dt.hour
    D['month'] = D.t.dt.strftime('%Y-%m')
    D.to_csv(f'{H.TB}/r2_match_S1_detail.csv', index=False)
    pd.set_option('display.width', 250)
    print(M.round(3).to_string(index=False))
    print(D.groupby(['period', 'matched30']).cc_mad_ours.describe()[['count', '25%', '50%', '75%']].round(2))
    # random-time baseline of our CC (MAD units)
    rng = np.random.default_rng(H.SEED + 41)
    base = []
    for K in range(1, 6):
        pi, pv, med, mad = P[K]
        o = np.argsort(pi); pi, pv = pi[o], pv[o]
        t = rng.integers(H.tidx('2011-07-01'), H.tidx('2012-12-01'), 3000)
        lo = np.searchsorted(pi, t - 45); hi = np.searchsorted(pi, t + 45, side='right')
        base.append((np.array([pv[a:b].max() if b > a else med for a, b in zip(lo, hi)]) - med) / mad)
    base = np.concatenate(base)
    print('random-time CC (MAD units) median', np.median(base))
    # figure
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.2))
    x = np.arange(6)
    for j, p in enumerate(['pre2011', 'ETS2011', 'interETS', 'ETS2012', 'post2012']):
        m = M[(M.period == p)].set_index('template').reindex([1, 2, 3, 4, 5, 'all'])
        ax[0].bar(x + (j - 2) * 0.16, m.recall_30, 0.16, color=PCOL[p], label=p)
    ax[0].set_xticks(x); ax[0].set_xticklabels(['T1', 'T2', 'T3', 'T4', 'T5', 'all'])
    ax[0].set_ylabel('fraction of Table-S1 detections we match (±30 s)'); ax[0].legend(fontsize=7)
    ax[0].axvspan(-0.5, 5.5, ymin=0, ymax=0, color='w')
    bins = np.linspace(-2, 14, 65)
    ax[1].hist(base, bins, histtype='step', color='0.5', density=True, label='random times')
    for mt, ls in [(True, '-'), (False, '--')]:
        for p in ['ETS2011', 'interETS']:
            v = D[(D.matched30 == mt) & (D.period == p)].cc_mad_ours
            ax[1].hist(v, bins, histtype='step', color=PCOL[p], ls=ls, density=True,
                       label=f'Table S1 {p}, {"matched" if mt else "unmatched"} (n={len(v)})')
    ax[1].set_xlabel('our summed CC at their time (MAD above median)'); ax[1].legend(fontsize=6)
    for mt, c in [(True, 'k'), (False, 'orange')]:
        v = D[D.matched30 == mt].hour_PT
        ax[2].plot(np.arange(24), np.bincount(v, minlength=24) / max(len(v), 1), color=c, marker='o', ms=3,
                   label=f'Table S1 {"matched" if mt else "unmatched"} (n={len(v)})')
    ax[2].set_xlabel('hour of day (PT)'); ax[2].set_ylabel('fraction'); ax[2].legend(fontsize=7)
    fig.suptitle(f'Our calibrated catalogue ({ch["chosen_variant"]}: {ch["chosen_spec"]}) vs Table S1, detection by detection')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2_match.png', dpi=130); plt.close(fig)


if __name__ == '__main__':
    main()
