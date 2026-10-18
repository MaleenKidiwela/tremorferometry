#!/usr/bin/env python
"""Amendment R2.2: calibrate per-template thresholds on the 2011 ETS count and choose the processing variant on the
calibration period (1 Jul - 6 Sep 2011) ONLY. Writes tables/r2_calibration.csv, tables/r2_thresholds.csv,
tables/r2_choice.json and appends the choice to NOTES.md. Does not count anything after 6 Sep 2011.
usage: hg_r2_calib.py [extra.npz-spec ...]   (extra = V7 combination files data/r2/real_<tag>_T*.npz, named V7)"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
import hg_common as H

R2 = f'{H.HG}/data/r2'
SPEC = {'V0': 'disp,0,sum,all,90', 'V1': 'disp,0,mean,all,90', 'V2': 'disp,0,net,all,90', 'V3a': 'disp,0,sum,snr15,90',
        'V3b': 'disp,0,sum,near5,90', 'V4': 'vel,0,sum,all,90', 'V5': 'disp,-30,sum,all,90', 'V5b': 'disp,-9,sum,all,90',
        'V6': 'disp,0,sum,all,30', 'V6b': 'disp,0,sum,all,20',
        'V8': 'disp,0,sum,all,90+mask'}
GROUP = {'V1': 'norm', 'V2': 'norm', 'V3a': 'chset', 'V3b': 'chset', 'V4': 'other', 'V5': 'other', 'V5b': 'other', 'V6': 'other', 'V6b': 'other',
         'V8': 'other'}
CAL_END = H.tidx('2011-09-06')
_M = {}


def eqmask(i):
    """True where a detection window [t, t+90] touches our v2 earthquake masks (2011-2012 only)."""
    import hg_controls as C
    i = np.asarray(i, np.int64)
    return C.mask(i, 'v2') | C.mask(i + 45, 'v2') | C.mask(np.minimum(i + 90, H.N - 1), 'v2')


def load_variant(v, extra=None):
    """dict K -> (pi, pv, med, mad), self-detections removed (origin-aligned indices)."""
    out = {}
    for K in range(1, 6):
        if extra:
            z = np.load(extra.split('+')[0].replace('TK', f'T{K}'))
            pi, pv, med, mad, i0 = z['real_pi'], z['real_pv'], float(z['real_med']), float(z['real_mad']), int(z['i0'])
        else:
            z = np.load(f'{R2}/var_T{K}.npz')
            t = 'V0' if v == 'V8' else v
            pi, pv, med, mad, i0 = z[f'{t}_pi'], z[f'{t}_pv'], float(z[f'{t}_med']), float(z[f'{t}_mad']), int(z['i0'])
        keep = np.abs(pi - i0) > H.L
        if v == 'V8' or (extra and extra.endswith('+mask')):
            keep &= ~eqmask(pi)
        out[K] = (pi[keep], pv[keep], med, mad)
    return out


def thr_for_count(vals, n):
    v = np.sort(vals)[::-1]
    n = int(round(n))
    if n <= 0:
        return float(v[0]) + 1e-6
    if n >= v.size:
        return float(v[-1]) - 1e-6
    return 0.5 * (float(v[n - 1]) + float(v[n]))


def thresholds(P, tset, TG):
    a, b = H.per_idx('ETS2011')
    if tset in ('S1', 'C1_equal'):
        return {K: thr_for_count(P[K][1][(P[K][0] >= a) & (P[K][0] < b)], TG[(TG.target_set == tset) & (TG.template == str(K))
                                                                            & (TG.period == 'ETS2011')].n_target.iloc[0])
                for K in P}
    # C0: common k in MAD units, pooled 2011-ETS count = Table S1 pooled (485)
    ntot = int(TG[(TG.target_set == 'pooled_S1') & (TG.period == 'ETS2011')].n_target.iloc[0])
    def cnt(k):
        return sum(((P[K][1] > P[K][2] + k * P[K][3]) & (P[K][0] >= a) & (P[K][0] < b)).sum() for K in P)
    lo, hi = 0.0, 40.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if cnt(mid) > ntot else (lo, mid)
    return {K: P[K][2] + hi * P[K][3] for K in P}


def s1():
    S = pd.read_csv(f'{H.TB}/hg2019_tableS1.csv')
    return {K: np.sort(S[S.template == K].i.values) for K in range(1, 6)}


def near_any(x, ref, w):
    """bool per x: any ref within +-w (ref sorted)."""
    x = np.asarray(x, np.int64)
    if ref.size == 0 or x.size == 0:
        return np.zeros(x.size, bool)
    return (np.searchsorted(ref, x + w, side='right') - np.searchsorted(ref, x - w, side='left')) > 0


def calib_metrics(P, TH, SS):
    days = pd.date_range('2011-07-01', '2011-09-06', inclusive='left')
    a0 = H.tidx('2011-07-01')
    daily = np.zeros(len(days)); dref = np.zeros(len(days))
    pre_a, pre_b = H.per_idx('pre2011')
    npre = 0; per_t = {}; nth = nth_m = nour = nour_m = 0
    for K, (pi, pv, med, mad) in P.items():
        s = (pv > TH[K]) & (pi >= a0) & (pi < CAL_END)
        daily += np.bincount((pi[s] - a0) // 86400, minlength=len(days))[:len(days)]
        th = SS[K][(SS[K] >= a0) & (SS[K] < CAL_END)]
        dref += np.bincount((th - a0) // 86400, minlength=len(days))[:len(days)]
        ours = np.sort(pi[s])
        nth += th.size; nth_m += near_any(th, ours, 30).sum(); nour += ours.size; nour_m += near_any(ours, th, 30).sum()
        n = int(((pv > TH[K]) & (pi >= pre_a) & (pi < pre_b)).sum()); npre += n; per_t[K] = n
    r = float(np.corrcoef(daily, dref)[0, 1])
    return npre, r, per_t, nth_m / max(nth, 1), nour_m / max(nour, 1)


def main():
    TG = pd.read_csv(f'{H.TB}/r2_targets_S1.csv'); TG['template'] = TG.template.astype(str)
    SS = s1()
    npre_ref = int(TG[(TG.target_set == 'pooled_S1') & (TG.period == 'pre2011')].n_target.iloc[0])
    extra = [a for a in sys.argv[1:] if not a.startswith('--')]
    if extra:
        SPEC['V7'] = v7spec()[0] + ('+mask' if extra[0].endswith('+mask') else '')
    variants = list(SPEC)
    rows, thr_rows = [], []
    for v in variants:
        P = load_variant(v, extra[0] if v == 'V7' else None)
        for tset in ['S1', 'C1_equal', 'C0_commonMAD']:
            TH = thresholds(P, tset, TG)
            npre, r, per_t, rec, prec = calib_metrics(P, TH, SS)
            rows.append({'variant': v, 'spec': SPEC[v], 'target_set': tset, 'recall30_cal': rec, 'precision30_cal': prec,
                         'n_pre_ours': npre, 'n_pre_S1': npre_ref, 'score_abs_ln_pre': abs(np.log(max(npre, 1) / npre_ref)),
                         'r_daily_S1_cal': r, **{f'pre_T{K}': per_t[K] for K in per_t}})
            for K in P:
                thr_rows.append({'variant': v, 'target_set': tset, 'template': K, 'thr': TH[K],
                                 'thr_in_MAD_above_median': (TH[K] - P[K][2]) / P[K][3] if P[K][3] > 0 else np.nan,
                                 'median': P[K][2], 'MAD': P[K][3]})
    C = pd.DataFrame(rows); TT = pd.DataFrame(thr_rows)
    C.to_csv(f'{H.TB}/r2_calibration.csv', index=False); TT.to_csv(f'{H.TB}/r2_thresholds.csv', index=False)
    pd.set_option('display.width', 250)
    print(C.round(3).to_string(index=False))
    print(TT[TT.target_set == 'S1'].pivot(index='variant', columns='template', values='thr_in_MAD_above_median').round(2))
    c = C[C.target_set == 'S1'].set_index('variant')
    best = c.recall30_cal.max()
    near = c[c.recall30_cal >= best - 0.03]
    m = near.score_abs_ln_pre.min()
    near2 = near[near.score_abs_ln_pre <= m + 1e-9]
    choice = near2.r_daily_S1_cal.idxmax()
    # V7 rule: groups whose best single beats V0 (on the primary criterion, recall at +-30 s)
    v0 = c.loc['V0', 'recall30_cal']
    better = {}
    for v, g in GROUP.items():
        if v in c.index and c.loc[v, 'recall30_cal'] > v0:
            if g not in better or c.loc[v, 'recall30_cal'] > c.loc[better[g], 'recall30_cal']:
                better[g] = v
    out = {'chosen_variant': choice, 'chosen_spec': SPEC[choice], 'target_set': 'S1',
           'recall30_cal': float(c.loc[choice, 'recall30_cal']), 'score_pre': float(c.loc[choice, 'score_abs_ln_pre']),
           'r_daily': float(c.loc[choice, 'r_daily_S1_cal']), 'groups_beating_V0': better,
           'V7_needed': len(better) >= 2 and 'V7' not in c.index,
           'V7_file': extra[0] if extra else None,
           'made': pd.Timestamp.now(tz='America/Los_Angeles').strftime('%Y-%m-%d %H:%M PT')}
    json.dump(out, open(f'{H.TB}/r2_choice.json', 'w'), indent=1)
    print(out)
    line = (f"- {out['made']}: R2a calibration choice (1 Jul-6 Sep 2011 only, Table-S1 targets, before any held-out "
            f"count): {choice} ({out['chosen_spec']}), recall of Table-S1 calibration detections within 30 s "
            f"{out['recall30_cal']:.3f}, pre-ETS |ln| {out['score_pre']:.3f}, daily r {out['r_daily']:.2f}; groups beating V0: "
            f"{better}; V7 needed: {out['V7_needed']}. Table tables/r2_calibration.csv.\n")
    open(f'{H.HG}/NOTES.md', 'a').write(line)


def v7spec():
    """combination of the best single of each group that beat V0 (R2.3 V7 rule) -> (spec, file pattern)."""
    ch = json.load(open(f'{H.TB}/r2_choice.json'))
    kind, off, norm, chset, hw, mask = 'disp', '0', 'sum', 'all', '90', ''
    for g, v in ch['groups_beating_V0'].items():
        s = SPEC[v]
        if v == 'V8':
            mask = '+mask'
            continue
        k_, o_, n_, c_, h_ = s.split(',')
        if g == 'norm':
            norm = n_
        elif g == 'chset':
            chset = c_
        else:
            kind, off, hw = (k_ if k_ != 'disp' else kind), (o_ if o_ != '0' else off), (h_ if h_ != '90' else hw)
    spec = f'{kind},{off},{norm},{chset},{hw}'
    import hg_r2
    return spec, f'{R2}/real_{hg_r2.tag(spec)}_TK.npz{mask}'


if __name__ == '__main__':
    if sys.argv[1:] == ['--v7spec']:
        print(' '.join(v7spec()))
    else:
        main()
