#!/usr/bin/env python
"""Amendment R2b (PLAN.md): debug the low overlap between our calibrated catalogue and Table S1.
Modes
  self              R2b.1 self-detection check (V0 displacement and V4 velocity templates) -> tables/r2b_self.csv
  stack K           R2b.2 stacks at Table-S1 times vs controls for template K -> data/r2b/stack_TK.parquet + stack_TK.npz
  stacksum          collect R2b.2 -> tables/r2b_stack.csv, r2b_pertrace.csv, figures/r2b_stack.png
  overlap           R2b.3 rank / tolerance / lag histogram / self-consistency -> tables/r2b_overlap_*.csv, figures
  grid              R2b.3 processing grid on the calibration period only -> tables/r2b_grid.csv
Conventions: template window [origin + off, origin + off + 90]; our detection index = origin time (window start - off);
Table-S1 time t ~ origin - 9..11 s, so their window start corresponds to t + delta, delta in {0, 9, 11}."""
import os, sys, json, glob
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.signal import oaconvolve, butter, sosfiltfilt, sosfilt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hg_common as H
import hg_r2 as R
import hg_mf

OUT = f'{H.HG}/data/r2b'
os.makedirs(OUT, exist_ok=True)
LAGS = 15
DELTAS = [0, 9, 11]
NDR = 5
PERIODS = {'ETS2011': H.PER['ETS2011'], 'interETS': H.PER['interETS'], 'ETS2012': H.PER['ETS2012'],
           'inter_SON2011': ('2011-09-06', '2011-12-01'), 'inter_DJF': ('2011-12-01', '2012-03-01'),
           'inter_MAM': ('2012-03-01', '2012-06-01'), 'inter_JJA': ('2012-06-01', '2012-08-30')}


def log(msg):
    line = f'{pd.Timestamp.now(tz="America/Los_Angeles"):%Y-%m-%d %H:%M:%S} PT r2b {msg}'
    print(line, flush=True)
    open(f'{H.LOGD}/r2b.log', 'a').write(line + '\n')


def s1():
    return pd.read_csv(f'{H.TB}/hg2019_tableS1.csv')


def tmpl(K, kind, off=0):
    io = H.tidx(H.templates().set_index('template').loc[K].time)
    T = R.cut(io + off, kind)
    chs = R.chans_for(K, T, 'all')
    return io, chs, np.array([T[c][0] for c in chs]), np.array([T[c][1] for c in chs])


def seg_cc(K, kind, a, b, chs, W):
    """summed CC (sum over channels) for window starts a..b-1 (global indices)."""
    S = np.zeros(b - a); NL = np.zeros(b - a, int)
    for ch, w in zip(chs, W):
        s, comp = ch.rsplit('.', 1)
        X, D = R.cont(s, kind)
        c = H.COMPS.index(comp)
        x = np.nan_to_num(np.asarray(X[c, a:b + H.L - 1], np.float64))
        num = np.correlate(x, w, mode='valid')
        d = np.asarray(D[c, a:b], np.float64)
        live = d > 0
        S += np.where(live, num / np.where(live, d, 1), 0.0); NL += live
    return S, NL


# ---------------------------------------------------------------- R2b.1
def run_self():
    S1 = s1(); TT = H.templates()
    rows = []
    for kind in ['disp', 'vel']:
        for K in range(1, 6):
            io, chs, W, _ = tmpl(K, kind)
            day0 = (io // 86400) * 86400
            a, b = day0 - 86400, day0 + 2 * 86400
            S, NL = seg_cc(K, kind, a, b, chs, W)
            mf = pd.Series(S).rolling(2 * H.L + 1, center=True, min_periods=1).max().values
            pk = np.where((S == mf) & (S > 0))[0]
            dsel = pk[(pk + a >= day0) & (pk + a < day0 + 86400)]
            vals = np.sort(S[dsel])[::-1]
            r = {'kind': kind, 'template': K, 'n_ch': len(chs), 'S_at_origin': S[io - a],
                 'self_is_day_max': bool(S[io - a] >= vals[0] - 1e-6), 'day_max_excl_self': float(
                     np.max(S[dsel[np.abs(dsel + a - io) > H.L]])) if dsel.size else np.nan}
            th = S1[(S1.template == K) & (np.abs(S1.i - io) <= 300)].i.values
            if th.size:
                t = int(th[np.argmin(np.abs(th - io))])
                r['S1_time_minus_origin_s'] = t - io
                r['S_at_S1_time'] = S[t - a]
                w = S[t - a - LAGS: t - a + LAGS + 1]
                r['S_best_within15_of_S1'] = float(w.max()); r['lag_best_vs_S1_s'] = int(np.argmax(w) - LAGS)
            for K2 in range(1, 6):
                if K2 == K:
                    continue
                io2 = H.tidx(TT.set_index('template').loc[K2].time)
                S2, _ = seg_cc(K, kind, io2 - 60, io2 + 61, chs, W)
                r[f'S_at_T{K2}_origin_best60'] = float(S2.max())
            rows.append(r)
            log(f'self {kind} T{K} done')
    D = pd.DataFrame(rows)
    D.to_csv(f'{H.TB}/r2b_self.csv', index=False)
    pd.set_option('display.width', 250)
    print(D.round(3).to_string(index=False))


# ---------------------------------------------------------------- R2b.2
class Data:
    """all channels of the configuration in RAM (float32), for fast window cuts."""
    def __init__(self, chs, kind):
        self.X = {}
        for ch in chs:
            s, comp = ch.rsplit('.', 1)
            self.X[ch] = np.asarray(R.cont(s, kind)[0][H.COMPS.index(comp)], np.float32)

    def cut(self, ch, starts, n):
        idx = np.asarray(starts, np.int64)[:, None] + np.arange(n)[None, :]
        idx = np.clip(idx, 0, H.N - 1)
        return self.X[ch][idx]


def metrics(data, chs, W, starts, shifts=None, signs=None):
    """per-trace mean CC (lag 0 and best within +-LAGS) and stack CC / semblance for windows starting at `starts`
    (window start = template-window start). shifts: per-station day shift (samples); signs: per-channel template sign."""
    n = len(starts)
    L = H.L + 2 * LAGS
    ccl = np.zeros((n, 2 * LAGS + 1)); nl = np.zeros(n)
    stacks = []
    for k, (ch, w) in enumerate(zip(chs, W)):
        si = H.STAS.index(ch.rsplit('.', 1)[0])
        st = np.asarray(starts, np.int64) - LAGS + (0 if shifts is None else shifts[si])
        st = np.mod(st, H.N - L)
        Xw = data.cut(ch, st, L).astype(np.float64)                     # (n, L)
        ok = np.all(np.isfinite(Xw), 1)
        Xw = np.where(np.isfinite(Xw), Xw, 0.0)
        ww = w * (1.0 if signs is None else signs[k])
        # sliding CC over lags
        win = np.lib.stride_tricks.sliding_window_view(Xw, H.L, axis=1)  # (n, 2*LAGS+1, 90)
        m = win.mean(-1, keepdims=True)
        den = np.sqrt(((win - m) ** 2).sum(-1))
        num = (win * ww).sum(-1)
        c = np.where(den > 0, num / np.where(den > 0, den, 1), 0.0)
        live = ok & (den[:, LAGS] > 0)
        ccl += np.where(live[:, None], c, 0.0); nl += live
        # stack at lag 0 (fixed alignment), unit-normalised demeaned traces
        x0 = win[:, LAGS, :] - m[:, LAGS, :]
        nx = np.linalg.norm(x0, axis=1)
        good = live & (nx > 0)
        if good.sum() >= 2:
            xn = x0[good] / nx[good, None]
            sk = xn.sum(0)
            sc = float(np.dot(sk, ww) / np.linalg.norm(sk)) if np.linalg.norm(sk) > 0 else np.nan
            semb = float(np.dot(sk, sk) / good.sum() ** 2)
            stacks.append((sc, semb, good.sum()))
    ccl = ccl / np.maximum(nl, 1)[:, None]
    sc = np.array([s[0] for s in stacks]); sb = np.array([s[1] for s in stacks]); ng = np.array([s[2] for s in stacks])
    return {'cc0': ccl[:, LAGS], 'ccbest': ccl.max(1), 'lagbest': ccl.argmax(1) - LAGS,
            'stack_cc': float(np.nanmean(sc)) if sc.size else np.nan, 'semblance': float(np.mean(sb)) if sb.size else np.nan,
            'semblance_chance': float(np.mean(1.0 / ng)) if ng.size else np.nan}


def run_stack(K):
    S1 = s1(); S1 = S1[S1.template == K]
    det = pd.read_csv(f'{H.TB}/r2_match_S1_detail.csv')
    det = det[det.template == K].set_index('i')
    ours = pd.read_csv(f'{H.TB}/r2_catalogue_chosen.csv'); ours = ours[ours.template == K]
    offs, signs = hg_mf.draws()
    rng = np.random.default_rng(H.SEED + 60 + K)
    rows, trace_rows = [], []
    for kind in ['disp', 'vel']:
        io, chs, W, Wn = tmpl(K, kind)
        data = Data(chs, kind)
        sgn = [np.array([signs[r, H.CHANS.index(c)] for c in chs]) for r in range(NDR)]
        # tercile of our CC at their times (all Table-S1 detections of this template)
        cc_ours = det.reindex(S1.i.values).cc_mad_ours.values
        q = np.nanpercentile(cc_ours, [100 / 3, 200 / 3])
        terc = np.where(cc_ours <= q[0], 'low', np.where(cc_ours <= q[1], 'mid', 'high'))
        for pname, (pa, pb) in PERIODS.items():
            a, b = H.tidx(pa), H.tidx(pb)
            msk = (S1.i.values >= a) & (S1.i.values < b)
            if msk.sum() < 3:
                continue
            groups = [('all', msk)] + ([(f'tercile_{t}', msk & (terc == t)) for t in ['low', 'mid', 'high']]
                                       if pname in ('ETS2011', 'interETS', 'ETS2012') else [])
            for gname, gm in groups:
                t = S1.i.values[gm]
                if t.size < 3:
                    continue
                for dl in DELTAS:
                    base = {'template': K, 'kind': kind, 'period': pname, 'group': gname, 'delta_s': dl, 'n': t.size}
                    m = metrics(data, chs, W, t + dl)
                    rows.append({**base, 'set': 'Table S1', 'draw': 0, 'stack_cc': m['stack_cc'], 'semblance': m['semblance'],
                                 'semblance_chance': m['semblance_chance'], 'cc0_median': np.median(m['cc0']),
                                 'ccbest_median': np.median(m['ccbest'])})
                    if gname == 'all':
                        trace_rows.append(pd.DataFrame({'template': K, 'kind': kind, 'period': pname, 'delta_s': dl,
                                                        'set': 'Table S1', 'i': t, 'cc0': m['cc0'], 'ccbest': m['ccbest'],
                                                        'lagbest': m['lagbest']}))
                    if dl != 9 and gname != 'all':
                        continue
                    for r in range(NDR):
                        day = (t // 86400) * 86400
                        tr_ = day + rng.integers(0, 86400 - 200, t.size)
                        ctl = {'random same day': metrics(data, chs, W, tr_),
                               'station-day shift': metrics(data, chs, W, t + dl, shifts=offs[r] * 86400),
                               'polarity-scrambled template': metrics(data, chs, W, t + dl, signs=sgn[r])}
                        for nm, mm in ctl.items():
                            rows.append({**base, 'set': nm, 'draw': r, 'stack_cc': mm['stack_cc'], 'semblance': mm['semblance'],
                                         'semblance_chance': mm['semblance_chance'], 'cc0_median': np.median(mm['cc0']),
                                         'ccbest_median': np.median(mm['ccbest'])})
                            if gname == 'all' and r == 0 and dl == 9:
                                trace_rows.append(pd.DataFrame({'template': K, 'kind': kind, 'period': pname, 'delta_s': dl,
                                                                'set': nm, 'i': -1, 'cc0': mm['cc0'], 'ccbest': mm['ccbest'],
                                                                'lagbest': mm['lagbest']}))
            # our own detections (selected by our template) as a reference, and selection-matched noise templates
            if pname in ('ETS2011', 'interETS', 'ETS2012'):
                o = ours[(ours.i >= a) & (ours.i < b)].i.values
                if o.size >= 3:
                    mm = metrics(data, chs, W, o)
                    rows.append({'template': K, 'kind': kind, 'period': pname, 'group': 'all', 'delta_s': 0, 'n': o.size,
                                 'set': 'ours (V4 calibrated)', 'draw': 0, 'stack_cc': mm['stack_cc'],
                                 'semblance': mm['semblance'], 'semblance_chance': mm['semblance_chance'],
                                 'cc0_median': np.median(mm['cc0']), 'ccbest_median': np.median(mm['ccbest'])})
                    trace_rows.append(pd.DataFrame({'template': K, 'kind': kind, 'period': pname, 'delta_s': 0,
                                                    'set': 'ours (V4 calibrated)', 'i': o, 'cc0': mm['cc0'],
                                                    'ccbest': mm['ccbest'], 'lagbest': mm['lagbest']}))
                for j, (ni0, npi) in enumerate(noise_top(K, kind, a, b, int(msk.sum()))):
                    if npi.size < 3:
                        continue
                    T = R.cut(ni0, kind)
                    Wn_ = np.array([T[c][0] if T[c] is not None else np.zeros(H.L) for c in chs])
                    mm = metrics(data, chs, Wn_, npi)
                    rows.append({'template': K, 'kind': kind, 'period': pname, 'group': 'all', 'delta_s': 0, 'n': npi.size,
                                 'set': 'selection-matched noise template (vs own window)', 'draw': j,
                                 'stack_cc': mm['stack_cc'], 'semblance': mm['semblance'],
                                 'semblance_chance': mm['semblance_chance'], 'cc0_median': np.median(mm['cc0']),
                                 'ccbest_median': np.median(mm['ccbest'])})
                    if j == 0:
                        trace_rows.append(pd.DataFrame({'template': K, 'kind': kind, 'period': pname, 'delta_s': 0,
                                                        'set': 'selection-matched noise template', 'i': npi,
                                                        'cc0': mm['cc0'], 'ccbest': mm['ccbest'], 'lagbest': mm['lagbest']}))
            log(f'stack T{K} {kind} {pname} done')
        del data
    pd.DataFrame(rows).to_parquet(f'{OUT}/stack_T{K}.parquet')
    pd.concat(trace_rows, ignore_index=True).to_parquet(f'{OUT}/pertrace_T{K}.parquet')
    return f'T{K} stack done'


def noise_top(K, kind, a, b, n, nj=10):
    """(noise window start, its top-n detection starts in [a, b)) for the first nj noise templates, configuration of T{K}.
    vel: chosen-setting runs data/r2/noise_vel_0_sum_all_90_*.npz; disp: R1 runs data/mf/noise_*.npz (T{K}_prim)."""
    out = []
    for j in range(nj):
        if kind == 'vel':
            f = f'{R.R2}/noise_{R.tag("vel,0,sum,all,90")}_{j:03d}.npz'
            if not os.path.isfile(f):
                continue
            z = np.load(f, allow_pickle=True)
            pi = pv = None
            for q in range(int(z['ncfg'])):
                if K in list(z[f'c{q}_templates']):
                    pi, pv = z[f'c{q}_pi'], z[f'c{q}_pv']
        else:
            import hg_analyze as A
            z = np.load(f'{H.MF}/noise_{j:03d}.npz', allow_pickle=True)
            pi, pv, _ = A.noise_peaks(z, f'T{K}_prim')
        i0 = int(z['i0'])
        s = (pi >= a) & (pi < b) & (np.abs(pi - i0) > H.L)
        pi, pv = pi[s], pv[s]
        out.append((i0, np.sort(pi[np.argsort(pv)[::-1][:n]])))
    return out


def run_stacksum():
    D = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f'{OUT}/stack_T?.parquet'))], ignore_index=True)
    P = pd.concat([pd.read_parquet(f) for f in sorted(glob.glob(f'{OUT}/pertrace_T?.parquet'))], ignore_index=True)
    G = D.groupby(['template', 'kind', 'period', 'group', 'delta_s', 'set']).agg(
        n=('n', 'first'), stack_cc=('stack_cc', 'mean'), stack_cc_max=('stack_cc', 'max'), semblance=('semblance', 'mean'),
        semblance_chance=('semblance_chance', 'mean'), cc0_median=('cc0_median', 'mean'),
        ccbest_median=('ccbest_median', 'mean')).reset_index()
    G.to_csv(f'{H.TB}/r2b_stack.csv', index=False)
    P.groupby(['template', 'kind', 'period', 'delta_s', 'set']).agg(n=('cc0', 'size'), cc0_median=('cc0', 'median'),
                                                                     ccbest_median=('ccbest', 'median'),
                                                                     ccbest_p90=('ccbest', lambda x: np.percentile(x, 90)),
                                                                     lag_abs_median=('lagbest', lambda x: np.median(np.abs(x)))
                                                                     ).reset_index().to_csv(f'{H.TB}/r2b_pertrace.csv', index=False)
    pd.set_option('display.width', 250)
    g = G[(G.group == 'all') & (G.period.isin(['ETS2011', 'interETS', 'ETS2012'])) & (G.kind == 'vel') &
          ((G.delta_s == 9) | (G.set.str.startswith('ours')) | (G.set.str.startswith('selection')))]
    print(g.pivot_table(index=['template', 'period'], columns='set', values='stack_cc').round(3).to_string())
    # figure: per-trace best CC distributions (velocity, delta 9), per template, ETS2011 and interETS
    fig, ax = plt.subplots(2, 5, figsize=(18, 6.5), sharex=True)
    bins = np.linspace(-0.3, 1.0, 66)
    cols = {'Table S1': 'k', 'random same day': '0.6', 'station-day shift': 'g', 'polarity-scrambled template': 'purple',
            'ours (V4 calibrated)': 'r', 'selection-matched noise template': 'c'}
    for row, pname in enumerate(['ETS2011', 'interETS']):
        for K in range(1, 6):
            g_ = ax[row, K - 1]
            for nm, c in cols.items():
                x = P[(P.template == K) & (P.kind == 'vel') & (P.period == pname) & (P.set == nm) &
                      ((P.delta_s == 9) | P.set.isin(['ours (V4 calibrated)', 'selection-matched noise template']))]
                if len(x):
                    g_.hist(x.ccbest, bins, histtype='step', color=c, density=True, label=f'{nm} ({len(x)})')
            g_.set_title(f'T{K} {pname}', fontsize=9)
            if row == 1:
                g_.set_xlabel('mean per-channel CC with our template (best lag ±15 s)')
        ax[row, 0].legend(fontsize=6)
    fig.suptitle('Do Table-S1 windows look like our templates? velocity, their time + 9 s (selection-matched = top-N '
                 'detections of a noise template vs its own window)', fontsize=10)
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2b_pertrace.png', dpi=120); plt.close(fig)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4), sharey=True)
    for g_, pname in zip(ax, ['ETS2011', 'interETS', 'ETS2012']):
        x = G[(G.kind == 'vel') & (G.period == pname) & (G.group == 'all')]
        for k, (nm, c) in enumerate(cols.items()):
            y = x[(x.set == nm) & ((x.delta_s == 9) | x.set.isin(['ours (V4 calibrated)']))]
            if nm.startswith('selection'):
                y = x[x.set.str.startswith('selection')]
            g_.plot(y.template + (k - 2.5) * 0.1, y.stack_cc, 'o', color=c, label=nm)
        g_.set_title(pname); g_.set_xlabel('template'); g_.set_xticks(range(1, 6))
    ax[0].set_ylabel('stack CC with template (velocity)'); ax[0].legend(fontsize=6)
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2b_stack.png', dpi=120); plt.close(fig)


# ---------------------------------------------------------------- R2b.3 overlap
def within(x, ref, lo, hi):
    """bool per x: any ref in [x + lo, x + hi] (ref sorted)."""
    if ref.size == 0 or len(x) == 0:
        return np.zeros(len(x), bool)
    x = np.asarray(x, np.int64)
    return (np.searchsorted(ref, x + hi, side='right') - np.searchsorted(ref, x + lo, side='left')) > 0


def run_overlap():
    import hg_r2_calib as RC
    S1 = s1()
    C = pd.read_csv(f'{H.TB}/r2_catalogue_chosen.csv')
    det = pd.read_csv(f'{H.TB}/r2_match_S1_detail.csv')
    OFF = 11
    rows = []
    for pname in ['ETS2011', 'interETS', 'ETS2012']:
        a, b = H.per_idx(pname)
        for K in list(range(1, 6)) + ['all']:
            ks = range(1, 6) if K == 'all' else [K]
            for frac in [0.10, 0.25, 0.50, 1.0]:
                m30 = m5 = n = 0; r30 = r5 = nr = 0
                for k in ks:
                    th = np.sort(S1[(S1.template == k)].i.values)
                    o = C[(C.template == k) & (C.i >= a) & (C.i < b)].sort_values('S', ascending=False)
                    o = o.iloc[:max(1, int(round(frac * len(o))))] if len(o) else o
                    # ours - theirs in [-30, 30]  <=> theirs in [ours-30, ours+30]
                    m30 += within(o.i.values, th, -30, 30).sum(); m5 += within(o.i.values, th, -OFF - 5, -OFF + 5).sum()
                    n += len(o)
                    d = det[(det.template == k) & (det.i >= a) & (det.i < b)].sort_values('cc_mad_ours', ascending=False)
                    d = d.iloc[:max(1, int(round(frac * len(d))))] if len(d) else d
                    oi = np.sort(C[C.template == k].i.values)
                    r30 += within(d.i.values, oi, -30, 30).sum(); r5 += within(d.i.values, oi, OFF - 5, OFF + 5).sum()
                    nr += len(d)
                rows.append({'period': pname, 'template': K, 'top_frac': frac, 'n_ours': n,
                             'ours_in_S1_30': m30 / max(n, 1), 'ours_in_S1_5off': m5 / max(n, 1), 'n_S1': nr,
                             'S1_in_ours_30': r30 / max(nr, 1), 'S1_in_ours_5off': r5 / max(nr, 1)})
    RK = pd.DataFrame(rows); RK.to_csv(f'{H.TB}/r2b_overlap_rank.csv', index=False)
    pd.set_option('display.width', 250)
    print(RK[RK.template == 'all'].round(3).to_string(index=False))
    # tolerance curve + chance
    tl = []
    for pname in ['ETS2011', 'interETS']:
        a, b = H.per_idx(pname); nd = (b - a) / 86400
        for centre in [0, OFF]:
            for w in [2, 5, 10, 20, 30, 60]:
                hit = n = exp = 0.0
                for k in range(1, 6):
                    th = S1[(S1.template == k) & (S1.i >= a) & (S1.i < b)].i.values
                    oi = np.sort(C[C.template == k].i.values)
                    hit += within(th, oi, centre - w, centre + w).sum(); n += th.size
                    rate = ((oi >= a) & (oi < b)).sum() / (nd * 86400)
                    exp += th.size * (1 - np.exp(-rate * (2 * w + 1)))
                tl.append({'period': pname, 'centre_s': centre, 'halfwidth_s': w, 'recall': hit / n, 'chance': exp / n})
    TL = pd.DataFrame(tl); TL.to_csv(f'{H.TB}/r2b_overlap_tolerance.csv', index=False); print(TL.round(4).to_string(index=False))
    # lag histogram ours - theirs within +-600 s
    lags = {p: [] for p in ['ETS2011', 'interETS', 'ETS2012']}
    for p in lags:
        a, b = H.per_idx(p)
        for k in range(1, 6):
            th = S1[(S1.template == k) & (S1.i >= a) & (S1.i < b)].i.values
            oi = np.sort(C[C.template == k].i.values)
            for t in th:
                lo, hi = np.searchsorted(oi, t - 600), np.searchsorted(oi, t + 600, side='right')
                lags[p] += list(oi[lo:hi] - t)
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.6))
    for g_, (p, v) in zip(ax, lags.items()):
        g_.hist(v, bins=np.arange(-600, 601, 10), color='k'); g_.set_title(f'{p}: ours − Table S1 (s), same template')
        g_.set_xlabel('lag (s)')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2b_lag_hist.png', dpi=120); plt.close(fig)
    pd.DataFrame([{'period': p, 'lag_s': x} for p, v in lags.items() for x in v]).to_csv(
        f'{H.TB}/r2b_lag_pairs.csv', index=False)
    # self-consistency: our calibrated variants (Table-S1 thresholds) vs each other and vs Table S1, 2011 ETS
    TH = pd.read_csv(f'{H.TB}/r2_thresholds.csv'); TH = TH[TH.target_set == 'S1']
    a, b = H.per_idx('ETS2011')
    cats = {}
    for v in ['V0', 'V4', 'V5b', 'V6b', 'V3b', 'V1']:
        P = RC.load_variant(v)
        th = TH[TH.variant == v].set_index('template').thr
        cats[v] = {K: np.sort(P[K][0][(P[K][1] > th[K]) & (P[K][0] >= a) & (P[K][0] < b)]) for K in range(1, 6)}
    cats['Table S1 (+11 s)'] = {K: np.sort(S1[(S1.template == K) & (S1.i >= a) & (S1.i < b)].i.values + OFF)
                                for K in range(1, 6)}
    names = list(cats)
    M = np.zeros((len(names), len(names)))
    for x, A_ in enumerate(names):
        for y, B_ in enumerate(names):
            hit = sum(within(cats[A_][K], cats[B_][K], -30, 30).sum() for K in range(1, 6))
            M[x, y] = hit / max(sum(len(cats[A_][K]) for K in range(1, 6)), 1)
    SC = pd.DataFrame(M, index=names, columns=names)
    SC.to_csv(f'{H.TB}/r2b_selfconsistency.csv'); print(SC.round(2).to_string())
    per = []
    for K in range(1, 6):
        for A_, B_ in [('V0', 'V4'), ('V0', 'V6b'), ('V4', 'V5b'), ('V0', 'V3b'), ('V4', 'Table S1 (+11 s)'), ('V0', 'Table S1 (+11 s)')]:
            per.append({'template': K, 'A': A_, 'B': B_, 'frac_A_in_B_30': within(cats[A_][K], cats[B_][K], -30, 30).mean()
                        if len(cats[A_][K]) else np.nan})
    pd.DataFrame(per).to_csv(f'{H.TB}/r2b_selfconsistency_per_template.csv', index=False)
    print(pd.DataFrame(per).pivot_table(index=['A', 'B'], columns='template', values='frac_A_in_B_30').round(2))


# ---------------------------------------------------------------- R2b.3 grid (calibration period only)
GA, GB = H.tidx('2011-06-24'), H.tidx('2011-09-09')          # data span (margins for templates/taper)
CA, CB = H.tidx('2011-07-01'), H.tidx('2011-09-06')
GRID = {'G0 V4 as is': {}, 'G1 window origin-11': {'off': -11}, 'G2 Z only': {'chset': 'Z'},
        'G3 horizontals only': {'chset': 'H'}, 'G4 5 nearest stations': {'chset': 'near5'}, 'G5 network norm': {'norm': 'net'},
        'G6 causal 4-pole': {'filt': 'causal4'}, 'G7 zero-phase 2-pole': {'filt': 'zp2'}, 'G8 corners 0.01-0.05': {'filt': 'wide'},
        'G9 decluster 20 s': {'hw': 20}, 'G10 displacement': {'kind': 'disp'}}


def filt_variant(x, how):
    if how == 'zp4':
        sos, f = butter(2, H.BAND, btype='band', fs=1.0, output='sos'), sosfiltfilt
    elif how == 'causal4':
        sos, f = butter(2, H.BAND, btype='band', fs=1.0, output='sos'), sosfilt
    elif how == 'zp2':
        sos, f = butter(1, H.BAND, btype='band', fs=1.0, output='sos'), sosfiltfilt
    elif how == 'wide':
        sos, f = butter(2, (0.01, 0.05), btype='band', fs=1.0, output='sos'), sosfiltfilt
    import hg_cont
    y = np.full(x.shape, np.nan, np.float64)
    a_, b_ = hg_cont.runs(np.isfinite(x))
    w = np.hanning(240)
    for i, j in zip(a_, b_):
        if j - i < 1800:
            continue
        s = x[i:j].astype(np.float64); t = np.arange(s.size)
        s = s - np.polyval(np.polyfit(t, s, 1), t)
        s[:120] *= w[:120]; s[-120:] *= w[120:]
        g = f(sos, s); g[:120] = np.nan; g[-120:] = np.nan
        y[i:j] = g
    return y


_G = {}


def gdata(sta, filt, kind):
    k = (sta, filt, kind)
    if k not in _G:
        nd = (GB - GA) // 86400
        X = np.full((3, GB - GA), np.nan)
        for q in range(nd):
            p = H.day_file(sta, H.T0 + pd.Timedelta(seconds=GA + q * 86400))
            if p is not None:
                X[:, q * 86400:(q + 1) * 86400] = np.load(p)['d']
        Y = np.array([filt_variant(X[c], filt) for c in range(3)])
        if kind == 'vel':
            V = np.full(Y.shape, np.nan); V[:, 1:-1] = 0.5 * (Y[:, 2:] - Y[:, :-2]); Y = V
        _G[k] = Y
    return _G[k]


def gcc(x, w):
    """normalised CC of unit template w along x (NaN -> dead)."""
    bad = ~np.isfinite(x)
    x0 = np.where(bad, 0.0, x)
    num = oaconvolve(x0, w[::-1], mode='valid')
    s1 = oaconvolve(x0, np.ones(H.L), mode='valid'); s2 = oaconvolve(x0 * x0, np.ones(H.L), mode='valid')
    nb = oaconvolve(bad.astype(float), np.ones(H.L), mode='valid')
    den = np.sqrt(np.maximum(s2 - s1 * s1 / H.L, 0))
    live = (nb < 0.5) & (den > 1e-12 * np.nanmax(den))
    return np.where(live, num / np.where(live, den, 1), 0.0), live, np.where(live, den, 0.0)


def run_grid():
    S1 = s1(); TT = H.templates().set_index('template')
    import hg_r2_calib as RC
    out = []
    base = {'kind': 'vel', 'off': 0, 'chset': 'all', 'norm': 'sum', 'filt': 'zp4', 'hw': 90}
    for gname, ch in GRID.items():
        P = {**base, **ch}
        res = {'grid': gname}
        n_th = n_hit = n_hit5 = n_our = n_top = n_top_hit = 0
        lags = []
        for K in range(1, 6):
            io = H.tidx(TT.loc[K].time)
            i0 = io + P['off'] - GA
            chs, Ws, Wn = [], [], []
            for s in H.STAS:
                Y = gdata(s, P['filt'], P['kind'])
                for c in range(3):
                    w = Y[c, i0:i0 + H.L]
                    if not np.all(np.isfinite(w)):
                        continue
                    chs.append((s, c)); w = w - w.mean(); Wn.append(np.linalg.norm(w)); Ws.append(w / np.linalg.norm(w))
            if P['chset'] == 'Z':
                keep = [q for q, (s, c) in enumerate(chs) if c == 2]
            elif P['chset'] == 'H':
                keep = [q for q, (s, c) in enumerate(chs) if c != 2]
            elif P['chset'] == 'near5':
                r = TT.loc[K]
                d = {s: np.hypot(R.COORD[s][0] - r.lat, (R.COORD[s][1] - r.lon) * np.cos(np.radians(r.lat))) for s in H.STAS}
                near = sorted(d, key=d.get)[:5]
                keep = [q for q, (s, c) in enumerate(chs) if s in near]
            else:
                keep = list(range(len(chs)))
            S = None; D2 = None; W2 = None
            for q in keep:
                s, c = chs[q]
                cc, live, den = gcc(gdata(s, P['filt'], P['kind'])[c], Ws[q])
                if S is None:
                    S = np.zeros(cc.size); D2 = np.zeros(cc.size); W2 = np.zeros(cc.size)
                if P['norm'] == 'net':
                    S += cc * den * Wn[q]; D2 += den ** 2; W2 += live * Wn[q] ** 2
                else:
                    S += cc
            if P['norm'] == 'net':
                dd = np.sqrt(D2 * W2); S = np.where(dd > 0, S / np.where(dd > 0, dd, 1), 0.0)
            NL = np.ones(S.size, int)
            pi, pv, _, med, mad = R.peaks(S.astype(np.float32), NL, P['hw'], P['off'])
            pi = pi + GA
            keepp = (np.abs(pi - io) > H.L) & (pi >= CA) & (pi < CB)
            pi, pv = pi[keepp], pv[keepp]
            ea, eb = H.per_idx('ETS2011')
            nt = int(((S1.template == K) & (S1.i >= ea) & (S1.i < eb)).sum())
            thr = RC.thr_for_count(pv[(pi >= ea) & (pi < eb)], nt)
            sel = pv > thr
            ours = np.sort(pi[sel]); ov = pv[sel][np.argsort(pi[sel])]
            th = np.sort(S1[(S1.template == K) & (S1.i >= CA) & (S1.i < CB)].i.values)
            hit = within(th, ours, -30, 30)
            n_th += th.size; n_hit += hit.sum(); n_our += ours.size
            # offset-aware +-5 s
            for t in th[hit]:
                j = np.argmin(np.abs(ours - t)); lags.append(ours[j] - t)
            res[f'T{K}_recall30'] = hit.mean() if th.size else np.nan
            res[f'T{K}_thr_MAD'] = (thr - med) / mad
            top = ours[ov >= np.percentile(ov, 75)] if ours.size else ours
            n_top += top.size; n_top_hit += within(top, th, -30, 30).sum()
            res[f'_ours_T{K}'] = ours; res[f'_th_T{K}'] = th
        off = float(np.median(lags)) if lags else 0.0
        h5 = sum(within(res[f'_th_T{K}'], res[f'_ours_T{K}'], off - 5, off + 5).sum() for K in range(1, 6))
        res.update({'recall30_cal': n_hit / max(n_th, 1), 'median_offset_s': off, 'recall5_off_cal': h5 / max(n_th, 1),
                    'ours_top25_in_S1_30': n_top_hit / max(n_top, 1), 'n_theirs_cal': n_th, 'n_ours_cal': n_our})
        out.append({k: v for k, v in res.items() if not k.startswith('_')})
        log(f'grid {gname}: recall30 {res["recall30_cal"]:.3f}')
    G = pd.DataFrame(out); G.to_csv(f'{H.TB}/r2b_grid.csv', index=False)
    pd.set_option('display.width', 250)
    print(G.round(3).to_string(index=False))


if __name__ == '__main__':
    m = sys.argv[1]
    if m == 'self':
        run_self()
    elif m == 'stack':
        log(run_stack(int(sys.argv[2])))
    elif m == 'stacksum':
        run_stacksum()
    elif m == 'overlap':
        run_overlap()
    elif m == 'grid':
        run_grid()
