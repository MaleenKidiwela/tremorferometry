#!/usr/bin/env python
"""Amendment R2c: moment budget of Table-S1 detections per period (2011 ETS, inter-ETS, 2012 ETS) by least-squares
amplitude scaling of each template (displacement, 0.02-0.05 Hz, raw window from the Table-1 origin) onto the data at
t + delta_K (template self-offset), cross-checked with our whole-space thrust synthetic; noise contribution from random
same-day times (20 draws) and from selection-matched noise-template detections (10 noise templates, top-N).
Outputs tables/r2c_moment_template.csv, tables/r2c_moment_pooled.csv, tables/r2c_threeway.csv, tables/r2c_selfcal.csv,
figures/r2c_moment.png."""
import os, sys
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/scripts')
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hg_common as H
import hg_r2 as R
import hg_r2b as B
import s13_green as G

DELTA = {1: 9, 2: 12, 3: 9, 4: 9, 5: 9}
PERS = ['ETS2011', 'interETS', 'ETS2012']
NRND, NJ = 20, 10
SOS = butter(2, H.BAND, btype='band', fs=2.0, output='sos')


def synth_template(K, chs):
    r = H.templates().set_index('template').loc[K]
    t = np.arange(-200, 400, 0.5)
    M = G.mt_ned(315, 12, 90, 1.0)
    out = []
    for ch in chs:
        s, comp = ch.rsplit('.', 1)
        n, e = G.ned_offset(r.lat, r.lon, *R.COORD[s])
        u = sosfiltfilt(SOS, G.synth(M, (0, 0, r.depth_km * 1e3), (n, e, 0.0), t, tau=10.0), axis=-1)
        g = u[H.COMPS.index(comp), 400:400 + 2 * H.L:2]
        out.append(g - g.mean())
    return np.array(out)


def scale(data, chs, W, starts, lags=0):
    """LS scale of W (nch, 90) onto windows at starts (fixed lag, or best |CC| lag within +-lags -> returns scale at the
    lag of maximum positive projection)."""
    st = np.asarray(starts, np.int64)
    best = np.full(st.size, -np.inf)
    for L in range(-lags, lags + 1):
        num = np.zeros(st.size); den = np.zeros(st.size)
        for k, ch in enumerate(chs):
            X = data.cut(ch, np.mod(st + L, H.N - H.L), H.L).astype(np.float64)
            ok = np.all(np.isfinite(X), 1)
            X = np.where(np.isfinite(X), X, 0.0)
            X = X - X.mean(1, keepdims=True)
            num += np.where(ok, X @ W[k], 0.0); den += np.where(ok, W[k] @ W[k], 0.0)
        a = np.where(den > 0, num / np.where(den > 0, den, 1), np.nan)
        best = np.where(np.isnan(best) | (a > best), a, best) if lags else a
    return best


def tsum(a, M0):
    """post-hoc robust total: N x 5-95 % trimmed mean of the scales x M0 (the plain sum is dominated by a few
    very large windows - earthquakes / glitches)."""
    a = a[np.isfinite(a)]
    if a.size == 0:
        return np.nan
    lo, hi = np.percentile(a, [5, 95])
    return a.size * a[(a >= lo) & (a <= hi)].mean() * M0


def main():
    S1 = pd.read_csv(f'{H.TB}/hg2019_tableS1.csv')
    TT = H.templates().set_index('template')
    rng = np.random.default_rng(H.SEED + 70)
    rows, cal = [], []
    noise_tot = {p: np.zeros(NJ) for p in PERS}; rnd_tot = {p: np.zeros(NRND) for p in PERS}
    boot_tot = {p: np.zeros(500) for p in PERS}
    rndt_tot = {p: np.zeros(NRND) for p in PERS}; noiset_tot = {p: np.zeros(NJ) for p in PERS}
    boott_tot = {p: np.zeros(500) for p in PERS}
    for K in range(1, 6):
        io = H.tidx(TT.loc[K].time); M0K = float(TT.loc[K].M0_Nm)
        T = R.cut(io, 'disp')
        chs = R.chans_for(K, T, 'all')
        data = B.Data(chs, 'disp')
        W = np.array([data.X[c][io:io + H.L] for c in chs], np.float64); W = W - W.mean(1, keepdims=True)
        Gs = synth_template(K, chs)
        a_self = scale(data, chs, W, [io])[0]
        msyn = scale(data, chs, Gs, [io])[0]
        cal.append({'template': K, 'M0_table1': M0K, 'M0_synthetic_fit': msyn, 'synthetic_over_table1': msyn / M0K,
                    'template_self_scale': a_self})
        wn = np.sqrt((W ** 2).sum())
        for p in PERS:
            a, b = H.per_idx(p)
            t = S1[(S1.template == K) & (S1.i >= a) & (S1.i < b)].i.values
            n = t.size
            r = {'template': K, 'period': p, 'n_S1': n, 'delta_s': DELTA[K]}
            if n == 0:
                rows.append(r); continue
            av = scale(data, chs, W, t + DELTA[K])
            ab = scale(data, chs, W, t + DELTA[K], lags=15)
            ms = scale(data, chs, Gs, t + DELTA[K])
            r.update({'sum_M0_template': np.nansum(av) * M0K, 'median_scale': np.nanmedian(av),
                      'trim_M0_template': tsum(av, M0K), 'trim_M0_synthetic': tsum(ms, 1.0),
                      'frac_positive': float(np.mean(av > 0)),
                      'sum_M0_template_bestlag15': np.nansum(ab) * M0K, 'sum_M0_synthetic': np.nansum(ms)})
            bi = [rng.integers(0, n, n) for _ in range(500)]
            bs = np.array([np.nansum(av[i]) * M0K for i in bi]); boot_tot[p] += bs
            bt = np.array([tsum(av[i], M0K) for i in bi]); boott_tot[p] += np.nan_to_num(bt)
            r['boot_trim_p16'], r['boot_trim_p84'] = np.nanpercentile(bt, [16, 84])
            r['boot_p16'], r['boot_p84'] = np.percentile(bs, [16, 84])
            # (i) random same-day times
            rs, rt = [], []
            for q in range(NRND):
                tr = (t // 86400) * 86400 + rng.integers(0, 86400 - 200, n)
                ar = scale(data, chs, W, tr)
                rs.append(np.nansum(ar) * M0K); rt.append(tsum(ar, M0K))
            rs = np.array(rs); rnd_tot[p] += rs; rt = np.array(rt); rndt_tot[p] += rt
            r.update({'random_trim_mean': rt.mean(), 'random_trim_p16': np.percentile(rt, 16),
                      'random_trim_p84': np.percentile(rt, 84)})
            r.update({'random_mean': rs.mean(), 'random_p16': np.percentile(rs, 16), 'random_p84': np.percentile(rs, 84)})
            # (ii) selection-matched noise templates (R1 displacement runs, T{K} configuration)
            ns, nt_ = [], []
            for j, (ni0, npi) in enumerate(B.noise_top(K, 'disp', a, b, n, nj=NJ)):
                Wn = np.array([data.X[c][ni0:ni0 + H.L] for c in chs], np.float64)
                if not np.all(np.isfinite(Wn)):
                    ns.append(np.nan); nt_.append(np.nan); continue
                Wn = Wn - Wn.mean(1, keepdims=True)
                an = scale(data, chs, Wn, npi)
                fac = np.sqrt((Wn ** 2).sum()) / wn * M0K * (n / max(npi.size, 1))
                ns.append(np.nansum(an) * fac); nt_.append(tsum(an, 1.0) * fac)
            ns = np.array(ns + [np.nan] * (NJ - len(ns))); nt_ = np.array(nt_ + [np.nan] * (NJ - len(nt_)))
            noise_tot[p] += np.nan_to_num(ns); noiset_tot[p] += np.nan_to_num(nt_)
            r.update({'selnoise_trim_mean': np.nanmean(nt_), 'selnoise_trim_p16': np.nanpercentile(nt_, 16),
                      'selnoise_trim_p84': np.nanpercentile(nt_, 84), 'net_trim_M0': r['trim_M0_template'] - np.nanmean(nt_)})
            r.update({'selnoise_mean': np.nanmean(ns), 'selnoise_p16': np.nanpercentile(ns, 16),
                      'selnoise_p84': np.nanpercentile(ns, 84),
                      'net_M0': r['sum_M0_template'] - np.nanmean(ns)})
            rows.append(r)
            print(f'T{K} {p} n={n} sum {r["sum_M0_template"]:.3g} sel-noise {r["selnoise_mean"]:.3g} '
                  f'random {r["random_mean"]:.3g} syn {r["sum_M0_synthetic"]:.3g}', flush=True)
        del data
    D = pd.DataFrame(rows); C = pd.DataFrame(cal)
    D.to_csv(f'{H.TB}/r2c_moment_template.csv', index=False); C.to_csv(f'{H.TB}/r2c_selfcal.csv', index=False)
    pool = []
    for p in PERS:
        d = D[D.period == p]
        pool.append({'period': p, 'n_S1': d.n_S1.sum(), 'sum_M0_template': d.sum_M0_template.sum(),
                     'boot_p16': np.percentile(boot_tot[p], 16), 'boot_p84': np.percentile(boot_tot[p], 84),
                     'sum_M0_template_bestlag15': d.sum_M0_template_bestlag15.sum(), 'sum_M0_synthetic': d.sum_M0_synthetic.sum(),
                     'random_mean': rnd_tot[p].mean(), 'random_p16': np.percentile(rnd_tot[p], 16),
                     'random_p84': np.percentile(rnd_tot[p], 84),
                     'selnoise_mean': noise_tot[p].mean(), 'selnoise_p16': np.percentile(noise_tot[p], 16),
                     'selnoise_p84': np.percentile(noise_tot[p], 84),
                     'net_M0': d.sum_M0_template.sum() - noise_tot[p].mean(),
                     'net_p16': d.sum_M0_template.sum() - np.percentile(noise_tot[p], 84),
                     'net_p84': d.sum_M0_template.sum() - np.percentile(noise_tot[p], 16),
                     'trim_M0_template': d.trim_M0_template.sum(), 'trim_boot_p16': np.percentile(boott_tot[p], 16),
                     'trim_boot_p84': np.percentile(boott_tot[p], 84), 'trim_M0_synthetic': d.trim_M0_synthetic.sum(),
                     'random_trim_mean': rndt_tot[p].mean(), 'random_trim_p16': np.percentile(rndt_tot[p], 16),
                     'random_trim_p84': np.percentile(rndt_tot[p], 84),
                     'selnoise_trim_mean': noiset_tot[p].mean(), 'selnoise_trim_p16': np.percentile(noiset_tot[p], 16),
                     'selnoise_trim_p84': np.percentile(noiset_tot[p], 84),
                     'net_trim_M0': d.trim_M0_template.sum() - noiset_tot[p].mean(),
                     'net_trim_p16': d.trim_M0_template.sum() - np.percentile(noiset_tot[p], 84),
                     'net_trim_p84': d.trim_M0_template.sum() - np.percentile(noiset_tot[p], 16)})
    P = pd.DataFrame(pool); P.to_csv(f'{H.TB}/r2c_moment_pooled.csv', index=False)
    # three-way split (SVI LFE numbers from S13; GNSS = literature range, not ours)
    sc = pd.read_csv(f'{H.S13}/tables/scale_svi_episodes.csv').set_index('episode')
    tw = []
    for p, ep in [('ETS2011', '2011-08'), ('ETS2012', '2012-08')]:
        n_lfe = int(sc.loc[ep, 'n']); m_lock = float(sc.loc[ep, 'thrust_m0'])
        pr = P.set_index('period').loc[p]
        tw.append({'ETS': ep, 'n_LFE_SVI': n_lfe, 'LFE_moment': n_lfe * 4.5e11,
                   'LFE_locked_slow_moment': n_lfe * m_lock, 'LFE_locked_sd': n_lfe * float(sc.loc[ep, 'thrust_ctrl_sd']),
                   'discrete_VLFE_S1_sum': pr.sum_M0_template, 'discrete_VLFE_net': pr.net_M0,
                   'discrete_VLFE_net_p16': pr.net_p16, 'discrete_VLFE_net_p84': pr.net_p84,
                   'discrete_VLFE_synthetic_sum': pr.sum_M0_synthetic,
                   'discrete_VLFE_trim': pr.trim_M0_template, 'discrete_VLFE_trim_selnoise': pr.selnoise_trim_mean,
                   'discrete_VLFE_net_trim': pr.net_trim_M0, 'discrete_VLFE_net_trim_p16': pr.net_trim_p16,
                   'discrete_VLFE_net_trim_p84': pr.net_trim_p84, 'discrete_VLFE_trim_synthetic': pr.trim_M0_synthetic,
                   'GNSS_lit_M0_low_Mw6.5': 10 ** (1.5 * 6.5 + 9.1), 'GNSS_lit_M0_high_Mw6.8': 10 ** (1.5 * 6.8 + 9.1)})
    TW = pd.DataFrame(tw); TW.to_csv(f'{H.TB}/r2c_threeway.csv', index=False)
    pd.set_option('display.width', 250)
    print(C.to_string(index=False)); print(D.round(3).to_string(index=False)); print(P.to_string(index=False))
    print(TW.T.to_string())
    fig, ax = plt.subplots(1, 1, figsize=(9, 4))
    x = np.arange(len(PERS))
    ax.bar(x - 0.3, P.sum_M0_template, 0.2, color='k', label='Table S1, template scaling (sum)')
    ax.bar(x - 0.1, P.selnoise_mean, 0.2, color='c', label='selection-matched noise (mean of 10)')
    ax.bar(x + 0.1, P.random_mean, 0.2, color='0.6', label='random same-day times (mean of 20)')
    ax.bar(x + 0.3, P.sum_M0_synthetic, 0.2, color='orange', label='Table S1, thrust-synthetic scaling (sum)')
    ax.set_xticks(x); ax.set_xticklabels(PERS); ax.set_ylabel('summed moment (N m)'); ax.axhline(0, color='k', lw=0.5)
    ax.legend(fontsize=7); ax.set_title('Moment of Table-S1 detections vs noise (all five templates)')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2c_moment.png', dpi=130); plt.close(fig)


if __name__ == '__main__':
    main()
