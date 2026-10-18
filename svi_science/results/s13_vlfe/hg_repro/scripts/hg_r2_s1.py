#!/usr/bin/env python
"""Amendment R2a.5: tests on Table S1's OWN detections, using the chosen setting's calibrated thresholds and runs.
(a) Table-S1 counts vs station-day-shift / polarity-scrambled controls and their-style noise templates (same thresholds);
(b) far/near flag and ratio at Table-S1 times vs random times; (c) LFE coincidence +-60 s vs random times;
(d) seasonality: daily/monthly counts vs a microseism/storm proxy (daily median of 10-min RMS, Z, CN.PGC and UW.LRIV,
5-10 s and 20-50 s, from the 1-sps displacement day files); Spearman rho on inter-ETS days for Table S1, our catalogue and
the noise-template detections; (e) clustering between ETS (Fano, gaps < 1 h, rho with tremor) vs our scrambled controls.
Outputs tables/r2_s1_controls.csv, r2_s1_farnear.csv, r2_s1_lfe.csv, r2_s1_season.csv, r2_s1_clustering.csv,
r2_proxy_daily.csv; figures/r2_s1_controls.png, r2_s1_season.png."""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.signal import butter, sosfiltfilt
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hg_common as H
import hg_analyze as A
import hg_controls as C
import hg_lfe as LF
import hg_r2_stage2 as S2

NDRAW = 8
PERS = ['ETS2011', 'interETS', 'ETS2012']
DAYS = pd.date_range('2011-07-01', '2012-12-01', inclusive='left')


def proxy():
    f = f'{H.TB}/r2_proxy_daily.csv'
    if os.path.isfile(f):
        return pd.read_csv(f, index_col=0, parse_dates=True)
    bands = {'5-10s': butter(4, (0.1, 0.2), btype='band', fs=1.0, output='sos'),
             '20-50s': butter(2, (0.02, 0.05), btype='band', fs=1.0, output='sos')}
    out = {}
    for sta in ['CN.PGC', 'UW.LRIV']:
        for b in bands:
            out[f'{sta}_{b}'] = np.full(len(DAYS), np.nan)
        for k, d in enumerate(DAYS):
            p = H.day_file(sta, d)
            if p is None:
                continue
            x = np.load(p)['d'][2].astype(np.float64)
            ok = np.isfinite(x)
            if ok.mean() < 0.5:
                continue
            x = np.where(ok, x - np.nanmean(x), 0.0)
            for b, sos in bands.items():
                y = sosfiltfilt(sos, x)
                w = y.reshape(144, 600); okw = ok.reshape(144, 600).all(1)
                if okw.sum() >= 36:
                    out[f'{sta}_{b}'][k] = np.median(np.sqrt((w[okw] ** 2).mean(1)))
    P = pd.DataFrame(out, index=DAYS)
    P.to_csv(f)
    return P


def daily(i):
    a0 = H.tidx('2011-07-01')
    i = np.asarray(i, np.int64); i = i[(i >= a0) & (i < H.tidx('2012-12-01'))]
    return np.bincount((i - a0) // 86400, minlength=len(DAYS))[:len(DAYS)].astype(float)


def main():
    ch = json.load(open(f'{H.TB}/r2_choice.json'))
    TH = pd.read_csv(f'{H.TB}/r2_thresholds.csv')
    th = TH[(TH.variant == ch['chosen_variant']) & (TH.target_set == 'S1')].set_index('template').thr
    S = pd.read_csv(f'{H.TB}/hg2019_tableS1.csv')
    AC = pd.read_parquet(f'{H.HG}/data/r2/r2_controls_catalogue.parquet')
    NC = pd.read_csv(f'{H.TB}/r2_noise_counts.csv')
    real, NZ, mask, src = S2.load(ch)
    keepm = (lambda i: ~S2.RC.eqmask(i)) if mask else (lambda i: np.ones(len(i), bool))
    # (a) counts vs controls
    rows = []
    for p in PERS + ['interETS_b']:
        a, b = H.per_idx(p); nd = (b - a) / 86400
        for K in list(range(1, 6)) + ['all']:
            sk = S if K == 'all' else S[S.template == K]
            ns = int(((sk.i >= a) & (sk.i < b)).sum())
            d = AC[(AC.i >= a) & (AC.i < b)] if K == 'all' else AC[(AC.template == K) & (AC.i >= a) & (AC.i < b)]
            nours = int((d.set == 'real').sum())
            a1 = np.array([(d.set == f'A1_{r}').sum() for r in range(NDRAW)])
            a2 = np.array([(d.set == f'A2_{r}').sum() for r in range(NDRAW)])
            nn = NC[(NC.period == p) & ((NC.template == K) if K != 'all' else True)].groupby('j').n.sum().values
            ntpl = 5 if K == 'all' else 1
            e = (ns - nn) / ns if ns else np.full(nn.shape, np.nan)
            rows.append({'period': p, 'template': K, 'days': nd, 'N_S1': ns, 'N_ours': nours,
                         'S1_per_tpl_day': ns / nd / ntpl, 'noise_per_tpl_day': nn.mean() / nd / ntpl,
                         'A1_per_tpl_day': a1.mean() / nd / ntpl, 'A2_per_tpl_day': a2.mean() / nd / ntpl,
                         'excess_vs_noise_mean': (ns - nn.mean()) / ns if ns else np.nan,
                         'excess_vs_noise_p16': np.percentile(e, 16), 'excess_vs_noise_p84': np.percentile(e, 84),
                         'excess_vs_scrambled': (ns - max(a1.mean(), a2.mean())) / ns if ns else np.nan,
                         'scrambled_range': f'A1 {a1.min()}-{a1.max()}, A2 {a2.min()}-{a2.max()}'})
    RA = pd.DataFrame(rows); RA.to_csv(f'{H.TB}/r2_s1_controls.csv', index=False)
    pd.set_option('display.width', 250)
    print(RA.round(3).to_string(index=False))
    # (b) far/near at Table-S1 times vs random
    rng = np.random.default_rng(H.SEED + 51)
    fr = []
    for K in range(1, 6):
        chans = real[K]['chans'] if 'chans' in real[K] else real[K]['prim']
        z = C.zstas_of(chans)
        for kind, t in [('Table S1', S[S.template == K].i.values),
                        ('random times', rng.integers(H.tidx('2011-07-01'), H.tidx('2012-12-01'), 3000))]:
            f = C.farnear(t, z); f['template'] = K; f['kind'] = kind; f['period'] = __import__('hg_r2_match').per_name(np.asarray(t))
            fr.append(f)
    FR = pd.concat(fr, ignore_index=True)
    FS = FR[FR.refs_valid].groupby(['kind', 'period']).agg(n=('flag_far', 'size'), flag_frac=('flag_far', 'mean'),
                                                          excess_med=('ref_excess_med', 'median'),
                                                          ratio_med=('ref_ratio', 'median')).reset_index()
    FS.to_csv(f'{H.TB}/r2_s1_farnear.csv', index=False); print(FS.round(3).to_string(index=False))
    # (c) LFE coincidence
    L = LF.lfe_times()
    lr = []
    for p in PERS:
        a, b = H.per_idx(p)
        rnd = rng.integers(a, b, 20000)
        for K in list(range(1, 6)) + ['all']:
            sk = S if K == 'all' else S[S.template == K]
            t = sk[(sk.i >= a) & (sk.i < b)].i.values
            for cat in ['both', 'SVI', 'WA_NW']:
                c, cr = LF.counts(L[cat], t), LF.counts(L[cat], rnd)
                lr.append({'period': p, 'template': K, 'lfe_cat': cat, 'n_S1': len(t),
                           'frac_S1_with_lfe': (c > 0).mean() if len(c) else np.nan, 'frac_random_with_lfe': (cr > 0).mean(),
                           'mean_lfe_S1': c.mean() if len(c) else np.nan, 'mean_lfe_random': cr.mean(),
                           'ratio': c.mean() / cr.mean() if len(c) and cr.mean() > 0 else np.nan})
    LR = pd.DataFrame(lr); LR.to_csv(f'{H.TB}/r2_s1_lfe.csv', index=False)
    print(LR[LR.lfe_cat == 'both'].round(4).to_string(index=False))
    # (d) seasonality
    PX = proxy()
    ours = pd.read_csv(f'{H.TB}/r2_catalogue_chosen.csv')
    nz_daily = {K: np.zeros(len(DAYS)) for K in range(1, 6)}
    for dz in NZ:
        for K in range(1, 6):
            pi, pv = dz[K]
            s = (pv > th[K]) & (np.abs(pi - dz['i0']) > H.L)
            s &= keepm(pi)
            nz_daily[K] += daily(pi[s]) / len(NZ)
    a, b = H.per_idx('interETS')
    di = (DAYS >= H.PER['interETS'][0]) & (DAYS < H.PER['interETS'][1])
    sr = []
    series = {}
    for K in list(range(1, 6)) + ['all']:
        series[('Table S1', K)] = daily(S.i.values if K == 'all' else S[S.template == K].i.values)
        series[('ours', K)] = daily(ours.i.values if K == 'all' else ours[ours.template == K].i.values)
        series[('noise templates (mean)', K)] = sum(nz_daily.values()) if K == 'all' else nz_daily[K]
    for (kind, K), y in series.items():
        for col in PX.columns:
            x = PX[col].values
            ok = di & np.isfinite(x)
            rho = spearmanr(y[ok], x[ok]).correlation if y[ok].std() > 0 else np.nan
            sr.append({'series': kind, 'template': K, 'proxy': col, 'rho_interETS': rho, 'n_days': int(ok.sum())})
    SR = pd.DataFrame(sr); SR.to_csv(f'{H.TB}/r2_s1_season.csv', index=False)
    print(SR.pivot_table(index=['series', 'template'], columns='proxy', values='rho_interETS').round(2).to_string())
    # (e) clustering between ETS: Table S1 vs our scrambled controls
    TT = H.templates()
    cl = []
    for K in range(1, 6):
        r = TT[TT.template == K].iloc[0]
        trem = A.tremor_daily(r.lat, r.lon, 100)
        trem = trem[(trem.index >= H.PER['interETS'][0]) & (trem.index < H.PER['interETS'][1])].values
        sets = {'Table S1': S[S.template == K].i.values, 'ours': ours[ours.template == K].i.values}
        for q in range(NDRAW):
            sets[f'A1_{q}'] = AC[(AC.template == K) & (AC.set == f'A1_{q}')].i.values
            sets[f'A2_{q}'] = AC[(AC.template == K) & (AC.set == f'A2_{q}')].i.values
        for name, x in sets.items():
            o = C.clustering(np.asarray(x, np.int64), a, b, trem)
            o.update({'template': K, 'set': name, 'kind': name.split('_')[0] if name.startswith('A') else name})
            cl.append(o)
    CL = pd.DataFrame(cl); CL.to_csv(f'{H.TB}/r2_s1_clustering.csv', index=False)
    print(CL.groupby(['template', 'kind'])[['n', 'fano', 'frac_iet_lt1h', 'poisson_frac_lt1h', 'rho_tremor100']].mean().round(3))
    figs(series, PX, RA, ch)


def figs(series, PX, RA, ch):
    tr = A.tremor_daily(H.templates().lat.mean(), H.templates().lon.mean(), 100).reindex(DAYS, fill_value=0)
    fig, ax = plt.subplots(4, 1, figsize=(12, 11), sharex=True)
    for K, c in zip(range(1, 6), ['k', 'r', 'b', 'g', 'orange']):
        y = pd.Series(series[('Table S1', K)], index=DAYS).resample('MS').sum()
        ax[0].plot(y.index, y.values, color=c, marker='o', ms=3, label=f'T{K}')
    ax[0].set_ylabel('Table S1 detections / month'); ax[0].legend(fontsize=7, ncol=5)
    for K, c in zip(range(1, 6), ['k', 'r', 'b', 'g', 'orange']):
        y = pd.Series(series[('ours', K)], index=DAYS).resample('MS').sum()
        ax[1].plot(y.index, y.values, color=c, marker='o', ms=3, label=f'T{K}')
    yn = pd.Series(series[('noise templates (mean)', 'all')], index=DAYS).resample('MS').sum()
    ax[1].plot(yn.index, yn.values, color='c', lw=2, label='noise templates, 5 configs (mean of 100)')
    ax[1].set_ylabel(f'ours ({ch["chosen_variant"]}) / month'); ax[1].legend(fontsize=7, ncol=3)
    for col, c in zip(PX.columns, ['navy', 'skyblue', 'darkred', 'salmon']):
        ax[2].semilogy(PX.index, PX[col], color=c, lw=0.7, label=col)
    ax[2].set_ylabel('daily median 10-min RMS (m)'); ax[2].legend(fontsize=7, ncol=2)
    ax[3].bar(tr.index, tr.values, width=1, color='#3b4cc0'); ax[3].set_ylabel('PNSN tremor/day ≤100 km')
    for g in ax:
        for p in ['ETS2011', 'ETS2012']:
            g.axvspan(*pd.to_datetime(H.PER[p]), color='0.88', zorder=0)
    ax[3].set_xlim(DAYS[0], DAYS[-1])
    fig.suptitle('Table S1 vs ours vs noise templates, monthly, with a microseism/storm proxy (PGC, LRIV Z)')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2_s1_season.png', dpi=120); plt.close(fig)
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.8), sharey=True)
    for g, p in zip(ax, PERS):
        x = RA[(RA.period == p) & (RA.template != 'all')]
        k = np.arange(len(x))
        for j, (col, c, lab) in enumerate([('S1_per_tpl_day', 'k', 'Table S1'), ('noise_per_tpl_day', 'c', 'noise templates'),
                                           ('A1_per_tpl_day', 'g', 'station-day-shift'), ('A2_per_tpl_day', 'purple', 'polarity-scrambled')]):
            g.bar(k + (j - 1.5) * 0.2, x[col], 0.2, color=c, label=lab)
        g.set_xticks(k); g.set_xticklabels([f'T{v}' for v in x.template]); g.set_title(p); g.set_yscale('log')
    ax[0].set_ylabel('detections per template per day'); ax[0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{H.FIG}/r2_s1_controls.png', dpi=130); plt.close(fig)


if __name__ == '__main__':
    main()
