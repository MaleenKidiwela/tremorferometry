#!/usr/bin/env python
"""Step 5 (PLAN.md rows A1-A5, §3): our controls applied to the reproduced H&G catalogue.
- A1 station-day-shift and A2 polarity-scrambled detections: peaks of the A1_r / A2_r sums (hg_mf.py real) above the
  same H&G threshold (thrA, primary config), counted by output time per period.
- A3 far/near: refs data/mf/ref_*.npy (LHZ displacement, 0.02-0.05 Hz). Per time t: per ref, RMS over [t-60, t+330) and
  over [t-3600, t+3600) minus that window (live if >= 90 % finite in both); ref_excess = ratio; FLAG if >= 2 live refs
  have excess >= 3. ref_ratio = median ref window RMS / median near-Z RMS on [t, t+90) of the config. GNW excess column.
  Only t < 2013-01-01 (refs end 2012-12-31).
- A4 clustering (inter-ETS): Fano factor of daily counts, fraction of inter-event times < 1 h (and the Poisson value),
  Spearman rho of daily counts with PNSN tremor within 100 km of the template.
- A5 masks v2 / v1 (../data/scan*/masks2_*.npz, True = masked) as columns.
Survival (PLAN §3): surv = max(0, N_real,unflagged - max(mean A1, mean A2)_unflagged) / N_real.
Outputs tables/controls_counts.csv, controls_survival.csv, controls_clustering.csv, catalogue_prim_flags.csv,
controls_calibration.csv; figures/controls_daily.png, controls_separation.png, controls_farnear.png."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import hg_common as H
import hg_analyze as A

REFS = ['IU.COR', 'US.NEW', 'US.HAWA', 'UW.OMAK', 'TA.I02D']
W0, W1, BG = -60, 330, 3600
IEND = H.tidx('2013-01-01')
NDRAW = 8


class Win:
    """fast window RMS on a 1-sps series with NaN."""
    def __init__(self, x):
        x = np.asarray(x)
        f = np.isfinite(x)
        y = np.where(f, x.astype(np.float64), 0.0)
        cap = 50 * np.percentile(np.abs(y[f]), 99.9)    # glitch cap (0.8 m TA glitches would kill cumsum precision)
        np.clip(y, -cap, cap, out=y)
        self.c2 = np.concatenate([[0.0], np.cumsum(y * y)])
        self.cn = np.concatenate([[0], np.cumsum(f)])
        self.n = x.size

    def sums(self, a, b):
        a = np.clip(a, 0, self.n); b = np.clip(b, 0, self.n)
        return self.c2[b] - self.c2[a], self.cn[b] - self.cn[a]


_REF = None
_NEAR = {}


def refs():
    global _REF
    if _REF is None:
        _REF = {s: Win(np.load(f'{H.MF}/ref_{s}.npy')) for s in REFS}
    return _REF


def near(sta):
    if sta not in _NEAR:
        _NEAR[sta] = Win(np.load(f'{H.MF}/cont_{sta}.npy', mmap_mode='r')[2])
    return _NEAR[sta]


def excess(W, t):
    sw, nw = W.sums(t + W0, t + W1)
    sb, nb = W.sums(t - BG, t + BG)
    sb, nb = sb - sw, nb - nw
    ok = (nw >= 0.9 * (W1 - W0)) & (nb >= 0.9 * (2 * BG - (W1 - W0)))
    rw = np.sqrt(sw / np.maximum(nw, 1)); rb = np.sqrt(sb / np.maximum(nb, 1))
    ex = np.where(ok & (rb > 0), rw / np.where(rb > 0, rb, 1), np.nan)
    return ex, np.where(ok, rw, np.nan)


def farnear(t, zstas):
    t = np.asarray(t, np.int64)
    EX, RW = [], []
    for s in REFS:
        e, r = excess(refs()[s], t)
        EX.append(e); RW.append(r)
    EX, RW = np.array(EX), np.array(RW)
    nlive = np.isfinite(EX).sum(0)
    n3 = (EX >= 3).sum(0)
    NR = []
    for s in zstas:
        sw, nw = near(s).sums(t, t + H.L)
        NR.append(np.where(nw >= 0.9 * H.L, np.sqrt(sw / np.maximum(nw, 1)), np.nan))
    NR = np.array(NR)
    with np.errstate(all='ignore'):
        ratio = np.nanmedian(RW, 0) / np.nanmedian(NR, 0)
        exg, _ = excess(near('UW.GNW'), t)
    valid = t < IEND
    out = pd.DataFrame({'ref_nlive': nlive, 'ref_n_excess3': n3, 'ref_excess_med': np.nanmedian(EX, 0),
                        'ref_ratio': ratio, 'gnw_excess': exg})
    out['flag_far'] = (out.ref_n_excess3 >= 2) & valid
    out.loc[~valid, ['ref_excess_med', 'ref_ratio']] = np.nan
    out['refs_valid'] = valid & (nlive >= 2)
    return out


_MASK = {}


def mask(t, key='v2'):
    out = np.zeros(len(t), bool)
    for y in (2011, 2012):
        if y not in _MASK:
            f = f'{H.S13}/data/scan/masks2_2011.npz' if y == 2011 else f'{H.S13}/data/scan_2012/masks2_2012.npz'
            _MASK[y] = dict(np.load(f))
        sec = (H.T0 - pd.Timestamp(f'{y}-01-01')).total_seconds()
        s = np.asarray(t, np.int64) + int(sec)
        m = (s >= 0) & (s < _MASK[y][key].size)
        out[m] = _MASK[y][key][s[m]]
    return out


def zstas_of(chans):
    return sorted({c.rsplit('.', 1)[0] for c in chans if c.endswith('.Z')})


def ctrl_cat(R, TH, K, name, thk):
    thr = float(TH[(TH.template == K) & (TH.variant == 'prim')][thk].iloc[0])
    pi, pv = R[K][f'{name}_pi'], R[K][f'{name}_pv']
    s = pv > thr
    return pd.DataFrame({'template': K, 'set': name, 'i': pi[s], 'S': pv[s]})


def clustering(i, a, b, trem):
    i = np.sort(i[(i >= a) & (i < b)])
    nd = (b - a) // 86400
    daily = np.bincount((i - a) // 86400, minlength=nd)[:nd]
    lam = len(i) / (b - a)
    iet = np.diff(i)
    rho = spearmanr(daily, trem[:nd]).correlation if daily.std() > 0 else np.nan
    return {'n': len(i), 'fano': daily.var() / daily.mean() if daily.mean() > 0 else np.nan,
            'frac_iet_lt1h': (iet < 3600).mean() if iet.size else np.nan,
            'poisson_frac_lt1h': 1 - np.exp(-lam * 3600), 'rho_tremor100': rho}


def main(thk, R, NZ):
    TH = pd.read_csv(f'{H.TB}/thresholds.csv')
    CAT = pd.read_csv(f'{H.TB}/catalogue_prim_{thk}.csv')
    TT = H.templates()
    # --- real + control catalogues with flags
    allc = [CAT.assign(set='real')[['template', 'set', 'i', 'S']]]
    for K in range(1, 6):
        for r in range(NDRAW):
            allc += [ctrl_cat(R, TH, K, f'A1_{r}', thk), ctrl_cat(R, TH, K, f'A2_{r}', thk)]
    AC = pd.concat(allc, ignore_index=True)
    fl = []
    for K in range(1, 6):
        m = AC.template == K
        f = farnear(AC[m].i.values, zstas_of(R[K]['prim']))
        f.index = AC[m].index
        fl.append(f)
    AC = AC.join(pd.concat(fl))
    AC['mask_v2'] = mask(AC.i.values, 'v2'); AC['mask_v1'] = mask(AC.i.values, 'v1')
    AC['time_utc'] = H.itime(AC.i.values)
    AC[AC.set == 'real'].to_csv(f'{H.TB}/catalogue_prim_flags_{thk}.csv', index=False)
    AC.to_parquet(f'{H.MF}/controls_catalogue_{thk}.parquet')
    # --- calibration: templates themselves, random times, noise-template self times
    rng = np.random.default_rng(H.SEED + 11)
    cal = []
    for K in range(1, 6):
        z = zstas_of(R[K]['prim'])
        f = farnear([int(R[K]['i0'])], z).assign(kind='template event', template=K)
        cal.append(f)
        rt = rng.integers(H.tidx('2011-07-01'), IEND - 400, 3000)
        cal.append(farnear(rt, z).assign(kind='random times', template=K))
    CAL = pd.concat(cal, ignore_index=True)
    CS = CAL[CAL.refs_valid].groupby(['kind', 'template']).agg(n=('flag_far', 'size'), flag_frac=('flag_far', 'mean'),
                                                               excess_med=('ref_excess_med', 'median'),
                                                               ratio_med=('ref_ratio', 'median')).reset_index()
    CS.to_csv(f'{H.TB}/controls_calibration.csv', index=False)
    print(CS.round(3).to_string(index=False))
    # --- counts per period
    rows = []
    for K in range(1, 6):
        for p in ['ETS2011', 'interETS', 'interETS_b', 'ETS2012', 'all18']:
            a, b = H.per_idx(p)
            d = AC[(AC.template == K) & (AC.i >= a) & (AC.i < b)]
            nd = (b - a) / 86400
            def cnt(sets, unfl=False, nomask=False):
                v = []
                for s in sets:
                    x = d[d.set == s]
                    if unfl:
                        x = x[~x.flag_far]
                    if nomask:
                        x = x[~x.mask_v2]
                    v.append(len(x))
                return np.array(v, float)
            r_ = cnt(['real'])[0]; ru = cnt(['real'], True)[0]; rm = cnt(['real'], True, True)[0]
            a1 = cnt([f'A1_{q}' for q in range(NDRAW)]); a2 = cnt([f'A2_{q}' for q in range(NDRAW)])
            a1u = cnt([f'A1_{q}' for q in range(NDRAW)], True); a2u = cnt([f'A2_{q}' for q in range(NDRAW)], True)
            a1m = cnt([f'A1_{q}' for q in range(NDRAW)], True, True); a2m = cnt([f'A2_{q}' for q in range(NDRAW)], True, True)
            ctrl = max(a1u.mean(), a2u.mean())
            surv = max(0.0, ru - ctrl)
            ctrl_m = max(a1m.mean(), a2m.mean())
            rows.append({'template': K, 'period': p, 'days': nd, 'N_real': r_, 'N_real_unflagged': ru,
                         'A1_mean': a1.mean(), 'A1_range': f'{a1.min():.0f}-{a1.max():.0f}',
                         'A2_mean': a2.mean(), 'A2_range': f'{a2.min():.0f}-{a2.max():.0f}',
                         'A1_unfl_mean': a1u.mean(), 'A2_unfl_mean': a2u.mean(),
                         'survivors': surv, 'surv_frac': surv / r_ if r_ else np.nan,
                         'surv_frac_noA3': max(0.0, r_ - max(a1.mean(), a2.mean())) / r_ if r_ else np.nan,
                         'surv_frac_unmasked_v2': max(0.0, rm - ctrl_m) / rm if rm else np.nan,
                         'N_real_unmasked_v2_unfl': rm,
                         'real_flag_frac': 1 - ru / r_ if r_ else np.nan})
    CC = pd.DataFrame(rows)
    pooled = CC.groupby('period')[['days', 'N_real', 'N_real_unflagged', 'A1_mean', 'A2_mean', 'A1_unfl_mean', 'A2_unfl_mean',
                                   'N_real_unmasked_v2_unfl']].sum()
    pooled['days'] /= 5
    # pooled survival = sum over templates of per-template survivors (each template has its own control)
    pooled['survivors'] = CC.groupby('period').survivors.sum()
    pooled['surv_frac'] = pooled.survivors / pooled.N_real
    pooled['rate_real_per_tpl_day'] = pooled.N_real / pooled.days / 5
    pooled['rate_A1_per_tpl_day'] = pooled.A1_mean / pooled.days / 5
    pooled['rate_A2_per_tpl_day'] = pooled.A2_mean / pooled.days / 5
    CC.to_csv(f'{H.TB}/controls_counts_{thk}.csv', index=False)
    pooled.reset_index().to_csv(f'{H.TB}/controls_survival_{thk}.csv', index=False)
    print(CC.round(3).to_string(index=False)); print(pooled.round(3).to_string())
    # --- clustering
    cl = []
    a, b = H.per_idx('interETS')
    for K in range(1, 6):
        r = TT[TT.template == K].iloc[0]
        trem = A.tremor_daily(r.lat, r.lon, 100)
        trem = trem[(trem.index >= H.PER['interETS'][0]) & (trem.index < H.PER['interETS'][1])].values
        for s in ['real'] + [f'A1_{q}' for q in range(NDRAW)] + [f'A2_{q}' for q in range(NDRAW)]:
            x = AC[(AC.template == K) & (AC.set == s)].i.values
            o = clustering(x, a, b, trem); o.update({'template': K, 'set': s, 'kind': s.split('_')[0]})
            cl.append(o)
    CL = pd.DataFrame(cl)
    CL.to_csv(f'{H.TB}/controls_clustering_{thk}.csv', index=False)
    print(CL.groupby(['template', 'kind'])[['n', 'fano', 'frac_iet_lt1h', 'poisson_frac_lt1h', 'rho_tremor100']].mean().round(3))
    figs(AC, NZ, TH, R, CAL, thk)


def figs(AC, NZ, TH, R, CAL, thk):
    days = pd.date_range('2011-07-01', '2012-12-01', inclusive='left')
    a18, b18 = H.per_idx('all18')
    TT = H.templates()
    tr = A.tremor_daily(TT.lat.mean(), TT.lon.mean(), 100)
    fig, ax = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
    def daily(sets):
        x = AC[AC.set.isin(sets) & (AC.i >= a18) & (AC.i < b18)]
        return x.groupby(x.time_utc.dt.floor('D')).size().reindex(days, fill_value=0) / len(sets)   # per draw, 5 templates
    real = daily(['real'])
    a1 = daily([f'A1_{q}' for q in range(NDRAW)]); a2 = daily([f'A2_{q}' for q in range(NDRAW)])
    NR = pd.read_csv(f'{H.TB}/noise_rate.csv')
    bgn = NR[(NR.variant == 'prim') & (NR.thr_kind == thk) & (NR.period == 'all18')].noise_rate_per_day.sum()
    ax[0].plot(real.index, real.values, 'r', lw=0.8, label='real templates (H&G-style catalogue)')
    ax[0].plot(a2.index, a2.values, color='purple', lw=0.8, label='polarity-scrambled templates (mean of 8)')
    ax[0].plot(a1.index, a1.values, color='green', lw=0.8, label='station-day-shift (mean of 8)')
    ax[0].axhline(bgn, color='c', lw=1.5, label=f'H&G-style noise templates {bgn:.1f}/day')
    ax[0].set_ylabel('detections/day (5 templates)'); ax[0].legend(fontsize=8, ncol=2)
    ru = AC[(AC.set == 'real') & ~AC.flag_far & (AC.i >= a18) & (AC.i < b18)]
    ru = ru.groupby(ru.time_utc.dt.floor('D')).size().reindex(days, fill_value=0)
    x = AC[AC.set.str.startswith('A2') & ~AC.flag_far & (AC.i >= a18) & (AC.i < b18)]
    a2u = x.groupby(x.time_utc.dt.floor('D')).size().reindex(days, fill_value=0) / NDRAW
    ax[1].plot(ru.index, ru.rolling(7, center=True).mean(), 'r', label='real, not seen far away (7-day mean)')
    ax[1].plot(a2u.index, a2u.rolling(7, center=True).mean(), color='purple', label='polarity-scrambled, not seen far away')
    ax[1].set_ylabel('detections/day'); ax[1].legend(fontsize=8)
    ax[2].bar(tr.index, tr.values, width=1, color='#3b4cc0')
    ax[2].set_ylabel('PNSN tremor/day ≤100 km')
    for a in ax:
        for p in ['ETS2011', 'ETS2012']:
            a.axvspan(*pd.to_datetime(H.PER[p]), color='0.88', zorder=0)
    ax[2].set_xlim(days[0], days[-1])
    fig.suptitle(f'Reproduced H&G catalogue vs our scrambled controls (same matched filter, same threshold: {thk})')
    fig.tight_layout(); fig.savefig(f'{H.FIG}/controls_daily_{thk}.png', dpi=130); plt.close(fig)
    # separation curves, ETS vs inter-ETS
    fig, ax = plt.subplots(2, 5, figsize=(17, 6.5), sharey='row')
    for row, p in enumerate(['ETS2011', 'interETS']):
        a, b = H.per_idx(p); nd = (b - a) / 86400
        for K in range(1, 6):
            n = len(R[K]['prim']); x = np.linspace(0.05, 0.6, 111)
            def curve(name):
                pi, pv = R[K][f'{name}_pi'], R[K][f'{name}_pv']
                m = (pi >= a) & (pi < b)
                if name == 'real':
                    m &= np.abs(pi - int(R[K]['i0'])) > H.L
                return np.array([(pv[m] / n > t).sum() for t in x]) / nd
            g = ax[row, K - 1]
            for q in range(NDRAW):
                g.semilogy(x, curve(f'A1_{q}'), color='green', lw=0.5, alpha=0.6)
                g.semilogy(x, curve(f'A2_{q}'), color='purple', lw=0.5, alpha=0.6)
            g.semilogy(x, curve('real'), 'r', lw=2)
            thr = float(TH[(TH.template == K) & (TH.variant == 'prim')][thk].iloc[0]) / n
            g.axvline(thr, color='k', ls=':')
            g.set_title(f'T{K} {p}', fontsize=9)
            if row == 1:
                g.set_xlabel('summed CC / n_chan')
        ax[row, 0].set_ylabel('detections per day')
    ax[0, 0].plot([], [], 'r', label='real'); ax[0, 0].plot([], [], 'g', label='station-day-shift')
    ax[0, 0].plot([], [], color='purple', label='polarity-scrambled'); ax[0, 0].plot([], [], 'k:', label='H&G threshold')
    ax[0, 0].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{H.FIG}/controls_separation_{thk}.png', dpi=120); plt.close(fig)
    # far/near
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    v = AC[AC.refs_valid]
    bins = np.logspace(-0.5, 1.5, 50)
    for s, c, lab in [('real', 'r', 'real detections'), ('A2', 'purple', 'polarity-scrambled'), ('A1', 'green', 'day-shift')]:
        x = v[v.set == 'real'] if s == 'real' else v[v.set.str.startswith(s)]
        ax[0].hist(x.ref_excess_med.dropna(), bins=bins, histtype='step', color=c, density=True, label=lab)
        ax[1].hist(np.log10(x.ref_ratio.dropna()), bins=50, histtype='step', color=c, density=True, label=lab)
    rt = CAL[(CAL.kind == 'random times') & CAL.refs_valid]
    ax[0].hist(rt.ref_excess_med.dropna(), bins=bins, histtype='step', color='0.5', density=True, label='random times')
    ax[1].hist(np.log10(rt.ref_ratio.dropna()), bins=50, histtype='step', color='0.5', density=True, label='random times')
    te = CAL[CAL.kind == 'template event']
    for _, r in te.iterrows():
        ax[0].axvline(r.ref_excess_med, color='k', lw=0.8); ax[1].axvline(np.log10(r.ref_ratio), color='k', lw=0.8)
    ax[0].set_xscale('log'); ax[0].set_xlabel('median far-reference excess (window / surrounding hour)')
    ax[1].set_xlabel('log10 far/near amplitude ratio'); ax[0].legend(fontsize=8)
    ax[0].set_title('black lines: the five template events', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{H.FIG}/controls_farnear_{thk}.png', dpi=120); plt.close(fig)


if __name__ == '__main__':
    R, NZ = A.load()
    for thk in (sys.argv[1:] or A.THK):
        print('=====', thk, flush=True)
        main(thk, R, NZ)
