#!/usr/bin/env python
"""VLFE scan pilot raw figures (re-drawn): z(t) 2011 and 2011-08 with masks/ETS/H&G ticks; separation curves incl. the
post-hoc leave-dominant-station-out variant. Points coloured by network CC; no error bars."""
import os, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scan_common as S
R_ = f'{S.SD}/runs'; TB = S.C.TB; FG = S.C.FIG
zmax = np.load(f'{R_}/real/zmax.npy'); mh = S.mask_vector('head')
det = pd.read_csv(f'{TB}/scan_detections_2011.csv'); det = det[det.sign > 0]
lo = pd.read_parquet(f'{R_}/real/loso.parquet'); lo = lo[(lo.variant == 0)]
fig, ax = plt.subplots(2, 1, figsize=(13, 7.5))
for a_, (t0, t1), n, ttl in ((ax[0], (S.Y0, S.Y0 + pd.Timedelta(days=365)), 3600, '2011 (hourly max)'),
                             (ax[1], (pd.Timestamp('2011-08-10'), pd.Timestamp('2011-09-06')), 600, '2011-08 ETS (10-min max)')):
    z = np.where(np.isfinite(zmax), zmax, -np.inf)
    zu = np.where(mh, -np.inf, z).reshape(-1, n).max(axis=1); zm = np.where(~mh, -np.inf, z).reshape(-1, n).max(axis=1)
    tt = S.Y0 + pd.to_timedelta(np.arange(len(zu)) * n, unit='s'); sel = (tt >= t0) & (tt < t1)
    a_.plot(tt[sel], np.clip(zm[sel], 0, 60), color='tab:orange', lw=0.4, alpha=0.6, label='max z, inside earthquake mask')
    a_.plot(tt[sel], np.clip(zu[sel], 0, 60), color='0.35', lw=0.4, label='max z, unmasked')
    a_.axvspan(S.ETS[0], S.ETS[1], color='tab:blue', alpha=0.08, label='SVI ETS 2011-08-12…09-03')
    td = pd.to_datetime(det.time_utc); sd = (td >= t0) & (td < t1)
    sc = a_.scatter(td[sd], np.clip(det.z[sd], 0, 60), c=det.cc[sd], cmap='viridis', vmin=0.3, vmax=1, s=16, zorder=3,
                    label='headline detections (z ≥ 7)')
    tl = pd.to_datetime(lo.time_utc); sl = (tl >= t0) & (tl < t1) & (lo.z_loso >= 7)
    a_.scatter(tl[sl], np.clip(lo.z[sl], 0, 60), facecolor='none', edgecolor='r', s=40, lw=0.8, zorder=4,
               label='…also z ≥ 7 without the dominant station (post hoc)')
    a_.axhline(7, color='k', lw=0.6, ls='--')
    if t1 - t0 < pd.Timedelta(days=60):
        for k, ts_, la, lo_, dep, mw in S.HG:
            a_.axvline(pd.Timestamp(ts_), color='r', lw=1.0, ymax=0.12)
            a_.annotate(f'H&G{k} Mw{mw}', (pd.Timestamp(ts_), 52 - 5 * int(k)), color='r', fontsize=7)
    a_.set_ylim(0, 60); a_.set_ylabel('z, max over grid (clipped 60)'); a_.set_title(f'VLFE scan, S Vancouver Island, {ttl}', fontsize=9)
    a_.legend(fontsize=7, loc='upper left', ncol=3)
fig.colorbar(sc, ax=ax, label='network CC'); fig.savefig(f'{FG}/scan_z_2011.png', dpi=130); plt.close(fig)
# separation curves
SEP = pd.read_csv(f'{TB}/scan_separation.csv'); TL = pd.read_csv(f'{TB}/scan_separation_tail.csv'); LS = pd.read_csv(f'{TB}/scan_loso_separation.csv')
g = SEP[(SEP.tier == 'head') & SEP.filtered]; g = pd.concat([g, TL[TL.thr > 7]], ignore_index=True)
fig, ax = plt.subplots(1, 2, figsize=(12, 4.4))
for a_, d, x, ttl in ((ax[0], g, 'thr', 'pre-registered statistic (headline filter, headline mask)'),
                      (ax[1], LS, 'thr_z_loso', 'POST HOC: z without the dominant station (z ≥ 7 also required)')):
    a_.semilogy(d[x], d.real_pos_per_yr.clip(lower=0.1), 'o-', color='tab:red', label='real, positive')
    if 'real_neg_per_yr' in d and d.real_neg_per_yr.notna().any():
        a_.semilogy(d[x], d.real_neg_per_yr.clip(lower=0.1), 's-', color='tab:blue', label='real, negative')
    rb = 'real_pos_between_per_yr' if 'real_pos_between_per_yr' in d else 'real_between_per_yr'
    a_.semilogy(d[x], d[rb].clip(lower=0.1), 'o--', color='tab:purple', label='real positive, between ETS')
    a_.semilogy(d[x], d.ctrlA_per_yr.clip(lower=0.1), '^-', color='0.3', label='(a) station-day shifted (mean of 5)')
    a_.semilogy(d[x], d.ctrlB_between_per_yr.clip(lower=0.1), 'v-', color='tab:green', label='(b) polarity scrambled, between ETS')
    a_.axhline(1, color='k', lw=0.5, ls=':'); a_.axhline(2, color='k', lw=0.5, ls=':')
    a_.set_xlabel('threshold'); a_.set_ylabel('clusters per live year'); a_.set_title(ttl, fontsize=9); a_.legend(fontsize=7)
fig.tight_layout(); fig.savefig(f'{FG}/scan_separation.png', dpi=130); plt.close(fig)
