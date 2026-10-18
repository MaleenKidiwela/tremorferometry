#!/usr/bin/env python
"""Amendment Y9b part 1: the 5 between-ETS candidates. -> tables/y9b_candidates.csv, tables/y9b_hf_<k>.csv,
figures/y9b_cand_<k>.png; plus corrected patch-LFE counts per 2020-2021 episode -> tables/y9b_patch_lfe_2020_2021.csv"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe/scripts')
import numpy as np, pandas as pd
import cat_common as K
import s13_common as C

CANDS = ['2014-05-05 08:30:28', '2015-04-13 06:06:11', '2016-06-13 07:41:24', '2016-06-13 19:58:13', '2017-01-31 04:22:00']
HFSTA = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV', 'CN.LZB', 'UW.DOSE', 'UW.GNW', 'UW.SP2', 'UW.WISH']
PATCH = pd.read_csv(f'{K.TB}/inj_sources.csv').set_index('source').loc[['C0', 'C1'], ['lat', 'lon']]


def dk(a, b, c, d):
    return float(K.dkm(a, b, c, d))


def hf(sta, ep):
    """RMS in 2-8 Hz and 1-15 Hz (Z) in [t0-30, t0+150] over the median of 3-min windows in +-30 min (raw counts)."""
    from obspy import read, UTCDateTime
    t0 = UTCDateTime(ep); d = pd.Timestamp(ep, unit='s')
    fs = C.day_files(sta, d.strftime('%Y-%m-%d'))
    if not fs:
        return None
    try:
        st = read(fs[0], starttime=t0 - 1900, endtime=t0 + 1900).select(component='Z')
    except Exception:
        return None
    if not len(st):
        return None
    st.merge(fill_value=0); tr = st[0]; sr = tr.stats.sampling_rate
    if sr < 20:
        return None
    out = {}
    for lab, (f1, f2) in (('2_8', (2, 8)), ('1_15', (1, min(15, 0.45 * sr)))):
        x = tr.copy(); x.detrend('demean'); x.taper(0.02); x.filter('bandpass', freqmin=f1, freqmax=f2, corners=4, zerophase=True)
        y = x.data.astype(float); n0 = int((t0 - 1800 - x.stats.starttime) * sr); w = int(180 * sr)
        segs = [np.sqrt(np.mean(y[n0 + k * w:n0 + (k + 1) * w] ** 2)) for k in range(20) if n0 + (k + 1) * w <= len(y)]
        a = int((t0 - 30 - x.stats.starttime) * sr); win = np.sqrt(np.mean(y[a:a + w] ** 2))
        out[lab] = float(win / np.median(segs)) if segs and np.median(segs) > 0 else np.nan
    return out


def main():
    import scan_refs as SR
    cc = SR.all_coords()
    R = pd.read_csv(f'{K.TB}/y9_level5_between.csv')
    EQ = pd.concat([K.eq_table(y) for y in range(2014, 2018)]).drop_duplicates('id')
    Ptr = pd.read_csv(f'{K.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv', usecols=['time', 'lat', 'lon'])
    Ptr['t'] = pd.to_datetime(Ptr.time).astype('datetime64[ns]').astype('int64') / 1e9
    Rw = C.svi_rows(); Rw = Rw[(Rw.series == 'data') & Rw.counted_once]
    F = C.svi_families()
    Rw = Rw.join(F[['lat', 'lon']], on='fam')
    rows = []
    for k, ts in enumerate(CANDS):
        e = R[R.time_pt.str.startswith(ts)].iloc[0]; ep = float(e.epoch)
        q = EQ[(EQ.t >= ep - 300) & (EQ.t <= ep + 300)]
        qn = EQ[(EQ.t >= ep - 3600) & (EQ.t <= ep + 300)].sort_values('t')
        H = {}
        for s in HFSTA:
            h = hf(s, ep)
            if h:
                H[s] = h
        pd.DataFrame(H).T.to_csv(f'{K.TB}/y9b_hf_{k}.csv')
        r28 = [v['2_8'] for v in H.values() if np.isfinite(v['2_8'])]; r115 = [v['1_15'] for v in H.values() if np.isfinite(v['1_15'])]
        tr = Ptr[(abs(Ptr.t - ep) <= 3600)]
        tr30 = tr[[dk(a, b, e.lat_best, e.lon_best) <= 30 for a, b in zip(tr.lat, tr.lon)]] if len(tr) else tr
        lf = Rw[abs(Rw.t - ep) <= 3600]
        lfd = np.array([dk(a, b, e.lat_best, e.lon_best) if np.isfinite(a) else np.nan for a, b in zip(lf.lat, lf.lon)]) if len(lf) else np.array([])
        o = dict(k=k, time_pt=e.time_pt, lat=e.lat_best, lon=e.lon_best, dist_C0_km=dk(e.lat_best, e.lon_best, *PATCH.loc['C0']),
                 dist_C1_km=dk(e.lat_best, e.lon_best, *PATCH.loc['C1']), nearest_family=e.nearest_family, dist_family_km=e.dist_family_km,
                 eq_pm5min=';'.join(f"{x.time_pt} M{x.mag:.1f} {x.dist_km:.0f}km {x.src}" for x in q.itertuples()),
                 eq_prev1h=';'.join(f"{x.time_pt} M{x.mag:.1f} {x.dist_km:.0f}km" for x in qn.tail(3).itertuples()),
                 hf28_n_ge3=int(sum(np.array(r28) >= 3)), hf28_med=float(np.median(r28)) if r28 else np.nan, hf28_max=float(np.max(r28)) if r28 else np.nan,
                 hf115_n_ge3=int(sum(np.array(r115) >= 3)), hf115_med=float(np.median(r115)) if r115 else np.nan, n_hf_sta=len(r28),
                 tremor_pm1h_all=len(tr), tremor_pm1h_30km=len(tr30),
                 lfe075_pm1h=int((lf.P >= 0.75).sum()), lfe075_pm1h_30km=int(((lf.P >= 0.75).values & (lfd <= 30)).sum()) if len(lf) else 0,
                 lfecand_pm1h=int(((lf.P >= 0.1) & (lf.P < 0.75)).sum()),
                 lfecand_pm1h_30km=int((((lf.P >= 0.1) & (lf.P < 0.75)).values & (lfd <= 30)).sum()) if len(lf) else 0,
                 emerg_med_near=e.emerg_med_near, emerg_far_med=e.emerg_far_med, t3_spearman=e.t3_spearman, ratio_far=e.ratio_far,
                 label_eye_y9=e.label_eye)
        cls = 'regular earthquake' if (o['eq_pm5min'] or o['hf28_n_ge3'] >= 2) else \
              ('VLFE-like' if (o['hf28_n_ge3'] == 0 and e.t3_spearman <= -0.4 and e.emerg_far_med < e.emerg_med_near) else 'unclear')
        o['class_rule'] = cls
        rows.append(o)
        fig_cand(k, e, H, tr30, lf[lfd <= 30] if len(lf) else lf, cc)
    T = pd.DataFrame(rows); T.to_csv(f'{K.TB}/y9b_candidates.csv', index=False)
    pd.set_option('display.width', 250); pd.set_option('display.max_columns', 60)
    print(T.drop(columns=['eq_prev1h']).round(2).T.to_string())
    print(T[['time_pt', 'eq_prev1h']].to_string())
    patch_counts(Rw, F)


def fig_cand(k, e, H, tr30, lf30, cc):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from obspy import read, UTCDateTime
    ep = float(e.epoch); t0 = UTCDateTime(ep); d = pd.Timestamp(ep, unit='s')
    fig, axs = plt.subplots(1, 3, figsize=(16, 7), gridspec_kw=dict(width_ratios=[1.2, 1.2, 0.8]))
    sts = sorted(H, key=lambda s: dk(cc.loc[s, 'lat'], cc.loc[s, 'lon'], e.lat_best, e.lon_best) if s in cc.index else 999)
    for j, s in enumerate(sts):
        fs = C.day_files(s, d.strftime('%Y-%m-%d'))
        try:
            st = read(fs[0], starttime=t0 - 600, endtime=t0 + 600).select(component='Z'); st.merge(fill_value=0); tr = st[0]
        except Exception:
            continue
        for ax, (f1, f2) in zip(axs[:2], ((0.02, 0.05), (2, 8))):
            x = tr.copy(); x.detrend('demean'); x.taper(0.05); x.filter('bandpass', freqmin=f1, freqmax=f2, corners=4, zerophase=True)
            y = x.data.astype(float); t = np.arange(len(y)) / x.stats.sampling_rate - 600 + (x.stats.starttime - (t0 - 600))
            if f1 >= 1:
                from scipy.signal import hilbert
                y = np.abs(hilbert(y))
            ax.plot(t, y / (2.2 * np.nanmax(np.abs(y))) - j, 'k', lw=0.5)
            ax.text(605, -j, f"{s} {dk(cc.loc[s, 'lat'], cc.loc[s, 'lon'], e.lat_best, e.lon_best):.0f} km" + (f"  HF {H[s]['2_8']:.1f}" if f1 >= 1 else ''), fontsize=7)
    for ax, lab in zip(axs[:2], ('Z 0.02-0.05 Hz (counts, own-peak scaled)', 'Z 2-8 Hz envelope (own-peak scaled)')):
        ax.axvspan(0, 180, color='0.9', zorder=0); ax.set_xlim(-600, 760); ax.set_yticks([]); ax.set_title(lab, fontsize=9); ax.set_xlabel('s after t0')
    ax = axs[2]
    ax.scatter(PATCH.lon, PATCH.lat, marker='s', c='tab:blue', s=60, label='C0/C1 patches')
    if len(tr30):
        ax.scatter(tr30.lon, tr30.lat, c='orange', s=8, label=f'PNSN tremor ±1 h ≤30 km ({len(tr30)})')
    if len(lf30):
        ax.scatter(lf30.lon, lf30.lat, c='green', s=8, label=f'SVI LFE det. ±1 h ≤30 km ({len(lf30)})')
    ax.scatter([e.lon_best], [e.lat_best], marker='*', c='r', s=200, label='scan location')
    for s in sts:
        if s in cc.index:
            ax.plot(cc.loc[s, 'lon'], cc.loc[s, 'lat'], 'k^', ms=5); ax.text(cc.loc[s, 'lon'], cc.loc[s, 'lat'], s.split('.')[1], fontsize=6)
    ax.set_xlim(-125.2, -122.3); ax.set_ylim(47.4, 49.4); ax.legend(fontsize=7, loc='lower left'); ax.set_aspect(1 / np.cos(np.radians(48.3)))
    fig.suptitle(f'Candidate {k}: {e.time_pt}  MF {e.zmad:.1f}  scan z {e.z_best:.1f}  emergence near {e.emerg_med_near:.1f} far {e.emerg_far_med:.1f}  '
                 f'Spearman {e.t3_spearman:.2f}', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{K.FIG}/y9b_cand_{k}.png', dpi=95); plt.close(fig)


def patch_counts(Rw, F):
    """corrected patch-LFE counts: families within 10 km of the C0 or C1 centroid, P >= 0.75 (and candidates 0.1-0.75)."""
    near = [f for f, r in F.iterrows() if min(dk(r.lat, r.lon, *PATCH.loc['C0']), dk(r.lat, r.lon, *PATCH.loc['C1'])) <= 10]
    E = pd.read_csv(f'{K.ROOT}/svi_science/results/b_ets_record/episodes_b.csv'); E = E[(E.region == 'SVI') & (E.start >= '2020') & (E.start < '2022')]
    x = Rw[Rw.fam.isin(near)]; rows = []
    for r in E.itertuples():
        a = pd.Timestamp(r.start).value / 1e9 - 2 * 86400; b = pd.Timestamp(r.end).value / 1e9 + 3 * 86400
        w = x[(x.t >= a) & (x.t < b)]; wa = Rw[(Rw.t >= a) & (Rw.t < b)]
        rows.append(dict(episode=r.episode, major=r.major, patch_lfe075=int((w.P >= 0.75).sum()), patch_cand=int(((w.P >= 0.1) & (w.P < 0.75)).sum()),
                         all_svi_lfe075=int((wa.P >= 0.75).sum()), all_svi_cand=int(((wa.P >= 0.1) & (wa.P < 0.75)).sum())))
    T = pd.DataFrame(rows); T.to_csv(f'{K.TB}/y9b_patch_lfe_2020_2021.csv', index=False)
    print('patch families (<= 10 km of C0/C1):', len(near)); print(T.to_string(index=False))


if __name__ == '__main__':
    main()
