#!/usr/bin/env python
"""newvlfe stage 4 (PLAN.md §4): every surviving 'new' candidate (real, > 15 km from all template locations) of a tier
that passes; if no tier passes, the real 'new' B2 rows (z >= 5) are listed as candidates. The same columns are computed
for the control rows of the same selection (for comparison).
Columns: location, nearest LFE family, LFE status (+-10 min), patch LFEs (families <= 20 km, SVI P >= 0.75 + WA-NW
selected) vs 1,000 random same-UTC-day times, PNSN tremor <= 100 km +-10 min, storm percentile (PGC Z 5-10 s), production
template detection (bank, MF >= 5.5 MAD within +-30 s).
-> tables/nv_new_candidates.csv, nv_new_controls.csv, figures/nv_wave_<time>.png (real rows; raw 0.02-0.05 Hz Z)"""
import os, sys, json
os.environ['OMP_NUM_THREADS'] = '1'
S13 = '/home/jovyan/tremorferometry/svi_science/results/s13_vlfe'
NV = f'{S13}/newvlfe'; CAT = f'{S13}/catalogue'
for p in (f'{S13}/scripts', f'{S13}/hg_repro/scripts', f'{S13}/catalogue/scripts'):
    sys.path.insert(0, p)
import numpy as np, pandas as pd
import cat_common as K

NEAR = ['CN.PGC', 'CN.VGZ', 'CN.SNB', 'CN.NLLB', 'CN.PFB', 'CN.YOUB', 'UW.LRIV', 'CN.LZB']
REG = ['UW.GNW', 'UW.DOSE', 'UW.SP2', 'TA.D03D', 'TA.B05D', 'UW.WISH', 'US.NLWA']
FAR = ['IU.COR', 'US.NEW', 'US.HAWA', 'US.MSO', 'US.WVOR']
MAXFIG = 25


def storm(year):
    f = f'{NV}/data/storm_{year}.csv'
    if os.path.isfile(f):
        return pd.read_csv(f)
    from scipy.signal import butter, sosfiltfilt
    sos = butter(2, (0.1, 0.2), btype='band', fs=1.0, output='sos')
    rows = []
    for k in range(366 if year % 4 == 0 else 365):
        d = pd.Timestamp(f'{year}-01-01') + pd.Timedelta(days=k)
        fn = f'{S13}/data/disp1sps_v3/CN.PGC/{d.year}/{d.dayofyear:03d}.npz'
        v = np.nan
        if os.path.isfile(fn):
            x = np.load(fn)['d'][2].astype(float)
            if np.isfinite(x).mean() > 0.8:
                x = np.where(np.isfinite(x), x, 0.0); x = x - x.mean()
                y = sosfiltfilt(sos, x) * 1e9
                v = float(np.median([np.sqrt(np.mean(y[i:i + 600] ** 2)) for i in range(0, 86400, 600)]))
        rows.append(dict(day=str(d.date()), pgc_5_10s_rms_nm=v))
    S = pd.DataFrame(rows); S['pct'] = S.pgc_5_10s_rms_nm.rank(pct=True)
    S.to_csv(f, index=False)
    return S


def lfe_tables():
    import hg_r3 as X
    import s13_common as C
    LF = X.lfe_lists()
    R = C.svi_rows(); R = R[(R.series == 'data') & R.counted_once & (R.P >= 0.75)][['t', 'fam']]
    c = pd.read_parquet(f'{K.ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet', columns=['time', 'region', 'selected_v2', 'family_uid'])
    c = c[c.selected_v2 & (c.region == 'NW')].copy()
    c['t'] = c.time.astype('datetime64[ns]').astype('int64').values / 1e9
    c['fam'] = c.family_uid.str.split(':').str[-1]
    F = K.families()
    P = pd.read_csv(f'{K.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv', usecols=['time', 'lat', 'lon'])
    P['t'] = pd.to_datetime(P.time).astype('datetime64[ns]').astype('int64') / 1e9
    return X, LF, R, c, F, P


def wave_fig(e, fn, cc):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    y = int(e.year); SD = f'{S13}/data/scan_{y}'
    y0 = (pd.Timestamp(f'{y}-01-01') - pd.Timestamp('1970-01-01')).total_seconds()
    ts = int(round(e.epoch - y0)); w0, w1 = -120, 300
    tr = []
    for s in NEAR + REG + FAR:
        f = f'{SD}/bp_{y}v2/{s}.npy' if os.path.isfile(f'{SD}/bp_{y}v2/{s}.npy') else f'{SD}/bp2_ref/{s}.npy'
        if not os.path.isfile(f) or s not in cc.index:
            continue
        X = np.load(f, mmap_mode='r')
        if ts + w1 > X.shape[1] or ts + w0 < 0:
            continue
        x = np.asarray(X[2, ts + w0:ts + w1], float) * 1e9
        if np.isfinite(x).sum() < 0.8 * len(x):
            continue
        dist = float(K.dkm(cc.loc[s, 'lat'], cc.loc[s, 'lon'], e.lat, e.lon))
        tr.append((dist, s, x))
    tr.sort()
    if not tr:
        return False
    pks = {s: np.nanmax(np.abs(x[120:300])) for _, s, x in tr if s not in FAR and s != 'CN.LZB'}   # LZB: never used by the scan
    medpk = np.nanmedian(list(pks.values())) if pks else 1.0
    amp = np.nanmax([v for v in pks.values() if v <= 10 * medpk] or [1.0])        # stations > 10x the median peak are clipped (e.g. broken LRIV 2018)
    fig, ax = plt.subplots(figsize=(9, 0.42 * len(tr) + 2))
    t = np.arange(w0, w1)
    for k, (dist, s, x) in enumerate(tr):
        col = 'tab:red' if dist <= 40 else ('tab:blue' if dist <= 150 else 'k')
        ax.plot(t, k + 0.45 * np.clip(x / amp, -2.5, 2.5), color=col if s != 'CN.LZB' else '0.7', lw=0.8)
        pk = np.nanmax(np.abs(x[120:300]))
        ax.text(w1 + 3, k, f'{s.split(".")[1]} {dist:.0f} km  {pk:.0f} nm' + (' (not used)' if s == 'CN.LZB' else ''), fontsize=7, va='center')
    ax.axvline(0, color='0.5', lw=0.5); ax.axvspan(-60, 120, color='0.93', zorder=0)
    ax.set_yticks([]); ax.set_xlim(w0, w1 + 90); ax.set_xlabel('s after scan time (fit window shaded)')
    ax.set_title(f"{e.time_pt}  {e.lat:.2f} N {e.lon:.2f} W  z {e.z_best:.1f}  z_loso {e.z_loso:.1f}  ratio_far {e.ratio_far:.2f}\n"
                 f"vr_thrust {e.vr_thrust:.2f}  {e.lfe_status}  patch LFEs {e.patch_lfe_10m} (p {e.patch_p:.3f})  storm pct {e.storm_pct:.2f}\n"
                 f"0.02–0.05 Hz Z displacement, common scale = {amp:.0f} nm (max near/mid); red ≤ 40 km, blue ≤ 150 km, black far", fontsize=8)
    fig.tight_layout(); fig.savefig(fn, dpi=110); plt.close(fig)
    return True


def hf_ratio(sta, ep):
    """2-8 Hz 1-s RMS envelope (scan_hfenv recipe, raw counts), max in [t-60, t+60] / day median; cached in newvlfe/data/hfenv."""
    import s13_common as C
    day = pd.Timestamp(ep, unit='s').normalize()
    f = f'{NV}/data/hfenv/{sta}/{day.year}/{day.dayofyear:03d}.npy'
    if not os.path.isfile(f):
        from obspy import read, UTCDateTime
        fs = C.day_files(sta, day.strftime('%Y-%m-%d'))
        E = np.full(86400, np.nan, np.float32)
        if fs:
            t0 = UTCDateTime(str(day.date()))
            st = read(fs[0], starttime=t0 - 60, endtime=t0 + 86460)
            st = st.select(channel='BHZ') or st.select(channel='HHZ')
            for tr in st:
                sr = tr.stats.sampling_rate
                if sr > 1.5 and abs(sr - round(sr)) < 1e-3 * sr:
                    tr.stats.sampling_rate = float(round(sr))
            st.merge(method=1, fill_value=None); st = st.split()
            for tr in st:
                if tr.stats.npts * tr.stats.delta < 120:
                    continue
                tr.data = tr.data.astype(np.float64)
                tr.detrend('demean'); tr.taper(max_percentage=0.5, max_length=5.0)
                tr.filter('bandpass', freqmin=2.0, freqmax=8.0, corners=4, zerophase=True)
                sr = int(round(tr.stats.sampling_rate))
                i0 = int(np.ceil(tr.stats.starttime - t0)); k0 = int(round((t0 + i0 - tr.stats.starttime) * sr))
                x = tr.data[k0:]; n = len(x) // sr
                if n <= 0:
                    continue
                r = np.sqrt(np.mean(x[:n * sr].reshape(n, sr) ** 2, axis=1))
                a, b = max(i0, 0), min(i0 + n, 86400)
                if b > a:
                    E[a:b] = r[a - i0:b - i0]
        os.makedirs(os.path.dirname(f), exist_ok=True); np.save(f, E)
    E = np.load(f)
    u = int(ep - day.timestamp())
    seg = E[max(u - 60, 0):u + 61]
    if not np.isfinite(seg).any() or not np.isfinite(E).any():
        return np.nan
    return float(np.nanmax(seg) / np.nanmedian(E))


def mf_detect(year):
    d = {2013: 'prodq2013', 2018: 'prodr2018'}[year]
    R = pd.read_parquet(f'{CAT}/data/{d}/rows_scored_v2.parquet', columns=['template', 'kind', 'source', 'epoch', 'zmad', 'tier'])
    return R[(R.source == 'real') & ~R.template.astype(str).str.startswith('N')]


def main():
    D = pd.read_parquet(f'{NV}/tables/nv_rows_scored.parquet')
    V = json.load(open(f'{NV}/tables/nv_verdict.json'))
    sel = []
    for per in ('ETS', 'between'):
        for tz in V[f'new_{per}']['passing_tiers']:
            tier, z = tz.split(' z>=')
            sel.append((per, tier, float(z)))
    if sel:
        m = np.zeros(len(D), bool)
        for per, tier, z in sel:
            m |= (D.period == per) & D[tier] & (D.z_best >= z)
        how = 'passing tiers: ' + '; '.join(f'{p} {t} z>={z:g}' for p, t, z in sel)
    else:
        m = D.B2.values & (D.period.isin(['ETS', 'between', 'other'])).values
        how = 'no tier passes: all B2 rows (z >= 5)'
    X_ = D[m & (D['loc'] == 'new')].copy()
    if os.environ.get('NV_TEST'):
        X_ = X_.head(int(os.environ['NV_TEST']))
    X, LF, R, WA, F, P = lfe_tables()
    rng = np.random.default_rng(20261016)
    import scan_refs as SR
    CC = SR.all_coords()
    rows = []
    stm = {y: storm(y) for y in (2013, 2018)}
    MF = {y: mf_detect(y) for y in (2013, 2018)}
    for e in X_.itertuples():
        ep = float(e.epoch)
        o = dict(idx=e.Index)
        n_svi = int(X.cnt(LF['svi075'], [ep], 600)[0]); n_wa = int(X.cnt(LF['wa_nw'], [ep], 600)[0])
        n_c = int(X.cnt(LF['svi_c010_075'], [ep], 600)[0])
        dP = K.dkm(P.lat.values, P.lon.values, e.lat, e.lon)
        n_tr = int(((dP <= 100) & (np.abs(P.t.values - ep) <= 600)).sum())
        o.update(n_lfe_svi075_10m=n_svi, n_lfe_wanw_10m=n_wa, n_svi_cand_10m=n_c, n_pnsn_tremor_100km_10m=n_tr,
                 lfe_status='with LFEs' if n_svi + n_wa else ('candidates/tremor only' if n_c + n_tr else 'quiet'))
        Fn = F[K.dkm(F.lat.values, F.lon.values, e.lat, e.lon) <= 20]
        sv = Fn[Fn.region == 'S_Vancouver_Is'].family_id; sv = set(pd.to_numeric(sv, errors='coerce').dropna().astype(int))
        wa = set(Fn[Fn.region == 'Royer_NW_WA'].family_id.astype(str))
        tt = np.sort(np.r_[R[R.fam.isin(sv)].t.values, WA[WA.fam.isin(wa)].t.values])
        obs = int(X.cnt(tt, [ep], 600)[0])
        day0 = np.floor(ep / 86400) * 86400
        null = X.cnt(tt, day0 + rng.uniform(0, 86400, 1000), 600)
        o.update(n_patch_families=len(sv) + len(wa), patch_lfe_10m=obs, patch_null_mean=float(np.mean(null)),
                 patch_p=float(np.mean(null >= obs)) if obs > 0 else 1.0)
        day = str(pd.Timestamp(ep, unit='s').date())
        S = stm[int(e.year)]; q = S[S.day == day]
        o['storm_pct'] = float(q.pct.iloc[0]) if len(q) else np.nan
        hf = []
        for st_ in sorted(NEAR[:7], key=lambda s_: float(K.dkm(CC.loc[s_, 'lat'], CC.loc[s_, 'lon'], e.lat, e.lon))):
            try:
                v_ = hf_ratio(st_, ep)
            except Exception:  # noqa
                v_ = np.nan
            if np.isfinite(v_):
                hf.append((st_, v_))
            if len(hf) == 2:
                break
        o.update(hf_sta=','.join(h[0] for h in hf), hf1_ratio=hf[0][1] if hf else np.nan, hf2_ratio=hf[1][1] if len(hf) > 1 else np.nan,
                 hf_class='HF-loud' if len(hf) == 2 and min(h[1] for h in hf) >= 10 else ('HF-quiet' if len(hf) == 2 else 'HF unknown'))
        M = MF[int(e.year)]
        h = M[np.abs(M.epoch - ep) <= 30]
        o.update(mf_best_z=float(h.zmad.max()) if len(h) else np.nan, mf_templates=','.join(sorted(set(h.template.astype(str)))) if len(h) else '',
                 mf_tier=bool(h.tier.any()) if len(h) else False)
        rows.append(o)
    A = X_.join(pd.DataFrame(rows).set_index('idx')) if rows else X_
    cols = ['time_pt', 'kind', 'run', 'variant', 'year', 'period', 'lat', 'lon', 'dep', 'z_best', 'z_loso', 'ratio_far', 'n_far', 'vr_thrust',
            'share_excess', 'pw2_n', 'pw2_vr', 'hf_class', 'hf1_ratio', 'hf2_ratio', 'hf_sta', 'B0', 'B1', 'B2', 'dist_tpl_km', 'nearest_tpl', 'nearest_family', 'dist_family_km',
            'lfe_status', 'n_lfe_svi075_10m', 'n_lfe_wanw_10m', 'n_svi_cand_10m', 'n_pnsn_tremor_100km_10m', 'n_patch_families',
            'patch_lfe_10m', 'patch_null_mean', 'patch_p', 'storm_pct', 'mf_best_z', 'mf_templates', 'mf_tier', 'epoch']
    cols = [c for c in cols if c in A.columns]
    real = A[A.kind == 'real'].sort_values('epoch'); ctrl = A[A.kind == 'control'].sort_values('epoch')
    real[cols].to_csv(f'{NV}/tables/nv_new_candidates.csv', index=False)
    ctrl[cols].to_csv(f'{NV}/tables/nv_new_controls.csv', index=False)
    print(how)
    pd.set_option('display.width', 250)
    show = ['time_pt', 'period', 'lat', 'lon', 'z_best', 'z_loso', 'ratio_far', 'vr_thrust', 'hf_class', 'nearest_family', 'lfe_status',
            'patch_lfe_10m', 'patch_p', 'storm_pct', 'mf_best_z', 'mf_tier']
    print(real[[c for c in show if c in real]].round(2).to_string())
    for k in ('real', 'control'):
        x = A[A.kind == k]
        if len(x):
            print(k, len(x), 'with LFEs', float((x.lfe_status == 'with LFEs').mean()), 'patch p<0.01', float((x.patch_p < 0.01).mean()),
                  'storm>=0.85', float((x.storm_pct >= 0.85).mean()), 'HF-loud', float((x.hf_class == 'HF-loud').mean()))
    cc = CC
    for e in real.sort_values('z_best', ascending=False).head(MAXFIG).itertuples():
        tag = e.time_pt[:19].replace(' ', '_').replace(':', '')
        wave_fig(e, f'{NV}/figures/nv_wave_{tag}.png', cc)
    json.dump(dict(selection=how, n_real=len(real), n_control=len(ctrl)), open(f'{NV}/tables/nv_new_selection.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
