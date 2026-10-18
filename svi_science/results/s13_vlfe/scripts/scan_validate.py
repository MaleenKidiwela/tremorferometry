#!/usr/bin/env python
"""VLFE scan pilot: controls, headline threshold, validation 6a-6e, go/no-go, figures (PLAN_SCAN.md §5-7, S1-S5).
Writes tables/scan_separation.csv, tables/scan_filter_validation.csv, tables/scan_injection.csv,
tables/scan_injection_summary.csv, tables/scan_hg_refind.csv, tables/scan_detections_2011.csv,
tables/scan_verdict.json, figures/scan_*.png."""
import os, sys, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import scan_common as S

R_ = f'{S.SD}/runs'
INJ = os.environ.get('SCAN_INJ_RUN', 'inj2')
THRS = [4, 4.5, 5, 5.5, 6, 7]
SHIFTS = [3, 5, 7, 11, 13]
CRIT = ['n_sta', 'n_Z_agree', 'max_share', 'vr_ratio', 'vr_gt_ss', 'vr_gt_steep', 'isolation']


def crit_pass(D, share_rule='max_share'):
    P = pd.DataFrame(index=D.index)
    P['n_sta'] = D.n_sta >= 3
    P['n_Z_agree'] = D.n_Z_agree >= 2
    P['max_share'] = D.max_share <= 0.5
    P['share_excess'] = D.share_excess <= 0.25
    P['vr_ratio'] = D.vr_ratio >= 0.8
    P['vr_gt_ss'] = D.vr_thrust > D.vr_strikeslip
    P['vr_gt_steep'] = D.vr_thrust > D.vr_steepdip
    P['isolation'] = D.isolation < 2
    P['not_hf_loud'] = ~D.class_eq.fillna(False).astype(bool)
    return P


def headline_ok(D, use):
    P = crit_pass(D)
    ok = np.ones(len(D), bool)
    for c in use:
        ok &= P[c].values
    return ok


def between_mask():
    ep = S.episodes_svi(); ep = ep[ep.major.astype(bool)]
    m = np.zeros(S.NDAY * 86400, bool)
    y0 = S.Y0
    for r in ep.itertuples():
        a = int((pd.Timestamp(r.start) - pd.Timedelta(days=2) - y0).total_seconds())
        b = int((pd.Timestamp(r.end) + pd.Timedelta(days=3) - y0).total_seconds())
        a, b = max(a, 0), min(b, len(m))
        if b > a:
            m[a:b] = True
    return ~m                                                 # True = between-ETS


def ets_vec():
    m = np.zeros(S.NDAY * 86400, bool)
    a = int((S.ETS[0] - S.Y0).total_seconds()); b = int((S.ETS[1] - S.Y0).total_seconds())
    m[a:b] = True
    return m


def live_years(navail, maskv, sel=None):
    ok = (navail >= 3) & ~maskv
    if sel is not None:
        ok &= sel
    return ok.sum() / (365.25 * 86400)


def main():
    V = {}
    masks = {t: S.mask_vector(t) for t in ('head', 'm45', 't58')}
    btw = between_mask(); etsv = ets_vec()
    real = pd.read_parquet(f'{R_}/real/cand_cols.parquet')
    nav_real = np.load(f'{R_}/real/navail.npy')
    inj = pd.read_parquet(f'{R_}/{INJ}/cand_cols.parquet')
    J = pd.read_parquet(f'{S.SD}/inject_2011.parquet')
    nodes = pd.read_csv(f'{C.TB}/scan_grid_nodes.csv')
    # ---------------- S5.3 filter validation on injections ----------------
    pos = inj[(inj.variant == 0) & (inj.sign > 0)].sort_values('tsec')
    t = pos.tsec.values
    mt = []
    for r in J.itertuples():
        lo, hi = np.searchsorted(t, r.tsec - 60), np.searchsorted(t, r.tsec + 60, 'right')
        if hi > lo:
            k = pos.index[lo + int(np.argmax(pos.z.values[lo:hi]))]
            mt.append((r.inj, k))
        else:
            mt.append((r.inj, -1))
    J['cand'] = [k for _, k in mt]
    rec = J[(J.cand >= 0)]
    Rc = pos.loc[rec.cand.values].reset_index().assign(inj=rec.inj.values)
    JJ = J.merge(Rc[['inj', 'z', 'm0', 'node', 'lat', 'lon'] + [c for c in Rc.columns if c in (
        'n_sta', 'n_Z_agree', 'max_share', 'share_excess', 'vr_ratio', 'vr_thrust', 'vr_strikeslip', 'vr_steepdip',
        'isolation', 'class_eq', 'hf1', 'hf2', 'in_mask_head', 'region_km', 'zZ', 'm0_20s')]].rename(
        columns=dict(lat='lat_det', lon='lon_det', m0='m0_det')), on='inj', how='left')
    P = crit_pass(JJ.set_index('inj').loc[:, [c for c in JJ.columns if c not in ('inj',)]])
    vsel = (JJ.gf == 'standard') & (JJ.mw >= 3.4) & JJ.z.notna() & (JJ.z >= 4)
    fv = []
    use = []
    for c in CRIT:
        frac = float(P[c].values[vsel.values].mean())
        row = dict(criterion=c, pass_frac_inj_mw34=frac, n=int(vsel.sum()), kept=frac >= 0.85)
        if c == 'max_share' and frac < 0.85:
            f2 = float(P['share_excess'].values[vsel.values].mean())
            row.update(fallback='share_excess<=0.25', fallback_pass_frac=f2, kept=False)
            if f2 >= 0.85:
                use.append('share_excess')
        elif frac >= 0.85:
            use.append(c)
        fv.append(row)
    use.append('not_hf_loud')
    fv.append(dict(criterion='not_hf_loud', pass_frac_inj_mw34=float(P['not_hf_loud'].values[vsel.values].mean()),
                   n=int(vsel.sum()), kept=True, fallback='kept by rule (not checkable on injections)'))
    FV = pd.DataFrame(fv); FV.to_csv(f'{C.TB}/scan_filter_validation.csv', index=False)
    V['headline_filter'] = use
    print(FV.round(3).to_string())
    # ---------------- separation curves ----------------
    def tier_counts(D, nav, tier, sel_time=None, filt=True):
        mv = masks[tier]
        d = D[~mv[D.tsec.values]]
        if sel_time is not None:
            d = d[sel_time[d.tsec.values]]
        ok = headline_ok(d, use) if filt else np.ones(len(d), bool)
        return d[ok]
    rows = []
    ctrlA = {s: pd.read_parquet(f'{R_}/ctrlA_s{s}/cand_cols.parquet') for s in SHIFTS}
    navA = {s: np.load(f'{R_}/ctrlA_s{s}/navail.npy') for s in SHIFTS}
    nd = int(real.variant.max())
    for tier in ('head', 'm45', 't58'):
        for filt in (True, False):
            ly = live_years(nav_real, masks[tier]); lyb = live_years(nav_real, masks[tier], btw)
            r0 = tier_counts(real[real.variant == 0], nav_real, tier, filt=filt)
            rb = tier_counts(real[real.variant > 0], nav_real, tier, filt=filt)
            rbb = rb[btw[rb.tsec.values]]
            for thr in THRS:
                rp = r0[(r0.sign > 0) & (r0.z >= thr)]; rn = r0[(r0.sign < 0) & (-r0.z >= thr)]
                a_rates = []
                for s in SHIFTS:
                    da = tier_counts(ctrlA[s], navA[s], tier, filt=filt)
                    a_rates.append(((da.sign > 0) & (da.z >= thr)).sum() / live_years(navA[s], masks[tier]))
                bpos = rb[(rb.sign > 0) & (rb.z >= thr)]; bposb = rbb[(rbb.sign > 0) & (rbb.z >= thr)]
                rows.append(dict(tier=tier, filtered=filt, thr=thr, live_yr=ly, live_yr_between=lyb,
                                 real_pos=len(rp), real_neg=len(rn), real_pos_per_yr=len(rp) / ly,
                                 real_neg_per_yr=len(rn) / ly,
                                 real_pos_between=int(btw[rp.tsec.values].sum()),
                                 real_pos_between_per_yr=btw[rp.tsec.values].sum() / lyb,
                                 real_pos_ets=int(etsv[rp.tsec.values].sum()), real_neg_ets=int(etsv[rn.tsec.values].sum()),
                                 ctrlA_per_yr=float(np.mean(a_rates)), ctrlA_per_yr_each=';'.join(f'{x:.2f}' for x in a_rates),
                                 ctrlB_per_yr=len(bpos) / nd / ly, ctrlB_between_per_yr=len(bposb) / nd / lyb))
    SEP = pd.DataFrame(rows); SEP.to_csv(f'{C.TB}/scan_separation.csv', index=False)
    H = SEP[(SEP.tier == 'head') & SEP.filtered]
    okthr = H[(H.ctrlA_per_yr <= 1) & (H.ctrlB_between_per_yr <= 2)]
    thr = float(okthr.thr.min()) if len(okthr) else np.nan
    V['headline_threshold'] = thr
    print(SEP[SEP.tier == 'head'].round(2).to_string())
    # ---------------- detections at the headline threshold ----------------
    thr_use = thr if np.isfinite(thr) else 7.0
    V['threshold_used_for_validation'] = thr_use
    r0 = real[real.variant == 0].copy()
    r0['headline'] = (~masks['head'][r0.tsec.values]) & headline_ok(r0, use) & (r0.z.abs() >= thr_use)
    r0['between'] = btw[r0.tsec.values]; r0['ets'] = etsv[r0.tsec.values]
    r0['zZ_agree'] = np.where(r0.z >= 1.2 * thr_use, r0.zZ >= thr_use, np.nan)
    det = r0[r0.headline]
    det.to_csv(f'{C.TB}/scan_detections_2011.csv', index=False)
    hp = det[det.sign > 0]; hn = det[det.sign < 0]
    V['n_headline_pos'] = int(len(hp)); V['n_headline_neg'] = int(len(hn))
    V['ets_pos'] = int(hp.ets.sum()); V['ets_neg'] = int(hn.ets.sum())
    V['ets_sign_accept'] = bool(V['ets_pos'] >= 10 * max(V['ets_neg'], 0) and V['ets_pos'] > 0)
    V['between_pos'] = int(hp.between.sum())
    za = hp[hp.z >= 1.2 * thr_use]
    V['zonly_agree_frac'] = float((za.zZ >= thr_use).mean()) if len(za) else np.nan
    # 6a
    live_ets = ((nav_real >= 3) & ~masks['head'] & etsv).sum() / 86400
    live_rest = ((nav_real >= 3) & ~masks['head'] & ~etsv).sum() / 86400
    dens_ets = V['ets_pos'] / live_ets if live_ets > 0 else np.nan
    dens_rest = (len(hp) - V['ets_pos']) / live_rest if live_rest > 0 else np.nan
    V['6a_ets_per_day'] = dens_ets; V['6a_rest_per_day'] = dens_rest
    V['6a_ratio'] = dens_ets / dens_rest if dens_rest > 0 else np.inf
    V['6a_pass'] = bool(V['6a_ratio'] >= 10)
    V['6a_sum_m0_ets'] = float(hp[hp.ets].m0.sum()); V['6a_frac_of_S13_7.7e15'] = V['6a_sum_m0_ets'] / 7.7e15
    # 6b
    V.update(check_6b(nodes))
    # 6c
    allc = pd.read_parquet(f'{R_}/real/cand.parquet')
    allc = allc[(allc.variant == 0) & (allc.sign > 0)]
    hg = []
    for k, ts, la, lo, dep, mw in S.HG:
        tt = (pd.Timestamp(ts) - S.Y0).total_seconds()
        m = allc[(allc.tsec - tt).abs() <= 60]
        if len(m):
            b = m.loc[m.z.idxmax()]
            dkm = float(S.dist_km(nodes.lat[b.node], nodes.lon[b.node], la, lo))
            mwo = (np.log10(b.m0) - 9.1) / 1.5
            rr = real[(real.variant == 0) & ((real.tsec - tt).abs() <= 150)]
            cl = rr.loc[rr.z.abs().idxmax()] if len(rr) else None
            hg.append(dict(template=k, time=ts, hg_mw=mw, z=b.z, dt_s=int(b.tsec - tt), node_lat=nodes.lat[b.node],
                           node_lon=nodes.lon[b.node], dist_km=dkm, mw_scan=mwo, found=bool(b.z >= thr_use and dkm <= 40),
                           cluster_sign=int(cl.sign) if cl is not None else 0, cluster_z=float(cl.z) if cl is not None else np.nan,
                           cluster_headline=bool((~masks['head'][int(cl.tsec)]) and headline_ok(cl.to_frame().T, use)[0])
                           if cl is not None and 'n_sta' in cl.index and pd.notna(cl.n_sta) else False))
        else:
            hg.append(dict(template=k, time=ts, hg_mw=mw, found=False))
    HG = pd.DataFrame(hg); HG.to_csv(f'{C.TB}/scan_hg_refind.csv', index=False)
    V['6c_found'] = int(HG.found.sum()); V['6c_pass'] = bool(V['6c_found'] >= 3)
    print(HG.round(2).to_string())
    # 6d
    JJ['recovered'] = JJ.z.notna() & (JJ.z >= thr_use) & ~JJ.in_mask_head.fillna(True).astype(bool)
    okf = np.zeros(len(JJ), bool)
    hasz = JJ.z.notna().values
    okf[hasz] = headline_ok(JJ[hasz], use)
    JJ['recovered'] &= okf
    JJ['loc_err_km'] = S.dist_km(JJ.lat_det, JJ.lon_det, JJ.lat, JJ.lon)
    JJ['m0_bias_log10'] = np.log10(JJ.m0_det / JJ.m0)
    JJ.to_csv(f'{C.TB}/scan_injection.csv', index=False)
    rows = []
    for (gf, season, tau), g in JJ.groupby(['gf', 'season', 'dur_s']):
        for mw, gg in g.groupby('mw'):
            rows.append(dict(gf=gf, season=season, tau=tau, mw=mw, n=len(gg), rec=gg.recovered.mean()))
    for (gf, season), g in JJ.groupby(['gf', 'season']):
        for mw, gg in g.groupby('mw'):
            rows.append(dict(gf=gf, season=season, tau='both', mw=mw, n=len(gg), rec=gg.recovered.mean(),
                             loc_med=gg[gg.recovered].loc_err_km.median(), loc_p95=gg[gg.recovered].loc_err_km.quantile(.95),
                             bias_med=gg[gg.recovered].m0_bias_log10.median()))
    for gf, g in JJ.groupby('gf'):
        for mw, gg in g.groupby('mw'):
            rows.append(dict(gf=gf, season='all', tau='both', mw=mw, n=len(gg), rec=gg.recovered.mean(),
                             loc_med=gg[gg.recovered].loc_err_km.median(), loc_p95=gg[gg.recovered].loc_err_km.quantile(.95),
                             bias_med=gg[gg.recovered].m0_bias_log10.median()))
    IS = pd.DataFrame(rows); IS.to_csv(f'{C.TB}/scan_injection_summary.csv', index=False)

    def mwc(season, level):
        g = IS[(IS.gf == 'standard') & (IS.season == season) & (IS.tau == 'both')].sort_values('mw')
        x, y = g.mw.values, g.rec.values
        for i in range(len(x)):
            if y[i] >= level:
                if i == 0:
                    return float(x[0])
                return float(x[i - 1] + (level - y[i - 1]) * (x[i] - x[i - 1]) / (y[i] - y[i - 1]))
        return np.nan
    for se in ('summer', 'winter', 'all'):
        V[f'6d_mw50_{se}'] = mwc(se, 0.5); V[f'6d_mw80_{se}'] = mwc(se, 0.8)
    gs = IS[(IS.gf == 'standard') & (IS.season == 'summer') & (IS.tau == 'both')].set_index('mw')
    r36 = float(gs.rec.get(3.6, np.nan)); r33 = float(np.interp(3.3, gs.index.values, gs.rec.values))
    ga = IS[(IS.gf == 'standard') & (IS.season == 'all') & (IS.tau == 'both')].set_index('mw')
    loc36 = float(ga.loc_med.get(3.6, np.nan))
    V.update({'6d_rec_mw36_summer': r36, '6d_rec_mw33_summer': r33, '6d_locmed_mw36_km': loc36,
              '6d_loc95_mw36_km': float(ga.loc_p95.get(3.6, np.nan)),
              '6d_perturbed_rec_mw36': float(JJ[JJ.gf == 'perturbed'].recovered.mean())})
    V['6d_pass'] = bool(r36 >= 0.8 and r33 >= 0.5 and loc36 <= 30)
    # 6e
    E6 = pd.read_csv(f'{C.TB}/scan_loc6e_summary.csv')
    V['6e_dist_km'] = E6.set_index('episode').dist_to_centroid_km.to_dict()
    V['6e_pass'] = bool(E6.pass_20km.all())
    # controls at headline
    hrow = H[H.thr == thr_use].iloc[0]
    V['ctrlA_per_yr'] = float(hrow.ctrlA_per_yr); V['ctrlB_between_per_yr'] = float(hrow.ctrlB_between_per_yr)
    V['ctrlB_all_per_yr'] = float(hrow.ctrlB_per_yr)
    V['controls_pass'] = bool(V['ctrlA_per_yr'] <= 1 and V['ctrlB_between_per_yr'] <= 2)
    V['GO'] = bool(V['6a_pass'] and V.get('6b_pass', False) and V['6c_pass'] and V['6d_pass'] and V['controls_pass']
                   and V['6e_pass'])
    json.dump(V, open(f'{C.TB}/scan_verdict.json', 'w'), indent=1, default=float)
    print(json.dumps(V, indent=1, default=float))
    figures(real, r0, SEP, JJ, IS, masks, nav_real, thr_use, use, HG)


def check_6b(nodes):
    """6b: mean M0 at the S13 2011-08 centroid node at the 1,445 P >= 0.75 LFE origins (S13 origins + lag 1 s), dropping
    LFEs that fail S13 window QC at a station-component carrying >= 10 % of the template weight."""
    out = {}
    m0c = np.load(f'{R_}/real/m0_centroid.npy')
    R = C.svi_rows()
    R = R[(R.series == 'data') & R.counted_once & (R.date >= '2011-08-12') & (R.date <= '2011-09-03') & (R.P >= 0.75)]
    sh = C.origin_shift_svi(sorted(R.fam.unique()))
    to = R.t.values + R.fam.map(sh).fillna(0).values + 1.0
    ts = np.round(to - S.Y0.value / 1e9).astype(int)
    out['6b_n_lfe'] = int(len(ts))
    TP = np.load(f'{S.SD}/templates_grid.npz'); stas = list(TP['stations'])
    Tc = TP['thrust'][-1] * 1e21
    cfg = pd.read_parquet(f'{R_}/real/daycfg.parquet')
    BP = {s: np.load(f'{S.SD}/bp_2011g/{s}.npy', mmap_mode='r') for s in stas}
    # QC thresholds per station-component over the episode windows (-150..+250 s)
    W = {}
    for s in stas:
        X = np.stack([np.asarray(BP[s][:, t - 150:t + 251]) for t in ts]) * 1e9           # (n, 3, 401)
        nanf = np.isnan(X).mean(axis=2)
        rms = np.sqrt(np.nanmean(X ** 2, axis=2))
        thr = 3 * np.nanmedian(rms, axis=0)
        W[s] = (nanf <= 0.05) & (rms <= thr[None, :])
    m0g, m0u, keep = [], [], []
    for j, t in enumerate(ts):
        d = t // 86400
        c = cfg[cfg.day == d]
        num = den = numu = denu = 0.0
        shares = {}
        tot = 0.0
        for i, s in enumerate(stas):
            x = np.asarray(BP[s][:, t + S.FIT0:t + S.FIT1 + 1]) * 1e9
            for k, cn in enumerate('NEZ'):
                r = c[(c.sta == s) & (c.comp == cn)].iloc[0]
                if not r.w_on or not np.isfinite(x[k]).all():
                    continue
                w = 1.0 / r.sigma_nm ** 2
                e = w * r.g ** 2 * (Tc[i, k] ** 2).sum()
                shares[(s, k)] = e; tot += e
                num += w * r.g * (Tc[i, k] * x[k]).sum(); den += e
                numu += w * (Tc[i, k] * x[k]).sum(); denu += w * (Tc[i, k] ** 2).sum()
        okq = all(W[s][j, k] for (s, k), e in shares.items() if tot > 0 and e / tot >= 0.10)
        keep.append(okq and den > 0)
        m0g.append(num / den * 1e12 if den > 0 else np.nan); m0u.append(numu / denu * 1e12 if denu > 0 else np.nan)
    keep = np.array(keep); m0g = np.array(m0g); m0u = np.array(m0u)
    gpu = m0c[ts]
    out['6b_n_kept'] = int(keep.sum())
    out['6b_mean_m0_gpu_centroid'] = float(np.nanmean(gpu[keep]))
    out['6b_mean_m0_cpu_gains'] = float(np.nanmean(m0g[keep]))
    out['6b_mean_m0_cpu_g1'] = float(np.nanmean(m0u[keep]))
    out['6b_ratio_to_5.3e12'] = out['6b_mean_m0_gpu_centroid'] / 5.3e12
    out['6b_pass'] = bool(abs(out['6b_ratio_to_5.3e12'] - 1) <= 0.2)
    out['6b_bug_flag'] = bool(out['6b_ratio_to_5.3e12'] > 2 or out['6b_ratio_to_5.3e12'] < 0.5)
    # standard error from 200 random same-day draws (scrambled control)
    rng = np.random.default_rng(5)
    dr = []
    for _ in range(200):
        tt = (ts // 86400) * 86400 + rng.integers(0, 86400, len(ts))
        dr.append(np.nanmean(m0c[tt]))
    out['6b_ctrl_mean'] = float(np.mean(dr)); out['6b_ctrl_sd'] = float(np.std(dr))
    return out


def figures(real, r0, SEP, JJ, IS, masks, nav, thr, use, HG):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    zmax = np.load(f'{R_}/real/zmax.npy')
    # 1) z(t) for 2011 with masks and ETS
    n = 600
    zr = np.nanmax(np.where(np.isfinite(zmax), zmax, -np.inf).reshape(-1, n), axis=1)
    tt = S.Y0 + pd.to_timedelta(np.arange(len(zr)) * n, unit='s')
    mk = masks['head'].reshape(-1, n).mean(axis=1) > 0.5
    fig, ax = plt.subplots(2, 1, figsize=(13, 7), sharex=False)
    for a_, (t0, t1), ttl in ((ax[0], (S.Y0, S.Y0 + pd.Timedelta(days=365)), '2011'),
                              (ax[1], (pd.Timestamp('2011-08-10'), pd.Timestamp('2011-09-06')), '2011-08 ETS')):
        sel = (tt >= t0) & (tt < t1)
        a_.plot(tt[sel], np.clip(zr[sel], 0, 40), color='0.6', lw=0.4, label='max z over grid (10-min max)')
        a_.fill_between(tt[sel], 0, 40, where=mk[sel], color='tab:orange', alpha=0.15, step='mid', label='earthquake mask (headline)')
        a_.axvspan(S.ETS[0], S.ETS[1], color='tab:blue', alpha=0.08, label='SVI ETS 2011-08-12…09-03')
        d = r0[r0.headline & (r0.sign > 0)]
        tdet = S.Y0 + pd.to_timedelta(d.tsec, unit='s'); sd = (tdet >= t0) & (tdet < t1)
        sc = a_.scatter(tdet[sd], np.clip(d.z[sd], 0, 40), c=d.cc[sd], cmap='viridis', vmin=0, vmax=1, s=18, zorder=3,
                        label=f'headline detections (z ≥ {thr:g})')
        a_.axhline(thr, color='k', lw=0.6, ls='--')
        if ttl != '2011':
            for k, ts_, la, lo, dep, mw in S.HG:
                a_.axvline(pd.Timestamp(ts_), color='r', lw=0.8, ymax=0.15)
                a_.annotate(f'H&G {k}', (pd.Timestamp(ts_), 36), color='r', fontsize=7, rotation=90)
        a_.set_ylim(0, 40); a_.set_ylabel('z (clipped at 40)'); a_.set_title(f'VLFE scan, S Vancouver Island, {ttl}', fontsize=9)
        a_.legend(fontsize=7, loc='upper left', ncol=2)
    fig.colorbar(sc, ax=ax, label='network CC')
    fig.savefig(f'{C.FIG}/scan_z_2011.png', dpi=130); plt.close(fig)
    # 2) separation curves
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for a_, filt in zip(ax, (False, True)):
        g = SEP[(SEP.tier == 'head') & (SEP.filtered == filt)]
        a_.semilogy(g.thr, g.real_pos_per_yr.clip(lower=0.05), 'o-', color='tab:red', label='real, positive')
        a_.semilogy(g.thr, g.real_neg_per_yr.clip(lower=0.05), 's-', color='tab:blue', label='real, negative (sign control)')
        a_.semilogy(g.thr, g.ctrlA_per_yr.clip(lower=0.05), '^-', color='0.3', label='(a) station-day shifted, mean of 5')
        a_.semilogy(g.thr, g.ctrlB_between_per_yr.clip(lower=0.05), 'v-', color='tab:green', label='(b) polarity scrambled, between-ETS')
        a_.semilogy(g.thr, g.real_pos_between_per_yr.clip(lower=0.05), 'o--', color='tab:purple', label='real positive, between-ETS')
        a_.axhline(1, color='k', lw=0.5, ls=':'); a_.axhline(2, color='k', lw=0.5, ls=':')
        a_.set_xlabel('z threshold'); a_.set_ylabel('clusters per live year (floor 0.05)')
        a_.set_title(('headline filter' if filt else 'raw clusters, no filter') + ' (headline mask)', fontsize=9)
        a_.legend(fontsize=7)
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan_separation.png', dpi=130); plt.close(fig)
    # 3) injection completeness
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for se, mkr in (('summer', 'o-'), ('winter', 's-')):
        for tau, col in ((1.0, 'tab:red'), (20.0, 'tab:blue')):
            g = IS[(IS.gf == 'standard') & (IS.season == se) & (IS.tau.astype(str) == str(tau))]
            ax[0].plot(g.mw, g.rec, mkr, color=col, label=f'{se}, {tau:g}-s source', alpha=0.8)
    ax[0].axhline(0.5, color='k', lw=0.5, ls=':'); ax[0].axhline(0.8, color='k', lw=0.5, ls=':')
    ax[0].set_xlabel('injected Mw'); ax[0].set_ylabel(f'recovered (headline, z ≥ {thr:g})'); ax[0].legend(fontsize=7)
    r = JJ[JJ.recovered & (JJ.gf == 'standard')]
    sc = ax[1].scatter(r.mw + np.random.default_rng(1).uniform(-0.05, 0.05, len(r)), r.loc_err_km, c=r.z.clip(0, 30), s=14, cmap='viridis')
    ax[1].set_xlabel('injected Mw'); ax[1].set_ylabel('location error (km, best node vs truth)'); fig.colorbar(sc, ax=ax[1], label='z')
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan_injection.png', dpi=130); plt.close(fig)
    # 4) detection map
    P = pd.read_csv(f'{C.TB}/scan_footprint_polygon.csv'); F = pd.read_csv(f'{C.TB}/scan_footprint_families.csv')
    cc = C.station_coords()
    fig, a_ = plt.subplots(figsize=(6.5, 6))
    a_.plot(P.lon, P.lat, 'k-', lw=0.8, label='LFE footprint (10-km buffer)')
    a_.scatter(F.lon, F.lat, s=4, c='0.6', label='SVI families')
    d = r0[r0.headline & (r0.sign > 0)]
    for lab, sel, mk_ in (('ETS', d.ets, 'o'), ('between ETS', d.between, '*'), ('other', ~(d.ets | d.between), 's')):
        if sel.any():
            sc = a_.scatter(d.lon[sel], d.lat[sel], c=d.z[sel].clip(0, 30), cmap='viridis', vmin=thr, vmax=30, s=40, marker=mk_,
                            edgecolor='k', lw=0.3, label=lab)
    for k, ts_, la, lo, dep, mw in S.HG:
        a_.plot(lo, la, 'rx', ms=8); a_.annotate(f'H&G{k}', (lo, la), color='r', fontsize=7)
    for s in S.PILOT:
        a_.plot(cc.loc[s, 'lon'], cc.loc[s, 'lat'], 'k^', ms=6); a_.annotate(s.split('.')[1], (cc.loc[s, 'lon'], cc.loc[s, 'lat']), fontsize=7)
    a_.set_xlim(-125, -122.4); a_.set_ylim(47.7, 49.4); a_.legend(fontsize=7, loc='lower left')
    a_.set_title(f'2011 headline detections (z ≥ {thr:g}); colour = z', fontsize=9)
    fig.tight_layout(); fig.savefig(f'{C.FIG}/scan_map_2011.png', dpi=130); plt.close(fig)


if __name__ == '__main__':
    main()
