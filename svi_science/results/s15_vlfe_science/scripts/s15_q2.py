"""S15 Q2 — tidal modulation of the 540 ETS VLFEs vs matched LFE draws (PLAN §2 Q2, §3.10–14)."""
import json
import numpy as np, pandas as pd
import s15_common as C
import s15_q1 as Q
import s5_lib as SL

COMPS = ['tau_updip_Pa', 'tau_Pa', 'coulomb_0p4', 'sigma_n_Pa']
CONST = ['m2', 'k1']
NDRAW = 1000
rng = np.random.default_rng(2015)
T2014 = (pd.Timestamp('2014-01-01') - pd.Timestamp('1970-01-01')).total_seconds()


def load_stress():
    Sb = SL.Stress(pd.read_parquet(f'{C.RES}/s5/tidal_stress_bostock_hourly.parquet', columns=['fid', 'time'] + COMPS), COMPS)
    Ss = SL.Stress(pd.read_parquet(f'{C.RES}/s5_setup/tidal_stress_svi_hourly.parquet', columns=['fid', 'time'] + COMPS), COMPS)
    return Sb, Ss


def quantities(t, fid, Sb, Ss):
    """t: UTC epoch s (any shape), fid: SVI family ids (same shape). -> dict of arrays (same shape):
    th_<comp>_<const> (radians, 0 = peak positive stress of that constituent) and z_<comp> (stress / family std)."""
    shp = np.shape(t); t = np.asarray(t, float).ravel(); fid = np.asarray(fid).ravel().astype(int)
    out = {}
    tns = (t * 1e9).astype(np.int64)
    th = SL.astro_th(tns.astype('datetime64[ns]'), CONST)              # (n, 2)
    for c in COMPS:
        for j, k in enumerate(CONST):
            out[f'th_{c}_{k}'] = np.full(len(t), np.nan)
        out[f'z_{c}'] = np.full(len(t), np.nan)
    for S, m in ((Sb, t < T2014), (Ss, t >= T2014)):
        if not m.any():
            continue
        fi = np.array([S.fidx[f] for f in fid[m]])
        for c in COMPS:
            for j, k in enumerate(CONST):
                out[f'th_{c}_{k}'][m] = th[m, j] + np.angle(S.hc[c][k])[fi]
            out[f'z_{c}'][m] = S.value(c, fi, tns[m]) / S.std[c][fi]
    return {k: v.reshape(shp) for k, v in out.items()}


def circ(th, axis=0):
    z = np.nanmean(np.exp(1j * th), axis=axis)
    return z


def summ(qv, qc, qu, rows, label, comp, k):
    """observed set vs matched-LFE control draws vs random-time draws, rows = boolean mask over the 540."""
    th = qv[f'th_{comp}_{k}'][rows]; n = int(np.isfinite(th).sum())
    z = circ(th); R = abs(z)
    p, _ = SL.rayleigh_p(np.cos(th[np.isfinite(th)]), np.sin(th[np.isfinite(th)]))
    zc = circ(qc[f'th_{comp}_{k}'][rows], 0); zu = circ(qu[f'th_{comp}_{k}'][rows], 0)    # (NDRAW,)
    ac, au = 2 * abs(zc), 2 * abs(zu)
    zcm = zc.mean()
    o = dict(set=label, comp=comp, const=k, n=n, R=R, rayleigh_z=n * R * R, rayleigh_p=p, a=2 * R, theta0_deg=np.degrees(np.angle(z)),
             C=z.real, a_rand_p99=float(np.percentile(au, 99)), a_rand_mean=float(au.mean()), pct_vs_rand=Q.pct(2 * R, au),
             a_lfe_median=float(np.median(ac)), a_lfe_p01=float(np.percentile(ac, 1)), a_lfe_p99=float(np.percentile(ac, 99)),
             theta0_lfe_deg=float(np.degrees(np.angle(zcm))), a_lfe_of_mean=float(2 * abs(zcm)), pct_vs_lfe=Q.pct(2 * R, ac),
             C_lfe_median=float(np.median(zc.real)), C_lfe_p01=float(np.percentile(zc.real, 1)), C_lfe_p99=float(np.percentile(zc.real, 99)),
             pct_C_vs_lfe=Q.pct(z.real, zc.real), sd_a_lfe=float(ac.std()))
    o['a_minus_lfe'] = o['a'] - o['a_lfe_median']
    o['zstress_mean'] = float(np.nanmean(qv[f'z_{comp}'][rows]))
    zs = np.nanmean(qc[f'z_{comp}'][rows], 0)
    o['zstress_lfe_median'] = float(np.median(zs)); o['pct_zstress_vs_lfe'] = Q.pct(o['zstress_mean'], zs)
    o['zstress_rand_mean'] = float(np.nanmean(qu[f'z_{comp}'][rows]))
    detected = (p < 0.01) and (o['pct_vs_rand'] > 99)
    if detected:
        rel = 'stronger' if o['pct_vs_lfe'] > 99 else ('weaker' if o['pct_vs_lfe'] < 1 else 'same')
    else:
        rel = 'none' if o['pct_vs_lfe'] < 1 else 'not resolved at this N'
    o['detected'] = detected; o['verdict'] = rel
    return o


def main():
    V, F, L, D, pools, pools_loc, finfo = Q.build()
    Fi = F.set_index('key')
    S115 = F[F.src == 'svi']
    sx, sy, sf = S115.x.values, S115.y.values, S115.fid.values
    with_lfe = V.n_lfe_10m.values >= 1
    ax = np.where(with_lfe, V.anchor_x.astype(float), V.vx.astype(float)).astype(float); ay = np.where(with_lfe, V.anchor_y.astype(float), V.vy.astype(float)).astype(float)
    V['stress_fid'] = sf[np.argmin(np.hypot(ax[:, None] - sx[None], ay[:, None] - sy[None]), 1)]
    V['stress_fid_dist_km'] = np.min(np.hypot(ax[:, None] - sx[None], ay[:, None] - sy[None]), 1)
    n = len(V)
    Ct = np.full((n, NDRAW), np.nan); Cf = np.zeros((n, NDRAW), int); widen = np.full(n, np.nan)
    for i in range(n):
        p = pools[i] if with_lfe[i] else pools_loc[i]
        if p is None:
            continue
        t0 = V.epoch.values[i]
        for w in (1.0, 2.0, None):
            m = np.abs(p[0] - t0) <= w * C.DAY if w else np.ones(len(p[0]), bool)
            if m.any():
                break
        widen[i] = w if w else 99
        j = rng.choice(np.where(m)[0], NDRAW)
        Ct[i] = p[0][j]; Cf[i] = Fi.svi_fid.reindex(p[1][j]).values.astype(int)
    okc = np.isfinite(Ct[:, 0])
    # rows without any pool: fall back to random times only (control NaN)
    Ut = V.epoch.values[:, None] + rng.uniform(-C.DAY, C.DAY, (n, NDRAW))
    Uf = np.repeat(V.stress_fid.values[:, None], NDRAW, 1)
    Sb, Ss = load_stress()
    qv = quantities(V.epoch.values, V.stress_fid.values, Sb, Ss)
    Cf2 = np.where(okc[:, None], Cf, Uf)
    qc = quantities(np.where(okc[:, None], Ct, np.nan), Cf2, Sb, Ss)
    qu = quantities(Ut, Uf, Sb, Ss)
    # declustering: first VLFE of each 1-h cluster
    gap = np.r_[np.inf, np.diff(V.epoch.values)]
    first_in_cluster = gap > 3600
    V['matched_window_d'] = widen
    V['has_lfe_control'] = okc
    sets = {'all_540': np.ones(n, bool), 'with_LFE_431': with_lfe, 'no_LFE_109': ~with_lfe,
            'no_2016_17': ~V.flag_2016_17.values, 'declustered_1h': first_in_cluster,
            'era_2010_13': V.year.values <= 2013, 'era_2014_26': V.year.values >= 2014,
            'big9_episodes': V.episode.isin(C.BIG9).values}
    out = []
    for lab, rows in sets.items():
        rows = rows & okc if lab != 'all_540_incl_no_control' else rows
        for comp in COMPS:
            for k in CONST:
                out.append(summ(qv, qc, qu, rows, lab, comp, k))
    T = pd.DataFrame(out)
    T.to_csv(f'{C.OUT_T}/q2_summary.csv', index=False)
    # per-event table
    for c in COMPS:
        for k in CONST:
            V[f'th_{c}_{k}_deg'] = np.degrees(np.angle(np.exp(1j * qv[f'th_{c}_{k}'])))
        V[f'z_{c}'] = qv[f'z_{c}']
    V['hour_PT'] = pd.to_datetime(V.epoch, unit='s', utc=True).dt.tz_convert('America/Los_Angeles').dt.hour
    V['first_in_1h_cluster'] = first_in_cluster
    keep = ['time_pt', 'epoch', 'run', 'episode', 'n_lfe_10m', 'flag_2016_17', 'stress_fid', 'stress_fid_dist_km', 'matched_window_d',
            'has_lfe_control', 'first_in_1h_cluster', 'hour_PT'] + [c for c in V.columns if c.startswith('th_') or c.startswith('z_')]
    V[keep].to_csv(f'{C.OUT_T}/q2_events.csv', index=False)
    np.savez_compressed(f'{C.OUT_T}/q2_phases.npz', vlfe=qv['th_tau_updip_Pa_m2'], lfe=qc['th_tau_updip_Pa_m2'].astype(np.float32),
                        rand=qu['th_tau_updip_Pa_m2'].astype(np.float32), vlfe_k1=qv['th_tau_updip_Pa_k1'],
                        lfe_k1=qc['th_tau_updip_Pa_k1'].astype(np.float32), rand_k1=qu['th_tau_updip_Pa_k1'].astype(np.float32),
                        okc=okc, with_lfe=with_lfe, flag=V.flag_2016_17.values, zv=qv['z_tau_updip_Pa'],
                        zc=qc['z_tau_updip_Pa'].astype(np.float32), zu=qu['z_tau_updip_Pa'].astype(np.float32))
    # time-of-day diagnostic (K1 confound): solar-day phase Rayleigh for VLFEs and matched LFEs
    sol = lambda t: 2 * np.pi * ((np.asarray(t, float) % C.DAY) / C.DAY)
    zt = circ(sol(V.epoch.values[okc])); ztc = circ(sol(Ct[okc]), 0)
    tod = dict(vlfe_R_solar=float(abs(zt)), vlfe_rayleigh_p=SL.rayleigh_p(np.cos(sol(V.epoch.values[okc])), np.sin(sol(V.epoch.values[okc])))[0],
               lfe_R_solar_median=float(np.median(abs(ztc))), lfe_R_solar_p99=float(np.percentile(abs(ztc), 99)),
               vlfe_peak_hour_utc=float((np.angle(zt) % (2 * np.pi)) / (2 * np.pi) * 24))
    info = dict(n_no_control=int((~okc).sum()), matched_window_counts=pd.Series(widen).value_counts(dropna=False).to_dict(),
                stress_fid_dist_km_median=float(V.stress_fid_dist_km.median()), stress_fid_dist_km_max=float(V.stress_fid_dist_km.max()),
                n_declustered=int(first_in_cluster.sum()), time_of_day=tod)
    json.dump(info, open(f'{C.OUT_T}/q2_info.json', 'w'), indent=1, default=str)
    pd.set_option('display.width', 300)
    cols = ['set', 'comp', 'const', 'n', 'a', 'theta0_deg', 'rayleigh_z', 'rayleigh_p', 'pct_vs_rand', 'a_lfe_median', 'a_lfe_p01', 'a_lfe_p99',
            'theta0_lfe_deg', 'pct_vs_lfe', 'C', 'C_lfe_median', 'zstress_mean', 'zstress_lfe_median', 'verdict']
    print(T[cols].round(3).to_string())
    print(json.dumps(info, indent=1, default=str))


if __name__ == '__main__':
    main()
