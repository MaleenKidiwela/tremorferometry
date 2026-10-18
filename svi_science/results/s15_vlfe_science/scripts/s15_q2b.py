"""S15 Q2 robustness columns (PLAN Amendment 2): cluster bootstrap, 6-h declustering, detection vs tidal phase,
per-episode modulation, stress-level curve."""
import json
import numpy as np, pandas as pd
import s15_common as C
import s15_q2 as Q2
import cat_common as K

rng = np.random.default_rng(3015)
SRC = {'C0': (48.141757352279875, -123.52741640459564), 'C1': (48.29464879935246, -123.61541420640553)}


def a_of(th):
    th = th[np.isfinite(th)]
    z = np.mean(np.exp(1j * th)); return 2 * abs(z), np.degrees(np.angle(z))


def main():
    Z = np.load(f'{C.OUT_T}/q2_phases.npz')
    Ev = pd.read_csv(f'{C.OUT_T}/q2_events.csv')
    th = Z['vlfe']; okc = Z['okc']; thc = Z['lfe']; thu = Z['rand']
    t = Ev.epoch.values
    res = {}
    # (i) cluster bootstrap
    for gap_h in (1, 6):
        cid = np.cumsum(np.r_[True, np.diff(t) > gap_h * 3600])
        m = okc
        ids = np.unique(cid[m])
        A = []
        for _ in range(2000):
            pick = rng.choice(ids, len(ids))
            sel = np.concatenate([np.where((cid == p) & m)[0] for p in pick])
            A.append(a_of(th[sel])[0])
        A = np.array(A)
        first = np.r_[True, np.diff(t) > gap_h * 3600] & m
        a_dec, th_dec = a_of(th[first])
        ac_dec = np.array([a_of(thc[first, k])[0] for k in range(thc.shape[1])])
        res[f'cluster_{gap_h}h'] = dict(n_clusters=int(len(ids)), a_all=a_of(th[m])[0], a_boot_p025=float(np.percentile(A, 2.5)),
                                        a_boot_p975=float(np.percentile(A, 97.5)), a_boot_p01=float(np.percentile(A, 1)),
                                        lfe_a_median=float(np.median([a_of(thc[m, k])[0] for k in range(200)])),
                                        lfe_a_p99=float(np.percentile([a_of(thc[m, k])[0] for k in range(1000)], 99)),
                                        declustered_n=int(first.sum()), declustered_a=a_dec, declustered_theta0=th_dec,
                                        declustered_lfe_a_median=float(np.median(ac_dec)), declustered_lfe_a_p99=float(np.percentile(ac_dec, 99)),
                                        declustered_pct_vs_lfe=float(np.mean(ac_dec < a_dec) * 100))
    # (iii) per-episode
    rows = []
    for e in C.BIG9 + ['all']:
        m = okc & ((Ev.episode == e).values if e != 'all' else True)
        if m.sum() < 3:
            continue
        a, t0 = a_of(th[m]); ac = np.array([a_of(thc[m, k])[0] for k in range(thc.shape[1])])
        zc = np.nanmean(np.exp(1j * thc[m]), 0).mean()
        rows.append(dict(episode=e, n=int(m.sum()), a_vlfe=a, theta0_vlfe=t0, a_lfe_median=float(np.median(ac)),
                         a_lfe_p01=float(np.percentile(ac, 1)), a_lfe_p99=float(np.percentile(ac, 99)),
                         theta0_lfe=float(np.degrees(np.angle(zc))), pct_vs_lfe=float(np.mean(ac < a) * 100)))
    PE = pd.DataFrame(rows); PE.to_csv(f'{C.OUT_T}/q2_per_episode.csv', index=False)
    # (iv) stress-level curve: deciles of z_tau_updip at random times (pooled), rate ratio VLFE / random and LFE / random
    zv, zc, zu = Z['zv'], Z['zc'], Z['zu']
    edges = np.nanpercentile(zu[okc].ravel(), np.linspace(0, 100, 11)); edges[0] -= 1e9; edges[-1] += 1e9
    hv = np.histogram(zv[okc], edges)[0] / okc.sum()
    hc = np.histogram(zc[okc].ravel(), edges)[0] / np.isfinite(zc[okc]).sum()
    hu = np.histogram(zu[okc].ravel(), edges)[0] / np.isfinite(zu[okc]).sum()
    SC = pd.DataFrame(dict(decile=np.arange(1, 11), z_lo=edges[:-1], z_hi=edges[1:], vlfe_over_random=hv / hu, lfe_over_random=hc / hu,
                           n_vlfe=np.histogram(zv[okc], edges)[0]))
    SC.to_csv(f'{C.OUT_T}/q2_stress_deciles.csv', index=False)
    # (ii) detection vs tidal phase (injections)
    F, _ = C.families()
    S115 = F[F.src == 'svi']
    fid = {}
    for s, (la, lo) in SRC.items():
        x, y = C.xy(la, lo); fid[s] = int(S115.fid.values[np.argmin(np.hypot(S115.x - x, S115.y - y))])
    rows = []
    for r in C.RUNS:
        J = pd.read_parquet(f'{C.CAT}/data/{r}/inj_list.parquet').merge(pd.read_parquet(f'{C.CAT}/data/{r}/inj_v2.parquet'), on='inj_id')
        J = K.rescore_inj(J)
        J = J[J.source.isin(['C0', 'C1'])].copy()
        J['rec'] = ((J.mf_best_z >= 5.5) & (J.mf_second_z >= 5.5) & J.battery_e2.fillna(False).astype(bool)).astype(int)
        J['mf_only'] = ((J.mf_best_z >= 5.5) & (J.mf_second_z >= 5.5)).astype(int)
        J['run'] = r
        rows.append(J[['run', 'source', 'Mw', 'period', 'epoch', 'rec', 'mf_only']])
    J = pd.concat(rows, ignore_index=True)
    Sb, Ss = Q2.load_stress()
    q = Q2.quantities(J.epoch.values.astype(float), J.source.map(fid).values, Sb, Ss)
    J['th_m2'] = q['th_tau_updip_Pa_m2']; J['th_k1'] = q['th_tau_updip_Pa_k1']; J['z_tau_updip'] = q['z_tau_updip_Pa']
    J.to_csv(f'{C.OUT_T}/q2_injections_phase.csv', index=False)
    det = []
    for lab, sub in [('Mw3.6', J[J.Mw == 3.6]), ('Mw3.8', J[J.Mw == 3.8]), ('Mw3.4-4.0', J), ('Mw3.6_ETS', J[(J.Mw == 3.6) & (J.period == 'ETS')])]:
        for col in ('rec', 'mf_only'):
            for k in ('m2', 'k1'):
                ph = sub[f'th_{k}'].values; y = sub[col].values.astype(float)
                X = np.c_[np.ones(len(y)), np.cos(ph), np.sin(ph)]
                c, *_ = np.linalg.lstsq(X, y, rcond=None)
                b = np.hypot(c[1], c[2]) / c[0]
                # bootstrap sd of b
                bs = []
                for _ in range(300):
                    i = rng.integers(0, len(y), len(y)); cc, *_ = np.linalg.lstsq(X[i], y[i], rcond=None); bs.append(np.hypot(cc[1], cc[2]) / cc[0])
                bins = np.histogram(np.angle(np.exp(1j * ph)), np.linspace(-np.pi, np.pi, 9))[0]
                rec_bins = [float(y[(np.angle(np.exp(1j * ph)) >= a0) & (np.angle(np.exp(1j * ph)) < a0 + np.pi / 4)].mean())
                            for a0 in np.linspace(-np.pi, np.pi, 9)[:-1]]
                det.append(dict(set=lab, outcome=col, const=k, n=len(y), mean=float(y.mean()), modulation_b=float(b),
                                b_boot_p95=float(np.percentile(bs, 95)), phase_deg=float(np.degrees(np.arctan2(c[2], c[1]))),
                                rec_by_45deg_bin=';'.join(f'{x:.3f}' for x in rec_bins)))
    DT = pd.DataFrame(det); DT.to_csv(f'{C.OUT_T}/q2_detection_vs_tide.csv', index=False)
    json.dump(res, open(f'{C.OUT_T}/q2b_clusters.json', 'w'), indent=1, default=float)
    pd.set_option('display.width', 250)
    print(json.dumps(res, indent=1, default=float)); print(PE.round(3).to_string()); print(SC.round(3).to_string()); print(DT.round(3).to_string())


if __name__ == '__main__':
    main()
