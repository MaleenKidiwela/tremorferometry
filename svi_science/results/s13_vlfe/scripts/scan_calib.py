#!/usr/bin/env python
"""VLFE scan step 3 (PLAN_SCAN.md §3): grid templates, station gain/polarity calibration from the guarded S13
per-episode pooled stacks, noise table.
Writes data/scan/templates_grid.npz, tables/scan_grid_nodes.csv, tables/sc_gain.csv, tables/scan_calib_episodes.csv,
tables/scan_noise_month.csv, tables/scan_noise_ref.csv, tables/scan_calib_verdict.json."""
import os, sys, json, pickle, glob
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_green as G
import scan_common as S
from s13_stack import TREL, FITM

STAS = S.PILOT + S.NOISE_ONLY


def node_templates(args):
    i, src, coords = args
    out = {}
    for m, mech in S.MECHS.items():
        M = G.mt_ned(*mech, 1.0)
        out[m] = np.stack([S.template(M, src, coords[s]) for s in STAS])
    out['thrust20'] = np.stack([S.template(G.mt_ned(*S.MECHS['thrust'], 1.0), src, coords[s], tau=20.0) for s in STAS])
    out['elem'] = np.stack([np.stack([S.template(E, src, coords[s]) for s in STAS]) for E in S.ELEM])
    return i, out


def build_templates():
    N = S.grid()
    N.to_csv(f'{C.TB}/scan_grid_nodes.csv', index=False)
    cc = C.station_coords()
    coords = {s: (cc.loc[s, 'lat'], cc.loc[s, 'lon']) for s in STAS}
    from multiprocessing import Pool
    res = {}
    with Pool(8) as p:
        for i, o in p.imap_unordered(node_templates, [(i, (r.lat, r.lon, r.dep), coords) for i, r in N.iterrows()]):
            res[i] = o
    keys = ['thrust', 'strikeslip', 'steepdip', 'thrust20', 'elem']
    A = {k: np.stack([res[i][k] for i in range(len(N))]).astype(np.float64) for k in keys}
    np.savez_compressed(f'{S.SD}/templates_grid.npz', stations=np.array(STAS), tau=S.TAU, lat=N.lat.values,
                        lon=N.lon.values, dep=N.dep.values, **A)
    return N, A


def epoch_of(E, netsta, date):
    e = E[(E.netsta == netsta) & E.cha.str.endswith('Z')]
    d = str(date)
    for r in e.itertuples():
        if r.start[:19] <= d and (not isinstance(r.end, str) or r.end == '' or d <= r.end[:19]):
            return f'{r.cha[:2]}:{r.start[:10]}', r.start[:19], (r.end[:19] if isinstance(r.end, str) else '9999')
    return None, None, None


def pick(b):
    for pre in ('scalef', 'refit', 'scale'):
        f = f'{C.TB}/{pre}_svi_{b}_units.csv'
        if os.path.isfile(f):
            return pre
    return None


def calibrate():
    Ep = pd.read_csv(f'{C.TB}/station_epochs.csv')
    Ep['end'] = Ep['end'].fillna('')
    EPS = pd.read_csv(f'{C.TB}/svi_episodes.csv')
    EPS = EPS[EPS.built]
    # the 2011 epochs to calibrate (SNB has two in 2011)
    targets = []
    for s in STAS:
        seen = set()
        for d in ('2011-01-15', '2011-06-15'):
            lab, a, b = epoch_of(Ep, s, d + 'T00:00:00')
            if lab and lab not in seen:
                seen.add(lab); targets.append((s, lab, a, b))
    rows, erows = [], []
    acc = {}
    for e in EPS.itertuples():
        pre = pick(e.episode)
        if pre is None:
            continue
        U = pd.read_csv(f'{C.TB}/{pre}_svi_{e.episode}_units.csv')
        u = U[U.unit == 'pooled'].iloc[0]
        Sc = pd.read_csv(f'{C.TB}/{pre}_svi_{e.episode}_stacomp.csv')
        Sc = Sc[Sc.unit == 'pooled']
        st = pickle.load(open(f'{C.DATA}/stacks/{pre}_svi_{e.episode}.pkl', 'rb'))['units']['pooled']
        lag = int(st['lag']); n = int(u.n)
        for s, lab, a, b in targets:
            if s not in st['stk'] or not (a[:10] <= e.start and e.end <= b[:10]):
                continue
            stk = st['stk'][s]; tpl = np.roll(st['tpl'][s], lag, axis=1)
            for c, cn in enumerate('NEZ'):
                r = Sc[(Sc.sta == s) & (Sc.comp == cn)]
                if r.empty:
                    continue
                r = r.iloc[0]
                # Amendment S2: ignore guard2 drops whose only reason is '<0.05xsta' (healthy verticals, tilt-noisy horizontals)
                reason = str(r.drop_reason) if 'drop_reason' in r.index and isinstance(r.drop_reason, str) else ''
                drop = bool(r.guard_dropped) if 'guard_dropped' in r.index else False
                if drop and reason == '<0.05xsta':
                    drop = False
                sig = r.sigma_ctrl_nm * 1e-9
                if drop or not np.isfinite(sig) or sig <= 0 or r.n_used < 5:
                    continue
                k = (s, lab, c)
                a_ = acc.setdefault(k, dict(S=0, T=0, n=0, v=0, eps=[]))
                w = 1.0 / sig ** 2                                      # Amendment S2: inverse-variance weights
                a_['S'] = a_['S'] + w * stk[c]; a_['T'] = a_['T'] + w * tpl[c]; a_['n'] += n; a_['v'] += w
                a_['eps'].append(e.episode)
            erows.append(dict(episode=e.episode, table=pre, n=n, lag=lag, m0=u.thrust_m0, sta=s, epoch=lab))
    pd.DataFrame(erows).to_csv(f'{C.TB}/scan_calib_episodes.csv', index=False)
    # joint grand fit
    num = den = 0.0
    G_ = {}
    for k, a_ in acc.items():
        Sg = a_['S'] / a_['v']; Tg = a_['T'] / a_['v']; sg = 1.0 / np.sqrt(a_['v'])
        ts = (Tg[FITM] * Sg[FITM]).sum(); tt = (Tg[FITM] ** 2).sum()
        G_[k] = (ts, tt, sg, len(a_['eps']), a_['n'])
        num += ts / sg ** 2; den += tt / sg ** 2
    M0g = num / den
    for (s, lab, c), (ts, tt, sg, ne, nn) in G_.items():
        g = ts / (M0g * tt)
        rows.append(dict(sta=s, epoch=lab, comp='NEZ'[c], n_episodes=ne, n_lfe=nn, sigma_grand_nm=sg * 1e9,
                         tpl_rms_per1e12_nm=np.sqrt(tt / FITM.sum()) * 1e21, g=g, g_snr=ts / (sg * np.sqrt(tt)),
                         keep=bool(g > 0 and abs(np.log(g)) <= np.log(2))))
    R = pd.DataFrame(rows).sort_values(['sta', 'epoch', 'comp'])
    R['M0_grand_per_lfe'] = M0g
    for s, lab, a, b in targets:
        R.loc[(R.sta == s) & (R.epoch == lab), 'epoch_start'] = a
        R.loc[(R.sta == s) & (R.epoch == lab), 'epoch_end'] = b
    R.to_csv(f'{C.TB}/sc_gain.csv', index=False)
    return R, M0g


def day_noise(args):
    f, = args
    d = np.load(f)['d']
    out = []
    for c in range(3):
        b = S.bandpass_runs(d[c].astype(np.float64))
        w = b[:86400 // 400 * 400].reshape(-1, 400)
        ok = np.isfinite(w).mean(axis=1) >= 0.95
        out.append(np.nanmedian(np.sqrt(np.nanmean(w[ok] ** 2, axis=1))) if ok.any() else np.nan)
    return f, out


def noise():
    from multiprocessing import Pool
    files = []
    for s in STAS:
        files += sorted(glob.glob(f'{C.DATA}/disp1sps_v2/{s}/*/*.npz')) + sorted(glob.glob(f'{S.V3}/{s}/2011/*.npz'))
    rows = []
    with Pool(8) as p:
        for f, o in p.imap_unordered(day_noise, [(f,) for f in files], chunksize=8):
            parts = f.split('/')
            rows.append(dict(sta=parts[-3], year=int(parts[-2]), doy=int(parts[-1][:3]), root=parts[-4],
                             N=o[0], E=o[1], Z=o[2]))
    D = pd.DataFrame(rows)
    D['date'] = pd.to_datetime(D.year.astype(str) + '-' + D.doy.astype(str), format='%Y-%j')
    D.to_csv(f'{S.SD}/noise_days.csv', index=False)
    ref = D.drop_duplicates(['sta', 'date']).groupby('sta').Z.median().rename('Z_ref_all_days')
    ref.to_csv(f'{C.TB}/scan_noise_ref.csv')
    V = D[(D.root == 'disp1sps_v3') & (D.year == 2011)].copy()
    V['month'] = V.date.dt.month
    M = V.groupby(['sta', 'month'])[['N', 'E', 'Z']].median().reset_index()
    netmed = M[M.sta.isin(S.PILOT + S.NOISE_ONLY)].groupby('month').Z.median().rename('Z_net_median')
    M = M.join(ref, on='sta').join(netmed, on='month')
    M['n_days'] = V.groupby(['sta', 'month']).size().values
    M['exclude'] = (M.Z > 3 * M.Z_ref_all_days) | (M.Z > 5 * M.Z_net_median)
    for c in ('N', 'E', 'Z', 'Z_ref_all_days', 'Z_net_median'):
        M[c] = M[c] * 1e9
    M.to_csv(f'{C.TB}/scan_noise_month.csv', index=False)
    return M


if __name__ == '__main__':
    what = sys.argv[1] if len(sys.argv) > 1 else 'all'
    V = {}
    if what in ('all', 'tpl'):
        N, A = build_templates(); V['n_nodes'] = len(N)
    if what in ('all', 'gain'):
        R, M0g = calibrate()
        print(R.round(3).to_string())
        P = R[R.sta.isin(S.PILOT)]
        P11 = P[~((P.sta == 'CN.SNB') & P.epoch.str.startswith('BH'))]      # count the 21 = 7 x 3 main 2011 epochs
        V.update(M0_grand_per_lfe=M0g, n_keep_of_21=int(P11.keep.sum()), n_21=int(len(P11)),
                 z_gains={r.sta: float(r.g) for r in P11[P11.comp == 'Z'].itertuples()})
        zg = V['z_gains']
        V['accept_12of21'] = V['n_keep_of_21'] >= 12
        V['accept_Z_PGC_VGZ_LRIV_0.7_1.4'] = all(0.7 <= zg.get(s, np.nan) <= 1.4 for s in ('CN.PGC', 'CN.VGZ', 'UW.LRIV'))
    if what in ('all', 'noise'):
        M = noise(); print(M.round(1).to_string()); V['noise_months_excluded'] = M[M.exclude][['sta', 'month']].values.tolist()
    fj = f'{C.TB}/scan_calib_verdict.json'
    old = json.load(open(fj)) if os.path.isfile(fj) else {}
    old.update(V)
    json.dump(old, open(fj, 'w'), indent=1, default=float)
    print(json.dumps(V, indent=1, default=float))
