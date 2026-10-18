#!/usr/bin/env python
"""S13 stack + template fit (rules: ../PLAN.md §3).

usage: s13_stack.py EPISODE START END STATIONS [--fams 10,56,...] [--ndraw 500] [--tag pilot]
For units = pooled (all SVI families) + each listed family (+ pooled P>=0.5, pooled odd/even-day halves):
stack 0.02-0.05 Hz displacement (window -150..+250 s about origin, 1 sps) per station-component; scrambled controls
(a) random time within the same UTC day (ndraw draws, same QC) and (b) +-1 day shifts; joint thrust-template fit.
Writes tables/<tag>_<EPISODE>_units.csv, tables/<tag>_<EPISODE>_stacomp.csv, data/stacks/<tag>_<EPISODE>.npz."""
import os, sys, argparse, json
os.environ.setdefault('OMP_NUM_THREADS', '1')
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s13_common as C
import s13_green as G

W0, W1 = -150, 250                  # window about origin (s)
FIT0, FIT1 = -60, 120               # template-fit window (s)
LAGS = np.arange(-8, 9)
NWIN = W1 - W0 + 1
TREL = np.arange(W0, W1 + 1)
MECHS = {'thrust': (315, 12, 90), 'strikeslip': (315, 90, 0), 'steepdip': (315, 80, 90)}
DIL = 150                           # mask dilation (s) after the band-pass


def load_station(netsta, days):
    """Continuous band-passed N,E,Z over the days (NaN where missing), t0 = first day 00:00 UTC (epoch s)."""
    X = np.full((3, len(days) * 86400), np.nan, np.float64)
    for i, d in enumerate(days):
        f = f'{C.DISP}/{netsta}/{d.year}/{d.dayofyear:03d}.npz'
        if os.path.isfile(f):
            X[:, i * 86400:(i + 1) * 86400] = np.load(f)['d']
    bad = ~np.isfinite(X)
    B = G.bandpass(np.where(bad, 0.0, X), 1.0)
    # dilate the mask
    from scipy.ndimage import binary_dilation
    for c in range(3):
        if bad[c].any():
            m = binary_dilation(bad[c], structure=np.ones(2 * DIL + 1, bool))
            B[c, m] = np.nan
    return B


def gather(B, idx):
    """windows (n, 3, NWIN) for start indices idx (may be out of range -> NaN)."""
    n = len(idx)
    out = np.full((n, 3, NWIN), np.nan)
    ok = (idx >= 0) & (idx + NWIN <= B.shape[1])
    if ok.any():
        ii = idx[ok][:, None] + np.arange(NWIN)[None, :]
        out[ok] = B[:, ii].transpose(1, 0, 2)
    return out


def stack_windows(W, thr, w=None):
    """QC'd mean over windows. W (n,3,NWIN); thr (3,) RMS thresholds. Returns stack (3,NWIN), n_used (3,)."""
    nanfrac = np.isnan(W).mean(axis=2)
    rms = np.sqrt(np.nanmean(W ** 2, axis=2))
    keep = (nanfrac <= 0.05) & (rms <= thr[None, :])
    Wf = np.where(np.isnan(W), 0.0, W)
    ww = np.ones(W.shape[0]) if w is None else np.asarray(w, float)
    k = keep * ww[:, None]
    s = (Wf * k[:, :, None]).sum(0) / np.maximum(k.sum(0), 1e-12)[:, None]
    return s, keep.sum(0)


def templates(src_ll_dep, stations, coords, mech):
    lat0, lon0, dep = src_ll_dep
    M = G.mt_ned(*mech, 1.0)
    t = np.arange(-400, 600, 0.25)
    T = {}
    for s in stations:
        n, e = G.ned_offset(lat0, lon0, coords.loc[s, 'lat'], coords.loc[s, 'lon'])
        u = G.bandpass(G.synth(M, (0, 0, dep * 1e3), (n, e, 0.0), t), 4.0)
        T[s] = np.vstack([np.interp(TREL, t, c) for c in u])
    return T


def elem_templates(src_ll_dep, stations, coords):
    """5 deviatoric elementary moment tensors (unit)."""
    lat0, lon0, dep = src_ll_dep
    E = [np.diag([1, -1, 0]), np.array([[0, 1, 0], [1, 0, 0], [0, 0, 0]]), np.array([[0, 0, 1], [0, 0, 0], [1, 0, 0]]),
         np.array([[0, 0, 0], [0, 0, 1], [0, 1, 0]]), np.diag([-1, -1, 2]) / np.sqrt(3)]
    t = np.arange(-400, 600, 0.25)
    out = []
    for M in E:
        T = {}
        for s in stations:
            n, e = G.ned_offset(lat0, lon0, coords.loc[s, 'lat'], coords.loc[s, 'lon'])
            u = G.bandpass(G.synth(np.asarray(M, float), (0, 0, dep * 1e3), (n, e, 0.0), t), 4.0)
            T[s] = np.vstack([np.interp(TREL, t, c) for c in u])
        out.append(T)
    return out


FITM = (TREL >= FIT0) & (TREL <= FIT1)


def fit(stacks, T, wts, lags=LAGS):
    """Joint weighted fit over station-components. stacks/T: dict sta->(3,NWIN); wts: dict sta->(3,) 1/sigma^2.
    Returns M0 (signed, at the lag maximising the explained energy), lag, VR."""
    best = None
    num_d = sum((wts[s][:, None] * stacks[s][:, FITM] ** 2).sum() for s in stacks)
    for L in lags:
        sd = ss = 0.0
        for s in stacks:
            tt = np.roll(T[s], L, axis=1)[:, FITM]
            sd += (wts[s][:, None] * tt * stacks[s][:, FITM]).sum()
            ss += (wts[s][:, None] * tt * tt).sum()
        expl = sd * sd / ss
        if best is None or expl > best[0]:
            best = (expl, sd / ss, L)
    expl, m0, L = best
    return m0, L, expl / num_d if num_d > 0 else np.nan


def fit_free(stacks, Ts, wts, L):
    rows, rhs = [], []
    for s in stacks:
        w = np.sqrt(wts[s])[:, None]
        rhs.append((w * stacks[s][:, FITM]).ravel())
        rows.append(np.vstack([(w * np.roll(T[s], L, axis=1)[:, FITM]).ravel() for T in Ts]).T)
    A = np.vstack(rows); b = np.concatenate(rhs)
    m, *_ = np.linalg.lstsq(A, b, rcond=None)
    vr = 1 - ((b - A @ m) ** 2).sum() / (b ** 2).sum()
    return m, vr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('episode'); ap.add_argument('start'); ap.add_argument('end'); ap.add_argument('stations')
    ap.add_argument('--fams', default=''); ap.add_argument('--ndraw', type=int, default=500)
    ap.add_argument('--tag', default='pilot'); ap.add_argument('--pcut', type=float, default=0.75)
    a = ap.parse_args()
    lf = f'{C.LOGD}/stack_{a.tag}.log'
    rng = np.random.default_rng(C.SEED)
    stations = a.stations.split(',')
    d0 = pd.Timestamp(a.start) - pd.Timedelta(days=2); d1 = pd.Timestamp(a.end) + pd.Timedelta(days=2)
    days = pd.date_range(d0, d1, freq='D')
    T0 = d0.value / 1e9
    R = C.svi_rows()
    R = R[(R.series == 'data') & R.counted_once & (R.date >= d0.strftime('%Y-%m-%d')) & (R.date <= d1.strftime('%Y-%m-%d'))]
    F = C.svi_families()
    shift = C.origin_shift_svi(sorted(R.fam.unique()))
    R = R.assign(to=R.t + R.fam.map(shift).fillna(0.0))
    coords = C.station_coords()
    # neighbour count (any family, P>=0.5) within +-10 s
    allt = np.sort(R[R.P >= 0.5].to.values)
    fams = [int(x) for x in a.fams.split(',') if x]
    units = {'pooled': R[R.P >= a.pcut], 'pooled_P050': R[R.P >= 0.5]}
    pr = units['pooled']
    dd = pd.to_datetime(pr.date).dt.day
    units['pooled_odd'] = pr[dd % 2 == 1]; units['pooled_even'] = pr[dd % 2 == 0]
    for f in fams:
        units[f'fam{f}'] = pr[pr.fam == f]
    C.log(f'{a.episode}: units ' + ', '.join(f'{k}:{len(v)}' for k, v in units.items()), lf)
    # load stations
    B = {}
    for s in stations:
        B[s] = load_station(s, days)
        C.log(f'{s} loaded, finite Z frac {np.isfinite(B[s][2]).mean():.2f}', lf)
    stations = [s for s in stations if np.isfinite(B[s][2]).mean() > 0.2]
    rows, srows, store = [], [], {}
    for uname, U in units.items():
        if len(U) < 10:
            continue
        to = U.to.values; P = U.P.values
        nb = (np.searchsorted(allt, to + 10) - np.searchsorted(allt, to - 10) - 1).mean()
        # location: count-weighted mean of the unit's families
        fl = F.loc[U.fam.values]
        src = (fl.lat.mean(), fl.lon.mean(), fl.depth_km.mean())
        idx0 = np.round(to - T0 + W0).astype(np.int64)
        dayi = np.floor((to - T0) / 86400).astype(np.int64)
        draws = [np.round(dayi * 86400 + rng.uniform(0, 86400, len(to)) + W0).astype(np.int64) for _ in range(a.ndraw)]
        stk, stkP, ctrl, shp, shm, nused = {}, {}, {}, {}, {}, {}
        thr = {}
        for s in stations:
            Wd = gather(B[s], idx0)
            Wc = np.concatenate([gather(B[s], draws[i]) for i in range(min(20, a.ndraw))])
            allw = np.concatenate([Wd, Wc])
            rms = np.sqrt(np.nanmean(allw ** 2, axis=2))
            thr[s] = 3 * np.nanmedian(rms, axis=0)
            stk[s], nused[s] = stack_windows(Wd, thr[s])
            stkP[s], _ = stack_windows(Wd, thr[s], P)
            ctrl[s] = np.stack([stack_windows(gather(B[s], dr), thr[s])[0] for dr in draws])
            shp[s], _ = stack_windows(gather(B[s], idx0 + 86400), thr[s])
            shm[s], _ = stack_windows(gather(B[s], idx0 - 86400), thr[s])
        sig = {s: ctrl[s][:, :, FITM].std(axis=0).mean(axis=1) for s in stations}   # (3,)
        wts = {s: 1.0 / np.maximum(sig[s], 1e-15) ** 2 for s in stations}
        res = {}
        for mname, mech in MECHS.items():
            T = templates(src, stations, coords, mech)
            m0, L, vr = fit(stk, T, wts)
            cm = np.array([fit({s: ctrl[s][i] for s in stations}, T, wts)[0] for i in range(a.ndraw)])
            mp = fit(shp, T, wts)[0]; mm = fit(shm, T, wts)[0]
            mP = fit(stkP, T, wts)[0]
            res[mname] = dict(m0=m0, lag=L, vr=vr, ctrl_mean=cm.mean(), ctrl_sd=cm.std(), ctrl_p99=np.quantile(cm, .99),
                              z=(m0 - cm.mean()) / cm.std(), p=(1 + (cm >= m0).sum()) / (1 + len(cm)),
                              m0_shift_p1=mp, m0_shift_m1=mm, m0_Pw=mP)
            if mname == 'thrust':
                Tthr, Lthr = T, L
                cvr = np.array([fit({s: ctrl[s][i] for s in stations}, T, wts)[2] for i in range(min(200, a.ndraw))])
                res[mname]['ctrl_vr_p99'] = np.quantile(cvr, .99)
                # sign agreement among the larger template station-components
                amps, agree = [], []
                for s in stations:
                    tt = np.roll(T[s], L, axis=1)[:, FITM]
                    for c in range(3):
                        amps.append(np.abs(tt[c]).max())
                        agree.append(np.sign((tt[c] * stk[s][c, FITM]).sum()))
                amps = np.array(amps); agree = np.array(agree)
                big = amps >= np.median(amps)
                res[mname]['sign_agree_big'] = float((agree[big] > 0).mean())
        Ts = elem_templates(src, stations, coords)
        mfree, vrfree = fit_free(stk, Ts, wts, Lthr)
        row = dict(episode=a.episode, unit=uname, n=len(U), n_nb10=nb, src_lat=src[0], src_lon=src[1], src_dep=src[2],
                   n_sta=len(stations), stations=','.join(stations), vr_free=vrfree, mt_free=json.dumps(list(np.round(mfree, 1))))
        for mname, r in res.items():
            for k, v in r.items():
                row[f'{mname}_{k}'] = v
        row['thrust_m0_per_lfe_nbcorr'] = row['thrust_m0'] / (1 + nb)
        rows.append(row)
        for s in stations:
            q = np.quantile(np.abs(ctrl[s][:, :, (TREL >= -30) & (TREL <= 90)]).max(axis=2), 0.99, axis=0)
            pk = np.abs(stk[s][:, (TREL >= -30) & (TREL <= 90)]).max(axis=1)
            for c, cn in enumerate('NEZ'):
                srows.append(dict(episode=a.episode, unit=uname, sta=s, comp=cn, n_used=int(nused[s][c]),
                                  sigma_ctrl_nm=sig[s][c] * 1e9, thr_rms_nm=thr[s][c] * 1e9 / 3,
                                  peak_nm=pk[c] * 1e9, ctrl_peak_p99_nm=q[c] * 1e9, det=bool(pk[c] > q[c]),
                                  tpl_amp_per1e12_nm=np.abs(Tthr[s][c]).max() * 1e21))
        store[uname] = dict(stk={s: stk[s] for s in stations}, stkP={s: stkP[s] for s in stations},
                            shp={s: shp[s] for s in stations}, shm={s: shm[s] for s in stations},
                            c_lo={s: np.quantile(ctrl[s], .025, axis=0) for s in stations},
                            c_hi={s: np.quantile(ctrl[s], .975, axis=0) for s in stations},
                            c_ex={s: ctrl[s][:5] for s in stations},
                            tpl={s: Tthr[s] for s in stations}, lag=Lthr)
        C.log(f'{a.episode} {uname} n={len(U)} thrust M0={res["thrust"]["m0"]:.3g} z={res["thrust"]["z"]:.2f} '
              f'p={res["thrust"]["p"]:.3f} VR={res["thrust"]["vr"]:.2f} lag={Lthr}', lf)
    pd.DataFrame(rows).to_csv(f'{C.TB}/{a.tag}_{a.episode}_units.csv', index=False)
    pd.DataFrame(srows).to_csv(f'{C.TB}/{a.tag}_{a.episode}_stacomp.csv', index=False)
    os.makedirs(f'{C.DATA}/stacks', exist_ok=True)
    import pickle
    with open(f'{C.DATA}/stacks/{a.tag}_{a.episode}.pkl', 'wb') as fh:
        pickle.dump(dict(trel=TREL, units=store), fh)
    open(f'{C.OUT}/STACK_{a.tag}_{a.episode}_DONE', 'w').write('ok\n')


if __name__ == '__main__':
    main()
