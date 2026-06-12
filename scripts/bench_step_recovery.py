#!/usr/bin/env python
"""Step-injection benchmark: how well does each dv/v estimator recover a sudden coseismic step?

Injects a known stretch step at --t0 into REAL daily coda stacks (post-t0 traces are time-stretched),
then recovers dv/v(t) with several estimators and compares step shape / rise time / leakage:
  prod30_svd : production method — 30-cal-day trailing stack, full-record SVD-Wiener, stretch (dvv_roll30cal.py)
  prod30     : same without SVD-Wiener (causal)
  roll10     : 10-cal-day trailing stack (min 3), no SVD
  daily      : per-day stack vs all-time ref (no time smoothing) — the raw measurement
  daily_tv   : cross-family-median daily series + total-variation denoise (step-preserving, acausal)
  daily_kal  : cross-family-median daily series + causal robust Kalman with jump acceptance
Cross-family per-day MEDIAN (spatial, allowed) is applied to every method before metrics.
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, sys
import numpy as np, pandas as pd
sys.path.insert(0, "src")
from tremorferometry.dvv import stretch_dvv  # noqa: E402

_C = {}


def _init(t, fs, w1, w2, t0, inj_a):
    _C.update(t=t, fs=fs, tmin=w1 - t[0], tmax=w2 - t[0], t0=np.datetime64(t0), inj=inj_a,
              prenoise=(t < -1.0))


def _inject(S, Dn):
    """Time-stretch post-t0 traces about the LFE origin (t=0 of the npz time axis)."""
    t = _C["t"]; a = _C["inj"]
    S = S.copy()
    post = Dn > _C["t0"]
    for i in np.where(post)[0]:
        S[i] = np.interp(t * (1.0 + a), t, S[i])
    return S


def _stretch_rows(R, ref):
    out = np.full(len(R), np.nan); cc = np.full(len(R), np.nan)
    for i in range(len(R)):
        try:
            r = stretch_dvv(ref, R[i], fs=_C["fs"], t_min=_C["tmin"], t_max=_C["tmax"],
                            eps_max=0.02, n_eps=201)
            out[i] = r.dvv; cc[i] = r.cc_max
        except Exception:
            pass
    return out, cc


def _roll(S, Dn, days, minst, svd):
    """Trailing calendar-day rolling stack + (optional) SVD-Wiener + stretch — mirrors dvv_roll30cal.py."""
    n = len(S); idx = np.arange(n)
    cs = np.vstack([np.zeros((1, S.shape[1])), np.cumsum(S, axis=0)])
    los = np.searchsorted(Dn, Dn - np.timedelta64(days, 'D'), side='right')
    counts = idx + 1 - los
    keep = counts >= minst
    if keep.sum() < 8:
        return None
    R = (cs[idx + 1] - cs[los]) / counts[:, None]
    R = R[keep]; Dr = Dn[keep]
    pk = np.max(np.abs(R), axis=1, keepdims=True); pk[pk == 0] = 1.0; Rn = R / pk
    if svd:
        ne, ns = Rn.shape; sig = Rn[:, _C["prenoise"]].std()
        try:
            U, sv, Vt = np.linalg.svd(Rn, full_matrices=False)
            nmf = sig * (np.sqrt(ne) + np.sqrt(ns)); w = sv ** 2 / (sv ** 2 + nmf ** 2)
            Rn = (U * (sv * w)) @ Vt
        except np.linalg.LinAlgError:
            pass
    ref = Rn.mean(0)
    dvv, cc = _stretch_rows(Rn, ref)
    return Dr, dvv, cc


def _one_arm(S, Dn):
    res = {}
    r = _roll(S, Dn, 30, 5, svd=True)
    if r: res['prod30_svd'] = r
    r = _roll(S, Dn, 30, 5, svd=False)
    if r: res['prod30'] = r
    r = _roll(S, Dn, 10, 3, svd=False)
    if r: res['roll10'] = r
    pk = np.max(np.abs(S), axis=1, keepdims=True); pk[pk == 0] = 1.0; Sn = S / pk
    ref = Sn.mean(0)
    dvv, cc = _stretch_rows(Sn, ref)
    res['daily'] = (Dn, dvv, cc)
    return res


def _patch(arg):
    patch, S, dstr = arg
    order = np.argsort(dstr)
    S = S[order].astype(np.float64)
    Dn = pd.to_datetime(np.asarray(dstr)[order]).values.astype('datetime64[D]')
    return patch, {'ctrl': _one_arm(S, Dn), 'inj': _one_arm(_inject(S, Dn), Dn)}


def tv_denoise(y, lam, n_iter=200, rho=2.0):
    """1-D total-variation denoise via ADMM (step-preserving)."""
    from scipy.sparse import diags
    from scipy.sparse.linalg import factorized
    n = len(y)
    D = diags([-np.ones(n - 1), np.ones(n - 1)], [0, 1], shape=(n - 1, n)).tocsc()
    A = (diags([np.ones(n)], [0]) + rho * (D.T @ D)).tocsc()
    solve = factorized(A)
    x = y.copy(); z = np.zeros(n - 1); u = np.zeros(n - 1)
    for _ in range(n_iter):
        x = solve(y + rho * (D.T @ (z - u)))
        Dx = D @ x
        z = np.sign(Dx + u) * np.maximum(np.abs(Dx + u) - lam / rho, 0)
        u += Dx - z
    return x


def kalman_jump(y, sig_meas, q_scale=1 / 30., k_gate=4.0):
    """Causal random-walk Kalman; a gated innovation confirmed by the NEXT sample = accepted jump."""
    R = sig_meas ** 2; Q = (sig_meas * q_scale) ** 2
    x = y[0]; P = R; out = np.empty(len(y))
    for i, yi in enumerate(y):
        P = P + Q
        S = P + R; innov = yi - x
        if abs(innov) > k_gate * np.sqrt(S) and i + 1 < len(y) and abs(y[i + 1] - x) > k_gate * np.sqrt(S):
            x = yi; P = R          # accept as jump
        elif abs(innov) > k_gate * np.sqrt(S):
            pass                   # isolated outlier: ignore
        else:
            K = P / S; x = x + K * innov; P = (1 - K) * P
        out[i] = x
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--npz", required=True); ap.add_argument("--station", required=True)
    ap.add_argument("--t0", default="2016-07-01"); ap.add_argument("--amp", type=float, default=-0.0015)
    ap.add_argument("--window", nargs=2, type=float, default=[2.0, 4.0])
    ap.add_argument("--n-patches", type=int, default=12); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--out-prefix", default="bench_step")
    args = ap.parse_args()

    d = np.load(args.npz, allow_pickle=True)
    t = d["t"]; fs = float(d["fs"]); stacks = d["stacks"]; pat = d["patches"].astype(str)
    dates = pd.to_datetime(d["dates"]).values.astype('datetime64[D]')
    cnt = pd.Series(pat).value_counts()
    top = cnt.index[:args.n_patches].tolist()

    # ---- calibrate injection: what does the estimator read for a stretch about the LFE origin?
    _init(t, fs, *args.window, args.t0, 0.0)
    m = pat == top[0]
    ref0 = (stacks[m].astype(np.float64) / np.max(np.abs(stacks[m]), axis=1, keepdims=True)).mean(0)
    probe = np.interp(t * (1.0 + args.amp), t, ref0)
    r = stretch_dvv(ref0, probe, fs=fs, t_min=args.window[0] - t[0], t_max=args.window[1] - t[0],
                    eps_max=0.02, n_eps=801)
    meas_per_inj = r.dvv / args.amp
    print(f"calibration: injected stretch {args.amp:+.4%} reads as dvv {r.dvv:+.4%} "
          f"(measured/injected = {meas_per_inj:.2f}) [origin-offset scale factor]", flush=True)
    target = r.dvv  # what a perfect estimator should plateau at

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed
    tasks = [(p, stacks[pat == p], dates[pat == p]) for p in top]
    ctx = mp.get_context("spawn")
    per = {}
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx, initializer=_init,
                             initargs=(t, fs, args.window[0], args.window[1], args.t0, args.amp)) as ex:
        futs = [ex.submit(_patch, tk) for tk in tasks]
        for f in as_completed(futs):
            patch, res = f.result()
            per[patch] = res
    print(f"{len(per)} patches processed", flush=True)

    # ---- cross-patch daily median per method, per arm (ctrl / inj)
    series = {'ctrl': {}, 'inj': {}}
    for arm in ['ctrl', 'inj']:
        for meth in ['prod30_svd', 'prod30', 'roll10', 'daily']:
            frames = []
            for patch, res in per.items():
                if meth not in res[arm]: continue
                Dn, dvv, cc = res[arm][meth]
                frames.append(pd.DataFrame({'date': Dn, 'dvv': dvv}).set_index('date').dvv.rename(patch))
            M = pd.concat(frames, axis=1, sort=True)
            series[arm][meth] = M.median(axis=1).dropna()
        # post-hoc filters applied independently per arm (nonlinear — must see each arm's own data)
        sd = series[arm]['daily']
        sig = sd.diff().std() / np.sqrt(2)
        series[arm]['daily_tv'] = pd.Series(tv_denoise(sd.values, lam=4 * sig), index=sd.index)
        series[arm]['daily_kal'] = pd.Series(kalman_jump(sd.values, sig), index=sd.index)

    # ---- DIFFERENTIAL step response: inj − ctrl (natural variability cancels exactly)
    diff = {}
    for meth in series['inj']:
        a, b = series['inj'][meth].align(series['ctrl'][meth], join='inner')
        diff[meth] = (a - b).dropna()

    # ---- metrics on the differential response
    t0 = np.datetime64(args.t0)
    rows = []
    for meth, s in diff.items():
        pre = s[s.index < t0 - np.timedelta64(2, 'D')]
        base = pre.median(); noise = pre.std()
        post = s[(s.index > t0 + np.timedelta64(35, 'D')) & (s.index < t0 + np.timedelta64(150, 'D'))]
        amp = post.median() - base
        lead = s[(s.index >= t0 - np.timedelta64(15, 'D')) & (s.index < t0)]
        leak = (lead.median() - base) / target if len(lead) else np.nan
        rel = (s[s.index > t0] - base) / target
        rise = np.nan
        vals = rel.values; idxs = rel.index
        for i in range(len(vals) - 4):
            if np.all(vals[i:i + 5] >= 0.8):
                rise = (idxs[i] - t0) / np.timedelta64(1, 'D'); break
        # estimator's own stability: pre-step noise of the RAW (inj-arm) series
        raw_pre = series['inj'][meth]
        raw_pre = raw_pre[(raw_pre.index >= t0 - np.timedelta64(400, 'D')) & (raw_pre.index < t0 - np.timedelta64(10, 'D'))]
        rows.append(dict(method=meth, series_noise_pct=100 * raw_pre.std(), diff_pre_noise_pct=100 * noise,
                         amp_frac=amp / target, rise80_days=rise, pre_leak_frac=leak))
    met = pd.DataFrame(rows)
    print('\ninjected step reads as %+.4f%% at t0=%s\n' % (100 * target, args.t0))
    print(met.to_string(index=False, float_format=lambda x: f'{x:.3f}'))
    met.to_csv(f'data/{args.out_prefix}_{args.station}_metrics.csv', index=False)

    # ---- figure
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 1, figsize=(13, 9), sharex=False)
    zoom = (t0 - np.timedelta64(180, 'D'), t0 + np.timedelta64(240, 'D'))
    colors = {'prod30_svd': 'tab:red', 'prod30': 'tab:orange', 'roll10': 'tab:green',
              'daily': '0.8', 'daily_tv': 'tab:blue', 'daily_kal': 'tab:purple'}
    ax = axes[0]
    for meth, s in diff.items():
        ss = s[(s.index >= zoom[0]) & (s.index <= zoom[1])]
        lw = 0.8 if meth == 'daily' else 1.6
        ax.plot(ss.index, 100 * ss.values, color=colors[meth], lw=lw, label=meth,
                alpha=0.65 if meth == 'daily' else 1.0)
    ax.plot([zoom[0], t0, t0, zoom[1]], 100 * np.array([0, 0, target, target]),
            'k-', lw=2.2, alpha=0.5, label='truth')
    ylo, yhi = ax.get_ylim()
    ax.vlines(t0, ylo, yhi, color='k', ls='--', lw=1)
    ax.set_xlim(zoom); ax.legend(ncol=4, fontsize=9); ax.set_ylabel('Δdv/v (inj − ctrl, %)')
    ax.set_title(f'{args.station}: DIFFERENTIAL step response, injected {100*target:+.3f}% at {args.t0}')
    ax = axes[1]
    for meth, s in series['ctrl'].items():
        ax.plot(s.index, 100 * s.values, color=colors[meth], lw=0.8 if meth == 'daily' else 1.4,
                label=meth, alpha=0.65 if meth == 'daily' else 1.0)
    ylo, yhi = ax.get_ylim()
    ax.vlines(t0, ylo, yhi, color='k', ls='--', lw=1)
    ax.set_ylabel('dv/v (%)'); ax.set_title('control arm (no injection) — natural variability, full record')
    fig.tight_layout()
    out = f'figures/{args.out_prefix}_{args.station}.png'
    fig.savefig(out, dpi=140)
    print('->', out)


if __name__ == "__main__":
    main()
