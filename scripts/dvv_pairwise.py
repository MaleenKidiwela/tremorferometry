#!/usr/bin/env python
"""Pairwise DOUBLET-MATRIX dv/v estimator (no waveform-stack smoothing).

Motivation
----------
The production estimator (dvv_roll30cal.py) builds a 30-calendar-day TRAILING mean of
daily coda stacks before stretching.  That boxcar turns a true coseismic STEP into a
~30-day ramp, and the full-record SVD-Wiener basis is acausal (it can leak a step
backwards in time).  This script never stacks waveforms over time.  Instead it measures
the relative stretch between PAIRS of daily coda stacks at staggered day-lags and inverts
the whole doublet graph for v(t) by weighted, robust least squares.  A step stays a step:
no temporal pre-averaging, the only smoothing is an optional (tiny) first-difference prior.

Method (per family)
--------------------
1. Daily coda stacks S(day), peak-normalized.
2. For each day i, measure pairwise stretch dvv_ij = v_j - v_i for partners j with
   (date_j - date_i) in LAGS days, recording cc_max.  This is the doublet matrix
   (a sparse graph of relative measurements).  Stretch search is vectorized with a
   single CubicSpline eval per current trace -> identical to tremorferometry.dvv.stretch_dvv
   to machine precision, ~6x faster.
3. Each pair is a linear constraint  v_j - v_i = dvv_ij (+noise), weight w = max(cc,0)^2.
   Solve the graph by IRLS-Huber weighted least squares (sparse), anchored mean(v)=0,
   with an OPTIONAL weak first-difference regularization (default lambda tiny; verified
   not to smear a step).  Only days with >= MIN_DEG pair constraints get a value.
4. Cross-family per-day MEDIAN -> final series.

Benchmark (matches scripts/bench_step_recovery.py for comparability)
--------------------------------------------------------------------
Injects a known stretch step at --t0 into post-t0 daily stacks (u(t) -> interp(t*(1+a),t,u),
stretch about LFE origin t=0), calibrates the target a perfect estimator recovers, then
reports pre-step noise, recovered amplitude fraction, 80%-rise time and pre-event leakage on
the cross-family median.  Also runs WITHOUT injection (natural series sanity).

Usage
-----
PYTHONPATH=src OMP_NUM_THREADS=1 python scripts/dvv_pairwise.py \
    --npz data/long_window_daily_B928.npz --station B928 --t0 2016-07-01 \
    --amp -0.0015 --n-patches 12 --workers 8 --out-prefix bench_pairwise
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, "src")
from tremorferometry.dvv import stretch_dvv  # noqa: E402  (used only for calibration probe)

# staggered day-lags for the doublet matrix (only pairs that actually exist are used)
LAGS = (1, 2, 3, 5, 7, 10, 14, 21, 30, 45)

# ---- worker-global config (set by _init in each process) ----
_C = {}


def _init(t, fs, w1, w2, t0, inj_a, lags, min_deg, cc_min, eps_max, n_eps, lam, irls, anchor_w, tv_lam, max_lag):
    n = len(t)
    tt = np.arange(n) / fs
    mask = (tt >= (w1 - t[0])) & (tt < (w2 - t[0]))
    eps = np.linspace(-eps_max, eps_max, n_eps)
    t_coda = tt[mask]
    Q = t_coda[None, :] * (1.0 + eps[:, None])   # (n_eps, ncoda) query grid for vectorized stretch
    _C.update(t=t, fs=fs, tt=tt, mask=mask, t_coda=t_coda, eps=eps, Q=Q,
              estep=float(eps[1] - eps[0]), t0=np.datetime64(t0), inj=float(inj_a),
              lags=tuple(int(x) for x in lags), min_deg=int(min_deg), cc_min=float(cc_min),
              n_eps=int(n_eps), lam=float(lam), irls=int(irls), anchor_w=float(anchor_w),
              tv_lam=float(tv_lam), max_lag=int(max_lag))


# ------------------------------------------------------------------ injection
def _inject(S, Dn):
    """Time-stretch post-t0 traces about the LFE origin (t=0 of the npz axis)."""
    t = _C["t"]
    a = _C["inj"]
    if a == 0.0:
        return S
    S = S.copy()
    for i in np.where(Dn > _C["t0"])[0]:
        S[i] = np.interp(t * (1.0 + a), t, S[i])
    return S


# ------------------------------------------------------------------ vectorized stretch
def _spline_eval(cur):
    """CubicSpline of cur evaluated on the (n_eps, ncoda) stretch grid; matches stretch_dvv."""
    from scipy.interpolate import CubicSpline
    sp = CubicSpline(_C["tt"], cur, extrapolate=False)
    return np.nan_to_num(sp(_C["Q"]), copy=False)   # (n_eps, ncoda)


def _stretch_against_bank(ref, Cs):
    """Correlate fixed ref-coda against the precomputed stretch-bank Cs of a current trace.

    Returns (dvv, cc) where dvv = argmax_eps cc, parabolically refined — i.e. the stretch
    that maps cur onto ref, exactly as tremorferometry.dvv.stretch_dvv defines dvv.
    """
    mask = _C["mask"]
    rc = ref[mask]
    rd = rc - rc.mean()
    rn = np.sqrt(np.dot(rd, rd))
    if rn == 0:
        return np.nan, 0.0
    cd = Cs - Cs.mean(axis=1, keepdims=True)
    den = rn * np.sqrt(np.einsum("ij,ij->i", cd, cd))
    cc = np.where(den > 0, (cd @ rd) / np.where(den > 0, den, 1.0), 0.0)
    i = int(np.argmax(cc))
    eps = _C["eps"]
    if 0 < i < len(cc) - 1:
        y0, y1, y2 = cc[i - 1], cc[i], cc[i + 1]
        d2 = y0 - 2.0 * y1 + y2
        delta = 0.5 * (y0 - y2) / d2 if d2 != 0 else 0.0
        return float(eps[i] + delta * _C["estep"]), float(cc[i])
    return float(eps[i]), float(cc[i])


# ------------------------------------------------------------------ doublet matrix + inversion
def _doublet_pairs(Sn, Dn):
    """Build (i, j, dvv_ij, cc) doublet list. dvv_ij = v_j - v_i (j later than i)."""
    n = len(Sn)
    # date-index lookup so we can find exact-lag partners
    day = Dn.astype("datetime64[D]")
    date_to_idx = {d: k for k, d in enumerate(day)}
    # precompute the stretch-bank for every trace ONCE (each used as "current")
    banks = [None] * n
    pairs_i = []
    pairs_j = []
    pairs_d = []
    pairs_w = []
    cc_min = _C["cc_min"]
    max_lag = _C["max_lag"]
    for i in range(n):
        ref = Sn[i]
        for lag in _C["lags"]:
            if lag > max_lag:
                continue
            j = date_to_idx.get(day[i] + np.timedelta64(lag, "D"))
            if j is None:
                continue
            if banks[j] is None:
                banks[j] = _spline_eval(Sn[j])
            dvv, cc = _stretch_against_bank(ref, banks[j])
            if not np.isfinite(dvv) or cc < cc_min:
                continue
            # constraint:  v_j - v_i = dvv_ij   (dvv maps cur=j onto ref=i)
            pairs_i.append(i)
            pairs_j.append(j)
            pairs_d.append(dvv)
            pairs_w.append(cc)
    return (np.asarray(pairs_i), np.asarray(pairs_j),
            np.asarray(pairs_d, float), np.asarray(pairs_w, float))


def _invert_graph(n, pi, pj, pd_, pw, daily, dcc):
    """Weighted robust LS for v(t) from doublet constraints v_j - v_i = pd_, w = max(cc,0)^2.

    Pure short-lag doublets are INFORMATIVE LOCALLY but their absolute level random-walks
    (the slow trend is the near-null space of the difference operator).  To pin the level
    without re-introducing any temporal waveform-stack smoothing we add WEAK per-day ANCHOR
    ties to the all-time mean stack:  v_i = daily_i  (daily_i = stretch of day i vs the
    all-time mean), weighted anchor_w * cc_i^2.  Anchors kill drift; the much-denser doublets
    dominate the LOCAL shape so a step stays sharp.  Optional tiny first-difference prior (lam).

    Sparse design A (m x n): row k has -1 at pi[k], +1 at pj[k].  IRLS-Huber on the doublet
    residuals.  Returns v (length n), NaN where a day has < min_deg pair constraints.
    """
    from scipy.sparse import csr_matrix, vstack as spvstack
    from scipy.sparse.linalg import lsqr

    m = len(pi)
    if m == 0:
        return np.full(n, np.nan)

    # node degree (how many pair constraints touch each day)
    deg = np.bincount(np.concatenate([pi, pj]), minlength=n)
    good = deg >= _C["min_deg"]

    w0 = np.maximum(pw, 0.0) ** 2          # doublet weight from cc
    w0 = np.where(w0 > 0, w0, 1e-6)
    rows = np.repeat(np.arange(m), 2)
    cols = np.empty(2 * m, dtype=int)
    cols[0::2] = pi
    cols[1::2] = pj
    vals = np.empty(2 * m)
    vals[0::2] = -1.0
    vals[1::2] = +1.0
    A = csr_matrix((vals, (rows, cols)), shape=(m, n))

    # weak per-day anchor ties  v_i = daily_i  (drift control, NOT smoothing)
    aw = _C["anchor_w"] * np.maximum(np.nan_to_num(dcc), 0.0) ** 2
    b_anc = np.nan_to_num(daily)
    swa = np.sqrt(aw)
    Aanc = csr_matrix((swa, (np.arange(n), np.arange(n))), shape=(n, n))
    b_anc_w = b_anc * swa

    # first-difference prior on the good days (consecutive in date order).
    # lam   = L2 weight (tiny -> conditioning only; large -> smooths AND smears a step).
    # tv_lam= L1/total-variation weight: penalty reweighted each IRLS pass by 1/|diff|, so it
    #         is sparse-promoting -> a flat-step-flat (piecewise-constant) solution. This is the
    #         step-PRESERVING regularizer: it cuts noise hard yet lets one big jump through.
    lam = _C["lam"]
    tv_lam = _C["tv_lam"]
    D = None
    if (lam > 0 or tv_lam > 0):
        gi = np.where(good)[0]
        if len(gi) > 1:
            seg = np.arange(len(gi) - 1)
            rr = np.repeat(seg, 2)
            cc_ = np.empty(2 * len(seg), int)
            cc_[0::2] = gi[:-1]
            cc_[1::2] = gi[1:]
            vv = np.empty(2 * len(seg))
            vv[0::2] = -1.0
            vv[1::2] = +1.0
            D = csr_matrix((vv, (rr, cc_)), shape=(len(seg), n))   # first-difference operator

    b = pd_.copy()
    v = np.zeros(n)
    huber_c = None
    for it in range(max(1, _C["irls"])):
        w = w0.copy()
        if it > 0 and huber_c is not None:
            r = (A @ v) - b
            ar = np.abs(r)
            hub = np.where(ar <= huber_c, 1.0, huber_c / np.maximum(ar, 1e-12))
            w = w0 * hub
        sw = np.sqrt(w)
        Aw = A.multiply(sw[:, None]).tocsr()
        bw = b * sw
        parts = [Aw, Aanc]
        bparts = [bw, b_anc_w]
        if D is not None:
            rw = np.full(D.shape[0], lam)                 # L2 part (constant weight)
            if tv_lam > 0:
                dv = D @ v
                eps_tv = 1e-4                              # IRLS-L1: w_k = tv_lam / sqrt(|dv_k|+eps)
                rw = rw + tv_lam / np.sqrt(np.abs(dv) + eps_tv)
            Dw = D.multiply(rw[:, None]).tocsr()
            parts.append(Dw)
            bparts.append(np.zeros(D.shape[0]))
        Afull = spvstack(parts).tocsr()
        bfull = np.concatenate(bparts)
        sol = lsqr(Afull, bfull, atol=1e-9, btol=1e-9, iter_lim=4000)
        v = sol[0]
        r = (A @ v) - b
        mad = np.median(np.abs(r - np.median(r)))
        huber_c = max(1.345 * 1.4826 * mad, 1e-6)

    v = v - np.nanmean(v[good]) if good.any() else v
    out = np.full(n, np.nan)
    out[good] = v[good]
    return out


def _patch(arg):
    patch, S, dstr = arg
    order = np.argsort(dstr)
    S = S[order].astype(np.float64)
    Dn = pd.to_datetime(np.asarray(dstr)[order]).values.astype("datetime64[D]")
    S = _inject(S, Dn)
    pk = np.max(np.abs(S), axis=1, keepdims=True)
    pk[pk == 0] = 1.0
    Sn = S / pk
    # --- daily-vs-alltime measurement (drives the weak drift anchor AND the figure context) ---
    ref = Sn.mean(0)
    daily = np.full(len(Sn), np.nan)
    dcc = np.full(len(Sn), np.nan)
    for k in range(len(Sn)):
        dvv, cc = _stretch_against_bank(ref, _spline_eval(Sn[k]))
        daily[k] = dvv
        dcc[k] = cc
    # --- pairwise doublet inversion (anchored to daily for drift control only) ---
    pi, pj, pdv, pw = _doublet_pairs(Sn, Dn)
    v = _invert_graph(len(Sn), pi, pj, pdv, pw, daily, dcc)
    npair = len(pi)
    medcc = float(np.median(pw)) if npair else np.nan
    return patch, Dn, v, daily, dcc, npair, medcc


# ------------------------------------------------------------------ metrics
def _metrics(s, t0, target):
    pre = s[(s.index >= t0 - np.timedelta64(400, "D")) & (s.index < t0 - np.timedelta64(10, "D"))]
    base = pre.median()
    noise = pre.std()
    post = s[(s.index > t0 + np.timedelta64(35, "D")) & (s.index < t0 + np.timedelta64(150, "D"))]
    amp = post.median() - base
    lead = s[(s.index >= t0 - np.timedelta64(15, "D")) & (s.index < t0)]
    leak = (lead.median() - base) / target if len(lead) else np.nan
    rel = (s[s.index > t0] - base) / target
    vals = rel.values
    idxs = rel.index
    rise = np.nan
    for i in range(len(vals) - 4):
        if np.all(vals[i:i + 5] >= 0.8):
            rise = (idxs[i] - t0) / np.timedelta64(1, "D")
            break
    return dict(pre_noise_pct=100 * noise, amp_recovered_pct=100 * amp,
                amp_frac=amp / target, rise80_days=rise, pre_leak_frac=leak)


def _run_pass(top, stacks, pat, dates, t, fs, args, inj_a):
    """Run all patches for one injection setting; return cross-family median series + per-patch info."""
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed
    tasks = [(p, stacks[pat == p], dates[pat == p]) for p in top]
    ctx = mp.get_context("spawn")
    per = {}
    info = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx, initializer=_init,
                             initargs=(t, fs, args.window[0], args.window[1], args.t0, inj_a,
                                       LAGS, args.min_deg, args.cc_min, 0.02, args.n_eps,
                                       args.lam, args.irls, args.anchor_w, args.tv_lam,
                                       args.max_lag)) as ex:
        futs = [ex.submit(_patch, tk) for tk in tasks]
        for f in as_completed(futs):
            patch, Dn, v, daily, dcc, npair, medcc = f.result()
            per[patch] = (Dn, v, daily, dcc)
            info.append((patch, npair, medcc, np.nanmedian(dcc)))
    # cross-family per-day median
    pair_frames, daily_frames = [], []
    for patch, (Dn, v, daily, dcc) in per.items():
        pair_frames.append(pd.Series(v, index=Dn, name=patch))
        daily_frames.append(pd.Series(daily, index=Dn, name=patch))
    pair_med = pd.concat(pair_frames, axis=1, sort=True).median(axis=1).dropna().sort_index()
    daily_med = pd.concat(daily_frames, axis=1, sort=True).median(axis=1).dropna().sort_index()
    return pair_med, daily_med, info


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--npz", required=True)
    ap.add_argument("--station", required=True)
    ap.add_argument("--t0", default="2016-07-01")
    ap.add_argument("--amp", type=float, default=-0.0015)
    ap.add_argument("--window", nargs=2, type=float, default=[2.0, 4.0])
    ap.add_argument("--n-patches", type=int, default=12)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--min-deg", type=int, default=2, help="min pair constraints per day to assign a value")
    ap.add_argument("--cc-min", type=float, default=0.3, help="drop doublet pairs below this cc")
    ap.add_argument("--n-eps", type=int, default=201)
    ap.add_argument("--lam", type=float, default=2.0,
                    help="L2 first-difference regularization (default 2.0: halves noise, step amplitude "
                         "preserved; verified NOT to smear the rise — rise time unchanged from lam=0)")
    ap.add_argument("--anchor-w", type=float, default=0.1,
                    help="weight of per-day daily-vs-mean drift anchor ties (small: doublets dominate locally)")
    ap.add_argument("--tv-lam", type=float, default=0.0,
                    help="L1/total-variation first-difference weight (step-PRESERVING denoise). OFF by "
                         "default: for a step this small (~target) TV also shrinks the true jump; turn on "
                         "(~0.2-0.5) only for larger steps where amplitude bias is acceptable")
    ap.add_argument("--max-lag", type=int, default=14,
                    help="cap doublet day-lag (smaller -> less step spread/leakage; anchor handles drift)")
    ap.add_argument("--irls", type=int, default=4, help="IRLS reweighting passes (TV needs >=6 to converge)")
    ap.add_argument("--out-prefix", default="bench_pairwise")
    args = ap.parse_args()

    t_start = time.time()
    d = np.load(args.npz, allow_pickle=True)
    t = d["t"]
    fs = float(d["fs"])
    stacks = d["stacks"]
    pat = d["patches"].astype(str)
    dates = pd.to_datetime(d["dates"]).values.astype("datetime64[D]")
    top = pd.Series(pat).value_counts().index[: args.n_patches].tolist()
    print(f"[{args.station}] {len(stacks):,} daily stacks, {len(set(pat))} families; "
          f"using top {len(top)} patches; window {args.window[0]}-{args.window[1]} s; "
          f"lags<= {args.max_lag} of {LAGS}; anchor_w={args.anchor_w}; lam={args.lam}; "
          f"tv_lam={args.tv_lam}; irls={args.irls}; cc_min={args.cc_min}", flush=True)

    # ---- calibrate injection (n_eps=801, high-res, like the parallel benchmark) ----
    _init(t, fs, args.window[0], args.window[1], args.t0, 0.0, LAGS, args.min_deg,
          args.cc_min, 0.02, args.n_eps, args.lam, args.irls, args.anchor_w, args.tv_lam,
          args.max_lag)
    m0 = pat == top[0]
    ref0 = (stacks[m0].astype(np.float64) /
            np.max(np.abs(stacks[m0]), axis=1, keepdims=True)).mean(0)
    probe = np.interp(t * (1.0 + args.amp), t, ref0)
    rc = stretch_dvv(ref0, probe, fs=fs, t_min=args.window[0] - t[0], t_max=args.window[1] - t[0],
                     eps_max=0.02, n_eps=801)
    target = rc.dvv
    print(f"calibration: injected {args.amp:+.4%} reads as dvv {target:+.4%} "
          f"(measured/injected = {target / args.amp:.3f}; origin-offset scale) -> this is TRUTH",
          flush=True)

    # ---- injected pass (step recovery) ----
    pair_inj, daily_inj, info = _run_pass(top, stacks, pat, dates, t, fs, args, args.amp)
    tot_pairs = sum(i[1] for i in info)
    print(f"injected pass: {len(info)} patches, {tot_pairs:,} doublet pairs total, "
          f"median pair cc {np.median([i[2] for i in info if np.isfinite(i[2])]):.3f}", flush=True)

    # ---- natural pass (no injection) for sanity ----
    pair_nat, daily_nat, info_nat = _run_pass(top, stacks, pat, dates, t, fs, args, 0.0)
    print(f"natural pass: median daily cc {np.nanmedian([i[3] for i in info_nat]):.3f}; "
          f"natural pairwise std {100 * pair_nat.std():.4f}%", flush=True)

    # ---- metrics ----
    t0 = np.datetime64(args.t0)
    rows = []
    rows.append(dict(method="pairwise", **_metrics(pair_inj, t0, target)))
    rows.append(dict(method="daily", **_metrics(daily_inj, t0, target)))
    met = pd.DataFrame(rows)
    print(f"\ninjected step reads as {100 * target:+.4f}% at t0={args.t0}\n")
    print(met.to_string(index=False, float_format=lambda x: f"{x:.4f}"))
    csv = f"data/{args.out_prefix}_{args.station}_metrics.csv"
    met.to_csv(csv, index=False)
    # per-day series dump (for differential ctrl/inj comparison across runs)
    pd.DataFrame({"pairwise": pair_inj, "daily": daily_inj}).rename_axis("date").to_csv(
        f"data/{args.out_prefix}_{args.station}_series.csv")

    runtime = time.time() - t_start
    print(f"\nruntime {runtime:.1f} s ({args.workers} workers)", flush=True)

    # ---- figure ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(3, 1, figsize=(13, 12))
    zoom = (t0 - np.timedelta64(180, "D"), t0 + np.timedelta64(240, "D"))
    base = daily_inj[daily_inj.index < t0].median()

    # panel 0: injected step zoom
    ax = axes[0]
    di = daily_inj[(daily_inj.index >= zoom[0]) & (daily_inj.index <= zoom[1])]
    ax.plot(di.index, 100 * di.values, color="0.7", lw=0.8, alpha=0.7, label="daily (vs all-time ref)")
    pi = pair_inj[(pair_inj.index >= zoom[0]) & (pair_inj.index <= zoom[1])]
    ax.plot(pi.index, 100 * pi.values, color="tab:blue", lw=1.7, label="pairwise (doublet inversion)")
    ax.plot([zoom[0], t0, t0, zoom[1]],
            100 * np.array([base, base, base + target, base + target]),
            "k-", lw=2.2, alpha=0.55, label="truth")
    ylo, yhi = ax.get_ylim()
    ax.vlines(t0, ylo, yhi, color="k", ls="--", lw=1)
    ax.set_xlim(zoom)
    ax.set_ylabel("dv/v (%)")
    ax.legend(ncol=3, fontsize=9)
    ax.set_title(f"{args.station}: injected {100 * target:+.4f}% step at {args.t0} "
                 f"(pairwise vs daily) — zoom ±180/240 d")

    # panel 1: full record (injected)
    ax = axes[1]
    ax.plot(daily_inj.index, 100 * daily_inj.values, color="0.8", lw=0.5, alpha=0.6, label="daily")
    ax.plot(pair_inj.index, 100 * pair_inj.values, color="tab:blue", lw=0.9, label="pairwise")
    ax.vlines(t0, *ax.get_ylim(), color="k", ls="--", lw=1)
    ax.set_ylabel("dv/v (%)")
    ax.legend(fontsize=9)
    ax.set_title("full record (injected step)")

    # panel 2: natural (no injection) sanity
    ax = axes[2]
    ax.plot(daily_nat.index, 100 * daily_nat.values, color="0.8", lw=0.5, alpha=0.6, label="daily (natural)")
    ax.plot(pair_nat.index, 100 * pair_nat.values, color="tab:green", lw=0.9, label="pairwise (natural)")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_ylabel("dv/v (%)")
    ax.legend(fontsize=9)
    ax.set_title(f"natural series (NO injection) — pairwise std {100 * pair_nat.std():.4f}% "
                 f"(sanity: flat, like production 2–4 s)")

    fig.tight_layout()
    out = f"figures/{args.out_prefix}_{args.station}.png"
    fig.savefig(out, dpi=140)
    print("->", out)
    print("->", csv)


if __name__ == "__main__":
    main()
