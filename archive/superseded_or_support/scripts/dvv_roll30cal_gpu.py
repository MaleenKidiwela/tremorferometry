#!/usr/bin/env python
"""GPU dv/v: identical 30-cal-day rolling-stack + SVD-Wiener prep as dvv_roll30cal.py, but the expensive
per-day stretch (ne days x n_eps stretches x cubic resample + correlation) is GPU-batched with cupy.
Cubic spline = not-a-knot (matches scipy CubicSpline, extrapolate=False -> 0). Same CLI + output columns.
  python dvv_roll30cal_gpu.py --station B040 --npz ... --window 2 4 --origin-anchor --eps-max 0.04 --n-eps 401 --out ...
  add --validate N to compare N families against the CPU stretch_dvv.
"""
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, sys, time
import numpy as np, pandas as pd
sys.path.insert(0, "src")
import cupy as cp

# ---------------- family prep (rolling stack + SVD-Wiener) : identical to CPU dvv_roll30cal._family ----------------
def prep_family(S, dstr, days, minst, mindays, prenoise_mask):
    if len(S) < mindays:
        return None
    order = np.argsort(dstr); S = S[order].astype(np.float64); dstr = np.asarray(dstr)[order]
    Dn = pd.to_datetime(dstr).values.astype("datetime64[D]")
    n = len(S); idx = np.arange(n)
    cs = np.vstack([np.zeros((1, S.shape[1])), np.cumsum(S, axis=0)])
    los = np.searchsorted(Dn, Dn - np.timedelta64(days, "D"), side="right")
    counts = idx + 1 - los; keep = counts >= minst
    if int(keep.sum()) < 8:
        return None
    R = (cs[idx + 1] - cs[los]) / counts[:, None]
    R = R[keep]; Dr = dstr[keep]; ck = counts[keep]
    pk = np.max(np.abs(R), axis=1, keepdims=True); pk[pk == 0] = 1.0; Rn = R / pk
    ne, ns = Rn.shape; sig = Rn[:, prenoise_mask].std()
    try:
        U, sv, Vt = np.linalg.svd(Rn, full_matrices=False)
        nmf = sig * (np.sqrt(ne) + np.sqrt(ns)); w = sv**2 / (sv**2 + nmf**2); Rf = (U * (sv * w)) @ Vt
    except np.linalg.LinAlgError:
        Rf = Rn
    return Rf, Dr, ck

# ---------------- GPU cubic-spline stretch prep (fixed per station: depends on ns,fs,window,eps grid) ---------
def stretch_prep(ns, fs, tmin, tmax, eps_max, n_eps):
    h = 1.0 / fs; t = np.arange(ns) / fs
    mask = (t >= tmin) & (t < tmax)
    if mask.sum() < 8:
        raise ValueError("coda window < 8 samples")
    tc = t[mask]; eps = np.linspace(-eps_max, eps_max, n_eps)
    Q = tc[None, :] * (1.0 + eps[:, None])                       # (n_eps, ncoda) query times
    k = np.clip((Q / h).astype(np.int64), 0, ns - 2)            # interval index
    xk = k * h; a = (xk + h - Q) / h; b = (Q - xk) / h
    oob = (Q < 0) | (Q > (ns - 1) * h)
    ac = (a**3 - a) * h * h / 6.0; bc = (b**3 - b) * h * h / 6.0
    A = np.zeros((ns, ns))                                       # not-a-knot second-derivative system (uniform h)
    A[0, 0], A[0, 1], A[0, 2] = 1.0, -2.0, 1.0
    A[-1, -3], A[-1, -2], A[-1, -1] = 1.0, -2.0, 1.0
    ii = np.arange(1, ns - 1); A[ii, ii - 1] = 1.0; A[ii, ii] = 4.0; A[ii, ii + 1] = 1.0
    return dict(h=h, mask=cp.asarray(mask), Ainv=cp.asarray(np.linalg.inv(A)),
                k=cp.asarray(k), kp=cp.asarray(k + 1), a=cp.asarray(a), b=cp.asarray(b),
                ac=cp.asarray(ac), bc=cp.asarray(bc), oob=cp.asarray(oob),
                eps=cp.asarray(eps), eps_np=eps, ncoda=int(mask.sum()), n_eps=n_eps, ns=ns)

def gpu_cc(Rf, i0, P, day_batch=256):
    """Return cc (ne, n_eps) on CPU: normalized coda correlation of each day's stretched stack vs the mean ref."""
    cur_full = cp.asarray(Rf[:, i0:])                            # (ne, ns)
    ne, ns = cur_full.shape; h = P["h"]
    ref = cur_full.mean(0)
    rc = ref[P["mask"]]; rc = rc - rc.mean(); rn = cp.sqrt((rc * rc).sum())
    if float(rn) == 0:
        return None
    out = np.empty((ne, P["n_eps"]), np.float64)
    k, kp, a, b, ac, bc, oob = P["k"], P["kp"], P["a"], P["b"], P["ac"], P["bc"], P["oob"]
    for s in range(0, ne, day_batch):
        cur = cur_full[s:s + day_batch]                          # (nb, ns)
        RHS = cp.zeros_like(cur); RHS[:, 1:-1] = 6.0 / (h * h) * (cur[:, :-2] - 2 * cur[:, 1:-1] + cur[:, 2:])
        M = RHS @ P["Ainv"].T                                    # (nb, ns) second derivatives
        S = a * cur[:, k] + b * cur[:, kp] + ac * M[:, k] + bc * M[:, kp]   # (nb, n_eps, ncoda)
        S = cp.where(oob, 0.0, S)
        Sd = S - S.mean(-1, keepdims=True)
        num = (Sd * rc).sum(-1); den = rn * cp.sqrt((Sd * Sd).sum(-1))
        out[s:s + day_batch] = cp.asnumpy(cp.where(den > 0, num / den, 0.0))
    return out

def refine(cc, eps):                                            # per-day argmax + parabolic (matches _parabolic_refine)
    n_eps = len(eps); step = eps[1] - eps[0]
    im = cc.argmax(1); dvv = np.empty(len(cc)); ccm = cc[np.arange(len(cc)), im]
    for r in range(len(cc)):
        i = im[r]
        if i == 0 or i == n_eps - 1: dvv[r] = eps[i]; continue
        y0, y1, y2 = cc[r, i - 1], cc[r, i], cc[r, i + 1]; den = y0 - 2 * y1 + y2
        dvv[r] = eps[i] + (0.5 * (y0 - y2) / den * step if den != 0 else 0.0)
    return dvv, ccm

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--station", required=True); p.add_argument("--npz", required=True)
    p.add_argument("--window", nargs=2, type=float, default=[2.0, 4.0])
    p.add_argument("--days", type=int, default=30); p.add_argument("--min-stacks", type=int, default=5)
    p.add_argument("--mindays", type=int, default=60); p.add_argument("--n-eps", type=int, default=201)
    p.add_argument("--eps-max", type=float, default=0.02); p.add_argument("--out", required=True)
    p.add_argument("--origin-anchor", action="store_true"); p.add_argument("--day-batch", type=int, default=256)
    p.add_argument("--validate", type=int, default=0, help="compare N families vs CPU stretch_dvv and exit")
    a = p.parse_args(); t_s = 1.0
    d = np.load(a.npz, allow_pickle=True)
    t = d["t"]; fs = float(d["fs"]); stacks = d["stacks"]; pat = d["patches"].astype(str)
    dates = pd.to_datetime(d["dates"]).strftime("%Y-%m-%d").values
    i0 = int(np.searchsorted(t, t_s)) if a.origin_anchor else 0
    tmin = (a.window[0] - t_s) if a.origin_anchor else (a.window[0] - t[0])
    tmax = (a.window[1] - t_s) if a.origin_anchor else (a.window[1] - t[0])
    prenoise = (t < -1.0)
    ns = len(t) - i0
    P = stretch_prep(ns, fs, tmin, tmax, a.eps_max, a.n_eps)
    uniq = pd.unique(pat)
    print(f"[{a.station}] {len(stacks)} stacks, {len(uniq)} families | ns={ns} ncoda={P['ncoda']} n_eps={a.n_eps} eps_max={a.eps_max} | GPU", flush=True)
    if a.validate:
        from tremorferometry.dvv import stretch_dvv
        errs = []
        for fam in uniq[:a.validate]:
            r = prep_family(stacks[pat == fam], dates[pat == fam], a.days, a.min_stacks, a.mindays, prenoise)
            if r is None: continue
            Rf, Dr, ck = r; cc = gpu_cc(Rf, i0, P, a.day_batch)
            if cc is None: continue
            dvv_g, ccm_g = refine(cc, P["eps_np"]); ref = Rf[:, i0:].mean(0)
            for i in range(min(30, len(Rf))):
                rc = stretch_dvv(ref, Rf[i, i0:], fs=fs, t_min=tmin, t_max=tmax, eps_max=a.eps_max, n_eps=a.n_eps)
                errs.append((abs(rc.dvv - dvv_g[i]), abs(rc.cc_max - ccm_g[i])))
        e = np.array(errs)
        print(f"VALIDATE {len(e)} day-samples: max|ddvv| {e[:,0].max():.2e} median {np.median(e[:,0]):.2e} | max|dcc| {e[:,1].max():.2e}")
        return
    t0 = time.time(); out = []
    for fam in uniq:
        r = prep_family(stacks[pat == fam], dates[pat == fam], a.days, a.min_stacks, a.mindays, prenoise)
        if r is None: continue
        Rf, Dr, ck = r; cc = gpu_cc(Rf, i0, P, a.day_batch)
        if cc is None: continue
        dvv, ccm = refine(cc, P["eps_np"])
        out.extend([(fam, Dr[i], float(dvv[i]), float(ccm[i]), int(ck[i])) for i in range(len(Rf))])
    df = pd.DataFrame(out, columns=["patch", "date", "dvv", "cc_max", "n_stack"]); df.to_csv(a.out, index=False)
    print(f"[{a.station}] DONE {len(df):,} rows, {df.patch.nunique()} families, mean cc {df.cc_max.mean():.3f}, "
          f"{time.time()-t0:.0f}s -> {a.out}", flush=True)

if __name__ == "__main__":
    main()
