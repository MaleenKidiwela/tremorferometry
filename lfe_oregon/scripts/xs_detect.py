#!/usr/bin/env python
"""Cross-station LFE detector: a line-by-line port of Savard & Bostock (2015, BSSA 105:2128) as coded
in their MATLAB (notes/oregon_sources/CS-quake-detection/matlab_original: main_scc.m, cleanwf.m,
lfedetect.m, mccc2.m, ccsearch.m, gettwin.m, shift.m).  Where code and paper differ the CODE is
followed (see METHOD_FROM_PAPER.md sec.4); every constant below cites its MATLAB line.

Thomas & Bostock (2015) sec.2: "We applied the method of Savard and Bostock (2014) to stations ALP0,
ALS0, and BWO0 ... 3 days of continuous data (09/03-09/05)".

Only deviation forced by that sentence: the code requires >= 4 stations at three places
(main_scc.m:348 `size(newuk,1)<8`, cleanwf.m:163 `ns < 4`, lfedetect.m:76 `<4`, lfedetect.m:237
`nu >= 4`).  With three stations that can never be met, so --min-sta (default 3 = all three) replaces
the 4 at all four places.  Everything else is literal, including two code quirks kept on purpose
(and counted in the day stats):
  Q1 ccsearch never accepts the set containing EVERY channel (ccsearch.m:47);
  Q2 the double-count check (lfedetect.m:151-206) runs BEFORE ccsearch, so a better-envelope window
     that later fails ccsearch still erases the earlier detection.
The one bug fixed: history times are absolute (MATLAB mixes hour-relative `tt` across hours).

Masks: a station whose 18 s window touches an invalid cache sample is treated as absent
(= MATLAB's all-zero/NaN trace, cleanwf.m:21-37).

  python scripts/xs_detect.py --dates 2009-09-03,2009-09-04,2009-09-05 [--min-sta 3] [--cenv 0.43]
-> results/xs/<date>_det.parquet      one row per detection (absolute S times per station)
   results/xs/<date>_stage1.npz       per-window envelope cc of every channel pair (sweep column)
   results/xs/<date>_stats.json
"""
from __future__ import annotations

import argparse
import json
import os
import time
from collections import deque
from itertools import combinations

import numpy as np
import pandas as pd
from scipy.fft import irfft, rfft
from scipy.signal import detrend, hilbert

import xs_common as C

NDT = 0.025                 # main_scc.m:10
TLEN = 18.0                 # main_scc.m:7   window length (paper: 18 s, 17.5 s overlap)
TINC = 0.5                  # main_scc.m:8
CENV = 0.43                 # lfedetect.m:28 ThreshCC (paper: 0.4-0.5)
SIGR = 0.3                  # lfedetect.m:24 sigrthresh (paper sigma_max 0.3 s)
CS_TWIN = 2.0               # lfedetect.m:215 ccsearch twin=2 -> +/-1 s lag (paper says +/-2 s)
DC = 0.01                   # ccsearch.m:22 Cmin decrement
MASK_S = 4.0                # lfedetect.m:137-138 +/-4 s around the envelope-stack maximum
DUP_OVERLAP = 4.0           # lfedetect.m:180
T_FRAC, TF_HIGH, TF_LOW, TF_MAX, T_LOW = 70.0, 6.0, 2.0, 70.0, 60.0   # cleanwf.m:8-12
SAME_STA_TWIN = 0.3         # main_scc.m:357 gettwin(newlst, 0.3, ...)
V_S, MINZ, MAXANG = 3.6, 25.0, np.radians(45.0)                      # gettwin.m:37-41
NW = int(round(TLEN / NDT))                                           # 720 samples
NFFT = 1 << int(np.ceil(np.log2(2 * NW - 1)))                         # mccc2.m:20 -> 2048


# ------------------------------------------------------------------ gettwin.m ----------------
def gettwin(stas):
    """(2ns, 2ns) twin matrix exactly as gettwin(..., 0.3, 'double'): 2 x max lag in s."""
    ns = len(stas)
    xy = np.array([C.utm_km(s) for s in stas])
    tw = np.zeros((ns, ns))
    for i, j in combinations(range(ns), 2):
        D = float(np.hypot(*(xy[i] - xy[j])))
        t = 1.5 * D / V_S
        if D < MINZ * np.tan(MAXANG):
            t = ((MINZ / np.cos(MAXANG)) - np.sqrt((MINZ * np.tan(MAXANG) - D) ** 2 + MINZ ** 2)) / V_S
        tw[i, j] = tw[j, i] = 2 * t
    tw += np.eye(ns) * SAME_STA_TWIN
    return np.tile(tw, (2, 2))


# ------------------------------------------------------------------ mccc2.m ------------------
def mccc2(seis, twin):
    """(tdel, rmean, sigr, r) -- VanDecar & Crosson MCCC exactly as mccc2.m (weight=0).
    seis (n, T); twin scalar or (n, n) (2 x max lag, s)."""
    seis = seis - seis.mean(1, keepdims=True)
    n, T = seis.shape
    nt = 1 << int(np.ceil(np.log2(2 * T - 1)))
    if np.isscalar(twin):
        twin = np.full((n, n), float(twin))
    F = rfft(seis, nt, axis=1)
    sigt = np.sqrt(np.einsum("ij,ij->i", seis, seis))                       # sqrt(max acf)
    r = np.zeros((n, n))
    tcc = np.zeros((n, n))
    iu, ju = np.triu_indices(n, 1)
    if iu.size:
        cc = irfft(np.conj(F[iu]) * F[ju], nt, axis=1)                        # (p, nt)
        itw = np.minimum(nt // 2, np.round(twin[iu, ju] / (2 * NDT)).astype(int))
        lag = np.arange(nt)
        ok = (lag[None, :] <= itw[:, None]) | (lag[None, :] >= nt - itw[:, None])
        ccm = np.where(ok, cc, 0.0)                                           # mask multiplies
        k = ccm.argmax(1)
        cmax = ccm[np.arange(iu.size), k]
        with np.errstate(invalid="ignore", divide="ignore"):
            r[iu, ju] = cmax / (sigt[iu] * sigt[ju])
        kk = np.where(k > nt / 2, k - nt, k)                                  # 0-based lag
        tcc[iu, ju] = kk * NDT
    with np.errstate(invalid="ignore", divide="ignore"):
        z = 0.5 * np.log((1 + r) / (1 - r))
    zmean = (z.sum(1) + z.sum(0)) / max(n - 1, 1)
    rmean = np.tanh(zmean)
    tdel = np.array([(-tcc[:i, i].sum() + tcc[i, i + 1:].sum()) / n for i in range(n)])
    res = np.zeros((n, n))
    res[iu, ju] = tcc[iu, ju] - (tdel[iu] - tdel[ju])
    with np.errstate(invalid="ignore", divide="ignore"):
        sigr = np.sqrt(((res ** 2).sum(1) + (res ** 2).sum(0)) / (n - 2)) if n > 2 else np.full(n, np.inf)
    return tdel, rmean, sigr, r


def fshift(x, t0):
    """shift.m: FFT (sinc) shift of each row of x by t0[i] seconds (positive = delay)."""
    n = x.shape[1]
    om = np.fft.rfftfreq(n, NDT) * 2 * np.pi
    F = np.fft.rfft(x, axis=1) * np.exp(-1j * om[None, :] * np.asarray(t0)[:, None])
    if n % 2 == 0:
        F[:, -1] = np.fft.rfft(x, axis=1)[:, -1] * np.cos(np.pi * np.asarray(t0) / NDT)
    return np.fft.irfft(F, n, axis=1)


def ccsearch(seis, twin, smax, allow_full=False):
    """ccsearch.m (weight=0). -> (ix, tdel, rmean, sigr, cmin_stop).
    allow_full=True is the PAPER's reading ("the largest subset S3"): the full set may be accepted."""
    _, rmean0, sigr0, cc = mccc2(seis, twin)
    cmax = np.nanmax(cc) if cc.size else 0.0
    ix, tdel, rmean, sigr = None, None, None, None
    pflag = True
    while pflag and cmax > 0:
        ir, ic = np.nonzero(cc > cmax)
        ix2 = np.unique(np.concatenate([ir, ic]))
        if ix2.size > 2:
            td2, rm2, sg2, _ = mccc2(seis[ix2], twin)
            if sg2.max() > smax or (ix2.size == seis.shape[0] and not allow_full):   # Q1
                pflag = False
            elif ix2.size == seis.shape[0]:
                ix, tdel, rmean, sigr = ix2, td2, rm2, sg2
                pflag = False
            else:
                ix, tdel, rmean, sigr = ix2, td2, rm2, sg2
        cmax -= DC
    return ix, tdel, rmean, sigr, cmax


# ------------------------------------------------------------------ cleanwf.m (vectorised) --
def cleanwf_batch(N, E, valid):
    """N, E (W, ns, T) hour-normalised windows; valid (W, ns) mask-clean & present.
    -> keep (W, ns) bool, Nc, Ec detrended/demeaned/unit-max (cleanwf.m:57-68)."""
    keep = valid & ((np.abs(N.sum(2)) > 0) | (np.abs(E.sum(2)) > 0))          # cleanwf.m:31
    ne = np.abs(hilbert(N, axis=2))
    ee = np.abs(hilbert(E, axis=2))
    with np.errstate(invalid="ignore", divide="ignore"):
        lown = ne.sum(2) / ne.max(2) < T_LOW                                  # cleanwf.m:50-51
        lowe = ee.sum(2) / ee.max(2) < T_LOW
    keep &= ~(lown | lowe)
    Nc = detrend(N, axis=2)
    Ec = detrend(E, axis=2)
    Nc -= Nc.mean(2, keepdims=True)
    Ec -= Ec.mean(2, keepdims=True)
    am = np.maximum(np.abs(Nc).max(2), np.abs(Ec).max(2))
    am = np.where(am > 0, am, 1.0)[..., None]
    Nc /= am
    Ec /= am
    T = N.shape[2]
    f = np.arange(T // 2 + 1) * C.FS / T
    bad = np.zeros(keep.shape, bool)
    for X in (Nc, Ec):
        S = np.abs(np.fft.rfft(X, axis=2)) ** 2 / (C.FS * T)
        S[..., 1:-1] *= 2
        Sf = S * f
        tot = Sf.sum(2)
        with np.errstate(invalid="ignore", divide="ignore"):
            flo = Sf[..., f < TF_LOW].sum(2) / tot * 100
            fhi = Sf[..., f > TF_HIGH].sum(2) / tot * 100
        I = S.argmax(2)
        cs = np.concatenate([np.zeros(S.shape[:2] + (1,)), np.cumsum(S, 2)], 2)
        lo = np.clip(I - 5, 0, None)
        hi = np.clip(I + 5, None, S.shape[2] - 1)
        pk = (np.take_along_axis(cs, (hi + 1)[..., None], 2) - np.take_along_axis(cs, lo[..., None], 2))[..., 0]
        with np.errstate(invalid="ignore", divide="ignore"):
            areapk = pk / S.sum(2) * 100
        bad |= (flo > T_FRAC) | (fhi > T_FRAC) | (areapk > TF_MAX)
    keep &= ~bad
    return keep, Nc, Ec


def stage1_batch(Nc, Ec, keep, twin, chunk=900):
    """Envelope cc (mccc2 r) for every channel pair of every window; -inf where a channel's
    station is not kept.  -> R (W, P) float32, pairs (P, 2) in [N stations ; E stations] order."""
    W, ns, T = Nc.shape
    env = np.concatenate([np.abs(hilbert(Nc, axis=2)), np.abs(hilbert(Ec, axis=2))], 1)  # (W,2ns,T)
    env -= env.mean(2, keepdims=True)
    pairs = np.array(list(combinations(range(2 * ns), 2)))
    itw = np.minimum(NFFT // 2, np.round(twin[pairs[:, 0], pairs[:, 1]] / (2 * NDT)).astype(int))
    lag = np.arange(NFFT)
    okl = (lag[None, :] <= itw[:, None]) | (lag[None, :] >= NFFT - itw[:, None])
    R = np.empty((W, len(pairs)), np.float32)
    for a in range(0, W, chunk):
        e = env[a:a + chunk]
        F = rfft(e, NFFT, axis=2)
        nrm = np.sqrt((e ** 2).sum(2))
        cc = irfft(np.conj(F[:, pairs[:, 0]]) * F[:, pairs[:, 1]], NFFT, axis=2)
        cc = np.where(okl[None], cc, 0.0).max(2)
        with np.errstate(invalid="ignore", divide="ignore"):
            R[a:a + chunk] = cc / (nrm[:, pairs[:, 0]] * nrm[:, pairs[:, 1]])
    kc = np.concatenate([keep, keep], 1)
    R[~(kc[:, pairs[:, 0]] & kc[:, pairs[:, 1]])] = -np.inf
    R[~np.isfinite(R) & (R != -np.inf)] = -np.inf
    return R, pairs, env


# ------------------------------------------------------------------ lfedetect.m (per window) -
def lfedetect_one(nc, ec, keep_s, twin_full, cenv, min_sta, rpairs, pairs):
    """Stages 2-4 for one window.  nc, ec (ns, T) cleaned.  Returns dict with env-stage results
    (always, if stage 1 passes) and a 'det' sub-dict if the detection passes."""
    ns = nc.shape[0]
    hit = rpairs > cenv
    if not hit.any():
        return None
    ix2 = np.unique(pairs[hit].ravel())
    ian = ix2[ix2 < ns]
    iae = ix2[ix2 >= ns] - ns
    if np.unique(np.concatenate([ian, iae])).size < min_sta:            # lfedetect.m:76
        return None
    ne = np.abs(hilbert(nc, axis=1))
    ee = np.abs(hilbert(ec, axis=1))
    sel = np.concatenate([ian, iae + ns])
    tdel0, _, _, _ = mccc2(np.vstack([ne[ian], ee[iae]]), twin_full[np.ix_(sel, sel)])
    nc, ne, ec, ee = nc.copy(), ne.copy(), ec.copy(), ee.copy()
    if ian.size:
        nc[ian] = fshift(nc[ian], tdel0[:ian.size])
        ne[ian] = fshift(ne[ian], tdel0[:ian.size])
    if iae.size:
        ec[iae] = fshift(ec[iae], tdel0[ian.size:])
        ee[iae] = fshift(ee[iae], tdel0[ian.size:])
    esum = np.vstack([ne[ian], ee[iae]]).mean(0)
    ie = int(esum.argmax())
    maxwin = float(esum[ie])
    tabs = ie * NDT
    i1 = max(0, ie - int(MASK_S / NDT))
    i2 = min(nc.shape[1] - 1, ie + int(MASK_S / NDT))
    out = dict(ian=ian, iae=iae, tdel0=tdel0, maxwin=maxwin, tabs=tabs, i1=i1, i2=i2,
               nc=nc, ec=ec)
    return out


def finish_detect(o, ns, min_sta, sigr_max, allow_full=False, why=None):
    """ccsearch + station times (lfedetect.m:213-296)."""
    ian, iae = o["ian"], o["iae"]
    seis = np.vstack([o["nc"][ian, o["i1"]:o["i2"] + 1], o["ec"][iae, o["i1"]:o["i2"] + 1]])
    iy, tdel, rmean, sigr, cstop = ccsearch(seis, CS_TWIN, sigr_max, allow_full)
    if iy is None:
        if why is not None:
            why.update(reason="ccsearch_none", n_cs=0)
        return None
    iw = np.concatenate([ian, iae + ns])
    ix = iw[iy]
    isn = ix < ns
    ixn, ixe = ix[isn], ix[~isn] - ns
    ixz = np.unique(np.concatenate([ixn, ixe]))
    if ixz.size < min_sta:
        if why is not None:
            why.update(reason="too_few_sta", n_cs=int(iy.size), n_cs_sta=int(ixz.size))
        return None
    sta_rows = []
    both = np.concatenate([ixn, ixe])
    iae_ns = np.concatenate([ian, iae])
    for s in ixz:
        iq = np.flatnonzero(both == s)
        comp = "B" if iq.size == 2 else ("N" if iq[0] < ixn.size else "E")
        tdels = float(tdel[iq].mean())
        ir = iae_ns == s
        ts = o["tabs"] - tdels - float(o["tdel0"][ir].mean())
        wt = int(np.fix(sigr[iq].mean() / 0.1))
        sta_rows.append((int(s), comp, ts, wt, float(sigr[iq].mean()), float(rmean[iq].mean())))
    return dict(sta=sta_rows, n_ch=int(iy.size), cmin_stop=float(cstop + DC),
                sigr_max=float(sigr.max()), rmean_mean=float(rmean.mean()))


# ------------------------------------------------------------------ day driver ---------------
def run_day(date, stas, min_sta, cenv, sigr_max, out_dir, hours=24, allow_full=False):
    t_start = time.time()
    ns = len(stas)
    D = {s: C.load_day(s, date) for s in stas}
    t_day = pd.Timestamp(date).timestamp()
    nday = int(86400 * C.FS)
    twin = gettwin(stas)
    step = int(round(TINC / NDT))
    # continuous day arrays with 18 s of zero/invalid tail (MATLAB: hour 23 is zero padded)
    X = np.zeros((ns, 2, nday + NW), np.float32)
    M = np.zeros((ns, nday + NW), bool)
    for i, s in enumerate(stas):
        if D[s] is None:
            continue
        t0, x, m = D[s]
        assert abs(t0 - t_day) < 1e-6
        X[i, :, :nday] = x[1:3, :nday]
        M[i, :nday] = m[:nday]
    bad_cum = np.concatenate([np.zeros((ns, 1), np.int64), np.cumsum(~M, axis=1)], 1)
    dets, stage1, fails = [], [], []
    hist = deque(maxlen=int(np.ceil(TLEN / TINC)))                  # main_scc.m:394 numMem = 36
    st = dict(windows=0, passed_clean=np.zeros(ns, int), stage1_pass=0, dc_suppressed=0,
              dc_erased_prev=0, ccsearch_fail=0, detections=0, q2_lost=0)
    pend = []                                                         # detections awaiting exit

    def emit(h):
        if h["detect"]:
            dets.append(h["row"])

    for hr in range(hours):
        h0 = hr * 3600 * int(C.FS)
        seg = X[:, :, h0:h0 + int((3600 + TLEN) * C.FS)].astype(np.float64)
        umax = np.abs(seg).reshape(ns, -1).max(1)                     # main_scc.m:285
        umax = np.where(umax > 0, umax, np.finfo(float).eps)
        seg = seg / umax[:, None, None]
        nwin = int(3600 / TINC)
        starts = np.arange(nwin) * step
        idx = starts[:, None] + np.arange(NW)[None, :]
        N = np.transpose(seg[:, 0, :][:, idx], (1, 0, 2))             # (W, ns, T)
        E = np.transpose(seg[:, 1, :][:, idx], (1, 0, 2))
        a = h0 + starts
        valid = (bad_cum[:, a + NW] - bad_cum[:, a] == 0).T           # (W, ns)
        keep, Nc, Ec = cleanwf_batch(N, E, valid)
        st["passed_clean"] += keep.sum(0)
        ok_w = keep.sum(1) >= min_sta                                 # cleanwf.m:163 / main_scc:348
        R, pairs, _ = stage1_batch(Nc, Ec, keep, twin)
        stage1.append(R.astype(np.float16))
        st["windows"] += nwin
        for w in range(nwin):
            t_w = t_day + (h0 + starts[w]) / C.FS
            rec = dict(detect=0, i1s=0.0, i2s=0.0, maxwin=0.0, row=None)
            if ok_w[w]:
                ks = np.flatnonzero(keep[w])
                # restrict to kept stations (cleanwf removes rows); pair indices remapped
                sub = np.concatenate([ks, ks + ns])
                pos = -np.ones(2 * ns, int)
                pos[sub] = np.arange(sub.size)
                pm = (pos[pairs[:, 0]] >= 0) & (pos[pairs[:, 1]] >= 0)
                o = lfedetect_one(Nc[w, ks], Ec[w, ks], None, twin[np.ix_(sub, sub)], cenv, min_sta,
                                  R[w, pm].astype(float), np.stack([pos[pairs[pm, 0]], pos[pairs[pm, 1]]], 1))
                if o is not None:
                    st["stage1_pass"] += 1
                    i1s = t_w + o["i1"] * NDT                          # lfedetect.m:143 (absolute)
                    i2s = t_w + o["i2"] * NDT
                    rec.update(i1s=i1s, i2s=i2s, maxwin=o["maxwin"])
                    blocked = False
                    erased = []
                    for h in hist:
                        if not h["detect"]:
                            continue
                        ov = min(i2s, h["i2s"]) - max(i1s, h["i1s"])
                        if ov > DUP_OVERLAP:
                            if h["maxwin"] > o["maxwin"]:
                                blocked = True
                                break
                            erased.append(h)
                    if blocked:
                        st["dc_suppressed"] += 1
                    else:
                        for h in erased:
                            h["detect"] = 0
                        st["dc_erased_prev"] += len(erased)
                        why = dict(t_win=t_w, n_env_ch=int(o["ian"].size + o["iae"].size),
                                   n_env_sta=int(np.unique(np.concatenate([o["ian"], o["iae"]])).size),
                                   maxwin=o["maxwin"])
                        fd = finish_detect(o, len(ks), min_sta, sigr_max, allow_full, why)
                        if fd is None:
                            fails.append(why)
                            st["ccsearch_fail"] += 1
                            st["q2_lost"] += len(erased)
                        else:
                            row = dict(t_win=t_w, t_ref=t_w + o["tabs"], maxwin=o["maxwin"],
                                       n_env_ch=int(o["ian"].size + o["iae"].size),
                                       n_ch=fd["n_ch"], n_sta=len(fd["sta"]), cmin_stop=fd["cmin_stop"],
                                       sigr_max=fd["sigr_max"], rmean=fd["rmean_mean"],
                                       env_cc_max=float(R[w, pm][np.isfinite(R[w, pm])].max()))
                            for (si, comp, ts, wt, sg, rm) in fd["sta"]:
                                sname = stas[ks[si]]
                                row[f"tS_{sname}"] = t_w + ts
                                row[f"comp_{sname}"] = comp
                                row[f"wt_{sname}"] = wt
                                row[f"sigr_{sname}"] = sg
                            rec.update(detect=1, row=row)
            if len(hist) == hist.maxlen:
                emit(hist[0])                                          # main_scc.m:412 exiting window
            hist.append(rec)
    for h in hist:
        emit(h)
    Dt = pd.DataFrame(dets)
    st["detections"] = int(len(Dt))
    st["passed_clean"] = dict(zip(stas, st["passed_clean"].tolist()))
    st["secs"] = round(time.time() - t_start, 1)
    st.update(date=date, stations=stas, min_sta=min_sta, cenv=cenv, sigr_max=sigr_max, allow_full=allow_full)
    os.makedirs(out_dir, exist_ok=True)
    Dt.to_parquet(f"{out_dir}/{date}_det.parquet", index=False)
    pd.DataFrame(fails).to_parquet(f"{out_dir}/{date}_fails.parquet", index=False)
    np.savez_compressed(f"{out_dir}/{date}_stage1.npz", R=np.concatenate(stage1), pairs=pairs,
                        t0=t_day, tinc=TINC, stations=np.array(stas))
    json.dump(st, open(f"{out_dir}/{date}_stats.json", "w"), indent=1, default=int)
    return Dt, st


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dates", default=",".join(C.DETECT_DAYS))
    ap.add_argument("--stations", default=",".join(C.XS_STATIONS))
    ap.add_argument("--min-sta", type=int, default=3)
    ap.add_argument("--cenv", type=float, default=CENV)
    ap.add_argument("--sigr", type=float, default=SIGR)
    ap.add_argument("--allow-full-set", action="store_true", help="paper reading of S3 (ablation of code quirk Q1)")
    ap.add_argument("--hours", type=int, default=24, help="testing: only the first N hours")
    ap.add_argument("--out", default=f"{C.RES}/xs")
    a = ap.parse_args()
    for d in a.dates.split(","):
        Dt, st = run_day(d, a.stations.split(","), a.min_sta, a.cenv, a.sigr, a.out, a.hours, a.allow_full_set)
        print(json.dumps({k: v for k, v in st.items() if k != "stations"}, default=int), flush=True)


if __name__ == "__main__":
    main()
