#!/usr/bin/env python
"""STAGE 5 (FAST) -- ab_grow.py with the CPU side of the matched filter moved off the critical path.

SAME METHOD, SAME CLI, SAME OUTPUTS as scripts/ab_grow.py (which is left untouched so the provenance
of runs made with it stays exact).  Every scientific function (cut_templates, merge_duplicates,
load_stacks, unit_max3, _sliding_max, the GPU correlation / network-function / 8 x MAD / decluster
arithmetic) is ab_grow's own, imported or copied verbatim.  Only the scheduling differs:

  * profile (2026-09-26, NW 2008-05-08, 2,695 tpl_it0 templates, 15 stations): ~22 s/day, of which
    the serial CPU day preparation -- DayCache load 0.8 s, ab_common.day_refs 1.6 s, the 15
    window_stats_guarded calls 2.3 s, float64 casts/uploads -- stalls the GPU;
  * here a pool of --prep-workers SPAWNED processes (default 6, <= 8) prepares the next days ahead
    (load the station-days, day_refs over the whole loaded pool, window_stats_guarded per template
    station, and the float32 denominator DEN = where(okc, sqrt(max(e, 0)), inf) -- the very
    expression ab_grow evaluates on the GPU, done in float64 then cast, so it is bit-identical)
    while the GPU correlates the current day.  okc is recovered on the GPU as isfinite(DEN);
  * the restack runs over date chunks in --restack-workers processes via
    ab_grow_cont.restack_parallel (verified: identical sums up to float summation order).

Equivalence test: scripts/ab_grow_fast_equiv.py (3 NW days x ~300 templates vs ab_grow.matched_filter).

  python scripts/ab_grow_fast.py --sub NW --pool ... --out ... --max-iters 14 --tol 0.01

--merge fair (2026-09-28, user-approved): see merge_fair / site_scores below.  Default --merge most
is the old merge_duplicates, unchanged.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import os
import sys
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402
import ab_grow as G                                                         # noqa: E402
from ab_grow import (LONG_PRE, LONG_LEN, TPL_LEN, MAX_STA, MIN_STA_VALID,  # noqa: E402,F401
                     cut_templates, merge_duplicates, load_stacks, unit_max3, _sliding_max)

PREP_WORKERS = 6
PREP_THREADS = 3          # threads per prep worker (numpy releases the GIL in the heavy ops)
RESTACK_WORKERS = 8
POW2_NORM_FOLD = True     # see the template-FFT block
COUNT_GEMM = True         # see the count-path block in matched_filter


# ------------------------------------------------------------------ day preparation (CPU) -------
def _day_refs(D, n, threads, net_rel=A.GUARD_NET_REL):
    """ab_common.day_refs verbatim, with a bounded thread count."""
    with ThreadPoolExecutor(max(min(len(D), threads), 1)) as ex:
        R = dict(zip(D, ex.map(lambda s: A.channel_refs(D[s][1], D[s][2], n), D)))
    tops = np.array([r.max() for r in R.values() if r.max() > 0])
    if tops.size:
        floor = net_rel * float(np.median(tops))
        for s in R:
            if R[s].max() < floor:
                R[s] = np.zeros_like(R[s])
    return R


def prepare_day(args):
    """Everything ab_grow.matched_filter computes on the CPU for one date, in a worker process.

    -> dict(date, have, N, T0 {sta: t0}, X {sta: (3, N) float32}, DEN {sta: (3, NV) float32})
       or dict(date, have, skip=True).  DEN is +inf where the guarded channel-window is invalid."""
    date, pool, stas, nt, min_sta_valid, threads = args
    dc = A.DayCache(pool)
    D = dc.day(date)
    have = [s for s in stas if s in D]
    if len(have) < min_sta_valid:
        return dict(date=date, have=have, skip=True)
    N = min(D[s][1].shape[1] for s in have)
    REF = _day_refs(D, nt, threads)                           # whole loaded pool, network floor

    def one(s):
        ew, okw = A.window_stats_guarded(D[s][1][:, :N], D[s][2][:N], nt, 1, ref=REF[s])
        # == ab_grow: torch.where(okc, e.clamp_min(0).sqrt(), inf).float()  (float64 sqrt, cast)
        den = np.where(okw, np.sqrt(np.maximum(ew, 0.0)), np.inf).astype(np.float32)
        return den

    with ThreadPoolExecutor(max(min(len(have), threads), 1)) as ex:
        DEN = dict(zip(have, ex.map(one, have)))
    X = {s: np.ascontiguousarray(D[s][1][:, :N]) for s in have}          # float32, as cached
    T0 = {s: float(D[s][0]) for s in have}
    return dict(date=date, have=have, N=N, T0=T0, X=X, DEN=DEN, skip=False)


class DayPrefetcher:
    """Ordered iterator over prepare_day results, keeping `depth` dates in flight."""

    def __init__(self, dates, pool, stas, nt, min_sta_valid, workers=PREP_WORKERS,
                 threads=PREP_THREADS, depth=None):
        self.dates = list(dates)
        self.args = [(d, list(pool), list(stas), nt, min_sta_valid, threads) for d in self.dates]
        self.workers = max(1, min(int(workers), 8))
        self.depth = depth or self.workers + 1
        # start at construction, so the first days are prepared while the caller builds templates
        self.ex = ProcessPoolExecutor(self.workers, mp_context=mp.get_context("spawn"))
        self.futs, self.nxt = [], 0
        while self.nxt < len(self.args) and len(self.futs) < self.depth:
            self.futs.append(self.ex.submit(prepare_day, self.args[self.nxt]))
            self.nxt += 1

    def __iter__(self):
        futs = self.futs
        nxt = self.nxt
        try:
            while futs:
                r = futs.pop(0).result()
                if nxt < len(self.args):
                    futs.append(self.ex.submit(prepare_day, self.args[nxt]))
                    nxt += 1
                yield r
        finally:
            self.ex.shutdown(wait=True, cancel_futures=True)


# ------------------------------------------------------------------ matched filter (GPU) --------
def matched_filter(dc, TPL, dates, nmad=A.NMAD, net="mean", min_sta_valid=MIN_STA_VALID,
                   decl_s=TPL_LEN, tpl_len=TPL_LEN, dev=None, L=8192, gchunk=256, fbatch=64,
                   _dump=None, prep_workers=PREP_WORKERS, prep_threads=PREP_THREADS, _prof=None):
    """ab_grow.matched_filter with the per-date CPU preparation prefetched by DayPrefetcher.
    `dc` only supplies the station pool (dc.stas).  The GPU arithmetic below is ab_grow's,
    line for line, from the per-station upload onwards."""
    import time
    import torch
    dev = dev or A.torch_dev()
    torch.backends.cuda.matmul.allow_tf32 = False
    nt = int(round(tpl_len * A.FS))
    S = L - nt + 1
    assert S > 0
    gs_all = sorted(TPL)
    G_all = len(gs_all)
    stas = sorted({s for g in TPL for s in TPL[g]})
    PF = DayPrefetcher(dates, dc.stas, stas, nt, min_sta_valid, prep_workers, prep_threads)
    # ---- template FFTs, once per call (verbatim ab_grow) ------------------------------------
    TS = {}
    for s in stas:
        gi = [i for i, g in enumerate(gs_all) if s in TPL[g]]
        M = np.zeros((len(gi), 3, L), np.float64)
        nrm = np.zeros((len(gi), 3), np.float64)
        ks = np.zeros(len(gi), np.int64)
        for r, i in enumerate(gi):
            tpl, off, _ = TPL[gs_all[i]][s]
            M[r, :, :nt] = tpl
            nrm[r] = np.sqrt((tpl ** 2).sum(1))
            ks[r] = int(round(off * A.FS))
        inv = np.where(nrm > 0, 1.0 / np.where(nrm > 0, nrm, 1.0), 0.0)
        Tf = (torch.fft.rfft(torch.from_numpy(M).to(dev), n=L).conj()
              * torch.from_numpy(inv).to(dev)[:, :, None]).to(torch.complex64)
        # FAST: irfft's 1/L normalisation folded into the template spectrum.  L is a power of 2, so
        # x 2^-log2(L) is EXACT in floating point and commutes exactly with every multiply/add of
        # the product and the c2r transform -> bit-identical cc, one full-size pass fewer per batch
        # (irfft is then called with norm="forward", i.e. unscaled inverse).
        fold = bool(POW2_NORM_FOLD and (L & (L - 1)) == 0)
        if fold:
            Tf = Tf * (1.0 / L)
        # FAST: index tensors per (chunk, batch, offset) built once per call and kept on the GPU
        # (were host->device copies inside the date loop); sel=None when a batch's rows all share
        # one offset (always true for cut="common") -> no gather copy before index_add_.
        gi_np = np.asarray(gi, np.int64)
        IDX = {}
        for c0 in range(0, G_all, gchunk):
            c1 = min(c0 + gchunk, G_all)
            a0, a1 = np.searchsorted(gi_np, [c0, c1])
            for b0 in range(a0, a1, fbatch):
                b1 = min(b0 + fbatch, a1)
                loc = gi_np[b0:b1] - c0
                kb = ks[b0:b1]
                lst = []
                for k in np.unique(kb):
                    sel = np.flatnonzero(kb == k)
                    lst.append((int(k), torch.from_numpy(loc[sel]).to(dev),
                                None if sel.size == b1 - b0 else torch.from_numpy(sel).to(dev)))
                IDX[(c0, b0)] = lst
        TS[s] = dict(gi=gi_np, Tf=Tf, IDX=IDX, norm="forward" if fold else "backward",
                     live=torch.tensor(nrm > 0, device=dev, dtype=torch.float32), ks=ks)
        del M
    # FAST (count path): the valid-channel and valid-station counts are small INTEGERS, so they can
    # be formed by one exact GEMM per chunk instead of per-batch matmul + two index_add_ passes:
    #   cch[g]  = sum_s sum_c live[g,s,c] OK[s,c]            (ab_grow: live @ OKf, index_add_)
    #   csta[g] = sum_s [ live[g,s,.] . OK[s,.] > 0 ]         (ab_grow: (nch > 0), index_add_)
    # with the station indicator taken per live-component PATTERN p (bitmask over Z,N,E):
    #   Ind[s,p] = any_c in p OK[s,c].  Integer sums < 2^24 are exact in float32 in any order.
    # Valid only when every template offset is 0 (cut="common"); otherwise the ref path is used.
    cnt_fast = bool(COUNT_GEMM) and all((TS[s_]["ks"] == 0).all() for s_ in stas)
    if cnt_fast:
        si = {s_: j for j, s_ in enumerate(stas)}
        def _pat(tpl):
            return int(sum(int(v) << c for c, v in enumerate(np.sqrt((tpl ** 2).sum(1)) > 0)))
        pats = {s_: sorted({_pat(TPL[g][s_][0]) for g in gs_all if s_ in TPL[g]} - {0})
                for s_ in stas}
        prow, r_ = {}, 0
        for s_ in stas:
            for p_ in pats[s_]:
                prow[(s_, p_)] = r_
                r_ += 1
        NPR = r_
        CM = {}
        for c0 in range(0, G_all, gchunk):
            gs_c = gs_all[c0:c0 + gchunk]
            ml = np.zeros((len(gs_c), 3 * len(stas)), np.float32)
            mp_ = np.zeros((len(gs_c), NPR), np.float32)
            for i, g in enumerate(gs_c):
                for s_ in TPL[g]:
                    lv = np.sqrt((TPL[g][s_][0] ** 2).sum(1)) > 0          # == TS live (nrm > 0)
                    ml[i, 3 * si[s_]:3 * si[s_] + 3] = lv
                    p_ = int(sum(int(v) << c for c, v in enumerate(lv)))
                    if p_:
                        mp_[i, prow[(s_, p_)]] = 1.0
            CM[c0] = (torch.from_numpy(ml).to(dev), torch.from_numpy(mp_).to(dev))
    out = []
    tw = time.time()
    for P_ in PF:
        date, have = P_["date"], P_["have"]
        t_got = time.time()
        if P_["skip"]:
            A.log(f"[mf] {date} only {len(have)} stations -> skipped")
            continue
        N = P_["N"]
        NV = N - nt + 1
        B = -(-NV // S)
        NP = B * S + nt - 1
        Xf, IDEN, OKf, T0 = {}, {}, {}, {}
        for s in have:
            X = P_["X"].pop(s)
            Xp = torch.zeros((3, NP), dtype=torch.float64, device=dev)
            Xp[:, :N] = torch.from_numpy(X).to(dev).double()      # == float64 cast, then upload
            Xf[s] = torch.fft.rfft(Xp.unfold(-1, L, S), n=L).to(torch.complex64)   # (3,B,F)
            den = torch.from_numpy(P_["DEN"].pop(s)).to(dev)
            okc = torch.isfinite(den)
            iden = torch.zeros((3, B * S), device=dev)
            iden[:, :NV] = torch.where(okc, 1.0 / den, torch.zeros_like(den))  # 0 where invalid
            IDEN[s] = iden.view(3, B, S)
            OKf[s] = okc.float()
            T0[s] = P_["T0"][s]
            del Xp, den, okc, X
        if cnt_fast:
            OKall = torch.zeros((3 * len(stas), NV), device=dev)
            IND = torch.zeros((NPR, NV), device=dev)
            for s_ in have:
                j = si[s_]
                OKall[3 * j:3 * j + 3] = OKf[s_]
                okb = OKf[s_] > 0
                for p_ in pats[s_]:
                    m_ = torch.zeros_like(okb[0])
                    for c in range(3):
                        if p_ >> c & 1:
                            m_ |= okb[c]
                    IND[prow[(s_, p_)]] = m_.float()
        del P_
        ndet = 0
        for c0 in range(0, G_all, gchunk):
            c1 = min(c0 + gchunk, G_all)
            gs = gs_all[c0:c1]
            G = len(gs)
            acc = torch.zeros((G, NV), device=dev)
            if cnt_fast:
                cch = CM[c0][0] @ OKall                  # exact integer GEMM (allow_tf32 False)
                csta = CM[c0][1] @ IND
            else:
                cch = torch.zeros((G, NV), device=dev)
                csta = torch.zeros((G, NV), device=dev)
            for s in have:
                T = TS[s]
                a0, a1 = np.searchsorted(T["gi"], [c0, c1])
                for b0 in range(a0, a1, fbatch):
                    b1 = min(b0 + fbatch, a1)
                    P = Xf[s][None] * T["Tf"][b0:b1, :, None, :]            # (r,3,B,F)
                    cc = torch.fft.irfft(P, n=L, norm=T["norm"])[..., :S]    # (r,3,B,S) cc/nrm
                    del P
                    w = IDEN[s]
                    ccs = cc[:, 0] * w[0]
                    ccs.addcmul_(cc[:, 1], w[1])
                    ccs.addcmul_(cc[:, 2], w[2])
                    del cc
                    ccs = ccs.reshape(b1 - b0, B * S)[:, :NV]
                    nch = None if cnt_fast else T["live"][b0:b1] @ OKf[s]
                    for k, rs, ss in T["IDX"][(c0, b0)]:
                        lo, hi = max(-k, 0), min(NV - k, NV)
                        if hi <= lo:
                            continue
                        if ss is None:                    # all rows: views, no gather copy
                            src_c = ccs[:, lo + k:hi + k]
                        else:
                            src_c = ccs[ss, lo + k:hi + k]
                        acc[:, lo:hi].index_add_(0, rs, src_c)
                        del src_c
                        if cnt_fast:
                            continue
                        src_n = nch[:, lo + k:hi + k] if ss is None else nch[ss, lo + k:hi + k]
                        cch[:, lo:hi].index_add_(0, rs, src_n)
                        csta[:, lo:hi].index_add_(0, rs, (src_n > 0).float())
                        del src_n
                    del ccs, nch
            enough = csta >= float(min_sta_valid)
            fn = acc if net == "sum" else acc / cch.clamp_min(1.0)
            fn = torch.where(enough, fn, torch.full_like(fn, float("nan")))
            if _dump is not None:
                for i, g in enumerate(gs):
                    _dump[(date, g)] = fn[i].cpu().numpy()
            med = torch.nanmedian(fn, dim=1, keepdim=True).values
            madv = torch.nanmedian((fn - med).abs(), dim=1, keepdim=True).values
            thr = med + nmad * madv.clamp_min(1e-9)
            dec = int(round(decl_s * A.FS))
            hit = torch.nan_to_num(fn, nan=-1e9) >= thr
            vals = torch.where(hit, fn, torch.full_like(fn, -1e9))
            pooled = _sliding_max(vals, dec)
            keep = hit & (vals >= pooled - 1e-9)
            ij = torch.nonzero(keep, as_tuple=False)
            if ij.numel():
                gi = ij[:, 0].cpu().numpy()
                ti = ij[:, 1].cpu().numpy()
                v = fn[ij[:, 0], ij[:, 1]].cpu().numpy()
                nm = ((fn[ij[:, 0], ij[:, 1]] - med[ij[:, 0], 0])
                      / madv[ij[:, 0], 0].clamp_min(1e-9)).cpu().numpy()
                ns = csta[ij[:, 0], ij[:, 1]].cpu().numpy().astype(np.int16)
                nc = cch[ij[:, 0], ij[:, 1]].cpu().numpy().astype(np.int16)
                ccm = (acc[ij[:, 0], ij[:, 1]] / cch[ij[:, 0], ij[:, 1]].clamp_min(1.0)).cpu().numpy()
                out.append(pd.DataFrame(dict(grp=[gs[i] for i in gi], t=T0[have[0]] + ti / A.FS,
                                             cc=v, cc_mean=ccm, nmad=nm, n_sta=ns, n_ch=nc,
                                             thr=thr[ij[:, 0], 0].cpu().numpy(), date=date)))
                ndet += len(out[-1])
            del acc, cch, csta, fn, hit, vals, pooled, keep, ij
        t_end = time.time()
        if _prof is not None:
            _prof.append(dict(date=date, wait_s=t_got - tw, gpu_s=t_end - t_got))
        tw = t_end
        A.log(f"[mf] {date} {len(have)} sta -> {ndet}")
        del Xf, IDEN, OKf
        if cnt_fast:
            del OKall, IND
        if dev == "cuda":
            torch.cuda.empty_cache()
    del TS
    return (pd.concat(out, ignore_index=True) if out
            else pd.DataFrame(columns=["grp", "t", "cc", "cc_mean", "nmad", "n_sta", "n_ch",
                                       "thr", "date"]))


def restack(dc, stas, det, max_stack=0, long_pre=LONG_PRE, long_len=LONG_LEN,
            workers=RESTACK_WORKERS):
    """ab_grow.restack over date chunks in parallel (ab_grow_cont.restack_parallel)."""
    import ab_grow_cont as GC
    if max_stack:
        det = det.groupby("grp", group_keys=False).apply(lambda g: g.nlargest(max_stack, "cc"))
    if not len(det):
        return {}
    return GC.restack_parallel(stas, det, workers=max(1, min(int(workers), 8)),
                               long_pre=long_pre, long_len=long_len)


# ------------------------------------------------------------------ FAIR merge (2026-09-28) -----
# User-approved replacement for merge_duplicates' representative choice (Merlin 09-28: "most
# detections" = noise catch-alls swallow real seeds).  Selected with --merge fair; default unchanged.
#  (a) duplicate test = scripts/ab_final_dedupe.py: for every detection of a, the b-detections within
#      +-FAIR_WIN s give lags (FAIR_BIN bins); modal lag = the most populated bin after a +-1-bin box
#      (+-0.75 s); fraction = count in that box / min(n_a, n_b).  Control: the same with b shifted by
#      +-FAIR_SHIFT s (lags in 45..75 s, both signs, the larger kept).  Duplicate if excess > FAIR_THR.
#  (b) representative = the higher site-count score (ab_famscore logic on the member's CURRENT stack
#      band-passed 2-8 Hz; sites = stations < 3 km apart, ab_mute.sites); tie -> more detections.
#  (c) greedy STAR: in descending (score, n_det) a template joins the representative it shares the
#      largest excess with, else becomes a representative (no chaining).
#  (d) a pair with one member scoring >= FAIR_REAL sites and the other < FAIR_REAL is never merged
#      (kept separate, logged in fair_straddle_it<N>.csv).
#  Union of detections: member times aligned onto the representative by the modal lag, then one
#  detection per TPL_LEN s per family (strongest cc) = the matched filter's own declustering.
FAIR_WIN, FAIR_SHIFT, FAIR_BIN, FAIR_NB, FAIR_THR, FAIR_REAL = 15.0, 60.0, 0.5, 1, 0.10, 3


def _famscore_one(args):
    """ab_famscore.score for one family, on a dict {sta: (3, n) 2-8 Hz stack}; site = {sta: site id}."""
    from scipy.signal import hilbert
    g, st, site = args
    good = []
    for s, x in st.items():
        e = np.abs(hilbert(np.asarray(x, float), axis=1)).sum(0)
        med = np.median(e)
        imp = e.max() / med if med > 0 else 0.0
        dur = (e > e.max() / 2).sum() / A.FS
        if imp > 5 and dur < 3:
            good.append(s)
    return g, len({site.get(s, s) for s in good})


def _bp28_score(args):
    from scipy.signal import butter, sosfiltfilt
    sos = butter(4, [2.0 / (A.FS / 2), 8.0 / (A.FS / 2)], "band", output="sos")
    out = []
    for g, st, site in args:
        out.append(_famscore_one((g, {s: sosfiltfilt(sos, v, axis=-1) for s, v in st.items()}, site)))
    return out


def site_scores(ST, gids, stas, workers=16):
    """{grp: n impulsive sites} of each template's current stack, band-passed 2-8 Hz (linear, so
    = the mean of 2-8 Hz-filtered unit-max(1-8 Hz) segments; no data re-read)."""
    import ab_mute as MU
    sid = MU.sites(list(stas))
    site = {s: int(k) for s, k in zip(stas, sid)}
    gl = [g for g in gids if g in ST]
    chunks = [[(g, ST[g], site) for g in gl[i::workers]] for i in range(workers)]
    res = {}
    with ProcessPoolExecutor(workers, mp_context=mp.get_context("fork")) as ex:
        for part in ex.map(_bp28_score, chunks):
            res.update(dict(part))
    return {g: res.get(g, 0) for g in gids}


def _fair_pairs(det, gids, gi_all):
    """-> DataFrame lo, hi (template indices, lo < hi), cnt_r, lag_r (t_hi - t_lo, s), cnt_s."""
    pairs, _ = G._numba()
    nbin = int(round(FAIR_WIN / FAIR_BIN))
    K = 2 * nbin + 1
    Gn = gids.size
    n = np.bincount(gi_all, minlength=Gn).astype(np.float64)
    tt = det.t.to_numpy(np.float64)
    codes = {"r": [], "p": [], "m": []}          # real, shift +60, shift -60
    for _, idx in det.groupby("date").indices.items():
        t = tt[idx]; g = gi_all[idx]
        o = np.argsort(t, kind="stable"); t, g = t[o], g[o]
        ia, gb, d = pairs(t, g, FAIR_SHIFT + FAIR_WIN + FAIR_BIN)
        if ia.size == 0:
            continue
        ga = g[ia]
        keep = ga < gb                                  # one orientation; d = t_b - t_a
        ga, gb, d = ga[keep], gb[keep], d[keep]
        pc = ga.astype(np.int64) * Gn + gb
        for lab, dd in (("r", d), ("p", d - FAIR_SHIFT), ("m", d + FAIR_SHIFT)):
            m = np.abs(dd) <= FAIR_WIN
            k = np.round(dd[m] / FAIR_BIN).astype(np.int64) + nbin
            codes[lab].append(pc[m] * K + k)
    out = {}
    for lab in codes:
        c = np.concatenate(codes[lab]) if codes[lab] else np.zeros(0, np.int64)
        u, cnt = np.unique(c, return_counts=True)
        out[lab] = (u // K, u % K, cnt)
    # prune: a pair can only be a duplicate if its TOTAL real count within +-WIN > THR * min(n)
    pr, kr, cr = out["r"]
    up, inv = np.unique(pr, return_inverse=True)
    tot = np.bincount(inv, weights=cr)
    lo_, hi_ = up // Gn, up % Gn
    mn = np.minimum(n[lo_], n[hi_])
    cand = tot > FAIR_THR * mn
    up = up[cand]
    if not up.size:
        return pd.DataFrame(columns=["lo", "hi", "cnt_r", "lag_r", "cnt_s"])
    ker = np.ones(2 * FAIR_NB + 1)

    def modal(lab):
        p, k, c = out[lab]
        sel = np.isin(p, up)
        row = np.searchsorted(up, p[sel])
        H = np.zeros((up.size, K))
        np.add.at(H, (row, k[sel]), c[sel])
        Hs = np.apply_along_axis(lambda r: np.convolve(r, ker, "same"), 1, H)
        kk = Hs.argmax(1)
        return Hs[np.arange(up.size), kk], (kk - nbin) * FAIR_BIN

    cnt_r, lag_r = modal("r")
    cnt_p, _ = modal("p")
    cnt_m, _ = modal("m")
    return pd.DataFrame(dict(lo=up // Gn, hi=up % Gn, cnt_r=cnt_r, lag_r=lag_r,
                             cnt_s=np.maximum(cnt_p, cnt_m)))


def merge_fair(det, score, tpl_len=TPL_LEN):
    """-> (merged det [grp = representative, t aligned to it, src_grp, t_raw], family table,
    straddle table, stats).  See the block comment above."""
    _, nms = G._numba()
    gids = np.sort(det.grp.unique()).astype(np.int64)
    Gn = gids.size
    gi_all = np.searchsorted(gids, det.grp.to_numpy(np.int64))
    n = np.bincount(gi_all, minlength=Gn).astype(np.float64)
    s = np.array([score.get(int(g), -1) for g in gids], float)
    P = _fair_pairs(det, gids, gi_all)
    stats = dict(n_families_in=int(Gn), n_det_in=int(len(det)), n_pairs_tested=int(len(P)))
    if len(P):
        mn = np.minimum(n[P.lo], n[P.hi])
        P["excess"] = (P.cnt_r - P.cnt_s) / mn
        D = P[P.excess > FAIR_THR].copy()
    else:
        D = P.assign(excess=[])
    real = s >= FAIR_REAL
    strad = (real[D.lo.to_numpy(int)] != real[D.hi.to_numpy(int)]) if len(D) else np.zeros(0, bool)
    S_ = D[strad].copy()
    D = D[~strad]
    stats.update(n_dup_edges=int(len(D)), n_straddle_kept_separate=int(len(S_)))
    nb = {}
    for lo, hi, e, lag in zip(D.lo.to_numpy(int), D.hi.to_numpy(int), D.excess.to_numpy(), D.lag_r.to_numpy()):
        nb.setdefault(lo, []).append((hi, e, lag))       # lag = t_hi - t_lo
        nb.setdefault(hi, []).append((lo, e, -lag))      # lag = t_lo - t_hi
    order = np.lexsort((np.arange(Gn), -n, -s))          # score desc, then n_det desc
    is_rep = np.zeros(Gn, bool)
    rep_of = np.arange(Gn)
    off = np.zeros(Gn)                                   # t_rep - t_member at the modal lag
    for i in order:
        best, be, bl = -1, FAIR_THR, 0.0
        for j, e, lag in nb.get(i, ()):
            if is_rep[j] and e > be:
                best, be, bl = j, e, lag
        if best < 0:
            is_rep[i] = True
        else:
            rep_of[i] = best
            off[i] = bl                                  # t_rep = t_member + lag
    out = det.copy()
    out["src_grp"] = out.grp.astype(np.int64)
    out["t_raw"] = out.t
    out["grp"] = gids[rep_of[gi_all]]
    out["t"] = out.t.to_numpy(np.float64) + off[gi_all]
    out = out.sort_values(["grp", "t"], kind="stable").reset_index(drop=True)
    keep = nms(out.grp.to_numpy(np.int64), out.t.to_numpy(np.float64),
               np.argsort(-out.cc.to_numpy(np.float64), kind="stable").astype(np.int64), tpl_len)
    out = out[keep].reset_index(drop=True)
    fam = pd.DataFrame(dict(grp=gids, n_det=n.astype(int), score=s.astype(int), rep=gids[rep_of],
                            is_rep=is_rep, lag_to_rep_s=off))
    S_ = S_.assign(a=gids[S_.lo.to_numpy(int)], b=gids[S_.hi.to_numpy(int)],
                   score_a=s[S_.lo.to_numpy(int)], score_b=s[S_.hi.to_numpy(int)],
                   n_a=n[S_.lo.to_numpy(int)].astype(int), n_b=n[S_.hi.to_numpy(int)].astype(int))
    stats.update(n_components_merged=int(np.unique(rep_of[~is_rep]).size),
                 n_templates_merged_away=int((~is_rep).sum()), n_families_out=int(is_rep.sum()),
                 n_det_out=int(len(out)),
                 n_reps_ge3=int((real & is_rep).sum()), n_tpl_ge3=int(real.sum()))
    return out, fam, S_[["a", "b", "score_a", "score_b", "n_a", "n_b", "lag_r", "excess"]], stats


# ------------------------------------------------------------------ driver (ab_grow.main) -------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", required=True, choices=["NW", "SW", "NVI", "MENDO"])
    ap.add_argument("--tpl", default="", help="cluster stacks npz from ab_cluster (tpl_it0.npz)")
    ap.add_argument("--out", default="")
    ap.add_argument("--pool", default="", help="default = the 25-station pool A.POOL[sub]")
    ap.add_argument("--dates", default="", help="default = A.scan_days(sub)")
    ap.add_argument("--tpl-len", type=float, default=TPL_LEN)
    ap.add_argument("--net", default="mean", choices=["mean", "sum"])
    ap.add_argument("--nmad", type=float, default=A.NMAD)
    ap.add_argument("--min-sta-valid", type=int, default=MIN_STA_VALID)
    ap.add_argument("--cut", default="common", choices=["common", "envelope"])
    ap.add_argument("--max-sta", type=int, default=MAX_STA)
    ap.add_argument("--max-stack", type=int, default=0, help="0 = all detections, as the paper says")
    ap.add_argument("--max-iters", type=int, default=8)
    ap.add_argument("--tol", type=float, default=0.01, help="converged when |dN|/N < tol")
    ap.add_argument("--prep-workers", type=int, default=PREP_WORKERS, help="<= 8")
    ap.add_argument("--prep-threads", type=int, default=PREP_THREADS)
    ap.add_argument("--restack-workers", type=int, default=RESTACK_WORKERS, help="<= 8")
    ap.add_argument("--no-merge", action="store_true",
                    help="skip merge_duplicates (2026-09-28 Merlin step 1); the matched filter's "
                         "tpl-len NMS still declusters within a template. Default: merge as before.")
    ap.add_argument("--merge", default="most", choices=["most", "fair"],
                    help="most = merge_duplicates (representative = most detections; default, as before); "
                         "fair = merge_fair (2026-09-28: +-15 s modal lag with +60 s control, representative "
                         "= higher 2-8 Hz site-count score, star, never across the 3-site line)")
    ap.add_argument("--score-workers", type=int, default=16)
    a = ap.parse_args()
    out = a.out or f"{A.RES}/{a.sub}/grow"
    os.makedirs(out, exist_ok=True)
    pub = A.SUBARRAY[a.sub]
    stas = A.pool_from_arg(a.pool, a.sub)
    dates = a.dates.split(",") if a.dates else A.scan_days(a.sub)
    tpl_npz = a.tpl or f"{A.RES}/{a.sub}/clusters/tpl_it0.npz"
    ST = load_stacks(tpl_npz)
    A.log(f"[grow-fast {a.sub}] {len(ST)} cluster stacks, pool {stas}, {len(dates)} dates")
    dc = A.DayCache(stas)
    rows, prev = [], None
    for it in range(a.max_iters + 1):
        TPL, info = cut_templates(ST, tpl_len=a.tpl_len, max_sta=a.max_sta, cut=a.cut,
                                  min_sta_valid=a.min_sta_valid)
        info.to_csv(f"{out}/tpl_it{it}_stations.csv", index=False)
        if not TPL:
            A.log("[grow] no scannable template")
            break
        np.savez_compressed(f"{out}/tpl_it{it}.npz",
                            **{f"{g}_{s}": TPL[g][s][0].astype(np.float32)
                               for g in TPL for s in TPL[g]},
                            meta=np.array(json.dumps(
                                {f"{g}_{s}": [TPL[g][s][1], TPL[g][s][2]]
                                 for g in TPL for s in TPL[g]})))
        raw = matched_filter(dc, TPL, dates, nmad=a.nmad, net=a.net, min_sta_valid=a.min_sta_valid,
                             decl_s=a.tpl_len, tpl_len=a.tpl_len,
                             prep_workers=a.prep_workers, prep_threads=a.prep_threads)
        raw.to_parquet(f"{out}/detraw_it{it}.parquet", index=False)
        import time as _time
        _t0 = _time.time()
        if a.no_merge:
            det, mst = raw, dict(n_templates_merged_away=0, n_components_merged=0, n_dup_edges=0,
                                 no_merge=True)
        elif len(raw) and a.merge == "fair":
            sc = site_scores(ST, sorted(TPL), stas, workers=a.score_workers)
            _ts = _time.time() - _t0
            pd.DataFrame(dict(grp=list(sc), n_sites_3km=list(sc.values()))).to_csv(
                f"{out}/fair_scores_it{it}.csv", index=False)
            det, fam, strad, mst = merge_fair(raw, sc, tpl_len=a.tpl_len)
            mst["score_s"] = round(_ts, 1)
            fam.to_csv(f"{out}/families_it{it}.csv", index=False)
            strad.to_csv(f"{out}/fair_straddle_it{it}.csv", index=False)
        elif len(raw):
            det, fam, mst = merge_duplicates(raw)
            fam.to_csv(f"{out}/families_it{it}.csv", index=False)
        else:
            det, mst = raw, dict(n_templates_merged_away=0, n_components_merged=0, n_dup_edges=0)
        A.log(f"[merge it{it}] {mst} ({_time.time() - _t0:.0f} s)")
        det.to_parquet(f"{out}/det_it{it}.parquet", index=False)
        per = det.groupby("grp").size() if len(det) else pd.Series(dtype=int)
        rows.append(dict(iter=it, n_templates=len(TPL), n_det_raw=int(len(raw)),
                         n_dup_edges=int(mst.get("n_dup_edges", 0)),
                         n_merge_groups=int(mst.get("n_components_merged", 0)),
                         n_templates_merged_away=int(mst.get("n_templates_merged_away", 0)),
                         n_straddle=int(mst.get("n_straddle_kept_separate", 0)),
                         n_tpl_ge3=int(mst.get("n_tpl_ge3", -1)), n_reps_ge3=int(mst.get("n_reps_ge3", -1)),
                         n_families=int(per.size), n_detections=int(len(det)),
                         n_templates_with_det=int(per.size),
                         per_tpl_median=float(per.median()) if per.size else 0.0,
                         per_tpl_max=int(per.max()) if per.size else 0,
                         tpl_snr_median=float(info.snr_med.median()) if len(info) else np.nan))
        pd.DataFrame(rows).to_csv(f"{out}/iters.csv", index=False)
        A.log(f"[iter {it}] {len(TPL)} templates -> {len(raw)} raw -> {per.size} families, "
              f"{len(det)} detections; per family median {rows[-1]['per_tpl_median']:.0f} "
              f"max {rows[-1]['per_tpl_max']}; template SNR median {rows[-1]['tpl_snr_median']:.2f}")
        if prev is not None and prev > 0 and abs(len(det) - prev) / prev < a.tol:
            A.log(f"[grow] converged at iteration {it}")
            break
        prev = len(det)
        if it == a.max_iters or len(det) == 0:
            break
        _t0 = _time.time()
        ST = restack(dc, stas, det, max_stack=a.max_stack, workers=a.restack_workers)
        A.log(f"[restack it{it}] {len(ST)} families ({_time.time() - _t0:.0f} s)")
    it_final = rows[-1]["iter"] if rows else -1
    json.dump(dict(subarray=a.sub, pool=stas, dates=len(dates), iterations=rows,
                   n_final_templates=rows[-1]["n_templates_with_det"] if rows else 0,
                   n_final_templates_published=pub["n_tpl_pub"],
                   tpl_len_s=a.tpl_len, net=a.net, nmad=a.nmad, min_sta_valid=a.min_sta_valid, cut=a.cut,
                   max_sta=a.max_sta, max_stack=a.max_stack, final_iter=it_final,
                   no_merge=bool(a.no_merge), merge=("none" if a.no_merge else a.merge),
                   code="ab_grow_fast.py"),
              open(f"{out}/summary.json", "w"), indent=1)
    if rows:
        A.log(f"[grow-fast {a.sub}] final {rows[-1]['n_templates_with_det']} templates "
              f"(Royer: {pub['n_tpl_pub']})")


if __name__ == "__main__":
    main()
