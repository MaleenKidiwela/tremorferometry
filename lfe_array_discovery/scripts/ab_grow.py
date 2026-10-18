#!/usr/bin/env python
"""STAGE 5 -- Royer & Bostock (2014) sec.3.3 NETWORK CROSS-CORRELATION, literally.

  "Iterative network cross-correlation and stacking was used to register further detections and
   improve template SNR.  We performed network cross-correlation by choosing subsets of 7 to 10 high
   SNR, three-component stations available within a template and scanning through all available
   tremor dates.  As before, when the summed network cross-correlation coefficient exceeded 8 times
   the median absolute deviation, now for a 24-hour period, we logged detections.  New templates were
   formed by stacking waveforms (normalized to unit maximum amplitude across three components) for
   all newly registered detections.  After several iterations of network cross-correlation and
   stacking, we obtained final sets of templates suitable for location and waveform analysis.  Our
   final suites of templates number 122, 54 and 47 ... Each template possesses 100's to 1000's of
   independent detections."

Published yields to check against: 122 final templates (NW), 54 (SW).

Literal elements that differ from `k_grow.py`:
  * MAD over the 24-HOUR period, not per group over a 12 h block;
  * the restack normalisation is "unit maximum amplitude across three components" -- max|x| over the
    three components jointly, NOT the rms normalisation of `k_common.norm3`;
  * "all newly registered detections" enter the stack, so there is no MAX_STACK = 200 cap
    (--max-stack 0 is the default; a cap is available only as a reported ablation);
  * iteration runs to convergence (--tol on the change in detection count), not a fixed 3.

MASK DISCIPLINE: a detection window that touches an invalid sample is not scored; the network
function divides by the per-sample count of contributing stations rather than a constant; the MAD
pool contains only samples where at least --min-sta stations are valid.

GAP GUARD (2026-09-25, matched_filter only; the _ref keeps the old `e > 0`): a channel-window
enters the per-channel cc only if it passes ab_common's relative guard (exact energies; e > 1e-3 x
the channel's truly-non-zero daily median; <= 1 % zero/constant samples; dead-channel and
network floors).  `e > 0` let residue windows through -> network cc up to 106 on 7 NW days.

  python scripts/ab_grow.py --sub NW --tpl results/ab_royer/NW/clusters/tpl_it0.npz \
      --out results/ab_royer/NW/grow
-> <out>/tpl_it<K>.npz      templates after iteration K
-> <out>/detraw_it<K>.parquet raw detector output: grp, t, cc, nmad, n_sta, date
-> <out>/det_it<K>.parquet  after the duplicate-family merge (+ src_grp, t_raw)
-> <out>/families_it<K>.csv template -> representative, lag
-> <out>/iters.csv          per iteration: templates, raw/merged detections, merges, families,
                            median/max per family, median template SNR (Bostock 2012: iterate
                            until template SNR and detection count stop improving; --tol/--max-iters)
-> <out>/summary.json       final counts beside Royer's published yields
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

LONG_PRE, LONG_LEN = 12.0, 44.0      # the cluster-stack window written by ab_cluster
TPL_LEN = 15.0                       # SILENT in sec.3.3; 15 s is sec.3.1's own window (swept)
TPL_PRE = 2.0                        # --cut envelope only: cut this far before the envelope peak
MIN_STA = 7                          # sec.3.3 "subsets of 7 to 10" -- the SIZE of the chosen subset
MAX_STA = 10
MIN_STA_VALID = 4                    # availability floor per sample, SAME as ab_auto --min-sta
                                     # (paper SILENT; AUDIT FIX A6 -- was 7 = every station present)


def envelope(x, fs=A.FS, smooth_s=0.4):
    from scipy.signal import hilbert
    e = np.sqrt((np.abs(hilbert(x, axis=-1)) ** 2).sum(0))
    k = max(int(round(smooth_s * fs)), 1)
    return np.convolve(e, np.ones(k) / k, mode="same")


def unit_max3(seg):
    """Royer sec.3.3's stacking normalisation: unit MAXIMUM amplitude across three components."""
    m = float(np.abs(seg).max())
    return seg / m if m > 0 else seg


def cut_templates(ST, tpl_len=TPL_LEN, tpl_pre=TPL_PRE, max_sta=MAX_STA, cut="common",
                  min_sta_valid=MIN_STA_VALID):
    """{cid: {sta: (3,n)}} long stacks -> {cid: {sta: (tpl, off_s, snr)}}.

    sec.3.3: "choosing subsets of 7 to 10 high SNR, three-component stations available within a
    template".  AUDIT FIX A4: this is a RANKING -- stations are ordered by the stack's SNR (peak over
    median envelope) and the top `max_sta` are used; there is no SNR cut (the old MIN_PEAK_MAD=3
    gate is gone).  The SNR is kept as a column.  A template is scannable when it has at least
    `min_sta_valid` stations (the same availability floor as stage 3); how many it used is a column.

    AUDIT FIX A5: `cut="common"` (default) takes the sec.3.1/3.2 15 s detection window itself,
    [LONG_PRE, LONG_PRE + tpl_len] of the long stack, identically for every station, so the network
    moveout is preserved and P and S both stay inside.  `cut="envelope"` (the old per-station cut
    tpl_pre before the envelope peak) is an ablation."""
    nt = int(round(tpl_len * A.FS))
    nL = int(round(LONG_LEN * A.FS))
    out, info = {}, []
    for cid, T in ST.items():
        cand = []
        for s, st in T.items():
            if not np.isfinite(st).all():
                continue
            # "three-component stations": a vertical-only station-day stack (zero horizontals)
            # is not eligible for the growth subset
            if not ((st[1] ** 2).sum() > 0 and (st[2] ** 2).sum() > 0):
                continue
            env = envelope(st)
            if cut == "common":
                i0 = int(round(LONG_PRE * A.FS))
                k = i0 + int(np.argmax(env[i0:i0 + nt]))
            else:
                k = int(np.argmax(env))
                i0 = k - int(round(tpl_pre * A.FS))
            if i0 < 0 or i0 + nt > min(nL, st.shape[1]):
                continue
            med = float(np.median(env))
            h = float(env[k] / med) if med > 0 else 0.0
            seg = st[:, i0:i0 + nt].copy()
            if not (seg ** 2).sum() > 0:
                continue
            cand.append((h, s, seg, i0 / A.FS - LONG_PRE))
        cand.sort(key=lambda r: -r[0])
        sel = cand[:max_sta]
        info.append(dict(cid=int(cid), n_sta_stack=len(T), n_sta_usable=len(cand),
                         n_sta_used=len(sel), scannable=bool(len(sel) >= min_sta_valid),
                         snr_med=float(np.median([r[0] for r in sel])) if sel else np.nan,
                         stations=",".join(r[1] for r in sel)))
        if len(sel) >= min_sta_valid:
            out[int(cid)] = {r[1]: (r[2], r[3], r[0]) for r in sel}
    return out, pd.DataFrame(info)


def _sliding_max(v, dec):
    """Exactly torch max_pool1d(v[None], 2*dec+1, stride=1, padding=dec)[0] (implicit -inf padding)
    by log-doubling (van Herk style): O(log k) passes instead of O(k) per sample."""
    import torch
    n = v.shape[1]
    w = 2 * dec + 1
    m = torch.nn.functional.pad(v, (dec, dec), value=float("-inf"))
    p = 1
    while 2 * p <= w:
        m = torch.maximum(m[:, :-p], m[:, p:])        # m[i] = max x[i .. i+2p-1]
        p *= 2
    return torch.maximum(m[:, :n], m[:, w - p:w - p + n])


def matched_filter(dc, TPL, dates, nmad=A.NMAD, net="mean", min_sta_valid=MIN_STA_VALID,
                   decl_s=TPL_LEN, tpl_len=TPL_LEN, dev=None, L=8192, gchunk=256, fbatch=64,
                   _dump=None):
    """Network cross-correlation over `dates`, 8 x MAD over each 24-hour period.

    SAME method and outputs as `matched_filter_ref` (the pre-2026-09-25 code); only the arithmetic
    is reorganised for speed:
      * overlap-save block correlation at block length L (valid outputs per block S = L - nt + 1)
        instead of one day-length FFT per template channel;
      * every template channel's L-point FFT is computed ONCE per call (resident on the GPU) and
        reused for every date; each station-day is FFT'd once per date;
      * complex64 products / float32 inverse FFTs (the ref cast cc to float32 after a float64 FFT);
      * the per-template offset shift-and-add is vectorised per unique offset.
    Unchanged semantics: per-channel normalised cc with the validity-masked window energies
    (DEN = inf where invalid, no clamping on zero-filled gaps), per-station offset shift-and-add in
    the same station order, --net mean/sum, per-sample valid-channel and valid-station counts,
    min_sta_valid, 8 x MAD over the 24 h, declustering at decl_s.

    AUDIT FIX A7: the network function is the per-CHANNEL normalised cc summed over every
    (station, component) channel of the template (Gibbons & Ringdal 2006; Shelly et al. 2006),
    `--net sum`, or that sum divided by the per-sample count of valid channels, `--net mean`.
    A sample enters the day's median/MAD pool when at least `min_sta_valid` template stations
    are valid there (AUDIT FIX A6)."""
    import torch
    dev = dev or A.torch_dev()
    torch.backends.cuda.matmul.allow_tf32 = False
    nt = int(round(tpl_len * A.FS))
    S = L - nt + 1
    assert S > 0
    gs_all = sorted(TPL)
    G_all = len(gs_all)
    stas = sorted({s for g in TPL for s in TPL[g]})
    # ---- template FFTs, once per call: per station, rows sorted by global template index -----
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
        # the per-channel template norm is folded into the template spectrum (1/nrm; 0 where the
        # component is dead, which the ref masks out via tok = nrm > 0)
        inv = np.where(nrm > 0, 1.0 / np.where(nrm > 0, nrm, 1.0), 0.0)
        Tf = (torch.fft.rfft(torch.from_numpy(M).to(dev), n=L).conj()
              * torch.from_numpy(inv).to(dev)[:, :, None]).to(torch.complex64)
        TS[s] = dict(gi=np.asarray(gi, np.int64), Tf=Tf,
                     live=torch.tensor(nrm > 0, device=dev, dtype=torch.float32), ks=ks)
        del M
    out = []
    for date in dates:
        D = dc.day(date)
        have = [s for s in stas if s in D]
        if len(have) < min_sta_valid:
            A.log(f"[mf] {date} only {len(have)} stations -> skipped")
            continue
        N = min(D[s][1].shape[1] for s in have)
        NV = N - nt + 1
        B = -(-NV // S)
        NP = B * S + nt - 1
        Xf, IDEN, OKf, T0 = {}, {}, {}, {}
        from concurrent.futures import ThreadPoolExecutor
        # 2026-09-25 GAP GUARD (replaces `e > 0`): exact window energies + relative guard against
        # the channel's own truly-non-zero median, flat-sample fraction, and the network floor
        # (A.window_stats_guarded / A.day_refs).  Refs use the whole station-day.
        REF = A.day_refs(D, nt)                               # whole pool, network floor
        with ThreadPoolExecutor(min(len(have), 16)) as ex:
            WS = dict(zip(have, ex.map(lambda s: A.window_stats_guarded(
                D[s][1][:, :N], D[s][2][:N], nt, 1, ref=REF[s]), have)))
        for s in have:
            t0, X, m = D[s]
            X = X[:, :N]
            ew, okw = WS.pop(s)                                # (3, NV), (3, NV) guarded
            Xp = torch.zeros((3, NP), dtype=torch.float64, device=dev)
            Xp[:, :N] = torch.from_numpy(np.ascontiguousarray(X, np.float64)).to(dev)
            Xf[s] = torch.fft.rfft(Xp.unfold(-1, L, S), n=L).to(torch.complex64)   # (3,B,F)
            e = torch.from_numpy(ew).to(dev)
            okc = torch.from_numpy(okw).to(dev)
            den = torch.where(okc, e.clamp_min(0).sqrt(),
                              torch.full_like(e, float("inf"))).float()       # = ref DEN
            iden = torch.zeros((3, B * S), device=dev)
            iden[:, :NV] = torch.where(okc, 1.0 / den, torch.zeros_like(den))  # 0 where invalid
            IDEN[s] = iden.view(3, B, S)
            OKf[s] = okc.float()
            T0[s] = float(t0)
            del Xp, e, den, okc
        ndet = 0
        for c0 in range(0, G_all, gchunk):
            c1 = min(c0 + gchunk, G_all)
            gs = gs_all[c0:c1]
            G = len(gs)
            acc = torch.zeros((G, NV), device=dev)
            cch = torch.zeros((G, NV), device=dev)
            csta = torch.zeros((G, NV), device=dev)
            for s in have:
                T = TS[s]
                a0, a1 = np.searchsorted(T["gi"], [c0, c1])
                for b0 in range(a0, a1, fbatch):
                    b1 = min(b0 + fbatch, a1)
                    loc = T["gi"][b0:b1] - c0                               # rows in this chunk
                    P = Xf[s][None] * T["Tf"][b0:b1, :, None, :]            # (r,3,B,F)
                    cc = torch.fft.irfft(P, n=L)[..., :S]                   # (r,3,B,S) cc/nrm
                    del P
                    w = IDEN[s]
                    ccs = cc[:, 0] * w[0]                                    # sum_c cc/(nrm DEN)
                    ccs.addcmul_(cc[:, 1], w[1])
                    ccs.addcmul_(cc[:, 2], w[2])
                    del cc
                    ccs = ccs.reshape(b1 - b0, B * S)[:, :NV]
                    nch = T["live"][b0:b1] @ OKf[s]                          # (r,NV) valid channels
                    ks = T["ks"][b0:b1]
                    for k in np.unique(ks):
                        k = int(k)
                        sel = np.flatnonzero(ks == k)
                        lo, hi = max(-k, 0), min(NV - k, NV)
                        if hi <= lo:
                            continue
                        rs = torch.from_numpy(loc[sel]).to(dev)
                        ss = torch.from_numpy(sel).to(dev)
                        src_n = nch[ss, lo + k:hi + k]
                        acc[:, lo:hi].index_add_(0, rs, ccs[ss, lo + k:hi + k])
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
            pooled = _sliding_max(vals, dec)          # == max_pool1d(k=2dec+1, pad=dec), exact
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
        A.log(f"[mf] {date} {len(have)} sta -> {ndet}")
        del Xf, IDEN, OKf
        if dev == "cuda":
            torch.cuda.empty_cache()
    del TS
    return (pd.concat(out, ignore_index=True) if out
            else pd.DataFrame(columns=["grp", "t", "cc", "cc_mean", "nmad", "n_sta", "n_ch",
                                       "thr", "date"]))


def matched_filter_ref(dc, TPL, dates, nmad=A.NMAD, net="mean", min_sta_valid=MIN_STA_VALID,
                       decl_s=TPL_LEN, tpl_len=TPL_LEN, dev=None, batch=24, chunk=64, _dump=None):
    """REFERENCE (pre-2026-09-25) implementation, kept for the equivalence test.  Network cross-correlation over `dates`, 8 x MAD over each 24-hour period.

    AUDIT FIX A7: the network function is the per-CHANNEL normalised cc summed over every
    (station, component) channel of the template (Gibbons & Ringdal 2006; Shelly et al. 2006),
    `--net sum`, or that sum divided by the per-sample count of valid channels, `--net mean`.
    A sample enters the day's median/MAD pool when at least `min_sta_valid` template stations
    are valid there (AUDIT FIX A6)."""
    import torch
    dev = dev or A.torch_dev()
    torch.backends.cuda.matmul.allow_tf32 = False
    nt = int(round(tpl_len * A.FS))
    gs_all = sorted(TPL)
    stas = sorted({s for g in TPL for s in TPL[g]})
    out = []
    for date in dates:
        D = dc.day(date)
        have = [s for s in stas if s in D]
        if len(have) < min_sta_valid:
            A.log(f"[mf] {date} only {len(have)} stations -> skipped")
            continue
        N = min(D[s][1].shape[1] for s in have)
        NV = N - nt + 1
        Xf, DEN, OK, T0 = {}, {}, {}, {}
        for s in have:
            t0, X, m = D[s]
            X = X[:, :N]
            m = m[:N]
            ew, vw = A.window_stats(X, m, nt, 1)              # (3, NV), (NV,)
            Xt = torch.from_numpy(np.ascontiguousarray(X, np.float64)).to(dev)
            Xf[s] = torch.fft.rfft(Xt, n=N)
            e = torch.from_numpy(ew).to(dev)                  # (3, NV) float64
            okc = torch.from_numpy(vw).to(dev)[None] & (e > 0)  # (3, NV)
            DEN[s] = torch.where(okc, e.clamp_min(0).sqrt(),
                                 torch.full_like(e, float("inf"))).float()
            OK[s] = okc
            T0[s] = float(t0)
            del Xt, e
        ndet = 0
        for c0 in range(0, len(gs_all), chunk):     # template chunks bound GPU memory
            gs = gs_all[c0:c0 + chunk]
            acc = torch.zeros((len(gs), NV), device=dev)
            cch = torch.zeros((len(gs), NV), device=dev)        # valid channels
            csta = torch.zeros((len(gs), NV), device=dev)       # valid stations
            for b0 in range(0, len(gs), batch):
                sub = gs[b0:b0 + batch]
                for s in have:
                    idx = [i for i, g in enumerate(sub) if s in TPL[g]]
                    if not idx:
                        continue
                    M = np.zeros((len(idx), 3, N), np.float64)
                    offs, nrm = [], []
                    for r, i in enumerate(idx):
                        tpl, off, _ = TPL[sub[i]][s]
                        M[r, :, :nt] = tpl
                        offs.append(off)
                        nrm.append(np.sqrt((tpl ** 2).sum(1)))        # per component
                    Tf = torch.fft.rfft(torch.from_numpy(M).to(dev), n=N)
                    cc = torch.fft.irfft(Xf[s][None] * Tf.conj(), n=N)[..., :NV].float()  # (r,3,NV)
                    del Tf
                    nrmt = torch.tensor(np.array(nrm), device=dev, dtype=torch.float32)  # (r,3)
                    tok = (nrmt > 0)[:, :, None] & OK[s][None]                          # (r,3,NV)
                    cc = cc / (nrmt[:, :, None] * DEN[s][None]).clamp_min(1e-30)
                    cc = torch.where(tok, cc, torch.zeros_like(cc))
                    ccs = cc.sum(1)                                                      # (r,NV)
                    nch = tok.float().sum(1)
                    del cc, tok
                    for r, i in enumerate(idx):
                        k = int(round(offs[r] * A.FS))
                        lo, hi = max(-k, 0), min(NV - k, NV)
                        if hi <= lo:
                            continue
                        acc[b0 + i, lo:hi] += ccs[r, lo + k:hi + k]
                        cch[b0 + i, lo:hi] += nch[r, lo + k:hi + k]
                        csta[b0 + i, lo:hi] += (nch[r, lo + k:hi + k] > 0).float()
                    del ccs, nch
                if dev == "cuda":
                    torch.cuda.empty_cache()
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
            pooled = torch.nn.functional.max_pool1d(vals[None], 2 * dec + 1, stride=1, padding=dec)[0]
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
        A.log(f"[mf] {date} {len(have)} sta -> {ndet}")
        del Xf, DEN, OK
        if dev == "cuda":
            torch.cuda.empty_cache()
    return (pd.concat(out, ignore_index=True) if out
            else pd.DataFrame(columns=["grp", "t", "cc", "cc_mean", "nmad", "n_sta", "n_ch",
                                       "thr", "date"]))


def restack(dc, stas, det, max_stack=0, long_pre=LONG_PRE, long_len=LONG_LEN):
    """"stacking waveforms (normalized to unit maximum amplitude across three components) for all
    newly registered detections".

    SMOKE-TEST FIX: iterate DATE-outer (the DayCache holds one date), so each day is loaded once per
    iteration instead of once per (template, date) -- same sums, ~#templates fewer day loads."""
    nL = int(round(long_len * A.FS))
    if max_stack:
        det = det.groupby("grp", group_keys=False).apply(lambda g: g.nlargest(max_stack, "cc"))
    acc, cnt = {}, {}
    nref = int(round(TPL_LEN * A.FS))
    rc = A.RefCache(dc, nref)                  # 2026-09-25 gap guard on the stacked channels too
    for date, dd in det.sort_values("t").groupby("date"):
        refs = {s: rc.get(s, date) for s in stas}
        for r in dd.itertuples():
            g = int(r.grp)
            for s in stas:
                seg = dc.cut(s, r.t - long_pre, long_len, date)
                if seg is None:
                    continue
                okc = refs[s] > 0                  # dead channel that day (guard ref 0) -> 0
                if not okc.any():
                    continue
                if not okc.all():
                    seg = np.where(okc[:, None], seg, 0.0)
                key = (g, s)
                if key not in acc:
                    acc[key] = np.zeros((3, nL))
                    cnt[key] = 0
                acc[key] += unit_max3(seg[:, :nL])
                cnt[key] += 1
    ST = {}
    for (g, s), v in acc.items():
        if cnt[(g, s)] >= 2:
            ST.setdefault(g, {})[s] = v / cnt[(g, s)]
    return ST


# ------------------------------------------------------------------ duplicate-family merge -------
# 2026-09-25 (v4b runaway: 2,695 templates -> a median of 3 templates per 4-s bin, NW 3.56 M ->
# 19.3 M detections).  The Royer/Bostock group's rule (Thomas & Bostock 2015, notes/royer_sources/
# SOURCES.md S5): "no family shares more than 10 % of detection times, defined as a difference
# between detection times of <5 s, with any other family."  Applied after every iteration's
# detection and before the restack: families whose shared fraction (detections of the SMALLER
# family with a detection of the other within 5 s / its detection count) exceeds 10 % are linked;
# connected components merge into one family represented by the member with the most detections;
# the union of detections (each member's times shifted onto the representative's alignment by
# the modal pair lag), deduplicated to one per 5 s (strongest cc kept), is the family's stack set.
DUP_WIN_S = 5.0
DUP_FRAC = 0.10


def _numba():
    import numba

    @numba.njit(cache=False)
    def pairs(t, g, win):
        n = t.size
        c = 0
        for i in range(n):
            j = i + 1
            while j < n and t[j] - t[i] < win:
                if g[j] != g[i]:
                    c += 2
                j += 1
        A_ = np.empty(c, np.int64)
        B_ = np.empty(c, np.int64)
        D_ = np.empty(c, np.float64)
        k = 0
        for i in range(n):
            j = i + 1
            while j < n and t[j] - t[i] < win:
                if g[j] != g[i]:
                    A_[k] = i; B_[k] = g[j]; D_[k] = t[j] - t[i]; k += 1
                    A_[k] = j; B_[k] = g[i]; D_[k] = t[i] - t[j]; k += 1
                j += 1
        return A_, B_, D_

    @numba.njit(cache=False)
    def nms(fam, t, order, win):
        # fam, t sorted by (fam, t); order = indices by descending score
        n = t.size
        sup = np.zeros(n, np.bool_)
        keep = np.zeros(n, np.bool_)
        for q in range(n):
            i = order[q]
            if sup[i]:
                continue
            keep[i] = True
            j = i - 1
            while j >= 0 and fam[j] == fam[i] and t[i] - t[j] < win:
                sup[j] = True; j -= 1
            j = i + 1
            while j < n and fam[j] == fam[i] and t[j] - t[i] < win:
                sup[j] = True; j += 1
        return keep

    return pairs, nms


def merge_duplicates(det, win=DUP_WIN_S, frac=DUP_FRAC, lag_bin=0.1):
    """-> (merged det [grp = representative, t aligned to it, src_grp, t_raw], family table,
    stats dict).  See the block comment above."""
    from scipy import sparse
    pairs, nms = _numba()
    gids = np.sort(det.grp.unique()).astype(np.int64)
    G = gids.size
    gi_all = np.searchsorted(gids, det.grp.to_numpy(np.int64))
    n_fam = np.bincount(gi_all, minlength=G).astype(np.float64)
    M = sparse.csr_matrix((G, G))
    H = []
    for date, idx in det.groupby("date").indices.items():
        t = det.t.to_numpy(np.float64)[idx]
        g = gi_all[idx]
        o = np.argsort(t, kind="stable")
        t, g = t[o], g[o]
        a, b, d = pairs(t, g, win)
        if a.size == 0:
            continue
        # nearest partner per (detection, other family)
        key = a * G + b
        oo = np.lexsort((np.abs(d), key))
        first = np.ones(oo.size, bool)
        first[1:] = key[oo][1:] != key[oo][:-1]
        sel = oo[first]
        ga, gb, dd = g[a[sel]], b[sel], d[sel]
        M = M + sparse.csr_matrix((np.ones(ga.size), (ga, gb)), shape=(G, G))
        H.append(pd.DataFrame(dict(a=ga, b=gb, k=np.round(dd / lag_bin).astype(np.int32)))
                 .value_counts().rename("n").reset_index())
    M = M.tocoo()
    stats = dict(n_families_in=int(G), n_det_in=int(len(det)))
    # m[a, b] = detections of a with a b-detection within win; pair fraction = the SMALLER family's
    # m / n (ties: both directions), symmetrised
    Fr_small = sparse.csr_matrix((G, G))
    if M.nnz:
        sm = n_fam[M.row] <= n_fam[M.col]
        Fr_small = sparse.coo_matrix((M.data[sm] / n_fam[M.row[sm]], (M.row[sm], M.col[sm])),
                                     shape=(G, G)).tocsr()
        Fr_small = Fr_small.maximum(Fr_small.T).tocsr()
    stats["n_pairs_sharing"] = int(Fr_small.nnz // 2)
    stats["n_dup_edges"] = int((Fr_small.data > frac).sum() // 2)
    # GREEDY STAR merge (not connected components: on NW it0 single-linkage chaining put 1,763 of
    # 2,695 templates into ONE component with lags up to 26 s).  Families in descending detection
    # count; a family joins the representative it shares the largest fraction with, if that
    # fraction > frac; otherwise it becomes a representative.  So every merged member is a DIRECT
    # duplicate of its representative, and no two representatives share > frac.
    order = np.lexsort((np.arange(G), -n_fam))
    is_rep = np.zeros(G, bool)
    rep_of = np.arange(G)
    for i in order:
        row = Fr_small.getrow(i)
        best, bf = -1, frac
        for j, f in zip(row.indices, row.data):
            if is_rep[j] and f > bf:
                best, bf = j, f
        if best < 0:
            is_rep[i] = True
        else:
            rep_of[i] = best
    lab = rep_of                                  # component label = representative index
    rep = np.arange(G)                            # rep[lab] = lab
    ncomp = int(is_rep.sum())
    off = np.zeros(G)
    members = np.flatnonzero(~is_rep)
    if members.size and H:
        Hd = pd.concat(H, ignore_index=True)
        # modal lag t_rep - t_member... keep rows (a = member, b = its rep): k = (t_b - t_a)
        Hd = Hd[rep_of[Hd.a.to_numpy()] == Hd.b.to_numpy()]
        Hd = Hd[~is_rep[Hd.a.to_numpy()]]
        Hd = Hd.groupby(["a", "k"], as_index=False)["n"].sum()
        mode = Hd.sort_values("n").groupby("a").tail(1)
        # t_rep ~ t_member + k*lag_bin  ->  member time aligned to rep = t + k*lag_bin
        off[mode.a.to_numpy()] = -mode.k.to_numpy() * lag_bin
    multi = np.unique(rep_of[members]) if members.size else np.array([], int)
    out = det.copy()
    out["src_grp"] = out.grp.astype(np.int64)
    out["t_raw"] = out.t
    out["grp"] = gids[rep[lab[gi_all]]]
    out["t"] = out.t.to_numpy(np.float64) - off[gi_all]
    # one detection per win per family (strongest cc)
    out = out.sort_values(["grp", "t"], kind="stable").reset_index(drop=True)
    score = out.cc.to_numpy(np.float64)
    keep = nms(out.grp.to_numpy(np.int64), out.t.to_numpy(np.float64),
               np.argsort(-score, kind="stable").astype(np.int64), win)
    out = out[keep].reset_index(drop=True)
    fam = pd.DataFrame(dict(grp=gids, n_det=n_fam.astype(int), rep=gids[rep[lab]],
                            comp_size=np.bincount(lab)[lab], lag_to_rep_s=off))
    stats.update(n_components_merged=int(multi.size),
                 n_templates_merged_away=int(G - ncomp), n_families_out=int(ncomp),
                 n_det_out=int(len(out)))
    return out, fam, stats


def load_stacks(npz):
    z = np.load(npz, allow_pickle=True)
    ST = {}
    for k in z.files:
        if k == "meta":
            continue
        cid, sta = k.split("_", 1)
        ST.setdefault(int(cid), {})[sta] = z[k].astype(np.float64)
    return ST


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", required=True, choices=["NW", "SW", "NVI"])
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
    a = ap.parse_args()
    out = a.out or f"{A.RES}/{a.sub}/grow"
    os.makedirs(out, exist_ok=True)
    pub = A.SUBARRAY[a.sub]
    # sec.3.3: the 7-10 stations are chosen from the stations "available within a template", i.e.
    # the subarray's 25-station array of sec.2.1 / Fig.1D (A.POOL), scanning "all available tremor
    # dates" (A.scan_days, ab_tremor_days.py)
    stas = A.pool_from_arg(a.pool, a.sub)      # '' literal pool; 'ext:<R>' EXTENSION pool
    dates = a.dates.split(",") if a.dates else A.scan_days(a.sub)
    tpl_npz = a.tpl or f"{A.RES}/{a.sub}/clusters/tpl_it0.npz"
    ST = load_stacks(tpl_npz)
    A.log(f"[grow {a.sub}] {len(ST)} cluster stacks, pool {stas}, {len(dates)} dates")
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
                             decl_s=a.tpl_len, tpl_len=a.tpl_len)
        # raw detector output kept; det_it<K> is the MERGED catalogue (duplicate families merged,
        # grp = representative template of tpl_it<K>, t on its alignment) -> used by restack/locate
        raw.to_parquet(f"{out}/detraw_it{it}.parquet", index=False)
        import time as _time
        _t0 = _time.time()
        if len(raw):
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
        ST = restack(dc, stas, det, max_stack=a.max_stack)
    it_final = rows[-1]["iter"] if rows else -1
    json.dump(dict(subarray=a.sub, pool=stas, dates=len(dates), iterations=rows,
                   n_final_templates=rows[-1]["n_templates_with_det"] if rows else 0,
                   n_final_templates_published=pub["n_tpl_pub"],
                   tpl_len_s=a.tpl_len, net=a.net, nmad=a.nmad, min_sta_valid=a.min_sta_valid, cut=a.cut,
                   max_sta=a.max_sta, max_stack=a.max_stack, final_iter=it_final),
              open(f"{out}/summary.json", "w"), indent=1)
    if rows:
        A.log(f"[grow {a.sub}] final {rows[-1]['n_templates_with_det']} templates "
              f"(Royer: {pub['n_tpl_pub']})")


if __name__ == "__main__":
    main()
