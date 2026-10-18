#!/usr/bin/env python
# FROZEN COPY for avenue A4 (2026-09-25) of the file named without the a4f_ prefix (md5 at copy time in results/a4_seed/A4_LOG.md); do not edit.
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
pool contains only samples where at least --min-sta stations are valid.  No relative-energy guard.

  python scripts/ab_grow.py --sub NW --tpl results/ab_royer/NW/clusters/tpl_it0.npz \
      --out results/ab_royer/NW/grow
-> <out>/tpl_it<K>.npz      templates after iteration K
-> <out>/det_it<K>.parquet  grp, t, cc, nmad, n_sta, date
-> <out>/iters.csv          per iteration: n_templates, n_detections, median/max per template
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
import a4f_common as A                                                       # noqa: E402

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


def matched_filter(dc, TPL, dates, nmad=A.NMAD, net="mean", min_sta_valid=MIN_STA_VALID,
                   decl_s=TPL_LEN, tpl_len=TPL_LEN, dev=None, batch=24, chunk=64):
    """Network cross-correlation over `dates`, 8 x MAD over each 24-hour period.

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
    for date, dd in det.sort_values("t").groupby("date"):
        for r in dd.itertuples():
            g = int(r.grp)
            for s in stas:
                seg = dc.cut(s, r.t - long_pre, long_len, date)
                if seg is None:
                    continue
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
        det = matched_filter(dc, TPL, dates, nmad=a.nmad, net=a.net, min_sta_valid=a.min_sta_valid,
                             decl_s=a.tpl_len, tpl_len=a.tpl_len)
        det.to_parquet(f"{out}/det_it{it}.parquet", index=False)
        per = det.groupby("grp").size() if len(det) else pd.Series(dtype=int)
        rows.append(dict(iter=it, n_templates=len(TPL), n_detections=int(len(det)),
                         n_templates_with_det=int(per.size),
                         per_tpl_median=float(per.median()) if per.size else 0.0,
                         per_tpl_max=int(per.max()) if per.size else 0))
        pd.DataFrame(rows).to_csv(f"{out}/iters.csv", index=False)
        A.log(f"[iter {it}] {len(TPL)} templates -> {len(det)} detections; "
              f"per template median {rows[-1]['per_tpl_median']:.0f} max {rows[-1]['per_tpl_max']}")
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
