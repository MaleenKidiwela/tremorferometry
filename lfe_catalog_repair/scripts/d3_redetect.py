#!/usr/bin/env python
"""D3 stage 2 — GPU re-detection with 3-C cc stacking and per-template MAD thresholds.

Reuses the `scripts/densify_gnw_gpu.py` architecture verbatim where possible (CPU pool does
read -> resample_poly(40 Hz) -> demean -> 2-8 Hz bandpass -> despike-MAD; cupy lives only in the main
process; one rfft(day) per component, templates batched, GPU sliding norm, sparse crossings only over
PCIe; `_dedup_min_gap` for declustering) and ADDS what D1 showed was missing:

  1. all three components (Z + the two horizontals), cc_stack = mean of the per-component cc
  2. per-template, per-day robust noise statistics of the STACKED cc series:
         med_day  = median(cc_stack[::m])        mad_day = median(|cc_stack[::m] - med_day|)
     (stride-m subsample = 43,200 ~independent values/day; MAD is the raw median-absolute-deviation,
      NOT scaled by 1.4826 — the Shelly/EQcorrscan convention for `median + k*MAD` thresholds)
  3. an extraction floor of  min(med_day + K_EXTRACT*mad_day, OLD_THR)  so that BOTH the new
     MAD rules (k = 8/9/10, recomputed per 30-d window in post-processing from the stored daily
     stats) AND the old fixed rule (cc_stack >= 0.80) are exact subsets of what is stored
  4. declustering at min-separation = template length (2.0 s) — monotone-safe, so applying a higher
     threshold afterwards gives exactly the same peaks as declustering at that threshold
  5. per detection: cc_stack, cc per component, n_comp, SNR (peak |Z| in the 2 s template window /
     RMS of the preceding 10 s of Z, same band), the floor used and the daily median/MAD.

Outputs
  <out_prefix>_detections.parquet   one row per surviving peak
  <out_prefix>_daily.parquet        (template, day) x (med, mad, n_det, capped, n_valid)
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import sys
import time
import warnings
from datetime import datetime, timedelta
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, "src")
sys.path.insert(0, "scripts")
from scipy.signal import resample_poly   # noqa: E402

_G = {}


def _init_worker(root, station, fs, bandpass, despike_mad, comps):
    _G.update(root=Path(root), station=station, fs=fs, bandpass=bandpass,
              despike_mad=despike_mad, comps=comps)


def _preprocess(day):
    """read the day file, preprocess EVERY component onto ONE common [day, day+86400] grid.

    Returns (day, {comp: float32 array of n_fix_real}, day_start_epoch).
    """
    from obspy import read, UTCDateTime
    root, station, fs, bp = _G["root"], _G["station"], _G["fs"], _G["bandpass"]
    comps = _G["comps"]
    cands = list(root.glob(f"*.{station}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"))
    if not cands:
        return day, {}, 0.0
    try:
        st_all = read(str(cands[0]))
    except Exception:
        return day, {}, 0.0
    t0 = UTCDateTime(datetime(day.year, day.month, day.day))
    out = {}
    for c in comps:
        st = st_all.select(component=c).copy()
        if len(st) == 0:
            continue
        bad = False
        for tr in st:
            fs0 = float(tr.stats.sampling_rate)
            if abs(fs0 - fs) > 1e-6:
                fr = Fraction(int(round(fs)), int(round(fs0))).limit_denominator(1000)
                try:
                    tr.data = resample_poly(tr.data.astype(np.float64),
                                            fr.numerator, fr.denominator).astype(np.float32)
                    tr.stats.sampling_rate = float(fs)
                except Exception:
                    bad = True
                    break
        if bad:
            continue
        try:
            st.merge(fill_value=0)
        except Exception:
            continue
        if len(st) == 0:
            continue
        st.detrend("demean")
        try:
            st.filter("bandpass", freqmin=bp[0], freqmax=bp[1], corners=4, zerophase=True)
        except Exception:
            continue
        try:
            st.trim(t0, t0 + 86400.0, pad=True, fill_value=0, nearest_sample=False)
        except Exception:
            continue
        tr = st[0]
        data = tr.data.astype(np.float32)
        dm = _G.get("despike_mad", 0.0)
        if dm and dm > 0:
            x = data.astype(np.float64)
            med = np.median(x)
            mad = np.median(np.abs(x - med)) * 1.4826
            if mad > 0:
                np.clip(x, med - dm * mad, med + dm * mad, out=x)
                data = x.astype(np.float32)
        out[c] = data
    return day, out, float(t0.timestamp)


def _dedup_min_gap(ti, si, vals, n_gap):
    """densify_gnw_gpu._dedup_min_gap, verbatim (greedy windowed argmax per template)."""
    out = {}
    if ti.size == 0:
        return out
    btmpl = np.unique(ti)
    starts = np.searchsorted(ti, btmpl, side="left")
    ends = np.searchsorted(ti, btmpl, side="right")
    for t, a, b in zip(btmpl.tolist(), starts.tolist(), ends.tolist()):
        s = si[a:b]
        c = vals[a:b]
        picks = []
        i = 0
        L = s.size
        while i < L:
            j = s[i]
            w_end = j + n_gap
            best = i
            k = i + 1
            while k < L and s[k] < w_end:
                if c[k] > c[best]:
                    best = k
                k += 1
            picks.append(best)
            nxt = s[best] + n_gap
            while i < L and s[i] < nxt:
                i += 1
        pi = np.asarray(picks)
        out[int(t)] = (s[pi], c[pi])
    return out


def load_templates(npz3c, npzz, summary, comps, min_snr=0.0, limit=0, family_list=None):
    """Return (fam_ids, {comp: (F,m) float32 normalised+reversed}, has[comp] (F,) float32)."""
    dz = np.load(npzz, allow_pickle=True)
    d3 = np.load(npz3c, allow_pickle=True) if npz3c and Path(npz3c).exists() else None
    s = pd.read_csv(summary)
    if min_snr > 0:
        s = s[s.snr >= min_snr]
    if family_list:
        keep = [ln.strip() for ln in open(family_list) if ln.strip()]
        s = s[s.family_id.isin(set(keep))]
    fams = [f for f in s.family_id.tolist() if f in dz.files]
    if limit:
        fams = fams[:limit]
    m = int(np.asarray(dz[fams[0]]).size)
    T = {c: np.zeros((len(fams), m), np.float32) for c in comps}
    has = {c: np.zeros(len(fams), np.float32) for c in comps}

    def _norm_rev(x):
        x = np.asarray(x, np.float64)
        x = x - x.mean()
        n = np.linalg.norm(x)
        return ((x / n) if n > 0 else x)[::-1].astype(np.float32)

    for i, f in enumerate(fams):
        for c in comps:
            arr = None
            if c == "Z":
                # ALWAYS the project's canonical stored Z template (built from ALL the family's
                # members): it keeps the old-vs-new comparison like-for-like on the Z channel, and
                # the rebuilt Z is only needed where no stored one exists. Horizontals come from
                # the 3-C rebuild.
                arr = dz[f] if f in dz.files else (
                    d3[f"{f}|{c}"] if (d3 is not None and f"{f}|{c}" in d3.files) else None)
            elif d3 is not None and f"{f}|{c}" in d3.files:
                arr = d3[f"{f}|{c}"]
            if arr is None or np.asarray(arr).size != m:
                continue
            v = _norm_rev(arr)
            if not np.all(np.isfinite(v)):
                continue
            T[c][i] = v
            has[c][i] = 1.0
    return fams, T, has, m


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--net", required=True)
    p.add_argument("--sta", required=True)
    p.add_argument("--templates-z", required=True)
    p.add_argument("--templates-3c", default="")
    p.add_argument("--summary", required=True)
    p.add_argument("--wfdir", default="lfe_catalog_repair/traces")
    p.add_argument("--comps", default="Z,N,E", help="component letters, Z first (B943: 'Z,1,2')")
    p.add_argument("--fs", type=float, default=40.0)
    p.add_argument("--fmin", type=float, default=2.0)
    p.add_argument("--fmax", type=float, default=8.0)
    p.add_argument("--start", required=True)
    p.add_argument("--end", required=True)
    p.add_argument("--k-extract", type=float, default=8.0)
    p.add_argument("--old-thr", type=float, default=0.80)
    p.add_argument("--abs-floor", type=float, default=0.0,
                   help="never extract below this cc_stack (0 = no clamp)")
    p.add_argument("--min-gap-s", type=float, default=2.0)
    p.add_argument("--despike-mad", type=float, default=8.0)
    p.add_argument("--max-raw-det", type=int, default=4_000_000)
    p.add_argument("--max-det-per-tmpl-day", type=int, default=2000)
    p.add_argument("--tb", type=int, default=64, help="template batch size on the GPU")
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--min-snr", type=float, default=0.0)
    p.add_argument("--limit-templates", type=int, default=0)
    p.add_argument("--family-list", default="", help="file with one family_id per line (subset)")
    p.add_argument("--limit-days", type=int, default=0)
    p.add_argument("--out-prefix", required=True)
    args = p.parse_args()

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    import cupy as cp

    comps = [c.strip() for c in args.comps.split(",") if c.strip()]
    print(f"[gpu] {cp.cuda.runtime.getDeviceProperties(0)['name'].decode()}", flush=True)
    fams, T, has, m = load_templates(args.templates_3c, args.templates_z, args.summary,
                                     comps, args.min_snr, args.limit_templates,
                                     args.family_list or None)
    F = len(fams)
    ncomp_tmpl = np.sum([has[c] for c in comps], axis=0)
    print(f"[load] {F} templates, m={m}; per-component coverage: "
          + ", ".join(f"{c}={int(has[c].sum())}" for c in comps)
          + f"; templates with all {len(comps)} comps: {int((ncomp_tmpl == len(comps)).sum())}", flush=True)
    Tg = {c: cp.asarray(T[c]) for c in comps}
    has_g = {c: cp.asarray(has[c]) for c in comps}

    import scipy.fft as _sf
    fs = args.fs
    N_FIX = int(round(fs)) * 86400 + 16
    valid_fix = N_FIX - m + 1
    N = _sf.next_fast_len(N_FIX + m - 1)
    n_gap = max(1, int(round(args.min_gap_s * fs)))
    pre_n = int(round(10.0 * fs))
    try:
        cp.fft.config.get_plan_cache().set_size(4)   # bound the cuFFT plan cache (each plan for
    except Exception:                                 # N~3.5e6 x batch carries a large workspace)
        pass
    print(f"[gpu] N_FIX={N_FIX:,} N={N:,} tb={args.tb} n_gap={n_gap} "
          f"k_extract={args.k_extract} old_thr={args.old_thr}", flush=True)

    d0 = datetime.fromisoformat(args.start)
    d1 = datetime.fromisoformat(args.end)
    root = Path(args.wfdir)
    days = []
    d = d0
    while d < d1:
        f = root / f"{args.net}.{args.sta}" / f"{d.year}" / f"{d.timetuple().tm_yday:03d}.mseed"
        if f.exists():
            days.append(d)
        d += timedelta(days=1)
    if args.limit_days:
        days = days[:args.limit_days]
    print(f"[plan] {len(days)} day-files {args.start}..{args.end}", flush=True)

    det_rows = []      # arrays appended per day
    daily_rows = []
    t_start = time.time()
    ctx = mp.get_context("spawn")
    pool = ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx,
                               initializer=_init_worker,
                               initargs=(str(root), args.sta, fs, (args.fmin, args.fmax),
                                         args.despike_mad, comps))
    n_days_done = 0
    n_skipped = 0
    tb_state = [args.tb]      # shrinks permanently on a GPU OOM (see the retry below)
    try:
        for day, data, t0ep in pool.map(_preprocess, days, chunksize=1):
            if not data or "Z" not in data:
                n_skipped += 1
                continue
            comps_day = [c for c in comps if c in data]
            n_real = min(int(data[c].size) for c in comps_day)
            if n_real < m + pre_n:
                n_skipped += 1
                continue
            real_valid = max(0, min(n_real, N_FIX) - m + 1)
            D = {}
            ws = {}
            for c in comps_day:
                arr = data[c]
                if arr.size >= N_FIX:
                    dg = cp.asarray(arr[:N_FIX])
                else:
                    dg = cp.zeros(N_FIX, dtype=cp.float32)
                    dg[:arr.size] = cp.asarray(arr)
                D[c] = cp.fft.rfft(dg, n=N).astype(cp.complex64)
                x = dg.astype(cp.float64)
                c1 = cp.concatenate((cp.zeros(1), cp.cumsum(x)))
                wm = (c1[m:] - c1[:-m]) / m
                c2 = cp.concatenate((cp.zeros(1), cp.cumsum(x * x)))
                ws[c] = cp.sqrt(cp.maximum((c2[m:] - c2[:-m]) - m * wm * wm, 1e-12)).astype(cp.float32)
                del dg, x, c1, c2, wm

            zdat = data["Z"].astype(np.float64)
            zc2 = np.concatenate(([0.0], np.cumsum(zdat * zdat)))
            zabs = np.abs(data["Z"])

            tb_cur = tb_state[0]
            while True:
              try:
                med_all = np.zeros(F, np.float32)
                mad_all = np.zeros(F, np.float32)
                nc_all = np.zeros(F, np.float32)
                ti_a, si_a, v_a, cc_a, nc_a = [], [], [], {c: [] for c in comps}, []
                raw_cross = 0
                for b0 in range(0, F, tb_cur):
                    b1 = min(F, b0 + tb_cur)
                    nb = b1 - b0
                    cc_sum = cp.zeros((nb, valid_fix), cp.float32)
                    cnt = cp.zeros(nb, cp.float32)
                    ccs = {}
                    for c in comps_day:
                        if float(has[c][b0:b1].sum()) == 0.0:
                            continue
                        Tk = cp.fft.rfft(Tg[c][b0:b1], n=N, axis=1).astype(cp.complex64)
                        Tk *= D[c][None, :]
                        conv = cp.fft.irfft(Tk, n=N, axis=1)
                        del Tk
                        cc = conv[:, m - 1:m - 1 + valid_fix]          # view into conv
                        cc /= ws[c][None, :]
                        cp.copyto(cc, cp.float32(0.0),
                                  where=~(cp.isfinite(cc) & (cp.abs(cc) <= 1.0)))
                        w = has_g[c][b0:b1]
                        cc *= w[:, None]
                        cc_sum += cc
                        cnt += w
                        ccs[c] = (cc, conv, w)
                        del conv
                    cnt_safe = cp.maximum(cnt, 1.0)
                    cc_stack = cc_sum / cnt_safe[:, None]
                    del cc_sum
                    # robust noise stats on a stride-m subsample of the REAL part of the day
                    sub = cc_stack[:, :real_valid:m]
                    med = cp.median(sub, axis=1)
                    mad = cp.median(cp.abs(sub - med[:, None]), axis=1)
                    thr = med + cp.float32(args.k_extract) * mad
                    thr = cp.minimum(thr, cp.float32(args.old_thr))
                    if args.abs_floor > 0:
                        thr = cp.maximum(thr, cp.float32(args.abs_floor))
                    thr = cp.where(cnt > 0, thr, cp.float32(2.0))
                    mask = cc_stack[:, :real_valid] >= thr[:, None]
                    ti, si = cp.nonzero(mask)
                    raw_cross += int(ti.size)
                    vals = cc_stack[ti, si]
                    per = {c: (ccs[c][0][ti, si] if c in ccs else None) for c in comps}
                    ti_h = cp.asnumpy(ti).astype(np.int64) + b0
                    si_h = cp.asnumpy(si).astype(np.int64)
                    ti_a.append(ti_h)
                    si_a.append(si_h)
                    v_a.append(cp.asnumpy(vals).astype(np.float32))
                    for c in comps:
                        if per[c] is None:
                            cc_a[c].append(np.full(ti_h.size, np.nan, np.float32))
                        else:
                            cc_a[c].append(cp.asnumpy(per[c]).astype(np.float32))
                    nc_a.append(cp.asnumpy(cnt)[cp.asnumpy(ti).astype(np.int64)].astype(np.int8)
                                if ti.size else np.empty(0, np.int8))
                    med_all[b0:b1] = cp.asnumpy(med)
                    mad_all[b0:b1] = cp.asnumpy(mad)
                    nc_all[b0:b1] = cp.asnumpy(cnt)
                    del cc_stack, sub, med, mad, thr, mask, ti, si, vals, ccs, cnt, cnt_safe
                for c in comps_day:
                    del D[c], ws[c]
                cp.get_default_memory_pool().free_all_blocks()

                break
              except cp.cuda.memory.OutOfMemoryError:
                  # drop every big device array the failed batch may still be holding via the
                  # frame locals, THEN release the pool, THEN halve the batch.
                  for _n in ("cc_sum", "cc", "conv", "Tk", "cc_stack", "sub", "mask", "ti", "si",
                             "vals", "per", "ccs", "cnt", "cnt_safe", "med", "mad", "thr"):
                      if _n in dir():
                          try:
                              del locals()[_n]
                          except Exception:
                              pass
                  ccs = per = None
                  cc_sum = cc = conv = Tk = cc_stack = sub = mask = None
                  ti = si = vals = cnt = cnt_safe = med = mad = thr = None
                  import gc as _gc
                  _gc.collect()
                  cp.get_default_memory_pool().free_all_blocks()
                  try:
                      cp.fft.config.get_plan_cache().clear()
                  except Exception:
                      pass
                  if tb_cur <= 8:
                      raise
                  tb_cur = max(8, tb_cur // 2)
                  tb_state[0] = tb_cur
                  print(f"  [oom] {day.date()}: retrying with tb={tb_cur}", flush=True)
            if args.max_raw_det and raw_cross > args.max_raw_det:
                daily_rows.append(pd.DataFrame(dict(
                    day=np.full(F, day.strftime("%Y-%m-%d")), tmpl=np.arange(F),
                    med=med_all, mad=mad_all, n_comp=nc_all,
                    n_det=np.zeros(F, np.int32), capped=np.zeros(F, np.int8),
                    n_valid=np.full(F, real_valid, np.int64),
                    qc_drop=np.ones(F, np.int8))))
                n_skipped += 1
                n_days_done += 1
                continue

            ti = np.concatenate(ti_a) if ti_a else np.empty(0, np.int64)
            si = np.concatenate(si_a) if si_a else np.empty(0, np.int64)
            vv = np.concatenate(v_a) if v_a else np.empty(0, np.float32)
            ccc = {c: (np.concatenate(cc_a[c]) if cc_a[c] else np.empty(0, np.float32)) for c in comps}
            ncv = np.concatenate(nc_a) if nc_a else np.empty(0, np.int8)
            res = _dedup_min_gap(ti, si, vv, n_gap)
            # map declustered picks back to the concatenated arrays: rebuild index
            n_det = np.zeros(F, np.int32)
            capped = np.zeros(F, np.int8)
            keep_idx = []
            order = np.lexsort((si, ti))
            ti_s, si_s = ti[order], si[order]
            for t, (samp, ccv) in res.items():
                a = np.searchsorted(ti_s, t, side="left")
                b = np.searchsorted(ti_s, t, side="right")
                loc = a + np.searchsorted(si_s[a:b], samp)
                if args.max_det_per_tmpl_day and loc.size > args.max_det_per_tmpl_day:
                    sel = np.argpartition(ccv, -args.max_det_per_tmpl_day)[-args.max_det_per_tmpl_day:]
                    loc = loc[sel]
                    capped[t] = 1
                n_det[t] = loc.size
                keep_idx.append(order[loc])
            if keep_idx:
                ki = np.concatenate(keep_idx)
                s_i = si[ki]
                lo = np.maximum(s_i - pre_n, 0)
                cnt_pre = (s_i - lo).astype(np.float64)
                pre_ms = np.where(cnt_pre > 0, (zc2[s_i] - zc2[lo]) / np.maximum(cnt_pre, 1), np.nan)
                pre_rms = np.sqrt(pre_ms)
                gather = s_i[:, None] + np.arange(m)[None, :]
                np.clip(gather, 0, zabs.size - 1, out=gather)
                peak = zabs[gather].max(axis=1)
                snr = np.where(pre_rms > 0, peak / pre_rms, np.nan).astype(np.float32)
                det_rows.append(pd.DataFrame(dict(
                    tmpl=ti[ki].astype(np.int32),
                    t=(t0ep + s_i / fs).astype(np.float64),
                    samp=s_i.astype(np.int64),
                    cc_stack=vv[ki],
                    **{f"cc_{c}": ccc[c][ki] for c in comps},
                    n_comp=ncv[ki],
                    snr=snr,
                    med=med_all[ti[ki]],
                    mad=mad_all[ti[ki]],
                )))
            daily_rows.append(pd.DataFrame(dict(
                day=np.full(F, day.strftime("%Y-%m-%d")), tmpl=np.arange(F),
                med=med_all, mad=mad_all, n_comp=nc_all, n_det=n_det, capped=capped,
                n_valid=np.full(F, real_valid, np.int64), qc_drop=np.zeros(F, np.int8))))
            n_days_done += 1
            if n_days_done % 25 == 0:
                el = time.time() - t_start
                nd = sum(len(x) for x in det_rows)
                print(f"  day {n_days_done}/{len(days)} ({el/60:.1f} min, "
                      f"{n_days_done/max(el,1e-9):.2f} d/s) dets so far {nd:,}", flush=True)
    finally:
        pool.shutdown(wait=True)

    out = Path(args.out_prefix)
    out.parent.mkdir(parents=True, exist_ok=True)
    det = (pd.concat(det_rows, ignore_index=True) if det_rows
           else pd.DataFrame(columns=["tmpl", "t", "cc_stack"]))
    det["template"] = pd.Categorical.from_codes(det["tmpl"].values, categories=fams) if len(det) else []
    det.to_parquet(f"{out}_detections.parquet", index=False)
    dly = pd.concat(daily_rows, ignore_index=True)
    dly["template"] = pd.Categorical.from_codes(dly["tmpl"].values, categories=fams)
    dly.to_parquet(f"{out}_daily.parquet", index=False)
    el = time.time() - t_start
    print(f"[DONE] {n_days_done} days processed, {n_skipped} skipped, {len(det):,} detections in "
          f"{el/60:.1f} min ({n_days_done/max(el/3600,1e-9)/365:.2f} station-years/hour)", flush=True)
    print(f"  -> {out}_detections.parquet ({os.path.getsize(f'{out}_detections.parquet')/1e6:.1f} MB)")
    print(f"  -> {out}_daily.parquet")


if __name__ == "__main__":
    main()
