#!/usr/bin/env python
"""D3 stage 5 — NETWORK matched filter (the fix the sensitivity budget calls for).

Stage 3 measured the budget: a single median Lin/Bostock LFE at PGC has amplitude s/n = 0.078 of the
ambient 1-8 Hz noise, while the single-channel matched-filter cc noise is sigma = 1/sqrt(2*B*L)
(0.189 at L = 2 s, 0.109 at L = 6 s). One 3-component station therefore sits at ~1-2 MAD of signal -
undetectable at any threshold. Summing the cc over Nch coherently aligned channels divides sigma by
sqrt(Nch), so the network sum is the only way in.

TEMPLATE TIMES. The stage-4 pre-check showed that stacking at the PGC *re-detections* yields nothing
(they are noise), so network templates cannot be built from them. This script builds them from the
project's **refined Lin/Bostock picks** (`lfe_features/data/lin_precise_b011.csv`, cross-correlation
refined at B011, 0.3 km from PGC) grouped by source cell - the only times in this pilot known to be
real LFEs. Half of the picks build the templates, the other half is held out to score recall.

WHAT IT DOES
  1. per (cell, station, component): stack the band-passed 40 Hz waveform at the FIT-half refined
     arrival times, corrected for the station's predicted S moveout from the cell (1-D model), then
     cut a template of length L starting 0.5 s before the S arrival on the stack;
     stations are kept on stack SNR against a same-n random-time null;
  2. GPU matched filter for every kept (cell, station, component) over the whole window;
  3. shift each cc series back by that station's moveout and SUM -> the network cc, then
     median + k*MAD per 30-day window of the NETWORK sum, k scanned over 4,5,6,8,9,10;
  4. decluster at the template length; store network cc, contributing channel count and per-station
     cc for every detection;
  5. score recall/precision against the HELD-OUT refined picks (+-0.5 s), against a +600 s shifted
     control, and report the fraction of detections seen at >= 3 stations.
Repeat for L = 2 / 4 / 6 s to test template length.

Outputs
  results/d3_stage5_templates_L<L>.npz / .csv
  results/d3_stage5_detections_L<L>.parquet
  results/d3_stage5_summary.csv
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
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, str(Path(__file__).resolve().parent))
import d3_vel1d as V                                              # noqa: E402
from d3_stage4_stacks import _init, _load_day, align_stack, STA_LL, FS   # noqa: E402
from d3_redetect import _dedup_min_gap                            # noqa: E402

FMIN, FMAX = 1.0, 8.0
# 22/12 s deliberately matches d3_stage4_stacks so that d3_stage4_locate.py can consume these
# stacks unchanged (it imports NPRE from there and expects sample NPRE = the S arrival).
PRE_S, POST_S = 22.0, 12.0          # stacking window around the (moveout-corrected) S arrival
NPRE = int(round(PRE_S * FS))
NWIN = int(round((PRE_S + POST_S) * FS))
NOISE = int(round(8.0 * FS))
TPL_LEAD = 0.5                      # template starts this long before the picked S
KS = (4.0, 5.0, 6.0, 8.0, 9.0, 10.0)


def snr_of(arr, i0, i1):
    from scipy.signal import hilbert
    env = np.abs(hilbert(arr - arr.mean()))
    rms = float(np.sqrt(np.mean(arr[:NOISE] ** 2)))
    i0, i1 = max(i0, NOISE), min(i1, len(arr) - 1)
    if not (rms > 0 and i1 > i0):
        return np.nan, -1
    return float(env[i0:i1].max() / rms), int(i0 + np.argmax(env[i0:i1]))


def cell_key(lat, lon, deg):
    return f"{round(lat/deg)*deg:.2f}_{round(lon/deg)*deg:.2f}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--picks", default="lfe_catalog_repair/results/d3_lin_refined_window.csv")
    ap.add_argument("--stations", required=True, help="NET:STA,...")
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--cell-deg", type=float, default=0.1)
    ap.add_argument("--min-picks", type=int, default=40)
    ap.add_argument("--lengths", default="2,4,6")
    ap.add_argument("--min-stack-snr", type=float, default=0.0,
                    help="absolute stack-SNR floor; 0 = use the null-ratio gate only")
    ap.add_argument("--snr-ratio", type=float, default=1.3, help="stack SNR / null SNR gate")
    ap.add_argument("--single-cell", action="store_true",
                    help="treat ALL picks as ONE source at their centroid, correcting each event "
                         "by ITS OWN predicted moveout. The 0.1-deg cells give only ~35 picks in "
                         "the fit half, and at the measured per-event s/n = 0.078 a 35-event stack "
                         "reaches just sqrt(1+35*0.0061) = 1.10x the noise -- far too shallow. "
                         "Per-event correction lets the whole fit half (388 picks -> 1.83x) build "
                         "one template per station/component.")
    ap.add_argument("--depth", type=float, default=35.0,
                    help="assumed source depth (km) for the moveout; Lin family depths are unreliable")
    ap.add_argument("--tol", type=float, default=1.0,
                    help="recall match tolerance (s); 1.0 in --single-cell mode absorbs the "
                         "moveout error of events away from the centroid")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--tb", type=int, default=64)
    ap.add_argument("--use-all", action="store_true",
                    help="ignore the per-station stack gate and sum EVERY station. The gate exists "
                         "to keep pure-noise channels out, but the measured stack ratios (1.05-1.19 "
                         "at the neighbours vs 2.3 at PGC) say they carry weak but non-zero signal, "
                         "and the point of the network sum is exactly to add those up.")
    ap.add_argument("--min-net-sta", type=int, default=3)
    ap.add_argument("--outdir", default="lfe_catalog_repair/results")
    a = ap.parse_args()
    rng = np.random.default_rng(0)

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    import cupy as cp
    import scipy.fft as _sf

    stns = [tuple(s.split(":")) for s in a.stations.split(",") if s.strip()]
    d0, d1 = datetime.fromisoformat(a.start), datetime.fromisoformat(a.end)
    P = pd.read_csv(a.picks)
    P["cell"] = [cell_key(la, lo, a.cell_deg) for la, lo in zip(P.lat, P.lon)]
    cnt = P.groupby("cell").size().sort_values(ascending=False)
    cells = cnt[cnt >= a.min_picks].index.tolist()
    print(f"[s5] {len(P)} refined picks -> {len(cells)} cells with >= {a.min_picks} "
          f"({cnt[cells].sum()} picks); stations {[s for _, s in stns]}", flush=True)
    if not cells:
        print("[s5] no cell reaches the pick minimum -> stop")
        return
    if a.single_cell:
        cells = ["ALL"]
        P = P.copy()
        P["cell"] = "ALL"
    fit, held, fit_ll, held_ll = {}, {}, {}, {}
    for c in cells:
        g = P[P.cell == c].copy()
        g = g[(g.arrival > d0.timestamp() + 60) & (g.arrival < d1.timestamp() - 60)]
        g = g.sort_values("arrival")
        fit[c] = g.arrival.values[0::2].astype(float)
        held[c] = g.arrival.values[1::2].astype(float)
        fit_ll[c] = (g.lat.values[0::2], g.lon.values[0::2])
        held_ll[c] = (g.lat.values[1::2], g.lon.values[1::2])
    cell_ll = {c: (float(P.loc[P.cell == c, "lat"].mean()),
                   float(P.loc[P.cell == c, "lon"].mean())) for c in cells}
    print(f"[s5] fit/held per cell: "
          + ", ".join(f"{c}:{len(fit[c])}/{len(held[c])}" for c in cells), flush=True)

    def moveout_of(sta, la, lo, ref_sta, depth):
        """S-time difference (station - reference) for a source at (la, lo, depth), per EVENT."""
        la = np.atleast_1d(np.asarray(la, float))
        lo = np.atleast_1d(np.asarray(lo, float))
        out = np.empty(la.size)
        for i in range(la.size):
            xr = float(np.hypot((la[i] - STA_LL[ref_sta][0]) * 111.2,
                                (lo[i] - STA_LL[ref_sta][1]) * 111.2 * np.cos(np.radians(la[i]))))
            x = float(np.hypot((la[i] - STA_LL[sta][0]) * 111.2,
                               (lo[i] - STA_LL[sta][1]) * 111.2 * np.cos(np.radians(la[i]))))
            out[i] = V.tt(x, depth, "S") - V.tt(xr, depth, "S")
        return out
    # station moveout relative to the reference station (first in the list)
    ref = stns[0][1]
    mo = {}
    for c in cells:
        la, lo = cell_ll[c]
        xr = float(np.hypot((la - STA_LL[ref][0]) * 111.2,
                            (lo - STA_LL[ref][1]) * 111.2 * np.cos(np.radians(la))))
        tref = V.tt(xr, a.depth, "S")
        for _, sta in stns:
            x = float(np.hypot((la - STA_LL[sta][0]) * 111.2,
                               (lo - STA_LL[sta][1]) * 111.2 * np.cos(np.radians(la))))
            mo[(c, sta)] = V.tt(x, a.depth, "S") - tref
    print(f"[s5] moveout vs {ref}: "
          + ", ".join(f"{s}={np.median([mo[(c,s)] for c in cells]):+.2f}s" for _, s in stns), flush=True)

    days = []
    d = d0
    while d < d1:
        days.append((d.year, d.month, d.day))
        d += timedelta(days=1)

    # ---------------- pass 1: build the network templates ----------------
    print("[s5] pass 1: stacking the network templates at the FIT-half refined picks", flush=True)
    acc, accn, ncomp_of, ev_off = {}, {}, {}, {}
    items = [(net, sta, y, m, dd) for (net, sta) in stns for (y, m, dd) in days]
    ctx = mp.get_context("spawn")
    t0 = time.time()
    with ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init,
                             initargs=("lfe_catalog_repair/traces", FS, (FMIN, FMAX))) as ex:
        for (net, sta, ymd), data in ex.map(_load_day, items, chunksize=1):
            if not data:
                continue
            comps = [c for c in ["Z"] + [x for x in sorted(data) if x != "Z"] if c in data]
            ncomp_of[sta] = len(comps)
            n = min(data[c].size for c in comps)
            X = np.stack([data[c][:n] for c in comps], axis=1)
            day_ep = datetime(*ymd).timestamp()
            for c in cells:
                key = (c, sta)
                if key not in ev_off:
                    ev_off[key] = moveout_of(sta, fit_ll[c][0], fit_ll[c][1], ref, a.depth)
                for src, off, store in ((fit[c], ev_off[key], acc),
                                        (None, mo[(c, sta)], accn)):
                    if src is None:
                        tt_ = rng.uniform(day_ep + 100, day_ep + 86300,
                                          max(len(fit[c]) // len(days), 2))
                        oo = off
                    else:
                        tt_, oo = src, off
                    idx = np.round((tt_ + oo - day_ep) * FS).astype(np.int64) - NPRE
                    ok = (idx >= 0) & (idx + NWIN <= n)
                    if not np.any(ok):
                        continue
                    g = idx[ok][:, None] + np.arange(NWIN)[None, :]
                    store.setdefault(key, []).append(X[g])
    print(f"[s5]   stacked in {(time.time()-t0)/60:.1f} min", flush=True)

    stack_rows, stacks = [], {}
    for (c, sta), parts in sorted(acc.items()):
        cuts = np.concatenate(parts, axis=0)
        if cuts.shape[0] < 10:
            continue
        amp = np.sqrt((cuts ** 2).mean(axis=(1, 2)))
        amp[amp <= 0] = 1.0
        st = (cuts / amp[:, None, None]).mean(0)           # RAW stack: no alignment refinement
        h = np.sqrt((st[:, 1:] ** 2).sum(1)) if st.shape[1] > 1 else st[:, 0]
        s_h, i_h = snr_of(h, NPRE - int(2 * FS), NPRE + int(4 * FS))
        nsh = np.nan
        if (c, sta) in accn:
            nc = np.concatenate(accn[(c, sta)], axis=0)[:len(cuts)]
            amp = np.sqrt((nc ** 2).mean(axis=(1, 2)))
            amp[amp <= 0] = 1.0
            nst = (nc / amp[:, None, None]).mean(0)
            nh = np.sqrt((nst[:, 1:] ** 2).sum(1)) if nst.shape[1] > 1 else nst[:, 0]
            nsh = snr_of(nh, NPRE - int(2 * FS), NPRE + int(4 * FS))[0]
        keep = bool(s_h >= a.min_stack_snr
                    and (not np.isfinite(nsh) or s_h > a.snr_ratio * nsh))
        s_z, _iz = snr_of(st[:, 0], NPRE - int(6 * FS), NPRE + int(4 * FS))
        stacks[(c, sta)] = (st, i_h)
        stack_rows.append(dict(family=c, cell=c, station=sta, n_cuts=int(cuts.shape[0]),
                               n_comp=st.shape[1], snr_h=s_h, snr_z=s_z, null_snr_h=nsh,
                               null_snr_z=np.nan, cc_mean=np.nan, i_S=i_h,
                               pick_offset_s=(i_h - NPRE) / FS, moveout_s=mo[(c, sta)],
                               dt_pred_s=mo[(c, sta)], keep=int(keep)))
    S = pd.DataFrame(stack_rows)
    S.to_csv(f"{a.outdir}/d3_stage5_stacks.csv", index=False)
    # stage-4-compatible dump so d3_stage4_locate.py can pick and locate these (real) stacks
    _out = {}
    for (c, sta), (st, iS) in stacks.items():
        for j, cn in enumerate(["Z", "H1", "H2"][:st.shape[1]]):
            _out[f"{c}|{sta}|{cn}"] = st[:, j]
    np.savez(f"{a.outdir}/d3_stage5_stacks.npz", **_out)
    print(f"[s5] {len(S)} cell-station stacks, {int(S.keep.sum())} kept "
          f"(snr_h >= {a.min_stack_snr} and > {a.snr_ratio}x null)")
    for sta, g in S.groupby("station"):
        print(f"    {sta:5s} keep {int(g.keep.sum()):2d}/{len(g):2d}  snr_h med {g.snr_h.median():.2f} "
              f"(null {g.null_snr_h.median():.2f})  pick offset med {g.pick_offset_s.median():+.2f} s")
    kept = S if a.use_all else S[S.keep == 1]
    if a.use_all:
        print(f"[s5] --use-all: summing all {len(kept)} cell-station stacks regardless of the gate")
    ncell = kept.groupby("cell").station.nunique()
    use_cells = ncell[ncell >= a.min_net_sta].index.tolist()
    print(f"[s5] cells with >= {a.min_net_sta} usable stations: {len(use_cells)} / {len(cells)}")
    if not use_cells:
        print("[s5] no cell has a usable network template -> stop (reported as such)")
        S.to_csv(f"{a.outdir}/d3_stage5_stacks.csv", index=False)
        return

    # ---------------- pass 2: network matched filter, per template length ----------------
    summary = []
    N_FIX = int(round(FS)) * 86400 + 16
    for L in [float(x) for x in a.lengths.split(",")]:
        m = int(round(L * FS))
        lead = int(round(TPL_LEAD * FS))
        chans, T = [], []
        for _, r in kept.iterrows():
            if r.cell not in use_cells:
                continue
            st, iS = stacks[(r.cell, r.station)]
            i0 = int(iS) - lead
            if i0 < 0 or i0 + m > st.shape[0]:
                continue
            for j in range(st.shape[1]):
                v = st[i0:i0 + m, j].astype(np.float64)
                v = v - v.mean()
                nn = np.linalg.norm(v)
                if nn <= 0:
                    continue
                chans.append((r.cell, r.station, j, mo[(r.cell, r.station)]))
                T.append((v / nn)[::-1].astype(np.float32))
        if not chans:
            print(f"[s5] L={L}: no channel templates -> skip")
            continue
        T = np.asarray(T, np.float32)
        Tg = cp.asarray(T)
        valid = N_FIX - m + 1
        N = _sf.next_fast_len(N_FIX + m - 1)
        by_sta = {}
        for i, (c, sta, j, off) in enumerate(chans):
            by_sta.setdefault(sta, []).append(i)
        cell_idx = {}
        for i, (c, sta, j, off) in enumerate(chans):
            cell_idx.setdefault(c, []).append(i)
        print(f"[s5] L={L:.0f} s (m={m}): {len(chans)} channel templates over "
              f"{len(cell_idx)} cells, {len(by_sta)} stations; "
              f"predicted sigma_1ch = {1/np.sqrt(2*(FMAX-FMIN)*L):.3f}", flush=True)

        det_rows, daily_rows = [], []
        ctx = mp.get_context("spawn")
        with ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init,
                                 initargs=("lfe_catalog_repair/traces", FS, (FMIN, FMAX))) as ex:
            per_day = {}
            for (net, sta, ymd), data in ex.map(_load_day, items, chunksize=1):
                per_day.setdefault(ymd, {})[sta] = data
                if len(per_day[ymd]) < len(stns):
                    continue
                dd = per_day.pop(ymd)
                day_ep = datetime(*ymd).timestamp()
                # per-channel cc on the GPU
                cc_all = {}
                for sta, data in dd.items():
                    if not data or sta not in by_sta:
                        continue
                    comps = [c for c in ["Z"] + [x for x in sorted(data) if x != "Z"] if c in data]
                    for j, c in enumerate(comps):
                        arr = data[c]
                        if arr.size >= N_FIX:
                            dg = cp.asarray(arr[:N_FIX])
                        else:
                            dg = cp.zeros(N_FIX, cp.float32)
                            dg[:arr.size] = cp.asarray(arr)
                        D = cp.fft.rfft(dg, n=N).astype(cp.complex64)
                        x = dg.astype(cp.float64)
                        c1 = cp.concatenate((cp.zeros(1), cp.cumsum(x)))
                        wm = (c1[m:] - c1[:-m]) / m
                        c2 = cp.concatenate((cp.zeros(1), cp.cumsum(x * x)))
                        ws = cp.sqrt(cp.maximum((c2[m:] - c2[:-m]) - m * wm * wm, 1e-12)).astype(cp.float32)
                        rows = [i for i in by_sta[sta] if chans[i][2] == j]
                        for b0 in range(0, len(rows), a.tb):
                            rr = rows[b0:b0 + a.tb]
                            Tk = cp.fft.rfft(Tg[rr], n=N, axis=1).astype(cp.complex64)
                            Tk *= D[None, :]
                            conv = cp.fft.irfft(Tk, n=N, axis=1)
                            del Tk
                            cc = conv[:, m - 1:m - 1 + valid]
                            cc /= ws[None, :]
                            cp.copyto(cc, cp.float32(0.0),
                                      where=~(cp.isfinite(cc) & (cp.abs(cc) <= 1.0)))
                            for kk, i in enumerate(rr):
                                cc_all[i] = cc[kk].copy()
                            del conv, cc
                        del dg, D, x, c1, c2, ws, wm
                        cp.get_default_memory_pool().free_all_blocks()
                # network sum per cell, with the moveout shift
                for c, idxs in cell_idx.items():
                    have = [i for i in idxs if i in cc_all]
                    if len(have) < 3:
                        continue
                    net_sum = cp.zeros(valid, cp.float32)
                    nch = 0
                    for i in have:
                        off = chans[i][3]
                        sh = int(round(off * FS))
                        src = cc_all[i]
                        if sh >= 0:
                            net_sum[:valid - sh] += src[sh:]
                        else:
                            net_sum[-sh:] += src[:valid + sh]
                        nch += 1
                    net_sum /= nch
                    sub = net_sum[::m]
                    med = float(cp.median(sub))
                    mad = float(cp.median(cp.abs(sub - med)))
                    thr = med + min(KS) * mad
                    mask = net_sum >= thr
                    ti = cp.nonzero(mask)[0]
                    if ti.size:
                        si = cp.asnumpy(ti).astype(np.int64)
                        vv = cp.asnumpy(net_sum[ti]).astype(np.float32)
                        res = _dedup_min_gap(np.zeros(si.size, np.int64), si, vv, m)
                        if 0 in res:
                            samp, ccv = res[0]
                            per = {}
                            for i in have:
                                sta = chans[i][1]
                                sh = int(round(chans[i][3] * FS))
                                ii = np.clip(samp + sh, 0, valid - 1)
                                per.setdefault(sta, []).append(cp.asnumpy(cc_all[i][cp.asarray(ii)]))
                            pst = {s_: np.mean(v, axis=0) for s_, v in per.items()}
                            det_rows.append(pd.DataFrame(dict(
                                cell=c, t=day_ep + samp / FS, cc_net=ccv, n_chan=nch,
                                med=med, mad=mad,
                                n_sta_above=np.sum([(pst[s_] > 0.25).astype(int) for s_ in pst], axis=0),
                                **{f"cc_{s_}": pst[s_] for s_ in pst})))
                    daily_rows.append(dict(cell=c, day=datetime(*ymd).strftime("%Y-%m-%d"),
                                           med=med, mad=mad, n_chan=nch))
                    del net_sum, sub
                for i in list(cc_all):
                    del cc_all[i]
                cp.get_default_memory_pool().free_all_blocks()

        if not det_rows:
            print(f"[s5] L={L}: no detections")
            continue
        Dt = pd.concat(det_rows, ignore_index=True)
        Dl = pd.DataFrame(daily_rows)
        Dl.to_csv(f"{a.outdir}/d3_stage5_daily_L{int(L)}.csv", index=False)
        # window thresholds (the window here is the whole run: 34 days < 30-day cadence x 2)
        w = Dl.groupby("cell").agg(medw=("med", "median"), madw=("mad", "median")).reset_index()
        Dt = Dt.merge(w, on="cell", how="left")
        for k in KS:
            Dt[f"k{int(k)}"] = Dt.cc_net >= (Dt.medw + k * Dt.madw)
        Dt.to_parquet(f"{a.outdir}/d3_stage5_detections_L{int(L)}.parquet", index=False)
        sig_meas = float(Dl.mad.median() / 0.6745)
        print(f"[s5] L={L:.0f}: {len(Dt):,} network detections; measured network cc sigma "
              f"= {sig_meas:.4f} (MAD {Dl.mad.median():.4f}), "
              f"predicted {1/np.sqrt(2*(FMAX-FMIN)*L)/np.sqrt(Dl.n_chan.median()):.4f} "
              f"at {int(Dl.n_chan.median())} channels")
        for k in KS:
            s = Dt[Dt[f"k{int(k)}"]]
            rec = prec = recs = np.nan
            if len(s):
                allheld = np.concatenate([held[c] for c in use_cells if len(held[c])])
                td = s.t.values
                # per-cell matching (each cell has its own held-out picks)
                nm = nd = 0
                nl = 0
                for c, g in s.groupby("cell"):
                    hv = np.sort(held.get(c, np.array([])))
                    if not len(hv):
                        continue
                    x = np.sort(g.t.values)
                    i = np.searchsorted(hv, x)
                    dd_ = np.full(len(x), np.inf)
                    ok = i < len(hv)
                    dd_[ok] = np.abs(hv[i[ok]] - x[ok])
                    ok2 = i > 0
                    dd_[ok2] = np.minimum(dd_[ok2], np.abs(hv[i[ok2] - 1] - x[ok2]))
                    nd += len(x)
                    nm += int((dd_ <= a.tol).sum())
                    j = np.searchsorted(x, hv)
                    e = np.full(len(hv), np.inf)
                    ok = j < len(x)
                    e[ok] = np.abs(x[j[ok]] - hv[ok])
                    ok2 = j > 0
                    e[ok2] = np.minimum(e[ok2], np.abs(x[j[ok2] - 1] - hv[ok2]))
                    nl += int((e <= a.tol).sum())
                rec = nl / max(sum(len(held.get(c, [])) for c in s.cell.unique()), 1)
                prec = nm / max(nd, 1)
                # shifted control
                nm2 = 0
                for c, g in s.groupby("cell"):
                    hv = np.sort(held.get(c, np.array([])))
                    if not len(hv):
                        continue
                    x = np.sort(g.t.values) + 600.0
                    i = np.searchsorted(hv, x)
                    dd_ = np.full(len(x), np.inf)
                    ok = i < len(hv)
                    dd_[ok] = np.abs(hv[i[ok]] - x[ok])
                    ok2 = i > 0
                    dd_[ok2] = np.minimum(dd_[ok2], np.abs(hv[i[ok2] - 1] - x[ok2]))
                    nm2 += int((dd_ <= a.tol).sum())
                recs = nm2 / max(nd, 1)
            summary.append(dict(L=L, k=k, n_det=len(s), n_cells=int(s.cell.nunique()) if len(s) else 0,
                                sigma_net=sig_meas, recall=rec, precision=prec,
                                precision_shift=recs,
                                frac_ge3sta=float((s.n_sta_above >= 3).mean()) if len(s) else np.nan))
            print(f"    k={k:4.1f}: {len(s):8,} det  recall {rec if rec==rec else float('nan'):.4f}  "
                  f"precision {prec if prec==prec else float('nan'):.4f}  "
                  f"(shifted {recs if recs==recs else float('nan'):.4f})  "
                  f"frac >=3 stations {summary[-1]['frac_ge3sta']:.3f}")
    if summary:
        pd.DataFrame(summary).to_csv(f"{a.outdir}/d3_stage5_summary.csv", index=False)
        print(f"  -> {a.outdir}/d3_stage5_summary.csv")


if __name__ == "__main__":
    main()
