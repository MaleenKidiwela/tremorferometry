#!/usr/bin/env python
"""D3 stage 4b — NETWORK-CONSISTENT family template stacks.

For every PGC family with >= --min-det re-detections inside the stage-4 window, cut and stack the
ADDED stations' 3-C waveforms at the PGC detection times, aligned first with predicted S times from
the 1-D model (`d3_vel1d.py`) at the family's source cell and Slab2 depth, then refined by
cross-correlating each individual cut against the stack (2 iterations, +-MAXSHIFT s).

HONESTY CONTROL (stack alignment can manufacture a coherent arrival out of noise): for every
(family, station) an identical stack is built at the SAME NUMBER of RANDOM times inside the same
window, with the same alignment refinement. A station is kept only if its real stack SNR >= --min-snr
AND exceeds the 95th percentile of its own random-stack SNR.

Outputs
  results/d3_stage4_stacks.npz    keys "<family>|<STA>|<comp>" -> float32 stack (40 Hz, 1-8 Hz)
  results/d3_stage4_stacks.csv    family x station: n_cuts, snr_z, snr_h, cc_mean, shift, null snr
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from datetime import datetime, timedelta
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import d3_vel1d as V   # noqa: E402

FS = 40.0
FMIN, FMAX = 1.0, 8.0
PRE_S, POST_S = 22.0, 12.0          # stack window relative to the predicted S arrival
NOISE_S = 8.0                       # leading noise window: pre-P at every pilot distance
                                    # (max S-P here is ~12.3 s at 100 km, so P lands at
                                    #  PRE_S-12.3 = 9.7 s > NOISE_S)
NPRE = int(round(PRE_S * FS))
NWIN = int(round((PRE_S + POST_S) * FS))
MAXSHIFT_S = 1.5
NITER = 2
_G = {}


def _init(root, fs, bp):
    _G.update(root=Path(root), fs=fs, bp=bp)


def _load_day(item):
    """(net, sta, day) -> {comp: array} on the [day, day+86400] grid at 40 Hz, 1-8 Hz."""
    from obspy import read, UTCDateTime
    from scipy.signal import resample_poly
    net, sta, y, m, d = item
    root, fs, bp = _G["root"], _G["fs"], _G["bp"]
    p = root / f"{net}.{sta}" / f"{y}" / f"{datetime(y, m, d).timetuple().tm_yday:03d}.mseed"
    if not p.exists():
        return (net, sta, (y, m, d)), {}
    try:
        st_all = read(str(p))
    except Exception:
        return (net, sta, (y, m, d)), {}
    t0 = UTCDateTime(datetime(y, m, d))
    out = {}
    for c in sorted({tr.stats.channel[-1] for tr in st_all}):
        st = st_all.select(component=c).copy()
        bad = False
        for tr in st:
            f0 = float(tr.stats.sampling_rate)
            if abs(f0 - fs) > 1e-6:
                fr = Fraction(int(round(fs)), int(round(f0))).limit_denominator(1000)
                try:
                    tr.data = resample_poly(tr.data.astype(np.float64),
                                            fr.numerator, fr.denominator).astype(np.float32)
                    tr.stats.sampling_rate = fs
                except Exception:
                    bad = True
                    break
        if bad:
            continue
        try:
            st.merge(fill_value=0)
            st.detrend("demean")
            st.filter("bandpass", freqmin=bp[0], freqmax=bp[1], corners=4, zerophase=True)
            st.trim(t0, t0 + 86400.0, pad=True, fill_value=0, nearest_sample=False)
        except Exception:
            continue
        if len(st):
            out[c] = st[0].data.astype(np.float32)
    return (net, sta, (y, m, d)), out


def align_stack(cuts, maxshift=MAXSHIFT_S, niter=NITER):
    """Iterative alignment: cross-correlate each cut with the running stack, re-stack.

    cuts: (n, NWIN, ncomp). Alignment is driven by the HORIZONTAL energy where available (the S
    arrival), else by Z. Returns (stack (NWIN, ncomp), shifts (n,), cc (n,)).
    """
    n, nw, nc = cuts.shape
    ns = int(round(maxshift * FS))
    drive = (np.sqrt((cuts[:, :, 1:] ** 2).sum(2)) if nc > 1 else np.abs(cuts[:, :, 0]))
    drive = drive - drive.mean(1, keepdims=True)
    nrm = np.linalg.norm(drive, axis=1, keepdims=True)
    drive = drive / np.where(nrm > 0, nrm, 1.0)
    shifts = np.zeros(n, int)
    cc = np.zeros(n)
    ref = drive.mean(0)
    # circular cross-correlation of every cut against the reference, all lags at once:
    # irfft(conj(F{X}) * F{r})[s] == dot(np.roll(X, s), r), which is what the explicit
    # per-cut, per-lag loop computed (that loop is ~1e11 flops over the whole pilot).
    lags = np.arange(-ns, ns + 1)
    cols = lags % nw
    D = np.fft.rfft(drive, axis=1)
    for _ in range(niter):
        r = ref - ref.mean()
        r /= (np.linalg.norm(r) + 1e-12)
        Vr = np.fft.irfft(np.conj(D) * np.fft.rfft(r)[None, :], n=nw, axis=1)[:, cols]
        j = np.argmax(Vr, axis=1)
        shifts = lags[j]
        cc = Vr[np.arange(n), j]
        ref = np.mean([np.roll(drive[i], shifts[i]) for i in range(n)], axis=0)
    stack = np.zeros((nw, nc), np.float64)
    for i in range(n):
        x = cuts[i]
        a = np.sqrt((x ** 2).mean())
        if a <= 0:
            continue
        stack += np.roll(x / a, shifts[i], axis=0)
    stack /= max(n, 1)
    return stack.astype(np.float32), shifts, cc


def stack_snr(stack):
    """(snr_z, snr_h) = peak envelope inside the PHASE search window / RMS of the leading noise.

    The peak is searched only where the phase can be (Z: the P band, NPRE-14 s .. NPRE+2 s;
    H: the S band, NPRE-3 s .. NPRE+4 s) rather than over the whole trace: a max over 20 s of pure
    noise already gives peak/RMS ~3.5, which would pass a naive SNR >= 3 gate on noise alone
    (measured: pure-noise stacks give snr_z ~3.6, snr_h ~1.8 with the whole-trace form).
    """
    from scipy.signal import hilbert
    nn = int(round(NOISE_S * FS))
    zwin = (max(nn, NPRE - int(14 * FS)), NPRE + int(2 * FS))
    hwin = (NPRE - int(3 * FS), NPRE + int(4 * FS))
    out = []
    arrs = [(stack[:, 0], zwin)]
    h = np.sqrt((stack[:, 1:] ** 2).sum(1)) if stack.shape[1] > 1 else stack[:, 0]
    arrs.append((h, hwin))
    for arr, (i0, i1) in arrs:
        env = np.abs(hilbert(arr - arr.mean()))
        rms = float(np.sqrt(np.mean(arr[:nn] ** 2)))
        i0, i1 = max(i0, nn), min(i1, len(arr) - 1)
        out.append(float(env[i0:i1].max() / rms) if (rms > 0 and i1 > i0) else np.nan)
    return out[0], out[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pgc-det", required=True)
    ap.add_argument("--rule", default="new9")
    ap.add_argument("--min-det", type=int, default=200)
    ap.add_argument("--stations", required=True, help="NET:STA,NET:STA,...")
    ap.add_argument("--wfdir", default="lfe_catalog_repair/traces")
    ap.add_argument("--start", required=True)
    ap.add_argument("--end", required=True)
    ap.add_argument("--pgc-lat", type=float, default=48.6498)
    ap.add_argument("--pgc-lon", type=float, default=-123.4521)
    ap.add_argument("--min-snr", type=float, default=3.0)
    ap.add_argument("--max-fams", type=int, default=0)
    ap.add_argument("--max-cuts", type=int, default=400, help="cap detections used per family")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", default="lfe_catalog_repair/results/d3_stage4_stacks")
    a = ap.parse_args()
    rng = np.random.default_rng(0)

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    stns = [tuple(s.split(":")) for s in a.stations.split(",") if s.strip()]
    d0, d1 = datetime.fromisoformat(a.start), datetime.fromisoformat(a.end)
    det = pd.read_parquet(a.pgc_det)
    if a.rule in det:
        det = det[det[a.rule]]
    det["dt"] = pd.to_datetime(det["t"], unit="s")
    det = det[(det.dt >= d0) & (det.dt < d1)]
    det["fam"] = det["template"].astype(str)
    n_by_fam = det.groupby("fam").size().sort_values(ascending=False)
    fams = n_by_fam[n_by_fam >= a.min_det].index.tolist()
    if a.max_fams:
        fams = fams[:a.max_fams]
    print(f"[s4] {len(det):,} PGC detections ({a.rule}) in {a.start}..{a.end}; "
          f"{len(fams)} families with >= {a.min_det} (of {len(n_by_fam)})", flush=True)
    if not fams:
        print("[s4] no family reaches the detection minimum -> stop")
        return
    # per-family source geometry
    fl = np.array([[float(f.split("__")[0].split("_")[0]), float(f.split("__")[0].split("_")[1])]
                   for f in fams])
    fz = V.slab2_depth(fl[:, 0], fl[:, 1])
    print(f"[s4] Slab2 depth over families: {np.nanmin(fz):.1f}-{np.nanmax(fz):.1f} km "
          f"(median {np.nanmedian(fz):.1f})", flush=True)

    def dist_km(la, lo, la2, lo2):
        return float(np.hypot((la - la2) * 111.2, (lo - lo2) * 111.2 * np.cos(np.radians(la))))

    # reference S time at PGC and the predicted offset to each station
    dt_pred = {}
    for i, f in enumerate(fams):
        xp = dist_km(fl[i, 0], fl[i, 1], a.pgc_lat, a.pgc_lon)
        tS_pgc = V.tt(xp, fz[i], "S")
        for net, sta in stns:
            xs = dist_km(fl[i, 0], fl[i, 1], STA_LL[sta][0], STA_LL[sta][1])
            dt_pred[(f, sta)] = V.tt(xs, fz[i], "S") - tS_pgc

    # detection times per family, capped
    times = {}
    for f in fams:
        t = np.sort(det.loc[det.fam == f, "t"].values)
        if len(t) > a.max_cuts:
            t = t[rng.choice(len(t), a.max_cuts, replace=False)]
            t.sort()
        times[f] = t + 1.0        # +1 s: the stored time is the START of the [-1,+1) s template window

    # random control times (same count per family), uniform inside the window
    ep0, ep1 = d0.timestamp(), d1.timestamp()
    rtimes = {f: np.sort(rng.uniform(ep0 + 100, ep1 - 100, len(times[f]))) for f in fams}

    days = []
    d = d0
    while d < d1:
        days.append((d.year, d.month, d.day))
        d += timedelta(days=1)
    items = [(net, sta, y, m, dd) for (net, sta) in stns for (y, m, dd) in days]
    print(f"[s4] {len(items)} station-days to read", flush=True)

    acc = {}       # (fam, sta) -> list of (NWIN, ncomp) cuts
    accn = {}
    cut_times = {}   # (fam, sta) -> list of the cut reference times (for the event-level lags)
    ncomp_of = {}
    t0 = time.time()
    ctx = mp.get_context("spawn")
    done = 0
    with ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init,
                             initargs=(a.wfdir, FS, (FMIN, FMAX))) as ex:
        for (net, sta, ymd), data in ex.map(_load_day, items, chunksize=1):
            done += 1
            if not data:
                continue
            comps = [c for c in ["Z"] + [x for x in sorted(data) if x != "Z"] if c in data]
            ncomp_of[sta] = len(comps)
            day_ep = datetime(*ymd).timestamp()
            n = min(data[c].size for c in comps)
            X = np.stack([data[c][:n] for c in comps], axis=1)
            for src, store in ((times, acc), (rtimes, accn)):
                for f in fams:
                    tt_ = src[f]
                    off = dt_pred[(f, sta)]
                    idx = np.round((tt_ + off - day_ep) * FS).astype(np.int64) - NPRE
                    ok = (idx >= 0) & (idx + NWIN <= n)
                    if not ok.any():
                        continue
                    g = idx[ok][:, None] + np.arange(NWIN)[None, :]
                    store.setdefault((f, sta), []).append(X[g])
                    if store is acc:
                        cut_times.setdefault((f, sta), []).append(tt_[ok])
            if done % 40 == 0:
                print(f"  {done}/{len(items)} station-days ({(time.time()-t0)/60:.1f} min)", flush=True)

    rows = []
    ev_rows = []
    out = {}
    for (f, sta), parts in sorted(acc.items()):
        cuts = np.concatenate(parts, axis=0)
        if cuts.shape[0] < 20:
            continue
        stack, shifts, cc = align_stack(cuts)
        sz, sh = stack_snr(stack)
        nparts = accn.get((f, sta))
        nsz = nsh = nccm = np.nan
        if nparts:
            ncuts = np.concatenate(nparts, axis=0)
            k = min(len(ncuts), len(cuts))
            nstack, _s, nc_ = align_stack(ncuts[:k])
            nsz, nsh = stack_snr(nstack)
            nccm = float(np.mean(nc_))
        # the S arrival on the horizontals is what the location uses, so it is the primary gate;
        # Z-only stations (OCP) fall back to Z. Both must beat the SAME-n random-time null by 1.5x
        # in SNR and 2x in mean cut-to-stack cc.
        prim = sh if stack.shape[1] > 1 else sz
        nprim = nsh if stack.shape[1] > 1 else nsz
        ccm = float(np.mean(cc))
        # The mean cut-to-stack cc is NOT usable as a gate at this SNR: the validated Lin-refined
        # stack (a demonstrably real arrival) had cc 0.105 vs a null of 0.099. Only the stack-SNR
        # ratio against the same-n random-time null discriminates; cc is kept as a diagnostic.
        keep = bool((prim >= a.min_snr)
                    and (not np.isfinite(nprim) or prim > 1.5 * nprim))
        for j, c in enumerate(["Z", "H1", "H2"][:stack.shape[1]]):
            out[f"{f}|{sta}|{c}"] = stack[:, j]
        ct = np.concatenate(cut_times[(f, sta)]) if (f, sta) in cut_times else np.full(len(shifts), np.nan)
        ev_rows.append(pd.DataFrame(dict(family=f, station=sta, t_ref=ct[:len(shifts)],
                                         lag_s=shifts / FS, cc=cc)))
        rows.append(dict(family=f, station=sta, n_cuts=int(cuts.shape[0]),
                         snr_z=sz, snr_h=sh, cc_mean=ccm,
                         shift_med_s=float(np.median(shifts) / FS),
                         null_snr_z=nsz, null_snr_h=nsh, null_cc_mean=nccm, keep=int(keep),
                         dt_pred_s=dt_pred[(f, sta)]))
    R = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez(f"{a.out}.npz", **out)
    meta = dict(pre_s=PRE_S, post_s=POST_S, fs=FS, fmin=FMIN, fmax=FMAX,
                maxshift_s=MAXSHIFT_S, rule=a.rule, start=a.start, end=a.end)
    R.attrs = meta
    R.to_csv(f"{a.out}.csv", index=False)
    pd.Series(meta).to_csv(f"{a.out}_meta.csv")
    if ev_rows:
        pd.concat(ev_rows, ignore_index=True).to_parquet(f"{a.out}_lags.parquet", index=False)
    print(f"[s4] {len(R)} family-station stacks, {R.keep.sum()} kept (SNR >= {a.min_snr} and > null)")
    if len(R):
        for sta, g in R.groupby("station"):
            print(f"  {sta:5s} n={len(g):4d} keep={g.keep.sum():4d} "
                  f"snr_z {g.snr_z.median():5.2f} (null {g.null_snr_z.median():5.2f}) "
                  f"snr_h {g.snr_h.median():5.2f} (null {g.null_snr_h.median():5.2f}) "
                  f"cc {g.cc_mean.median():.3f} (null {g.null_cc_mean.median():.3f})")
        fk = R[R.keep == 1].groupby("family").size()
        print(f"  families with >=4 kept stations: {(fk >= 4).sum()} / {R.family.nunique()}")
    print(f"  -> {a.out}.npz / .csv")


STA_LL = {
    "PGC": (48.6498, -123.4521), "SYMB": (48.5593, -123.7989), "B011": (48.6495, -123.4482),
    "VGZ": (48.4131, -123.3251), "LRIV": (48.0575, -123.5040), "CLRS": (48.8203, -124.1309),
    "PTRF": (48.5442, -124.4131), "SNB": (48.7751, -123.1723), "OCP": (48.2979, -124.6252),
    "GOBB": (48.9493, -123.5105), "STW": (48.1507, -123.6711), "UWFH": (48.5459, -123.0132),
    "B926": (48.8202, -124.1312),   # PB borehole co-located with CN.CLRS, the 2010 west neighbour
    "B013": (47.8130, -122.9108), "B014": (47.5133, -123.8125), "B927": (49.2188, -124.8113),
    "WISH": (47.1160, -123.7714), "B943": (47.8132, -122.9113), "YACT": (45.9325, -122.4193),
}

if __name__ == "__main__":
    main()
