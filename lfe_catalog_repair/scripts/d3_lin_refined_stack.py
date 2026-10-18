#!/usr/bin/env python
"""D3 stage 3f — the pipeline audit that decides the whole D3 verdict.

The Lin ORIGIN times in `data/detections_lin_ets_2010_vi.parquet` sit **11.28 s before the S
arrival with a scatter of +-2.2 s** (measured here from the project's own refined picks in
`lfe_features/data/lin_precise_b011.csv`: median arrival-OT = 11.28 s, p05 8.0, p95 14.0,
std 2.22 s). At 1-8 Hz a +-2.2 s label error destroys stack coherence outright, so "stack at
OT + 11 s" CANNOT work and its null result says nothing about the data.

This script therefore stacks PGC at the project's **refined** arrival times (776 curated Lin events
inside 2010-08-10..2010-09-13, aligned by cross-correlation at B011, which is 0.3 km from PGC).
  * If a coherent S arrival appears above a random-time null -> the LFE signal IS present at PGC,
    the D3 processing chain is sound, and the failure of the single-station matched filter is a
    SENSITIVITY problem (the Stage 5 hypothesis).
  * If it does not -> the ground truth is not visible at this station and the pilot cannot use it.
Reported for template lengths 2 / 4 / 6 s as well, which is the input Stage 5 needs.

Output: results/d3_lin_refined_stack.txt + _stack.npz
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from d3_stage4_stacks import _init, _load_day, align_stack, FS   # noqa: E402

PRE, POST = 12.0, 12.0
NPRE = int(round(PRE * FS))
NWIN = int(round((PRE + POST) * FS))
NOISE = int(round(8.0 * FS))


def snr_of(arr, i0, i1):
    from scipy.signal import hilbert
    env = np.abs(hilbert(arr - arr.mean()))
    rms = float(np.sqrt(np.mean(arr[:NOISE] ** 2)))
    i0, i1 = max(i0, NOISE), min(i1, len(arr) - 1)
    if not (rms > 0 and i1 > i0):
        return np.nan, -1
    return float(env[i0:i1].max() / rms), int(i0 + np.argmax(env[i0:i1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--times", default="lfe_catalog_repair/results/d3_lin_refined_window.csv")
    ap.add_argument("--start", default="2010-08-10")
    ap.add_argument("--end", default="2010-09-13")
    ap.add_argument("--net", default="CN")
    ap.add_argument("--sta", default="PGC")
    ap.add_argument("--wfdir", default="lfe_catalog_repair/traces")
    ap.add_argument("--band", default="1,8")
    ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--out", default="lfe_catalog_repair/results/d3_lin_refined_stack")
    a = ap.parse_args()
    rng = np.random.default_rng(0)
    fmin, fmax = [float(x) for x in a.band.split(",")]

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    T = pd.read_csv(a.times)
    t_arr = np.sort(T["arrival"].values.astype(float))
    d0, d1 = datetime.fromisoformat(a.start), datetime.fromisoformat(a.end)
    ep0, ep1 = d0.timestamp(), d1.timestamp()
    t_arr = t_arr[(t_arr > ep0 + 60) & (t_arr < ep1 - 60)]
    t_rnd = np.sort(rng.uniform(ep0 + 100, ep1 - 100, len(t_arr)))
    lines = [f"[refined] {len(t_arr)} refined Lin arrivals at {a.net}.{a.sta} "
             f"{a.start}..{a.end}, band {fmin}-{fmax} Hz"]
    print(lines[-1], flush=True)

    days = []
    d = d0
    while d < d1:
        days.append((a.net, a.sta, d.year, d.month, d.day))
        d += timedelta(days=1)
    acc, accn = [], []
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init,
                             initargs=(a.wfdir, FS, (fmin, fmax))) as ex:
        for (net, sta, ymd), data in ex.map(_load_day, days, chunksize=1):
            if not data:
                continue
            comps = [c for c in ["Z"] + [x for x in sorted(data) if x != "Z"] if c in data]
            n = min(data[c].size for c in comps)
            X = np.stack([data[c][:n] for c in comps], axis=1)
            day_ep = datetime(*ymd).timestamp()
            for src, store in ((t_arr, acc), (t_rnd, accn)):
                idx = np.round((src - day_ep) * FS).astype(np.int64) - NPRE
                ok = (idx >= 0) & (idx + NWIN <= n)
                if ok.any():
                    g = idx[ok][:, None] + np.arange(NWIN)[None, :]
                    store.append(X[g])
    if not acc:
        print("no cuts")
        return
    cuts = np.concatenate(acc, axis=0)
    ncuts = np.concatenate(accn, axis=0)
    lines.append(f"[refined] cut {cuts.shape[0]} real / {ncuts.shape[0]} random windows, "
                 f"{cuts.shape[2]} components")
    print(lines[-1], flush=True)

    def stack_raw(c):
        amp = np.sqrt((c ** 2).mean(axis=(1, 2)))
        amp[amp <= 0] = 1.0
        return (c / amp[:, None, None]).mean(0)

    for tag, fn in (("RAW (no alignment)", stack_raw),
                    ("cc-REFINED (2 iterations)", lambda c: align_stack(c, niter=2)[0])):
        st = fn(cuts)
        nst = fn(ncuts[:len(cuts)])
        h = np.sqrt((st[:, 1:] ** 2).sum(1)) if st.shape[1] > 1 else st[:, 0]
        nh = np.sqrt((nst[:, 1:] ** 2).sum(1)) if nst.shape[1] > 1 else nst[:, 0]
        sh, ih = snr_of(h, NPRE - int(2 * FS), NPRE + int(4 * FS))
        nsh, _ = snr_of(nh, NPRE - int(2 * FS), NPRE + int(4 * FS))
        sz, iz = snr_of(st[:, 0], NPRE - int(2 * FS), NPRE + int(4 * FS))
        nsz, _ = snr_of(nst[:, 0], NPRE - int(2 * FS), NPRE + int(4 * FS))
        lines.append(f"  {tag}: snr_h {sh:.2f} (null {nsh:.2f}, ratio {sh/max(nsh,1e-9):.2f}x)  "
                     f"snr_z {sz:.2f} (null {nsz:.2f}, ratio {sz/max(nsz,1e-9):.2f}x)  "
                     f"peak_h at {(ih-NPRE)/FS:+.2f} s from the refined arrival")
        print(lines[-1], flush=True)
        if tag.startswith("RAW"):
            np.savez(f"{a.out}.npz", stack=st, null=nst)
            # template-length sweep on the RAW stack: cc-noise sigma expected ~ 1/sqrt(2*B*L)
            for L in (2.0, 4.0, 6.0):
                nL = int(round(L * FS))
                i0 = NPRE - int(round(0.5 * FS))
                seg = h[i0:i0 + nL]
                nseg = nh[i0:i0 + nL]
                e = float(np.sqrt((seg ** 2).mean()))
                ne = float(np.sqrt((nseg ** 2).mean()))
                sig = 1.0 / np.sqrt(2 * (fmax - fmin) * L)
                lines.append(f"    template L={L:.0f} s: rms(signal window)/rms(null window) "
                             f"= {e/max(ne,1e-12):.2f}; expected single-channel cc sigma "
                             f"~ 1/sqrt(2*B*L) = {sig:.3f}")
                print(lines[-1], flush=True)
    Path(f"{a.out}.txt").write_text("\n".join(lines) + "\n")
    print(f"  -> {a.out}.txt / .npz")


if __name__ == "__main__":
    main()
