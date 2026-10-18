#!/usr/bin/env python
"""D3 stage 3e — is the LFE SIGNAL there at all, or are the discovery TEMPLATES the problem?

The single-station 3-C matched filter recovers Lin/Bostock LFEs at chance level (recall ~0.001,
identical to a +600 s shifted control). Two explanations must be separated:
  (A) a 2 s single-station template is simply not sensitive enough (sigma_cc ~0.17/sqrt(3)), or
  (B) the project's discovery-stage "families" are not LFE waveforms, so no detector could work.

This script settles it WITHOUT any detector: it stacks PGC's 3-C waveforms directly at the
Lin/Bostock event times of each Lin family, at the S time PREDICTED by the 1-D model from that
family's own catalogue location and depth. If a coherent S arrival appears above a random-time null,
the LFE signal IS present at PGC and (A) is the live explanation; if it does not, the ground truth
itself is not visible at this station and neither explanation can be tested here.

It also reports, for each Lin family:
  * how far the observed stack peak sits from the model-predicted S time (a check on the velocity
    model AND on the +11 s ANCHOR the project already uses for B011),
  * the stack SNR as a function of template length (2 / 4 / 6 s), which is the input Stage 5 needs,
  * the same numbers for the NEAREST discovery template's stack, side by side.

Output: results/d3_lin_template_test.csv  +  figures/d3_lin_template_test.png
"""
from __future__ import annotations

import argparse
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import d3_vel1d as V                                        # noqa: E402
from d3_stage4_stacks import _init, _load_day, align_stack, FS   # noqa: E402

LIN_DET = "data/detections_lin_ets_2010_vi.parquet"
LIN_FAM = "catalogs/lin_families_ets_2010_vi.csv"
PGC = (48.6498, -123.4521)
# ORIGIN-TIME referenced, with a WIDE search. The first version of this test stacked at the S time
# PREDICTED by the 1-D model from the Lin family depth (3-14 s) and searched +-4 s around it, and
# found nothing -- but the project's own calibration (`refine_lin_picks.py`, ANCHOR = 11.0 s) puts
# the S arrival at B011 (0.3 km from PGC) at OT + ~11 s, i.e. up to 8 s outside that search window.
# The Lin family DEPTHS in `lin_families_ets_2010_vi.csv` are also suspect (12-35 km, where the
# interface is 30-45 km) because `lin_catalog.py` sign-flips them and averages over a 0.1 deg cell.
# So the window is now referenced to the ORIGIN TIME and the arrival is searched over the whole
# plausible travel-time range, which MEASURES the offset instead of assuming it.
PRE, POST = 12.0, 30.0
NPRE = int(round(PRE * FS))          # sample of the ORIGIN time
NWIN = int(round((PRE + POST) * FS))
NOISE = int(round(8.0 * FS))         # OT-12 s .. OT-4 s
SEARCH_LO, SEARCH_HI = 0.0, 29.0     # s after OT


def snr_of(arr, i0, i1):
    from scipy.signal import hilbert
    env = np.abs(hilbert(arr - arr.mean()))
    rms = float(np.sqrt(np.mean(arr[:NOISE] ** 2)))
    i0, i1 = max(i0, NOISE), min(i1, len(arr) - 1)
    return (float(env[i0:i1].max() / rms) if (rms > 0 and i1 > i0) else np.nan,
            int(i0 + np.argmax(env[i0:i1])) if i1 > i0 else -1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="2010-08-10")
    ap.add_argument("--end", default="2010-09-13")
    ap.add_argument("--wfdir", default="lfe_catalog_repair/traces")
    ap.add_argument("--net", default="CN")
    ap.add_argument("--sta", default="PGC")
    ap.add_argument("--top", type=int, default=25, help="Lin families by event count in the window")
    ap.add_argument("--max-cuts", type=int, default=800)
    ap.add_argument("--refine", action="store_true",
                    help="also report the cc-refined stack (default reports BOTH)")
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--out", default="lfe_catalog_repair/results/d3_lin_template_test.csv")
    a = ap.parse_args()
    rng = np.random.default_rng(0)

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor

    lin = pd.read_parquet(LIN_DET)
    lin["time"] = pd.to_datetime(lin["time"])
    lin["ep"] = lin["time"].astype("datetime64[ns]").astype("int64") / 1e9
    lf = pd.read_csv(LIN_FAM).set_index("family_id")
    d0, d1 = datetime.fromisoformat(a.start), datetime.fromisoformat(a.end)
    sel = lin[(lin.time >= d0) & (lin.time < d1)]
    cnt = sel.groupby("family_id").size().sort_values(ascending=False)
    fams = cnt.index[:a.top].tolist()
    print(f"[lintpl] {len(sel):,} Lin events in {a.start}..{a.end}; "
          f"top {len(fams)} families hold {int(cnt[:a.top].sum()):,}", flush=True)

    geom = {}
    for f in fams:
        la, lo, dep = float(lf.lat[f]), float(lf.lon[f]), float(lf.depth_km[f])
        x = float(np.hypot((la - PGC[0]) * 111.2,
                           (lo - PGC[1]) * 111.2 * np.cos(np.radians(la))))
        geom[f] = (la, lo, dep, x, V.tt(x, dep, "S"), V.tt(x, dep, "P"))
    print("[lintpl] family geometry (dist km, Lin depth km, predicted tS, tS-tP):")
    for f in fams[:8]:
        la, lo, dep, x, tS, tP = geom[f]
        print(f"   {f}  x={x:5.1f} z={dep:5.1f}  tS={tS:5.2f} s  S-P={tS-tP:4.2f} s  n={cnt[f]}")

    times = {f: np.sort(sel.loc[sel.family_id == f, "ep"].values) for f in fams}
    for f in fams:
        if len(times[f]) > a.max_cuts:
            t = times[f][rng.choice(len(times[f]), a.max_cuts, replace=False)]
            times[f] = np.sort(t)
    ep0, ep1 = d0.timestamp(), d1.timestamp()
    rtimes = {f: np.sort(rng.uniform(ep0 + 100, ep1 - 100, len(times[f]))) for f in fams}

    days = []
    d = d0
    while d < d1:
        days.append((a.net, a.sta, d.year, d.month, d.day))
        d += timedelta(days=1)
    acc, accn = {}, {}
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx, initializer=_init,
                             initargs=(a.wfdir, FS, (1.0, 8.0))) as ex:
        for (net, sta, ymd), data in ex.map(_load_day, days, chunksize=1):
            if not data:
                continue
            comps = [c for c in ["Z"] + [x for x in sorted(data) if x != "Z"] if c in data]
            n = min(data[c].size for c in comps)
            X = np.stack([data[c][:n] for c in comps], axis=1)
            day_ep = datetime(*ymd).timestamp()
            for src, store in ((times, acc), (rtimes, accn)):
                for f in fams:
                    idx = np.round((src[f] - day_ep) * FS).astype(np.int64) - NPRE
                    ok = (idx >= 0) & (idx + NWIN <= n)
                    if not ok.any():
                        continue
                    g = idx[ok][:, None] + np.arange(NWIN)[None, :]
                    store.setdefault(f, []).append(X[g])

    rows = []
    stacks = {}
    for f in fams:
        if f not in acc:
            continue
        cuts = np.concatenate(acc[f], axis=0)
        ncuts = np.concatenate(accn[f], axis=0) if f in accn else None
        r = dict(family=f, n_events=int(cuts.shape[0]), dist_km=geom[f][3],
                 lin_depth=geom[f][2], tS_pred=geom[f][4])
        for tag, refine in (("raw", 0), ("refined", 2)):
            if refine:
                st, sh, cc = align_stack(cuts, niter=refine)
                r[f"cc_mean_{tag}"] = float(np.mean(cc))
            else:
                # NO alignment: a straight stack at the predicted time cannot manufacture coherence
                amp = np.sqrt((cuts ** 2).mean(axis=(1, 2)))
                amp[amp <= 0] = 1.0
                st = (cuts / amp[:, None, None]).mean(0)
            h = np.sqrt((st[:, 1:] ** 2).sum(1)) if st.shape[1] > 1 else st[:, 0]
            lo = NPRE + int(SEARCH_LO * FS)
            hi = NPRE + int(SEARCH_HI * FS)
            s_h, i_h = snr_of(h, lo, hi)
            s_z, i_z = snr_of(st[:, 0], lo, hi)
            r[f"snr_h_{tag}"] = s_h
            r[f"snr_z_{tag}"] = s_z
            r[f"peak_after_OT_s_{tag}"] = (i_h - NPRE) / FS
            r[f"peakZ_after_OT_s_{tag}"] = (i_z - NPRE) / FS
            if ncuts is not None:
                k = min(len(ncuts), len(cuts))
                if refine:
                    nst, _s, ncc = align_stack(ncuts[:k], niter=refine)
                    r[f"null_cc_mean_{tag}"] = float(np.mean(ncc))
                else:
                    amp = np.sqrt((ncuts[:k] ** 2).mean(axis=(1, 2)))
                    amp[amp <= 0] = 1.0
                    nst = (ncuts[:k] / amp[:, None, None]).mean(0)
                nh = np.sqrt((nst[:, 1:] ** 2).sum(1)) if nst.shape[1] > 1 else nst[:, 0]
                r[f"null_snr_h_{tag}"] = snr_of(nh, NPRE + int(SEARCH_LO * FS),
                                                NPRE + int(SEARCH_HI * FS))[0]
                r[f"null_snr_z_{tag}"] = snr_of(nst[:, 0], NPRE + int(SEARCH_LO * FS),
                                                NPRE + int(SEARCH_HI * FS))[0]
            if tag == "raw":
                stacks[f] = st
        rows.append(r)
    T = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    T.to_csv(a.out, index=False)
    np.savez(str(a.out).replace(".csv", "_stacks.npz"), **{f: v for f, v in stacks.items()})
    print("\n family        n_ev  dist  RAW snr_h (null)  snr_z (null)  peakH-OT  peakZ-OT   tS_pred")
    for _, r in T.iterrows():
        print(f" {r.family:8s} {int(r.n_events):6d} {r.dist_km:5.1f}  "
              f"{r.snr_h_raw:6.2f} ({r.null_snr_h_raw:5.2f})  "
              f"{r.snr_z_raw:6.2f} ({r.null_snr_z_raw:5.2f})   "
              f"{r.peak_after_OT_s_raw:+7.2f}  {r.peakZ_after_OT_s_raw:+7.2f}  {r.tS_pred:6.2f}")
    print(f"\n RAW stack (no alignment, cannot manufacture coherence):")
    print(f"   snr_h median {T.snr_h_raw.median():.2f} vs null {T.null_snr_h_raw.median():.2f}; "
          f"families with snr_h > 2x null: {int((T.snr_h_raw > 2*T.null_snr_h_raw).sum())}/{len(T)}")
    print(f"   snr_z median {T.snr_z_raw.median():.2f} vs null {T.null_snr_z_raw.median():.2f}")
    print(f"   horizontal peak AFTER the Lin origin time: median "
          f"{T.peak_after_OT_s_raw.median():+.2f} s [p10 {T.peak_after_OT_s_raw.quantile(.1):+.2f}, "
          f"p90 {T.peak_after_OT_s_raw.quantile(.9):+.2f}]  "
          f"(project ANCHOR = +11.0 s; 1-D model from the Lin depths gives "
          f"{T.tS_pred.median():.1f} s)")
    print(f" REFINED: snr_h median {T.snr_h_refined.median():.2f} vs null "
          f"{T.null_snr_h_refined.median():.2f}; cc {T.cc_mean_refined.median():.3f} vs null "
          f"{T.null_cc_mean_refined.median():.3f}")
    print(f"  -> {a.out}")


if __name__ == "__main__":
    main()
