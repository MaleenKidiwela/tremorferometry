#!/usr/bin/env python
"""GPU-accelerated NLLB-style family discovery -- stage 2 (clustering).

Three speedups over discover_nllb_pnsn_driven.py, all validated equivalent:
  1. LOAD EACH DAY ONCE.  The CPU script re-loads a day file for every spatial bin
     that has a candidate on it (HDW: ~100k day-loads).  Here we group all candidate
     cuts by calendar day and load each day exactly once (~11k loads).
  2. resample_poly instead of obspy's FFT resample (0.12 s vs 1.1 s per 100 Hz day,
     and no spurious dv/v drift -- same lesson as the densify pipeline).
  3. GPU all-pairs max-shifted CC (repeater_gpu): 11.6 s -> 0.06 s for a 2000-window
     bin, max|d|=1.8e-7 vs CPU, ZERO cluster-decision flips at cc>=0.8.

Single GPU process; CPU workers only load+cut waveforms (KB each), never the N^2 CC
matrix -- so memory stays bounded and we run at the FULL 2000-cap (no OOM, no 50-cap
undersample).

Output schema is IDENTICAL to discover_nllb_pnsn_driven.py AND adds the `snr` column
(template SNR = env peak / pre-pulse RMS) that select_coverage_families.py /
plot_gnw_family_map.py require -- the step that was silently missing.
  <out>.npz                : per-family template, keyed by family_id
  <out>.summary.csv        : family_id,n_members,n_years,first_year,last_year,lat,lon,snr
  <out>.members.parquet    : member detection times
"""
import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--station", default="HDW")
    p.add_argument("--wfdir", default="data/waveforms")
    p.add_argument("--candidates", default="data/hdw_pnsn_candidates.parquet")
    p.add_argument("--out", default="data/hdw_pnsn_families.npz")
    p.add_argument("--grid-deg", type=float, default=0.05)
    p.add_argument("--fs", type=float, default=40.0)
    p.add_argument("--fmin", type=float, default=2.0)
    p.add_argument("--fmax", type=float, default=8.0)
    p.add_argument("--template-pre", type=float, default=-1.0)
    p.add_argument("--template-post", type=float, default=1.0)
    p.add_argument("--max-shift-samples", type=int, default=20)
    p.add_argument("--cc-threshold", type=float, default=0.80)
    p.add_argument("--min-bin-candidates", type=int, default=10)
    p.add_argument("--max-bin-candidates", type=int, default=2000)
    p.add_argument("--min-family-members", type=int, default=3)
    p.add_argument("--min-years", type=int, default=3)
    p.add_argument("--workers", type=int, default=24)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--limit-bins", type=int, default=0,
                   help="debug: only process the first N bins (0 = all)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# CPU I/O worker -- load ONE day with resample_poly, cut+normalize its candidates.
# Returns list of (cand_idx, cut_vector). Never imports cupy.
# ---------------------------------------------------------------------------
def _cut_day(args_tuple):
    from obspy import read
    from scipy.signal import resample_poly

    (year, jday), items, cfg = args_tuple
    cands = list(Path(cfg["wfdir"]).glob(
        f"*.{cfg['station']}/{year}/{jday:03d}.mseed"))
    if not cands:
        return []
    try:
        st = read(str(cands[0]))
    except Exception:
        return []
    st = st.select(component="Z")
    if len(st) == 0:
        return []
    fs = cfg["fs"]
    # per-trace resample_poly to fs (also fixes FS jitter that breaks merge)
    for tr in st:
        nat = float(tr.stats.sampling_rate)
        if abs(nat - fs) > 1e-6:
            fr = Fraction(int(round(fs)), int(round(nat))).limit_denominator(1000)
            try:
                tr.data = resample_poly(
                    tr.data.astype(np.float64), fr.numerator, fr.denominator
                ).astype(np.float32)
                tr.stats.sampling_rate = fs
            except Exception:
                return []
    try:
        st.merge(fill_value=0)
    except Exception:
        return []
    if len(st) == 0:
        return []
    st.detrend("demean")
    st.filter("bandpass", freqmin=cfg["fmin"], freqmax=cfg["fmax"],
              corners=4, zerophase=True)
    tr = st[0]
    data = tr.data.astype(np.float32)
    start_ts = tr.stats.starttime.timestamp

    template_n = cfg["template_n"]
    pre_off = int(round(cfg["template_pre"] * fs))
    out = []
    for cand_idx, t_us in items:
        offset_s = (t_us / 1e6) - start_ts
        i_peak = int(round(offset_s * fs))
        i0 = i_peak + pre_off
        i1 = i0 + template_n
        if i0 < 0 or i1 > data.size:
            continue
        seg = data[i0:i1].astype(np.float64)
        seg = seg - seg.mean()
        nrm = float(np.linalg.norm(seg))
        if nrm == 0 or not np.isfinite(nrm):
            continue
        out.append((int(cand_idx), (seg / nrm).astype(np.float32)))
    return out


def _template_snr(T):
    """env peak / pre-pulse RMS -- reproduces GNW family SNR definition."""
    from scipy.signal import hilbert
    env = np.abs(hilbert(T))
    pre = T[:34]                      # ~1 s before the centered peak, guard=5 samples
    rms = float(np.sqrt(np.mean(pre ** 2)))
    return float(env.max() / rms) if rms > 0 else 0.0


def main():
    args = parse_args()
    t_start = time.time()
    template_n = int(round((args.template_post - args.template_pre) * args.fs))

    # ----- candidates + spatial binning (identical to discover_nllb_pnsn_driven) -----
    cand = pd.read_parquet(args.candidates)
    cand["time"] = pd.to_datetime(cand["time"])
    cand["lat_bin"] = np.round(cand["lat"] / args.grid_deg) * args.grid_deg
    cand["lon_bin"] = np.round(cand["lon"] / args.grid_deg) * args.grid_deg
    cand["bin_label"] = (cand["lat_bin"].map(lambda v: f"{v:.3f}") + "_"
                         + cand["lon_bin"].map(lambda v: f"{v:.3f}"))
    cand["time_us"] = cand["time"].astype("datetime64[ns]").astype("int64") // 1000

    bc = cand["bin_label"].value_counts()
    keep_bins = bc[bc >= args.min_bin_candidates].index.tolist()
    if args.limit_bins:
        keep_bins = keep_bins[:args.limit_bins]
    cand = cand[cand["bin_label"].isin(keep_bins)].copy()
    print(f"[1/3] {len(cand):,} candidates in {len(keep_bins)} bins "
          f">={args.min_bin_candidates} cands", flush=True)

    # cap per bin (subsample) -- same semantics as the CPU script's max_bin_candidates
    rng = np.random.default_rng(args.seed)
    parts = []
    for _, sub in cand.groupby("bin_label", sort=False):
        if len(sub) > args.max_bin_candidates:
            sub = sub.iloc[rng.choice(len(sub), args.max_bin_candidates, replace=False)]
        parts.append(sub)
    cand = pd.concat(parts).reset_index(drop=True)
    cand["cand_idx"] = np.arange(len(cand), dtype=np.int64)
    cand["year"] = cand["time"].dt.year
    cand["jday"] = cand["time"].dt.dayofyear

    # ----- PASS 1: load each day ONCE, cut+normalize all its candidate windows -----
    cfg = dict(wfdir=args.wfdir, station=args.station, fs=args.fs,
               fmin=args.fmin, fmax=args.fmax,
               template_n=template_n, template_pre=args.template_pre)
    day_tasks = []
    for (y, j), g in cand.groupby(["year", "jday"], sort=False):
        items = list(zip(g["cand_idx"].values.tolist(), g["time_us"].values.tolist()))
        day_tasks.append(((int(y), int(j)), items, cfg))
    print(f"[2/3] cutting windows: {len(day_tasks)} unique days, "
          f"{args.workers} workers (load-each-day-once + resample_poly)", flush=True)

    cuts = {}
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = [ex.submit(_cut_day, t) for t in day_tasks]
        for f in as_completed(futs):
            for ci, vec in f.result():
                cuts[ci] = vec
            done += 1
            if done % 500 == 0:
                print(f"  day {done}/{len(day_tasks)}, cuts {len(cuts):,} "
                      f"({time.time()-t_start:.0f}s)", flush=True)
    print(f"  cut {len(cuts):,}/{len(cand):,} candidates "
          f"({time.time()-t_start:.0f}s)", flush=True)

    # ----- PASS 2: per-bin GPU CC + complete-linkage clustering -----
    import cupy as cp  # noqa: F401  (GPU only here, never in workers)
    from tremorferometry.repeater_gpu import all_pairs_cc_max_shifted_gpu
    from scipy.cluster.hierarchy import linkage, fcluster

    print(f"[3/3] GPU clustering {len(keep_bins)} bins "
          f"(cc>={args.cc_threshold} complete-linkage)", flush=True)
    families = []
    nb = 0
    for bl, g in cand.groupby("bin_label", sort=False):
        rows = [(cuts[ci], pd.Timestamp(t).to_pydatetime())
                for ci, t in zip(g["cand_idx"].values, g["time"].values)
                if ci in cuts]
        nb += 1
        if len(rows) < args.min_family_members:
            continue
        X = np.array([r[0] for r in rows])
        times = [r[1] for r in rows]
        n = X.shape[0]
        cc = all_pairs_cc_max_shifted_gpu(X, args.max_shift_samples)
        iu = np.triu_indices(n, k=1)
        dist = 1.0 - np.clip(cc[iu], -1.0, 1.0)
        Z = linkage(dist, method="complete")
        labels = fcluster(Z, t=(1.0 - args.cc_threshold), criterion="distance")
        lat_c = float(g["lat"].mean())
        lon_c = float(g["lon"].mean())
        for c in np.unique(labels):
            mi = np.flatnonzero(labels == c)
            if mi.size < args.min_family_members:
                continue
            yrs = sorted({times[i].year for i in mi})
            if len(yrs) < args.min_years:
                continue
            T = X[mi].mean(axis=0)
            T = (T / (np.linalg.norm(T) + 1e-12)).astype(np.float32)
            families.append(dict(
                family_id=f"{bl}__c{int(c)}", n_members=int(mi.size),
                n_years=len(yrs), first_year=int(min(yrs)), last_year=int(max(yrs)),
                lat=lat_c, lon=lon_c, snr=_template_snr(T), template=T,
                member_times=[times[i] for i in mi],
            ))
        if nb % 50 == 0:
            print(f"  bin {nb}/{len(keep_bins)}, families {len(families)} "
                  f"({time.time()-t_start:.0f}s)", flush=True)

    print(f"\n[done] {len(families)} families in {time.time()-t_start:.0f}s", flush=True)

    # ----- save (identical schema + snr) -----
    npz_arrays, summary_rows, member_rows = {}, [], []
    for fam in families:
        npz_arrays[fam["family_id"]] = fam["template"]
        summary_rows.append({k: fam[k] for k in
                             ("family_id", "n_members", "n_years", "first_year",
                              "last_year", "lat", "lon", "snr")})
        for t in fam["member_times"]:
            member_rows.append(dict(family_id=fam["family_id"], time=t))
    np.savez(args.out, **npz_arrays)
    pd.DataFrame(summary_rows).to_csv(
        Path(args.out).with_suffix(".summary.csv"), index=False)
    pd.DataFrame(member_rows).to_parquet(
        Path(args.out).with_suffix(".members.parquet"), index=False)
    print(f"  saved -> {args.out} + .summary.csv (with snr) + .members.parquet")


if __name__ == "__main__":
    main()
