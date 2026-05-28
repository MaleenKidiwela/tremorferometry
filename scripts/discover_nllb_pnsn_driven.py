"""PNSN-driven family discovery at NLLB (Branch B of §8.2).

The exact analog of the procedure that yielded the 16 STRICT families at
PGC, but run at NLLB. Single-station, no Lin catalog needed.

Stage 1 -- candidate detection
    For every PNSN tremor window in the V.I. bbox, find envelope peaks at
    NLLB above an SNR threshold. Each peak inherits the tremor window's
    lat/lon (used later for spatial binning). Parallel over days.

Stage 2 -- per-bin clustering
    Bin candidates by 0.05 deg lat/lon. For each bin with enough
    candidates:
        - Cut tight 2-s windows centered on each peak at NLLB
          (bandpass 2-8 Hz, L2-normalize).
        - All-pairs max-shifted CC.
        - Complete-linkage cluster at CC >= 0.80 (every pair within a
          family must have CC >= 0.80).
        - Keep clusters with >=3 members spanning >=3 years.

Outputs:
    data/nllb_pnsn_families.npz : per-family NLLB template, keyed by family_id
    data/nllb_pnsn_families.summary.csv : centroid, n_members, year span
    data/nllb_pnsn_families.members.parquet : member detection times
"""

from __future__ import annotations

import argparse
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster

sys.path.insert(0, "src")
from tremorferometry.repeater import all_pairs_cc_max_shifted, _load_day_filt  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--pnsn", nargs="+",
                   default=["catalogs/pnsn_tremor_2010.csv",
                            "catalogs/pnsn_tremor_2014-2026.csv"])
    p.add_argument("--wfdir", default="data/waveforms")
    p.add_argument("--station", default="NLLB")
    p.add_argument("--bbox", nargs=4, type=float,
                   default=[48.0, 50.0, -125.5, -122.5],
                   help="lat_min lat_max lon_min lon_max")
    p.add_argument("--grid-deg", type=float, default=0.05)
    p.add_argument("--window-s", type=float, default=300.0,
                   help="PNSN tremor window length (seconds)")
    p.add_argument("--snr-threshold", type=float, default=3.0)
    p.add_argument("--peak-sep-s", type=float, default=6.0)
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
    p.add_argument("--workers", type=int, default=32)
    p.add_argument("--candidates-out", default="data/nllb_pnsn_candidates.parquet")
    p.add_argument("--use-cached-candidates", action="store_true",
                   help="Skip stage 1 if candidates parquet already exists")
    p.add_argument("--out", default="data/nllb_pnsn_families.npz")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


# Top-level worker for ProcessPool (must be pickleable).
def _candidates_one_day(args_tuple):
    """Run envelope-peak detection at one day at one station for given windows.

    args_tuple: (day_str, windows_us_pairs, lat_lons, cfg)
        windows_us_pairs: list of (t0_us, t1_us) microsecond timestamps
        lat_lons:         list of (lat, lon)
        cfg: dict of parameters
    """
    from obspy import UTCDateTime
    from scipy.signal import hilbert

    day_str, windows_us, latlons, cfg = args_tuple
    day = datetime.fromisoformat(day_str)
    loaded = _load_day_filt(
        Path(cfg["wfdir"]), cfg["station"], day,
        (cfg["fmin"], cfg["fmax"]), cfg["fs"],
    )
    if loaded is None:
        return day_str, []
    data, start = loaded
    fs = float(cfg["fs"])
    n_sep = int(cfg["peak_sep_s"] * fs)
    rows = []
    end_utc = start + (data.size - 1) / fs
    for (t0_us, t1_us), (lat, lon) in zip(windows_us, latlons):
        t0_utc = UTCDateTime(datetime.utcfromtimestamp(t0_us / 1e6))
        t1_utc = UTCDateTime(datetime.utcfromtimestamp(t1_us / 1e6))
        if t0_utc < start or t1_utc > end_utc:
            continue
        i0 = int(round((t0_utc - start) * fs))
        i1 = int(round((t1_utc - start) * fs))
        if i1 <= i0 + n_sep * 2:
            continue
        seg = data[i0:i1]
        env = np.abs(hilbert(seg))
        bg = float(np.median(env))
        if bg == 0:
            continue
        thresh = bg * cfg["snr_threshold"]
        if not (env > thresh).any():
            continue
        i = 0
        while i < env.size:
            if env[i] >= thresh:
                j_end = min(env.size, i + n_sep)
                j = int(np.argmax(env[i:j_end])) + i
                peak_time = (start + (i0 + j) / fs).datetime
                snr = float(env[j] / bg)
                rows.append(dict(time=peak_time, lat=lat, lon=lon,
                                 snr=snr, station=cfg["station"]))
                i = j + n_sep
            else:
                i += 1
    return day_str, rows


def _process_bin(args_tuple):
    """All-pairs CC + complete-linkage cluster for one spatial bin."""
    bin_label, peak_times_us, cfg = args_tuple
    from scipy.signal import hilbert

    fs = float(cfg["fs"])
    template_n = int(round((cfg["template_post"] - cfg["template_pre"]) * fs))
    half_template = template_n // 2

    # Cut tight 2-s windows centered on each peak time
    cuts = []
    times_kept = []
    times_dt = [datetime.utcfromtimestamp(t / 1e6) for t in peak_times_us]
    # Sort by day to amortize day loading
    order = np.argsort([t.timetuple().tm_yday + 1000 * t.year for t in times_dt])
    times_dt = [times_dt[i] for i in order]

    cur_day = None
    cur_loaded = None
    for dt in times_dt:
        day = datetime(dt.year, dt.month, dt.day)
        if day != cur_day:
            cur_loaded = _load_day_filt(
                Path(cfg["wfdir"]), cfg["station"], day,
                (cfg["fmin"], cfg["fmax"]), fs,
            )
            cur_day = day
        if cur_loaded is None:
            continue
        data, start = cur_loaded
        offset_s = (dt - start.datetime).total_seconds()
        i_peak = int(round(offset_s * fs))
        i0 = i_peak + int(round(cfg["template_pre"] * fs))
        i1 = i0 + template_n
        if i0 < 0 or i1 > data.size:
            continue
        seg = data[i0:i1].astype(np.float64)
        seg = seg - seg.mean()
        nrm = float(np.linalg.norm(seg))
        if nrm == 0 or not np.isfinite(nrm):
            continue
        cuts.append((seg / nrm).astype(np.float32))
        times_kept.append(dt)

    if len(cuts) < cfg["min_family_members"]:
        return bin_label, None
    X = np.array(cuts)

    cc = all_pairs_cc_max_shifted(X, max_shift_samples=cfg["max_shift_samples"])
    np.fill_diagonal(cc, 1.0)
    # Complete-linkage: every pair within a cluster must satisfy d <= 1 - cc_thresh.
    # Use scipy linkage with method='complete'.
    # Build a distance vector (condensed). Distance = 1 - cc.
    n = cc.shape[0]
    if n < 2:
        return bin_label, None
    # condensed upper triangle
    iu = np.triu_indices(n, k=1)
    # Clip CC to [-1, 1] before computing distance: FFT-based all-pairs CC
    # can produce values marginally above 1.0 due to numerical roundoff,
    # which translates to negative distances and breaks scipy.linkage.
    cc_clipped = np.clip(cc[iu], -1.0, 1.0)
    dist = 1.0 - cc_clipped
    Z = linkage(dist, method="complete")
    labels = fcluster(Z, t=(1.0 - cfg["cc_threshold"]), criterion="distance")
    # Filter clusters: at least min_family_members, spanning at least min_years
    out_clusters = []
    for c in np.unique(labels):
        idx = np.flatnonzero(labels == c)
        if idx.size < cfg["min_family_members"]:
            continue
        years = sorted({times_kept[i].year for i in idx})
        if len(years) < cfg["min_years"]:
            continue
        # Family template = mean of L2-normed members, re-L2-normed
        T = X[idx].mean(axis=0)
        T = (T / (np.linalg.norm(T) + 1e-12)).astype(np.float32)
        out_clusters.append(dict(
            family_id=f"{bin_label}__c{int(c)}",
            n_members=int(idx.size),
            first_year=int(min(years)),
            last_year=int(max(years)),
            n_years=int(len(years)),
            template=T,
            member_times=[times_kept[i] for i in idx],
        ))
    return bin_label, out_clusters


def main():
    args = parse_args()

    # Stage 1: candidate detection at NLLB
    print("[1/3] Loading PNSN catalogs and filtering to V.I. bbox...")
    cats = [pd.read_csv(p) for p in args.pnsn]
    pnsn = pd.concat(cats, ignore_index=True)
    pnsn["time"] = pd.to_datetime(pnsn["time"])
    lat_min, lat_max, lon_min, lon_max = args.bbox
    mask = ((pnsn["lat"] >= lat_min) & (pnsn["lat"] <= lat_max)
            & (pnsn["lon"] >= lon_min) & (pnsn["lon"] <= lon_max))
    pnsn = pnsn[mask].copy()
    print(f"  {len(pnsn):,} PNSN tremor windows in V.I. bbox "
          f"({pnsn.time.min()} -- {pnsn.time.max()})")

    # Skip stage 1 if cached candidates exist
    if args.use_cached_candidates and Path(args.candidates_out).exists():
        cand_df = pd.read_parquet(args.candidates_out)
        cand_df["time"] = pd.to_datetime(cand_df["time"])
        print(f"[2/3] Stage 1: SKIPPED -- loaded cached candidates: "
              f"{len(cand_df):,} rows from {args.candidates_out}")
        # jump to stage 2 below
    else:
        # Build per-day work units
        print(f"[2/3] Stage 1: envelope-peak detection at {args.station} "
              f"across PNSN windows (workers={args.workers})...")
        pnsn = pnsn.sort_values("time").reset_index(drop=True)
        pnsn["t0_us"] = pnsn["time"].astype("datetime64[ns]").astype("int64") // 1000
        pnsn["t1_us"] = pnsn["t0_us"] + int(args.window_s * 1_000_000)
        pnsn["date"] = pnsn["time"].dt.date.astype(str)

        work = {}
        for _, r in pnsn.iterrows():
            work.setdefault(r["date"], []).append(
                ((r["t0_us"], r["t1_us"]), (r["lat"], r["lon"]))
            )
        cfg = dict(
            wfdir=args.wfdir, station=args.station, fs=args.fs,
            fmin=args.fmin, fmax=args.fmax,
            snr_threshold=args.snr_threshold, peak_sep_s=args.peak_sep_s,
        )
        tasks = []
        for day_str, items in work.items():
            windows_us = [a for a, _ in items]
            latlons = [b for _, b in items]
            tasks.append((day_str, windows_us, latlons, cfg))

        candidates = []
        n_done = 0
        t0 = time.time()
        with ProcessPoolExecutor(max_workers=args.workers) as ex:
            futs = {ex.submit(_candidates_one_day, t): t[0] for t in tasks}
            for f in as_completed(futs):
                day_str, rows = f.result()
                candidates.extend(rows)
                n_done += 1
                if n_done % 200 == 0:
                    el = time.time() - t0
                    print(f"  day {n_done}/{len(tasks)}, "
                          f"cumulative candidates: {len(candidates):,}  ({el:.0f}s)")

        cand_df = pd.DataFrame(candidates)
        cand_df.to_parquet(args.candidates_out, index=False)
        print(f"  saved {len(cand_df):,} candidates -> {args.candidates_out}")

    # Stage 2: cluster per spatial bin
    print(f"[3/3] Stage 2: complete-linkage clustering per "
          f"{args.grid_deg} deg bin at CC>={args.cc_threshold}...")
    cand_df["lat_bin"] = np.round(cand_df["lat"] / args.grid_deg) * args.grid_deg
    cand_df["lon_bin"] = np.round(cand_df["lon"] / args.grid_deg) * args.grid_deg
    cand_df["bin_label"] = (cand_df["lat_bin"].map(lambda v: f"{v:.3f}")
                            + "_"
                            + cand_df["lon_bin"].map(lambda v: f"{v:.3f}"))
    cand_df["time_us"] = (cand_df["time"]
                          .astype("datetime64[ns]").astype("int64") // 1000)

    bin_counts = cand_df["bin_label"].value_counts()
    keep_bins = bin_counts[bin_counts >= args.min_bin_candidates].index.tolist()
    print(f"  {len(bin_counts)} bins with candidates, "
          f"{len(keep_bins)} with >={args.min_bin_candidates}")

    rng = np.random.default_rng(args.seed)
    cluster_cfg = dict(
        wfdir=args.wfdir, station=args.station, fs=args.fs,
        fmin=args.fmin, fmax=args.fmax,
        template_pre=args.template_pre, template_post=args.template_post,
        max_shift_samples=args.max_shift_samples,
        cc_threshold=args.cc_threshold,
        min_family_members=args.min_family_members,
        min_years=args.min_years,
    )
    bin_tasks = []
    for bin_label in keep_bins:
        sub = cand_df[cand_df["bin_label"] == bin_label]
        if len(sub) > args.max_bin_candidates:
            sub = sub.iloc[rng.choice(len(sub), args.max_bin_candidates, replace=False)]
        bin_tasks.append((bin_label,
                          np.array(sub["time_us"].values, dtype=np.int64),
                          cluster_cfg))

    families = []
    n_done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_process_bin, t): t[0] for t in bin_tasks}
        for f in as_completed(futs):
            bin_label, clusters = f.result()
            n_done += 1
            if clusters:
                for c in clusters:
                    c["lat"] = float(cand_df[cand_df["bin_label"] == bin_label]["lat"].mean())
                    c["lon"] = float(cand_df[cand_df["bin_label"] == bin_label]["lon"].mean())
                    families.append(c)
            if n_done % 20 == 0:
                print(f"  bin {n_done}/{len(bin_tasks)}, "
                      f"families so far: {len(families)}")

    print(f"\n[done] Discovered {len(families)} NLLB families at "
          f"CC>={args.cc_threshold} complete-linkage, "
          f">={args.min_family_members} members, "
          f">={args.min_years} years.")

    # Save
    npz_arrays = {}
    summary_rows = []
    member_rows = []
    for f in families:
        npz_arrays[f["family_id"]] = f["template"]
        summary_rows.append(dict(
            family_id=f["family_id"], n_members=f["n_members"],
            n_years=f["n_years"], first_year=f["first_year"],
            last_year=f["last_year"], lat=f["lat"], lon=f["lon"],
        ))
        for t in f["member_times"]:
            member_rows.append(dict(family_id=f["family_id"], time=t))
    np.savez(args.out, **npz_arrays)
    pd.DataFrame(summary_rows).to_csv(
        Path(args.out).with_suffix(".summary.csv"), index=False)
    pd.DataFrame(member_rows).to_parquet(
        Path(args.out).with_suffix(".members.parquet"), index=False)
    print(f"  saved -> {args.out} + summary + members")


if __name__ == "__main__":
    main()
