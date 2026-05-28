"""Extract per-family NLLB templates using PGC matched-filter detection times.

For each of the 51 PGC LFE families:
  1. Pull N high-CC PGC detection times from mf_pgc_all51_cc08.csv.
  2. At NLLB, cut a wide window (+/- 10 s) around each PGC time, bandpass
     2-8 Hz, find the Hilbert-envelope peak within +/- 5 s of the PGC
     time -- that's the NLLB direct-S arrival.
  3. Cut a tight 2-s window centered on that peak at NLLB.
  4. L2-normalize each cut.
  5. Reject cuts that look like band-limited noise (peak/median envelope
     ratio < 2) -- the family may not be visible at NLLB.
  6. Stack the surviving cuts -> NLLB template per family.
  7. Estimate the median NLLB-vs-PGC arrival offset per family (useful
     diagnostic + needed for downstream MF).

Output: data/nllb_templates.npz with arrays per family (template, offset_s,
n_used, n_attempted, peak_ratio_median).

Caveat: this finds the NLLB template that's coherent for sources known
to PGC. If a family is *not* visible at NLLB (poor SNR / geometry),
the stack will be flat noise and gets flagged. Independent PNSN-driven
discovery at NLLB is a separate step that catches NLLB-strong sources
that PGC missed.
"""

from __future__ import annotations

import argparse
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, "src")


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--mf-csv", default="data/mf_pgc_all51_cc08.csv")
    p.add_argument("--wfdir", default="data/waveforms")
    p.add_argument("--station", default="NLLB")
    p.add_argument("--network", default="CN")
    p.add_argument("--out", default="data/nllb_templates.npz")
    p.add_argument("--n-per-family", type=int, default=200,
                   help="Random sample of PGC detections per family")
    p.add_argument("--pre-pad", type=float, default=10.0,
                   help="Seconds before PGC time to load")
    p.add_argument("--post-pad", type=float, default=10.0,
                   help="Seconds after PGC time to load")
    p.add_argument("--peak-search-s", type=float, default=5.0,
                   help="Envelope-peak search half-width around PGC time")
    p.add_argument("--template-s", type=float, default=2.0,
                   help="Output template duration")
    p.add_argument("--fs", type=float, default=40.0)
    p.add_argument("--fmin", type=float, default=2.0)
    p.add_argument("--fmax", type=float, default=8.0)
    p.add_argument("--min-peak-ratio", type=float, default=2.0,
                   help="Reject cuts where peak/median envelope < this")
    p.add_argument("--workers", type=int, default=24)
    p.add_argument("--cc-high", type=float, default=0.85,
                   help="Only use PGC detections with cc>=this for templates")
    return p.parse_args()


# Globals for ProcessPool workers (avoid re-pickling per task).
_CFG: dict = {}


def _load_day_filt(root: Path, station: str, network: str,
                   day: datetime, fs: float, fmin: float, fmax: float):
    from obspy import read

    candidates = list(root.glob(
        f"{network}.{station}/{day.year}/{day.timetuple().tm_yday:03d}.mseed"
    ))
    if not candidates:
        return None
    try:
        st = read(str(candidates[0]))
    except Exception:
        return None
    st = st.select(component="Z")
    if len(st) == 0:
        return None
    for tr in st:
        if abs(tr.stats.sampling_rate - fs) > 1e-6:
            try:
                tr.resample(fs)
            except Exception:
                return None
    try:
        st.merge(fill_value=0)
    except Exception:
        return None
    if len(st) == 0:
        return None
    st.detrend("demean")
    st.filter("bandpass", freqmin=fmin, freqmax=fmax,
              corners=4, zerophase=True)
    tr = st[0]
    return tr.data.astype(np.float32), tr.stats.starttime.datetime


def _cut_one_family(args):
    family, times_us, cfg = args
    from scipy.signal import hilbert

    fs = cfg["fs"]
    template_n = int(round(cfg["template_s"] * fs))
    pre_n = int(round(cfg["pre_pad"] * fs))
    post_n = int(round(cfg["post_pad"] * fs))
    peak_half_n = int(round(cfg["peak_search_s"] * fs))
    template_half_n = template_n // 2

    cuts: list[np.ndarray] = []
    offsets_s: list[float] = []
    peak_ratios: list[float] = []
    n_attempted = 0

    # Sort by day so we load each day at most once.
    times_dt = [datetime.utcfromtimestamp(tus / 1e6) for tus in times_us]
    order = np.argsort([t.timetuple().tm_yday + 1000 * t.year for t in times_dt])
    times_dt = [times_dt[i] for i in order]

    cur_day = None
    cur_loaded = None

    for dt in times_dt:
        day = datetime(dt.year, dt.month, dt.day)
        if day != cur_day:
            cur_loaded = _load_day_filt(
                Path(cfg["wfdir"]), cfg["station"], cfg["network"],
                day, fs, cfg["fmin"], cfg["fmax"],
            )
            cur_day = day
        if cur_loaded is None:
            continue
        data, start = cur_loaded
        offset_s_pgc = (dt - start).total_seconds()
        i_pgc = int(round(offset_s_pgc * fs))
        i0 = i_pgc - pre_n
        i1 = i_pgc + post_n
        if i0 < 0 or i1 > data.size:
            continue
        seg = data[i0:i1]
        env = np.abs(hilbert(seg))
        if env.size == 0:
            continue
        # search for envelope peak in +/- peak_search_s around the PGC time
        center = pre_n  # index of PGC time inside seg
        lo = max(0, center - peak_half_n)
        hi = min(env.size, center + peak_half_n)
        if hi - lo < 1:
            continue
        local_peak_idx = lo + int(np.argmax(env[lo:hi]))
        local_peak_val = float(env[local_peak_idx])
        med = float(np.median(env)) + 1e-12
        peak_ratio = local_peak_val / med
        n_attempted += 1
        if peak_ratio < cfg["min_peak_ratio"]:
            continue
        # cut template_n window centered on the peak
        t0 = local_peak_idx - template_half_n
        t1 = t0 + template_n
        if t0 < 0 or t1 > seg.size:
            continue
        cut = seg[t0:t1].astype(np.float64)
        nrm = float(np.linalg.norm(cut))
        if nrm == 0 or not np.isfinite(nrm):
            continue
        cuts.append((cut / nrm).astype(np.float32))
        offsets_s.append((local_peak_idx - center) / fs)
        peak_ratios.append(peak_ratio)

    if not cuts:
        return family, None, 0, n_attempted, []

    C = np.array(cuts)
    # build template by averaging, then re-stack only the cuts that match
    # the average -- this rejects misaligned outliers.
    mean = C.mean(axis=0)
    mean = mean / (np.linalg.norm(mean) + 1e-12)
    ccs = C @ mean  # shape (n,)
    keep_mask = ccs > 0.5
    if keep_mask.sum() < 5:
        # too few coherent cuts -- give up for this family
        return family, None, int(keep_mask.sum()), n_attempted, peak_ratios
    template = C[keep_mask].mean(axis=0)
    template = (template / (np.linalg.norm(template) + 1e-12)).astype(np.float32)
    median_offset = float(np.median(np.array(offsets_s)[keep_mask]))
    return family, template, int(keep_mask.sum()), n_attempted, peak_ratios


def main():
    args = parse_args()

    print(f"[1/3] Loading PGC detections (cc>={args.cc_high})...")
    df = pd.read_csv(args.mf_csv)
    df = df[df["cc"] >= args.cc_high]
    df["time"] = pd.to_datetime(df["time"], format="mixed")
    df["t_us"] = df["time"].astype("datetime64[ns]").astype("int64") // 1000
    print(f"  {len(df):,} high-cc PGC detections across {df['template'].nunique()} families")

    rng = np.random.default_rng(0)
    cfg = dict(
        wfdir=args.wfdir, station=args.station, network=args.network,
        fs=args.fs, fmin=args.fmin, fmax=args.fmax,
        pre_pad=args.pre_pad, post_pad=args.post_pad,
        peak_search_s=args.peak_search_s, template_s=args.template_s,
        min_peak_ratio=args.min_peak_ratio,
    )

    tasks = []
    for fam, group in df.groupby("template"):
        if len(group) > args.n_per_family:
            sub = group.iloc[rng.choice(len(group), args.n_per_family, replace=False)]
        else:
            sub = group
        times_us = np.array(sub["t_us"].values, dtype=np.int64)
        tasks.append((fam, times_us, cfg))
    print(f"  prepared {len(tasks)} families, sampling up to {args.n_per_family} detections each")

    print(f"[2/3] Building NLLB templates (workers={args.workers})...")
    results: dict[str, dict] = {}
    n_done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_cut_one_family, t): t[0] for t in tasks}
        for f in as_completed(futs):
            fam, tmpl, n_kept, n_att, peak_ratios = f.result()
            n_done += 1
            results[fam] = dict(
                template=tmpl, n_kept=n_kept, n_attempted=n_att,
                peak_ratio_median=float(np.median(peak_ratios)) if peak_ratios else 0.0,
            )
            status = "OK" if tmpl is not None else "skip"
            print(f"  [{n_done}/{len(tasks)}] {fam}: {status}  "
                  f"kept={n_kept}/{n_att}  pr_med={results[fam]['peak_ratio_median']:.2f}")

    print(f"[3/3] Saving to {args.out}...")
    out_arrays: dict[str, np.ndarray] = {}
    meta_rows = []
    for fam in sorted(results):
        r = results[fam]
        if r["template"] is not None:
            out_arrays[fam] = r["template"]
        meta_rows.append(dict(
            family=fam, n_kept=r["n_kept"], n_attempted=r["n_attempted"],
            peak_ratio_median=r["peak_ratio_median"],
            template_present=r["template"] is not None,
        ))
    np.savez(args.out, **out_arrays)
    meta = pd.DataFrame(meta_rows)
    meta_csv = Path(args.out).with_suffix(".meta.csv")
    meta.to_csv(meta_csv, index=False)
    print(f"  {sum(1 for r in results.values() if r['template'] is not None)} "
          f"families with NLLB templates, {len(results)} total")
    print(f"  saved {args.out} + {meta_csv}")


if __name__ == "__main__":
    main()
