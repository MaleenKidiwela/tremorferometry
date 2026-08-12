"""Build long-window (-3 to +10 s) daily stacks for ALL 51 LFE families at PGC.

For each PGC day file:
  - Load + bandpass 2-8 Hz ONCE.
  - For each family that has any CC>=0.8 detections on that day:
      cut all detections, L2-norm each, sum, divide by count, L2-norm result.
  - Write one row per (family, day) where the count is >= min_det.

Templates with detections are read from `data/mf_pgc_2005-2026_cc08.csv`.
Outputs an .npz with stacks (n_rows, n_samp), patches (n_rows,), dates,
and n_det.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--mf-csv", default="data/mf_pgc_2005-2026_cc08.csv")
    p.add_argument("--wfdir", default="data/waveforms")
    p.add_argument("--network", default="CN")
    p.add_argument("--station", default="PGC")
    p.add_argument("--pre", type=float, default=3.0)
    p.add_argument("--post", type=float, default=10.0)
    p.add_argument("--fs", type=float, default=40.0)
    p.add_argument("--fmin", type=float, default=2.0)
    p.add_argument("--fmax", type=float, default=8.0)
    p.add_argument("--cc-min", type=float, default=0.80)
    p.add_argument("--min-det", type=int, default=20)
    p.add_argument("--workers", type=int, default=32)
    p.add_argument("--out", default="data/long_window_daily_all51.npz")
    return p.parse_args()


def _process_day(args_tuple):
    day_str, by_tmpl, cfg = args_tuple
    from obspy import read

    wfdir = Path(cfg["wfdir"])
    day = datetime.fromisoformat(day_str)
    p = wfdir / f"{cfg['network']}.{cfg['station']}" / f"{day.year}" \
        / f"{day.timetuple().tm_yday:03d}.mseed"
    if not p.exists():
        return day_str, []
    try:
        st = read(str(p))
    except Exception:
        return day_str, []
    st = st.select(component="Z")
    if len(st) == 0:
        return day_str, []
    for tr in st:
        if abs(tr.stats.sampling_rate - cfg["fs"]) > 1e-6:
            try:
                tr.resample(cfg["fs"])
            except Exception:
                return day_str, []
    try:
        st.merge(fill_value=0)
    except Exception:
        return day_str, []
    st.detrend("demean")
    st.filter("bandpass", freqmin=cfg["fmin"], freqmax=cfg["fmax"],
              corners=4, zerophase=True)
    tr = st[0]
    data = tr.data.astype(np.float32)
    start = tr.stats.starttime.datetime

    pre_n = cfg["pre_n"]
    post_n = cfg["post_n"]
    win_n = pre_n + post_n

    out_rows = []
    for tmpl, times_us in by_tmpl.items():
        sum_win = np.zeros(win_n, dtype=np.float64)
        n_used = 0
        for tus in times_us:
            dt = datetime.utcfromtimestamp(tus / 1e6)
            offset_s = (dt - start).total_seconds()
            i0 = int(round(offset_s * cfg["fs"])) - pre_n
            i1 = i0 + win_n
            if i0 < 0 or i1 > data.size:
                continue
            seg = data[i0:i1].astype(np.float64)
            nrm = float(np.linalg.norm(seg))
            if nrm == 0 or not np.isfinite(nrm):
                continue
            sum_win += seg / nrm
            n_used += 1
        if n_used < cfg["min_det"]:
            continue
        daily = sum_win / n_used
        nrm = float(np.linalg.norm(daily))
        if nrm <= 0:
            continue
        daily = (daily / nrm).astype(np.float32)
        out_rows.append((tmpl, daily, n_used))
    return day_str, out_rows


def main():
    args = parse_args()

    pre_n = int(round(args.pre * args.fs))
    post_n = int(round(args.post * args.fs))
    win_n = pre_n + post_n
    t_axis = (np.arange(win_n) - pre_n) / args.fs

    print(f"[1/3] Loading MF detections (cc>={args.cc_min})...")
    df = pd.read_csv(args.mf_csv)
    df = df[df["cc"] >= args.cc_min]
    df["time"] = pd.to_datetime(df["time"], format="mixed")
    df["date"] = df["time"].dt.date.astype(str)
    n_tmpl = df["template"].nunique()
    print(f"  {len(df):,} detections across {df['date'].nunique()} days "
          f"and {n_tmpl} templates")

    # Group: per day -> {template: array of microsecond timestamps}
    # NB: this is the data structure handed to each worker.
    print("[2/3] Grouping detections by (day, template)...")
    # pandas may return datetime64[us] from to_datetime; cast to ns explicitly
    df["t_us"] = df["time"].astype("datetime64[ns]").astype("int64") // 1000
    grouped = df.groupby(["date", "template"])["t_us"].apply(
        lambda s: np.array(s.values, dtype=np.int64)
    )
    day_list = sorted({d for d, _ in grouped.index})
    by_day_tmpl: dict[str, dict[str, np.ndarray]] = {d: {} for d in day_list}
    for (d, t), arr in grouped.items():
        by_day_tmpl[d][t] = arr
    print(f"  {len(day_list)} days to process")

    cfg = dict(
        wfdir=args.wfdir, network=args.network, station=args.station,
        fs=args.fs, fmin=args.fmin, fmax=args.fmax,
        pre_n=pre_n, post_n=post_n, min_det=args.min_det,
    )
    arg_iter = [(d, by_day_tmpl[d], cfg) for d in day_list]

    print(f"[3/3] Cutting daily long-window stacks (workers={args.workers})...")
    stacks: list[np.ndarray] = []
    patches: list[str] = []
    dates: list[str] = []
    n_det: list[int] = []

    n_done = 0
    n_kept = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_process_day, a): a[0] for a in arg_iter}
        for f in as_completed(futs):
            d, rows = f.result()
            n_done += 1
            if rows:
                for tmpl, w, nu in rows:
                    stacks.append(w)
                    patches.append(tmpl)
                    dates.append(d)
                    n_det.append(nu)
                    n_kept += 1
            if n_done % 200 == 0:
                print(f"  day {n_done}/{len(arg_iter)}, kept {n_kept} (template,day) rows")

    stacks_arr = np.array(stacks, dtype=np.float32)
    patches_arr = np.array(patches)
    dates_arr = np.array(dates)
    n_det_arr = np.array(n_det, dtype=np.int32)

    # Sort by (patch, date) for downstream convenience
    order = np.lexsort((dates_arr, patches_arr))
    stacks_arr = stacks_arr[order]
    patches_arr = patches_arr[order]
    dates_arr = dates_arr[order]
    n_det_arr = n_det_arr[order]

    np.savez(args.out,
             stacks=stacks_arr, patches=patches_arr, dates=dates_arr,
             n_det=n_det_arr, t=t_axis, fs=args.fs)
    print(f"  saved {args.out}  ({stacks_arr.shape}, {n_kept} rows)")


if __name__ == "__main__":
    main()
