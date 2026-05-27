"""Long-window daily stacks for one LFE family, plotted as an image.

For each day with >= min_det LFE detections (after CC filtering), cut a long
window (-pre to +post seconds) around each detection from the bandpassed
day waveform, sum & L2-normalize within the day, and store one row per day.

Output: image plot with
    x = trace time relative to LFE alignment (0 = detection time, same axis
        the stretching dv/v reference uses)
    y = date
    color = normalized amplitude

This is the LFE analog of the classic ambient-noise NCF image-of-time plot.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--template", default="PGC_79")
    p.add_argument("--station", default="PGC")
    p.add_argument("--network", default="CN")
    p.add_argument("--pre", type=float, default=3.0)
    p.add_argument("--post", type=float, default=10.0)
    p.add_argument("--fs", type=float, default=40.0)
    p.add_argument("--fmin", type=float, default=2.0)
    p.add_argument("--fmax", type=float, default=8.0)
    p.add_argument("--cc-min", type=float, default=0.80)
    p.add_argument("--min-det", type=int, default=20,
                   help="Min detections per day to keep that day")
    p.add_argument("--workers", type=int, default=24)
    p.add_argument("--mf-csv", default="data/mf_pgc_2005-2026_cc08.csv")
    p.add_argument("--wfdir", default="data/waveforms")
    p.add_argument("--out", default="figures/long_window_daily_image_PGC_79.png")
    p.add_argument("--npz", default="data/long_window_daily_PGC_79.npz")
    return p.parse_args()


# Workers must reach top-level args via globals (ProcessPool pickling).
_WORKER_CFG: dict = {}


def _set_cfg(cfg: dict):
    _WORKER_CFG.update(cfg)


def _process_day(args_tuple):
    day_str, times_us, cfg = args_tuple
    from obspy import read

    wfdir = Path(cfg["wfdir"])
    network = cfg["network"]
    station = cfg["station"]
    fs = cfg["fs"]
    fmin = cfg["fmin"]
    fmax = cfg["fmax"]
    pre_n = cfg["pre_n"]
    post_n = cfg["post_n"]
    win_n = pre_n + post_n

    day = datetime.fromisoformat(day_str)
    p = wfdir / f"{network}.{station}" / f"{day.year}" / f"{day.timetuple().tm_yday:03d}.mseed"
    if not p.exists():
        return day_str, None, 0
    try:
        st = read(str(p))
    except Exception:
        return day_str, None, 0
    st = st.select(component="Z")
    if len(st) == 0:
        return day_str, None, 0
    for tr in st:
        if abs(tr.stats.sampling_rate - fs) > 1e-6:
            try:
                tr.resample(fs)
            except Exception:
                return day_str, None, 0
    try:
        st.merge(fill_value=0)
    except Exception:
        return day_str, None, 0
    st.detrend("demean")
    st.filter("bandpass", freqmin=fmin, freqmax=fmax, corners=4, zerophase=True)
    tr = st[0]
    data = tr.data.astype(np.float32)
    start = tr.stats.starttime.datetime

    sum_win = np.zeros(win_n, dtype=np.float64)
    n_used = 0
    for tus in times_us:
        dt = datetime.utcfromtimestamp(tus / 1e6)
        offset_s = (dt - start).total_seconds()
        i0 = int(round(offset_s * fs)) - pre_n
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
        return day_str, None, n_used
    daily = sum_win / n_used
    nrm = float(np.linalg.norm(daily))
    if nrm > 0:
        daily = daily / nrm
    return day_str, daily.astype(np.float32), n_used


def main():
    args = parse_args()

    pre_n = int(round(args.pre * args.fs))
    post_n = int(round(args.post * args.fs))
    win_n = pre_n + post_n
    t_axis = (np.arange(win_n) - pre_n) / args.fs

    print(f"[1/4] Loading MF detections for {args.template} (cc>={args.cc_min})...")
    df = pd.read_csv(args.mf_csv)
    df = df[(df["template"] == args.template) & (df["cc"] >= args.cc_min)]
    df["time"] = pd.to_datetime(df["time"], format="mixed")
    df["date"] = df["time"].dt.date.astype(str)
    print(f"  {len(df):,} detections across {df['date'].nunique()} days")

    # group: list of microsecond timestamps per day
    by_day = df.groupby("date")["time"].apply(
        lambda s: np.array([t.value // 1000 for t in s], dtype=np.int64)
    ).to_dict()
    day_list = sorted(by_day.keys())
    print(f"[2/4] Building daily long-window stacks (workers={args.workers})...")

    cfg = dict(
        wfdir=args.wfdir, network=args.network, station=args.station,
        fs=args.fs, fmin=args.fmin, fmax=args.fmax,
        pre_n=pre_n, post_n=post_n, min_det=args.min_det,
    )
    arg_iter = [(d, by_day[d], cfg) for d in day_list]

    daily_arr = []
    daily_dates = []
    daily_n = []
    n_done = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futs = {ex.submit(_process_day, a): a[0] for a in arg_iter}
        for f in as_completed(futs):
            d, w, nu = f.result()
            n_done += 1
            if n_done % 200 == 0:
                print(f"  {n_done}/{len(arg_iter)} days, kept {len(daily_arr)}")
            if w is None:
                continue
            daily_arr.append(w)
            daily_dates.append(d)
            daily_n.append(nu)

    # sort by date
    order = np.argsort(daily_dates)
    daily_arr = np.array([daily_arr[i] for i in order])
    daily_dates = np.array([daily_dates[i] for i in order])
    daily_n = np.array([daily_n[i] for i in order])
    print(f"  kept {len(daily_arr)} day-stacks")

    np.savez(args.npz, stacks=daily_arr, dates=daily_dates, n_det=daily_n,
             t=t_axis, fs=args.fs, template=args.template)
    print(f"  saved {args.npz}")

    print("[3/4] All-time reference stack...")
    ns = daily_n.astype(np.float64)
    total = ns.sum()
    weighted = (daily_arr.T * ns).sum(axis=1) / total
    nrm = float(np.linalg.norm(weighted))
    ref = (weighted / nrm).astype(np.float32) if nrm > 0 else weighted.astype(np.float32)

    print("[4/4] Plotting...")
    dates = pd.to_datetime(daily_dates)
    fig = plt.figure(figsize=(11, 9))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 6], hspace=0.05)

    ax0 = fig.add_subplot(gs[0])
    ax0.plot(t_axis, ref, color="k", lw=0.8)
    ax0.axvline(0, color="C3", lw=0.5, alpha=0.5)
    ax0.axvspan(0, 2, color="grey", alpha=0.12)
    ax0.set_xlim(t_axis[0], t_axis[-1])
    ax0.set_ylabel("ref (a.u.)")
    ax0.set_title(
        f"{args.template} ({args.network}.{args.station}): {len(daily_arr)} daily L2-normed stacks "
        f"({df['time'].min().date()} – {df['time'].max().date()}), "
        f"bandpass {args.fmin}-{args.fmax} Hz, cc>={args.cc_min}"
    )
    plt.setp(ax0.get_xticklabels(), visible=False)

    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    # convert dates to ordinal floats for imshow extents
    date_nums = mdates.date2num(dates.to_pydatetime())
    # we plot with imshow using row index, then map ticks to dates
    vmax = float(np.percentile(np.abs(daily_arr), 99))
    im = ax1.imshow(
        daily_arr,
        aspect="auto",
        cmap="RdBu_r",
        vmin=-vmax, vmax=vmax,
        extent=[t_axis[0], t_axis[-1], date_nums[-1], date_nums[0]],
        interpolation="nearest",
    )
    ax1.yaxis_date()
    ax1.yaxis.set_major_locator(mdates.YearLocator(2))
    ax1.yaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax1.set_ylabel("date")
    ax1.set_xlabel("trace time relative to LFE alignment (s)  [stretching window: 0–2 s]")
    ax1.axvline(0, color="k", lw=0.4, alpha=0.4)
    cbar = fig.colorbar(im, ax=ax1, pad=0.02)
    cbar.set_label("amplitude (L2-normed)")

    fig.tight_layout()
    outp = Path(args.out)
    outp.parent.mkdir(exist_ok=True, parents=True)
    fig.savefig(outp, dpi=150)
    print(f"  wrote {outp}")


if __name__ == "__main__":
    main()
