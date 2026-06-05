"""Build long-window (-3..+10 s) daily stacks WITH instrument-response removal.

vs build_long_window_daily_all51.py:
  * remove_response (per trace, at native rate, output=VEL, pre_filt) using the saved StationXML
    -> corrects the BHZ->HHZ / sample-rate / response changes that caused the dv/v steps
  * resample_poly (fast) instead of obspy FFT resample
  * spawn workers (NOT fork) so they don't COW the parent's big detection DataFrame (the old
    fork stacker ballooned to ~135 GB; spawn keeps it ~parent + small per-worker)
Otherwise identical: per (family, day), cut all detections, L2-norm each, sum, normalize.
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from math import gcd
from pathlib import Path
import numpy as np, pandas as pd
from scipy.signal import resample_poly

_INV = None
def _inv(path):
    global _INV
    if _INV is None:
        from obspy import read_inventory
        _INV = read_inventory(path)
    return _INV

def _process_day(args_tuple):
    day_str, by_tmpl, cfg = args_tuple
    from obspy import read
    wf = Path(cfg["wfdir"]); day = datetime.fromisoformat(day_str)
    p = wf / f"{cfg['network']}.{cfg['station']}" / f"{day.year}" / f"{day.timetuple().tm_yday:03d}.mseed"
    if not p.exists():
        return day_str, []
    try:
        st = read(str(p)).select(component="Z")
    except Exception:
        return day_str, []
    if len(st) == 0:
        return day_str, []
    fs = cfg["fs"]
    def _resamp():
        for tr in st:
            fs0 = tr.stats.sampling_rate
            if abs(fs0 - fs) > 1e-6:
                r = int(round(fs0)); g = gcd(int(fs), r)
                tr.data = resample_poly(tr.data.astype(np.float64), int(fs)//g, r//g)
                tr.stats.sampling_rate = float(fs)
    try:
        if cfg.get("no_deconv"):
            # resample_poly ONLY (no response removal). Fine for WITHIN-era dv/v (the resample_poly
            # fix is what kills the obspy artifact; deconv only matters for CROSS-era comparison).
            _resamp()
        elif cfg.get("deconv_native"):
            # CORRECT order: deconvolve at NATIVE rate (true per-instrument phase), then resample.
            st.remove_response(inventory=_inv(cfg["inv_path"]), output="VEL", pre_filt=cfg["pre_filt"], water_level=60)
            _resamp()
        else:
            # FAST order: resample to 40 Hz first, then deconvolve (~2.75x faster).
            _resamp()
            st.remove_response(inventory=_inv(cfg["inv_path"]), output="VEL", pre_filt=cfg["pre_filt"], water_level=60)
    except Exception:
        return day_str, []
    try:
        st.merge(fill_value=0)
    except Exception:
        return day_str, []
    if len(st) == 0:
        return day_str, []
    st.detrend("demean")
    st.filter("bandpass", freqmin=cfg["fmin"], freqmax=cfg["fmax"], corners=4, zerophase=True)
    tr = st[0]; data = tr.data.astype(np.float64); start = tr.stats.starttime.datetime
    _dm = cfg.get("despike_mad", 0.0)
    if _dm and _dm > 0:
        _med = np.median(data); _mad = np.median(np.abs(data - _med)) * 1.4826
        if _mad > 0:
            np.clip(data, _med - _dm*_mad, _med + _dm*_mad, out=data)
    pre_n, post_n = cfg["pre_n"], cfg["post_n"]; win_n = pre_n + post_n
    out = []
    for tmpl, times_us in by_tmpl.items():
        acc = np.zeros(win_n); nu = 0
        for tus in times_us:
            off = (datetime.utcfromtimestamp(tus/1e6) - start).total_seconds()
            i0 = int(round(off*cfg["fs"])) - pre_n; i1 = i0 + win_n
            if i0 < 0 or i1 > data.size:
                continue
            seg = data[i0:i1]; nrm = np.linalg.norm(seg)
            if nrm == 0 or not np.isfinite(nrm):
                continue
            acc += seg/nrm; nu += 1
        if nu < cfg["min_det"]:
            continue
        daily = acc/nu; nrm = np.linalg.norm(daily)
        if nrm <= 0:
            continue
        out.append((tmpl, (daily/nrm).astype(np.float32), nu))
    return day_str, out

def main():
    import multiprocessing as mp
    p = argparse.ArgumentParser()
    p.add_argument("--mf-csv", default="data/mf_gnw_all.csv")
    p.add_argument("--wfdir", default="data/waveforms")
    p.add_argument("--network", default="UW"); p.add_argument("--station", default="GNW")
    p.add_argument("--inv", default="data/UW.GNW.response.xml")
    p.add_argument("--pre", type=float, default=3.0); p.add_argument("--post", type=float, default=10.0)
    p.add_argument("--fs", type=float, default=40.0)
    p.add_argument("--fmin", type=float, default=2.0); p.add_argument("--fmax", type=float, default=8.0)
    p.add_argument("--cc-min", type=float, default=0.80); p.add_argument("--min-det", type=int, default=20)
    p.add_argument("--workers", type=int, default=22)
    p.add_argument("--deconv-native", action="store_true",
                   help="deconvolve at native sample rate BEFORE resampling (correct, slower)")
    p.add_argument("--no-deconv", action="store_true",
                   help="skip response removal entirely (resample_poly only; OK for within-era dv/v)")
    p.add_argument("--start-year", type=int, default=0); p.add_argument("--end-year", type=int, default=9999)
    p.add_argument("--out", default="data/long_window_daily_GNW_resp.npz")
    p.add_argument("--despike-mad", type=float, default=0.0,
                   help="winsorize bandpassed samples beyond this many MAD/day (0=off; "
                        "match the densify --despike-mad so detection and stacking agree)")
    args = p.parse_args()

    pre_n = int(round(args.pre*args.fs)); post_n = int(round(args.post*args.fs)); win_n = pre_n+post_n
    t_axis = (np.arange(win_n)-pre_n)/args.fs
    print(f"[1/3] load detections {args.mf_csv} (cc>={args.cc_min})...", flush=True)
    df = pd.read_csv(args.mf_csv); df = df[df["cc"] >= args.cc_min]
    df["time"] = pd.to_datetime(df["time"], format="mixed")
    df = df[(df["time"].dt.year >= args.start_year) & (df["time"].dt.year <= args.end_year)]
    df["date"] = df["time"].dt.date.astype(str)
    df["t_us"] = df["time"].astype("datetime64[ns]").astype("int64")//1000
    grouped = df.groupby(["date","template"])["t_us"].apply(lambda s: np.array(s.values, dtype=np.int64))
    days = sorted({d for d,_ in grouped.index})
    bydt = {d: {} for d in days}
    for (d,t),arr in grouped.items(): bydt[d][t] = arr
    print(f"  {len(df):,} det, {len(days)} days", flush=True)
    cfg = dict(wfdir=args.wfdir, network=args.network, station=args.station, inv_path=args.inv,
               fs=args.fs, fmin=args.fmin, fmax=args.fmax, pre_n=pre_n, post_n=post_n,
               min_det=args.min_det, pre_filt=[0.5,1.0,15.0,18.0], deconv_native=args.deconv_native,
               no_deconv=args.no_deconv, despike_mad=args.despike_mad)
    arg_iter = [(d, bydt[d], cfg) for d in days]

    print(f"[2/3] stacking (resp-removed, spawn, workers={args.workers})...", flush=True)
    ctx = mp.get_context("spawn")
    stacks=[]; patches=[]; dates=[]; ndet=[]; done=0; kept=0
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx) as ex:
        futs = [ex.submit(_process_day, a) for a in arg_iter]
        for f in as_completed(futs):
            d, rows = f.result(); done += 1
            for tmpl,w,nu in rows:
                stacks.append(w); patches.append(tmpl); dates.append(d); ndet.append(nu); kept += 1
            if done % 500 == 0:
                print(f"  {done}/{len(arg_iter)} days, {kept} (fam,day) rows", flush=True)
    stacks=np.array(stacks,np.float32); patches=np.array(patches); dates=np.array(dates); ndet=np.array(ndet,np.int32)
    order = np.lexsort((dates, patches))
    np.savez(args.out, stacks=stacks[order], patches=patches[order], dates=dates[order],
             n_det=ndet[order], t=t_axis, fs=args.fs)
    print(f"[3/3] saved {args.out} ({stacks.shape}, {kept} rows)", flush=True)

if __name__ == "__main__":
    main()
