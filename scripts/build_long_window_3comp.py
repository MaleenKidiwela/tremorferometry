"""Build long-window (-3..+10 s) daily stacks for ALL THREE components (Z, EH1, EH2)
from a densified detection catalog. Same per-(family,day) L2-normalized stacking as
build_long_window_resp.py, but reads the day file once and stacks each component at the
Z-derived detection times. Writes THREE npz (one per component), each in the exact format
dvv_roll30cal.py expects (stacks/patches/dates/n_det/t/fs).

resample_poly only (--no-deconv equivalent): correct for single-era boreholes like B011;
avoids needing a StationXML. Horizontals aligned by the Z detection time (S is coherent
across components at the same instant).

Usage: python scripts/build_long_window_3comp.py --mf-csv-glob 'data/mf_b011p90_*.csv' \
         --network PB --station B011 --fs 100 --out-prefix data/long_window_daily_B011p90
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, glob
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
from math import gcd
from pathlib import Path
import numpy as np, pandas as pd
from scipy.signal import resample_poly

CHANS = ["Z", "H1", "H2"]           # logical names
SELECT = {"Z": ["Z"], "H1": ["1", "N"], "H2": ["2", "E"]}


def _get(st, keys):
    for k in keys:
        s = st.select(channel=f"*{k}")
        if len(s):
            return s
    return None


def _process_day(args_tuple):
    day_str, by_tmpl, cfg = args_tuple
    from obspy import read
    wf = Path(cfg["wfdir"]); day = datetime.fromisoformat(day_str)
    p = wf / f"{cfg['network']}.{cfg['station']}" / f"{day.year}" / f"{day.timetuple().tm_yday:03d}.mseed"
    if not p.exists():
        return day_str, {}
    try:
        st0 = read(str(p))
    except Exception:
        return day_str, {}
    fs = cfg["fs"]; pre_n, post_n = cfg["pre_n"], cfg["post_n"]; win_n = pre_n + post_n
    comp_data = {}
    for cname in CHANS:
        sc = _get(st0, SELECT[cname])
        if sc is None:
            continue
        sc = sc.copy()
        try:
            for tr in sc:
                fs0 = tr.stats.sampling_rate
                if abs(fs0 - fs) > 1e-6:
                    r = int(round(fs0)); g = gcd(int(fs), r)
                    tr.data = resample_poly(tr.data.astype(np.float64), int(fs)//g, r//g)
                    tr.stats.sampling_rate = float(fs)
            sc.merge(fill_value=0)
            if len(sc) == 0:
                continue
            sc.detrend("demean")
            sc.filter("bandpass", freqmin=cfg["fmin"], freqmax=cfg["fmax"], corners=4, zerophase=True)
        except Exception:
            continue
        tr = sc[0]; data = tr.data.astype(np.float64)
        _dm = cfg.get("despike_mad", 0.0)
        if _dm and _dm > 0:
            _med = np.median(data); _mad = np.median(np.abs(data - _med)) * 1.4826
            if _mad > 0:
                np.clip(data, _med - _dm*_mad, _med + _dm*_mad, out=data)
        comp_data[cname] = (data, tr.stats.starttime.datetime)
    if "Z" not in comp_data:
        return day_str, {}

    out = {c: [] for c in CHANS}
    for tmpl, times_us in by_tmpl.items():
        for cname, (data, start) in comp_data.items():
            acc = np.zeros(win_n); nu = 0
            for tus in times_us:
                off = (datetime.utcfromtimestamp(tus/1e6) - start).total_seconds()
                i0 = int(round(off*fs)) - pre_n; i1 = i0 + win_n
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
            out[cname].append((tmpl, (daily/nrm).astype(np.float32), nu))
    return day_str, out


def main():
    import multiprocessing as mp
    ap = argparse.ArgumentParser()
    ap.add_argument("--mf-csv-glob", required=True, help="glob for densify per-year CSVs")
    ap.add_argument("--wfdir", default="data/waveforms")
    ap.add_argument("--network", default="PB"); ap.add_argument("--station", default="B011")
    ap.add_argument("--pre", type=float, default=3.0); ap.add_argument("--post", type=float, default=10.0)
    ap.add_argument("--fs", type=float, default=100.0)
    ap.add_argument("--fmin", type=float, default=2.0); ap.add_argument("--fmax", type=float, default=8.0)
    ap.add_argument("--cc-min", type=float, default=0.80); ap.add_argument("--min-det", type=int, default=20)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--despike-mad", type=float, default=8.0)
    ap.add_argument("--out-prefix", default="data/long_window_daily_B011p90")
    a = ap.parse_args()

    pre_n = int(round(a.pre*a.fs)); post_n = int(round(a.post*a.fs)); win_n = pre_n+post_n
    t_axis = (np.arange(win_n)-pre_n)/a.fs
    files = sorted(glob.glob(a.mf_csv_glob))
    print(f"[1/3] {len(files)} detection files (cc>={a.cc_min}); STREAMING per-year "
          f"(bounded memory)", flush=True)
    cfg = dict(wfdir=a.wfdir, network=a.network, station=a.station, fs=a.fs,
               fmin=a.fmin, fmax=a.fmax, pre_n=pre_n, post_n=post_n,
               min_det=a.min_det, despike_mad=a.despike_mad)
    acc = {c: dict(stacks=[], patches=[], dates=[], ndet=[]) for c in CHANS}
    ctx = mp.get_context("spawn")
    print(f"[2/3] stacking 3-comp (spawn, workers={a.workers})...", flush=True)
    for fi, f in enumerate(files):
        df = pd.read_csv(f, usecols=["template", "time", "cc"])
        df = df[df["cc"] >= a.cc_min]
        df["time"] = pd.to_datetime(df["time"], format="mixed")
        df["date"] = df["time"].dt.date.astype(str)
        df["t_us"] = df["time"].astype("datetime64[ns]").astype("int64")//1000
        grouped = df.groupby(["date", "template"])["t_us"].apply(lambda s: np.array(s.values, dtype=np.int64))
        days = sorted({d for d, _ in grouped.index})
        bydt = {d: {} for d in days}
        for (d, t), arr in grouped.items():
            bydt[d][t] = arr
        nrows = len(df); del df, grouped
        arg_iter = [(d, bydt[d], cfg) for d in days]
        done = 0
        with ProcessPoolExecutor(max_workers=a.workers, mp_context=ctx) as ex:
            futs = [ex.submit(_process_day, x) for x in arg_iter]
            for fut in as_completed(futs):
                d, out = fut.result(); done += 1
                for cname, rows in out.items():
                    for tmpl, w, nu in rows:
                        acc[cname]["stacks"].append(w); acc[cname]["patches"].append(tmpl)
                        acc[cname]["dates"].append(d); acc[cname]["ndet"].append(nu)
        del bydt, arg_iter
        print(f"  [{fi+1}/{len(files)}] {os.path.basename(f)}: {nrows:,} det, {len(days)} days "
              f"| cumulative Z rows {len(acc['Z']['stacks'])}", flush=True)
    for cname in CHANS:
        A = acc[cname]
        if not A["stacks"]:
            print(f"  [{cname}] no rows — skipped"); continue
        stacks = np.array(A["stacks"], np.float32); patches = np.array(A["patches"])
        dates = np.array(A["dates"]); ndet = np.array(A["ndet"], np.int32)
        order = np.lexsort((dates, patches))
        out = f"{a.out_prefix}_{cname}.npz"
        np.savez(out, stacks=stacks[order], patches=patches[order], dates=dates[order],
                 n_det=ndet[order], t=t_axis, fs=a.fs)
        print(f"[3/3] saved {out} ({stacks.shape})", flush=True)


if __name__ == "__main__":
    main()
