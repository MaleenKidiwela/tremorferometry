#!/usr/bin/env python
"""E2 step 1 — run the SeisBench LFEDetect picker over the waveforms already on disk.

Contract (from `results/external/other/seisbench_lfedetect/*.json`, Muenchmeyer et al. 2024):
    3-component ZNE, 20 Hz, 1200 samples = 60 s windows, 1-8 Hz zerophase bandpass,
    norm="peak" over (components, time), sigmoid output, classes P and S, default overlap 600.
SeisBench's own `annotate` does the whole pre-processing chain (bandpass at native rate, then
lowpass+decimate to 20 Hz -- both integer factors here, 40/20 and 100/20), so no hand-rolled
pre-processing is used and the training-time operation order is preserved.

THRESHOLD SWEEP.  `picks_from_annotations` triggers on prob > t and off at t/2, so the trigger
intervals for a higher threshold are strictly nested inside those of a lower one: the pick set at
threshold t is exactly {picks with thr_first <= t <= prob}.  We therefore pick at every sweep
threshold, take the union keyed on (phase, peak sample), and store `thr_first` + `prob`.  Any sweep
threshold can then be reconstructed exactly from the stored table without re-annotating.

Usage
  python lfe_catalog_repair/scripts/e2_pick.py --set pgc2010 --model cascadia \
      --stations CN.PGC --out lfe_catalog_repair/results/e2_picks_pgc2010_cascadia.parquet
Environment (GPU): export LD_LIBRARY_PATH=$(ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib | tr '\n' ':')$LD_LIBRARY_PATH
                   export PYTHONPATH=/opt/conda/lib/python3.13/site-packages   # cu128 torch
"""
from __future__ import annotations

import argparse
import glob
import os
import time
from datetime import date, datetime, timedelta

import numpy as np
import pandas as pd

TRACES = "lfe_catalog_repair/traces"
RESULTS = "lfe_catalog_repair/results"
MODEL_CACHE_NAMES = ("cascadia", "jcms", "cms", "jcs", "jms")
SWEEP = (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)

# station -> (start, end) per named set; end exclusive
SETS = {
    # (a) PGC over the Lin/Bostock ground-truth window
    "pgc2010": {"CN.PGC": ("2010-05-01", "2010-11-01")},
    # (b) southern Vancouver Island network, 2017 ETS (3-component stations only)
    "vi2017": {s: ("2017-02-20", "2017-04-21") for s in
               ["CN.PGC", "CN.SYMB", "CN.SNB", "CN.VGZ", "CN.CLRS", "CN.PTRF", "PB.B011", "UW.LRIV"]},
    # (b2) the 2010 network for association scoring against Lin/Bostock (34 densest Lin days)
    "vi2010": {s: ("2010-08-10", "2010-09-14") for s in
               ["CN.PGC", "CN.SNB", "CN.VGZ", "CN.LZB", "CN.PFB", "CN.YOUB", "CN.NLLB",
                "PB.B011", "PB.B926", "UW.LRIV"]},
    # (c) off-ground-truth stations, two full years
    "offgt": {s: ("2017-01-01", "2019-01-01") for s in ["UW.WISH", "PB.B943", "UW.YACT"]},
}


def day_files(sta_dir: str, start: str, end: str):
    d0 = datetime.strptime(start, "%Y-%m-%d").date()
    d1 = datetime.strptime(end, "%Y-%m-%d").date()
    out = []
    d = d0
    while d < d1:
        p = f"{sta_dir}/{d.year}/{d.timetuple().tm_yday:03d}.mseed"
        if os.path.exists(p):
            out.append((d, p))
        d += timedelta(days=1)
    return out


def load_model(name: str, device: str):
    import seisbench.models as sbm
    m = sbm.LFEDetect.from_pretrained(name, version_str="1")
    m.to(device)
    m.eval()
    return m


def picks_for_day(model, st, batch_size: int, overlap: int):
    """Annotate one day-stream and return the exact threshold-sweep union of picks."""
    ann = model.annotate(st, batch_size=batch_size, overlap=overlap)
    rows = {}
    for phase in ("P", "S"):
        sel = ann.select(channel=f"LFEDetect_{phase}")
        if len(sel) == 0:
            continue
        for thr in SWEEP:                       # ascending: thr_first set on first appearance
            for pk in model.picks_from_annotations(sel, thr, phase):
                key = (phase, round(float(pk.peak_time.timestamp) * 20))
                if key in rows:
                    continue
                rows[key] = dict(
                    phase=phase,
                    time=pk.peak_time.datetime,
                    prob=float(pk.peak_value),
                    thr_first=thr,
                    t_start=pk.start_time.datetime,
                    t_end=pk.end_time.datetime,
                )
    return list(rows.values())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", dest="setname", required=True)
    ap.add_argument("--model", default="cascadia", choices=MODEL_CACHE_NAMES)
    ap.add_argument("--stations", default="", help="comma list to restrict the set")
    ap.add_argument("--out", default="")
    ap.add_argument("--batch-size", type=int, default=1024)
    ap.add_argument("--overlap", type=int, default=600)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--torch-threads", type=int, default=2)
    ap.add_argument("--allow-incomplete", action="store_true",
                    help="also run station-days with fewer than 3 components (zero-filled)")
    args = ap.parse_args()

    import torch
    torch.set_num_threads(args.torch_threads)
    from obspy import read

    dev = args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu"
    model = load_model(args.model, dev)
    gpu = torch.cuda.get_device_name(0) if dev == "cuda" else "cpu"

    conf = SETS[args.setname]
    stations = [s for s in conf] if not args.stations else [s for s in args.stations.split(",")]
    out = args.out or f"{RESULTS}/e2_picks_{args.setname}_{args.model}.parquet"

    frames, tp_rows = [], []
    for sta in stations:
        start, end = conf[sta]
        net, code = sta.split(".")
        files = day_files(f"{TRACES}/{sta}", start, end)
        t0 = time.time()
        n_days = n_skip = 0
        rows = []
        for d, p in files:
            try:
                st = read(p)
            except Exception as e:
                print(f"[warn] {sta} {d} read failed: {e}", flush=True)
                n_skip += 1
                continue
            ncomp = len({tr.stats.channel for tr in st})
            if ncomp < 3 and not args.allow_incomplete:
                n_skip += 1
                continue
            try:
                pk = picks_for_day(model, st, args.batch_size, args.overlap)
            except Exception as e:
                print(f"[warn] {sta} {d} annotate failed: {e}", flush=True)
                n_skip += 1
                continue
            for r in pk:
                r["net"], r["sta"], r["day"] = net, code, d
            rows.extend(pk)
            n_days += 1
            if n_days % 50 == 0:
                el = time.time() - t0
                print(f"  {sta} {n_days}/{len(files)} days  {n_days/el*60:.1f} sta-days/min  "
                      f"{len(rows)} picks", flush=True)
        el = max(time.time() - t0, 1e-6)
        df = pd.DataFrame(rows)
        if len(df):
            df["model"] = args.model
            df["set"] = args.setname
            frames.append(df)
        tp_rows.append(dict(set=args.setname, model=args.model, station=sta, n_days=n_days,
                            n_skipped=n_skip, n_picks=len(df), seconds=round(el, 1),
                            station_days_per_min=round(n_days / el * 60, 2), device=dev, gpu=gpu))
        print(f"[done] {sta} {n_days} days ({n_skip} skipped) {len(df)} picks "
              f"{n_days/el*60:.1f} sta-days/min", flush=True)

    if frames:
        allp = pd.concat(frames, ignore_index=True)
        for c in ("time", "t_start", "t_end"):
            allp[c] = pd.to_datetime(allp[c]).astype("datetime64[ns]")
        allp = allp[["set", "model", "net", "sta", "day", "phase", "time", "prob",
                     "thr_first", "t_start", "t_end"]]
        os.makedirs(os.path.dirname(out), exist_ok=True)
        allp.to_parquet(out, index=False)
        print(f"[write] {out}  {len(allp)} rows", flush=True)

    tp = pd.DataFrame(tp_rows)
    tpf = f"{RESULTS}/e2_throughput.csv"
    tp.to_csv(tpf, mode="a", header=not os.path.exists(tpf), index=False)
    print(f"[write] {tpf} (+{len(tp)} rows)", flush=True)


if __name__ == "__main__":
    main()
