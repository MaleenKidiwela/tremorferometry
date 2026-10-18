#!/usr/bin/env python
"""E2 step 6b (added 2026-09-10, orchestrator request) — is the pilot GPU-bound or CPU-bound?

Observation that prompted this: with one queue left the L40S sat at 0-13 % utilisation and 860 MB
while `e2_pick.py` used 101 % of one core.  The pod is capped at 32 CPUs (cgroup cpu.max
3200000/100000; `nproc` reports 176 and lies), so the "GPU-hours" in the fleet projection are
really miniSEED read + resample + filter time on ONE CORE PER QUEUE.

Two measurements, both re-using `e2_pick.py`'s own code path (`e2_pick.picks_for_day`, i.e.
`model.annotate` + the exact threshold-sweep post-processing) without modifying it:

 1. `--profile N`  per-day wall-time split on one process:
       read      obspy.read of the day file
       pre       SeisBench's own `annotate_stream_pre` (1-8 Hz zerophase bandpass + decimation
                 to 20 Hz) -- pure CPU
       annotate  the full `model.annotate` (windowing + GPU forward + stitching) MINUS `pre`
       post      `picks_from_annotations` over the 7 sweep thresholds x 2 phases
    plus the same day annotated on the CPU, to bound what the GPU forward actually costs.
 2. `--workers 1,4,8,12`  concurrency scaling: N independent processes (spawn, one CUDA context
    each, `--torch-threads` each) on DISJOINT day sets -- no day file is read twice anywhere in the
    test, so the NFS page cache cannot flatter the later configurations.  nvidia-smi is sampled
    every 5 s for utilisation and memory.

Writes results/e2_profile.csv and results/e2_scaling.csv.
"""
from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import subprocess
import sys
import threading
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
RES = "lfe_catalog_repair/results"


def _worker(args):
    """process `paths` with the real pick code path; returns (n_days, seconds, n_picks)."""
    paths, model_name, torch_threads, batch, overlap, device, omp = args
    # omp=None reproduces `e2_run_pick.sh`, which sets no thread limit at all: OpenMP then defaults
    # to os.cpu_count() = 176 on this pod (the cgroup allows 32), so four queues each spawn 176
    # BLAS threads.  omp=<n> is the fix under test.
    if omp:
        os.environ["OMP_NUM_THREADS"] = str(omp)
        os.environ["OPENBLAS_NUM_THREADS"] = str(omp)
        os.environ["MKL_NUM_THREADS"] = str(omp)
    else:
        for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
            os.environ.pop(k, None)
    import torch
    torch.set_num_threads(torch_threads)
    import e2_pick
    from obspy import read
    model = e2_pick.load_model(model_name, device)
    # one warm-up day so the CUDA context / cuDNN autotune is not charged to the measurement
    st = read(paths[0])
    e2_pick.picks_for_day(model, st, batch, overlap)
    t0 = time.time()
    n_days = n_picks = 0
    for p in paths[1:]:
        try:
            st = read(p)
        except Exception:
            continue
        if len({tr.stats.channel for tr in st}) < 3:
            continue
        n_picks += len(e2_pick.picks_for_day(model, st, batch, overlap))
        n_days += 1
    return n_days, time.time() - t0, n_picks


class GpuSampler(threading.Thread):
    def __init__(self, period=5.0):
        super().__init__(daemon=True)
        self.period = period
        self.util, self.mem = [], []
        self._stop = threading.Event()

    def run(self):
        while not self._stop.is_set():
            try:
                out = subprocess.run(
                    ["nvidia-smi", "--query-gpu=utilization.gpu,memory.used",
                     "--format=csv,noheader,nounits"], capture_output=True, text=True, timeout=10)
                u, m = out.stdout.strip().split(",")
                self.util.append(float(u)); self.mem.append(float(m))
            except Exception:
                pass
            self._stop.wait(self.period)

    def stop(self):
        self._stop.set()


def profile(paths, model_name, batch, overlap, torch_threads):
    import torch
    torch.set_num_threads(torch_threads)
    import e2_pick
    import seisbench.models as sbm  # noqa: F401
    from obspy import read
    model = e2_pick.load_model(model_name, "cuda")
    cpu_model = None
    rows = []
    for i, p in enumerate(paths):
        t = {}
        t0 = time.time(); st = read(p); t["read"] = time.time() - t0
        argdict = model.default_args.copy()
        argdict.update(dict(batch_size=batch, overlap=overlap, stride=1))
        t0 = time.time()
        pre = model.annotate_stream_pre(st.copy(), argdict)
        torch.cuda.synchronize()
        t["pre"] = time.time() - t0
        t0 = time.time()
        ann = model.annotate(st, batch_size=batch, overlap=overlap)
        torch.cuda.synchronize()
        t["annotate_total"] = time.time() - t0
        t0 = time.time()
        n = 0
        for phase in ("P", "S"):
            sel = ann.select(channel=f"LFEDetect_{phase}")
            for thr in e2_pick.SWEEP:
                n += len(model.picks_from_annotations(sel, thr, phase))
        t["post"] = time.time() - t0
        t["gpu_and_window"] = t["annotate_total"] - t["pre"]
        if i < 2:                       # bound the pure-CPU forward on the first two days
            if cpu_model is None:
                cpu_model = e2_pick.load_model(model_name, "cpu")
            t0 = time.time()
            cpu_model.annotate(st, batch_size=batch, overlap=overlap)
            t["annotate_total_cpu"] = time.time() - t0
        t["day"] = os.path.basename(p)
        t["n_picks_sweep"] = n
        t["total"] = t["read"] + t["annotate_total"] + t["post"]
        rows.append(t)
        print(f"  {t['day']}: read {t['read']:.1f}s  pre {t['pre']:.1f}s  "
              f"gpu+window {t['gpu_and_window']:.1f}s  post {t['post']:.1f}s  "
              f"total {t['total']:.1f}s"
              + (f"  [cpu annotate {t.get('annotate_total_cpu', float('nan')):.1f}s]"
                 if "annotate_total_cpu" in t else ""), flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(f"{RES}/e2_profile.csv", index=False)
    med = df.median(numeric_only=True)
    print("\n[profile medians per station-day]")
    for k in ("read", "pre", "gpu_and_window", "post", "annotate_total", "total"):
        print(f"  {k:16s} {med[k]:7.2f} s   {100*med[k]/med['total']:5.1f} % of total")
    if "annotate_total_cpu" in df:
        print(f"  annotate on CPU  {df.annotate_total_cpu.median():7.2f} s "
              f"(vs {med['annotate_total']:.2f} s with the L40S)")
    print(f"[write] {RES}/e2_profile.csv")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--station", default="CN.PGC")
    ap.add_argument("--set", dest="setname", default="pgc2010")
    ap.add_argument("--model", default="cascadia")
    ap.add_argument("--batch-size", type=int, default=1024)
    ap.add_argument("--overlap", type=int, default=600)
    ap.add_argument("--torch-threads", type=int, default=2)
    ap.add_argument("--profile", type=int, default=0)
    ap.add_argument("--omp", type=int, default=2,
                    help="OMP/OPENBLAS/MKL thread limit per worker; 0 = leave unset "
                         "(reproduces e2_run_pick.sh)")
    ap.add_argument("--workers", default="")
    ap.add_argument("--days-per-worker", type=int, default=6)
    args = ap.parse_args()

    import e2_pick
    start, end = e2_pick.SETS[args.setname][args.station]
    allf = [p for _, p in e2_pick.day_files(f"{e2_pick.TRACES}/{args.station}", start, end)]
    print(f"[days] {len(allf)} day files available for {args.station} {start}..{end} "
          f"(set {args.setname})")
    cursor = 0
    if args.profile:
        print(f"[profile] {args.profile} days, 1 process, torch-threads {args.torch_threads}")
        profile(allf[:args.profile], args.model, args.batch_size, args.overlap, args.torch_threads)
        cursor = args.profile
    if args.workers:
        ctx = mp.get_context("spawn")
        rows = []
        for nw in [int(x) for x in args.workers.split(",")]:
            need = nw * (args.days_per_worker + 1)          # +1 warm-up day per worker
            if cursor + need > len(allf):
                print(f"[skip] {nw} workers needs {need} unread days, only "
                      f"{len(allf)-cursor} left")
                continue
            chunks = [allf[cursor + k * (args.days_per_worker + 1):
                           cursor + (k + 1) * (args.days_per_worker + 1)] for k in range(nw)]
            cursor += need
            smp = GpuSampler(5.0); smp.start()
            t0 = time.time()
            with ctx.Pool(nw) as pool:
                out = pool.map(_worker, [(c, args.model, args.torch_threads, args.batch_size,
                                          args.overlap, "cuda", args.omp) for c in chunks])
            wall = time.time() - t0
            smp.stop(); smp.join(timeout=6)
            days = sum(o[0] for o in out)
            picks = sum(o[2] for o in out)
            per_q = np.mean([o[0] / o[1] * 60 for o in out])
            rows.append(dict(station=args.station, workers=nw,
                             torch_threads=args.torch_threads, omp_limit=args.omp,
                             days_per_worker=args.days_per_worker, station_days=days,
                             wall_s=round(wall, 1),
                             aggregate_sta_days_per_min=round(days / wall * 60, 2),
                             per_queue_sta_days_per_min=round(per_q, 2),
                             steady_state_aggregate_per_min=round(nw * per_q, 1),
                             sta_days_per_min_per_cpu=round(days / wall * 60
                                                            / (nw * args.torch_threads), 3),
                             gpu_util_mean=round(float(np.mean(smp.util)), 1) if smp.util else None,
                             gpu_util_max=round(float(np.max(smp.util)), 1) if smp.util else None,
                             gpu_mem_max_MiB=round(float(np.max(smp.mem)), 0) if smp.mem else None,
                             n_gpu_samples=len(smp.util), picks=picks))
            print(f"[scaling] {nw:2d} workers x {args.torch_threads} threads: {days} station-days "
                  f"in {wall:.1f} s = {days/wall*60:.1f} sta-days/min aggregate "
                  f"({per_q:.1f} per queue), GPU util mean "
                  f"{rows[-1]['gpu_util_mean']} % max {rows[-1]['gpu_util_max']} %, "
                  f"mem max {rows[-1]['gpu_mem_max_MiB']} MiB", flush=True)
        sc = pd.DataFrame(rows)
        out = f"{RES}/e2_scaling.csv"
        sc.to_csv(out, mode="a", header=not os.path.exists(out), index=False)
        print(sc.to_string(index=False))
        print(f"[write] {RES}/e2_scaling.csv")


if __name__ == "__main__":
    main()
