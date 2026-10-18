#!/usr/bin/env python
"""Track C step 5 — the resumable production watcher.

Loop:  oldest RUNNABLE, unprocessed month  ->  cache  ->  matched filter  ->  QC  ->  mark done.
Nothing runnable (the download has not reached it yet) -> sleep --sleep-s (default 600 s) and look again.
Resumable: state is results/trackC/months_done.csv plus the month's own output files; a month whose
outputs already exist is not recomputed unless --force.  Single writer: an flock on
results/trackC/.watch.lock — a second watcher exits immediately.

  nohup bash scripts/c_run_all.sh > /dev/null 2>&1 &     # the normal way to start it
  python scripts/c_watch.py --once                       # one month, then stop
  python scripts/c_watch.py --months 2015-12,2018-06     # a named list, in order
  kill $(cat logs/trackC_watch.pid)                      # stop; safe at any point
"""
from __future__ import annotations
import os, sys, time, json, fcntl, argparse, subprocess, datetime
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C
import c_month_ready as MR

PY = "/opt/conda/bin/python"
HERE = os.path.dirname(os.path.abspath(__file__))
LOGS = f"{C.EXT}/logs"


def gpu_env():
    env = dict(os.environ)
    libs = subprocess.run("ls -d /opt/conda/lib/python3.13/site-packages/nvidia/*/lib",
                          shell=True, capture_output=True, text=True).stdout.split()
    env["LD_LIBRARY_PATH"] = ":".join(libs + [env.get("LD_LIBRARY_PATH", "")])
    env["PYTHONPATH"] = "/opt/conda/lib/python3.13/site-packages"
    env["CUDA_PATH"] = "/opt/conda/targets/x86_64-linux"
    env["OMP_NUM_THREADS"] = "4"
    env["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
    return env


def run(cmd, log, env=None):
    with open(log, "a") as fh:
        fh.write(f"\n=== {' '.join(cmd)}  {datetime.datetime.now(datetime.timezone.utc).isoformat()}Z ===\n")
        fh.flush()
        return subprocess.call(cmd, stdout=fh, stderr=fh, env=env or os.environ, cwd=C.EXT)


def _derived_epochs():
    p = f"{C.TRC}/tiers_by_epoch.csv"
    if not os.path.exists(p):
        return set()
    T = pd.read_csv(p)
    return set(T[T.status.astype(str).str.startswith("derived")].epoch.astype(str))


def done_table():
    p = f"{C.TRC}/months_done.csv"
    return pd.read_csv(p) if os.path.exists(p) else pd.DataFrame(columns=["ym"])


def mark_done(ym, meta):
    p = f"{C.TRC}/months_done.csv"
    row = dict(ym=ym, finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
               n_detections=meta.get("n_detections"), n_det_ext=meta.get("n_det_ext"),
               n_det_core=meta.get("n_det_core"), n_det_core_hfpass=meta.get("n_det_core_hfpass"),
               n_family_days=meta.get("n_family_days"),
               n_family_days_unobservable=meta.get("n_family_days_unobservable"),
               template_set=meta.get("template_set"),
               epoch_template=("|".join(f"{k}:{v}" for k, v in (meta.get("epoch_template_counts") or {}).items())),
               n_cat_ext=meta.get("n_cat_ext"), n_cat_core=meta.get("n_cat_core"),
               runtime_s=meta.get("runtime_s"))
    D = done_table()
    D = D[D.ym != ym] if len(D) else D
    D = pd.concat([D, pd.DataFrame([row])], ignore_index=True).sort_values("ym")
    D.to_csv(p, index=False)


def do_month(ym, a, env):
    t0 = time.time()
    print(f"[watch] {ym} START {datetime.datetime.now(datetime.timezone.utc).isoformat()}Z", flush=True)
    if a.force or not os.path.exists(f"{C.TRC}/daily_noise_{ym}.parquet"):
        rc = run([PY, f"{HERE}/c_cache_month.py", ym, "--workers", str(a.cache_workers)]
                 + (["--force"] if a.force else []), f"{LOGS}/trackC_cache_{ym}.log")
        if rc != 0:
            print(f"[watch] {ym} CACHE FAILED rc={rc}", flush=True)
            return False
    if a.force or not os.path.exists(f"{C.TRC}/detections_{ym}.parquet"):
        for attempt in range(a.retries + 1):
            rc = run([PY, f"{HERE}/c_run_month.py", ym], f"{LOGS}/trackC_run_{ym}.log", env)
            if rc == 0:
                break
            print(f"[watch] {ym} RUN rc={rc} (attempt {attempt+1}/{a.retries+1})", flush=True)
            time.sleep(120)
        else:
            return False
    rc = run([PY, f"{HERE}/c_qc_month.py", ym], f"{LOGS}/trackC_qc_{ym}.log")
    if rc != 0:
        print(f"[watch] {ym} QC FAILED rc={rc}", flush=True)
    # Re-derive the per-network-epoch tiers.  Rows already marked `derived` are FROZEN, so this can
    # only ADD an epoch (when its quiet months finally arrive); when one appears, the whole record is
    # re-flagged once so an epoch never carries two different thresholds.
    derived_before = _derived_epochs()
    run([PY, f"{HERE}/c_tiers_by_epoch.py"], f"{LOGS}/trackC_tiers.log")
    new = _derived_epochs() - derived_before
    if new:
        print(f"[watch] network epoch(s) {sorted(new)} became derivable — re-flagging the record",
              flush=True)
        run([PY, f"{HERE}/c_reflag_tiers.py", "--all"], f"{LOGS}/trackC_tiers.log")
    meta = {}
    mp = f"{C.TRC}/run_{ym}.json"
    if os.path.exists(mp):
        meta = json.load(open(mp))
    mark_done(ym, meta)
    print(f"[watch] {ym} DONE in {time.time()-t0:.0f}s: {meta.get('n_det_ext')} T_ext / "
          f"{meta.get('n_det_core')} T_core detections", flush=True)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sleep-s", type=float, default=600.0)
    ap.add_argument("--cache-workers", type=int, default=12)
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--retries", type=int, default=2)
    ap.add_argument("--once", action="store_true")
    ap.add_argument("--months", default="", help="explicit comma-separated list, in order")
    ap.add_argument("--max-months", type=int, default=0)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--from-ym", default="", help="ignore months before this one (e.g. 2014-01)")
    ap.add_argument("--to-ym", default="", help="ignore months after this one")
    ap.add_argument("--allow-partial", action="store_true",
                    help="run a month even when some station-days are still unfinished (dry runs only)")
    a = ap.parse_args()
    os.makedirs(C.TRC, exist_ok=True); os.makedirs(LOGS, exist_ok=True)
    lock = open(f"{C.TRC}/.watch.lock", "w")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        print("[watch] another watcher holds the lock — exiting")
        return
    lock.write(f"{os.getpid()} {datetime.datetime.now(datetime.timezone.utc).isoformat()}Z\n"); lock.flush()
    open(f"{LOGS}/trackC_watch.pid", "w").write(str(os.getpid()) + "\n")
    env = gpu_env()
    n = 0
    todo = [m for m in a.months.split(",") if m] if a.months else None
    stop = f"{C.TRC}/STOP"
    while True:
        # GRACEFUL STOP: `touch results/trackC/STOP` makes the watcher finish the month it is on and
        # exit before starting another.  Remove the file before restarting.
        if os.path.exists(stop):
            print(f"[watch] {stop} present — stopping after the finished month", flush=True)
            try:
                os.remove(stop)
            except OSError:
                pass
            return
        done = set(done_table().ym.astype(str)) if not a.force else set()
        if todo is not None:
            pend = [m for m in todo if m not in done]
            ym = pend[0] if pend else None
        else:
            man = C.manifest(); ep = C.epochs()
            ym = None
            for cand in sorted(man.date.str[:7].unique()):
                if cand in done or (a.from_ym and cand < a.from_ym) or (a.to_ym and cand > a.to_ym):
                    continue
                st = MR.month_state(cand, man, ep, a.min_cov)
                if st["runnable"] or a.allow_partial:
                    if st["n_station_days"] == 0:
                        continue
                    ym = cand
                    break
        if ym is None:
            if a.once or todo is not None:
                print("[watch] nothing left to do", flush=True)
                return
            print(f"[watch] nothing runnable; sleeping {a.sleep_s:.0f}s", flush=True)
            time.sleep(a.sleep_s)
            continue
        do_month(ym, a, env)
        n += 1
        if a.once or (a.max_months and n >= a.max_months):
            print("[watch] stopping (--once / --max-months)", flush=True)
            return


if __name__ == "__main__":
    main()
