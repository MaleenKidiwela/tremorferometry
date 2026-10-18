#!/usr/bin/env python
"""P3D step 3 - resumable queue runner (<=8 concurrent FDSN requests, total).

Walks results/download_queue.csv in priority order (P0 -> P3, oldest day first, all
stations of a day together), appends one row per finished station-day to
results/download_manifest.csv (atomic: single os.write under a lock to an O_APPEND fd, so
a reader never sees a torn line), refreshes results/download_progress.json every 5 min,
and calls scripts/download_qc.py when every queued day of a calendar month is finished.

Resume: any (sta, date) already carrying a terminal status in the manifest is skipped, and
so is any day whose .mseed / .mseed.nd file is already on disk - restart never re-fetches.

Usage:  python scripts/download_worker.py [--workers 8] [--limit N] [--priorities P0,P1]
"""
import os
import csv
import sys
import json
import time
import signal
import argparse
import threading
import subprocess
from datetime import date, datetime, timezone
from collections import defaultdict, deque
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE, "scripts"))
from r_download_day import fetch_day                                    # noqa: E402

RESULTS = os.path.join(HERE, "results", "p3d")
TRACES = os.path.join(HERE, "traces")
QUEUE = os.path.join(RESULTS, "download_queue.csv")
MANIFEST = os.path.join(RESULTS, "download_manifest.csv")
PROGRESS = os.path.join(RESULTS, "download_progress.json")
QC = os.path.join(HERE, "scripts", "r_download_qc.py")
COLS = ["net", "sta", "band", "date", "doy", "status", "bytes", "n_channels", "fs",
        "coverage_frac", "n_frag", "updated_utc", "note"]
MAX_PER_PROVIDER = 8
PROGRESS_EVERY = 300.0
TERMINAL = {"data", "nodata", "no_epoch", "present_elsewhere"}

_stop = threading.Event()
_mlock = threading.Lock()
_mfd = None
_rate = deque()                      # (t, bytes) for the MB/s window
_counts = defaultdict(int)           # status -> n (this run)
_bytes_total = 0
_per_sy = defaultdict(int)           # "STA YYYY" -> n days with data
_month_left = defaultdict(int)
_month_done_fired = set()
_current = {"month": "", "priority": ""}
_t0 = time.time()


def log(msg):
    print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] {msg}", flush=True)


def manifest_open():
    global _mfd
    if not os.path.exists(MANIFEST):
        with open(MANIFEST, "w", newline="") as f:
            csv.writer(f).writerow(COLS)
    _mfd = os.open(MANIFEST, os.O_WRONLY | os.O_APPEND | os.O_CREAT, 0o644)


def manifest_append(rec):
    """One complete line, one os.write, under a lock -> readers never see a torn row."""
    import io
    buf = io.StringIO()
    csv.writer(buf, lineterminator="\n").writerow([rec.get(c, "") for c in COLS])
    data = buf.getvalue().encode()
    with _mlock:
        os.write(_mfd, data)


def load_done():
    done = {}
    if os.path.exists(MANIFEST):
        with open(MANIFEST) as f:
            for r in csv.DictReader(f):
                if r.get("status") in TERMINAL:
                    done[(r["sta"], r["date"])] = r["status"]
                elif (r["sta"], r["date"]) not in done:
                    done[(r["sta"], r["date"])] = r.get("status", "")
    return done


def load_queue(prios):
    rows = []
    with open(QUEUE) as f:
        for r in csv.DictReader(f):
            if r["status"] != "pending" or r["priority"] not in prios:
                continue
            rows.append(r)
    order = {"P0": 0, "P1": 1, "P2": 2}
    rows.sort(key=lambda r: (order.get(r["priority"], 9), r["date"], r["sta"]))
    return rows


def write_progress(total, done_now):
    now = time.time()
    with _mlock:
        pts = [(t, b) for t, b in _rate if now - t <= 1800]
    mbps = (sum(b for _, b in pts) / 1e6 / max(1e-9, (now - pts[0][0]))) if len(pts) > 1 else 0.0
    dph = (len(pts) / max(1e-9, (now - pts[0][0])) * 3600) if len(pts) > 1 else 0.0
    left = total - done_now
    eta_h = left / dph if dph > 0 else None
    out = dict(
        updated_utc=datetime.now(timezone.utc).isoformat(timespec="seconds"),
        elapsed_h=round((now - _t0) / 3600, 3),
        queue_total_this_run=total, finished_this_run=done_now, remaining=left,
        status_counts=dict(_counts), gb_done_this_run=round(_bytes_total / 1e9, 3),
        mb_per_s_last30min=round(mbps, 3), station_days_per_hour_last30min=round(dph, 1),
        eta_hours_remaining=(round(eta_h, 1) if eta_h else None),
        eta_finish_utc=(datetime.fromtimestamp(now + eta_h * 3600, timezone.utc)
                        .isoformat(timespec="seconds") if eta_h else None),
        current_month=_current["month"], current_priority=_current["priority"],
        months_complete=sorted(_month_done_fired),
        station_year_days=dict(sorted(_per_sy.items())),
    )
    tmp = PROGRESS + ".tmp"
    with open(tmp, "w") as f:
        json.dump(out, f, indent=1)
    os.replace(tmp, PROGRESS)
    return out


def run_qc(month):
    try:
        subprocess.run([sys.executable, QC, "--month", month], check=False,
                       stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, timeout=1800)
        log(f"QC written for {month}")
    except Exception as e:
        log(f"QC failed for {month}: {e}")


def main():
    global _bytes_total
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=MAX_PER_PROVIDER)
    ap.add_argument("--priorities", default="P0,P1,P2")
    ap.add_argument("--limit", type=int, default=0, help="stop after N station-days (test)")
    ap.add_argument("--no-qc", action="store_true", default=True)
    ap.add_argument("--qc", dest="no_qc", action="store_false", help="run the per-month QC hook")
    a = ap.parse_args()
    workers = min(a.workers, MAX_PER_PROVIDER)
    prios = set(a.priorities.split(","))

    signal.signal(signal.SIGTERM, lambda *_: _stop.set())
    signal.signal(signal.SIGINT, lambda *_: _stop.set())

    manifest_open()
    done = load_done()
    rows = load_queue(prios)
    # every queued month, so a month that is already fully on disk still fires its QC
    for r in rows:
        _month_left[r["date"][:7]] += 1
    todo = [r for r in rows if done.get((r["sta"], r["date"])) not in TERMINAL]
    already = len(rows) - len(todo)
    for r in rows:
        if done.get((r["sta"], r["date"])) in TERMINAL:
            _month_left[r["date"][:7]] -= 1
    for m, left in sorted(_month_left.items()):
        if left <= 0 and not os.path.exists(os.path.join(RESULTS, f"download_qc_{m}.csv")):
            _month_done_fired.add(m)
            if not a.no_qc:
                run_qc(m)                      # month already complete on disk (resume)
    if a.limit:
        todo = todo[:a.limit]
    log(f"queue {len(rows)} pending rows in {sorted(prios)}; {already} already done; "
        f"{len(todo)} to fetch; workers={workers}")

    total, done_now = len(todo), 0
    last_prog = 0.0
    lock = threading.Lock()

    attempted = set()

    def work(r):
        if _stop.is_set():
            return None
        attempted.add((r["sta"], r["date"]))
        rec = fetch_day(r["net"], r["sta"], r["band"], date.fromisoformat(r["date"]),
                        r.get("loc", ""), TRACES)
        rec["updated_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
        manifest_append(rec)
        return rec

    def account(rec, r):
        nonlocal done_now, last_prog
        global _bytes_total
        with lock:
            done_now += 1
            _counts[rec["status"]] += 1
            nb = int(rec["bytes"] or 0)
            _bytes_total += nb
            _rate.append((time.time(), nb))
            while _rate and time.time() - _rate[0][0] > 1800:
                _rate.popleft()
            if rec["status"] == "data":
                _per_sy[f"{rec['sta']} {rec['date'][:4]}"] += 1
            _current["month"], _current["priority"] = r["date"][:7], r["priority"]
            m = r["date"][:7]
            _month_left[m] -= 1
            fire = (_month_left[m] <= 0 and m not in _month_done_fired)
            if fire:
                _month_done_fired.add(m)
            if time.time() - last_prog > PROGRESS_EVERY or done_now == total:
                last_prog = time.time()
                p = write_progress(total, done_now)
                log(f"{done_now}/{total} {dict(_counts)} {p['gb_done_this_run']} GB "
                    f"{p['mb_per_s_last30min']} MB/s eta={p['eta_hours_remaining']} h "
                    f"month={p['current_month']}")
        if fire and not a.no_qc:
            run_qc(m)

    def drive(items, tag):
        nonlocal done_now
        with ThreadPoolExecutor(max_workers=workers) as ex:
            it = iter(items)
            fut = {}
            for r in [next(it, None) for _ in range(workers * 2)]:
                if r is not None:
                    fut[ex.submit(work, r)] = r
            while fut:
                for f in list(fut):
                    if not f.done():
                        continue
                    r = fut.pop(f)
                    try:
                        rec = f.result()
                    except Exception as e:                      # should not happen
                        rec = dict(net=r["net"], sta=r["sta"], band=r["band"],
                                   date=r["date"], doy=r["doy"], status="error",
                                   bytes="", n_channels="", fs="", coverage_frac="",
                                   n_frag="", note=f"worker:{type(e).__name__}",
                                   updated_utc=datetime.now(timezone.utc)
                                   .isoformat(timespec="seconds"))
                        manifest_append(rec)
                    if rec is not None:
                        account(rec, r)
                    if not _stop.is_set():
                        nxt = next(it, None)
                        if nxt is not None:
                            fut[ex.submit(work, nxt)] = nxt
                if fut:
                    time.sleep(0.2)
        log(f"{tag} finished: {dict(_counts)}")

    drive(todo, "pass 1")

    # second pass: retry everything that errored, once
    if not _stop.is_set():
        errs = []
        seen = load_done()
        for r in rows:
            if ((r["sta"], r["date"]) in attempted
                    and seen.get((r["sta"], r["date"])) not in TERMINAL):
                errs.append(r)
        if errs:
            log(f"pass 2: retrying {len(errs)} error rows once")
            total += len(errs)
            for r in errs:
                _month_left[r["date"][:7]] += 1
            drive(errs, "pass 2")

    write_progress(total, done_now)
    log("ALL DONE" if not _stop.is_set() else "STOPPED (resumable)")


if __name__ == "__main__":
    main()
