#!/usr/bin/env python
"""P4-0 step 2 - fetch the 8 probe days for every candidate (and the reference pair).

Uses the P3D downloader VERBATIM (`r_download_day.fetch_day`: retry-wrapped FDSN, `.nd` no-data
markers, native sample rate, nothing filtered or merged) behind a <= 8-thread pool, so IRIS never
sees more than MAX_CONCURRENT_FDSN concurrent requests.  Resumable: a day already on disk is only
re-verified.

Before any request, every queued station-day that already exists in the READ-ONLY
../lfe_extension/traces tree is SYMLINKED into traces/ (never copied, never modified there) as long
as its channels match the band we want - this is where CN.PGC and PB.B011 come from.

Writes results/p4/download_manifest.csv (one row per station-day, appended atomically).
"""
from __future__ import annotations
import os
import sys
import csv
import time
import argparse
import threading
from datetime import date, datetime, timezone
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import u_common as U                                                   # noqa: E402
from r_download_day import fetch_day, day_path, verify                 # noqa: E402

QUEUE = os.path.join(U.P4, "probe_queue.csv")
MANIFEST = os.path.join(U.P4, "download_manifest.csv")
COLS = ["net", "sta", "band", "loc", "date", "doy", "role", "src", "status", "bytes",
        "n_channels", "fs", "coverage_frac", "n_frag", "updated_utc", "note"]
_lock = threading.Lock()
_n = {"done": 0, "data": 0, "nodata": 0, "error": 0, "link": 0, "bytes": 0}
_t0 = time.time()


def log(m):
    print(f"[{datetime.now(timezone.utc).strftime('%H:%M:%S')}] {m}", flush=True)


def link_from_extension(net, sta, band, d):
    """Symlink an existing ../lfe_extension/traces day into traces/ if its band matches."""
    ns = f"{net}.{sta}"
    doy = d.timetuple().tm_yday
    dst = U.trace_path(ns, d.year, doy)
    if os.path.exists(dst) or os.path.exists(dst + ".nd"):
        return "present"
    src = os.path.join(U.EXT_TRACES, ns, str(d.year), f"{doy:03d}.mseed")
    if not os.path.exists(src):
        return ""
    try:
        from obspy import read
        chans = {tr.stats.channel for tr in read(src, headonly=True)}
    except Exception:
        return ""
    if not any(c.startswith(band) for c in chans):
        return ""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    try:
        os.symlink(os.path.realpath(src), dst)
    except FileExistsError:
        pass
    return "linked"


def one(row):
    net, sta, band = row["net"], row["sta"], row["band"]
    loc = (row.get("loc") or "").strip()
    d = date.fromisoformat(row["date"])
    src = ""
    try:
        src = link_from_extension(net, sta, band, d)
    except Exception as e:
        src = f"link_err:{type(e).__name__}"
    rec = fetch_day(net, sta, band, d, loc, U.TRACES, prov="IRIS", tries=3, timeout=180)
    rec["loc"] = loc
    rec["role"] = row.get("role", "")
    rec["src"] = src
    rec["updated_utc"] = datetime.now(timezone.utc).isoformat(timespec="seconds")
    with _lock:
        _n["done"] += 1
        _n[rec["status"]] = _n.get(rec["status"], 0) + 1
        if src == "linked":
            _n["link"] += 1
        try:
            _n["bytes"] += int(rec["bytes"] or 0)
        except (TypeError, ValueError):
            pass
        with open(MANIFEST, "a", newline="") as f:
            csv.DictWriter(f, COLS, extrasaction="ignore").writerow(rec)
        if _n["done"] % 25 == 0:
            el = time.time() - _t0
            log(f"  {_n['done']} done  data={_n.get('data',0)} nodata={_n.get('nodata',0)} "
                f"err={_n.get('error',0)} linked={_n['link']}  {_n['bytes']/2**30:.1f} GiB  "
                f"{el/60:.1f} min  ({_n['done']/max(el,1)*3600:.0f} day/h)")
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=U.MAX_CONCURRENT_FDSN)
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    with open(QUEUE) as f:
        rows = list(csv.DictReader(f))
    done = set()
    if os.path.exists(MANIFEST):
        with open(MANIFEST) as f:
            for r in csv.DictReader(f):
                if r.get("status") in ("data", "nodata"):
                    done.add((r["net"], r["sta"], r["date"]))
    else:
        with open(MANIFEST, "w", newline="") as f:
            csv.DictWriter(f, COLS).writeheader()
    todo = [r for r in rows if (r["net"], r["sta"], r["date"]) not in done]
    if a.limit:
        todo = todo[:a.limit]
    log(f"[download] {len(todo)} station-days to fetch ({len(done)} already terminal), "
        f"{a.workers} concurrent FDSN requests")
    with ThreadPoolExecutor(max_workers=min(a.workers, U.MAX_CONCURRENT_FDSN)) as ex:
        list(ex.map(one, todo))
    log(f"[download] finished: {_n} in {(time.time()-_t0)/60:.1f} min")


if __name__ == "__main__":
    main()
