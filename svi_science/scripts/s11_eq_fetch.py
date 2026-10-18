#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6.2): fetch ComCat phase-data QuakeML for the selected events.

  python s11_eq_fetch.py [--workers 2] [--tier 2]
Cache: data/eq_picks/raw/comcat/<cc_id>.geojson, <cc_id>.phase.xml, markers <cc_id>.done / .nophase / .fail.
Cached events are skipped; .fail events are retried on the next run. Writes data/eq_picks/fetch_log.csv (append).
"""
import argparse
import json
import os
import random
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
SEL = f"{ROOT}/svi_science/results/s11_depth/tomo/eq/selection.csv"
RAW = f"{ROOT}/data/eq_picks/raw/comcat"
LOG = f"{ROOT}/data/eq_picks/fetch_log.csv"
DETAIL = "https://earthquake.usgs.gov/fdsnws/event/1/query?eventid={}&format=geojson"
SUMQ = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=geojson&starttime={}-01-01&endtime={}-01-01"
        "&minlatitude=45.9&maxlatitude=50.5&minlongitude=-128&maxlongitude=-120&minmagnitude=0.9&limit=20000")
PHURL = "https://earthquake.usgs.gov/product/phase-data/{}/uw/{}/quakeml.xml"
SUMMARY = {}       # id -> (updated ms, has phase-data)
UA = "tremorferometry-s11-eq/1.0 (research; python-urllib)"
BACKOFF = (4, 8, 16, 32, 64, 128)
MIN_GAP = 0.5
lock = threading.Lock()
stats = dict(done=0, nophase=0, fail=0, skip=0, requests=0, bytes=0)


def get(url):
    """GET with retry and backoff; returns bytes, or raises after the last try. 404/204 -> None."""
    last = None
    for k in range(len(BACKOFF) + 1):
        time.sleep(MIN_GAP + random.random() * 0.2)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": UA})
            with urllib.request.urlopen(req, timeout=90) as r:
                if r.status == 204:
                    return None
                b = r.read()
            with lock:
                stats["requests"] += 1; stats["bytes"] += len(b)
            return b
        except urllib.error.HTTPError as e:
            with lock:
                stats["requests"] += 1
            if e.code in (404, 204, 403):
                return None
            last = f"HTTP {e.code}"
        except Exception as e:  # timeouts, resets
            last = f"{type(e).__name__}: {e}"
        if k < len(BACKOFF):
            time.sleep(BACKOFF[k] * (1.5 if last == "HTTP 429" else 1.0))
    raise RuntimeError(last)


def wr(path, b):
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        f.write(b)
    os.replace(tmp, path)


def one(cc):
    base = f"{RAW}/{cc}"
    if os.path.exists(base + ".done") or os.path.exists(base + ".nophase"):
        with lock:
            stats["skip"] += 1
        return cc, "skip", ""
    try:
        # fast path: a uw event whose last update is its uw product -> phase-data URL is known without the detail call
        sm = SUMMARY.get(cc)
        if cc.startswith("uw") and sm and sm[1] and not os.path.exists(base + ".geojson"):
            url = PHURL.format(cc, sm[0])
            try:
                x = get(url)
            except RuntimeError:
                x = None
            if x is not None and (b"quakeml" in x[:2000] or b"<pick" in x[:200000]):
                wr(base + ".phase.xml", x)
                wr(base + ".done", f"uw {url} (summary path)\n".encode())
                if os.path.exists(base + ".fail"):
                    os.remove(base + ".fail")
                with lock:
                    stats["done"] += 1
                return cc, "done", "uw-summary"
        if os.path.exists(base + ".geojson"):
            gj = open(base + ".geojson", "rb").read()
        else:
            gj = get(DETAIL.format(cc))
            if gj is None:
                wr(base + ".nophase", b"no event (404/204)\n")
                with lock:
                    stats["nophase"] += 1
                return cc, "noevent", ""
            wr(base + ".geojson", gj)
        props = json.loads(gj)["properties"]
        ph = props.get("products", {}).get("phase-data", [])
        ph = [p for p in ph if "quakeml.xml" in p.get("contents", {})]
        if not ph:
            wr(base + ".nophase", b"no phase-data product\n")
            with lock:
                stats["nophase"] += 1
            return cc, "nophase", ""
        ph.sort(key=lambda p: (p.get("source") != "uw", -int(p.get("preferredWeight", 0))))
        url = ph[0]["contents"]["quakeml.xml"]["url"]
        x = get(url)
        if x is None:
            raise RuntimeError("phase-data 404")
        wr(base + ".phase.xml", x)
        wr(base + ".done", f"{ph[0].get('source')} {url}\n".encode())
        if os.path.exists(base + ".fail"):
            os.remove(base + ".fail")          # our own marker, superseded by success
        with lock:
            stats["done"] += 1
        return cc, "done", ph[0].get("source")
    except Exception as e:
        wr(base + ".fail", f"{time.time():.0f} {e}\n".encode())
        with lock:
            stats["fail"] += 1
        return cc, "fail", str(e)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--tier", type=int, default=2, help="fetch tiers <= this")
    a = ap.parse_args()
    os.makedirs(RAW, exist_ok=True)
    S = pd.read_csv(SEL)
    S = S[(S.source == "comcat") & (S.tier <= a.tier)]
    ids = list(dict.fromkeys(S.cc_id))
    # yearly summary queries (cached under raw/comcat_summary/) give each event's last-update time
    SD = RAW + "_summary"
    os.makedirs(SD, exist_ok=True)
    for y in range(2003, 2027):
        f = f"{SD}/summary_{y}.geojson"
        if not os.path.exists(f):
            try:
                b = get(SUMQ.format(y, y + 1))
                if b:
                    wr(f, b)
            except RuntimeError as e:
                print("summary", y, "failed", e, flush=True)
                continue
        try:
            for ft in json.load(open(f))["features"]:
                p = ft["properties"]
                # fast path only for events whose only products are uw origin + phase-data (then the event's
                # last-update time is the phase-data product's time); others go through the detail call
                v = (p["updated"], (p.get("types") or "") == ",origin,phase-data," and (p.get("sources") or "") == ",uw,")
                for i in (p.get("ids") or "").strip(",").split(","):
                    if i:
                        SUMMARY[i] = v
        except Exception as e:
            print("summary", y, "unreadable", e, flush=True)
    print(f"summary: {len(SUMMARY)} ids; selected with phase-data flag "
          f"{sum(1 for i in ids if SUMMARY.get(i, (0, False))[1])}/{len(ids)}", flush=True)
    t0 = time.time()
    print(f"fetch start {time.strftime('%Y-%m-%d %H:%M:%S')} : {len(ids)} events, workers {a.workers}", flush=True)
    new = not os.path.exists(LOG)
    with open(LOG, "a") as lf, ThreadPoolExecutor(min(a.workers, 2)) as ex:
        if new:
            lf.write("t_unix,cc_id,status,info\n")
        for n, (cc, st, info) in enumerate(ex.map(one, ids), 1):
            if st != "skip":
                lf.write(f"{time.time():.0f},{cc},{st},\"{str(info)[:120]}\"\n"); lf.flush()
            if n % 50 == 0 or n == len(ids):
                el = time.time() - t0
                print(f"{n}/{len(ids)} {el / 60:.1f} min {stats}", flush=True)
    print("fetch end", stats, flush=True)
    return 0 if stats["fail"] == 0 else 3


if __name__ == "__main__":
    raise SystemExit(main())
