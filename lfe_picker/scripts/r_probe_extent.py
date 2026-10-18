#!/usr/bin/env python
"""P3D step 1c - map the ACTUAL data extent of a station month by month.

The FDSN availability service is retired (service.iris.edu and service.earthscope.org both answer
"This service has been retired"), and StationXML epochs lie: PQ.JRDRB's channel epoch opens
2021-04-01 but the archive holds nothing until 2025.  So: ask dataselect for a few minutes on
three days of every month of the epoch and record what comes back.  Cheap (a few hundred tiny
requests) and empirical, which is the only thing the GOBB lesson lets us trust.

  python scripts/r_probe_extent.py --stations PQ.JRDRB,PQ.LHLYB,... [--workers 8]
-> results/p3d/probe_month_extent.csv   net,sta,month,days_probed,days_with_data,bytes
"""
import os
import csv
import argparse
import urllib.request
import urllib.error
from datetime import date, timedelta
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results", "p3d")
IRIS = "https://service.iris.edu/fdsnws/dataselect/1/query"
PROBE_DAYS = (5, 15, 25)
PROBE_MIN = 10


def one(args):
    net, sta, band, d = args
    url = (f"{IRIS}?net={net}&sta={sta}&loc=*&cha={band}Z"
           f"&start={d.isoformat()}T12:00:00&end={d.isoformat()}T12:{PROBE_MIN:02d}:00")
    for _ in range(3):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                return (net, sta, d, len(r.read()))
        except urllib.error.HTTPError as e:
            if e.code in (204, 404):
                return (net, sta, d, 0)
        except Exception:
            pass
    return (net, sta, d, -1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", required=True)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--band", default="HH")
    ap.add_argument("--start", default="")
    ap.add_argument("--end", default="2026-09-15")
    a = ap.parse_args()
    ep = {}
    with open(os.path.join(RESULTS, "probe_summary.csv")) as f:
        for r in csv.DictReader(f):
            ep[f"{r['net']}.{r['sta']}"] = (r.get("epoch_start", ""), r.get("epoch_end", ""))
    jobs = []
    for k in a.stations.split(","):
        net, sta = k.split(".")
        s = date.fromisoformat(a.start or ep.get(k, ("2019-01-01",))[0] or "2019-01-01")
        e = min(date.fromisoformat(a.end),
                date.fromisoformat(ep.get(k, ("", a.end))[1] or a.end))
        d = date(s.year, s.month, 1)
        while d <= e:
            for dd in PROBE_DAYS:
                try:
                    t = d.replace(day=dd)
                except ValueError:
                    continue
                if s <= t <= e:
                    jobs.append((net, sta, a.band, t))
            d = (d.replace(day=28) + timedelta(days=5)).replace(day=1)
    print(f"{len(jobs)} probes, {a.workers} workers", flush=True)
    res = []
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        for i, r in enumerate(ex.map(one, jobs), 1):
            res.append(r)
            if i % 100 == 0:
                print(f"  {i}/{len(jobs)}", flush=True)
    agg = {}
    for net, sta, d, nb in res:
        k = (net, sta, d.strftime("%Y-%m"))
        a_ = agg.setdefault(k, [0, 0, 0, 0])
        a_[0] += 1
        if nb > 512:
            a_[1] += 1
            a_[2] += nb
        if nb < 0:
            a_[3] += 1
    p = os.path.join(RESULTS, "probe_month_extent.csv")
    with open(p, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["net", "sta", "month", "days_probed", "days_with_data", "bytes", "errors"])
        for (net, sta, m), v in sorted(agg.items()):
            w.writerow([net, sta, m, v[0], v[1], v[2], v[3]])
    print("->", p)
    # compact per-station picture
    for sta in sorted({k[1] for k in agg}):
        ms = sorted(k[2] for k in agg if k[1] == sta)
        line = "".join("#" if agg[(k[0], sta, m)][1] == agg[(k[0], sta, m)][0] else
                       ("+" if agg[(k[0], sta, m)][1] else ".")
                       for m in ms for k in [next(x for x in agg if x[1] == sta and x[2] == m)])
        live = [m for m in ms if agg[next(x for x in agg if x[1] == sta and x[2] == m)][1]]
        print(f"  {sta:<8} {ms[0]}..{ms[-1]}  {line}")
        print(f"  {'':<8} first data {live[0] if live else '-'}, last {live[-1] if live else '-'}, "
              f"{len(live)}/{len(ms)} months with data")


if __name__ == "__main__":
    main()
