#!/usr/bin/env python
"""Track C step 1 — is a month RUNNABLE, and which stations does it contribute?

RUNNABLE  every station-day of the month, for every station whose epoch (results/station_epochs.csv)
          covers that day, has a FINISHED manifest row: status in {data, nodata, error,
          present_elsewhere, no_epoch}.  `error` counts as finished (the downloader retries once and
          then gives up); `nodata` is a verified-empty day.
USABLE    a station contributes to the matched filter on a day when its row is `data` or
          `present_elsewhere` AND coverage_frac >= --min-cov (default 0.5).  present_elsewhere rows
          carry no coverage in the manifest (the file was written by another study); they are
          accepted here and re-checked when the day is cached.

  python scripts/c_month_ready.py 2015-12            # one month
  python scripts/c_month_ready.py --next             # oldest runnable, unprocessed month (for the watcher)
  python scripts/c_month_ready.py --list             # every month with its state
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C

FINISHED = {"data", "nodata", "error", "present_elsewhere", "no_epoch"}
USABLE = {"data", "present_elsewhere"}


def month_state(ym, man, ep, min_cov=0.5):
    days = C.month_days(ym)
    sub = man[man.date.str[:7] == ym]
    by = {(r.netsta, r.date): r for r in sub.itertuples()}
    missing, usable_days, per_sta = [], [], {}
    for d in days:
        ds = d.strftime("%Y-%m-%d")
        al = C.alive(ep, d)
        for r in al.itertuples():
            key = (r.netsta, ds)
            row = by.get(key)
            if row is None or str(row.status) not in FINISHED:
                missing.append((r.netsta, ds, "missing" if row is None else str(row.status)))
                continue
            if str(row.status) in USABLE:
                cov = row.coverage_frac
                ok = (not np.isfinite(cov)) or (cov >= min_cov)       # present_elsewhere: cov is NaN
                if ok:
                    per_sta.setdefault(r.sta, []).append(ds)
                    usable_days.append((r.sta, ds, r.band, r.epoch_id, float(cov) if np.isfinite(cov) else np.nan))
    return dict(ym=ym, n_days=len(days), runnable=len(missing) == 0, n_missing=len(missing),
                missing=missing[:20], stations=sorted(per_sta), n_station_days=len(usable_days),
                usable=usable_days)


def months_done():
    p = f"{C.TRC}/months_done.csv"
    if not os.path.exists(p):
        return set()
    return set(pd.read_csv(p).ym.astype(str))


def all_months(man):
    return sorted(man.date.str[:7].unique())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ym", nargs="?", default="")
    ap.add_argument("--next", action="store_true", help="print the oldest runnable, unprocessed month")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--min-cov", type=float, default=0.5)
    ap.add_argument("--json", action="store_true")
    a = ap.parse_args()
    man = C.manifest(); ep = C.epochs()
    if a.next or a.list:
        done = months_done()
        out = []
        for ym in all_months(man):
            st = month_state(ym, man, ep, a.min_cov)
            st["done"] = ym in done
            out.append(st)
            if a.next and st["runnable"] and not st["done"]:
                print(ym)
                return
        if a.next:
            print("")          # nothing runnable
            return
        for st in out:
            print(f"{st['ym']}  runnable={int(st['runnable'])} done={int(st['done'])} "
                  f"missing={st['n_missing']:4d} station-days={st['n_station_days']:4d} "
                  f"stations={len(st['stations'])}")
        return
    if not a.ym:
        ap.error("give a month (YYYY-MM) or --next / --list")
    st = month_state(a.ym, man, ep, a.min_cov)
    if a.json:
        print(json.dumps({k: v for k, v in st.items() if k != "usable"}, indent=1, default=str))
    else:
        print(f"[{st['ym']}] runnable={st['runnable']} ({st['n_missing']} unfinished station-days) "
              f"usable station-days={st['n_station_days']} over {len(st['stations'])} stations")
        print("  stations: " + ",".join(st["stations"]))
        for m in st["missing"]:
            print(f"  UNFINISHED {m[0]} {m[1]} {m[2]}")
    sys.exit(0 if st["runnable"] else 3)


if __name__ == "__main__":
    main()
