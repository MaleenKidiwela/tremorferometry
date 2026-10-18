#!/usr/bin/env python
"""P4-0a - the native-rate 2-8 Hz / 10-19 Hz pass on the NORTH calibration station-days.

Merlin's contrast bar is specified in 2-8 Hz.  The 20 Hz cache is band-passed 1-8 Hz, so the
cache-based contrast is really a 1-8 Hz number and is DILUTED by the 1-2 Hz microseism band.  This
runs u_fartremor.hf_one (native rate, four 1.5 h blocks, 6-s median window power) on the north
station-days so that the SAME native 2-8 Hz statistic exists on both sides of the calibration.
"""
from __future__ import annotations
import os
import sys
import argparse
import warnings
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import u_common as U                                                   # noqa: E402
from u_fartremor import hf_one                                         # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=24)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    P = pd.read_csv(os.path.join(U.P4, "north_stations.csv"))
    PD = pd.read_csv(os.path.join(U.P4, "north_probe_days.csv"))
    jobs = []
    for r in P.itertuples():
        for d in PD.date:
            dd = pd.Timestamp(d)
            if os.path.exists(U.trace_path(r.netsta, dd.year, int(dd.dayofyear))):
                jobs.append((r.netsta, d, dd.year, int(dd.dayofyear), float(r.zsign)))
    print(f"[hf-north] {len(jobs)} station-days", flush=True)
    import multiprocessing as mp
    with mp.Pool(a.workers, maxtasksperchild=8) as pool:
        H = pd.DataFrame(list(pool.imap_unordered(hf_one, jobs, chunksize=1)))
    H["day_type"] = H.date.map(dict(zip(PD.date, PD.day_type)))
    H.sort_values(["netsta", "date"]).to_csv(
        os.path.join(U.P4, "north_hf_station_day_stats.csv"), index=False)
    print("[hf-north] " + str(H.groupby("status").size().to_dict()), flush=True)


if __name__ == "__main__":
    main()
