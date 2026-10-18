#!/usr/bin/env python
"""S11 tomo (PLAN.md Amendment 8): SVI 2003-2006 era-detector family stacks incl. POLARIS.

  python s11_tomo_era_stacks.py --workers 6
writes tomo/tables/era_stacks/<NETSTA>.npz (keys fam -> (3, 480) mean stack, n), tables/era_stack_index.csv,
marker tomo/ERA_STACKS_DONE.  One worker per station; day files read once per station.
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
T = f"{ROOT}/svi_science/results/s11_depth/tomo"
OUT = f"{T}/tables/era_stacks"
CACHE = f"{ROOT}/lfe_extension/cache20hz_backward"
sys.path.insert(0, f"{ROOT}/lfe_catalog_repair/scripts")
import d3_vel1d as v1d  # noqa: E402

FS = 20.0
PRE, LEN = 16.0, 24.0
NS = int(LEN * FS)
STATIONS = ["C8.KHVB", "C8.MGCB", "C8.SHDB", "C8.SHVB", "C8.TWBB", "PO.KELB", "PO.KLNB", "PO.SILB", "PO.SSIB",
            "PO.TSJB", "PO.TWGB", "PO.TWKB", "CN.LZB", "CN.PGC", "CN.VGZ", "CN.SNB", "CN.NLLB", "PB.B001", "PB.B003",
            "PB.B004", "PB.B005", "PB.B006", "PB.B007", "PB.B009", "PB.B010", "PB.B011"]


def detections():
    C = pd.read_parquet(f"{ROOT}/data/catalogue/svi_catalogue_v1_4_backward.parquet",
                        columns=["fam", "t", "P", "series", "counted_once", "date"])
    c = C[(C.series == "data") & C.counted_once & (C.P >= 0.75) & (C.date <= "2006-06-30")].copy()
    c["fam"] = c.fam.astype(int)
    return c


def station_coords():
    I = pd.read_csv(f"{ROOT}/svi_science/results/backward_scope/inventory_stations.csv")
    I["netsta"] = I.net + "." + I.sta
    return I.drop_duplicates("netsta").set_index("netsta")[["lat", "lon"]]


def run_station(args):
    ns, D, F, xy = args
    lat, lon = xy
    ttS = {}
    for f, r in F.iterrows():
        dx = np.hypot((lon - r.lon) * 111.195 * np.cos(np.radians(r.lat)), (lat - r.lat) * 111.195)
        ttS[f] = float(v1d.tt(dx, r.dep, "S"))
    W = {f: [] for f in F.index}
    D = D[D.fam.isin(F.index)]
    nday = 0
    for date, g in D.groupby("date"):
        y, doy = int(date[:4]), pd.Timestamp(date).dayofyear
        p = f"{CACHE}/{ns}/{y}/{doy:03d}.npz"
        if not os.path.exists(p):
            continue
        z = np.load(p)
        x, t0 = z["x"], float(z["t0"])
        order = [str(o) for o in z["order"]]
        if len(order) < 3 or not order[0].endswith("Z"):
            continue
        nday += 1
        for r in g.itertuples():
            i0 = int(round((r.t + ttS[r.fam] - PRE - t0) * FS))
            if i0 < 0 or i0 + NS > x.shape[1]:
                continue
            seg = x[:, i0:i0 + NS].astype(np.float64)
            if not np.isfinite(seg).all():
                continue
            seg = seg - seg.mean(1, keepdims=True)
            rms = np.sqrt((seg ** 2).mean())
            if rms <= 0:
                continue
            W[r.fam].append((seg / rms, rms))
    out, idx = {}, []
    for f, L in W.items():
        if len(L) < 10:
            continue
        rr = np.array([b for _, b in L]); med = np.median(rr)
        A = np.array([a for a, b in L if b <= 5 * med])
        out[f"{f}"] = A.mean(0).astype(np.float32)
        idx.append(dict(netsta=ns, fam=f, n=len(A), n_drop=len(L) - len(A), offset=ttS[f] - PRE))
    np.savez_compressed(f"{OUT}/{ns}.npz", **out)
    return ns, nday, idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    D = detections()
    F = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv").set_index("bkey")
    F = F.loc[F.index.isin(D.fam.unique())]
    XY = station_coords()
    jobs = [(ns, D, F, (float(XY.loc[ns, "lat"]), float(XY.loc[ns, "lon"]))) for ns in STATIONS
            if not os.path.exists(f"{OUT}/{ns}.npz")]
    print(f"{len(D)} detections, {F.shape[0]} families, {len(jobs)} stations to do", flush=True)
    t0 = time.time(); rows = []
    import multiprocessing as mp
    with mp.get_context("fork").Pool(a.workers) as p:
        for ns, nday, idx in p.imap_unordered(run_station, jobs):
            rows += idx
            pd.DataFrame(idx).to_csv(f"{OUT}/{ns}_index.csv", index=False)
            print(f"{ns}: {nday} days, {len(idx)} family stacks  {time.time() - t0:.0f}s", flush=True)
    I = pd.concat([pd.read_csv(f"{OUT}/{ns}_index.csv") for ns in STATIONS if os.path.exists(f"{OUT}/{ns}_index.csv")
                   and os.path.getsize(f"{OUT}/{ns}_index.csv") > 1])
    I.to_csv(f"{T}/tables/era_stack_index.csv", index=False)
    open(f"{T}/ERA_STACKS_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
