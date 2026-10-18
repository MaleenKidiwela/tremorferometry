#!/usr/bin/env python
"""Build the gate VALIDATION proposal sets on the 2012 truth year.

positive controls
  pos_bost_N{20,50,100}  the 116 templated Bostock families, members = Bostock's own ETS-2012 event
                         times (the "known family" ceiling; also the source of the Gate-1 control table)
  pos_tpl_N{20,50,100}   the same families, members = the FROZEN template detections of the run with
                         VGZ / SNB / B926 EXCLUDED from the network sum (e5r_detections_a3b_noVSB,
                         12 s, T_ext 8.25, ETS days, own-matched) -- so the Gate-1 withheld stations
                         never contributed to the detections they are asked to confirm
negative controls
  neg_moveout_N50        pos_tpl times at a WRONG cell (S offset differs by >= 5 s at every withheld
                         station) -> must fail Gate 1
  neg_timeshift_N50      pos_tpl times shifted by a random 300-3600 s -> must fail Gate 1
  neg_dvv_N50            50 causality-certified dv/v "families" (../data/*_causality_cert.csv) with
                         their stored detection times (../data/mf_*.csv) at their tremor-cell
                         coordinates, Slab2 depth, with a +-30 s constant-offset search (the dv/v
                         detection time is not an origin time; the same search is applied to every null)
alias control
  alias_N50              a neighbouring family's own events proposed as a NEW family 1-2 km away
                         -> must fail Gate 0

Writes results/gates/proposals_<set>.parquet with columns pid, t, lat, lon, dep, offset_search_s, kind,
truth_bkey.
"""
from __future__ import annotations

import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import g_common as G                                             # noqa: E402

OUT = G.GATES
YEAR = 2012
WITHHELD = ["VGZ", "SNB", "B926"]


def ets_days(gt):
    gt = gt.copy()
    gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    n = gt.groupby("doy").size()
    return sorted(n[n >= 100].index)


def subsample(t, n, rng):
    t = np.asarray(t, float)
    if len(t) <= n:
        return t
    return np.sort(rng.choice(t, size=n, replace=False))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=20260917)
    ap.add_argument("--dvv-n", type=int, default=50)
    a = ap.parse_args()
    rng = np.random.default_rng(a.seed)
    os.makedirs(OUT, exist_ok=True)
    S = G.stations()
    F = pd.read_csv(f"{G.REP}/results/e5r_families_2012.csv").set_index("bkey")
    gt = pd.read_parquet(f"{G.REP}/results/e5r_gt_bostock2012.parquet")
    ETS = ets_days(gt)
    gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    gtE = gt[gt.doy.isin(ETS)]
    qc = pd.read_csv(f"{G.REP}/results/e5r_template_qc_2012.csv")
    templated = sorted(int(b) for b in qc[qc.usable & qc.template].bkey.unique())
    det = pd.read_parquet(f"{G.EXT}/results/e5r_detections_a3b_noVSB.parquet")
    det = det[(det.L == 12.0) & (det.net_mad >= 8.25) & det.doy.isin(ETS) & det.match_own]

    def write(rows, name):
        D = pd.DataFrame(rows)
        p = f"{OUT}/proposals_{name}.parquet"
        D.to_parquet(p, index=False)
        print(f"[{name}] {D.pid.nunique()} proposals, {len(D)} member rows -> {p}")
        return D

    # ---------------------------------------------------------------- positives
    for N in (20, 50, 100):
        rows = []
        for b in templated:
            t = gtE[gtE.bkey == b].t_origin.to_numpy(float)
            if len(t) < N:
                continue
            t = subsample(t, N, rng)
            fr = F.loc[b]
            rows += [dict(pid=f"posb{N}_f{b}", t=float(x), lat=float(fr.lat), lon=float(fr.lon),
                          dep=float(fr.dep), offset_search_s=0.0, kind="pos_bostock",
                          truth_bkey=int(b)) for x in t]
        write(rows, f"pos_bost_N{N}")
    for N in (20, 50, 100):
        rows = []
        for b, g in det.groupby("fam"):
            b = int(b)
            if b not in templated or len(g) < N:
                continue
            t = subsample(g.t.to_numpy(float), N, rng)
            fr = F.loc[b]
            rows += [dict(pid=f"post{N}_f{b}", t=float(x), lat=float(fr.lat), lon=float(fr.lon),
                          dep=float(fr.dep), offset_search_s=0.0, kind="pos_template",
                          truth_bkey=int(b)) for x in t]
        write(rows, f"pos_tpl_N{N}")

    # ---------------------------------------------------------------- negatives (i) wrong cell
    base = pd.read_parquet(f"{OUT}/proposals_pos_tpl_N50.parquet")
    glat = np.arange(48.0, 49.3001, 0.02)
    glon = np.arange(-124.8, -122.8001, 0.03)
    LA, LO = np.meshgrid(glat, glon, indexing="ij")
    LA, LO = LA.ravel(), LO.ravel()
    DZ = G.slab2_depth(LA, LO)
    wsta = [(s, float(S.loc[s, "lat"]), float(S.loc[s, "lon"])) for s in WITHHELD]
    ttS_cell = np.empty((len(LA), len(wsta)))
    for k, (s, la, lo) in enumerate(wsta):
        d = G.epi_km(LA, LO, la, lo)
        for i in range(len(LA)):
            ttS_cell[i, k] = G.ttSP(float(d[i]), float(DZ[i]))[0]
    rows = []
    for pid, g in base.groupby("pid"):
        r0 = g.iloc[0]
        tt0 = np.array([G.ttSP(float(G.epi_km(r0.lat, r0.lon, la, lo)), float(r0.dep))[0]
                        for _s, la, lo in wsta])
        ok = np.where((np.abs(ttS_cell - tt0[None, :]) >= 5.0).all(1))[0]
        j = int(rng.choice(ok))
        rows += [dict(pid=pid.replace("post50_", "negmo_"), t=float(x), lat=float(LA[j]),
                      lon=float(LO[j]), dep=float(DZ[j]), offset_search_s=0.0,
                      kind="neg_moveout", truth_bkey=int(r0.truth_bkey)) for x in g.t]
    write(rows, "neg_moveout_N50")

    # ---------------------------------------------------------------- negatives (ii) time shift
    rows = []
    for pid, g in base.groupby("pid"):
        r0 = g.iloc[0]
        sh = float(rng.uniform(300, 3600)) * rng.choice([-1.0, 1.0])
        rows += [dict(pid=pid.replace("post50_", "negts_"), t=float(x) + sh, lat=float(r0.lat),
                      lon=float(r0.lon), dep=float(r0.dep), offset_search_s=0.0,
                      kind="neg_timeshift", truth_bkey=int(r0.truth_bkey)) for x in g.t]
    write(rows, "neg_timeshift_N50")

    # ---------------------------------------------------------------- negatives (iii) dv/v families
    cells = {}
    for f in sorted(glob.glob(f"{G.ROOT}/data/*_causality_cert.csv")):
        sta = os.path.basename(f).split("_")[0]
        if not os.path.exists(f"{G.ROOT}/data/mf_{sta}p90f40_{YEAR}.csv"):
            continue
        C = pd.read_csv(f)
        C = C[C.reliable]
        for fam in C.fam:
            try:
                la, lo = float(fam.split("_")[0]), float(fam.split("_")[1])
            except Exception:
                continue
            if not (48.0 <= la <= 49.3 and -124.9 <= lo <= -122.8):
                continue
            dmin = min(G.epi_km(la, lo, float(S.loc[s, "lat"]), float(S.loc[s, "lon"]))
                       for s in WITHHELD)
            if dmin > 60.0:
                continue
            cells.setdefault((sta, fam), (la, lo, dmin))
    keys = sorted(cells)
    rng.shuffle(keys)
    rows, taken = [], 0
    bysta = {}
    for (sta, fam) in keys:
        bysta.setdefault(sta, []).append(fam)
    for sta in sorted(bysta):
        if taken >= a.dvv_n:
            break
        M = pd.read_csv(f"{G.ROOT}/data/mf_{sta}p90f40_{YEAR}.csv", usecols=["template", "time"])
        M["t"] = pd.to_datetime(M.time).astype("datetime64[ns]").astype("int64") / 1e9
        Mdoy = pd.to_datetime(M.time).dt.dayofyear
        M = M[Mdoy.isin(ETS)]
        for fam in bysta[sta]:
            if taken >= a.dvv_n:
                break
            t = M.t[M.template == fam].to_numpy(float)
            if len(t) < 50:
                continue
            la, lo, _ = cells[(sta, fam)]
            dz = float(G.slab2_depth([la], [lo])[0])
            t = subsample(t, 50, rng)
            rows += [dict(pid=f"dvv_{sta}_{fam}", t=float(x), lat=la, lon=lo, dep=dz,
                          offset_search_s=30.0, kind="neg_dvv", truth_bkey=-1) for x in t]
            taken += 1
    write(rows, "neg_dvv_N50")

    # ---------------------------------------------------------------- alias control
    rows = []
    lat = F.lat.values; lon = F.lon.values; bk = F.index.values
    for pid, g in base.groupby("pid"):
        r0 = g.iloc[0]
        b = int(r0.truth_bkey)
        # a plausible detector location error: 1-2 km from the true (accepted) family
        ang = rng.uniform(0, 2 * np.pi); rad = rng.uniform(1.0, 2.0)
        dlat = rad * np.sin(ang) / 111.2
        dlon = rad * np.cos(ang) / (111.2 * np.cos(np.radians(48.6)))
        rows += [dict(pid=pid.replace("post50_", "alias_"), t=float(x), lat=float(r0.lat + dlat),
                      lon=float(r0.lon + dlon), dep=float(r0.dep), offset_search_s=0.0,
                      kind="alias", truth_bkey=b) for x in g.t]
    write(rows, "alias_N50")


if __name__ == "__main__":
    main()
