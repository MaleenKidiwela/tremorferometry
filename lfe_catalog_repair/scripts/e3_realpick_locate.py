#!/usr/bin/env python
"""E3 step 2, config D — REAL LFEDetect picks at events whose answer is known.

Config C (e3_geom_synth.py) shows that at N >= 4 the S-only geometry is fine: with 0.3 s synthetic
pick noise the median horizontal error is 3.7 km, not the 25 km E2 measured on its N=4 associated
events.  Something other than geometry costs the remaining factor ~7, and there are two candidates:
(i) the real LFEDetect picks are not Gaussian at 0.3 s, or (ii) PyOcto assembled the wrong picks
into events.  Config D separates them: it uses the REAL picks but skips PyOcto entirely.

For every Bostock event in the window, the predicted S arrival at each station is
    (DB time) + (per-family offset, E3 step 5) + ttS(family hypocentre, Crosson 1-D).
Any real LFEDetect S pick within `--tol` of that prediction, at a station with epi < `--max-dist`,
is taken as that event's pick at that station.  Events with >= 3 such stations are written through
the UNMODIFIED `e2_locate_nll.py` machinery (S only, the same 0.4 s GAU card, the same grids) and
compared with the published family hypocentre.

C -> D is the cost of real (non-Gaussian) pick error; D -> E2 is the cost of PyOcto association.

Outputs results/e3_realpick_locate.csv, results/e3_realpick_summary.csv
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
BIN = "/home/jovyan/.local/bin"
OVERLAP = ["PGC", "LZB", "VGZ", "YOUB", "SNB", "PFB"]
R_KM_DEG = 111.195
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_locate_nll as NLL  # noqa: E402


def epoch(s):
    return pd.to_datetime(s).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default="cascadia")
    ap.add_argument("--threshold", type=float, default=0.3)
    ap.add_argument("--tol", type=float, default=1.0)
    ap.add_argument("--max-dist", type=float, default=40.0)
    ap.add_argument("--min-stations", type=int, default=3)
    args = ap.parse_args()

    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    fo = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv")
    acc = fo[fo.accept.astype(bool)]
    omap = dict(zip(acc.bkey.astype(int), acc.offset_s.astype(float)))
    bo = bo[bo.bkey.astype(int).isin(omap)].copy()
    bo["t"] = epoch(bo.time) + bo.bkey.astype(int).map(omap).to_numpy()

    pk = pd.read_parquet(f"{RES}/e2_picks_vi2010_{args.model}.parquet")
    pk = pk[(pk.phase == "S") & (pk.thr_first <= args.threshold) & (pk.prob >= args.threshold)]
    t0, t1 = epoch(pk.time).min(), epoch(pk.time).max()
    bo = bo[(bo.t >= t0 - 60) & (bo.t <= t1 + 60)].reset_index(drop=True)
    print(f"[gt] {len(bo)} Bostock events in the vi2010 window with an accepted family offset")

    rows = []
    for s in OVERLAP:
        g = pk[pk.sta == s]
        if len(g) == 0:
            continue
        ptimes = np.sort(epoch(g.time))
        pred = bo.t.to_numpy() + bo[f"ttS_{s}"].to_numpy()
        near = bo[f"epi_{s}"].to_numpy() <= args.max_dist
        i = np.searchsorted(ptimes, pred)
        best = np.full(len(pred), np.inf)
        bt = np.full(len(pred), np.nan)
        for off in (-1, 0):
            j = np.clip(i + off, 0, len(ptimes) - 1)
            d = ptimes[j] - pred
            m = np.abs(d) < np.abs(best)
            best = np.where(m, d, best); bt = np.where(m, ptimes[j], bt)
        ok = near & (np.abs(best) <= args.tol)
        rows.append(pd.DataFrame(dict(ev=np.where(ok)[0], station=s, t=bt[ok], resid=best[ok])))
        print(f"  {s}: {int(ok.sum())} of {int(near.sum())} near-station events matched "
              f"within {args.tol} s (residual sd {np.std(best[ok]):.2f} s)")
    A = pd.concat(rows, ignore_index=True)
    cnt = A.groupby("ev").station.nunique()
    keep = cnt[cnt >= args.min_stations].index
    print(f"[sel] events with >= 3 / 4 / 5 matched overlap stations: "
          f"{int((cnt >= 3).sum())} / {int((cnt >= 4).sum())} / {int((cnt >= 5).sum())}")
    A = A[A.ev.isin(keep)].copy()

    root = os.path.abspath(f"{RES}/nll_e3_realpick_t{args.threshold}")
    os.makedirs(root, exist_ok=True)
    stas = pd.read_csv(f"{RES}/e2_stations.csv")
    stas = stas[stas.sta.isin(OVERLAP)]
    NLL.write_control(root, stas, "obs.txt")
    NLL.build_grids(root)
    ids = sorted(A.ev.unique())
    with open(f"{root}/obs.txt", "w") as f:
        for eid in ids:
            for r in A[A.ev == eid].itertuples():
                tt = pd.Timestamp(r.t, unit="s")
                f.write(f"{r.station:<6s} ?    ?    ? {'S':<6s} ? "
                        f"{tt.strftime('%Y%m%d %H%M')} "
                        f"{tt.second + tt.microsecond/1e6:9.4f} GAU {0.4:9.2e} -1.00e+00 "
                        f"-1.00e+00 -1.00e+00 1.00\n")
            f.write("\n")
    print(f"[obs] {len(ids)} events written")
    r = subprocess.run([f"{BIN}/NLLoc", "nlloc.in"], capture_output=True, text=True,
                       errors="replace", cwd=root)
    print(f"[NLLoc] rc={r.returncode} " + " | ".join((r.stdout + r.stderr).strip().split("\n")[-2:]))
    loc = NLL.parse_sum(f"{root}/loc/lfe.sum.grid0.loc.hyp")
    for fn in os.listdir(f"{root}/loc"):
        if fn.endswith((".buf", ".hdr", ".scat", ".octree")) and ".sum." not in fn:
            os.remove(f"{root}/loc/{fn}")
    print(f"[loc] {len(loc)} of {len(ids)} located")
    if len(loc) != len(ids):
        lt = loc.time_utc.astype("datetime64[ns]").astype("int64").to_numpy() / 1e9
        et = bo.t.to_numpy()[ids]
        j = np.argmin(np.abs(lt[:, None] - et[None, :]), axis=1)
        loc["ev"] = np.asarray(ids)[j]
    else:
        loc["ev"] = ids
    truth = bo.iloc[loc.ev.to_numpy()]
    cos0 = np.cos(np.radians(48.6))
    loc["n_sta"] = cnt.reindex(loc.ev).to_numpy()
    loc["dh_km"] = np.hypot((loc.lon.to_numpy() - truth.lon.to_numpy()) * R_KM_DEG * cos0,
                            (loc.lat.to_numpy() - truth.lat.to_numpy()) * R_KM_DEG)
    loc["dz_km"] = loc.depth.to_numpy() - np.abs(truth.dep.to_numpy())
    loc["bkey"] = truth.bkey.to_numpy()
    loc.to_csv(f"{RES}/e3_realpick_locate_t{args.threshold}.csv", index=False)
    S = loc.groupby("n_sta").agg(n=("dh_km", "size"), dh_med=("dh_km", "median"),
                                 dh_q25=("dh_km", lambda v: v.quantile(.25)),
                                 dh_q75=("dh_km", lambda v: v.quantile(.75)),
                                 dz_med=("dz_km", "median"),
                                 frac_shallow=("depth", lambda v: float((v < 5).mean())),
                                 errh=("err_h_km", "median"), rms=("rms", "median")).reset_index()
    S["threshold"] = args.threshold
    S.to_csv(f"{RES}/e3_realpick_summary_t{args.threshold}.csv", index=False)
    print("\n[config D: REAL picks, known events, no PyOcto]")
    print(S.round(2).to_string(index=False))
    print(f"\n[write] {RES}/e3_realpick_locate_t{args.threshold}.csv, "
          f"e3_realpick_summary_t{args.threshold}.csv")


if __name__ == "__main__":
    main()
