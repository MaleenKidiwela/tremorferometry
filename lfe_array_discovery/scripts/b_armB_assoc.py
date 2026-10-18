#!/usr/bin/env python
"""ARM B stage 3 -- PyOcto association of the LFEDetect(jms) picks into events.  Prefix `b_`.

1-D velocity model: Crosson (1976) western-Washington layers with Vp/Vs = 1.78 -- the SAME model
`lfe_catalog_repair/scripts/e2_associate.py` used for the southern-Vancouver-Island LFEDetect
association (`vi2010`/`vi2017`), reused verbatim so Arm B is not a new methodological choice.

Depth is restricted to the Slab2 interface band 28-45 km as a PRIOR (PLAN sec.10 active-cell rule
ACTIVE_CELL_DEP).  The search box is derived from the DETECTION-STATION geometry (bbox + 0.3 deg),
never from the Bostock family locations: the hidden families must not touch the search.

S picks are the primary observable (Lin measured sigma_P = 4.2 s for LFE P picks); P picks are
supplied as well but no P-only event is accepted.

MINIMUM STATION COUNT = 3, which is e2_associate's documented rule for Vancouver Island (">= 3
stations inside VI, >= 4 elsewhere") and the count E2's own vi2010 association ran at (median
n_stations = 3 over its 9,515 events).  A 2-day pilot on this bed measured the cliff: at the
default threshold 0.30 the 4-station rule yields ZERO events and the 3-station rule 15, at 0.15
the 3-station rule yields 144 (n_stations histogram 113/28/3) -- so the 4-station variant is run
as a DIAGNOSTIC alongside and reported, never as the gate.

SITE-SHIFT NULL (--shift, stage 7).  `w_common.DayLoader` maps a NOMINAL day index i of a site
with offset k to the SOURCE day index (i + k) % n of the bed-day list, so a pick really recorded
on bed day j appears, under the shift, at nominal day (j - k) % n.  The picks are moved exactly
that way (whole days, members of a site together), which is bit-equivalent to re-running
LFEDetect on the shifted day-series -- the model is deterministic per station-day.

Usage
  python scripts/b_armB_assoc.py --picks "results/armB/picks_*.parquet" --thr 0.3 \
         --out results/armB/events_t30 [--shift '{"GLBC": 7, ...}']
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
os.environ.setdefault("W_BED", "c")
import w_common as W                                                        # noqa: E402

# Crosson (1976) Puget-lowland layers, verbatim from lfe_catalog_repair/scripts/e2_associate.py
VP = [5.40, 6.38, 6.59, 6.73, 6.86, 7.00, 7.90]
ZTOP = [0.0, 4.0, 9.0, 16.0, 20.0, 25.0, 32.0]
VPVS = 1.78
ZLIM = (28.0, 45.0)                       # Slab2 interface band prior (ACTIVE_CELL_DEP)
BOX_PAD_DEG = 0.3


def velocity_model(path, delta=1.0, xdist=260.0, zdist=80.0):
    import pyocto
    if not os.path.exists(path):
        z = np.arange(0.0, zdist + delta, delta)
        vp = np.array([VP[max(np.searchsorted(ZTOP, zi, side="right") - 1, 0)] for zi in z])
        prof = pd.DataFrame(dict(depth=z, vp=vp, vs=vp / VPVS))
        pyocto.VelocityModel1D.create_model(prof, delta, xdist, zdist, path)
    return pyocto.VelocityModel1D(path, tolerance=2.0, association_cutoff_distance=250.0)


def load_picks(pattern, thr, det, kinds=("ets",), shift=None):
    fr = [pd.read_parquet(p) for p in sorted(glob.glob(pattern))]
    D = pd.concat(fr, ignore_index=True) if fr else pd.DataFrame()
    if not len(D):
        raise SystemExit(f"no picks matched {pattern}")
    D = D[D.sta.isin(det) & D.kind.isin(kinds)]
    D = D[(D.thr_first <= thr) & (D.prob >= thr)].copy()
    t = pd.to_datetime(D.time).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9
    if shift:
        days = [pd.Timestamp(x) for x in W.bed_days()]
        base = np.array([d.timestamp() for d in days])
        doys = [int(d.dayofyear) for d in days]
        n = len(days)
        of = W.sites()[1]
        idx_of_doy = {d: i for i, d in enumerate(doys)}
        j = np.array([idx_of_doy.get(int(x), -1) for x in D.doy.to_numpy()])
        k = np.array([int(shift.get(of.get(s), 0)) for s in D.sta])
        ok = j >= 0
        i_new = np.where(ok, (j - k) % n, j)
        t = np.where(ok, t - base[np.clip(j, 0, n - 1)] + base[np.clip(i_new, 0, n - 1)], t)
        W.log(f"[shift] {int(ok.sum())}/{len(D)} picks moved by the site-day shift {shift}")
    out = pd.DataFrame(dict(station=D.sta.to_numpy(), phase=D.phase.to_numpy(), time=t,
                            probability=D.prob.to_numpy()))
    return out.sort_values("time").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--picks", default="results/armB/picks_*.parquet")
    ap.add_argument("--thr", type=float, default=0.3)
    ap.add_argument("--phases", default="PS", choices=["PS", "S"])
    ap.add_argument("--n-picks", type=int, default=3)
    ap.add_argument("--min-stations", type=int, default=3)
    ap.add_argument("--pick-tol", type=float, default=2.5)
    ap.add_argument("--pick-std", type=float, default=1.5)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--zlim", default="", help="override the depth prior, e.g. \"0,70\"")
    ap.add_argument("--time-before", type=float, default=120.0)
    ap.add_argument("--shift", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    import pyocto

    shift = json.loads(a.shift) if a.shift else None
    det = list(W.sites()[3])
    ST = W.station_table()
    st = pd.DataFrame([dict(id=s, latitude=float(ST.lat[s]), longitude=float(ST.lon[s]),
                            elevation=float(ST.elev_m[s])) for s in det if s in ST.index])
    picks = load_picks(a.picks, a.thr, det, shift=shift)
    if a.phases == "S":
        picks = picks[picks.phase == "S"]
    picks = picks.reset_index(drop=True)
    W.log(f"[picks] {len(picks)} picks >= {a.thr}, {picks.station.nunique()} stations, "
          f"{picks.groupby('phase').size().to_dict()}")

    lat = (float(st.latitude.min()) - BOX_PAD_DEG, float(st.latitude.max()) + BOX_PAD_DEG)
    lon = (float(st.longitude.min()) - BOX_PAD_DEG, float(st.longitude.max()) + BOX_PAD_DEG)
    zlim = tuple(float(v) for v in a.zlim.split(",")) if a.zlim else ZLIM
    W.log(f"[area] lat {lat} lon {lon} depth {zlim} km (Slab2 interface band prior)")
    vm = velocity_model(f"{W.RES}/armB/pyocto_crosson", 1.0, 260.0, 80.0)
    assoc = pyocto.OctoAssociator.from_area(
        lat=lat, lon=lon, zlim=zlim, velocity_model=vm, time_before=a.time_before,
        n_picks=a.n_picks, n_p_picks=0, n_s_picks=a.n_picks, n_p_and_s_picks=0,
        pick_match_tolerance=a.pick_tol, edt_pick_std=a.pick_std,
        min_interevent_time=3.0, n_threads=a.threads)
    assoc.transform_stations(st)
    ev, asg = assoc.associate(picks, st)
    W.log(f"[assoc] {len(ev)} events, {len(asg)} assigned picks")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    if len(ev) == 0:
        pd.DataFrame(columns=["idx", "time", "latitude", "longitude", "depth", "n_stations"]) \
            .to_csv(a.out + ".csv", index=False)
        json.dump(dict(n_events=0, n_assigned=0, n_picks=int(len(picks)), thr=a.thr,
                       shift=shift), open(a.out + ".json", "w"), indent=1)
        return
    assoc.transform_events(ev)
    asg = asg.rename(columns={"event_idx": "event"})
    cnt = asg.groupby("event").agg(
        n_picks=("pick_idx", "count"), n_stations=("station", "nunique"),
        n_S=("phase", lambda s: int((s == "S").sum())),
        n_P=("phase", lambda s: int((s == "P").sum())),
        mean_prob=("probability", "mean"),
        rms_resid=("residual", lambda r: float(np.sqrt(np.mean(np.square(r))))))
    ev = ev.join(cnt, on="idx")
    ev["time_utc"] = pd.to_datetime(ev.time, unit="s")
    keep = ev.n_stations >= a.min_stations
    W.log(f"[filter] {int(keep.sum())}/{len(ev)} events with >= {a.min_stations} distinct stations")
    ev = ev[keep].copy()
    asg = asg[asg.event.isin(set(ev.idx))].copy()
    ev.to_csv(a.out + ".csv", index=False)
    asg.to_parquet(a.out + "_assign.parquet", index=False)
    json.dump(dict(n_events=int(len(ev)), n_assigned=int(len(asg)), n_picks=int(len(picks)),
                   thr=float(a.thr), phases=a.phases, min_stations=int(a.min_stations),
                   area=dict(lat=lat, lon=lon, zlim=list(zlim)), shift=shift,
                   velocity_model="Crosson 1976 Puget layers, Vp/Vs 1.78 (e2_associate)",
                   median_n_stations=float(ev.n_stations.median()),
                   median_rms_resid=float(ev.rms_resid.median()),
                   depth_median=float(ev.depth.median())),
              open(a.out + ".json", "w"), indent=1)
    W.log(f"[write] {a.out}.csv  ({len(ev)} events); n_stations hist "
          f"{ev.n_stations.value_counts().sort_index().to_dict()}")


if __name__ == "__main__":
    main()
