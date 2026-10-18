#!/usr/bin/env python
"""E2 step 3a — PyOcto association of the LFEDetect picks into events.

1-D velocity model: Crosson (1976) western-Washington layers with Vp/Vs = 1.78, the same model used
for the D3 stage-4 locations (`d3_vel1d.py`) -- stated explicitly because the project has no
velocity model of its own.  casc1.6 is available but as a 17.7 GB ASCII 3-D block; a 3-D
association model is not what PyOcto wants and the 1-D model is what Muenchmeyer's released
workflow used, so the 3-D model is left to the location step.

Association rules (Lin's coincidence rule is the sanity floor, PyOcto adds travel-time consistency):
  * S picks are the primary observable (LFE P picks are known to be poor: Lin measured sigma_P =
    4.2 s and discarded them).  P picks are supplied to PyOcto as well but with a wider tolerance.
  * a candidate event needs >= 3 stations inside VI, >= 4 elsewhere; both are reported.

Usage:
  python lfe_catalog_repair/scripts/e2_associate.py --set vi2010 --model cascadia --threshold 0.3
"""
from __future__ import annotations

import argparse
import glob
import os

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
VP = [5.40, 6.38, 6.59, 6.73, 6.86, 7.00, 7.90]
ZTOP = [0.0, 4.0, 9.0, 16.0, 20.0, 25.0, 32.0]
VPVS = 1.78
AREAS = {
    "vi2010": dict(lat=(47.7, 49.6), lon=(-125.4, -122.3), zlim=(0.0, 70.0)),
    "vi2017": dict(lat=(47.7, 49.6), lon=(-125.4, -122.3), zlim=(0.0, 70.0)),
}


def velocity_model(path, delta=1.0, xdist=260.0, zdist=80.0):
    import pyocto
    if not os.path.exists(path):
        z = np.arange(0.0, zdist + delta, delta)
        vp = np.interp(z, ZTOP, VP, left=VP[0], right=VP[-1])
        # piecewise-constant layers, not interpolated ramps
        vp = np.array([VP[max(np.searchsorted(ZTOP, zi, side="right") - 1, 0)] for zi in z])
        prof = pd.DataFrame(dict(depth=z, vp=vp, vs=vp / VPVS))
        pyocto.VelocityModel1D.create_model(prof, delta, xdist, zdist, path)
    return pyocto.VelocityModel1D(path, tolerance=2.0, association_cutoff_distance=250.0)


def load_picks(setname, model, thr, stations, shift_null=None,
               shift_min=300.0, shift_max=3000.0, shift_mode="station"):
    """E3 step 1 (2026-09-11): `shift_null` = an integer SEED enables the CHANCE-ASSOCIATION NULL.

    Every pick at a given station is moved by ONE constant offset, drawn independently per station,
    uniform in +-[shift_min, shift_max] s (sign random, magnitude uniform).  A constant per-station
    shift preserves each station's pick count, its diurnal/day-to-day rate and therefore the ETS
    modulation exactly, while destroying real inter-station coincidence: any event PyOcto still
    finds is combinatorial.  Picks are NOT dropped when the shift moves them across a UTC day
    boundary (the brief: "picks that cross a day boundary stay in the set"), so the total pick count
    is identical to the observed run.  Default (None) = the original E2 behaviour.
    """
    fr = []
    for p in sorted(glob.glob(f"{RES}/e2_picks_{setname}_*{model}.parquet")):
        fr.append(pd.read_parquet(p))
    df = pd.concat(fr, ignore_index=True)
    df = df[(df.model == model) & (df.thr_first <= thr) & (df.prob >= thr)]
    df = df[df.sta.isin(stations)]
    out = pd.DataFrame(dict(
        station=df.sta.values,
        phase=df.phase.values,
        time=pd.to_datetime(df.time).astype("datetime64[ns]").astype("int64").to_numpy() / 1e9,
        probability=df.prob.values))
    if shift_null is not None:
        rng = np.random.default_rng(int(shift_null))
        stas_sorted = sorted(out.station.unique())          # seed -> same shifts for any threshold
        if shift_mode == "station":
            mag = rng.uniform(shift_min, shift_max, len(stas_sorted))
            sgn = rng.choice([-1.0, 1.0], len(stas_sorted))
            off = dict(zip(stas_sorted, mag * sgn))
            print("[shift-null] mode=station seed %d offsets (s): %s"
                  % (int(shift_null), {k: round(v, 1) for k, v in off.items()}))
            out["time"] = out.time.values + out.station.map(off).to_numpy()
        elif shift_mode == "station_hour":
            # Muenchmeyer et al. (2024) sec. 4 control: "constant shifts for each station for every
            # hour".  A per-station CONSTANT over the whole record also breaks the co-elevation of
            # the station rates inside a tremor burst, which DEFLATES the chance rate; re-drawing
            # the shift every hour keeps the burst envelope aligned at the few-minute scale, so
            # this variant is the CONSERVATIVE (higher-chance) end of the null.
            hr = np.floor(out.time.values / 3600.0).astype("int64")
            si = pd.Series(out.station.values).map(
                {s_: i for i, s_ in enumerate(stas_sorted)}).to_numpy()
            h0 = hr.min()
            nh = int(hr.max() - h0) + 1
            mag = rng.uniform(shift_min, shift_max, (len(stas_sorted), nh))
            sgn = rng.choice([-1.0, 1.0], (len(stas_sorted), nh))
            tab = mag * sgn
            print("[shift-null] mode=station_hour seed %d, %d stations x %d hours, |shift| in "
                  "[%.0f, %.0f] s" % (int(shift_null), len(stas_sorted), nh,
                                      shift_min, shift_max))
            out["time"] = out.time.values + tab[si, hr - h0]
        else:
            raise SystemExit(f"unknown shift_mode {shift_mode}")
    return out.sort_values("time").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--set", dest="setname", required=True)
    ap.add_argument("--model", default="cascadia")
    ap.add_argument("--threshold", type=float, default=0.3)
    ap.add_argument("--phases", default="PS", choices=["PS", "S"])
    ap.add_argument("--n-stations", type=int, default=3)
    ap.add_argument("--pick-tol", type=float, default=2.5)
    ap.add_argument("--pick-std", type=float, default=1.5)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min-stations", type=int, default=3,
                    help="post-filter: minimum number of DISTINCT stations per event")
    ap.add_argument("--out", default="")
    ap.add_argument("--shift-null", type=int, default=None,
                    help="E3 step 1: SEED for the chance-association null (per-station constant "
                         "time shift, uniform +-[300,3000] s).  Omit for the normal run.")
    ap.add_argument("--shift-mode", default="station", choices=["station", "station_hour"],
                    help="station = one constant offset per station (default); station_hour = a "
                         "fresh offset per station per hour (Muenchmeyer's own control)")
    ap.add_argument("--shift-min", type=float, default=300.0)
    ap.add_argument("--shift-max", type=float, default=3000.0)
    # E4 step 4 (2026-09-11, Merlin E3-audit prescription 3): the E2 tolerances were ~5x the
    # MEASURED S-pick residual (0.40-0.56 s) and the depth range was the whole crust.  These two
    # options let the association be run at tolerances matched to the residual and with the depth
    # restricted to the interface band; defaults reproduce E2/E3 exactly.
    ap.add_argument("--zlim", default="", help='e.g. "20,50" to override the area depth range')
    args = ap.parse_args()

    import pyocto
    stas = pd.read_csv(f"{RES}/e2_stations.csv")
    import e2_pick  # noqa
    stations_set = [s.split(".")[1] for s in e2_pick.SETS[args.setname]]
    stas = stas[stas.sta.isin(stations_set)].copy()
    st = pd.DataFrame(dict(id=stas.sta, latitude=stas.lat, longitude=stas.lon,
                           elevation=stas.elev_m))
    picks = load_picks(args.setname, args.model, args.threshold, set(stas.sta),
                       shift_null=args.shift_null, shift_min=args.shift_min,
                       shift_max=args.shift_max, shift_mode=args.shift_mode)
    if args.phases == "S":
        picks = picks[picks.phase == "S"]
    # PyOcto identifies picks by POSITION in the frame it is handed (_convert_picks enumerates),
    # so the index must be contiguous before the call or pick_idx cannot be joined back.
    picks = picks.reset_index(drop=True)
    print(f"[picks] {len(picks)} picks, {picks.station.nunique()} stations, "
          f"{picks.groupby('phase').size().to_dict()}")

    vm = velocity_model(f"{RES}/e2_pyocto_crosson", 1.0, 260.0, 80.0)
    area = dict(AREAS[args.setname])
    if args.zlim:
        area["zlim"] = tuple(float(v) for v in args.zlim.split(","))
        print(f"[area] depth range overridden to {area['zlim']} km")
    n_s = args.n_stations if args.phases == "S" else args.n_stations
    assoc = pyocto.OctoAssociator.from_area(
        lat=area["lat"], lon=area["lon"], zlim=area["zlim"],
        velocity_model=vm, time_before=120.0,
        n_picks=args.n_stations, n_p_picks=0, n_s_picks=n_s, n_p_and_s_picks=0,
        pick_match_tolerance=args.pick_tol, edt_pick_std=args.pick_std,
        min_interevent_time=3.0, n_threads=args.threads,
    )
    assoc.transform_stations(st)
    ev, asg = assoc.associate(picks, st)
    print(f"[assoc] {len(ev)} events, {len(asg)} assigned picks")
    if len(ev) == 0:
        return
    assoc.transform_events(ev)
    # BUGFIX 2026-09-10 (resumed run): pyocto 0.2.0 names the assignment key `event_idx` (the
    # script assumed `event`), and its events frame has (idx, time, x, y, z, picks) +
    # latitude/longitude/depth from transform_events -- no `residual` column (the residuals live
    # per assigned pick), so the RMS residual is aggregated here instead.
    asg = asg.rename(columns={"event_idx": "event"})
    cnt = asg.groupby("event").agg(
        n_picks=("pick_idx", "count"),
        n_stations=("station", "nunique"),
        n_S=("phase", lambda s: int((s == "S").sum())),
        n_P=("phase", lambda s: int((s == "P").sum())),
        mean_prob=("probability", "mean"),
        rms_resid=("residual", lambda r: float(np.sqrt(np.mean(np.square(r))))))
    ev = ev.join(cnt, on="idx")
    ev["time_utc"] = pd.to_datetime(ev.time, unit="s")
    # PyOcto's n_picks/n_s_picks are PICK counts; the >= N STATION rule of the plan is enforced here
    n_all = len(ev)
    keep = ev.n_stations >= args.min_stations
    print(f"[filter] {int(keep.sum())} of {n_all} events have >= {args.min_stations} distinct "
          f"stations")
    ev = ev[keep].copy()
    asg = asg[asg.event.isin(set(ev.idx))].copy()
    out = args.out or f"{RES}/e2_events_{args.setname}_{args.model}_t{args.threshold}.csv"
    ev.to_csv(out, index=False)
    asg.to_parquet(out.replace(".csv", "_assign.parquet"), index=False)
    print(f"[write] {out}  ({len(ev)} events, {len(asg)} assigned picks)")
    print(ev[["n_stations", "n_picks", "rms_resid", "depth", "mean_prob"]]
          .describe().round(2).to_string())
    print("n_stations histogram:", ev.n_stations.value_counts().sort_index().to_dict())


if __name__ == "__main__":
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    main()
