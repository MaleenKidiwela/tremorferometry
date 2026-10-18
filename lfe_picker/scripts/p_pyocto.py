#!/usr/bin/env python
"""SECONDARY association comparison: PyOcto on the learned S picks.

The primary association is the explicit physics back-projection (p2_backproj.py).  This is the
control Merlin asked for: the same per-station picks handed to a general-purpose associator with
the same Crosson 1-D velocities, so the network numbers are not an artefact of our own stacker.
S-only association, which PyOcto supports but which is weakly constrained; treat the output as a
cross-check on the ORDER of the detection count, not as a competing catalog.

  python scripts/p_pyocto.py --tag pnet_w05 --year 2012 --rate 3000 --doys 251,255
"""
from __future__ import annotations
import os, sys, json, argparse, tempfile
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

OUT = f"{P.RES}/p1"


def picks_from_scan(tag, year, rate, doys=None):
    z = np.load(f"{OUT}/scan_{tag}_{year}.npz")
    days = pd.read_csv(f"{OUT}/scan_{tag}_{year}_days.csv")
    if doys:
        days = days[days.doy.isin(doys)]
    rows = []
    for r in days.itertuples():
        k = f"{r.sta}_{r.doy}"
        if f"i_{k}" not in z:
            continue
        idx, val = z[f"i_{k}"], z[f"v_{k}"]
        n = max(1, int(round(rate * r.n / P.NDAY)))
        if len(val) == 0:
            continue
        sel = np.argsort(val)[::-1][:n]
        rows.append(pd.DataFrame(dict(station=r.sta, phase="S",
                                      time=r.t0 + idx[sel] / P.FS, probability=val[sel])))
    return pd.concat(rows, ignore_index=True).sort_values("time").reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--year", type=int, required=True)
    ap.add_argument("--rate", type=float, default=3000.0, help="picks per station-day fed in")
    ap.add_argument("--doys", default="")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    assert a.year != 2013, "2013 is sealed"
    import pyocto
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import p_vel1d as V

    doys = [int(x) for x in a.doys.split(",")] if a.doys else None
    pk = picks_from_scan(a.tag, a.year, a.rate, doys)
    ST = P.station_table().loc[P.S15].reset_index()
    stations = pd.DataFrame(dict(id=ST.sta, latitude=ST.lat, longitude=ST.lon,
                                 elevation=ST["elev_m"].astype(float)))

    tmp = tempfile.mkdtemp()
    depths = np.arange(0.0, 70.0, 1.0)
    model = pd.DataFrame(dict(depth=depths,
                              vp=np.interp(depths, V.ZTOP, V.VP, right=V.VP[-1]),
                              vs=np.interp(depths, V.ZTOP, V.VP, right=V.VP[-1]) / V.VPVS_DEFAULT))
    path = f"{tmp}/crosson"
    pyocto.VelocityModel1D.create_model(model, 1.0, 200.0, 70.0, path)
    vm = pyocto.VelocityModel1D(path, tolerance=1.0, association_cutoff_distance=150.0)
    assoc = pyocto.OctoAssociator.from_area(
        lat=P.GRID_LAT, lon=P.GRID_LON, zlim=(20.0, 50.0), velocity_model=vm,
        time_before=120.0, n_picks=4, n_p_picks=0, n_s_picks=4, n_p_and_s_picks=0,
        pick_match_tolerance=1.5)
    assoc.transform_stations(stations)
    ev, ap_ = assoc.associate(pk, stations)
    if len(ev):
        assoc.transform_events(ev)
    out = a.out or f"{OUT}/pyocto_{a.tag}_{a.year}.parquet"
    ev.to_parquet(out, index=False)
    ap_.to_parquet(out.replace(".parquet", "_picks.parquet"), index=False)
    print(f"[pyocto] {len(pk)} picks -> {len(ev)} events over "
          f"{pk.time.pipe(lambda t: (t.max()-t.min())/86400):.1f} days -> {out}")


if __name__ == "__main__":
    main()
