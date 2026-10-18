#!/usr/bin/env python
"""E2 step 3b — NonLinLoc location of the PyOcto-associated LFE events.

Velocity model: the casc1.6 3-D block is on disk but only as a 1.4 GB zip that expands to a 17.7 GB
ASCII table (E1 report §5); converting it into NonLinLoc 3-D grids for a pilot is not worth the I/O,
and Muenchmeyer's released LFE workflow used a 1-D model with NonLinLoc.  This step therefore uses
the same Crosson (1976) 1-D model as the association, as LAYER cards in a 2-D (GRID2D) model grid.
S-weighted: S picks get a 0.4 s Gaussian error, P picks 1.5 s (Lin measured sigma_P = 4.2 s for the
LFE P phase and discarded it entirely, so P is kept only as a weak constraint).

Writes `<events>.nll.csv` with lat/lon/depth, RMS, ERH/ERZ from the NonLinLoc statistics, plus the
Slab2 interface depth under each epicentre.
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
VP = [5.40, 6.38, 6.59, 6.73, 6.86, 7.00, 7.90]
ZTOP = [0.0, 4.0, 9.0, 16.0, 20.0, 25.0, 32.0]
VPVS = 1.78
LAT0, LON0 = 48.6, -123.8


def write_control(root, stas, obs_glob, nx=111, ny=111, nz=41, dloc=2.0, ndist=261,
                  ztop=-3.0, zmax=80.0):
    """BUGFIX 2026-09-10 (resumed run), two NonLinLoc 6.04 constraints found the hard way:
      * long ABSOLUTE paths in GTFILES corrupt the wave-type field (the model grid name comes out as
        `layer.P<binary garbage>.mod`), so every path in the control files is RELATIVE and the
        binaries are run with cwd=root;
      * the model/time grid must START ABOVE SEA LEVEL (ztop = -3 km): with a 0 km top every
        GTSRCE at positive elevation gave `Source point is not inside model grid` and no travel-time
        grid was written at all;
      * with `LOCGRID ... NO_SAVE` this build locates every event ("200 events located") but writes
        NO hypocentre file at all and reports "0 locations completed"; the flag has to be SAVE.  The
        location grid is therefore coarsened to 2 km (111 x 111 x 41 = 2.0 MB/event instead of
        15.8 MB) and main() deletes the grid/scatter files after parsing the summary.
    """
    os.makedirs(f"{root}/model", exist_ok=True)
    os.makedirs(f"{root}/time", exist_ok=True)
    os.makedirs(f"{root}/loc", exist_ok=True)
    lay = "\n".join(f"LAYER {z:6.2f} {vp:5.2f} 0.0 {vp/VPVS:5.3f} 0.0 2.7 0.0"
                    for z, vp in zip(ZTOP, VP))
    src = "\n".join(f"GTSRCE {r.sta} LATLON {r.lat:.5f} {r.lon:.5f} 0.0 {r.elev_m/1000.0:.4f}"
                    for r in stas.itertuples())
    common = f"""CONTROL 1 54321
TRANS SIMPLE {LAT0} {LON0} 0.0
"""
    # model/time grid: 1 km from ztop to zmax (independent of the coarser LOCATION grid --
    # if the time grid is shallower than LOCGRID, NLLoc silently skips every event:
    # "200 events read, 0 events located" with no error line)
    nzm = int(zmax - ztop) + 1
    for phase in ("P", "S"):
        with open(f"{root}/vel2grid_{phase}.in", "w") as f:
            f.write(common + f"""VGOUT model/layer
VGTYPE {phase}
VGGRID 2 {ndist} {nzm} 0.0 0.0 {ztop} 1.0 1.0 1.0 SLOW_LEN
{lay}
""")
        with open(f"{root}/grid2time_{phase}.in", "w") as f:
            f.write(common + f"""VGOUT model/layer
VGTYPE {phase}
VGGRID 2 {ndist} {nzm} 0.0 0.0 {ztop} 1.0 1.0 1.0 SLOW_LEN
{lay}
GTFILES model/layer time/layer {phase}
GTMODE GRID2D ANGLES_NO
{src}
GT_PLFD 1.0e-3 0
""")
    with open(f"{root}/nlloc.in", "w") as f:
        f.write(common + f"""LOCSIG e2 LFE pilot
LOCCOM LFEDetect + PyOcto + NonLinLoc (Crosson 1-D)
LOCFILES {obs_glob} NLLOC_OBS time/layer loc/lfe
LOCHYPOUT SAVE_NLLOC_SUM
LOCSEARCH OCT 20 20 8 0.01 50000 5000 0 0
LOCGRID {nx} {ny} {nz} -110.0 -110.0 0.0 {dloc} {dloc} {dloc} PROB_DENSITY SAVE
LOCMETH EDT_OT_WT 9999.0 3 -1 -1 -1 -1 -1 0
LOCGAU 0.2 0.0
LOCGAU2 0.02 0.05 2.0
LOCPHASEID P P Pn Pg
LOCPHASEID S S Sn Sg
LOCQUAL2ERR 0.1 0.5 1.0 2.0 99999.9
LOCANGLES ANGLES_NO 5
LOCPHSTAT 9999.0 -1 9999.0 1.0 1.0
""")
    return f"{root}/nlloc.in"


def build_grids(root):
    if os.path.exists(f"{root}/time/layer.S.PGC.time.hdr"):
        return
    for phase in ("P", "S"):
        for prog, ctl in (("Vel2Grid", f"vel2grid_{phase}.in"),
                          ("Grid2Time", f"grid2time_{phase}.in")):
            r = subprocess.run([f"{BIN}/{prog}", ctl], capture_output=True, text=True,
                               errors="replace", cwd=root)  # BUGFIX: NLL banners are not UTF-8
            tail = (r.stdout + r.stderr).strip().split("\n")[-3:]
            print(f"[{prog} {phase}] rc={r.returncode} " + " | ".join(tail))
            if r.returncode != 0 or "ERROR" in (r.stdout + r.stderr):
                bad = [l for l in (r.stdout + r.stderr).split("\n") if "ERROR" in l][:3]
                raise SystemExit(f"{prog} failed: " + " | ".join(bad))


def write_obs(asg, ev, path, err_S=0.4, err_P=1.5, max_events=None):
    ids = ev.idx.to_numpy()
    if max_events and len(ids) > max_events:
        rng = np.random.default_rng(0)
        ids = np.sort(rng.choice(ids, max_events, replace=False))
    g = asg[asg.event.isin(set(ids))].copy()
    g["t"] = pd.to_datetime(g.time, unit="s")
    n = 0
    with open(path, "w") as f:
        for eid, blk in g.groupby("event", sort=True):
            for r in blk.itertuples():
                err = err_S if r.phase == "S" else err_P
                f.write(f"{r.station:<6s} ?    ?    ? {r.phase:<6s} ? "
                        f"{r.t.strftime('%Y%m%d %H%M')} {r.t.second + r.t.microsecond/1e6:9.4f} "
                        f"GAU {err:9.2e} -1.00e+00 -1.00e+00 -1.00e+00 1.00\n")
            f.write("\n")
            n += 1
    return n, ids


def parse_sum(path):
    """Parse a NonLinLoc SAVE_NLLOC_SUM hypocentre file."""
    rows, cur = [], {}
    with open(path) as f:
        for ln in f:
            if ln.startswith("NLLOC "):
                cur = {}
            elif ln.startswith("GEOGRAPHIC "):
                p = ln.split()
                cur.update(year=int(p[2]), month=int(p[3]), day=int(p[4]),
                           hour=int(p[5]), minute=int(p[6]), sec=float(p[7]),
                           lat=float(p[9]), lon=float(p[11]), depth=float(p[13]))
            elif ln.startswith("QUALITY "):
                p = ln.split()
                cur.update(rms=float(p[8]), nphs=int(p[10]), gap=float(p[12]),
                           dmin=float(p[14]))
            elif ln.startswith("STATISTICS "):
                p = ln.split()
                cov = {p[i]: float(p[i + 1]) for i in range(1, len(p) - 1, 2)}
                exx, eyy, ezz = cov.get("CovXX", np.nan), cov.get("YY", np.nan), cov.get("ZZ", np.nan)
                cur.update(err_h_km=float(np.sqrt(max(exx, 0) + max(eyy, 0))),
                           err_z_km=float(np.sqrt(max(ezz, 0))))
            elif ln.startswith("END_NLLOC"):
                if cur:
                    rows.append(cur)
                cur = {}
    df = pd.DataFrame(rows)
    if len(df):
        df["time_utc"] = pd.to_datetime(dict(year=df.year, month=df.month, day=df.day,
                                             hour=df.hour, minute=df.minute,
                                             second=df.sec.astype(int))) \
            + pd.to_timedelta(df.sec % 1.0, unit="s")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", required=True)
    ap.add_argument("--max-events", type=int, default=0)
    ap.add_argument("--root", default="")
    args = ap.parse_args()

    ev = pd.read_csv(args.events)
    asg = pd.read_parquet(args.events.replace(".csv", "_assign.parquet"))
    stas = pd.read_csv(f"{RES}/e2_stations.csv")
    stas = stas[stas.sta.isin(set(asg.station))]
    root = args.root or os.path.abspath(f"{RES}/nll_" + os.path.basename(args.events)[:-4])
    os.makedirs(root, exist_ok=True)
    obs = f"{root}/obs.txt"
    ctl = write_control(root, stas, "obs.txt")   # relative: NLL is run with cwd=root
    build_grids(root)
    n, ids = write_obs(asg, ev, obs, max_events=args.max_events or None)
    print(f"[obs] {n} events written to {obs}")
    r = subprocess.run([f"{BIN}/NLLoc", os.path.basename(ctl)], capture_output=True, text=True,
                       errors="replace", cwd=root)  # BUGFIX: NLL banners are not UTF-8
    print(f"[NLLoc] rc={r.returncode}")
    print("\n".join((r.stdout + r.stderr).strip().split("\n")[-6:]))
    sumf = f"{root}/loc/lfe.sum.grid0.loc.hyp"
    if not os.path.exists(sumf):
        cands = [p for p in os.listdir(f"{root}/loc") if p.endswith(".hyp")]
        print("[warn] no summary file; candidates:", cands[:5])
        raise SystemExit(1)
    loc = parse_sum(sumf)
    print(f"[loc] {len(loc)} located of {n}")
    # keep the hypocentre files, drop the per-event probability grids (2 MB each)
    ndel = 0
    for f in os.listdir(f"{root}/loc"):
        if f.endswith((".buf", ".hdr", ".scat", ".octree")) and ".sum." not in f:
            os.remove(f"{root}/loc/{f}"); ndel += 1
    print(f"[clean] removed {ndel} per-event grid files")
    # attach the PyOcto event id: obs.txt is written in ascending event id and NLLoc preserves that
    # order in the summary file (verified: all events located, no skips).  If any event was skipped
    # the ids are matched on origin time instead.
    if len(loc) == len(ids):
        loc["event_idx"] = ids
    else:
        et = ev.set_index("idx").loc[ids, "time"].to_numpy()
        lt = loc.time_utc.astype("datetime64[ns]").astype("int64").to_numpy() / 1e9
        j = np.argmin(np.abs(lt[:, None] - et[None, :]), axis=1)
        loc["event_idx"] = np.asarray(ids)[j]
        print(f"[warn] {len(ids)-len(loc)} events not located; ids matched on origin time")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import d3_vel1d as v1d
    loc["slab2_km"] = v1d.slab2_depth(loc.lat.values, loc.lon.values)
    out = args.events.replace(".csv", ".nll.csv")
    loc.to_csv(out, index=False)
    print(f"[write] {out}")
    print(loc[["lat", "lon", "depth", "rms", "nphs", "err_h_km", "err_z_km", "slab2_km"]]
          .describe().round(2).to_string())


if __name__ == "__main__":
    main()
