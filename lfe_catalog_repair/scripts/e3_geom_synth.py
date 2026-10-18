#!/usr/bin/env python
"""E3 step 2, config C — the GEOMETRY + S-ONLY term, through the EXACT e2_locate_nll.py control.

The question config C answers: given S-only picks of a known quality on the six vi2010-overlap
stations, how far does the S-only origin-time/depth trade-off move the answer?  It is a synthetic
resolution test, and it is labelled as one: the SOURCES are the 129 published Bostock family
hypocentres and the arrival times are read out of NonLinLoc's OWN travel-time grids, so the model
term is identically zero by construction and the only error is the injected pick noise plus the
geometry.  It is NOT a measure of model error (that is config A, e3_family_relocate.py) nor of
picker error (that is config D, e3_realpick_locate.py).

Design (Merlin, 2026-09-11): sources = all 129 family hypocentres; station sets = every 3-, 4-, 5-
and 6-station subset of PGC/LZB/VGZ/YOUB/SNB/PFB (42 subsets); S picks only, with the same 0.4 s
GAU error card `e2_locate_nll.write_obs` writes; noise levels
  0.05 s  pure geometry + grid-top pinning,
  0.30 s  HEADLINE -- the MEASURED E2 S-pick scatter against the Bostock stack time base,
  0.50 s  bracket -- the measured scatter against Lin-hq (which includes Lin's own location error),
  heavy  85 % N(0, 0.3 s) + 15 % U(-2.5, 2.5) s -- PyOcto's pick-admission window,
with `--reps` realizations each.  Travel times come from `time/layer.S.<STA>.time.buf` in the run's
own NLL root (GRID2D: nx=1, ny = horizontal distance 0..260 km @1 km, nz = depth -3..80 km @1 km).

Usage:
  python lfe_catalog_repair/scripts/e3_geom_synth.py --sigma 0.3 --reps 3 --tag s030
"""
from __future__ import annotations

import argparse
import itertools
import os
import subprocess
import sys

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
BIN = "/home/jovyan/.local/bin"
MD = f"{RES}/external/other/seisbench_lfestacks_cascadia/metadata.csv"
OVERLAP = ["PGC", "LZB", "VGZ", "YOUB", "SNB", "PFB"]
R_KM_DEG = 111.195
LAT0, LON0 = 48.6, -123.8      # the SIMPLE transform origin of e2_locate_nll.write_control
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_locate_nll as NLL  # noqa: E402


def read_time_grid(root, phase, sta):
    """(x_km axis, z_km axis, T[nx, nz]) from a NonLinLoc GRID2D time grid."""
    h = open(f"{root}/time/layer.{phase}.{sta}.time.hdr").read().split("\n")
    p = h[0].split()
    nx, ny, nz = int(p[0]), int(p[1]), int(p[2])
    oy, oz = float(p[4]), float(p[5])
    dy, dz = float(p[7]), float(p[8])
    sp = h[1].split()
    sxyz = (float(sp[1]), float(sp[2]), float(sp[3]))
    T = np.fromfile(f"{root}/time/layer.{phase}.{sta}.time.buf", dtype="<f4")
    T = T.reshape(nx, ny, nz)[0]                   # GRID2D: nx == 1
    return oy + dy * np.arange(ny), oz + dz * np.arange(nz), T.astype(float), sxyz


def interp2(xa, za, T, xq, zq):
    xq = np.clip(xq, xa[0], xa[-1]); zq = np.clip(zq, za[0], za[-1])
    i = np.clip(np.searchsorted(xa, xq) - 1, 0, len(xa) - 2)
    j = np.clip(np.searchsorted(za, zq) - 1, 0, len(za) - 2)
    tx = (xq - xa[i]) / (xa[i + 1] - xa[i]); tz = (zq - za[j]) / (za[j + 1] - za[j])
    return ((1 - tx) * (1 - tz) * T[i, j] + tx * (1 - tz) * T[i + 1, j]
            + (1 - tx) * tz * T[i, j + 1] + tx * tz * T[i + 1, j + 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigma", default="0.3", help="0.05 / 0.3 / 0.5 / heavy")
    ap.add_argument("--reps", type=int, default=2)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-n", type=int, default=3)
    ap.add_argument("--max-n", type=int, default=6)
    args = ap.parse_args()

    root = os.path.abspath(f"{RES}/nll_e3_geom_{args.tag}")
    os.makedirs(root, exist_ok=True)
    stas = pd.read_csv(f"{RES}/e2_stations.csv")
    stas = stas[stas.sta.isin(OVERLAP)]
    NLL.write_control(root, stas, "obs.txt")       # the EXACT e2 control, unmodified
    NLL.build_grids(root)

    grids = {s: read_time_grid(root, "S", s) for s in OVERLAP}
    xs = {s: grids[s][3] for s in OVERLAP}         # station x, y in the SIMPLE transform

    m = pd.read_csv(MD)
    fam = m.groupby("source_family_id").first()[
        ["source_latitude_deg", "source_longitude_deg", "source_depth_km",
         "source_azimuthal_gap_deg"]].reset_index()
    cos0 = np.cos(np.radians(LAT0))
    fam["x"] = (fam.source_longitude_deg - LON0) * R_KM_DEG * cos0
    fam["y"] = (fam.source_latitude_deg - LAT0) * R_KM_DEG
    print(f"[src] {len(fam)} published family hypocentres; depth "
          f"{fam.source_depth_km.min():.1f}..{fam.source_depth_km.max():.1f} km")

    subsets = [c for n in range(args.min_n, args.max_n + 1)
               for c in itertools.combinations(OVERLAP, n)]
    print(f"[sets] {len(subsets)} station subsets of {OVERLAP}")

    rng = np.random.default_rng(args.seed)
    T0 = 1_281_398_400.0                 # arbitrary common origin (2010-08-10 00:00 UTC)
    rows, obs_lines, eid = [], [], 0
    for _, f in fam.iterrows():
        tts = {}
        for s in OVERLAP:
            xa, za, T, (sx, sy, sz) = grids[s]
            d = np.hypot(f.x - sx, f.y - sy)
            tts[s] = float(interp2(xa, za, T, d, f.source_depth_km))
        for sub in subsets:
            for rep in range(args.reps):
                if args.sigma == "heavy":
                    n = np.where(rng.random(len(sub)) < 0.15,
                                 rng.uniform(-2.5, 2.5, len(sub)),
                                 rng.normal(0.0, 0.3, len(sub)))
                else:
                    n = rng.normal(0.0, float(args.sigma), len(sub))
                t0 = T0 + eid * 600.0            # events far apart: NLL treats each obs block alone
                for s, ni in zip(sub, n):
                    tt = pd.Timestamp(t0 + tts[s] + ni, unit="s")
                    obs_lines.append(
                        f"{s:<6s} ?    ?    ? {'S':<6s} ? {tt.strftime('%Y%m%d %H%M')} "
                        f"{tt.second + tt.microsecond/1e6:9.4f} GAU {0.4:9.2e} -1.00e+00 "
                        f"-1.00e+00 -1.00e+00 1.00\n")
                obs_lines.append("\n")
                rows.append(dict(eid=eid, fam_id=float(f.source_family_id), n_sta=len(sub),
                                 subset="+".join(sub), rep=rep,
                                 lat_pub=f.source_latitude_deg, lon_pub=f.source_longitude_deg,
                                 dep_pub=f.source_depth_km, gap_pub=f.source_azimuthal_gap_deg))
                eid += 1
    truth = pd.DataFrame(rows)
    with open(f"{root}/obs.txt", "w") as fh:
        fh.writelines(obs_lines)
    print(f"[obs] {len(truth)} synthetic events written ({args.sigma} s noise, {args.reps} reps)")

    r = subprocess.run([f"{BIN}/NLLoc", "nlloc.in"], capture_output=True, text=True,
                       errors="replace", cwd=root)
    print(f"[NLLoc] rc={r.returncode} " + " | ".join((r.stdout + r.stderr).strip().split("\n")[-3:]))
    loc = NLL.parse_sum(f"{root}/loc/lfe.sum.grid0.loc.hyp")
    print(f"[loc] {len(loc)} of {len(truth)} located")
    for fn in os.listdir(f"{root}/loc"):
        if fn.endswith((".buf", ".hdr", ".scat", ".octree")) and ".sum." not in fn:
            os.remove(f"{root}/loc/{fn}")
    # obs blocks are written in ascending eid and NLLoc preserves that order
    if len(loc) == len(truth):
        loc["eid"] = truth.eid.to_numpy()
    else:
        lt = loc.time_utc.astype("datetime64[ns]").astype("int64").to_numpy() / 1e9
        loc["eid"] = np.round((lt - T0) / 600.0).astype(int)
        print(f"[warn] {len(truth)-len(loc)} events not located; eids recovered from origin time")
    out = truth.merge(loc, on="eid", how="inner")
    out["dh_km"] = np.hypot((out.lon - out.lon_pub) * R_KM_DEG * cos0,
                            (out.lat - out.lat_pub) * R_KM_DEG)
    out["dz_km"] = out.depth - out.dep_pub
    out["sigma"] = args.sigma
    out.to_csv(f"{RES}/e3_geom_synth_{args.tag}.csv", index=False)
    print()
    g = out.groupby("n_sta").agg(n=("dh_km", "size"), dh_med=("dh_km", "median"),
                                 dh_q25=("dh_km", lambda v: v.quantile(.25)),
                                 dh_q75=("dh_km", lambda v: v.quantile(.75)),
                                 dz_med=("dz_km", "median"),
                                 frac_shallow=("depth", lambda v: float((v < 5).mean())),
                                 errh=("err_h_km", "median"), errz=("err_z_km", "median"),
                                 rms=("rms", "median"))
    print(f"[config C, sigma={args.sigma}]")
    print(g.round(2).to_string())
    print(f"\n[write] {RES}/e3_geom_synth_{args.tag}.csv")


if __name__ == "__main__":
    main()
