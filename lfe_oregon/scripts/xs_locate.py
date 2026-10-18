#!/usr/bin/env python
"""Family templates -> P/S picks -> Hyp2000 (Thomas & Bostock 2015 sec.2-3).

  "the remaining LFE families are located using HYPOINVERSE (Klein, 2002) with two different
   velocity models. The first is the O0 model ... PNSN ... the second is a local 1D model employed
   by Williams et al. (2011)"  ;  "most locations have poor depth control so we fix their depths to
   the plate interface. Of the remaining four ... depths ranging from 31 to 34 km"  ;
  "we were unable to locate five families because the templates were not of sufficient quality to
   accurately pick arrivals".

Template = the linear stack of the family's final 3-day detections (xs_grow last iteration) at every
pool station (--stack-from det46 restacks the 46-day detections instead: column).
Picking, Hyp2000 I/O and parsing are lfe_array_discovery/scripts/ab_hypo.py's own functions
(make_picks, write_hyp_inputs, run_hyp), imported; only station tables/models are Oregon's.

NEEDED / SILENT (see METHOD_FROM_PAPER.md):
  * the O0 and Williams et al. (2011) layer tables are NOT on disk -> MODELS holds only what exists
    (GSC, Crosson 1976 'CRO', from ab_hypo); add 'O0' / 'WIL' rows when supplied;
  * "poor depth control" criterion SILENT -> both a free-depth and a fixed-to-interface solution
    are written for every family, with erz as the column a reader ranks on;
  * interface model for the fixed depth SILENT (Fig.1 contours are McCrory 2004) -> McCrory et al.
    (2012) grid on disk (data/slab_models/mccrory2012), evaluated at the free epicentre;
  * station elevations are set to 0 m (the station csv has none).

  python scripts/xs_locate.py --grow results/grow --scan results/scan46 --out results/hypo
-> <out>/stacks_loc.npz, picks.csv, hyp_<model>_<free|fixed>.csv, loc_summary.csv
"""
from __future__ import annotations

import argparse
import os
import re
import glob
import shutil

import numpy as np
import pandas as pd

import xs_common as C

A = C.patch_ab_common()
import ab_hypo as H                                                          # noqa: E402
import xs_grow as XG                                                         # noqa: E402

for r in C.STATIONS.itertuples():                   # Hyp2000 station codes are the bare STA
    A.COORD[r.sta] = (float(r.lat), float(r.lon), 0.0)
    A.NET[r.sta] = r.net


def iface_depth(lat, lon):
    g = pd.read_csv(f"{C.PROJ}/data/slab_models/mccrory2012/Modeled_slab_XYZ_Grid_GCS_WGS84.txt")
    d = np.hypot((g.Longitude - lon) * np.cos(np.radians(lat)), g.Latitude - lat)
    return float(-g.Depth_km.iloc[int(np.argmin(d))])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--grow", required=True)
    ap.add_argument("--scan", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--stack-from", default="det3", choices=["det3", "det46"])
    ap.add_argument("--models", default="GSC,CRO")
    ap.add_argument("--vpvs", type=float, default=1.73)
    ap.add_argument("--rmin", type=float, default=1.5)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    fam = pd.read_csv(f"{a.scan}/families.csv")
    fam = fam[fam.kept] if "kept" in fam else fam
    groups = set(fam.grp.astype(int))
    if a.stack_from == "det3":
        it = max(int(re.findall(r"det_it(\d+)", p)[0]) for p in glob.glob(f"{a.grow}/det_it*.parquet"))
        # each family's detections from the iteration it converged at (its final template's set)
        st = pd.read_csv(f"{a.grow}/state_it{it}.csv").set_index("grp")
        parts = []
        for g in groups:
            k = int(st.at[g, "it_conv"]) if st.at[g, "it_conv"] > 0 else it
            d = pd.read_parquet(f"{a.grow}/det_it{k}.parquet")
            parts.append(d[d.grp == g])
        det = pd.concat(parts, ignore_index=True)
        dates = C.DETECT_DAYS
    else:
        det = pd.read_parquet(f"{a.scan}/det46_families.parquet")
        det = det[det.grp.isin(groups)]
        dates = C.ALL_DAYS
    dc = A.DayCache(XG.POOL)
    acc = XG.restack(dc, det, XG.POOL, dates)
    ST = {}
    for (g, s), v in acc.items():
        if v["n"] >= 2:
            ST.setdefault(g, {})[s.split(".")[1]] = v["lin"] / v["n"]
    np.savez_compressed(f"{a.out}/stacks_loc.npz", **{f"{g}_{s}": v.astype(np.float32)
                                                      for g in ST for s, v in ST[g].items()})
    P = H.make_picks(ST)
    P.to_csv(f"{a.out}/picks.csv", index=False)
    stas = sorted({s for g in ST for s in ST[g]})
    rows = []
    for model in a.models.split(","):
        wd = f"{a.out}/wd_{model}_free"
        shutil.rmtree(wd, ignore_errors=True)
        H.write_hyp_inputs(wd, P, stas, a.vpvs, a.rmin, model=model)
        Hf, _ = H.run_hyp(wd)
        Hf.to_csv(f"{a.out}/hyp_{model}_free.csv", index=False)
        fx = []
        for r in Hf.itertuples():                  # fixed-to-interface solution per family
            z = iface_depth(r.lat, r.lon)
            wd2 = f"{a.out}/wd_{model}_fixed/{r.grp}"
            shutil.rmtree(wd2, ignore_errors=True)
            H.write_hyp_inputs(wd2, P[P.grp == r.grp], stas, a.vpvs, a.rmin, model=model)
            s = open(f"{wd2}/run.hyp").read().replace("ZTR 35", f"ZTR {z:.1f} T")
            open(f"{wd2}/run.hyp", "w").write(s)
            h2, _ = H.run_hyp(wd2)
            if len(h2):
                fx.append(h2.assign(z_iface=z))
        Hx = pd.concat(fx, ignore_index=True) if fx else pd.DataFrame()
        Hx.to_csv(f"{a.out}/hyp_{model}_fixed.csv", index=False)
        for h, mode in ((Hf, "free"), (Hx, "fixed")):
            if len(h):
                rows.append(h.assign(model=model, depth_mode=mode))
    S = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    npk = P.groupby("grp").apply(lambda d: int((d.get("rS", 0) >= a.rmin).sum() +
                                               (d.get("rP", 0) >= a.rmin).sum()))
    fam["n_picks"] = fam.grp.map(npk).fillna(0).astype(int)
    fam["located"] = fam.grp.isin(S.grp) if len(S) else False
    fam.to_csv(f"{a.out}/loc_summary_families.csv", index=False)
    S.to_csv(f"{a.out}/loc_summary.csv", index=False)
    A.log(f"{len(fam)} families, located (any model/mode) {int(fam.located.sum())}")
    open(f"{a.out}/LOC_DONE", "w").write("ok\n")


if __name__ == "__main__":
    main()
