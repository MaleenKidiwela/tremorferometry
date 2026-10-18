#!/usr/bin/env python
"""EXTENSION scorecard: literal vs extension (75/100/150 km pools), real beside every null, every
family keeps its row, every criterion a column (rank, don't discard).

Per family: detections, Hyp2000 depth / ERH / ERZ, station-bootstrap radius, hypoDD depth and
retention -- literal and each radius side by side -- plus the PERMANENT (non-temporary-network)
stations that carry the family, ranked by their stack SNR, for the later extension to 2026.

  python scripts/ab_ext_score.py  ->  results/ab_royer/ext/score/<SUB>_families.csv,
                                      <SUB>_permanent_carriers.csv, <SUB>_summary.csv, SUMMARY.md
"""
from __future__ import annotations

import glob
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

RADII = [75, 100, 150]
HEAD = "loc_headline.csv"
COLS = ["dep", "erh", "erz", "boot_km", "boot_dz", "dd_dep", "dd_kept", "num", "rms"]


def fam(grow, hypo):
    gd = sorted(glob.glob(f"{grow}/det_it*.parquet"),
                key=lambda p: int(p.split("it")[-1].split(".")[0]))
    if not gd:
        return pd.DataFrame()
    D = pd.read_parquet(gd[-1])
    F = D.groupby("grp").agg(n_det=("cc", "size"), n_days=("date", "nunique"),
                             nmad_med=("nmad", "median")).reset_index()
    p = f"{hypo}/{HEAD}"
    if os.path.exists(p):
        H = pd.read_csv(p)
        F = F.merge(H[["grp"] + [c for c in COLS if c in H]], on="grp", how="left")
    return F


def main():
    out = f"{A.EXT_DIR}/score"
    os.makedirs(out, exist_ok=True)
    md = ["# Royer & Bostock (2014) -- EXTENSION vs literal\n",
          "Extension = literal autocorrelation detections and initial templates, stations added "
          "ONLY at growth and location (public networks within R km).  Columns, never gates.\n"]
    for sub in ("NW", "SW", "NVI"):
        pops = [("real", f"{A.RES}/{sub}", "")] + [
            (os.path.basename(d), d, os.path.basename(d))
            for d in sorted(glob.glob(f"{A.RES}/{sub}/null[0-9]"))]
        allrows, summ, carriers = [], [], []
        for pop, lit, nk in pops:
            base = fam(f"{lit}/grow", f"{lit}/hypo").add_prefix("lit_").rename(
                columns={"lit_grp": "grp"})
            M = base
            for R in RADII:
                root = f"{A.EXT_DIR}/R{R}/{sub}" + (f"/{nk}" if nk else "")
                E = fam(f"{root}/grow", f"{root}/hypo").add_prefix(f"R{R}_").rename(
                    columns={f"R{R}_grp": "grp"})
                if len(E):
                    M = M.merge(E, on="grp", how="outer") if len(M) else E
                sp = f"{root}/hypo/station_snr.csv"
                if os.path.exists(sp):
                    Q = pd.read_csv(sp)
                    Q = Q[Q.permanent]
                    Q["rank"] = Q.groupby("grp").snr.rank(ascending=False)
                    Q["pop"], Q["radius_km"] = pop, R
                    carriers.append(Q)
                for c in ["dep", "dd_dep", "boot_km", "erh", "erz"]:
                    if f"R{R}_{c}" in M:
                        summ.append(dict(pop=pop, radius_km=R, criterion=c,
                                         ext_median=float(M[f"R{R}_{c}"].median()),
                                         lit_median=float(M[f"lit_{c}"].median())
                                         if f"lit_{c}" in M else np.nan,
                                         n_ext=int(M[f"R{R}_{c}"].notna().sum()),
                                         n_lit=int(M[f"lit_{c}"].notna().sum())
                                         if f"lit_{c}" in M else 0))
            if len(M):
                M.insert(0, "pop", pop)
                allrows.append(M)
        if allrows:
            pd.concat(allrows, ignore_index=True).to_csv(f"{out}/{sub}_families.csv",
                                                         index=False)
        if carriers:
            pd.concat(carriers, ignore_index=True).to_csv(
                f"{out}/{sub}_permanent_carriers.csv", index=False)
        S = pd.DataFrame(summ)
        S.to_csv(f"{out}/{sub}_summary.csv", index=False)
        md.append(f"\n## {sub}\n\n```\n{S.to_string(index=False) if len(S) else '(none)'}\n```\n")
    open(f"{out}/SUMMARY.md", "w").write("".join(md))
    A.log(f"[ext score] wrote {out}/SUMMARY.md")


if __name__ == "__main__":
    main()
