#!/usr/bin/env python
"""E3 step 6b — F1 readiness recomputed at a 40 km radius (Merlin F8).

F1 reported the fraction of source cells with >= 4 and >= 6 three-component stations within 60 km
and called 60 km "the radius where an LFE is plausibly above noise".  The MEASURED single-station
recall against Bostock is 16 % at < 25 km, 11 % at 25-40 km and 3.2 % at 40-60 km, so the usable
radius is ~40 km and the 60 km number overstates the usable geometry.  This recomputes the same
fractions at 40 km from F1's OWN tables -- no FDSN query -- and validates the method by
reproducing F1's stored n60 first.

Reads results/f1_cell_density.csv (cell lat/lon, region, year, roster, n60/n100/n150) and
results/f1_station_years.csv + f1_stations.csv (per-station-year 3-C availability + coordinates).
Writes results/f1_readiness_40km.csv.  F1's own report is NOT modified.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
R_KM_DEG = 111.195


def main():
    cells = pd.read_csv(f"{RES}/f1_cell_density.csv")
    sy = pd.read_csv(f"{RES}/f1_station_years.csv")
    st = pd.read_csv(f"{RES}/f1_stations.csv")[["net", "sta", "lat", "lon", "klass",
                                                "is_borehole_any"]]
    sy = sy.merge(st, on=["net", "sta"], how="left", suffixes=("", "_s"))
    sy = sy[sy.has_3c.astype(bool) & sy.lat.notna()]
    print(f"[in] {len(cells)} cell-year-roster rows ({cells.cell.nunique()} cells, "
          f"{cells.year.min()}-{cells.year.max()}, rosters {sorted(cells.roster.unique())}); "
          f"{len(sy)} station-years with 3-C data")

    out = []
    for (yr, roster), g in cells.groupby(["year", "roster"]):
        s = sy[sy.year == yr]
        if roster == "fleet":
            s = s[s.roster == "project"]
        if len(s) == 0:
            continue
        cl = g.cell_lat.to_numpy(); cn = g.cell_lon.to_numpy()
        sl = s.lat.to_numpy(); sn = s.lon.to_numpy()
        cos0 = np.cos(np.radians(cl.mean()))
        d = np.hypot((sn[None, :] - cn[:, None]) * R_KM_DEG * cos0,
                     (sl[None, :] - cl[:, None]) * R_KM_DEG)
        bb = (s.klass.to_numpy() == "broadband")
        r = g.copy()
        r["n40"] = (d <= 40).sum(1)
        r["n60_recomputed"] = (d <= 60).sum(1)
        r["n40_bb"] = (d[:, bb] <= 40).sum(1)
        out.append(r)
    C = pd.concat(out, ignore_index=True)
    # ---- validation: does the recomputed n60 reproduce F1's stored n60?
    for roster, g in C.groupby("roster"):
        dd = g.n60_recomputed - g.n60
        print(f"[check] roster={roster}: recomputed n60 - F1 n60: median {dd.median():+.1f}, "
              f"mean {dd.mean():+.2f}, exact match {float((dd == 0).mean()):.3f}, "
              f"|diff|<=1 {float((dd.abs() <= 1).mean()):.3f}")

    rows = []
    for (reg, yr, roster), g in C.groupby(["region", "year", "roster"]):
        rows.append(dict(region=reg, year=yr, roster=roster, n_cells=len(g),
                         med_n40=float(g.n40.median()), med_n60=float(g.n60.median()),
                         frac_ge3_40=float((g.n40 >= 3).mean()),
                         frac_ge4_40=float((g.n40 >= 4).mean()),
                         frac_ge6_40=float((g.n40 >= 6).mean()),
                         frac_ge4_40_bb=float((g.n40_bb >= 4).mean()),
                         frac_ge4_60=float((g.n60 >= 4).mean()),
                         frac_ge6_60=float((g.n60 >= 6).mean())))
    R = pd.DataFrame(rows).sort_values(["roster", "region", "year"])
    R.to_csv(f"{RES}/f1_readiness_40km.csv", index=False)

    print("\n[readiness, project fleet roster: fraction of cells with >= N 3-C stations]")
    f = R[R.roster == "fleet"]
    for yr in (2012, 2018, 2024):
        sub = f[f.year == yr]
        if not len(sub):
            continue
        print(f"  {yr}: " + "  ".join(
            f"{r.region.split('_')[0]} ge4@40 {r.frac_ge4_40:.2f}/ge4@60 {r.frac_ge4_60:.2f} | "
            f"ge6@40 {r.frac_ge6_40:.2f}/ge6@60 {r.frac_ge6_60:.2f}"
            for r in sub.itertuples()))
    print("\n[all regions pooled, project fleet]")
    p = f.groupby("year").apply(
        lambda g: pd.Series(dict(
            cells=int(g.n_cells.sum()),
            frac_ge4_40=float((g.frac_ge4_40 * g.n_cells).sum() / g.n_cells.sum()),
            frac_ge4_60=float((g.frac_ge4_60 * g.n_cells).sum() / g.n_cells.sum()),
            frac_ge6_40=float((g.frac_ge6_40 * g.n_cells).sum() / g.n_cells.sum()),
            frac_ge6_60=float((g.frac_ge6_60 * g.n_cells).sum() / g.n_cells.sum()))),
        include_groups=False)
    print(p.round(3).to_string())
    print(f"\n[write] {RES}/f1_readiness_40km.csv")


if __name__ == "__main__":
    main()
