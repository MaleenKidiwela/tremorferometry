#!/usr/bin/env python
"""P2 step 1 — the Slab2 back-projection grid and its S travel-time table.

Grid: 2-D cells of 3 km over 48.0-49.0 N / -124.5..-123.0 E, depth FIXED on Slab2
(../data/cas_slab2_input_04-18.csv, the same filtering as fault_tomography assemble_res_catalog).
Travel times: Crosson (1976) 1-D, Vp/Vs 1.78 (copied d3_vel1d.py -> p_vel1d.py) plus the a5
per-station S corrections.  NO learned location head anywhere (Merlin F3).

Diagnostics written here, all on DEV data only:
  * the Bostock family positions projected onto the grid, and the travel-time difference between the
    family's own (lat, lon, dep) and its nearest grid cell -- the grid's own discretisation error;
  * the moveout spread across sites, which is what makes back-projection resolve at all.

  python scripts/p2_tt.py
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

OUT = f"{P.RES}/p2"


def main():
    os.makedirs(OUT, exist_ok=True)
    g, shape, lats, lons = P.build_grid()
    print(f"[grid] {len(g)} cells, shape {shape}, depth {g.dep.min():.1f}-{g.dep.max():.1f} km")
    TT = P.tt_table(g, P.S15)
    np.savez_compressed(f"{OUT}/grid_tt.npz", cell=g.cell.values, lat=g.lat.values,
                        lon=g.lon.values, dep=g.dep.values, tt=TT,
                        stas=np.array(P.S15), shape=np.array(shape),
                        lats=lats, lons=lons)
    g.to_csv(f"{OUT}/grid_cells.csv", index=False)

    FAM = P.families()
    FT = P.fam_tt()
    rows = []
    for b in FAM.index:
        d = P.epi_km(FAM.lat[b], FAM.lon[b], g.lat.values, g.lon.values)
        c = int(np.argmin(d))
        dt = np.array([FT[(int(b), s)] - TT[c, k] for k, s in enumerate(P.S15)])
        rows.append(dict(bkey=int(b), lat=FAM.lat[b], lon=FAM.lon[b], dep=FAM.dep[b],
                         cell=c, cell_dist_km=float(d[c]), cell_dep=float(g.dep[c]),
                         dtt_med=float(np.median(dt)), dtt_rms=float(np.sqrt((dt ** 2).mean())),
                         dtt_demeaned_rms=float(np.sqrt(((dt - dt.mean()) ** 2).mean())),
                         inside=bool(P.GRID_LAT[0] <= FAM.lat[b] <= P.GRID_LAT[1] and
                                     P.GRID_LON[0] <= FAM.lon[b] <= P.GRID_LON[1])))
    F = pd.DataFrame(rows)
    F.to_csv(f"{OUT}/family_on_grid.csv", index=False)
    print(f"[fam] {F.inside.sum()}/{len(F)} of the 130 families inside the grid box; "
          f"cell distance median {F.cell_dist_km.median():.2f} km")
    print(f"[fam] |family tt - nearest-cell tt| after removing the common origin shift: "
          f"median RMS {F.dtt_demeaned_rms.median():.3f} s, p90 {F.dtt_demeaned_rms.quantile(.9):.3f} s")
    print(f"[fam] Slab2 depth at the family cell vs Bostock depth: median "
          f"{np.median(F.cell_dep - F.dep):+.1f} km")

    # moveout spread: the differential information the back-projection actually uses
    sm = np.array([P.SITES.index(P.SITE_OF[s]) for s in P.S15])
    B = np.zeros((len(g), len(P.SITES)))
    for j in range(len(P.SITES)):
        B[:, j] = TT[:, sm == j].mean(1)
    dmo = B - B.mean(1, keepdims=True)
    print(f"[moveout] differential moveout across the 9 sites: std over cells "
          f"{dmo.std(0).round(2)} s")
    json.dump(dict(ncell=len(g), shape=list(shape), grid_km=P.GRID_KM,
                   lat=list(P.GRID_LAT), lon=list(P.GRID_LON), stations=P.S15,
                   sites=P.SITES,
                   fam_inside=int(F.inside.sum()),
                   fam_cell_dist_km_median=float(F.cell_dist_km.median()),
                   dtt_demeaned_rms_median=float(F.dtt_demeaned_rms.median())),
              open(f"{OUT}/grid_meta.json", "w"), indent=1)


if __name__ == "__main__":
    main()
