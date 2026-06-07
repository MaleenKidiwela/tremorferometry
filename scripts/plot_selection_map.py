"""Coverage-SELECTION map in the same cartopy style as plot_family_map.py.

Shows the eligible family POOL (>= SNR floor) greyed out, with the FINAL selected
families (coverage-balanced + top-10% SNR add-back) highlighted as SNR-colored
diamonds -- so you can see what was kept vs. what was available, on a real
coastline/land basemap (matching the family map).

Example:
  python scripts/plot_selection_map.py --station B033 \
      --station-lat 43.2917 --station-lon -123.1245 \
      --summary  data/b033_pnsn_families_100km.summary.csv \
      --selected data/b033_coverage_selection.summary.csv \
      --min-snr 5 --out figures/smoke_b033_coverage_selection.png
"""
from __future__ import annotations

import argparse

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import pandas as pd

REF_STATIONS = {
    "GNW": (47.5641, -122.8250),
    "PGC": (48.6498, -123.4521),
    "NLLB": (49.2271, -123.9882),
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", required=True, help="FULL families summary csv")
    p.add_argument("--selected", required=True,
                   help="final selected subset csv (coverage + top-10%)")
    p.add_argument("--dvv", default=None,
                   help="dv/v csv; if given, colored highlight = families that actually "
                        "produced a dv/v curve, and selected-but-no-curve families are greyed")
    p.add_argument("--station", required=True)
    p.add_argument("--station-lat", type=float, required=True)
    p.add_argument("--station-lon", type=float, required=True)
    p.add_argument("--min-snr", type=float, default=5.0,
                   help="SNR floor defining the greyed eligible pool")
    p.add_argument("--out", required=True)
    p.add_argument("--pad", type=float, default=0.35)
    return p.parse_args()


def idcol(df):
    for c in ("family_id", "family", "id"):
        if c in df.columns:
            return c
    raise SystemExit("no family id column found")


def main():
    args = parse_args()
    full = pd.read_csv(args.summary)
    sel = pd.read_csv(args.selected)
    fid = idcol(full)

    pool = full[full["snr"] >= args.min_snr].copy()      # eligible pool
    sel_ids = set(sel[idcol(sel)])
    if args.dvv:
        prod_ids = set(pd.read_csv(args.dvv)["patch"].unique())   # families w/ a dv/v curve
    else:
        prod_ids = sel_ids
    hi = full[full[fid].isin(prod_ids)].copy()                    # colored: produced dv/v
    sel_nodvv = full[full[fid].isin(sel_ids - prod_ids)].copy()   # grey+edge: selected, no curve
    pool_only = pool[~pool[fid].isin(sel_ids | prod_ids)]         # light grey: rest of pool
    print(f"pool(SNR>={args.min_snr}):{len(pool)} selected:{len(sel_ids)} "
          f"produced-dv/v:{len(hi)} selected-no-dv/v:{len(sel_nodvv)}")

    lons = list(hi["lon"]) + [args.station_lon]
    lats = list(hi["lat"]) + [args.station_lat]
    ext = [min(lons) - args.pad, max(lons) + args.pad,
           min(lats) - args.pad, max(lats) + args.pad]

    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(10, 9))
    ax = plt.axes(projection=proj)
    ax.set_extent(ext, crs=proj)

    ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#cfe2f3")
    ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#f5ead6")
    ax.add_feature(cfeature.LAKES.with_scale("10m"), facecolor="#cfe2f3")
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), lw=0.6, edgecolor="#555")
    ax.add_feature(cfeature.BORDERS.with_scale("10m"), lw=0.5,
                   edgecolor="#888", linestyle=":")

    # light grey: eligible pool that was neither selected nor produced a curve
    ax.scatter(pool_only["lon"], pool_only["lat"], s=12, c="0.78", marker="o",
               edgecolor="none", alpha=0.5, zorder=4, transform=proj,
               label=f"eligible pool (n={len(pool)})")

    # grey with edge: SELECTED but produced NO dv/v curve
    if len(sel_nodvv):
        ax.scatter(sel_nodvv["lon"], sel_nodvv["lat"], s=34, c="0.55", marker="o",
                   edgecolor="k", linewidth=0.4, alpha=0.85, zorder=4.5, transform=proj,
                   label=f"selected, no dv/v (n={len(sel_nodvv)})")

    # colored by SNR (diamonds): families that actually produced a dv/v curve
    nmem = hi["n_members"].clip(upper=20) if "n_members" in hi.columns else 10
    sc = ax.scatter(
        hi["lon"], hi["lat"], c=hi["snr"], cmap="viridis",
        s=25 + 6 * nmem, marker="D", edgecolor="k", linewidth=0.4, zorder=5,
        transform=proj, label=f"produced dv/v (n={len(hi)})",
    )
    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("template SNR")

    ax.plot(args.station_lon, args.station_lat, marker="^", color="red",
            markersize=15, markeredgecolor="k", linestyle="none", zorder=6,
            transform=proj, label=f"{args.station} (target)")
    ax.text(args.station_lon + 0.03, args.station_lat + 0.02, args.station,
            fontsize=10, fontweight="bold", zorder=7, transform=proj)
    for name, (lat, lon) in REF_STATIONS.items():
        if name == args.station:
            continue
        if not (ext[0] <= lon <= ext[1] and ext[2] <= lat <= ext[3]):
            continue
        ax.plot(lon, lat, marker="^", color="#666", markersize=9,
                markeredgecolor="k", linestyle="none", zorder=6, transform=proj)
        ax.text(lon + 0.03, lat + 0.02, name, fontsize=8, zorder=7, transform=proj,
                clip_on=True)

    gl = ax.gridlines(draw_labels=True, lw=0.4, color="gray", alpha=0.4,
                      linestyle="--")
    gl.top_labels = False
    gl.right_labels = False

    ax.set_title(
        f"{args.station}: {len(hi)} families produced dv/v of {len(sel_ids)} selected "
        f"(eligible pool {len(pool)} @ SNR>={args.min_snr:.0f})\n"
        "colored = produced dv/v curve; grey+edge = selected but no curve; light grey = pool",
        fontsize=10,
    )
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
