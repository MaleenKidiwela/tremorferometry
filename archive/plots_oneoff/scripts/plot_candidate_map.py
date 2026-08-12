"""Map the raw stage-1 candidate density (envelope-peak detections) around a
target station -- shows the actual tremor footprint that families are drawn from.

Example:
  python scripts/plot_candidate_map.py --station HDW \
      --station-lat 47.649 --station-lon -123.053 \
      --candidates data/hdw_pnsn_candidates.parquet \
      --out figures/smoke_hdw_candidate_map.png
"""
from __future__ import annotations

import argparse

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm

REF_STATIONS = {
    "GNW": (47.5641, -122.8250),
    "PGC": (48.6498, -123.4521),
    "NLLB": (49.2271, -123.9882),
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--candidates", required=True)
    p.add_argument("--station", required=True)
    p.add_argument("--station-lat", type=float, required=True)
    p.add_argument("--station-lon", type=float, required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--gridsize", type=int, default=90)
    p.add_argument("--extent", nargs=4, type=float, default=None,
                   help="lon_min lon_max lat_min lat_max -- force extent (e.g. the "
                        "PNSN bbox) instead of auto-cropping to data")
    return p.parse_args()


def main():
    args = parse_args()
    c = pd.read_parquet(args.candidates, columns=["lat", "lon"])
    print(f"{len(c):,} candidates")

    if args.extent:
        ext = list(args.extent)
        bbox = args.extent  # draw the box edges
    else:
        pad = 0.25
        ext = [c["lon"].min() - pad, c["lon"].max() + pad,
               c["lat"].min() - pad, c["lat"].max() + pad]
        bbox = None

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

    hb = ax.hexbin(c["lon"], c["lat"], gridsize=args.gridsize, cmap="magma_r",
                   norm=LogNorm(vmin=1), mincnt=1, zorder=4, transform=proj,
                   linewidths=0.1)
    cb = fig.colorbar(hb, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("candidate count per cell (log)")

    if bbox is not None:
        lon0, lon1, lat0, lat1 = bbox
        ax.plot([lon0, lon1, lon1, lon0, lon0], [lat0, lat0, lat1, lat1, lat0],
                color="red", lw=1.4, ls="--", zorder=5, transform=proj,
                label="PNSN tremor bbox")
        # mark the actual tremor floor (southernmost candidate)
        tf = c["lat"].min()
        ax.plot([lon0, lon1], [tf, tf], color="darkorange", lw=1.2, ls=":",
                zorder=5, transform=proj)
        ax.text(lon0 + 0.05, tf + 0.02,
                f"tremor floor {c['lat'].min():.2f}N (no PNSN tremor south of here)",
                fontsize=7, color="darkorange", zorder=7, transform=proj)

    ax.plot(args.station_lon, args.station_lat, marker="^", color="red",
            markersize=15, markeredgecolor="k", linestyle="none", zorder=6,
            transform=proj, label=f"{args.station} (target)")
    ax.text(args.station_lon + 0.03, args.station_lat + 0.02, args.station,
            fontsize=10, fontweight="bold", zorder=7, transform=proj)
    for name, (lat, lon) in REF_STATIONS.items():
        if name == args.station:
            continue
        ax.plot(lon, lat, marker="^", color="#3a3a3a", markersize=9,
                markeredgecolor="w", linestyle="none", zorder=6, transform=proj)
        ax.text(lon + 0.03, lat + 0.02, name, fontsize=8, zorder=7,
                color="#222", transform=proj)

    gl = ax.gridlines(draw_labels=True, lw=0.4, color="gray", alpha=0.4,
                      linestyle="--")
    gl.top_labels = False
    gl.right_labels = False

    ax.set_title(
        f"{args.station}: {len(c):,} stage-1 LFE candidates "
        "(PNSN-tremor-driven envelope peaks)\ncolor ~ detection density",
        fontsize=10,
    )
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
