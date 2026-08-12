"""Map discovered LFE families above an SNR floor around a target station.

Reusable across stations (not GNW-hardcoded): pass --station / --station-lat /
--station-lon / --summary.  Styled like the GNW/PGC family maps.

Example:
  python scripts/plot_family_map.py --station HDW \
      --station-lat 47.649 --station-lon -123.053 \
      --summary data/hdw_pnsn_families.summary.csv --min-snr 10 \
      --out figures/smoke_hdw_family_map.png
"""
from __future__ import annotations

import argparse

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import pandas as pd

# Nearby reference stations for context
REF_STATIONS = {
    "GNW": (47.5641, -122.8250),
    "PGC": (48.6498, -123.4521),
    "NLLB": (49.2271, -123.9882),
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", required=True)
    p.add_argument("--station", required=True)
    p.add_argument("--station-lat", type=float, required=True)
    p.add_argument("--station-lon", type=float, required=True)
    p.add_argument("--min-snr", type=float, default=10.0)
    p.add_argument("--out", required=True)
    p.add_argument("--pad", type=float, default=0.35,
                   help="degrees of padding around the family/station extent")
    return p.parse_args()


def main():
    args = parse_args()
    s = pd.read_csv(args.summary)
    hi = s[s["snr"] >= args.min_snr].copy()
    print(f"{len(hi)} families with SNR >= {args.min_snr} (of {len(s)} total)")

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

    sc = ax.scatter(
        hi["lon"], hi["lat"], c=hi["snr"], cmap="viridis",
        s=25 + 6 * hi["n_members"].clip(upper=20), marker="D",
        edgecolor="k", linewidth=0.4, zorder=5,
        transform=proj, label=f"{args.station} family (n={len(hi)})",
    )
    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("template SNR")

    # target station (big red), reference stations (grey)
    ax.plot(args.station_lon, args.station_lat, marker="^", color="red",
            markersize=15, markeredgecolor="k", linestyle="none", zorder=6,
            transform=proj, label=f"{args.station} (target)")
    ax.text(args.station_lon + 0.03, args.station_lat + 0.02, args.station,
            fontsize=10, fontweight="bold", zorder=7, transform=proj)
    for name, (lat, lon) in REF_STATIONS.items():
        if name == args.station:
            continue
        # only show reference stations that fall within THIS station's map extent,
        # else bbox_inches='tight' stretches the figure to reach far-away labels.
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
        f"{args.station}: {len(hi)} LFE families with SNR >= {args.min_snr:.0f}\n"
        "marker size ~ member count; color ~ template SNR",
        fontsize=10,
    )
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
