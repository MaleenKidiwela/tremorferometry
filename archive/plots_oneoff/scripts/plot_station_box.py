"""Plot a station-centered discovery box over the master PNSN tremor catalog.

Shows the tremor density (full Cascadia catalog), the station, a +/-N/S,+/-E/W km
box centered on the station, and a radius circle -- so the box can be eyeballed
before running discovery. Reusable across stations.
"""
from __future__ import annotations

import argparse

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--catalog", default="catalogs/pnsn_tremor_cascadia_full.csv")
    p.add_argument("--station", required=True)
    p.add_argument("--station-lat", type=float, required=True)
    p.add_argument("--station-lon", type=float, required=True)
    p.add_argument("--ns-km", type=float, default=100.0, help="half-extent N/S (km)")
    p.add_argument("--ew-km", type=float, default=100.0, help="half-extent E/W (km)")
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    d = pd.read_csv(args.catalog, usecols=["lat", "lon"])

    coslat = np.cos(np.radians(args.station_lat))
    dlat = args.ns_km / 111.0
    dlon = args.ew_km / (111.0 * coslat)
    lat0, lat1 = args.station_lat - dlat, args.station_lat + dlat
    lon0, lon1 = args.station_lon - dlon, args.station_lon + dlon

    # tremor inside the box
    inb = d[(d.lat >= lat0) & (d.lat <= lat1) & (d.lon >= lon0) & (d.lon <= lon1)]
    nsouth = int((inb.lat < args.station_lat).sum())
    print(f"box lat {lat0:.3f}..{lat1:.3f}  lon {lon0:.3f}..{lon1:.3f}")
    print(f"tremor in box: {len(inb):,} | south of station: {nsouth:,} "
          f"({100*nsouth/max(len(inb),1):.0f}%)")

    pad = 0.45
    ext = [lon0 - pad, lon1 + pad, lat0 - pad, lat1 + pad]
    dview = d[(d.lat >= ext[2]) & (d.lat <= ext[3])
              & (d.lon >= ext[0]) & (d.lon <= ext[1])]

    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(10, 9))
    ax = plt.axes(projection=proj)
    ax.set_extent(ext, crs=proj)
    ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#cfe2f3")
    ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#f5ead6")
    ax.add_feature(cfeature.LAKES.with_scale("10m"), facecolor="#cfe2f3")
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), lw=0.6, edgecolor="#555")
    ax.add_feature(cfeature.BORDERS.with_scale("10m"), lw=0.5, edgecolor="#888",
                   linestyle=":")

    hb = ax.hexbin(dview.lon, dview.lat, gridsize=110, cmap="magma_r",
                   norm=LogNorm(vmin=1), mincnt=1, zorder=4, transform=proj,
                   linewidths=0.1)
    cb = fig.colorbar(hb, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("PNSN tremor count per cell (log)")

    # station-centered box
    ax.plot([lon0, lon1, lon1, lon0, lon0], [lat0, lat0, lat1, lat1, lat0],
            color="red", lw=2.0, zorder=6, transform=proj,
            label=f"±{args.ns_km:.0f}km N/S, ±{args.ew_km:.0f}km E/W box")
    # radius circle (use the N/S km as the radius reference)
    th = np.linspace(0, 2 * np.pi, 200)
    clat = args.station_lat + (args.ns_km / 111.0) * np.sin(th)
    clon = args.station_lon + (args.ns_km / (111.0 * coslat)) * np.cos(th)
    ax.plot(clon, clat, color="red", lw=1.0, ls="--", zorder=6, transform=proj,
            label=f"{args.ns_km:.0f} km radius")

    ax.plot(args.station_lon, args.station_lat, marker="^", color="lime",
            markersize=16, markeredgecolor="k", linestyle="none", zorder=7,
            transform=proj, label=f"{args.station}")
    ax.text(args.station_lon + 0.04, args.station_lat, args.station, fontsize=11,
            fontweight="bold", zorder=8, transform=proj)

    gl = ax.gridlines(draw_labels=True, lw=0.4, color="gray", alpha=0.4,
                      linestyle="--")
    gl.top_labels = gl.right_labels = False
    ax.set_title(
        f"{args.station}-centered discovery box: {len(inb):,} PNSN tremor "
        f"({100*nsouth/max(len(inb),1):.0f}% south of station)", fontsize=10)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
