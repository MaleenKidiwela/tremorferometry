"""Overlay the empirical downdip (eastern) edge of the ETS tremor band on a
tremor-density + family map around a station.

The downdip edge at each latitude is taken as the 92nd-percentile longitude of
tremor in that lat band (east = larger lon), i.e. where the density tapers out
to the east -- the surface trace of where the plate interface plunges below the
ETS depth window. Lightly smoothed across latitude.
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
    p.add_argument("--families", required=True)
    p.add_argument("--station", required=True)
    p.add_argument("--station-lat", type=float, required=True)
    p.add_argument("--station-lon", type=float, required=True)
    p.add_argument("--min-snr", type=float, default=9.0)
    p.add_argument("--box", nargs=4, type=float, required=True,
                   metavar=("LATMIN", "LATMAX", "LONMIN", "LONMAX"))
    p.add_argument("--pct", type=float, default=92.0,
                   help="percentile of lon per lat-bin defining the eastern edge")
    p.add_argument("--out", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    latmin, latmax, lonmin, lonmax = args.box
    t = pd.read_csv(args.catalog, usecols=["lat", "lon"])
    t = t[(t.lat >= latmin) & (t.lat <= latmax)
          & (t.lon >= lonmin) & (t.lon <= lonmax)]

    # downdip edge: per 0.1-deg lat bin, the pct-th percentile of lon (eastern envelope)
    edges = np.arange(latmin, latmax + 0.1, 0.1)
    elat, elon = [], []
    for a, b in zip(edges[:-1], edges[1:]):
        sub = t.lon[(t.lat >= a) & (t.lat < b)]
        if len(sub) < 30:
            continue
        elat.append((a + b) / 2)
        elon.append(np.percentile(sub, args.pct))
    elat, elon = np.array(elat), np.array(elon)
    # light smoothing (3-pt moving average)
    if len(elon) >= 3:
        elon = np.convolve(elon, np.ones(3) / 3, mode="same")
        elon[0], elon[-1] = elon[1], elon[-2]

    s = pd.read_csv(args.families)
    hi = s[s["snr"] >= args.min_snr]

    pad = 0.2
    ext = [lonmin - pad, lonmax + pad, latmin - pad, latmax + pad]
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

    hb = ax.hexbin(t.lon, t.lat, gridsize=90, cmap="Greys", norm=LogNorm(vmin=1),
                   mincnt=1, zorder=3, transform=proj, alpha=0.8, linewidths=0.1)
    cb = fig.colorbar(hb, ax=ax, shrink=0.55, pad=0.02)
    cb.set_label("PNSN tremor count (log)")

    # downdip edge
    ax.plot(elon, elat, color="red", lw=2.5, zorder=6, transform=proj,
            label=f"downdip (E) edge of tremor (p{args.pct:.0f})")

    sc = ax.scatter(hi["lon"], hi["lat"], c=hi["snr"], cmap="viridis",
                    s=45, marker="D", edgecolor="k", linewidth=0.4, zorder=5,
                    transform=proj, label=f"family SNR>={args.min_snr:g} (n={len(hi)})")
    fig.colorbar(sc, ax=ax, shrink=0.55, pad=0.08).set_label("template SNR")

    ax.plot(args.station_lon, args.station_lat, marker="^", color="red",
            markersize=15, markeredgecolor="k", linestyle="none", zorder=7,
            transform=proj, label=args.station)
    ax.text(args.station_lon + 0.03, args.station_lat + 0.02, args.station,
            fontsize=11, fontweight="bold", zorder=8, transform=proj)

    gl = ax.gridlines(draw_labels=True, lw=0.4, color="gray", alpha=0.4,
                      linestyle="--")
    gl.top_labels = gl.right_labels = False
    ax.set_title(f"{args.station}: ETS tremor band + downdip (eastern) edge\n"
                 "tremor terminates where the slab interface drops below the ETS depth window",
                 fontsize=10)
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)
    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote {args.out} | edge lon {elon.min():.2f}..{elon.max():.2f}")


if __name__ == "__main__":
    main()
