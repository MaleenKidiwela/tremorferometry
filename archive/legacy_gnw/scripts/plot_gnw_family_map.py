"""Map the high-SNR GNW densification seed families, styled like
figures/smoke_pgc_nllb_family_map.png (Salish Sea cartopy basemap with family
locations + station markers).

Input : data/gnw_pnsn_families.summary.csv (lat, lon, snr, n_members, n_years)
Output: figures/smoke_gnw_family_map.png
"""

from __future__ import annotations

import argparse

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import pandas as pd

# Station coords (FDSN, level=station)
STATIONS = {
    "GNW": (47.5641, -122.8250, "UW.GNW (densified)"),
    "PGC": (48.6498, -123.4521, "CN.PGC"),
    "NLLB": (49.2271, -123.9882, "CN.NLLB"),
}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", default="data/gnw_pnsn_families.summary.csv")
    p.add_argument("--min-snr", type=float, default=15.0)
    p.add_argument("--out", default="figures/smoke_gnw_family_map.png")
    return p.parse_args()


def main():
    args = parse_args()
    s = pd.read_csv(args.summary)
    hi = s[s["snr"] >= args.min_snr].copy()
    print(f"{len(hi)} seed families with SNR >= {args.min_snr} "
          f"(of {len(s)} total)")

    proj = ccrs.PlateCarree()
    fig = plt.figure(figsize=(10, 9))
    ax = plt.axes(projection=proj)
    ax.set_extent([-124.3, -122.3, 46.7, 49.4], crs=proj)

    # Basemap features (match reference style: tan land, light-blue ocean)
    ax.add_feature(cfeature.OCEAN.with_scale("10m"), facecolor="#cfe2f3")
    ax.add_feature(cfeature.LAND.with_scale("10m"), facecolor="#f5ead6")
    ax.add_feature(cfeature.LAKES.with_scale("10m"), facecolor="#cfe2f3")
    ax.add_feature(cfeature.COASTLINE.with_scale("10m"), lw=0.6, edgecolor="#555")
    ax.add_feature(cfeature.BORDERS.with_scale("10m"), lw=0.5,
                   edgecolor="#888", linestyle=":")

    # Families: diamonds colored by SNR, sized by member count
    sc = ax.scatter(
        hi["lon"], hi["lat"], c=hi["snr"], cmap="viridis",
        s=25 + 6 * hi["n_members"].clip(upper=20), marker="D",
        edgecolor="k", linewidth=0.4, zorder=5,
        transform=proj, label=f"GNW seed family (n={len(hi)})",
    )
    cb = fig.colorbar(sc, ax=ax, shrink=0.6, pad=0.02)
    cb.set_label("template SNR")

    # Stations
    for name, (lat, lon, lab) in STATIONS.items():
        big = name == "GNW"
        ax.plot(lon, lat, marker="^", color="red" if big else "#666",
                markersize=15 if big else 10, markeredgecolor="k",
                linestyle="none", zorder=6, transform=proj,
                label=lab)
        ax.text(lon + 0.04, lat + 0.02, name, fontsize=9 if big else 8,
                fontweight="bold" if big else "normal", zorder=7,
                transform=proj)

    gl = ax.gridlines(draw_labels=True, lw=0.4, color="gray",
                      alpha=0.4, linestyle="--")
    gl.top_labels = False
    gl.right_labels = False

    ax.set_title(
        f"GNW densification seeds: {len(hi)} high-SNR LFE families "
        f"(SNR >= {args.min_snr:.0f})\n"
        "marker size ~ member count; color ~ template SNR",
        fontsize=10,
    )
    ax.legend(loc="upper left", fontsize=8, framealpha=0.9)

    fig.savefig(args.out, dpi=150, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
