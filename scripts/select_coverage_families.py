"""Coverage-balanced, SNR-prioritized LFE-family selection for densification.

Instead of taking ALL families above an SNR cut (which over-samples dense regions and leaves
azimuthal/distance gaps), bin families by (azimuth sector x distance ring) from the station and
keep the top-K by SNR per cell, with an SNR floor. This gives good azimuth+distance coverage,
prioritizes the strongest templates, loosens to the floor only where nothing better exists, and
never piles 100s of redundant patches into one cell.

Writes a subset summary CSV (same schema as the input) -> feed to densify via --summary-csv.
"""
from __future__ import annotations
import argparse, re
import numpy as np, pandas as pd

def latlon(fid):
    m = re.match(r"(-?\d+\.\d+)_(-?\d+\.\d+)", str(fid))
    return (float(m.group(1)), float(m.group(2))) if m else (np.nan, np.nan)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", required=True, help="families summary csv (family_id, snr, ...)")
    p.add_argument("--station-lat", type=float, required=True)
    p.add_argument("--station-lon", type=float, required=True)
    p.add_argument("--min-snr", type=float, default=10.0, help="SNR floor (loosening threshold)")
    p.add_argument("--az-sectors", type=int, default=12)
    p.add_argument("--dist-rings", default="0,20,40,65,200", help="comma km ring edges")
    p.add_argument("--k", type=int, default=2, help="top-K families per (az,dist) cell")
    p.add_argument("--out", required=True, help="selected subset summary csv")
    p.add_argument("--out-fig", default=None)
    args = p.parse_args()

    s = pd.read_csv(args.summary)
    s["lat"] = s["family_id"].map(lambda f: latlon(f)[0])
    s["lon"] = s["family_id"].map(lambda f: latlon(f)[1])
    s = s.dropna(subset=["lat", "lon"])
    dlat = (s["lat"] - args.station_lat) * 111.0
    dlon = (s["lon"] - args.station_lon) * 111.0 * np.cos(np.radians(args.station_lat))
    s["dist_km"] = np.sqrt(dlat**2 + dlon**2)
    s["az_deg"] = (np.degrees(np.arctan2(dlon, dlat)) % 360)
    pool = s[s["snr"] >= args.min_snr].copy()
    rings = [float(x) for x in args.dist_rings.split(",")]
    pool["azbin"] = (pool["az_deg"] // (360.0 / args.az_sectors)).astype(int)
    pool["rbin"] = pd.cut(pool["dist_km"], rings, labels=False)
    pool = pool.dropna(subset=["rbin"])
    sel = pool.sort_values("snr", ascending=False).groupby(["azbin", "rbin"], sort=False).head(args.k)
    ncell = args.az_sectors * (len(rings) - 1)
    print(f"[select] pool SNR>={args.min_snr}: {len(pool)}  ->  selected {len(sel)} "
          f"(top-{args.k}/cell, {args.az_sectors}az x {len(rings)-1}dist)")
    print(f"  occupied cells {sel.groupby(['azbin','rbin']).ngroups}/{ncell}; "
          f"azimuth {sel['azbin'].nunique()}/{args.az_sectors} sectors; dist {sel['rbin'].nunique()}/{len(rings)-1} rings")
    print(f"  SNR: median {sel['snr'].median():.1f}, >=15: {(sel['snr']>=15).sum()}, "
          f"floor-{args.min_snr:g}..15 gap-fillers: {((sel['snr']>=args.min_snr)&(sel['snr']<15)).sum()}")
    # write subset summary (original columns only)
    orig_cols = [c for c in pd.read_csv(args.summary, nrows=0).columns]
    sel[orig_cols].to_csv(args.out, index=False)
    print(f"  wrote {args.out} ({len(sel)} families)")
    if args.out_fig:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(9, 8))
        ax.scatter(pool["lon"], pool["lat"], s=10, c="lightgray", label=f"pool ({len(pool)})")
        sc = ax.scatter(sel["lon"], sel["lat"], c=sel["snr"], cmap="viridis", vmin=args.min_snr, vmax=30,
                        s=70, edgecolor="k", lw=0.4, label=f"selected ({len(sel)})")
        ax.scatter([args.station_lon], [args.station_lat], marker="^", s=300, c="red", edgecolor="k", zorder=5, label="station")
        plt.colorbar(sc, label="SNR"); ax.legend(fontsize=8); ax.grid(alpha=0.3)
        ax.set_xlabel("lon"); ax.set_ylabel("lat"); ax.set_title(f"coverage-balanced selection ({len(sel)} families)")
        fig.tight_layout(); fig.savefig(args.out_fig, dpi=140); print(f"  wrote {args.out_fig}")

if __name__ == "__main__":
    main()
