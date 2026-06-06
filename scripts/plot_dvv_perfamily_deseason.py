"""Per-family DESEASONALIZED dv/v overlay (each family's own monthly climatology
removed, then 60-d median). Reproduces smoke_dvv_GNW_perfamily_deseason.png.

Deseasonalize = subtract each family's <clim-start>..<clim-end> monthly climatology
(the well-sampled modern era), so only non-seasonal change remains.
"""
from __future__ import annotations

import argparse

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--dvv-csv", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--clim-start", default="2011-01-01")
    p.add_argument("--clim-end", default="2026-12-31")
    p.add_argument("--smooth-days", type=int, default=60)
    p.add_argument("--min-fam-days", type=int, default=200,
                   help="skip a family with fewer than this many measurement days")
    p.add_argument("--eq-date", default="2001-02-28")
    p.add_argument("--label", default="response-removed, per-era, 1-4s")
    p.add_argument("--no-deseason", action="store_true",
                   help="do NOT remove each family's seasonal cycle; plot raw per-family dv/v")
    p.add_argument("--cc-min", type=float, default=0.95, help="drop measurements below this cc_max")
    p.add_argument("--err-max", type=float, default=0.30, help="drop measurements with dvv_err(%%) above this (0=off)")
    p.add_argument("--start-year", type=int, default=2000, help="x-axis / data start year")
    p.add_argument("--gap-days", type=int, default=45,
                   help="break a family's line where consecutive measurements are more than this apart")
    p.add_argument("--zoom-end-year", type=int, default=0,
                   help="if set, limit the x-axis end to this year (data/climatology unaffected)")
    return p.parse_args()


def main():
    args = parse_args()
    d = pd.read_csv(args.dvv_csv)
    d["date"] = pd.to_datetime(d["date"])
    d["dvv"] = d["dvv"] * 100.0

    fig, ax = plt.subplots(figsize=(14, 6))
    fams = sorted(d["patch"].unique())
    n_used = 0
    t0 = pd.Timestamp(f"{args.start_year}-01-01")
    for fam in fams:
        g = d[d["patch"] == fam].sort_values("date")
        # QC + start-year filter
        g = g[g["cc_max"] >= args.cc_min]
        if args.err_max > 0:
            g = g[g["dvv_err"] * 100.0 <= args.err_max]
        g = g[g["date"] >= t0]
        if g["date"].nunique() < args.min_fam_days:
            continue
        if args.no_deseason:
            vals = g["dvv"].values                         # raw per-family dv/v
        else:
            clim_win = g[(g["date"] >= args.clim_start) & (g["date"] <= args.clim_end)]
            if len(clim_win) < 50:
                continue
            clim = clim_win.groupby(clim_win["date"].dt.month)["dvv"].mean()
            vals = g["dvv"].values - g["date"].dt.month.map(clim).values
        s = pd.Series(vals, index=pd.DatetimeIndex(g["date"].values))
        s = s.groupby(level=0).median().sort_index()       # collapse same-day
        sm = s.rolling(f"{args.smooth_days}D").median().dropna()
        if len(sm) < 2:
            continue
        # break the line across real data gaps (don't plot through them)
        x = sm.index.to_numpy()
        y = sm.to_numpy().astype(float)
        gap = np.diff(x).astype("timedelta64[D]").astype(int) > args.gap_days
        y[1:][gap] = np.nan                                # NaN at first point after each gap
        ax.plot(x, y, lw=0.8, alpha=0.7)
        n_used += 1

    ax.axhline(0, color="k", lw=0.6)
    if args.eq_date:
        ax.axvline(pd.Timestamp(args.eq_date), color="blue", lw=1.2, ls="--",
                   zorder=6, label="2001 Nisqually")
    ax.set_ylabel("dv/v (%)" if args.no_deseason else "deseasonalized dv/v (%)")
    ax.set_xlabel("date")
    ax.set_ylim(-0.09, 0.09)
    x_end = pd.Timestamp(f"{args.zoom_end_year}-12-31") if args.zoom_end_year else d["date"].max()
    ax.set_xlim(t0, x_end)
    ax.xaxis.set_major_locator(mdates.YearLocator(2) if not args.zoom_end_year else mdates.MonthLocator(interval=3))
    if args.zoom_end_year:
        ax.xaxis.set_minor_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m" if args.zoom_end_year else "%Y"))
    if args.no_deseason:
        ax.set_title(f"GNW per-family dv/v ({args.label}, {args.smooth_days}-d median) "
                     f"- {n_used} families", fontsize=10)
    else:
        ax.set_title(f"GNW per-family DESEASONALIZED dv/v ({args.label}, each family's own "
                     f"seasonal cycle removed, {args.smooth_days}-d median) - {n_used} families",
                     fontsize=10)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out} ({n_used} families)")


if __name__ == "__main__":
    main()
