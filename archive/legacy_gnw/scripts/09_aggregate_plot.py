"""Aggregate dv/v across families + stations, Gaussian-smooth, and produce the
production-quality ETS dv/v figure with LFE rate panel and QC outputs.

This is what Phase C did inline; promoted here so the analysis is reproducible.

Usage:
    python scripts/09_aggregate_plot.py \\
        --config configs/ets_2010_vi.yaml \\
        --dvv data/dvv/phaseB_6sta.parquet \\
        --detections data/detections_lin_ets_2010_vi.parquet \\
        --out figures/ets_2010_vi_aggregate.png
"""

from __future__ import annotations

import argparse
from datetime import timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.dates import DateFormatter
from scipy.ndimage import gaussian_filter1d

from tremorferometry.config import load_config
from tremorferometry.qc import detection_count_independence, spatial_coherence


def _attach_n_det(df: pd.DataFrame, detections: pd.DataFrame, cfg, bin_days: int) -> pd.DataFrame:
    det = detections.copy()
    det["time"] = pd.to_datetime(det["time"])
    keep = set(df["family_id"].astype(str))
    det = det[det["family_id"].astype(str).isin(keep)]
    t0 = cfg.episode.t_start + timedelta(days=-90)
    t1 = cfg.episode.t_end + timedelta(days=30)
    edges = pd.date_range(t0, t1, freq=f"{bin_days}D")
    centers = [edges[i] + (edges[i + 1] - edges[i]) / 2 for i in range(len(edges) - 1)]
    det["bin_idx"] = pd.cut(det["time"], bins=edges, right=False, labels=False)
    n_per_fam_bin = det.groupby(["family_id", "bin_idx"]).size().reset_index(name="n_det")
    n_per_fam_bin["t_center"] = (
        n_per_fam_bin["bin_idx"]
        .astype("Int64")
        .map(lambda i: centers[i] if pd.notna(i) and 0 <= i < len(centers) else pd.NaT)
    )
    n_per_fam_bin = n_per_fam_bin.dropna(subset=["t_center"])
    n_per_fam_bin["t_center"] = pd.to_datetime(n_per_fam_bin["t_center"])
    df = df.drop(columns=["n_det"], errors="ignore").merge(
        n_per_fam_bin[["family_id", "t_center", "n_det"]],
        on=["family_id", "t_center"], how="left",
    )
    return df, n_per_fam_bin


def _wmean(g: pd.DataFrame) -> pd.Series:
    w = g["weight"].values
    if w.sum() == 0 or len(g) < 2:
        return pd.Series({"mean": np.nan, "se": np.nan, "n": len(g)})
    m = float(np.average(g["dvv"], weights=w))
    var = np.average((g["dvv"] - m) ** 2, weights=w)
    se = float(np.sqrt(var) / np.sqrt(len(g)))
    return pd.Series({"mean": m, "se": se, "n": len(g)})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--dvv", type=Path, required=True)
    ap.add_argument("--detections", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--smooth-sigma-days", type=float, default=3.0,
                    help="Gaussian sigma for time-axis smoothing of the aggregate")
    ap.add_argument("--y-limit", type=float, default=1.5, help="dv/v y-axis +/- limit (%)")
    args = ap.parse_args()

    cfg = load_config(args.config)
    df = pd.read_parquet(args.dvv)
    df["t_center"] = pd.to_datetime(df["t_center"])
    detections = pd.read_parquet(args.detections)

    df, n_per_fam_bin = _attach_n_det(df, detections, cfg, cfg.stack.bin_days)
    df["weight"] = df["cc_max"].clip(0, 1) * np.sqrt(df["n_det"].fillna(1))

    agg = df.groupby("t_center").apply(_wmean).reset_index().sort_values("t_center")
    agg = agg.dropna(subset=["mean"])

    sigma_bins = max(0.5, args.smooth_sigma_days / cfg.stack.bin_days)
    ts = agg.set_index("t_center").resample(f"{cfg.stack.bin_days}D").first()
    m = ts["mean"].values
    nan_mask = np.isnan(m)
    m_smooth = gaussian_filter1d(np.where(nan_mask, np.nanmedian(m), m), sigma=sigma_bins)
    ts["mean_smooth"] = m_smooth
    ts.loc[nan_mask, "mean_smooth"] = np.nan
    n_total = n_per_fam_bin.groupby("t_center")["n_det"].sum()

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(13, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    station_colors = dict(zip(
        sorted(df["station"].unique()), plt.cm.tab10.colors
    ))
    for sta, ssub in df.groupby("station"):
        ssub = ssub.sort_values("t_center")
        ax.scatter(ssub["t_center"], ssub["dvv"] * 100, s=5, alpha=0.18,
                   color=station_colors[sta])
    ax.errorbar(agg["t_center"], agg["mean"] * 100, yerr=agg["se"] * 100,
                fmt="o", ms=3, lw=1, color="black", alpha=0.6,
                label="per-bin weighted mean")
    ax.plot(ts.index, ts["mean_smooth"] * 100, lw=2.5, color="C3",
            label=f"Gaussian smoothed (sigma={args.smooth_sigma_days:.0f} d)", zorder=10)
    ax.axvspan(cfg.episode.t_start, cfg.episode.t_end, color="C3", alpha=0.10,
               label=f"{cfg.event_id}")
    ax.axhline(0, color="0.5", lw=0.5)
    for sta, col in station_colors.items():
        ax.scatter([], [], s=18, color=col, alpha=0.6, label=sta)
    ax.set_ylabel("dv/v (%)")
    ax.set_ylim(-args.y_limit, args.y_limit)
    ax.legend(loc="lower right", ncol=2, fontsize=8)
    ax.grid(alpha=0.3)
    ax.set_title(
        f"Cascadia dv/v aggregate: {df['family_id'].nunique()} LFE families x "
        f"{df['station'].nunique()} stations\n"
        f"coda {cfg.dvv.coda_window[0]:.0f}-{cfg.dvv.coda_window[1]:.0f} s, "
        f"{cfg.dvv.freq_band[0]:.0f}-{cfg.dvv.freq_band[1]:.0f} Hz, "
        f"{cfg.stack.bin_days}-day bins"
    )

    ax2.bar(n_total.index, n_total.values, width=cfg.stack.bin_days * 0.9,
            color="gray", alpha=0.7)
    ax2.axvspan(cfg.episode.t_start, cfg.episode.t_end, color="C3", alpha=0.10)
    ax2.set_ylabel(f"LFEs / {cfg.stack.bin_days}-d bin\n(all families)")
    ax2.set_xlabel("time")
    ax2.grid(alpha=0.3)
    for x in (ax, ax2):
        x.xaxis.set_major_formatter(DateFormatter("%b %Y"))

    plt.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)
    print(f"wrote {args.out}")

    sc_ok, sc_info = spatial_coherence(df)
    det_ok, det_info = detection_count_independence(df)
    print(f"QC: spatial_coherence pass={sc_ok}  {sc_info}")
    print(f"QC: n_det independence pass={det_ok}  {det_info}")
    # ETS-window summary
    inside = df[(df["t_center"] >= pd.Timestamp(cfg.episode.t_start))
              & (df["t_center"] <= pd.Timestamp(cfg.episode.t_end))]
    if not inside.empty:
        print(
            f"ETS dv/v summary: n={len(inside)}  "
            f"mean={inside['dvv'].mean() * 100:+.3f}%  "
            f"median={inside['dvv'].median() * 100:+.3f}%  "
            f"std={inside['dvv'].std() * 100:.3f}%"
        )


if __name__ == "__main__":
    main()
