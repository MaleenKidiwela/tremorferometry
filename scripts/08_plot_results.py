"""Plot dv/v(t) per family/station with the ETS interval shaded.

Usage:
    python scripts/08_plot_results.py \\
        --config configs/ets_2010_vi.yaml \\
        --dvv data/dvv/ets_2010_vi.parquet \\
        --out figures/ets_2010_vi.png
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from tremorferometry.config import load_config
from tremorferometry.qc import quality_filter


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--dvv", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    cfg = load_config(args.config)
    df = pd.read_parquet(args.dvv)
    df = quality_filter(df, min_cc=cfg.dvv.min_cc, min_n_det=cfg.dvv.min_detections)
    if df.empty:
        raise SystemExit("no dv/v rows survive quality filter")
    df["t_center"] = pd.to_datetime(df["t_center"])

    families = sorted(df["family_id"].unique())
    fig, axes = plt.subplots(
        len(families), 1, figsize=(11, 2.2 * len(families)), sharex=True
    )
    if len(families) == 1:
        axes = [axes]
    for ax, fid in zip(axes, families):
        sub = df[df["family_id"] == fid]
        for station, ssub in sub.groupby("station"):
            ssub = ssub.sort_values("t_center")
            ax.errorbar(
                ssub["t_center"], ssub["dvv"] * 100.0,
                yerr=ssub["dvv_err"] * 100.0, fmt="o-", ms=3, lw=0.7,
                label=str(station), alpha=0.8,
            )
        ax.axvspan(cfg.episode.t_start, cfg.episode.t_end, color="C3", alpha=0.15, label="ETS")
        ax.axhline(0.0, color="k", lw=0.5)
        ax.set_ylabel("dv/v (%)")
        ax.set_title(f"family {fid}")
        ax.legend(loc="upper right", fontsize=7, ncol=2)
    axes[-1].set_xlabel("time")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
