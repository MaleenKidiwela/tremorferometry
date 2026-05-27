"""End-to-end smoke test on synthetic data with a known dv/v pattern.

Generates per-family HDF5 stacks with an imposed Cascadia-ETS-like dv/v signal,
runs the measure step, and produces a figure overlaying imposed (truth) and
recovered dv/v.

Run after `pip install -e .` in the project env:

    python scripts/00_smoke_synthetic.py --out figures/smoke_synthetic.png
"""

from __future__ import annotations

import argparse
from datetime import datetime, timedelta
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from tremorferometry.measure import measure_many
from tremorferometry.synthetic import (
    dvv_ets_pattern,
    make_bin_edges,
    master_template,
    write_synthetic_family,
)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=Path("figures/smoke_synthetic.png"))
    ap.add_argument("--stacks-dir", type=Path, default=Path("data/synthetic_stacks"))
    ap.add_argument("--dvv-min", type=float, default=-0.005)
    ap.add_argument("--families", type=int, default=3)
    ap.add_argument("--stations", type=int, default=5)
    ap.add_argument("--bin-days", type=int, default=2)
    ap.add_argument("--noise", type=float, default=0.02)
    args = ap.parse_args()

    fs = 100.0
    ets_start = datetime(2026, 5, 27)
    ets_end = ets_start + timedelta(days=15)
    t_start = ets_start - timedelta(days=90)
    t_end = ets_end + timedelta(days=45)
    ref_window = (-90, -30)  # days relative to ets_start

    edges = make_bin_edges(t_start, t_end, args.bin_days)
    template = master_template(fs=fs)
    stations = [f"S{ii:02d}" for ii in range(args.stations)]

    args.stacks_dir.mkdir(parents=True, exist_ok=True)
    family_ids = [f"FAM{ii:03d}" for ii in range(args.families)]
    for k, fid in enumerate(family_ids):
        write_synthetic_family(
            out_path=args.stacks_dir / f"{fid}.h5",
            family_id=fid,
            stations=stations,
            bin_edges=edges,
            template=template,
            fs=fs,
            ets_start=ets_start,
            ets_end=ets_end,
            dvv_min=args.dvv_min,
            noise_level=args.noise,
            n_det_per_bin=50,
            seed=k,
        )

    ref_start = ets_start + timedelta(days=ref_window[0])
    ref_end = ets_start + timedelta(days=ref_window[1])
    h5_paths = sorted(args.stacks_dir.glob("*.h5"))
    df = measure_many(
        h5_paths=h5_paths,
        coda_window=(5.0, 25.0),
        fs=fs,
        eps_max=0.02,
        n_eps=401,
        ref_start=ref_start,
        ref_end=ref_end,
        min_cc=0.5,
        n_workers=min(8, max(1, len(h5_paths))),
    )

    truth_t = pd.date_range(t_start, t_end, freq="6h").to_pydatetime()
    truth_y = np.array([dvv_ets_pattern(t, ets_start, ets_end, dvv_min=args.dvv_min) for t in truth_t])

    fig, ax = plt.subplots(figsize=(11, 4.5))
    ax.plot(truth_t, truth_y * 100.0, "k-", lw=1.5, label="imposed (truth)", zorder=5)
    palette = plt.cm.viridis(np.linspace(0.1, 0.9, len(family_ids)))
    for c, fid in zip(palette, family_ids):
        sub = df[df["family_id"] == fid]
        for sta, ssub in sub.groupby("station"):
            ssub = ssub.sort_values("t_center")
            ax.plot(ssub["t_center"], ssub["dvv"] * 100.0, "-", color=c, lw=0.5, alpha=0.4)
            ax.scatter(ssub["t_center"], ssub["dvv"] * 100.0, s=8, color=c, alpha=0.5)
        ax.scatter([], [], color=c, label=fid)
    ax.axvspan(ets_start, ets_end, color="C3", alpha=0.12, label="ETS")
    ax.axhline(0.0, color="0.6", lw=0.5)
    ax.set_ylabel("dv/v (%)")
    ax.set_xlabel("time")
    ax.set_title(
        f"synthetic smoke: imposed dv/v_min = {args.dvv_min*100:.2f}%, "
        f"{args.families} families x {args.stations} stations, noise={args.noise}"
    )
    ax.legend(loc="lower right", ncol=2, fontsize=8)
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=140)
    print(f"wrote {args.out}")

    # quick recovery stat
    if not df.empty:
        df["t_center"] = pd.to_datetime(df["t_center"])
        df["dvv_true"] = df["t_center"].apply(
            lambda t: dvv_ets_pattern(t.to_pydatetime(), ets_start, ets_end, dvv_min=args.dvv_min)
        )
        residual = df["dvv"] - df["dvv_true"]
        print(
            f"recovery rms = {residual.std():.5f}, "
            f"bias = {residual.mean():+.5f}, "
            f"n = {len(df)}"
        )


if __name__ == "__main__":
    main()
