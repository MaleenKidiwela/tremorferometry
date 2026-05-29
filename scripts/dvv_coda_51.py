"""Compute coda-window (1-3 s) dv/v for all 51 LFE families at PGC.

Input: `data/long_window_daily_all51.npz` from
`build_long_window_daily_all51.py` (n_rows daily L2-normed stacks, each
on a -3 to +10 s axis, with patch + date labels).

Per family:
  - Build all-time reference = sum_d (daily_sum * n_det) / total_LFE_count,
    L2-normalize.
  - For each day: stretch_dvv vs ref on coda window 1-3 s
    (= 4-6 s in array-time from t[0]=-3 s).
  - Keep cc_max >= 0.8.

Outputs CSV + comparison figure vs the existing 0-2 s product.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, "src")
from tremorferometry.dvv import stretch_dvv  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--npz", default="data/long_window_daily_all51.npz")
    p.add_argument("--window", nargs=2, type=float, default=[1.0, 3.0],
                   help="Stretching window in seconds relative to LFE alignment")
    p.add_argument("--cc-min", type=float, default=0.80)
    p.add_argument("--eps-max", type=float, default=0.02)
    p.add_argument("--n-eps", type=int, default=401)
    p.add_argument("--station", default="PGC")
    p.add_argument("--out-csv", default="data/daily_dvv_51_coda_1to3.csv")
    p.add_argument("--out-fig", default="figures/smoke_dvv_51_coda_vs_direct.png")
    p.add_argument("--compare-csv", default="data/daily_dvv_51_alltime_ref.csv")
    return p.parse_args()


def main():
    args = parse_args()
    print(f"[1/3] Loading {args.npz}...")
    d = np.load(args.npz, allow_pickle=True)
    stacks = d["stacks"]
    patches = d["patches"]
    dates = pd.to_datetime(d["dates"])
    n_det = d["n_det"].astype(np.float64)
    t = d["t"]
    fs = float(d["fs"])
    print(f"  {len(stacks)} rows, {len(np.unique(patches))} unique patches, "
          f"t-axis: {t[0]:.2f} to {t[-1]:.2f} s")

    # Map LFE-aligned window to array-time (sample 0 corresponds to t[0])
    w_lo, w_hi = args.window
    t_min = w_lo - t[0]
    t_max = w_hi - t[0]
    print(f"  stretching window: {w_lo}-{w_hi} s LFE-aligned -> "
          f"{t_min:.2f}-{t_max:.2f} s array-time")

    print("[2/3] Per-family dv/v...")
    out_rows: list[dict] = []
    skipped = 0
    for tmpl in np.unique(patches):
        mask = patches == tmpl
        s = stacks[mask]
        n = n_det[mask]
        dd = dates[mask]
        # all-time reference: count-weighted mean LFE
        weighted = (s.T * n).sum(axis=1) / n.sum()
        nrm = float(np.linalg.norm(weighted))
        if nrm <= 0:
            skipped += 1
            continue
        ref = (weighted / nrm).astype(np.float64)
        for i, day in enumerate(dd):
            cur = s[i].astype(np.float64)
            try:
                r = stretch_dvv(ref, cur, fs=fs,
                                t_min=t_min, t_max=t_max,
                                eps_max=args.eps_max, n_eps=args.n_eps)
            except Exception:
                continue
            if r.cc_max < args.cc_min:
                continue
            out_rows.append(dict(
                patch=tmpl, date=day.strftime("%Y-%m-%d"),
                dvv=float(r.dvv), cc_max=float(r.cc_max),
                dvv_err=float(r.dvv_err), n_det=int(n[i]),
            ))
    df = pd.DataFrame(out_rows)
    print(f"  {len(df)} measurements across {df['patch'].nunique()} patches "
          f"(skipped {skipped} families with bad ref)")
    print(f"  mean cc {df['cc_max'].mean():.4f}, median dv/v {df['dvv'].median()*100:.4f}%, "
          f"5-95%ile {df['dvv'].quantile(0.05)*100:.3f} / "
          f"{df['dvv'].quantile(0.95)*100:.3f}%")
    df.to_csv(args.out_csv, index=False)
    print(f"  saved {args.out_csv}")

    print("[3/3] Plotting (coda only)...")
    df["date"] = pd.to_datetime(df["date"])
    df["dvv_pct"] = df["dvv"] * 100

    def cross_patch_median(d: pd.DataFrame, prefix: str | None) -> pd.Series:
        if prefix is not None:
            d = d[d["patch"].str.startswith(prefix)]
        s = d.set_index("date")["dvv_pct"].sort_index()
        return s.rolling("60D", min_periods=10).median()

    fig, ax = plt.subplots(1, 1, figsize=(13, 5.5))
    n_patches = 0
    for tmpl in sorted(df["patch"].unique()):
        sub = df[df["patch"] == tmpl].sort_values("date")
        if len(sub) < 30:
            continue
        roll = sub.set_index("date")["dvv_pct"].rolling("60D", min_periods=10).median()
        ax.plot(roll.index, roll.values, alpha=0.22, lw=0.6)
        n_patches += 1
    m = cross_patch_median(df, None)
    ax.plot(m.index, m.values, color="k", lw=1.8, label="cross-patch 60-d median")
    ax.axhline(0, color="r", lw=0.6, alpha=0.7)
    ax.set_ylabel("dv/v (%)")
    ax.set_ylim(-0.1, 0.1)
    ax.grid(True, alpha=0.3)
    ax.set_xlabel("date")
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.set_title(
        f"{args.station} dv/v on coda window {w_lo}-{w_hi} s, all-time-mean LFE reference  --  "
        f"{n_patches} patches, {len(df):,} measurements, mean cc {df['cc_max'].mean():.3f}"
    )
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout()
    fig.savefig(args.out_fig, dpi=150)
    print(f"  wrote {args.out_fig}")


if __name__ == "__main__":
    main()
