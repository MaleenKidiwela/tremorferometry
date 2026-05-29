"""Coda dv/v with TWO references: all-time-mean AND pre-Nisqually-EQ baseline.

For each family in a long-window-stacks npz, produce:
  - dv/v vs all-time-mean LFE stack (canonical, post-event-symmetric)
  - dv/v vs PRE-EQ-only stack (any pre-EQ trend not flattened by post-EQ data)

Saves two CSVs and a side-by-side comparison figure with the EQ date marked.

Designed for the Nisqually-2001 precursor question at GNW, but generic.
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
    p.add_argument("--npz", required=True,
                   help="long_window_daily_*.npz with stacks/patches/dates/n_det/t/fs")
    p.add_argument("--window", nargs=2, type=float, default=[1.0, 3.0],
                   help="coda window in seconds, LFE-aligned")
    p.add_argument("--cc-min", type=float, default=0.80)
    p.add_argument("--eps-max", type=float, default=0.02)
    p.add_argument("--n-eps", type=int, default=401)
    p.add_argument("--eq-date", required=True,
                   help="Earthquake date YYYY-MM-DD (e.g. 2001-02-28 for Nisqually)")
    p.add_argument("--pre-ref-start", default=None,
                   help="Start of pre-EQ reference window; default = earliest available")
    p.add_argument("--out-prefix", required=True,
                   help="Output prefix; will write {prefix}_alltime.csv, "
                        "{prefix}_pre.csv and {prefix}_dual.png")
    return p.parse_args()


def compute_dvv(stacks, dates, n_det, t_axis, fs, ref, t_min, t_max,
                cc_min, eps_max, n_eps):
    rows = []
    for i, day in enumerate(dates):
        cur = stacks[i].astype(np.float64)
        try:
            r = stretch_dvv(ref, cur, fs=fs, t_min=t_min, t_max=t_max,
                            eps_max=eps_max, n_eps=n_eps)
        except Exception:
            continue
        if r.cc_max < cc_min:
            continue
        rows.append(dict(date=day, dvv=float(r.dvv), cc_max=float(r.cc_max),
                         dvv_err=float(r.dvv_err), n_det=int(n_det[i])))
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    eq_date = pd.to_datetime(args.eq_date)
    pre_start = pd.to_datetime(args.pre_ref_start) if args.pre_ref_start else None

    print(f"[1/3] Loading {args.npz}...")
    d = np.load(args.npz, allow_pickle=True)
    stacks = d["stacks"]; patches = d["patches"]
    dates = pd.to_datetime(d["dates"]); n_det = d["n_det"]
    t_axis = d["t"]; fs = float(d["fs"])
    print(f"  rows={len(stacks)}, patches={len(np.unique(patches))}, "
          f"t-axis {t_axis[0]:.2f} to {t_axis[-1]:.2f} s, fs={fs}")

    w_lo, w_hi = args.window
    t_min = w_lo - t_axis[0]; t_max = w_hi - t_axis[0]
    print(f"  coda window {w_lo}-{w_hi} s LFE-aligned -> "
          f"{t_min:.2f}-{t_max:.2f} s array-time")

    print(f"[2/3] Per-family dv/v with TWO references...")
    rows_at: list[dict] = []
    rows_pre: list[dict] = []
    skipped_pre = 0
    for fam in sorted(np.unique(patches)):
        m = patches == fam
        s = stacks[m]; nd = n_det[m]; dd = dates[m]

        # All-time reference
        weighted = (s.T * nd.astype(np.float64)).sum(axis=1) / nd.sum()
        nrm = float(np.linalg.norm(weighted))
        if nrm <= 0: continue
        ref_at = (weighted / nrm).astype(np.float64)

        # Pre-EQ reference
        pre_mask = dd < eq_date
        if pre_start is not None:
            pre_mask &= dd >= pre_start
        if pre_mask.sum() < 10:
            ref_pre = None
            skipped_pre += 1
        else:
            s_pre = s[pre_mask]; nd_pre = nd[pre_mask]
            w_pre = (s_pre.T * nd_pre.astype(np.float64)).sum(axis=1) / nd_pre.sum()
            nrm_pre = float(np.linalg.norm(w_pre))
            if nrm_pre <= 0:
                ref_pre = None
                skipped_pre += 1
            else:
                ref_pre = (w_pre / nrm_pre).astype(np.float64)

        df_at = compute_dvv(s, dd, nd, t_axis, fs, ref_at,
                            t_min, t_max, args.cc_min, args.eps_max, args.n_eps)
        df_at["patch"] = fam
        rows_at.append(df_at)

        if ref_pre is not None:
            df_pre = compute_dvv(s, dd, nd, t_axis, fs, ref_pre,
                                 t_min, t_max, args.cc_min, args.eps_max, args.n_eps)
            df_pre["patch"] = fam
            rows_pre.append(df_pre)

    df_at = pd.concat(rows_at, ignore_index=True) if rows_at else pd.DataFrame()
    df_pre = pd.concat(rows_pre, ignore_index=True) if rows_pre else pd.DataFrame()
    print(f"  all-time:     {len(df_at):,} measurements across {df_at.patch.nunique()} patches")
    print(f"  pre-EQ ref:   {len(df_pre):,} measurements across {df_pre.patch.nunique()} patches "
          f"({skipped_pre} families had insufficient pre-EQ data)")

    out_at = f"{args.out_prefix}_alltime.csv"
    out_pre = f"{args.out_prefix}_pre.csv"
    df_at.to_csv(out_at, index=False)
    df_pre.to_csv(out_pre, index=False)
    print(f"  saved {out_at}, {out_pre}")

    print(f"[3/3] Plotting...")
    for df, label in [(df_at, "all-time"), (df_pre, "pre-EQ")]:
        df["date"] = pd.to_datetime(df["date"])
        df["dvv_pct"] = df["dvv"] * 100

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 8), sharex=True)

    for ax, df, label in [(ax1, df_at, "all-time-mean reference"),
                          (ax2, df_pre, f"pre-EQ reference (start: "
                                        f"{args.pre_ref_start or 'earliest available'})")]:
        if len(df) == 0:
            ax.set_title(f"{label} -- no measurements")
            continue
        # per-patch rolling medians thin
        for fam in sorted(df.patch.unique()):
            sub = df[df.patch == fam].sort_values("date")
            if len(sub) < 30:
                continue
            roll = sub.set_index("date")["dvv_pct"].rolling("60D", min_periods=10).median()
            ax.plot(roll.index, roll.values, alpha=0.20, lw=0.5)
        # cross-patch median
        s = df.set_index("date")["dvv_pct"].sort_index()
        m = s.rolling("60D", min_periods=10).median()
        ax.plot(m.index, m.values, color="k", lw=1.8, label="cross-patch 60-d median")
        ax.axhline(0, color="r", lw=0.6, alpha=0.6)
        ax.axvline(eq_date, color="purple", lw=1.5, alpha=0.8,
                   label=f"EQ {args.eq_date}")
        ax.set_ylabel("dv/v (%)")
        ax.set_ylim(-0.3, 0.3)
        ax.grid(True, alpha=0.3)
        ax.set_title(label)
        ax.legend(loc="upper right", fontsize=9)

    ax2.set_xlabel("date")
    ax2.xaxis.set_major_locator(mdates.YearLocator(2))
    ax2.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    fig.suptitle(f"GNW coda 1-3 s dv/v -- dual reference, EQ marker {args.eq_date}",
                 fontsize=11, y=0.995)
    fig.tight_layout()
    out_fig = f"{args.out_prefix}_dual.png"
    fig.savefig(out_fig, dpi=150)
    print(f"  wrote {out_fig}")


if __name__ == "__main__":
    main()
