"""Parallel per-patch coda dv/v — identical math to dvv_coda_51.py, but the per-patch
loop (the 382k serial stretch_dvv calls that made the single-threaded version take 40+ min)
runs across a ProcessPool. Each patch is independent: build its count-weighted all-time
reference, stretch every daily stack vs that ref on the coda window, keep cc>=cc_min.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

sys.path.insert(0, "src")
from tremorferometry.dvv import stretch_dvv  # noqa: E402

_C = {}


def _init(fs, t_min, t_max, eps_max, n_eps, cc_min):
    _C.update(fs=fs, t_min=t_min, t_max=t_max, eps_max=eps_max, n_eps=n_eps, cc_min=cc_min)


def _patch_worker(arg):
    """arg = (patch_name, s_patch float32 (k,nsamp), n_patch float64 (k,), date_strs list)."""
    patch, s, n, date_strs = arg
    weighted = (s.T * n).sum(axis=1) / n.sum()
    nrm = float(np.linalg.norm(weighted))
    if nrm <= 0:
        return []
    ref = (weighted / nrm).astype(np.float64)
    rows = []
    for i in range(len(date_strs)):
        cur = s[i].astype(np.float64)
        try:
            r = stretch_dvv(ref, cur, fs=_C["fs"], t_min=_C["t_min"], t_max=_C["t_max"],
                            eps_max=_C["eps_max"], n_eps=_C["n_eps"])
        except Exception:
            continue
        if r.cc_max < _C["cc_min"]:
            continue
        rows.append((patch, date_strs[i], float(r.dvv), float(r.cc_max),
                     float(r.dvv_err), int(n[i])))
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--npz", default="data/long_window_daily_GNW.npz")
    p.add_argument("--window", nargs=2, type=float, default=[1.0, 3.0])
    p.add_argument("--cc-min", type=float, default=0.80)
    p.add_argument("--eps-max", type=float, default=0.02)
    p.add_argument("--n-eps", type=int, default=401)
    p.add_argument("--station", default="GNW")
    p.add_argument("--workers", type=int, default=28)
    p.add_argument("--out-csv", default="data/daily_dvv_GNW_coda_1to3.csv")
    p.add_argument("--out-fig", default="figures/smoke_dvv_GNW_coda_1to3.png")
    args = p.parse_args()

    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed

    print(f"[1/3] Loading {args.npz}...", flush=True)
    d = np.load(args.npz, allow_pickle=True)
    stacks = d["stacks"]; patches = d["patches"]
    dates = pd.to_datetime(d["dates"]); n_det = d["n_det"].astype(np.float64)
    t = d["t"]; fs = float(d["fs"])
    w_lo, w_hi = args.window
    t_min = w_lo - t[0]; t_max = w_hi - t[0]
    uniq = np.unique(patches)
    print(f"  {len(stacks)} rows, {len(uniq)} patches; window {w_lo}-{w_hi}s -> "
          f"array {t_min:.2f}-{t_max:.2f}s; workers={args.workers}", flush=True)

    date_str_all = pd.Series(dates).dt.strftime("%Y-%m-%d").values
    tasks = []
    for tmpl in uniq:
        mask = patches == tmpl
        tasks.append((tmpl, stacks[mask], n_det[mask], list(date_str_all[mask])))

    print(f"[2/3] Per-patch dv/v across {len(tasks)} patches...", flush=True)
    ctx = mp.get_context("spawn")
    out_rows = []
    done = 0
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx, initializer=_init,
                             initargs=(fs, t_min, t_max, args.eps_max, args.n_eps, args.cc_min)) as ex:
        futs = [ex.submit(_patch_worker, tk) for tk in tasks]
        for f in as_completed(futs):
            out_rows.extend(f.result())
            done += 1
            if done % 10 == 0:
                print(f"  {done}/{len(tasks)} patches, {len(out_rows):,} measurements", flush=True)

    df = pd.DataFrame(out_rows, columns=["patch", "date", "dvv", "cc_max", "dvv_err", "n_det"])
    print(f"  {len(df):,} measurements across {df['patch'].nunique()} patches", flush=True)
    print(f"  mean cc {df['cc_max'].mean():.4f}, median dv/v {df['dvv'].median()*100:.4f}%", flush=True)
    df.to_csv(args.out_csv, index=False)
    print(f"  saved {args.out_csv}", flush=True)

    print("[3/3] Plotting...", flush=True)
    df["date"] = pd.to_datetime(df["date"]); df["dvv_pct"] = df["dvv"] * 100
    fig, ax = plt.subplots(1, 1, figsize=(13, 5.5))
    n_patches = 0
    for tmpl in sorted(df["patch"].unique()):
        sub = df[df["patch"] == tmpl].sort_values("date")
        if len(sub) < 30:
            continue
        roll = sub.set_index("date")["dvv_pct"].rolling("60D", min_periods=10).median()
        ax.plot(roll.index, roll.values, alpha=0.22, lw=0.6)
        n_patches += 1
    m = df.set_index("date")["dvv_pct"].sort_index().rolling("60D", min_periods=10).median()
    ax.plot(m.index, m.values, color="k", lw=1.8, label="cross-patch 60-d median")
    ax.axhline(0, color="r", lw=0.6, alpha=0.7)
    # mark the 2001 Nisqually M6.8
    nis = pd.Timestamp("2001-02-28")
    if df["date"].min() <= nis <= df["date"].max():
        ax.axvline(nis, color="b", lw=1.2, ls="--", alpha=0.8, label="2001 Nisqually M6.8")
    ax.set_ylabel("dv/v (%)"); ax.set_ylim(-0.1, 0.1); ax.grid(True, alpha=0.3)
    ax.set_xlabel("date")
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.set_title(f"{args.station} dv/v coda {w_lo}-{w_hi}s, all-time-mean ref  --  "
                 f"{n_patches} patches, {len(df):,} measurements, mean cc {df['cc_max'].mean():.3f}")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout(); fig.savefig(args.out_fig, dpi=150)
    print(f"  wrote {args.out_fig}", flush=True)


if __name__ == "__main__":
    main()
