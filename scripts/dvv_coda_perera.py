"""Per-ERA-reference coda dv/v — removes the instrument/sample-rate steps (GNW: 50Hz 1995-2010,
40Hz 2011-2018, 100Hz 2019-2026). Same stretch_dvv math as dvv_coda_parallel.py, but each daily
stack is referenced to a count-weighted mean built from ONLY its own sample-rate era, so the
resampling/response offset between eras cancels. Parallel across patches.

Tradeoff: each era is centered on its own reference (no absolute cross-era level comparison), but
WITHIN-era variation (seasonal, ETS, the 2000-2001 pre-Nisqually signal) is preserved artifact-free.
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.dates as mdates, matplotlib.pyplot as plt
sys.path.insert(0, "src")
from tremorferometry.dvv import stretch_dvv  # noqa: E402

_C = {}
def _init(fs, t_min, t_max, eps_max, n_eps, cc_min, bounds):
    _C.update(fs=fs, t_min=t_min, t_max=t_max, eps_max=eps_max, n_eps=n_eps, cc_min=cc_min, bounds=bounds)

def _era_of(date_strs, bounds):
    # date-based epochs (YYYY-MM-DD compares lexically); bounds are exact transition dates
    arr = np.asarray(date_strs)
    e = np.zeros(len(arr), int)
    for b in bounds:
        e += (arr >= b).astype(int)
    return e

def _patch_worker(arg):
    patch, s, n, date_strs = arg
    era = _era_of(date_strs, _C["bounds"])
    rows = []
    for e in np.unique(era):
        m = era == e
        se, ne = s[m], n[m]
        idx = np.where(m)[0]
        weighted = (se.T * ne).sum(axis=1) / ne.sum()
        nrm = float(np.linalg.norm(weighted))
        if nrm <= 0:
            continue
        ref = (weighted / nrm).astype(np.float64)
        for j, i in enumerate(idx):
            try:
                r = stretch_dvv(ref, se[j].astype(np.float64), fs=_C["fs"],
                                t_min=_C["t_min"], t_max=_C["t_max"],
                                eps_max=_C["eps_max"], n_eps=_C["n_eps"])
            except Exception:
                continue
            if r.cc_max < _C["cc_min"]:
                continue
            rows.append((patch, date_strs[i], float(r.dvv), float(r.cc_max), float(r.dvv_err), int(ne[j])))
    return rows

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--npz", default="data/long_window_daily_GNW.npz")
    p.add_argument("--window", nargs=2, type=float, default=[1.0, 3.0])
    p.add_argument("--cc-min", type=float, default=0.80)
    p.add_argument("--eps-max", type=float, default=0.02)
    p.add_argument("--n-eps", type=int, default=401)
    p.add_argument("--station", default="GNW")
    p.add_argument("--era-bounds", default="2010-09-10,2019-05-07",
                   help="exact instrument-change dates; GNW 50Hz<2010-09-10, 40Hz, 100Hz>=2019-05-07")
    p.add_argument("--workers", type=int, default=28)
    p.add_argument("--out-csv", default="data/daily_dvv_GNW_coda_1to3_perera.csv")
    p.add_argument("--out-fig", default="figures/smoke_dvv_GNW_coda_1to3_perera.png")
    args = p.parse_args()
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor, as_completed

    bounds = [x.strip() for x in args.era_bounds.split(",") if x.strip()]
    d = np.load(args.npz, allow_pickle=True)
    stacks, patches = d["stacks"], d["patches"]
    dates = pd.to_datetime(d["dates"]); n_det = d["n_det"].astype(np.float64)
    t = d["t"]; fs = float(d["fs"]); w_lo, w_hi = args.window
    t_min, t_max = w_lo - t[0], w_hi - t[0]
    date_str_all = pd.Series(dates).dt.strftime("%Y-%m-%d").values
    uniq = np.unique(patches)
    print(f"[1/3] {len(stacks)} stacks, {len(uniq)} patches; era bounds {bounds}", flush=True)
    tasks = [(tm, stacks[patches == tm], n_det[patches == tm], list(date_str_all[patches == tm])) for tm in uniq]

    print(f"[2/3] per-era dv/v across {len(tasks)} patches...", flush=True)
    ctx = mp.get_context("spawn")
    rows = []
    with ProcessPoolExecutor(max_workers=args.workers, mp_context=ctx, initializer=_init,
                             initargs=(fs, t_min, t_max, args.eps_max, args.n_eps, args.cc_min, bounds)) as ex:
        for f in as_completed([ex.submit(_patch_worker, tk) for tk in tasks]):
            rows.extend(f.result())
    df = pd.DataFrame(rows, columns=["patch", "date", "dvv", "cc_max", "dvv_err", "n_det"])
    df.to_csv(args.out_csv, index=False)
    print(f"  {len(df):,} measurements, mean cc {df['cc_max'].mean():.3f} -> {args.out_csv}", flush=True)

    print("[3/3] plot...", flush=True)
    df["date"] = pd.to_datetime(df["date"]); df["dvv_pct"] = df["dvv"] * 100
    fig, ax = plt.subplots(figsize=(13, 5.5))
    np_ = 0
    for tm in sorted(df["patch"].unique()):
        sub = df[df["patch"] == tm].sort_values("date")
        if len(sub) < 30:
            continue
        roll = sub.set_index("date")["dvv_pct"].rolling("60D", min_periods=10).median()
        ax.plot(roll.index, roll.values, alpha=0.18, lw=0.6); np_ += 1
    m = df.set_index("date")["dvv_pct"].sort_index().rolling("60D", min_periods=10).median()
    ax.plot(m.index, m.values, color="k", lw=1.8, label="cross-patch 60-d median")
    ax.axhline(0, color="r", lw=0.6, alpha=0.7)
    for b in bounds:
        ax.axvline(pd.Timestamp(b), color="green", ls=":", lw=1.0, alpha=0.6)
    ax.axvline(pd.Timestamp("2001-02-28"), color="b", ls="--", lw=1.2, alpha=0.8, label="2001 Nisqually M6.8")
    ax.set_ylim(-0.1, 0.1); ax.grid(True, alpha=0.3)
    ax.set_ylabel("dv/v (%)"); ax.set_xlabel("date")
    ax.xaxis.set_major_locator(mdates.YearLocator(2)); ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.set_title(f"{args.station} dv/v coda {w_lo}-{w_hi}s, PER-ERA ref (50/40/100 Hz) — "
                 f"{np_} patches, {len(df):,} meas, mean cc {df['cc_max'].mean():.3f} "
                 f"(green dotted = sample-rate change)")
    ax.legend(loc="upper right", fontsize=9)
    fig.tight_layout(); fig.savefig(args.out_fig, dpi=150)
    print(f"  wrote {args.out_fig}", flush=True)

if __name__ == "__main__":
    main()
