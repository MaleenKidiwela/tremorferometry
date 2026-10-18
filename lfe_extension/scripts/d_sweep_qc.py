#!/usr/bin/env python
"""Track D -- a read-only scan of the Track C sweep for the anomalies the brief asks about:
missing months, months with zero observable family-days or zero quiet family-days, sites flagged by
Track C's own QC (quiet-day rate above 2 lambda), HF-ratio shift alarms, station-day errors and any
month whose assembled catalog count disagrees with `months_done.csv`."""
from __future__ import annotations
import argparse, glob, os, subprocess, sys
import numpy as np, pandas as pd

EXT = "/home/jovyan/tremorferometry/lfe_extension"
TC = f"{EXT}/results/trackC"
TD = f"{EXT}/results/trackD"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check-ready", type=int, default=0,
                    help="run c_month_ready.py (read-only) on up to N missing months to see which are "
                         "permanently blocked and by which station-day")
    a = ap.parse_args()
    MD = pd.read_csv(f"{TC}/months_done.csv").drop_duplicates("ym", keep="last").sort_values("ym")
    want = pd.period_range("2014-01", "2026-09", freq="M").astype(str)
    have = set(MD.ym)
    missing = [m for m in want if m not in have]
    print(f"[months] {len(have)} finished of {len(want)} expected; missing {len(missing)}")
    if missing:
        print("  " + ", ".join(missing[:24]) + (" ..." if len(missing) > 24 else ""))
    if a.check_ready and missing:
        print(f"[readiness] probing {min(a.check_ready, len(missing))} missing month(s) with "
              f"c_month_ready.py (read-only)")
        env = dict(os.environ, PYTHONPATH="/opt/conda/lib/python3.13/site-packages")
        blocked = []
        for m in missing[:a.check_ready]:
            r = subprocess.run([sys.executable, f"{EXT}/scripts/c_month_ready.py", m],
                               capture_output=True, text=True, env=env, cwd=EXT)
            out = r.stdout
            if "runnable=False" in out:
                unf = [l.strip() for l in out.splitlines() if "UNFINISHED" in l]
                n = out.split("runnable=False (")[1].split(" ")[0] if "runnable=False (" in out else "?"
                blocked.append(dict(ym=m, n_unfinished=n,
                                    first_unfinished=unf[0] if unf else "",
                                    kind="EPOCH-BOUNDARY off-by-one" if n == "1" else "download incomplete"))
        B = pd.DataFrame(blocked)
        if len(B):
            print(B.to_string(index=False))
            B.to_csv(f"{TD}/sweep_blocked_months.csv", index=False)
            k = B[B.kind.str.startswith("EPOCH")]
            if len(k):
                print(f"[BLOCKER] {len(k)} month(s) wait on ONE station-day that the download queue marked "
                      f"`no_epoch` and never fetched -- they can never become runnable without "
                      f"`c_run_all.sh --allow-partial` or a one-line fix in c_month_ready.py")
    FD = pd.concat([pd.read_parquet(p, columns=["fam", "date", "observable", "quiet_flag"])
                    for p in sorted(glob.glob(f"{TC}/family_days_*.parquet"))], ignore_index=True)
    FD["ym"] = FD.date.astype(str).str[:7]
    g = FD.groupby("ym").agg(n=("fam", "size"), n_obs=("observable", "sum"),
                             n_quiet=("quiet_flag", "sum")).reset_index()
    g["frac_unobservable"] = 1 - g.n_obs / g.n
    print(f"[family-days] unobservable fraction: median {g.frac_unobservable.median():.4f}, "
          f"max {g.frac_unobservable.max():.4f} ({g.loc[g.frac_unobservable.idxmax(), 'ym']})")
    z = g[g.n_obs == 0]
    print(f"[ZERO observable family-days] {len(z)} month(s)" + ("" if not len(z) else ": " + ", ".join(z.ym)))
    zq = g[g.n_quiet == 0]
    print(f"[zero QUIET family-days] {len(zq)} month(s)"
          + ("" if not len(zq) else ": " + ", ".join(zq.ym) + "  (these months contribute no floor; "
             "lambda falls back through the family-epoch / epoch-season hierarchy)"))
    Q = pd.read_csv(f"{TC}/qc_summary.csv").drop_duplicates(["ym", "kind", "site"], keep="last")
    fs = Q[(Q.kind == "site") & (Q.flag_2lambda.astype(str) == "True")]
    print(f"[QC flag_2lambda] {len(fs)} site-month(s) above 2 lambda")
    if len(fs):
        print(fs.groupby("site").size().sort_values(ascending=False).to_string())
    fh = Q[Q.flag_hfr_shift.astype(str) == "True"]
    print(f"[QC flag_hfr_shift] {len(fh)} month(s)" + ("" if not len(fh) else ": " + ", ".join(sorted(set(fh.ym)))))
    ne = Q[(Q.kind == "network")] if "network" in set(Q.kind) else Q
    err = pd.to_numeric(ne.get("n_station_days_error", pd.Series(dtype=float)), errors="coerce").fillna(0)
    print(f"[station-days with a caching error] total {int(err.sum())}")
    x = f"{TD}/catalog_month_crosscheck.csv"
    if os.path.exists(x):
        X = pd.read_csv(x)
        bad = X[(X.d_ext.abs() > 0) | (X.d_core.abs() > 0)]
        print(f"[catalog vs months_done] {len(bad)} month(s) disagree")
        if len(bad):
            print(bad.to_string(index=False))
    g.to_csv(f"{TD}/sweep_qc_by_month.csv", index=False)
    print(f"[write] {TD}/sweep_qc_by_month.csv")


if __name__ == "__main__":
    main()
