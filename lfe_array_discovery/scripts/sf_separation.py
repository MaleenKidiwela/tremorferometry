#!/usr/bin/env python
"""REAL vs NULL as a function of every criterion -- the evidence, in place of a candidate count.

For each column of the per-template table (`sf_table.py`), this sweeps thresholds and reports how
many real and how many null templates sit beyond each one.  Nothing is cut: the sweep IS the
result.  Two summary statistics per column:

  * AUC -- the probability that a randomly drawn REAL template scores better than a randomly drawn
    NULL one (Mann-Whitney U / n_real n_null).  0.5 = the two populations are indistinguishable on
    that column.  Reported with the normal-approximation standard error, so "0.52 +- 0.02" reads
    as what it is.
  * LIFT -- at the strictest threshold that still leaves a usable null count, the ratio of the
    fraction of real templates beyond it to the fraction of null templates beyond it.  Lift 1.0 =
    no enrichment; lift 3.0 = the real chain puts three times the share of its population there.

Columns are then RANKED by |AUC - 0.5|, which is what "the strongest discriminators" means when
no threshold has been chosen in advance.

  python scripts/sf_separation.py --real table_X_real.csv --null table_X_null.csv --key X \
      --out results/sf_final/sep_X
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
sys.path.insert(0, HERE)
import sf_common as SF                                                      # noqa: E402

# column -> ("hi" if a LARGER value is the more family-like, "lo" otherwise)
COLS = [
    ("n_det", "hi"), ("n_days", "hi"), ("med_cc", "hi"), ("max_nmad", "hi"),
    ("rate_tremor", "hi"), ("fa_quiet", "lo"), ("tremor_quiet_ratio", "hi"),
    ("co_timing", "hi"), ("co_excess", "hi"), ("co_p", "lo"),
    ("n_sta", "hi"), ("rms", "lo"), ("ctrl_margin", "hi"), ("abs_dz", "lo"),
    ("dedupe_n_partners", "lo"), ("near_chestler_km", "lo"),
    ("restack_move_km", "lo"),
    ("sealed_n_det", "hi"), ("sealed_n_days", "hi"), ("sealed_rate", "hi"),
    ("sealed_fa", "lo"), ("sealed_ratio", "hi"), ("sealed_move_km", "lo"),
]
BOOL_COLS = ["dedupe_kept", "located", "ctrl_pass", "edge", "dep_rail",
             "depth_window_15_58", "tierA", "tierB", "tierC"]
MIN_NULL = 5


def auc(x_real, x_null, direction):
    """P(real beats null), ties counted as half.  NaNs are dropped and reported separately."""
    r = np.asarray(x_real, float)
    z = np.asarray(x_null, float)
    r = r[np.isfinite(r)]
    z = z[np.isfinite(z)]
    if len(r) < 3 or len(z) < 3:
        return np.nan, np.nan, len(r), len(z)
    if direction == "lo":
        r, z = -r, -z
    allv = np.concatenate([r, z])
    ranks = pd.Series(allv).rank(method="average").to_numpy()   # ties get the average rank
    U = ranks[:len(r)].sum() - len(r) * (len(r) + 1) / 2
    A = U / (len(r) * len(z))
    se = np.sqrt(A * (1 - A) / min(len(r), len(z)))
    return float(A), float(se), len(r), len(z)


def sweep(R, Z, col, direction, n=12):
    r = R[col].to_numpy(float)
    z = Z[col].to_numpy(float)
    pool = np.concatenate([r[np.isfinite(r)], z[np.isfinite(z)]])
    if len(pool) < 10:
        return pd.DataFrame()
    qs = np.unique(np.quantile(pool, np.linspace(0.0, 0.99, n)))
    rows = []
    nr, nz = int(np.isfinite(r).sum()), int(np.isfinite(z).sum())
    for t in qs:
        if direction == "hi":
            a, b = int(np.nansum(r >= t)), int(np.nansum(z >= t))
        else:
            a, b = int(np.nansum(r <= t)), int(np.nansum(z <= t))
        fr, fz = a / max(nr, 1), b / max(nz, 1)
        rows.append(dict(col=col, direction=direction, thr=float(t), n_real=a, n_null=b,
                         frac_real=fr, frac_null=fz,
                         lift=(fr / fz) if fz > 0 else np.inf,
                         lift_smoothed=((a + 0.5) / (nr + 1)) / ((b + 0.5) / (nz + 1))))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", required=True)
    ap.add_argument("--null", required=True)
    ap.add_argument("--key", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    R = pd.read_csv(a.real)
    Z = pd.read_csv(a.null)
    SF.log(f"[{a.key}] real {len(R)} rows, null {len(Z)} rows")
    sweeps, summ = [], []
    for col, d in COLS:
        if col not in R.columns or col not in Z.columns:
            continue
        S = sweep(R, Z, col, d)
        if len(S):
            sweeps.append(S)
        A, se, nr, nz = auc(R[col], Z[col], d)
        best = None
        if len(S):
            ok = S[S.n_null >= MIN_NULL]
            src = ok if len(ok) else S
            src = src.replace([np.inf], np.nan).dropna(subset=["lift"])
            if len(src):
                best = src.loc[src.lift.idxmax()]
        summ.append(dict(col=col, direction=d, auc=A, auc_se=se,
                         sep=(abs(A - 0.5) if A == A else np.nan),
                         n_real_finite=nr, n_null_finite=nz,
                         best_lift=(float(best.lift) if best is not None else np.nan),
                         best_thr=(float(best.thr) if best is not None else np.nan),
                         best_n_real=(int(best.n_real) if best is not None else 0),
                         best_n_null=(int(best.n_null) if best is not None else 0)))
    for col in BOOL_COLS:
        if col not in R.columns or col not in Z.columns:
            continue
        fr = float(R[col].astype(bool).mean())
        fz = float(Z[col].astype(bool).mean())
        summ.append(dict(col=col, direction="flag", auc=np.nan, auc_se=np.nan,
                         sep=abs(fr - fz), n_real_finite=len(R), n_null_finite=len(Z),
                         best_lift=(fr / fz) if fz > 0 else np.inf,
                         best_thr=np.nan,
                         best_n_real=int(R[col].astype(bool).sum()),
                         best_n_null=int(Z[col].astype(bool).sum())))
    SW = pd.concat(sweeps, ignore_index=True) if sweeps else pd.DataFrame()
    SU = pd.DataFrame(summ)
    SU["favours"] = ["--" if not (a == a) else ("real" if a > 0.5 else "null")
                     for a in SU.auc]
    SU.loc[SU.direction == "flag", "favours"] = [
        ("real" if r >= nn else "null")
        for r, nn in zip(SU.loc[SU.direction == "flag", "best_n_real"] / max(len(R), 1),
                         SU.loc[SU.direction == "flag", "best_n_null"] / max(len(Z), 1))]
    SU = SU.sort_values("sep", ascending=False)
    SW.to_csv(f"{a.out}_sweep.csv", index=False)
    SU.to_csv(f"{a.out}_summary.csv", index=False)

    # the ranked table: order the REAL rows by the strongest numeric discriminator
    num = SU[(SU.direction != "flag") & SU.auc.notna()]
    fav = num[num.auc > 0.5].sort_values("auc", ascending=False)   # rank only on a column where
    pick = fav.iloc[0] if len(fav) else None                       # the REAL chain scores better
    rank_col = str(pick.col) if pick is not None else "n_det"
    rank_dir = str(pick.direction) if pick is not None else "hi"
    RK = R.sort_values(rank_col, ascending=(rank_dir == "lo")).copy()
    zz = Z[rank_col].to_numpy(float)
    zz = zz[np.isfinite(zz)]

    def pct(v):
        if not np.isfinite(v) or not len(zz):
            return np.nan
        return float((zz <= v).mean() if rank_dir == "hi" else (zz >= v).mean())

    RK["null_percentile_on_rank_col"] = [pct(v) for v in RK[rank_col].to_numpy(float)]
    RK["rank_col"] = rank_col
    RK.to_csv(f"{a.out}_ranked.csv", index=False)
    out = dict(key=a.key, n_real=int(len(R)), n_null=int(len(Z)),
               rank_col=rank_col, rank_direction=rank_dir,
               top=[{k: (None if (isinstance(v, float) and not np.isfinite(v)) else v)
                     for k, v in r.items()}
                    for r in SU.head(8).to_dict("records")],
               max_sep=float(SU.sep.max()) if len(SU) else np.nan,
               best_favouring_real=(dict(col=str(pick.col), auc=float(pick.auc),
                                         auc_se=float(pick.auc_se), sep=float(pick.sep),
                                         best_lift=float(pick.best_lift),
                                         best_thr=float(pick.best_thr),
                                         n_real=int(pick.best_n_real),
                                         n_null=int(pick.best_n_null))
                                    if pick is not None else None),
               n_reverse_separations=int(((num.auc < 0.5) & (num.sep > 0.10)).sum()),
               any_separation=bool(pick is not None and float(pick.sep) > 0.10))
    json.dump(out, open(f"{a.out}.json", "w"), indent=1, default=float)
    print(json.dumps(out, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
