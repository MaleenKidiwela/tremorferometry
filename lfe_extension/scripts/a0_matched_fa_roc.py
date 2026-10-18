#!/usr/bin/env python
"""Track A step 0(c) -- the MATCHED-FALSE-ALARM ROC table (Merlin plan-review F3: E6 compared variants at fixed MAD, which is
a re-parameterisation, not a comparison; the only honest comparison is own recall at the SAME quiet-day false-alarm rate).

For each E6 pass and each target quiet-day FA in {1.0, 0.5, 0.1} detections per template per day:
  T*(FA)  = the threshold at which the pass reaches that FA, log-interpolated in FA between the tabulated 0.5-MAD grid points
            (the grid is {7, 7.5, 8, 8.5, 9, 10} from --budget-thr; FA falls near-exponentially in T, so log is the right axis);
  recall  = own-family recall (Bostock 2012 TEST half, +-2 s) at T*, reported BOTH ways:
            `recall_interp` linearly interpolated in T between the two bracketing grid recalls (= Merlin's recipe, the number
            the accept bar is checked against) and `recall_exact` recomputed at T* from the detections parquet (no interpolation);
  prec    = own precision on ETS days (>= 100 Bostock events) at T*, linearly interpolated in T.
Quiet days = 2012 doy 228-242 (no Bostock event, no PNSN tremor in the family box).  FA, matches and precision come from the
budget CSV (family-day-threshold counts); recall_exact from the declustered detections.  Matching: --one2one (greedy by time,
Merlin F11) or nearest-neighbour; for the OWN-family recall the two agree exactly (within-family declustering already makes the
detection side one-to-one), so the table is matcher-independent.
Inputs: lfe_catalog_repair/results (READ-ONLY).  Output: lfe_extension/results/e6_matched_fa_roc.csv (+ a per-T grid dump)."""
from __future__ import annotations
import argparse
import numpy as np, pandas as pd

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
OUT = "/home/jovyan/tremorferometry/lfe_extension/results"
QUIET = set(range(228, 243))
POST = set(range(270, 280))
TAGS = ["e6_sites", "e6_mergeA", "e6_mergeB", "e6_mergeB_cap25", "e6_mergeB_coinc", "e6_lev_a", "e6_lev_b", "e6_lev_c"]
TARGETS = [1.0, 0.5, 0.1]


def min_absdiff(a, b_sorted):
    a = np.asarray(a, float)
    if len(a) == 0:
        return np.zeros(0)
    if len(b_sorted) == 0:
        return np.full(len(a), np.inf)
    j = np.searchsorted(b_sorted, a)
    lo = b_sorted[np.clip(j - 1, 0, len(b_sorted) - 1)]; hi = b_sorted[np.clip(j, 0, len(b_sorted) - 1)]
    return np.minimum(np.abs(a - lo), np.abs(a - hi))


def one2one_absdiff(a, b_sorted, tol=2.0):
    a = np.asarray(a, float); out = np.full(len(a), np.inf)
    if len(a) == 0 or len(b_sorted) == 0:
        return out
    j = 0; nb = len(b_sorted)
    for i in np.argsort(a, kind="stable"):
        ti = a[i]
        while j < nb and b_sorted[j] < ti - tol:
            j += 1
        if j < nb and abs(b_sorted[j] - ti) <= tol:
            out[i] = abs(b_sorted[j] - ti); j += 1
    return out


def own_recall_at(det, test, T, match):
    hit = np.zeros(len(test), bool)
    dT = det[det.net_mad >= T]
    for fk, g in test.groupby("bkey"):
        td = np.sort(dT[dT.fam == float(fk)].t.to_numpy())
        hit[test.index.get_indexer(g.index)] = match(g.t_origin.to_numpy(), td, 2.0) <= 2.0
    return hit


def loginterp_T(fa_by_T, target):
    """smallest bracketing pair on the T grid with FA crossing `target`; log-linear in FA -> T*.  None if never reached."""
    Ts = sorted(fa_by_T); 
    for a, b in zip(Ts[:-1], Ts[1:]):
        fa, fb = fa_by_T[a], fa_by_T[b]
        if fa >= target >= fb and fa > 0 and fb > 0 and fa != fb:
            return a + (b - a) * (np.log(fa) - np.log(target)) / (np.log(fa) - np.log(fb)), a, b
    if fa_by_T[Ts[0]] <= target:
        return Ts[0], Ts[0], Ts[0]
    return None, None, None


def lin(x, x0, x1, y0, y1):
    return y0 if x1 == x0 else y0 + (x - x0) * (y1 - y0) / (x1 - x0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tags", default=",".join(TAGS))
    ap.add_argument("--L", type=float, default=8.0)
    ap.add_argument("--one2one", action="store_true", default=True)
    ap.add_argument("--nn", dest="one2one", action="store_false")
    ap.add_argument("--in-dir", default=IN)
    ap.add_argument("--gt-dir", default=IN)
    ap.add_argument("--out", default=f"{OUT}/e6_matched_fa_roc.csv")
    ap.add_argument("--grid-out", default=f"{OUT}/e6_fa_recall_grid.csv")
    a = ap.parse_args()
    match = one2one_absdiff if a.one2one else (lambda q, b, tol=2.0: min_absdiff(q, b))
    gt = pd.read_parquet(f"{a.gt_dir}/e5r_gt_bostock2012.parquet").reset_index(drop=True)
    gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    per_day = gt.groupby("doy").size(); ETS = set(per_day[per_day >= 100].index)
    test = gt[(gt.half == "test") & (~gt.guard)].copy().reset_index(drop=True)
    grid, roc = [], []
    for tag in a.tags.split(","):
        bud = pd.read_csv(f"{a.in_dir}/e5_budget_{tag}.csv"); bud = bud[bud.L == a.L]
        det = pd.read_parquet(f"{a.in_dir}/e5r_detections_{tag}.parquet"); det = det[det.L == a.L]
        fa_by_T, rec_by_T, prec_by_T, post_by_T = {}, {}, {}, {}
        for T in sorted(bud["T"].unique()):
            b = bud[bud["T"] == T]
            q = b[b.doy.isin(QUIET)]; e = b[b.doy.isin(ETS)]; p = b[b.doy.isin(POST)]
            fa_by_T[T] = float(q.n_det.mean()); post_by_T[T] = float(p.n_det.mean())
            prec_by_T[T] = float(e.n_match_own.sum() / max(e.n_det.sum(), 1))
            rec_by_T[T] = float(own_recall_at(det, test, T, match).mean())
            grid.append(dict(tag=tag, L=a.L, T=T, quiet_fa=fa_by_T[T], post_fa=post_by_T[T],
                             ets_precision_own=prec_by_T[T], recall_own=rec_by_T[T],
                             quiet_roll_site=float(q.n_det_roll_site.mean()) if "n_det_roll_site" in q else np.nan,
                             quiet_site_over_data=(float(q.n_det_roll_site.mean()) / max(fa_by_T[T], 1e-9)) if "n_det_roll_site" in q else np.nan))
        row = dict(tag=tag, L=a.L)
        for tgt in TARGETS:
            Tst, T0, T1 = loginterp_T(fa_by_T, tgt)
            if Tst is None:
                row[f"T_at_FA{tgt}"] = np.nan; row[f"recall_FA{tgt}"] = np.nan
                row[f"recall_exact_FA{tgt}"] = np.nan; row[f"precision_FA{tgt}"] = np.nan
                continue
            row[f"T_at_FA{tgt}"] = Tst
            row[f"recall_FA{tgt}"] = lin(Tst, T0, T1, rec_by_T[T0], rec_by_T[T1])
            row[f"recall_exact_FA{tgt}"] = float(own_recall_at(det, test, Tst, match).mean())
            row[f"precision_FA{tgt}"] = lin(Tst, T0, T1, prec_by_T[T0], prec_by_T[T1])
        roc.append(row)
    G = pd.DataFrame(grid); G.to_csv(a.grid_out, index=False)
    R = pd.DataFrame(roc); R.to_csv(a.out, index=False)
    pd.set_option("display.width", 250)
    print("[grid] quiet FA / own recall / ETS own precision on the budget threshold grid")
    print(G.round(4).to_string(index=False))
    print(f"\n[ROC] own recall of Bostock 2012 test-half events at matched quiet-day FA "
          f"(matching = {'one-to-one greedy' if a.one2one else 'nearest neighbour'})")
    print(R.round(4).to_string(index=False))
    print(f"[write] {a.out}\n[write] {a.grid_out}")


if __name__ == "__main__":
    main()
