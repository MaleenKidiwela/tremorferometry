#!/usr/bin/env python
"""Track A step 3, rejection test (b) -- HELD-OUT-YEAR recall at MATCHED false-alarm rate.
Compares two 2010 runs (frozen 2012 templates vs iterated templates, same 9 persistent stations, same MAD grid)
on own-family recall of Bostock's 2010 events at the SAME 2010 quiet-day FA per template-day.  2010 quiet days =
doy 222-226 / 250-256 with zero in-box PNSN tremor and zero Bostock events.  Adopt the iteration only if the
matched-FA recall rises by >= +0.01; reject if it falls by > 0.01 (Merlin step 3, test b).
Optionally restricts the comparison to a subset of families (--fams-csv with a bkey column), which is how the
"accepted families only" variant is scored without another GPU pass: a family's network sum is independent of
every other family's, so mixing per family is exact."""
from __future__ import annotations
import argparse
import numpy as np, pandas as pd

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
OUT = "/home/jovyan/tremorferometry/lfe_extension/results"
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"
TARGETS = [1.0, 0.5, 0.25, 0.1]


def one2one(a, b_sorted, tol=2.0):
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


def recall(det, ev, T):
    D = det[det.net_mad >= T]; h = np.zeros(len(ev), bool)
    for fk, g in ev.groupby("bkey"):
        td = np.sort(D[D.fam == float(fk)].t.to_numpy())
        h[ev.index.get_indexer(g.index)] = one2one(g.t_origin.to_numpy(), td) <= 2.0
    return h


def curve(tag, ev, quiet, Tgrid):
    det = pd.read_parquet(f"{OUT}/e5r_detections_{tag}.parquet"); det = det[det.L == 12.0]
    bud = pd.read_csv(f"{OUT}/e5_budget_{tag}.csv"); bud = bud[bud.L == 12.0]
    fd = bud[bud["T"] == sorted(bud["T"].unique())[0]][["doy", "fam"]].drop_duplicates()
    nq = max(int(fd.doy.isin(quiet).sum()), 1)
    rows = []
    for T in Tgrid:
        rows.append(dict(T=T, quiet_fa=len(det[(det.net_mad >= T) & det.doy.isin(quiet)]) / nq,
                         recall=float(recall(det, ev, T).mean()), n=int((det.net_mad >= T).sum())))
    return pd.DataFrame(rows), det


def at_fa(C, target):
    C = C.sort_values("T")
    for (t0, f0, r0), (t1, f1, r1) in zip(C[["T", "quiet_fa", "recall"]].values[:-1], C[["T", "quiet_fa", "recall"]].values[1:]):
        if f0 >= target >= f1 and f0 > 0 and f1 > 0 and f0 != f1:
            x = (np.log(f0) - np.log(target)) / (np.log(f0) - np.log(f1))
            return t0 + x * (t1 - t0), r0 + x * (r1 - r0)
    return np.nan, np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frozen", default="a2_xfer2010"); ap.add_argument("--iterated", required=True)
    ap.add_argument("--fams-csv", default="")
    ap.add_argument("--out", default=f"{OUT}/a3_heldout_matched_fa.csv")
    a = ap.parse_args()
    g10 = pd.read_parquet(f"{OUT}/a2_gt_bostock2010.parquet")
    off = pd.read_csv(f"{IN}/e3_bostock_family_offsets.csv"); off = off[off.accept.astype(bool)]
    omap = dict(zip(off.bkey.astype(int), off.offset_s.astype(float)))
    g10 = g10[g10.bkey.astype(int).isin(omap)].copy(); g10["time"] = g10.time.astype("datetime64[ns]")
    g10["t_origin"] = g10.time.astype("int64") / 1e9 + g10.bkey.astype(int).map(omap).to_numpy()
    g10["doy"] = pd.DatetimeIndex(g10.time).dayofyear
    if a.fams_csv:
        fams = set(pd.read_csv(a.fams_csv).bkey.astype(int))
        g10 = g10[g10.bkey.astype(int).isin(fams)]
    ev = g10.reset_index(drop=True)
    tr = pd.read_csv(PNSN); tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    F = pd.read_csv(f"{IN}/e5r_families_2012.csv")
    box = tr[(tr.time.dt.year == 2010) & tr.lat.between(F.lat.min() - .2, F.lat.max() + .2) & tr.lon.between(F.lon.min() - .2, F.lon.max() + .2)]
    tc = box.groupby(box.time.dt.dayofyear).size(); bc = g10.groupby("doy").size()
    quiet = [d for d in list(range(222, 227)) + list(range(250, 257)) if int(tc.get(d, 0)) == 0 and int(bc.get(d, 0)) == 0]
    Tgrid = [7.0, 7.5, 8.0, 8.25, 8.5, 8.75, 9.0, 10.0]
    Cf, _ = curve(a.frozen, ev, quiet, Tgrid); Ci, _ = curve(a.iterated, ev, quiet, Tgrid)
    Cf["set"] = "frozen"; Ci["set"] = "iterated"
    G = pd.concat([Cf, Ci]); G.to_csv(a.out.replace(".csv", "_grid.csv"), index=False)
    rows = []
    for t in TARGETS:
        Tf, rf = at_fa(Cf, t); Ti, ri = at_fa(Ci, t)
        rows.append(dict(target_fa=t, T_frozen=Tf, recall_frozen=rf, T_iterated=Ti, recall_iterated=ri,
                         delta=ri - rf, adopt=bool(np.isfinite(ri - rf) and ri - rf >= 0.01),
                         reject=bool(np.isfinite(ri - rf) and ri - rf < -0.01)))
    R = pd.DataFrame(rows); R.to_csv(a.out, index=False)
    pd.set_option("display.width", 250)
    print(f"[2010 quiet days] {quiet}; {len(ev)} Bostock 2010 events in {ev.bkey.nunique()} families"
          + (f" (restricted to {a.fams_csv})" if a.fams_csv else ""))
    print("[grid]"); print(G.round(4).to_string(index=False))
    print("\n[held-out recall at matched 2010 quiet-day FA]"); print(R.round(4).to_string(index=False))
    print(f"[write] {a.out}")


if __name__ == "__main__":
    main()
