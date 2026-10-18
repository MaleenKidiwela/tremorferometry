#!/usr/bin/env python
"""Track A step 2 -- score the FROZEN cross-year transfer.  Compares, at the frozen thresholds in MAD units:
  (a) the 2012 ablation on the persistent sites (the PREDICTION) against
  (b) the same 2012-built templates applied to 2010 (the TEST),
on own-family recall (+-2 s, one-to-one greedy by time) of Bostock's events and on quiet-day false alarms per
template-day.  2010 origin times = Bostock DB time + the ACCEPTED 2010 per-family offset (e3_bostock_family_offsets.csv,
accept == True, 90 families) -- exactly the table e4_template_test.py used to place the run.
Quiet days: 2012 doy 228-242; 2010 the doy 222-226 / 250-256 windows Merlin named, reported BOTH literally and restricted
to the days with ZERO in-box PNSN tremor AND zero Bostock events (2010 doy 251 carries 91 tremor detections and 21 Bostock
events and is not quiet).  ACCEPT: 2010 recall within 0.7-1.3x the 2012 ablation value AND 2010 quiet FA within 0.5-2x."""
from __future__ import annotations
import argparse, json
import numpy as np, pandas as pd

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
OUT = "/home/jovyan/tremorferometry/lfe_extension/results"
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"
MAGBINS = [(-9, 1.6, "M<1.6"), (1.6, 1.8, "M1.6-1.8"), (1.8, 2.0, "M1.8-2.0"), (2.0, 9, "M>=2.0")]


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
    D = det[det.net_mad >= T]
    hit = np.zeros(len(ev), bool)
    for fk, g in ev.groupby("bkey"):
        td = np.sort(D[D.fam == float(fk)].t.to_numpy())
        hit[ev.index.get_indexer(g.index)] = one2one(g.t_origin.to_numpy(), td) <= 2.0
    return hit


def load(tag, dirpath, L, hf_parquet="", hf_col="hfr_PGC", hf_cut=np.nan):
    det = pd.read_parquet(f"{dirpath}/e5r_detections_{tag}.parquet"); det = det[det.L == L].reset_index(drop=True)
    bud = pd.read_csv(f"{dirpath}/e5_budget_{tag}.csv"); bud = bud[bud.L == L]
    n_dropped = 0
    if hf_parquet:
        H = pd.read_parquet(hf_parquet); H = H[H.L == L][["doy", "fam", "t", hf_col]]
        det = det.merge(H, on=["doy", "fam", "t"], how="left")
        v = det[hf_col].to_numpy(); drop = np.isfinite(v) & (v > hf_cut); n_dropped = int(drop.sum())
        det = det[~drop].reset_index(drop=True)
    return det, bud, n_dropped


def famdays(bud):
    T0 = sorted(bud["T"].unique())[0]
    return bud[bud["T"] == T0][["doy", "fam"]].drop_duplicates()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag2012", default="a2_abl2012_7p"); ap.add_argument("--tag2010", default="a2_xfer2010")
    ap.add_argument("--dir", default=OUT); ap.add_argument("--L", type=float, default=12.0)
    ap.add_argument("--t-ext", type=float, default=8.25); ap.add_argument("--t-core", type=float, default=8.75)
    ap.add_argument("--hf2012", default=""); ap.add_argument("--hf2010", default="")
    ap.add_argument("--hf-cut", type=float, default=3.6384)
    ap.add_argument("--out", default=f"{OUT}/a2_transfer_report.json")
    a = ap.parse_args()
    d12, b12, n12d = load(a.tag2012, a.dir, a.L, a.hf2012, "hfr_PGC", a.hf_cut)
    d10, b10, n10d = load(a.tag2010, a.dir, a.L, a.hf2010, "hfr_PGC", a.hf_cut)
    # ---- ground truth
    gt12 = pd.read_parquet(f"{IN}/e5r_gt_bostock2012.parquet").reset_index(drop=True)
    ev12 = gt12[(gt12.half == "test") & (~gt12.guard)].copy().reset_index(drop=True)
    g10 = pd.read_parquet(f"{OUT}/a2_gt_bostock2010.parquet")
    off = pd.read_csv(f"{IN}/e3_bostock_family_offsets.csv"); off = off[off.accept.astype(bool)]
    omap = dict(zip(off.bkey.astype(int), off.offset_s.astype(float)))
    g10 = g10[g10.bkey.astype(int).isin(omap)].copy()
    g10["time"] = g10.time.astype("datetime64[ns]")
    g10["t_origin"] = g10.time.astype("int64") / 1e9 + g10.bkey.astype(int).map(omap).to_numpy()
    g10["doy"] = pd.DatetimeIndex(g10.time).dayofyear
    fams12 = set(d12.fam.unique()); fams10 = set(d10.fam.unique())
    ev10 = g10[g10.bkey.astype(float).isin(fams10)].reset_index(drop=True)
    # ---- quiet windows
    tr = pd.read_csv(PNSN); tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    F = pd.read_csv(f"{IN}/e5r_families_2012.csv")
    box = tr[tr.lat.between(F.lat.min() - .2, F.lat.max() + .2) & tr.lon.between(F.lon.min() - .2, F.lon.max() + .2)]
    trem10 = box[box.time.dt.year == 2010].groupby(box[box.time.dt.year == 2010].time.dt.dayofyear).size()
    bost10 = g10.groupby("doy").size()
    Q12 = list(range(228, 243))
    Q10_lit = [d for d in list(range(222, 227)) + list(range(250, 257))]
    Q10 = [d for d in Q10_lit if int(trem10.get(d, 0)) == 0 and int(bost10.get(d, 0)) == 0]
    fd12, fd10 = famdays(b12), famdays(b10)
    days10_run = sorted(fd10.doy.unique())
    Q10 = [d for d in Q10 if d in days10_run]; Q10_lit = [d for d in Q10_lit if d in days10_run]
    rep = dict(tag2012=a.tag2012, tag2010=a.tag2010, L=a.L, T_ext=a.t_ext, T_core=a.t_core,
               hf_veto_applied=bool(a.hf2012 and a.hf2010), hf_cut=a.hf_cut, n_dropped_2012=n12d, n_dropped_2010=n10d,
               days_2010=days10_run, quiet_2010_strict=Q10, quiet_2010_literal=Q10_lit, quiet_2012=Q12,
               tremor_2010_on_literal_quiet={int(d): int(trem10.get(d, 0)) for d in Q10_lit},
               bostock_2010_on_literal_quiet={int(d): int(bost10.get(d, 0)) for d in Q10_lit},
               n_families_2012=len(fams12), n_families_2010=len(fams10),
               n_events_2012_test=int(len(ev12)), n_events_2010=int(len(ev10)))
    rows = []
    for nm, T in (("T_ext", a.t_ext), ("T_core", a.t_core)):
        h12 = recall(d12, ev12, T); h10 = recall(d10, ev10, T)
        f12 = len(d12[(d12.net_mad >= T) & d12.doy.isin(Q12)]) / max(int(fd12.doy.isin(Q12).sum()), 1)
        f10 = len(d10[(d10.net_mad >= T) & d10.doy.isin(Q10)]) / max(int(fd10.doy.isin(Q10).sum()), 1)
        f10l = len(d10[(d10.net_mad >= T) & d10.doy.isin(Q10_lit)]) / max(int(fd10.doy.isin(Q10_lit).sum()), 1)
        r = dict(threshold=nm, T=T, recall_2012_ablation=float(h12.mean()), recall_2010_transfer=float(h10.mean()),
                 recall_ratio=float(h10.mean() / max(h12.mean(), 1e-9)),
                 quiet_fa_2012=f12, quiet_fa_2010_strict=f10, quiet_fa_2010_literal=f10l,
                 fa_ratio_strict=f10 / max(f12, 1e-9), fa_ratio_literal=f10l / max(f12, 1e-9),
                 n_catalog_2012=int((d12.net_mad >= T).sum()), n_catalog_2010=int((d10.net_mad >= T).sum()),
                 precision_own_2010=float(d10[d10.net_mad >= T].match_own.mean()) if (d10.net_mad >= T).any() else np.nan,
                 precision_own_2012=float(d12[d12.net_mad >= T].match_own.mean()) if (d12.net_mad >= T).any() else np.nan)
        # the 2010 strict-quiet set straddles the SSE (223/224/226 are PRE, 252-256 POST); the 2012 quiet window is pre-SSE
        # only, so the pre-SSE subset is the like-for-like comparison and the post subset is reported next to it.
        qpre = [d for d in Q10 if d < 240]; qpost = [d for d in Q10 if d >= 240]
        for lbl, dd in (("pre", qpre), ("post", qpost)):
            nd = max(int(fd10.doy.isin(dd).sum()), 1)
            r[f"quiet_fa_2010_{lbl}SSE"] = len(d10[(d10.net_mad >= T) & d10.doy.isin(dd)]) / nd
            r[f"fa_ratio_{lbl}SSE"] = r[f"quiet_fa_2010_{lbl}SSE"] / max(f12, 1e-9)
            r[f"quiet_days_2010_{lbl}SSE"] = dd
        r["accept_recall_0.7_1.3"] = bool(0.7 <= r["recall_ratio"] <= 1.3)
        r["accept_fa_0.5_2.0"] = bool(0.5 <= r["fa_ratio_strict"] <= 2.0)
        r["FAILURE_recall_below_half"] = bool(r["recall_ratio"] < 0.5)
        for lo, hi, lbl in MAGBINS:
            r[f"rec2012_{lbl}"] = float(h12[((ev12.mag >= lo) & (ev12.mag < hi)).to_numpy()].mean())
            r[f"rec2010_{lbl}"] = float(h10[((ev10.mag >= lo) & (ev10.mag < hi)).to_numpy()].mean())
        rows.append(r)
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/a2_transfer_table.csv", index=False)
    rep["rows"] = rows
    json.dump(rep, open(a.out, "w"), indent=1, default=float)
    pd.set_option("display.width", 260)
    print(f"[2010] days {days10_run[0]}-{days10_run[-1]}; strict-quiet days {Q10} (literal {Q10_lit}); "
          f"{len(ev10)} Bostock 2010 events in {ev10.bkey.nunique()} detected families")
    print(R[["threshold", "T", "recall_2012_ablation", "recall_2010_transfer", "recall_ratio", "accept_recall_0.7_1.3",
             "quiet_fa_2012", "quiet_fa_2010_strict", "fa_ratio_strict", "accept_fa_0.5_2.0",
             "quiet_fa_2010_preSSE", "fa_ratio_preSSE", "quiet_fa_2010_postSSE", "fa_ratio_postSSE",
             "precision_own_2012", "precision_own_2010"]].round(4).to_string(index=False))
    print(f"[write] {a.out}, {OUT}/a2_transfer_table.csv")


if __name__ == "__main__":
    main()
