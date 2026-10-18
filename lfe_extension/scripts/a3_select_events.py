#!/usr/bin/env python
"""Track A step 3 -- build the EVENT LISTS the iterated templates are stacked on (Merlin plan-review step 3, selection).

Selection = detections of the gated selection pass at T_core on TREMOR days only, HF-veto applied at the frozen cut,
cross-family assigned by raw weighted cc among families that share >= 50 % of their stations, already declustered
within the family at the template length by the matched filter itself.
  tremor day = >= 30 PNSN tremor detections in the family box that UTC day, OR >= 100 Bostock events (an ETS day).
Two lists are written because the PGC site carries > 40 % of most families' weight (53 % of family-days):
  a3_events_gated.parquet   -- from the full 15-station gated pass; used to re-stack every station OUTSIDE the PGC site
  a3_events_losoPGC.parquet -- from the pass with the PGC site excluded from the network sum; used to re-stack
                               PGC, B009, B010, B011, so those templates are never built on events those sensors selected.
Also reported: how many events each family gets (Merlin's fallback trigger is < 50 per family per SSE)."""
from __future__ import annotations
import argparse
import numpy as np, pandas as pd

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
OUT = "/home/jovyan/tremorferometry/lfe_extension/results"
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"


def cross_family_keep(t, v, fam, share_ok, w=2.0):
    o = np.argsort(t, kind="stable"); ts = t[o]; vs = v[o]; fs = fam[o]; n = len(ts)
    lo = np.searchsorted(ts, ts - w, "left"); hi = np.searchsorted(ts, ts + w, "right")
    keep = np.ones(n, bool)
    for i in range(n):
        for j in range(lo[i], hi[i]):
            if j == i or not share_ok[fs[i], fs[j]]:
                continue
            if vs[j] > vs[i] or (vs[j] == vs[i] and j < i):
                keep[i] = False
                break
    out = np.zeros(n, bool); out[o] = keep
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gated", default="a3_sel_gated"); ap.add_argument("--loso", default="a3_sel_loso_PGC")
    ap.add_argument("--L", type=float, default=12.0); ap.add_argument("--t-core", type=float, default=8.75)
    ap.add_argument("--hf-gated", default=""); ap.add_argument("--hf-loso", default="")
    ap.add_argument("--hf-cut", type=float, default=3.6384)
    ap.add_argument("--min-tremor", type=int, default=30)
    a = ap.parse_args()
    gt = pd.read_parquet(f"{IN}/e5r_gt_bostock2012.parquet"); gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    per_day = gt.groupby("doy").size(); ETS = set(per_day[per_day >= 100].index)
    F = pd.read_csv(f"{IN}/e5r_families_2012.csv")
    tr = pd.read_csv(PNSN); tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    box = tr[(tr.time.dt.year == 2012) & tr.lat.between(F.lat.min() - .2, F.lat.max() + .2) & tr.lon.between(F.lon.min() - .2, F.lon.max() + .2)]
    tcnt = box.groupby(box.time.dt.dayofyear).size()
    TREM = set(int(d) for d in tcnt[tcnt >= a.min_tremor].index) | set(int(d) for d in ETS)
    qc = pd.read_csv(f"{IN}/e5r_template_qc_2012.csv")
    sets = qc[qc.usable & qc.template].groupby("bkey").sta.apply(set).to_dict()
    nb = int(max(sets)) + 2
    share = np.zeros((nb, nb), bool)
    for fa_, A_ in sets.items():
        for fb_, B_ in sets.items():
            share[fa_, fb_] = len(A_ & B_) / max(min(len(A_), len(B_)), 1) >= 0.5
    print(f"[days] {len(TREM)} tremor days in 2012 doy 228-279: {sorted(TREM)}")
    for tag, hf, out in ((a.gated, a.hf_gated, "a3_events_gated"), (a.loso, a.hf_loso, "a3_events_losoPGC")):
        d = pd.read_parquet(f"{OUT}/e5r_detections_{tag}.parquet"); d = d[d.L == a.L]
        n0 = len(d)
        d = d[d.doy.isin(TREM) & (d.net_mad >= a.t_core)].reset_index(drop=True)
        n1 = len(d)
        if hf:
            H = pd.read_parquet(hf); H = H[H.L == a.L][["doy", "fam", "t", "hfr_PGC"]]
            d = d.merge(H, on=["doy", "fam", "t"], how="left")
            v = d.hfr_PGC.to_numpy(); d = d[~(np.isfinite(v) & (v > a.hf_cut))].reset_index(drop=True)
        n2 = len(d)
        keep = cross_family_keep(d.t.to_numpy(), d.net_cc.to_numpy().astype(float), d.fam.to_numpy().astype(int), share, 2.0)
        d = d[keep].reset_index(drop=True)
        E = pd.DataFrame(dict(bkey=d.fam.astype(int).to_numpy(), t_origin=d.t.to_numpy(),
                              net_mad=d.net_mad.to_numpy(), n_sites=d.n_site.to_numpy(), doy=d.doy.to_numpy()))
        E["half"] = "build"; E["guard"] = False        # every selected detection is used (--half all downstream)
        E.to_parquet(f"{OUT}/{out}.parquet", index=False)
        pf = E.groupby("bkey").size()
        print(f"[{tag}] {n0} detections -> tremor days & T>={a.t_core}: {n1} -> HF veto: {n2} -> cross-family: {len(E)}; "
              f"{E.bkey.nunique()} families; events/family p10/50/90 = {np.percentile(pf, [10,50,90]).round(0)}; "
              f"families with >= 30: {int((pf>=30).sum())}, >= 50: {int((pf>=50).sum())} -> {OUT}/{out}.parquet")


if __name__ == "__main__":
    main()
