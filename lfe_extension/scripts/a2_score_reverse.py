#!/usr/bin/env python
"""Track A step 2(c) -- the REVERSE frozen transfer: 2010-built templates applied to 2012, scored against the
2012-built templates on the same 9 stations (the step-2(a) ablation) at the frozen thresholds.  With 2010 templates
the whole 2012 catalog is held out, so recall is reported on ALL 2012 events and on the test half alone."""
from __future__ import annotations
import argparse, json
import numpy as np, pandas as pd

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
OUT = "/home/jovyan/tremorferometry/lfe_extension/results"
QUIET = list(range(228, 243))


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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fwd", default="a2_abl2012_7p"); ap.add_argument("--rev", default="a2_rev2012")
    ap.add_argument("--L", type=float, default=12.0)
    ap.add_argument("--t-ext", type=float, default=8.25); ap.add_argument("--t-core", type=float, default=8.75)
    a = ap.parse_args()
    gt = pd.read_parquet(f"{IN}/e5r_gt_bostock2012.parquet").reset_index(drop=True)
    ev_all = gt.reset_index(drop=True)
    ev_test = gt[(gt.half == "test") & (~gt.guard)].reset_index(drop=True)
    rows = []
    D = {}
    for tag in (a.fwd, a.rev):
        d = pd.read_parquet(f"{OUT}/e5r_detections_{tag}.parquet"); D[tag] = d[d.L == a.L].reset_index(drop=True)
        b = pd.read_csv(f"{OUT}/e5_budget_{tag}.csv"); b = b[b.L == a.L]
        T0 = sorted(b["T"].unique())[0]
        fd = b[b["T"] == T0][["doy", "fam"]].drop_duplicates()
        D[tag + "_nq"] = int(fd.doy.isin(QUIET).sum())
        D[tag + "_nfam"] = int(d.fam.nunique())
    for nm, T in (("T_ext", a.t_ext), ("T_core", a.t_core)):
        r = dict(threshold=nm, T=T)
        for lbl, tag in (("fwd_2012tpl", a.fwd), ("rev_2010tpl", a.rev)):
            d = D[tag]
            r[f"recall_all_{lbl}"] = float(recall(d, ev_all, T).mean())
            r[f"recall_test_{lbl}"] = float(recall(d, ev_test, T).mean())
            r[f"quiet_fa_{lbl}"] = len(d[(d.net_mad >= T) & d.doy.isin(QUIET)]) / max(D[tag + "_nq"], 1)
            r[f"n_catalog_{lbl}"] = int((d.net_mad >= T).sum())
            r[f"precision_own_{lbl}"] = float(d[d.net_mad >= T].match_own.mean())
            r[f"n_families_{lbl}"] = D[tag + "_nfam"]
        r["recall_ratio_rev_over_fwd_all"] = r["recall_all_rev_2010tpl"] / max(r["recall_all_fwd_2012tpl"], 1e-9)
        r["fa_ratio_rev_over_fwd"] = r["quiet_fa_rev_2010tpl"] / max(r["quiet_fa_fwd_2012tpl"], 1e-9)
        rows.append(r)
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/a2_reverse_table.csv", index=False)
    json.dump(rows, open(f"{OUT}/a2_reverse_report.json", "w"), indent=1, default=float)
    pd.set_option("display.width", 260)
    print(R.round(4).to_string(index=False))
    print(f"[note] the 2010-built set covers {rows[0]['n_families_rev_2010tpl']} families vs "
          f"{rows[0]['n_families_fwd_2012tpl']} for the 2012-built set (only the 90 families with an accepted 2010 offset "
          f"and >= 30 2010 events could be built), so recall over ALL 2012 events is bounded by that coverage.")
    print(f"[write] {OUT}/a2_reverse_table.csv")


if __name__ == "__main__":
    main()
