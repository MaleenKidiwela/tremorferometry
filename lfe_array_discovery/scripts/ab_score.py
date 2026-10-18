#!/usr/bin/env python
"""STAGE 8 -- score the Royer reproduction.  Existence-and-count, ranked, never gated.

Royer & Bostock released no catalogue, so there is no one-to-one match to compute.  The test is the
one the plan states for a region with no released event list: do we produce located, repeating
families above the whole-chain null, inside the published region and depth band, and how do our
stage-by-stage counts stand beside theirs?

REFERENCE LINES (their numbers, never gates):
  initial detections   NW 4915   SW 3306        (sec.3.1)
  initial templates    NW  224   SW   70        (sec.3.2)
  final templates      NW  122   SW   54        (sec.3.3)
  depth range          NW 31-46 km  SW 32-49 km (sec.5)
  Hyp2000 sigma_h/z    NW 2.0/1.3 km  SW 3.7/3.5 km
  hypoDD retention     NW 77 %   SW 72 %

RULES (memory `rank-dont-discard`): no criterion is ever applied as a cut.  Every family keeps its
row in `families_ranked.csv`; every criterion is a column; the null receives the identical table; the
evidence is the real-versus-null count as a function of each criterion (`separation_curves.csv`) plus
a ranking.  If a count comes in far below Royer's, that is a bug to find, not a finding
(memory `reproduce-known-answer-first`).

  python scripts/ab_score.py --sub NW
-> results/ab_royer/<SUB>/score/stage_counts.csv        real | null_k rows, our value | published
-> results/ab_royer/<SUB>/score/families_ranked.csv     every family, every criterion a column
-> results/ab_royer/<SUB>/score/separation_curves.csv   real-vs-null counts vs each criterion
-> results/ab_royer/<SUB>/score/SUMMARY.md
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_common as A                                                       # noqa: E402

CRITERIA = ["n_det", "n_sta", "nmad_med", "cc_med", "n_days", "boot_km", "dep", "dd_dep",
            "erh", "erz"]


def _roots(sub):
    out = [("real", f"{A.RES}/{sub}")]
    for d in sorted(glob.glob(f"{A.RES}/{sub}/null*")):
        out.append((os.path.basename(d), d))
    return out


def stage_counts(sub):
    pub = A.SUBARRAY[sub]
    rows = []
    for tag, root in _roots(sub):
        r = dict(pop=tag, subarray=sub)
        st = f"{root}/auto_window/_stats.csv"
        if os.path.exists(st):
            S = pd.read_csv(st)
            r["n_hours"] = int(len(S))
            r["n_pairs"] = int(S.n_det.sum()) if "n_det" in S else np.nan
            r["thr_med"] = float(S.thr.median()) if "thr" in S else np.nan
            r["nmad_in_sd_med"] = float(S.nmad_in_sd.median()) if "nmad_in_sd" in S else np.nan
        cj = f"{root}/clusters/summary.json"
        if os.path.exists(cj):
            j = json.load(open(cj))
            r["n_detections"] = j["n_detections"]
            r["n_initial_templates"] = j["n_initial_templates"]
            r["cluster_linkage"] = j["headline"]["linkage"]
            r["cluster_cut"] = j["headline"]["cut"]
            r["n_clusters_nearpub"] = j["nearest_published"]["n_clusters"]
            r["cluster_min_size"] = j["headline"].get("min_size")
            r["cluster_rule"] = j["headline"].get("how")
            r["n_pairs_strict"] = j.get("n_pairs_strict")
            r["n_pairs_pairnms"] = j.get("n_pairs_pairnms")
        gj = f"{root}/grow/summary.json"
        if os.path.exists(gj):
            j = json.load(open(gj))
            r["n_final_templates"] = j["n_final_templates"]
            r["n_iterations"] = len(j["iterations"])
            if j["iterations"]:
                r["n_grow_detections"] = j["iterations"][-1]["n_detections"]
        hs = f"{root}/hypo/loc_summary.csv"
        if os.path.exists(hs):
            H = pd.read_csv(hs)
            h = H[H.headline == True]                                          # noqa: E712
            if len(h):
                h = h.iloc[0]
                for c in ["n_hyp", "hyp_dep_med", "hyp_frac_in_band", "hyp_erh_mean",
                          "hyp_erz_mean", "n_dd", "dd_keep_frac", "dd_dep_med",
                          "dd_frac_in_band", "boot_km_med"]:
                    r[c] = h[c]
        ds = f"{root}/locate/depth_summary.csv"
        if os.path.exists(ds):
            D = pd.read_csv(ds)
            h = D[D.variant == "L9_xcdiff"]
            if len(h):
                h = h.iloc[0]
                r["L9_n_located"] = int(h.n_located)            # secondary locator
                r["L9_dep_med"] = float(h.dep_med)
                r["L9_frac_in_published_band"] = float(h.frac_in_published_band)
                r["L9_rail"] = float(h.rail)
                r["L9_boot_km_med"] = float(h.boot_km_med)
        rows.append(r)
    T = pd.DataFrame(rows)
    T["published_n_detections"] = pub["n_det_pub"]
    T["published_n_initial_templates"] = pub["n_tpl0_pub"]
    T["published_n_final_templates"] = pub["n_tpl_pub"]
    T["published_depth_band"] = f"{pub['depth_pub'][0]}-{pub['depth_pub'][1]}"
    return T


def families(root, sub, variant="L9_xcdiff"):
    """One row per grown template, every criterion a COLUMN.  Nothing is filtered."""
    gd = sorted(glob.glob(f"{root}/grow/det_it*.parquet"))
    if not gd:
        return pd.DataFrame()
    last = max(gd, key=lambda p: int(p.split("it")[-1].split(".")[0]))
    D = pd.read_parquet(last)
    F = (D.groupby("grp")
          .agg(n_det=("cc", "size"), cc_med=("cc", "median"), cc_max=("cc", "max"),
               nmad_med=("nmad", "median"), nmad_max=("nmad", "max"),
               n_sta=("n_sta", "median"), n_days=("date", "nunique"),
               t_first=("t", "min"), t_last=("t", "max")).reset_index())
    F["src_iter"] = int(last.split("it")[-1].split(".")[0])
    HL = f"{root}/hypo/loc_headline.csv"                   # headline Hyp2000 + hypoDD
    if os.path.exists(HL):
        Hc = pd.read_csv(HL)
        keep = [c for c in ["grp", "lat", "lon", "dep", "erh", "erz", "rms", "num", "gap",
                            "boot_km", "boot_dz", "dd_lat", "dd_lon", "dd_dep", "dd_kept"]
                if c in Hc]
        F = F.merge(Hc[keep], on="grp", how="left")
    L = f"{root}/locate/loc_{variant}.csv"
    if os.path.exists(L):                                      # secondary locator, prefixed
        Lc = pd.read_csv(L)
        Lc = Lc[[c for c in ["grp", "lat", "lon", "dep", "boot_km", "edge", "located"] if c in Lc]]
        F = F.merge(Lc.rename(columns={c: f"L9_{c}" for c in Lc if c != "grp"}), on="grp",
                    how="left")
    lo, hi = A.SUBARRAY[sub]["depth_pub"]
    if "dep" in F:
        F["in_published_depth_band"] = F.dep.between(lo, hi)
    cs = f"{root}/clusters/cluster_stacks.csv"
    if os.path.exists(cs):
        C = pd.read_csv(cs).rename(columns={"cid": "grp", "n_mem": "n_cluster_members"})
        F = F.merge(C[["grp", "n_cluster_members"]], on="grp", how="left")
    return F


def _md(df):
    """SMOKE-TEST FIX: DataFrame.to_markdown needs `tabulate`, which is not installed in this pod
    (ImportError killed the scorecard).  Fall back to a fixed-width block."""
    try:
        return df.to_markdown(index=False)
    except ImportError:
        return "```\n" + df.to_string(index=False) + "\n```"


def separation(real, nulls):
    """Real-vs-null counts as a function of each criterion -- the separation curve, not a cut."""
    rows = []
    for c in CRITERIA:
        if c not in real or real[c].dropna().empty:
            continue
        v = real[c].dropna().to_numpy(float)
        qs = np.unique(np.percentile(v, np.arange(0, 100, 5)))
        for q in qs:
            n_r = int((real[c] >= q).sum())
            n_n = [int((N[c] >= q).sum()) if c in N else 0 for N in nulls]
            rows.append(dict(criterion=c, threshold=float(q), n_real=n_r,
                             n_null_mean=float(np.mean(n_n)) if n_n else np.nan,
                             n_null_max=int(np.max(n_n)) if n_n else 0,
                             excess=n_r - (float(np.mean(n_n)) if n_n else 0.0)))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", default="NW,SW,NVI")
    ap.add_argument("--variant", default="L9_xcdiff")
    a = ap.parse_args()
    for sub in a.sub.split(","):
        out = f"{A.RES}/{sub}/score"
        os.makedirs(out, exist_ok=True)
        T = stage_counts(sub)
        T.to_csv(f"{out}/stage_counts.csv", index=False)
        rs = _roots(sub)
        real = families(rs[0][1], sub, a.variant)
        nulls = [families(r, sub, a.variant) for _, r in rs[1:]]
        for tag, F in [("real", real)] + list(zip([t for t, _ in rs[1:]], nulls)):
            if len(F):
                F.insert(0, "pop", tag)
        allF = pd.concat([f for f in [real] + nulls if len(f)], ignore_index=True) \
            if len(real) or any(len(f) for f in nulls) else pd.DataFrame()
        if len(allF):
            allF.to_csv(f"{out}/families_ranked.csv", index=False)
        S = separation(real, [f for f in nulls if len(f)]) if len(real) else pd.DataFrame()
        if len(S):
            S.to_csv(f"{out}/separation_curves.csv", index=False)
        pub = A.SUBARRAY[sub]
        with open(f"{out}/SUMMARY.md", "w") as fh:
            fh.write(f"# Royer & Bostock (2014) {sub} Washington -- reproduction scorecard\n\n")
            fh.write("Existence-and-count reproduction: they released no catalogue, so their "
                     "numbers are REFERENCE LINES and never gates.  Every family keeps its row; "
                     "every criterion is a column; the null gets the identical table.\n\n")
            fh.write("## Stage counts, real beside every null realisation\n\n")
            fh.write(_md(T) + "\n\n")
            fh.write(f"Royer: {pub['n_det_pub']} initial detections -> {pub['n_tpl0_pub']} initial "
                     f"templates -> {pub['n_tpl_pub']} final templates; depths "
                     f"{pub['depth_pub'][0]}-{pub['depth_pub'][1]} km; Hyp2000 sigma "
                     f"{pub['sig_h_pub']}/{pub['sig_z_pub']} km; hypoDD keeps "
                     f"{pub['dd_keep_pub']:.0%}.\n\n")
            if len(S):
                fh.write("## Separation curves (real vs the mean of the nulls)\n\n")
                fh.write(_md(S) + "\n")
        A.log(f"[score {sub}] wrote {out}/SUMMARY.md")


if __name__ == "__main__":
    main()
