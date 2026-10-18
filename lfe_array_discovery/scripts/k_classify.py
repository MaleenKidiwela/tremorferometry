#!/usr/bin/env python
"""Classify every grown template of the 2005 run: kept family / hidden family / candidate new
family / unlocatable burst."""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import k_common as K                                                        # noqa: E402
from k_dedupe2 import load_tpl, tpl_cc                                      # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loc", required=True)
    ap.add_argument("--score-kept", required=True)
    ap.add_argument("--score-hidden", default="")
    ap.add_argument("--tpl", required=True)
    ap.add_argument("--tpl-kept-ceiling", default="")
    ap.add_argument("--tpl-hidden-ceiling", default="")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    L = pd.read_csv(a.loc)
    SK = pd.read_csv(a.score_kept).set_index("grp")
    L["kept_bkey"] = [int(SK.bkey.get(g, -1)) if SK.frac.get(g, 0) >= 0.30 else -1 for g in L.grp]
    L["hid_bkey"] = -1
    if a.score_hidden and os.path.exists(a.score_hidden):
        SH = pd.read_csv(a.score_hidden).set_index("grp")
        L["hid_bkey"] = [int(SH.bkey.get(g, -1)) if SH.frac.get(g, 0) >= 0.30 else -1
                         for g in L.grp]
    T = load_tpl(a.tpl)
    CC = {}
    for tag, path in (("kept", a.tpl_kept_ceiling), ("hidden", a.tpl_hidden_ceiling)):
        if not path or not os.path.exists(path):
            continue
        C = load_tpl(path)
        best = {}
        for g in T:
            bb, bc = -1, -1.0
            for b, tb in C.items():
                cc, n = tpl_cc(T[g], tb)
                if np.isfinite(cc) and cc > bc:
                    bb, bc = b, cc
            best[g] = (bb, bc)
        CC[tag] = best
        L[f"cc_{tag}_best"] = [best.get(g, (-1, np.nan))[1] for g in L.grp]
        L[f"cc_{tag}_bkey"] = [best.get(g, (-1, np.nan))[0] for g in L.grp]
    def cls(r):
        if r.kept_bkey >= 0:
            return "kept_family"
        if r.hid_bkey >= 0:
            return "hidden_family"
        if not r.located:
            return "unlocatable"
        if bool(r.ctrl_pass) and np.isfinite(r.dz) and abs(r.dz) <= 8 and r.n_days >= 3:
            return "candidate_new"
        return "unlocatable"
    L["cls"] = [cls(r) for r in L.itertuples()]
    L.to_csv(a.out, index=False)
    K.log("\n" + str(L.cls.value_counts()))
    for c in ("kept_family", "hidden_family", "candidate_new"):
        sub = L[L.cls == c]
        K.log(f"{c}: {len(sub)} templates, "
              f"{sub.kept_bkey.nunique() if c=='kept_family' else (sub.hid_bkey.nunique() if c=='hidden_family' else len(sub))} distinct")
    cand = L[L.cls == "candidate_new"].sort_values("n_det", ascending=False)
    if len(cand):
        cols = [c for c in ["grp", "n_det", "n_days", "n_sta", "lat", "lon", "dep", "slab2",
                            "dz", "rms", "rms_shuf_p05", "cc_kept_best", "cc_hidden_best"]
                if c in cand.columns]
        K.log("\nCANDIDATE NEW FAMILIES (top 25)\n" + str(cand[cols].head(25).round(3)))


if __name__ == "__main__":
    main()
