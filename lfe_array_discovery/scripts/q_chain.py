#!/usr/bin/env python
"""FINAL SPEC -- the two SURVIVAL ORDERS, one implementation.  Prefix `q_`.

Order "frozen" (A) -- exactly the 1c' chain, unchanged except F5 (survival counts ETS-window days
only, behind W_SURV_ETS_ONLY):

    handoff -> quiet-day T_seed -> growth -> template-scramble twin -> SURVIVAL GATE
            -> iterate 1 -> quiet -> grow -> survive -> iterate 2 -> quiet -> grow -> survive
            -> dedupe -> locate/Gates 0-2

Order "staged" (B) -- Merlin's conditional alternative (PLAN sec.12 step 3), adopted ONLY if the
truth-seeded bootstrap test selects it.  EVERY threshold, clause and bar is the frozen one; only
the ORDER of frozen operations changes, and the twin goes through it identically:

    handoff -> quiet-day T_seed -> pass-0 growth -> CONTINUATION FILTER (a compute filter, not a
    claim: >= 5 independent hits on >= 2 ETS days AND >= 2 x the template-scramble twin)
            -> v_iterate (use_all, min_indep 5) -> growth -> iterate -> growth
            -> THE FROZEN SURVIVAL GATE with every clause, ETS days only
            -> dedupe -> locate/Gates 0-2

Both orders are driven by the same subprocess calls into w_grow / w_iterate / w_dedupe /
w_locate_gates, every stage guarded by its own output file so a relaunched queue resumes.
"""
from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import w_common as W                                                        # noqa: E402
import q_common as Q                                                        # noqa: E402

PY = sys.executable
G = f"{HERE}/w_grow.py"
IT = f"{HERE}/w_iterate.py"
DD = f"{HERE}/w_dedupe.py"
LG = f"{HERE}/w_locate_gates.py"
STUDY = os.path.dirname(HERE)


def sh(cmd, guard="", env=None):
    if guard and os.path.exists(guard):
        W.log(f"[chain] skip {os.path.basename(cmd[1])} (guard {os.path.basename(guard)})")
        return
    W.log("$ " + " ".join(os.path.basename(c) for c in cmd[:3]) + " ...")
    e = None
    if env:
        e = dict(os.environ)
        e.update({k: str(v) for k, v in env.items()})
    r = subprocess.run(cmd, cwd=STUDY, env=e)
    if r.returncode:
        raise SystemExit(f"chain stage failed ({r.returncode}): {' '.join(cmd)}")
    if guard and not os.path.exists(guard):
        raise SystemExit(f"chain stage produced no {guard}: {' '.join(cmd)}")


def _shift_args(js):
    return ["--shift", js] if js else []


# ------------------------------------------------------------------ continuation filter ---------
def continuation(d, tag, out_csv, twin=False):
    """Merlin's compute filter (order B).  Writes a survivor-shaped CSV holding ONLY the templates
    that continue, so `w_iterate --use-all --min-indep 5` selects exactly them.

    real chain: >= S3_CONT_MIN_INDEP independent hits on >= S3_CONT_MIN_DAYS ETS days AND
                n_indep >= S3_CONT_SCRAMBLE_RATIO x the template-scramble twin;
    twin chain: the count clause only (a twin has no twin of its own).
    """
    SV = W.read_csv_safe(f"{d}/{tag}_survival.csv", W.SURV_COLS)
    ip = f"{d}/{tag}_indep.parquet"
    H = (pd.read_parquet(ip) if os.path.exists(ip)
         else pd.DataFrame(columns=["tid", "t", "date"]))
    ETS = W.ets_date_set()
    rows = []
    for r in SV.itertuples():
        tid = str(r.tid)
        g = H[H.tid.astype(str) == tid] if len(H) else H
        ge = g[g.date.astype(str).isin(ETS)] if len(g) else g
        n_i = int(len(ge))
        n_d = int(ge.date.nunique()) if n_i else 0
        nm = float(getattr(r, "n_scramble", np.nan))
        ok_cnt = (n_i >= Q.S3_CONT_MIN_INDEP and n_d >= Q.S3_CONT_MIN_DAYS)
        if twin:
            ok_scr, ratio = True, np.nan
        else:
            ratio = n_i / max(nm, 1.0) if np.isfinite(nm) else np.inf
            ok_scr = bool(ratio >= Q.S3_CONT_SCRAMBLE_RATIO)
        rows.append(dict(tpl=r.tpl, tid=r.tid, n_hits=r.n_hits, n_indep=n_i, n_days=n_d,
                         frac_near_seed=r.frac_near_seed, n_scramble=nm,
                         scramble_ratio=ratio, survives=bool(ok_cnt and ok_scr)))
    R = (pd.DataFrame(rows, columns=W.SURV_COLS) if rows
         else pd.DataFrame(columns=W.SURV_COLS))
    R.to_csv(out_csv.replace(".csv", "_all.csv"), index=False)
    K = R[R.survives.astype(bool)] if len(R) else R
    K.to_csv(out_csv, index=False)
    W.log(f"[continuation{'/twin' if twin else ''}] {len(R)} templates -> {len(K)} continue "
          f"(>= {Q.S3_CONT_MIN_INDEP} indep on >= {Q.S3_CONT_MIN_DAYS} ETS days"
          f"{'' if twin else ' AND >= %gx the twin' % Q.S3_CONT_SCRAMBLE_RATIO})")
    return K


# ------------------------------------------------------------------ the chain -------------------
def run(d, seeds_csv, clusters_parquet, order="frozen", shift_json="", gates=True):
    """Run one whole chain in directory `d`.  Returns a dict of counts."""
    os.makedirs(d, exist_ok=True)
    js = shift_json or ""
    sa = _shift_args(js)

    sh([PY, G, "handoff", "--clusters", clusters_parquet, "--seeds", seeds_csv,
        "--out", f"{d}/bank0.npz"] + sa, f"{d}/bank0.npz")
    sh([PY, G, "scramble", "--clusters", clusters_parquet, "--seeds", seeds_csv,
        "--out", f"{d}/bank0scr.npz"] + sa, f"{d}/bank0scr.npz")
    sh([PY, G, "quiet", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0"] + sa,
       f"{d}/p0_thresholds.csv")
    sh([PY, G, "grow", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0",
        "--thr", f"{d}/p0_thresholds.csv"] + sa, f"{d}/p0_hits.parquet")
    # pass-0 twin: grown at the REAL pass-0 thresholds, the frozen convention (same seed, same
    # template index).  From iteration 1 the twin bank is a different template set and gets its own
    # quiet-day T.
    sh([PY, G, "grow", "--bank", f"{d}/bank0scr.npz", "--out", f"{d}/p0scr",
        "--thr", f"{d}/p0_thresholds.csv"] + sa, f"{d}/p0scr_hits.parquet")
    sh([PY, G, "survive", "--bank", f"{d}/bank0.npz", "--out", f"{d}/p0", "--seeds", seeds_csv,
        "--scramble-hits", f"{d}/p0scr_hits.parquet", "--scramble-bank", f"{d}/bank0scr.npz"],
       f"{d}/p0_survival.csv")

    if order == "frozen":
        prev, tag = f"{d}/bank0.npz", "p0"
        for it in (1, 2):
            sh([PY, IT, "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
                "--survivors", f"{d}/p0_survival.csv", "--out-bank", f"{d}/bank{it}.npz",
                "--out-diag", f"{d}/iter{it}_diag.csv"] + sa, f"{d}/bank{it}.npz")
            sh([PY, G, "quiet", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}"] + sa,
               f"{d}/p{it}_thresholds.csv")
            sh([PY, G, "grow", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}",
                "--thr", f"{d}/p{it}_thresholds.csv"] + sa, f"{d}/p{it}_hits.parquet")
            sh([PY, G, "survive", "--bank", f"{d}/bank{it}.npz", "--out", f"{d}/p{it}",
                "--seeds", seeds_csv], f"{d}/p{it}_survival.csv")
            prev, tag = f"{d}/bank{it}.npz", f"p{it}"
        final_bank, final_tag = prev, tag
    else:
        # ---- STAGED (order B): no gate at pass 0, only the compute filter --------------------
        sh([PY, G, "survive", "--bank", f"{d}/bank0scr.npz", "--out", f"{d}/p0scr",
            "--seeds", seeds_csv], f"{d}/p0scr_survival.csv")
        if not os.path.exists(f"{d}/p0_cont.csv"):
            continuation(d, "p0", f"{d}/p0_cont.csv", twin=False)
        if not os.path.exists(f"{d}/p0scr_cont.csv"):
            continuation(d, "p0scr", f"{d}/p0scr_cont.csv", twin=True)
        for chain, cont0 in (("", f"{d}/p0_cont.csv"), ("scr", f"{d}/p0scr_cont.csv")):
            prev = f"{d}/bank0{chain}.npz"
            src = cont0
            tag = f"p0{chain}"
            for it in (1, 2):
                sh([PY, IT, "--bank", prev, "--indep", f"{d}/{tag}_indep.parquet",
                    "--survivors", src, "--use-all", "--min-indep", str(Q.S3_ITER_MIN_INDEP),
                    "--out-bank", f"{d}/bank{it}{chain}.npz",
                    "--out-diag", f"{d}/iter{it}{chain}_diag.csv"] + sa,
                   f"{d}/bank{it}{chain}.npz")
                sh([PY, G, "quiet", "--bank", f"{d}/bank{it}{chain}.npz",
                    "--out", f"{d}/p{it}{chain}"] + sa, f"{d}/p{it}{chain}_thresholds.csv")
                sh([PY, G, "grow", "--bank", f"{d}/bank{it}{chain}.npz",
                    "--out", f"{d}/p{it}{chain}",
                    "--thr", f"{d}/p{it}{chain}_thresholds.csv"] + sa,
                   f"{d}/p{it}{chain}_hits.parquet")
                if it == 1 or chain == "scr":
                    # bookkeeping pass: independent hits for the next iteration / the twin count
                    sh([PY, G, "survive", "--bank", f"{d}/bank{it}{chain}.npz",
                        "--out", f"{d}/p{it}{chain}", "--seeds", seeds_csv],
                       f"{d}/p{it}{chain}_survival.csv")
                    src = f"{d}/p{it}{chain}_survival.csv"
                prev, tag = f"{d}/bank{it}{chain}.npz", f"p{it}{chain}"
        # ---- THE FROZEN SURVIVAL GATE, every clause, ETS days only, against the staged twin ---
        # a template whose twin did NOT continue has n_scramble = 0, not NaN: the twin produced
        # nothing, which is evidence FOR the template, not a reason to veto it.
        sh([PY, G, "survive", "--bank", f"{d}/bank2.npz", "--out", f"{d}/p2",
            "--seeds", seeds_csv, "--scramble-hits", f"{d}/p2scr_hits.parquet",
            "--scramble-bank", f"{d}/bank2scr.npz"], f"{d}/p2_survival.csv",
           env=dict(W_SCRAMBLE_MISSING_ZERO="1"))
        final_bank, final_tag = f"{d}/bank2.npz", "p2"

    sh([PY, DD, "--bank", final_bank, "--indep", f"{d}/{final_tag}_indep.parquet",
        "--out", f"{d}/dedupe.csv"], f"{d}/dedupe.csv")
    cmd = [PY, LG, "--bank", final_bank, "--indep", f"{d}/{final_tag}_indep.parquet",
           "--dedupe", f"{d}/dedupe.csv", "--survival", f"{d}/{final_tag}_survival.csv",
           "--out", f"{d}/lg"] + sa
    if not gates:
        cmd.append("--no-gates")
    sh(cmd, f"{d}/lg.json")

    j = Q.jload(f"{d}/lg.json", {}) or {}
    out = dict(order=order, dir=d,
               n_seeds=int(len(W.read_csv_safe(seeds_csv))),
               n_survivors_p0=Q.n_surv(f"{d}/p0_survival.csv"),
               n_survivors_final=Q.n_surv(f"{d}/{final_tag}_survival.csv"),
               n_located=j.get("n_located"), n_gate_pass=j.get("n_gate_pass"),
               final_bank=final_bank, final_tag=final_tag)
    if order == "staged":
        out["n_continue"] = int(len(W.read_csv_safe(f"{d}/p0_cont.csv")))
        out["n_continue_twin"] = int(len(W.read_csv_safe(f"{d}/p0scr_cont.csv")))
        out["n_twin_final"] = int(len(W.read_csv_safe(f"{d}/p2scr_survival.csv")))
    Q.jdump(out, f"{d}/chain.json")
    W.log(f"[chain/{order}] {d}: seeds {out['n_seeds']}, survivors p0 {out['n_survivors_p0']} / "
          f"final {out['n_survivors_final']}, located {out['n_located']}, "
          f"gated {out['n_gate_pass']}")
    return out
