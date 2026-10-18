#!/usr/bin/env python
"""s11_wa_kerr -- POST HOC (written 2026-10-12 ~3:55 pm PT, after the WA v2 coherent kernel-error injections failed:
two5_40-generated A Q2 recovery 0.656, B Q1 leak 0.411; ell80 B leak 0.268). Tests whether a different use of the
kernels makes Q1/Q2 robust to the kernel model. Same waveform-level injections as s11_v2.run_injections (A deep step +0.3 %
at 2015-01-01, B shallow annual 0.1 % peak day 233; generated with ell40 / two5_40 / ell80 coefficients), inverted with:
  M0 nominal      a_s, a_d of ell40 (= the pre-registered headline; reproduces v2_injection.csv)
  M1 site-const   a_s replaced by the mean ell40 a_s of the site's cells (no within-site shallow pattern)
  M2 kernel-avg   a_s, a_d = mean over ell20/ell40/ell80/two5_40 per series
  M3 shallow=1    a_s = 1 (shallow term per site is kernel-free; its scale is absorbed by S)
Also each model on the natural data (Q1, Q2). Columns only: no rule of PLAN §5 is replaced by this.
    S11_REGION=WA_NW python s11_wa_kerr.py --workers 8      -> tables/kerr_models.csv
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_v2 as V
import s11_paths as P

GEN = ("ell40", "two5_40", "ell80")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=8); a = ap.parse_args()
    log = C.Log(f"{C.LOGD}/kerr.log"); t0 = time.time()
    import multiprocessing as mp
    cells = pd.read_csv(f"{C.TB}/v2_cells.csv")
    key = (cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch).to_numpy()
    Kall = pd.read_parquet(f"{C.TB}/kernel_coef.parquet"); Kall = Kall[(Kall.comp == "3c") & (Kall.win == C.HEAD)]
    K = {v: Kall[Kall.variant == v].set_index("series") for v in ("ell20", "ell40", "ell80", "two5_40")}
    unit = P.unit_map(); need = set(key); vals = {}
    for ns in C.STATIONS:
        if not os.path.exists(f"{C.TB}/store/{ns}.npz"):
            continue
        V._G.clear(); V._G.update(store=P.load_store(ns), unit=unit, K=K)
        ks = [k for k in V._G["store"] if k in need]
        if not ks:
            continue
        with mp.get_context("fork").Pool(a.workers) as pool:
            for k, o in pool.imap_unordered(V._inj_series, ks):
                vals[k] = o

    def dvec(run):
        d = np.full(len(cells), np.nan)
        for i, (k, e) in enumerate(zip(key, cells.episode)):
            q = vals.get(k, {}).get(run)
            if q is None or not len(q):
                continue
            m = q.episode == e
            if m.any():
                d[i] = -100 * float(q[f"eps3c_{C.HEAD}"][m].iloc[0])
        return d
    dn = dvec("nat")
    models = {}
    c0 = cells.copy(); models["M0_nominal"] = c0
    c1 = cells.copy(); c1["a_s"] = c1.groupby("site").a_s.transform("mean"); models["M1_site_const_as"] = c1
    c2 = cells.copy()
    for col in ("a_s", "a_d"):
        c2[col] = np.nanmean(np.column_stack([pd.Series(key).map(K[v][col]).to_numpy() for v in K]), 1)
    models["M2_kernel_avg"] = c2
    c3 = cells.copy(); c3["a_s"] = 1.0; models["M3_shallow_1"] = c3
    rows = []
    for mn, cm in models.items():
        q = V.Q(V.fit(cm, "smooth"))
        r = dict(model=mn, Q1=q["Q1"], Q1_sig=q["Q1_sig"], Q1_peak=q["Q1_peak"], Q2=q["Q2"], Q2_sig=q["Q2_sig"], chi2=q["chi2"])
        for g in GEN:
            for run in ("A", "B"):
                di = dvec(f"{run}_{g}"); ok = np.isfinite(di) & np.isfinite(dn)
                cc = cm[ok].reset_index(drop=True)
                Fn, Fi = V.fit(cc, d=dn[ok]), V.fit(cc, d=di[ok])
                Fd = dict(Fn); Fd["th"] = Fi["th"] - Fn["th"]
                A1 = V.amp(Fd)[0]; q2 = float(V.q2vec(Fd) @ Fd["th"])
                if run == "A":
                    r[f"A_{g}_Q2rec"] = q2 / 0.3; r[f"A_{g}_Q1leak"] = A1
                else:
                    r[f"B_{g}_Q1leak"] = A1 / 0.1
        r["pass_all"] = bool(all(0.7 <= r[f"A_{g}_Q2rec"] <= 1.3 and r[f"A_{g}_Q1leak"] < 0.03 and r[f"B_{g}_Q1leak"] < 0.25 for g in GEN))
        rows.append(r)
        log(", ".join(f"{k} {v:.3f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    pd.DataFrame(rows).to_csv(f"{C.TB}/kerr_models.csv", index=False)
    log(f"kerr done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
