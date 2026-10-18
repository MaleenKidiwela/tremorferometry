#!/usr/bin/env python
"""s11_r7_hint -- S11 round 7, task 4 (PLAN.md §11.4): independent test of the round-6 WA top-third hint (longest /
largest episodes read ~0.07 % higher during the episode). Round-5 E4 L-SEL ΔD; top set fixed on the full WA fit; WA
family halves and station halves; SVI replication; stack-length column.
    S11_REGION=... python s11_r7_hint.py            (per region)
    python s11_r7_hint.py --stage verdict           (both regions -> round7/hint_verdict.json)
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_r7_common as Q
import s11_r5_common as R5
import s11_r4_common as R4

NPERM = 5000
MEAS = ["Du", "T"]


def fit_cells(cs):
    L = Q.e4_model(cs, [("D", "", cs.a_d.to_numpy(float))]); F = L.solve(cs.d.to_numpy(float))
    g = {}
    for e in L.eps:
        v = Q.dvec_g(F, L.eps, e, "D")
        if v is not None:
            g[e] = v
    return L, F, g


def contrast(F, g, top, eps_ok, rng, label, meas):
    ee = [e for e in eps_ok if e in g]
    if len(ee) < 6:
        return None
    th, Cm = F["th"], F["C"]
    val = np.array([g[e] @ th for e in ee]); sg = np.array([np.sqrt(g[e] @ Cm @ g[e]) for e in ee]); w = 1 / sg ** 2
    ist = np.array([e in top for e in ee])
    if ist.sum() < 2 or (~ist).sum() < 2:
        return None
    con = lambda m: np.sum(w[m] * val[m]) / np.sum(w[m]) - np.sum(w[~m] * val[~m]) / np.sum(w[~m])
    c = float(con(ist))
    gt = sum(w[i] * g[e] for i, e in enumerate(ee) if ist[i]) / w[ist].sum()
    gr = sum(w[i] * g[e] for i, e in enumerate(ee) if not ist[i]) / w[~ist].sum()
    s = float(np.sqrt((gt - gr) @ Cm @ (gt - gr)))
    cp = np.array([con(rng.permutation(ist)) for _ in range(NPERM)])
    p = float((np.sum(np.abs(cp) >= abs(c)) + 1) / (NPERM + 1))
    return dict(label=label, measure=meas, n=len(ee), n_top=int(ist.sum()), top=";".join(e for e in ee if e in top), contrast=c, sig=s,
                z=c / s, p_scr=p, scr_c95=float(np.percentile(np.abs(cp), 95)),
                top_mean=float(np.sum(w[ist] * val[ist]) / w[ist].sum()), rest_mean=float(np.sum(w[~ist] * val[~ist]) / w[~ist].sum()),
                g_diff=gt - gr)


def region_stage(lg, test=False):
    TB = Q.tbdir(test)
    H, inj = R5.load_H(lg)
    cs0, good = R5.setup_cells(H, Q.LSEL); good = good & np.isfinite(inj["nat"])
    cs = cs0[good].reset_index(drop=True)
    S = pd.read_csv(f"{Q.BASE}/round6/{Q.SUB}/tables/r6_sizes.csv").set_index("episode")
    if test:
        rt = np.random.default_rng(123)
        for k in MEAS:
            S[k] = rt.permutation(S[k].to_numpy())
    L, F, g = fit_cells(cs)
    # check against the round-6 headline ΔD
    r6 = pd.read_csv(f"{Q.BASE}/round6/{Q.SUB}/tables/r6_size_dD.csv")
    r6 = r6[(r6.source == "headline") & (r6.form == "E4") & (r6.setup == Q.LSEL)].set_index("episode")
    dd = [abs(float(g[e] @ F["th"]) - r6.dD[e]) for e in g if e in r6.index]
    lg(f"full fit reproduces round-6 headline ΔD: max |Δ| {max(dd) if dd else np.nan:.1e} over {len(dd)} episodes")
    tops = {}
    for m in MEAS:
        ok = [e for e in g if e in S.index and np.isfinite(S.loc[e, m])]
        k = int(np.ceil(len(ok) / 3)); tops[m] = (set(sorted(ok, key=lambda e: -S.loc[e, m])[:k]), ok)
        lg(f"{m}: {len(ok)} episodes with a size, top {k}: {sorted(tops[m][0])}")
    # splits
    w = 1 / cs.sig ** 2
    sw = w.groupby(cs.site).sum().sort_values(ascending=False)
    half_a = set(sw.index[0::2]); half_b = set(sw.index[1::2])
    splits = [("full", np.ones(len(cs), bool)), ("fam_odd", (cs.fam % 2 == 1).to_numpy()), ("fam_even", (cs.fam % 2 == 0).to_numpy()),
              ("sta_A", cs.site.isin(half_a).to_numpy()), ("sta_B", cs.site.isin(half_b).to_numpy())]
    lg(f"station halves: A {sorted(half_a)} (weight {float(sw[list(half_a)].sum() / sw.sum()):.2f}); B {sorted(half_b)}")
    rows, gd = [], {}
    for lab, m in splits:
        sub = cs[m].reset_index(drop=True)
        sub, nd = R4.identifiable(sub)
        Ls, Fs, gs = (L, F, g) if lab == "full" else fit_cells(sub)
        for meas in MEAS:
            top, ok = tops[meas]
            r = contrast(Fs, gs, top, ok, np.random.default_rng(Q.SEED), lab, meas)
            if r is None:
                lg(f"[{lab}|{meas}] too few episodes"); continue
            gd[(lab, meas)] = r.pop("g_diff")
            r.update(n_cells=len(sub), n_dropped_2site=nd, chi2_fit=Fs["chi2"]); rows.append(r)
            lg(f"[{lab}|{meas}] contrast {r['contrast']:+.4f} ± {r['sig']:.4f} % (z {r['z']:+.2f}, scrambled p {r['p_scr']:.3f}; "
               f"n {r['n']}, top {r['n_top']}; cells {len(sub)})")
    T = pd.DataFrame(rows)
    # differences between halves (independent cell sets -> σ_diff = hypot)
    for a_, b_ in (("fam_odd", "fam_even"), ("sta_A", "sta_B")):
        for meas in MEAS:
            x = T[(T.label == a_) & (T.measure == meas)]; y = T[(T.label == b_) & (T.measure == meas)]
            if len(x) and len(y):
                d = float(x.iloc[0].contrast - y.iloc[0].contrast); s = float(np.hypot(x.iloc[0].sig, y.iloc[0].sig))
                rows.append(dict(label=f"{a_}-{b_}", measure=meas, contrast=d, sig=s, z=d / s))
    T = pd.DataFrame(rows)
    T.to_csv(f"{TB}/r7_hint_contrasts.csv", index=False)
    # stack-length column
    nb = (w * cs.n).groupby(cs.episode).sum() / w.groupby(cs.episode).sum()
    srow = []
    for meas in MEAS:
        ok = [e for e in g if e in S.index and np.isfinite(S.loc[e, meas]) and e in nb.index]
        if len(ok) < 6:
            continue
        x1 = S.loc[ok, meas].to_numpy(float); x2 = np.log(nb[ok].to_numpy(float))
        y = np.array([float(g[e] @ F["th"]) for e in ok]); s = np.array([float(np.sqrt(g[e] @ F["C"] @ g[e])) for e in ok]); ww = 1 / s ** 2
        out = dict(measure=meas, n=len(ok), corr_size_logn=float(np.corrcoef(x1, x2)[0, 1]))
        for nm, X in (("size_only", np.column_stack([np.ones(len(ok)), x1])), ("size_and_logn", np.column_stack([np.ones(len(ok)), x1, x2]))):
            N = X.T @ (X * ww[:, None]); be = np.linalg.solve(N, X.T @ (ww * y)); r = y - X @ be
            chi2 = float(np.sum(ww * r ** 2) / max(len(y) - X.shape[1], 1)); Cv = np.linalg.inv(N) * max(chi2, 1.0)
            out[f"{nm}_b"] = float(be[1]); out[f"{nm}_sig"] = float(np.sqrt(Cv[1, 1])); out[f"{nm}_z"] = float(be[1] / np.sqrt(Cv[1, 1]))
            if X.shape[1] == 3:
                out["logn_b"] = float(be[2]); out["logn_sig"] = float(np.sqrt(Cv[2, 2]))
        srow.append(out)
        lg(f"stack length [{meas}]: corr(size, log n̄) {out['corr_size_logn']:+.2f}; size slope {out['size_only_b']:+.4f} ± {out['size_only_sig']:.4f} "
           f"-> with log n̄ {out['size_and_logn_b']:+.4f} ± {out['size_and_logn_sig']:.4f} (log n̄ {out['logn_b']:+.4f} ± {out['logn_sig']:.4f})")
    pd.DataFrame(srow).to_csv(f"{TB}/r7_hint_stack.csv", index=False)
    pd.DataFrame(dict(episode=list(g), dD=[float(g[e] @ F["th"]) for e in g], sig=[float(np.sqrt(g[e] @ F["C"] @ g[e])) for e in g],
                      nbar=[float(nb.get(e, np.nan)) for e in g], Du=[S.Du.get(e, np.nan) for e in g], T=[S["T"].get(e, np.nan) for e in g],
                      top_Du=[e in tops["Du"][0] for e in g], top_T=[e in tops["T"][0] for e in g])).to_csv(f"{TB}/r7_hint_dD.csv", index=False)


def verdict_stage(lg):
    out = dict(written=Q.pt_now())
    T = {}
    for sub in ("svi", "wa_nw"):
        p = f"{Q.BASE}/round7/{sub}/tables/r7_hint_contrasts.csv"
        if os.path.exists(p):
            T[sub] = pd.read_csv(p)
    st = {sub: pd.read_csv(f"{Q.BASE}/round7/{sub}/tables/r7_hint_stack.csv") for sub in T if os.path.exists(f"{Q.BASE}/round7/{sub}/tables/r7_hint_stack.csv")}
    for meas in MEAS:
        W = T.get("wa_nw"); Sv = T.get("svi")
        g = lambda D, lab: D[(D.label == lab) & (D.measure == meas)].iloc[0] if D is not None and len(D[(D.label == lab) & (D.measure == meas)]) else None
        cons = []
        for a_, b_ in (("fam_odd", "fam_even"), ("sta_A", "sta_B")):
            x, y, d = g(W, a_), g(W, b_), g(W, f"{a_}-{b_}")
            cons.append(bool(x is not None and y is not None and x.contrast > 0 and y.contrast > 0 and abs(d.z) < 2))
        wa_cons = all(cons)
        sv = g(Sv, "full")
        svi_rep = bool(sv is not None and sv.contrast > 0 and sv.z >= 2 and sv.p_scr < 0.05)
        v = "REPLICATED" if (wa_cons and svi_rep) else ("WA-ONLY" if wa_cons else "NOT ROBUST")
        sl = None
        if "wa_nw" in st:
            s_ = st["wa_nw"]; s_ = s_[s_.measure == meas]
            if len(s_):
                sl = bool(abs(s_.iloc[0].size_and_logn_z) < 1 and abs(s_.iloc[0].size_only_z) >= 1)
        wf = g(W, "full")
        out[meas] = dict(verdict=v, wa_consistent=wa_cons, fam_split_ok=cons[0], sta_split_ok=cons[1], svi_replicates=svi_rep,
                         wa_full=[float(wf.contrast), float(wf.sig), float(wf.p_scr)] if wf is not None else None,
                         svi_full=[float(sv.contrast), float(sv.sig), float(sv.p_scr)] if sv is not None else None,
                         stack_length_flag=sl)
        lg(f"[{meas}] {v}: " + json.dumps(out[meas], default=str))
    Q.dump(out, f"{Q.R7}/hint_verdict.json")
    open(f"{Q.R7}/HINT_VERDICT_DONE", "w").write(Q.pt_now() + "\n")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--stage", default="region"); ap.add_argument("--test", action="store_true")
    a = ap.parse_args()
    if a.stage == "verdict":
        verdict_stage(C.Log(f"{Q.R7}/hint_verdict.log")); return
    lg = Q.log("hint_test" if a.test else "hint"); t0 = time.time()
    region_stage(lg, a.test)
    lg(f"done {time.time() - t0:.0f}s")
    if not a.test:
        open(f"{Q.R7}/HINT_{C.REGION}_DONE", "w").write(Q.pt_now() + "\n")


if __name__ == "__main__":
    main()
