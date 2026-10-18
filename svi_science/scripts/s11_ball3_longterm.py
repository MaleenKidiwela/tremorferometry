#!/usr/bin/env python
"""s11_ball3_longterm -- S11 ballistic round 3, §1 pooled long-term constraints per station pair (round3/PLAN.md).
Items (YR headline, Y2 column for dS): dS and dP3 of every station X against each reference set (WA: G1BH, G2, LRIV,
G1all; SVI: BH, B011, B0910, PGC), clock-corrected (uncorrected column for the headline references); co-located dS inside
each group; S - P per station (P3 headline; P3s, S3c-P3 columns); scrambled control. Statistics as round 2
(s11_ball2_longterm.contrast_stats / series_frame, s11_ball_analysis.decomp). Per-pair table with 2-sigma bounds and the
detection / replication rule.

    python s11_ball3_longterm.py --region WA_NW|SVI [--tag _test] [--workers 8]
"""
from __future__ import annotations
import os, sys, json, argparse, zlib
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_ball3_lib as L
import s11_ball2_longterm as LT
import s11_ball_analysis as A1

SEED = 20261013


def run_item(args):
    name, d, conv, kind, bk, extra, nperm = args
    rng = np.random.default_rng(SEED + (zlib.crc32(name.encode()) % 100000))
    A1.NPERM = nperm
    r = dict(item=name, bin=bk, kind=kind, n_values=len(d), n_series=d.sid.nunique() if len(d) else 0)
    if len(d):
        o, _ = A1.decomp(d.assign(date=d.date), "sig", name, rng)
        ve = o.get("var_ep", np.nan)
        r.update(sd_year=np.sqrt(max(ve, 0)) if np.isfinite(ve) else np.nan, year_perm_p=o.get("perm_p", np.nan),
                 scatter_ms=o.get("scatter_ms", np.nan), scatter_rob_ms=o.get("scatter_robust_ms", np.nan))
        r.update(LT.contrast_stats(d, rng, nperm))
        r["cc_med"] = float(np.nanmedian(d.cc)) if "cc" in d else np.nan
    if extra:
        r.update(extra)
    for k_, T in conv.items():
        if T is not None and np.isfinite(T) and abs(T) > 0 and "decade_ms" in r and np.isfinite(r.get("decade_ms", np.nan)):
            sgn = 1 if k_.startswith("+") else -1
            r[f"dec_{k_}_pct"] = sgn * r["decade_ms"] / 1e3 / T * 100
            r[f"dec_{k_}_sig_pct"] = r["decade_sig_ms"] / 1e3 / abs(T) * 100
            r[f"dec_{k_}_bound2s_pct"] = (abs(r["decade_ms"]) + 2 * r["decade_sig_ms"]) / 1e3 / abs(T) * 100
            sgc = max(r["decade_sig_ms"], r.get("decade_perm_sd_ms", 0.0))
            r[f"dec_{k_}_bound2s_cons_pct"] = (abs(r["decade_ms"]) + 2 * sgc) / 1e3 / abs(T) * 100
    cy = None
    if kind == "data" and len(d):
        cy = d.groupby("lab").apply(lambda x: pd.Series(dict(common=np.average(x.r, weights=1 / x.sig ** 2), n=len(x),
                                                             cc=np.nanmedian(x.cc) if "cc" in x else np.nan))).reset_index()
        cy["item"], cy["bin"] = name, bk
    return r, cy


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--region", default="WA_NW"); ap.add_argument("--tag", default="")
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args(); reg, TB, _, _ = L.dirs(a.region); tag = a.tag
    st = L.load(reg, "lags_station", tag); pr = L.load(reg, "lags_pair", tag)
    Dc = pd.read_parquet(f"{TB}/pair_corrected{tag}.parquet"); DP = pd.read_parquet(f"{TB}/pairP_corrected{tag}.parquet")
    SP = L.sp_table(reg, tag).to_dict()
    refs = L.REFSETS[reg]; head = L.HEAD_REFS[reg]
    stations = sorted(set(st.netsta))
    tasks = []
    # ---- dS / dP3 per reference x station
    for obs, D_, nb in (("dS", Dc, ("YR", "Y2")), ("dP3", DP, ("YR",))):
        if obs == "dP3":
            D_ = D_.assign(ccS_A=D_.ccP3_A, ccS_B=D_.ccP3_B)
        Lg = L.oriented(D_[D_.bin.isin(nb)])
        for rname, rs in refs.items():
            for col, corr in (("rc", "corrected"), ("r", "uncorrected")):
                if corr == "uncorrected" and rname not in head:
                    continue
                F = L.ref_frame(Lg, rs, col)
                for bk in nb:
                    if bk == "Y2" and (corr == "uncorrected" or rname not in head):
                        continue
                    Fb = F[F.bin == bk]
                    for X, gX in Fb.groupby("X"):
                        d = LT.series_frame(gX, "v", "sig", ["fam", "X", "epX", "bin"], "ccX")
                        if d.sid.nunique() < 3:
                            continue
                        fams = d.fam.unique()
                        tsx = np.nanmean([L.TS(SP.get((f, X), np.nan)) for f in fams])
                        tsr = np.nanmean([np.nanmean([L.TS(SP.get((f, m), np.nan)) for m in rs]) for f in fams])
                        if obs == "dS":
                            conv = {"dVs_Xleg": tsx, "dVs_uniform": (tsx - tsr) if abs(tsx - tsr) >= 2 else np.nan}
                        else:
                            conv = {"dVp_Xleg": tsx / L.VPVS, "dVp_uniform": (tsx - tsr) / L.VPVS if abs(tsx - tsr) >= 2 else np.nan}
                        tasks.append((f"{obs} {rname}-{X} {corr}", d, conv, "data", bk,
                                      dict(obs=obs, ref=rname, X=X, corr=corr, TS_X=tsx, TS_ref=tsr), 10000 if rname in head else 2000))
    # ---- scrambled control for dS (uncorrected; same reference averaging), headline refs, YR
    pc = pr[(pr.kind == "ctrl") & (pr.bin == "YR")].copy(); pc.loc[pc.railS, "dS"] = np.nan
    pc["sig"] = np.maximum(pc.sig_dS, 0.5)
    pc = pc[np.isfinite(pc.dS)]
    pc["K"] = pc.groupby(["fam", "pair", "epA", "epB"]).dS.transform("size"); pc = pc[pc.K >= 3]
    pc["r"] = pc.dS - pc.groupby(["fam", "pair", "epA", "epB"]).dS.transform("mean")
    Lc = L.oriented(pc, cols=("r",))
    for rname in head:
        F = L.ref_frame(Lc, refs[rname], "r")
        for X, gX in F.groupby("X"):
            d = LT.series_frame(gX, "v", "sig", ["fam", "X", "epX", "bin"], "ccX")
            if d.sid.nunique() >= 3:
                tasks.append((f"dS {rname}-{X} corrected", d, {}, "ctrl", "YR", dict(obs="dS", ref=rname, X=X), 2000))
    # ---- co-located dS inside groups (corrected), YR
    ng = L.node_group(reg)
    for gname, mem in L.GROUPS[reg].items():
        W = Dc[(Dc.bin == "YR") & Dc.A.isin(mem) & Dc.B.isin(mem)]
        d = LT.series_frame(W.assign(v=W.rc, ccX=np.minimum(W.ccS_A, W.ccS_B)), "v", "sig", L.PKEY[:4] + ["bin"], "ccX")
        if len(d):
            tasks.append((f"coloc dS {gname} corrected", d, {}, "data", "YR", dict(obs="coloc", ref=gname), 10000))
        for pair, Wp in W.groupby("pair"):
            d = LT.series_frame(Wp.assign(v=Wp.rc, ccX=np.minimum(Wp.ccS_A, Wp.ccS_B)), "v", "sig", L.PKEY[:4] + ["bin"], "ccX")
            if d.sid.nunique() >= 3:
                tasks.append((f"coloc dS {pair} corrected", d, {}, "data", "YR", dict(obs="coloc", ref=gname, X=pair), 2000))
    # ---- S - P per station (YR): P3 headline, P3s and S3c-P3 columns; control
    for kind in ("data", "ctrl"):
        s = st[(st.kind == kind) & (st.bin == "YR")].copy()
        for est, (sk, pk, sg) in {"S-P3": ("S", "P3", "sigSmP3"), "S-P3s": ("S", "P3s", "sigSmP3s"),
                                  "S3c-P3": ("S3c", "P3", "sigS3cmP3")}.items():
            if kind == "ctrl" and est != "S-P3":
                continue
            s["SmP"] = np.where(s[f"rail{sk}"] | s[f"rail{pk}"], np.nan, s[f"lag{sk}"] - s[f"lag{pk}"])
            for ns in stations:
                d = LT.series_frame(s[s.netsta == ns], "SmP", sg, ["fam", "netsta", "epoch", "bin"], f"cc{pk}")
                if d.sid.nunique() < 3:
                    continue
                sp = float(np.nanmean([SP.get((f, ns), np.nan) for f in d.fam.unique()]))
                tasks.append((f"{est} {ns}", d, {"constVpVs": sp, "VsOnly": L.TS(sp), "+VpOnly": L.TP(sp)}, kind, "YR",
                              dict(obs=est, X=ns, SP_s=sp), 10000 if (est == "S-P3" and kind == "data") else 2000))
    print(f"{len(tasks)} items", flush=True)
    import multiprocessing as mp
    with mp.get_context("fork").Pool(a.workers) as p:
        out = p.map(run_item, tasks, chunksize=4)
    R = pd.DataFrame([o[0] for o in out]); Cv = pd.concat([o[1] for o in out if o[1] is not None], ignore_index=True)
    R.to_csv(f"{TB}/longterm{tag}.csv", index=False); Cv.to_csv(f"{TB}/longterm_curves{tag}.csv", index=False)
    calls(R, reg, TB, tag)


def calls(R, reg, TB, tag):
    head = L.HEAD_REFS[reg]
    D = R[(R.kind == "data") & (R.bin == "YR")]
    C = R[(R.kind == "ctrl") & (R.bin == "YR")].set_index("item")
    rows = []
    for _, r in D[(D.obs == "dS") & (D["corr"] == "corrected")].iterrows():
        cp = float(C.decade_perm_p.get(r["item"], np.nan)) if "decade_perm_p" in C else np.nan
        det = bool(np.isfinite(r.get("decade_ms", np.nan)) and abs(r.decade_ms) >= 3 * r.decade_sig_ms and r.decade_perm_p < 0.01
                   and not (np.isfinite(cp) and cp < 0.01))
        rows.append(dict(ref=r.ref, X=r.X, n_series=r.get("decade_n_series"), decade_ms=r.get("decade_ms"),
                         decade_sig_ms=r.get("decade_sig_ms"), perm_p=r.get("decade_perm_p"), perm_sd_ms=r.get("decade_perm_sd_ms"),
                         ctrl_perm_p=cp, chi2dof=r.get("decade_chi2dof"), TS_X=r.TS_X, TS_ref=r.TS_ref, cc_med=r.get("cc_med"),
                         leg_bound2s_pct=r.get("dec_dVs_Xleg_bound2s_pct"), leg_bound2s_cons_pct=r.get("dec_dVs_Xleg_bound2s_cons_pct"),
                         leg_dVs_pct=r.get("dec_dVs_Xleg_pct"), leg_sig_pct=r.get("dec_dVs_Xleg_sig_pct"),
                         uni_bound2s_pct=r.get("dec_dVs_uniform_bound2s_pct"), trend_ms_per_yr=r.get("trend_ms_per_yr"),
                         trend_sig=r.get("trend_sig"), year_perm_p=r.get("year_perm_p"), sd_year=r.get("sd_year"), det_single=det))
    P = pd.DataFrame(rows)
    # detection requires the same sign, p < 0.05 with all headline references
    if len(P):
        P["detection"] = False
        for X, g in P.groupby("X"):
            h = g[g.ref.isin(head)]
            for i, r in h.iterrows():
                oth = h[h.ref != r.ref]
                ok = r.det_single and len(oth) >= 1 and all((np.sign(o.decade_ms) == np.sign(r.decade_ms)) and o.perm_p < 0.05
                                                           for _, o in oth.iterrows())
                P.loc[i, "detection"] = bool(ok)
    P.to_csv(f"{TB}/pairs_decade{tag}.csv", index=False)
    if reg == "WA_NW" and len(P):
        P.to_csv(f"{L.ROOT3}/tables/wa_pairs_decade{tag}.csv", index=False)
    # rule 1 verdict
    v = {}
    co = D[D.obs == "coloc"]
    v["coloc"] = {r["item"]: dict(decade_ms=r.get("decade_ms"), sig=r.get("decade_sig_ms"), p=r.get("decade_perm_p"),
                                  within3s=bool(abs(r.get("decade_ms", 0)) < 3 * r.get("decade_sig_ms", np.inf)))
                  for _, r in co[co.X.isna()].iterrows()} if "X" in co else {}
    v["rule_a"] = bool(all(x["within3s"] for x in v["coloc"].values())) if v["coloc"] else None
    if len(P):
        H = P[P.ref.isin(head)]
        piv = H.pivot_table(index="X", columns="ref", values=["leg_bound2s_pct", "n_series"], aggfunc="first")
        good = []
        for X in piv.index:
            try:
                b = [piv.loc[X, ("leg_bound2s_pct", r)] for r in head]; n = [piv.loc[X, ("n_series", r)] for r in head]
            except KeyError:
                continue
            if all(np.isfinite(b)) and max(b) <= 0.03 and min(n) >= 20:
                good.append(dict(X=X, bounds=[float(x) for x in b], n=[int(x) for x in n]))
        v["rule_b_legs"] = good
        v["detections"] = P[P.detection][["ref", "X", "decade_ms", "decade_sig_ms", "perm_p"]].to_dict("records")
        v["rule_b"] = bool(len(good) > 0 or len(v["detections"]) > 0)
        v["PASS"] = bool(v["rule_a"] and v["rule_b"]) if v["rule_a"] is not None else v["rule_b"]
    json.dump(v, open(f"{TB}/rule1{tag}.json", "w"), indent=1, default=float)
    with pd.option_context("display.width", 250, "display.max_rows", 400):
        cols = ["ref", "X", "n_series", "decade_ms", "decade_sig_ms", "perm_p", "ctrl_perm_p", "TS_X", "leg_bound2s_pct",
                "leg_bound2s_cons_pct", "uni_bound2s_pct", "detection"]
        if len(P):
            print(P[P.ref.isin(head)].sort_values(["X", "ref"])[cols].round(4).to_string())
        print(json.dumps(v, indent=1, default=float)[:3000])


if __name__ == "__main__":
    main()
