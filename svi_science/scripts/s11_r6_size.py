#!/usr/bin/env python
"""s11_r6_size -- S11 round 6, task 2 (PLAN.md §10.2): does the per-episode deep-weighted change ΔD(e) scale with episode size?
Per region: ΔD(e) from the adopted form E4' (L-SEL headline; NOM, E4'+/-Scs2, round-5 E4 columns) regressed on size measures
(GNSS offset / moment / tremor / LFE count / duration / patch coverage / family-set activity), WLS 1/σ², GLS column,
5,000-draw scrambled control, leave-one-episode-out, waveform-level recovery of a size-proportional deep change, and the
top-third-vs-rest tests (during-episode contrast, before->after step).
    S11_REGION=SVI|WA_NW python s11_r6_size.py"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, json, time, argparse
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_r3_fit as LF
import s11_r4_common as R4
import s11_r5_common as R5
import s11_r6_model as M6

RES = f"{C.ROOT}/svi_science/results"
NPERM, SEED = 5000, 0
HEADLINE = {"SVI": "G", "WA_NW": "T"}
UNITS = {"G": "mm (GNSS Tier A)", "M0": "1e18 N·m (Tier B)", "T": "1,000 tremor epicentres", "L": "1,000 LFEs",
         "Du": "10 days", "Pc": "patch share (0-1)", "Fs": "1,000 sumP (family set)", "Fc": "family share (0-1)",
         "G_WA": "mm (WA stations, S1 window)"}


# ------------------------------------------------------------------------------ sizes
def sizes(H, lg):
    U = R5.U
    eps_ = sorted(H.episode.unique(), key=lambda e: U.mid[e])
    S = pd.DataFrame(index=pd.Index(eps_, name="episode"))
    S["major"] = [e in R4.majors() for e in eps_]
    S["mid_PT"] = [str(U.mid[e].date()) for e in eps_]
    b = pd.read_csv(f"{RES}/b_ets_record/episodes_b.csv")
    s1 = pd.read_csv(f"{RES}/s1/episode_geodetic.csv")
    s3 = pd.read_csv(f"{RES}/s3/episodes_all.csv")
    reg_b = "SVI" if C.REGION == "SVI" else "NW"
    bb = b[b.region == reg_b].drop_duplicates("episode").set_index("episode")
    S["T"] = bb.n_trem10.reindex(eps_) / 1000.0
    S["Du"] = bb.duration_d.reindex(eps_) / 10.0
    if C.REGION == "SVI":
        g = s1.drop_duplicates("episode").set_index("episode")
        S["G"] = g.D.reindex(eps_)
        S["M0"] = g.tierB_M0_Nm.reindex(eps_) / 1e18
        sv = s3[s3.region == "SVI"].set_index("episode")
        S["L"] = sv.n_headline.reindex(eps_) / 1000.0
        ev = pd.read_parquet(f"{RES}/s3/events_svi.parquet")
        ev = ev[(ev.series.astype(str) == "data") & (ev.P_same >= 0.75) & ev.keep_xfam8 & ev.counted_once]
        chk = ev.groupby("ep").size().reindex(sv.index).fillna(0) - sv.n_headline
        lg(f"S3 SVI headline events re-count vs n_headline: max |Δ| {float(np.abs(chk).max()):.0f}")
        npatch = 12
        cnt = ev.groupby(["ep", "group"]).size()
    else:
        S["L"] = bb.lfe_nw_once075.reindex(eps_) / 1000.0
        # G_WA: WA-station GNSS offset of the S1 (SVI) episode whose window overlaps this WA episode most
        s1w = s1.assign(st=pd.to_datetime(s1.start), en=pd.to_datetime(s1.end))
        gw = []
        for e in eps_:
            st, en = pd.to_datetime(U.start[e]), pd.to_datetime(U.end[e])
            ov = (np.minimum(s1w.en, en) - np.maximum(s1w.st, st)).dt.days + 1
            j = ov.idxmax() if (ov > 0).any() else None
            gw.append(float(s1w.D_WA[j]) if j is not None and np.isfinite(s1w.D_WA[j]) else np.nan)
        S["G_WA"] = gw
        nv = s3[s3.region == "NW"].set_index("episode")
        ev = pd.read_parquet(f"{RES}/s3/events_nw.parquet")
        ev = ev[(ev.series.astype(str) == "data") & ev.selected_v2 & ev.counted_once]
        chk = ev.groupby("ep").size().reindex(nv.index).fillna(0) - nv.n_headline
        lg(f"S3 NW headline events re-count vs n_headline: max |Δ| {float(np.abs(chk).max()):.0f}")
        ev = ev.assign(group=np.floor(ev.as_km / 10.0).astype(int))
        share = ev.group.value_counts() / len(ev)
        keep = share.index[share >= 0.01]
        ev = ev[ev.group.isin(keep)]; npatch = len(keep)
        cnt = ev.groupby(["ep", "group"]).size()
        sv = nv
    pc = []
    for e in eps_:
        if e not in sv.index:
            pc.append(np.nan); continue
        c = cnt.loc[e] if e in cnt.index.get_level_values(0) else pd.Series(dtype=float)
        pc.append(float((c >= 3).sum()) / npatch)
    S["Pc"] = pc
    S.attrs["npatch"] = npatch
    # family-set activity (the dv/v families themselves; full record)
    if C.REGION == "SVI":
        ad = pd.concat([pd.read_csv(f"{C.MTB}/active_days.csv"), pd.read_csv(f"{C.TB}/ext/active_days.csv")], ignore_index=True)
    else:
        ad = pd.read_csv(f"{C.MTB}/active_days.csv")
    fams = set(H.fam.unique()); ad = ad[ad.fam.isin(fams)]
    pe = ad.groupby(["episode", "fam"]).sumP.sum()
    S["Fs"] = [float(pe.loc[e].sum()) / 1000.0 if e in pe.index.get_level_values(0) else 0.0 for e in eps_]
    S["Fc"] = [float((pe.loc[e] >= 1).sum()) / len(fams) if e in pe.index.get_level_values(0) else 0.0 for e in eps_]
    lg(f"sizes: {len(eps_)} episodes; finite per measure " + ", ".join(f"{k} {int(S[k].notna().sum())}" for k in S.columns
                                                                       if k not in ("major", "mid_PT")) + f"; patches {npatch}; families {len(fams)}")
    return S


# ------------------------------------------------------------------------------ ΔD sources
def sources(adopted):
    kw, proj = M6.FORMS6[adopted]
    out = [("headline", adopted, "two3_20")]
    out.append(("nom", adopted, "ell40"))
    if "Scs2" in adopted:
        out.append(("noSc", adopted.replace("+Scs2", ""), "two3_20"))
    else:
        out.append(("Scs2", f"{adopted}+Scs2", "two3_20"))
    if adopted != "E4":
        out.append(("r5E4", "E4", "two3_20"))
    return out


def form_fac(name):
    if name in M6.FORMS6:
        kw, proj = M6.FORMS6[name]; base = name
    else:                                                  # e.g. E4l+Scs2 / E4h+Scs2
        base = name.replace("+Scs2", ""); kw, proj = M6.FORMS6[base]; kw = {**kw, "Sc": ("spline", 2)}
    return M6.factory6(kw, proj, M6.UTERMS.get(base))


def fit_source(H, inj, name, setup):
    cs, good = R5.setup_cells(H, setup)
    good = good & np.isfinite(inj["nat"])
    cs = cs[good].reset_index(drop=True)
    L = form_fac(name)(cs); F = L.solve()
    eps_ = list(L.eps)
    gv = {e: LF.dvec(F, e) for e in eps_}
    gv = {e: g for e, g in gv.items() if g is not None}
    dD = {e: float(g @ F["th"]) for e, g in gv.items()}
    sig = {e: float(np.sqrt(g @ F["C"] @ g)) for e, g in gv.items()}
    return dict(cs=cs, good=good, L=L, F=F, eps=eps_, g=gv, dD=dD, sig=sig)


# ------------------------------------------------------------------------------ statistics
def wls(x, y, w):
    X = np.column_stack([np.ones_like(x), x]); N = X.T @ (X * w[:, None])
    beta = np.linalg.solve(N, X.T @ (w * y)); r = y - X @ beta
    chi2 = float(np.sum(w * r ** 2) / max(len(y) - 2, 1))
    Cv = np.linalg.inv(N) * max(chi2, 1.0)
    return float(beta[1]), float(np.sqrt(Cv[1, 1])), float(beta[0]), chi2


def wls_slopes_perm(xp, y, w):
    """slopes for many permuted x rows (vectorised WLS)"""
    W = w.sum(); xm = (xp * w).sum(1) / W; ym = (w * y).sum() / W
    dx = xp - xm[:, None]
    return (w * dx * (y - ym)).sum(1) / (w * dx ** 2).sum(1)


def gls(x, y, Sig):
    X = np.column_stack([np.ones_like(x), x])
    Si = np.linalg.pinv(Sig, rcond=1e-10, hermitian=True)
    N = X.T @ Si @ X; beta = np.linalg.solve(N, X.T @ Si @ y); r = y - X @ beta
    chi2 = float(r @ Si @ r / max(len(y) - 2, 1))
    Cv = np.linalg.inv(N) * max(chi2, 1.0)
    return float(beta[1]), float(np.sqrt(Cv[1, 1]))


def regress(src, S, meas, rng, lg, recov=None, majors_only=False):
    ok = [e for e in src["dD"] if e in S.index and np.isfinite(S.loc[e, meas]) and np.isfinite(src["sig"][e]) and src["sig"][e] > 0]
    if majors_only:
        ok = [e for e in ok if bool(S.loc[e, "major"])]
    n = len(ok)
    if n < 5 or np.nanstd(S.loc[ok, meas].to_numpy(float)) == 0:
        return None
    x = S.loc[ok, meas].to_numpy(float); y = np.array([src["dD"][e] for e in ok]); s = np.array([src["sig"][e] for e in ok])
    w = 1 / s ** 2
    b, sb, a, chi2 = wls(x, y, w)
    Gm = np.column_stack([src["g"][e] for e in ok]); Sig = Gm.T @ src["F"]["C"] @ Gm
    bg, sbg = gls(x, y, Sig)
    xp = np.array([rng.permutation(x) for _ in range(NPERM)])
    bp = wls_slopes_perm(xp, y, w)
    p = float((np.sum(np.abs(bp) >= abs(b)) + 1) / (NPERM + 1))
    loo = []
    for i in range(n):
        m = np.arange(n) != i
        if np.std(x[m]) > 0:
            loo.append(wls(x[m], y[m], w[m])[0])
    loo = np.array(loo)
    q10, q90 = np.percentile(x, [10, 90])
    r = dict(measure=meas, unit=UNITS[meas], n=n, n_majors=int(S.loc[ok, "major"].sum()), b=b, sig_b=sb, z=b / sb, a=a, chi2red=chi2,
             p_scr=p, scr_b95=float(np.percentile(np.abs(bp), 95)), scr_b99=float(np.percentile(np.abs(bp), 99)),
             b_gls=bg, sig_b_gls=sbg, b_std=b * float(np.std(x)), sig_b_std=sb * float(np.std(x)),
             loo_max_rel=float(np.max(np.abs(loo - b)) / abs(b)) if abs(b) > 0 else np.nan,
             loo_sign_kept=bool(np.all(np.sign(loo) == np.sign(b))), size_q10=float(q10), size_q90=float(q90),
             size_min=float(x.min()), size_max=float(x.max()), b_det3=3 * sb, dD_det3_q90_q10=3 * sb * (q90 - q10),
             dD_scr99_q90_q10=float(np.percentile(np.abs(bp), 99)) * (q90 - q10), bound2=abs(b) + 2 * sb,
             majors_only=majors_only, episodes=";".join(ok))
    if recov is not None:
        r.update(recovery=recov)
        f = recov if (np.isfinite(recov) and not (0.9 <= recov <= 1.1)) else 1.0
        r.update(bound2_corr=(abs(b) + 2 * sb) / f, b_corr=b / f, dD_det3_corr=3 * sb * (q90 - q10) / f)
    return r


def recovery(src, inj, S, meas, setup):
    """waveform-level: superpose the round-4 per-major deep pulses (0.3 %) scaled by size/size_max"""
    good = src["good"]; dn = inj["nat"][good]; L = src["L"]
    mj = [e for e in R4.majors() if e in src["g"] and e in S.index and np.isfinite(S.loc[e, meas]) and f"Pd|{setup}|{e}" in inj]
    if len(mj) < 4:
        return np.nan, 0
    smax = float(np.max([S.loc[e, meas] for e in mj]))
    if smax <= 0:
        return np.nan, 0
    x = np.zeros(len(dn)); Dtrue = {e: 0.0 for e in src["eps"]}
    for e in mj:
        v = inj[f"Pd|{setup}|{e}"][good]; k = float(S.loc[e, meas]) / smax
        x += k * np.where(np.isfinite(v), v - dn, 0.0); Dtrue[e] = 0.3 * k
    th = L.Op @ x
    ok = [e for e in src["g"] if e in S.index and np.isfinite(S.loc[e, meas])]
    U = R5.U
    yt, yr, xs, ws = [], [], [], []
    for e in ok:
        nb = [z for z in src["eps"] if z != e and abs((U.mid[z] - U.mid[e]).days) <= LF.WIN3]
        yt.append(Dtrue[e] - np.mean([Dtrue[z] for z in nb])); yr.append(float(src["g"][e] @ th))
        xs.append(float(S.loc[e, meas])); ws.append(1 / src["sig"][e] ** 2)
    xs, ws = np.array(xs), np.array(ws)
    bt = wls(xs, np.array(yt), ws)[0]; br = wls(xs, np.array(yr), ws)[0]
    return (br / bt if bt != 0 else np.nan), len(mj)


def contrasts(src, S, meas, rng):
    """top third (by size) vs rest: (a) during-episode ΔD contrast; (b) before->after step at the episode"""
    U = R5.U; F = src["F"]; th, Cm = F["th"], F["C"]
    ok = [e for e in src["g"] if e in S.index and np.isfinite(S.loc[e, meas])]
    n = len(ok)
    if n < 6:
        return []
    k = int(np.ceil(n / 3))
    order = sorted(ok, key=lambda e: -S.loc[e, meas])
    top = set(order[:k])
    out = []
    # step functional: mean D of model episodes in (e, e+3yr] minus mean in [e-3yr, e)
    steps = {}
    for e in ok:
        aft = [z for z in src["eps"] if 0 < (U.mid[z] - U.mid[e]).days <= LF.WIN3]
        bef = [z for z in src["eps"] if -LF.WIN3 <= (U.mid[z] - U.mid[e]).days < 0]
        if not aft or not bef:
            continue
        g = np.zeros(len(th))
        for z in aft:
            g[F["idx"][("D", "", z)]] += 1.0 / len(aft)
        for z in bef:
            g[F["idx"][("D", "", z)]] -= 1.0 / len(bef)
        steps[e] = g
    for kind, gd in (("during", {e: src["g"][e] for e in ok}), ("step", steps)):
        ee = [e for e in ok if e in gd]
        if len(ee) < 6:
            continue
        val = np.array([gd[e] @ th for e in ee]); sg = np.array([np.sqrt(gd[e] @ Cm @ gd[e]) for e in ee]); w = 1 / sg ** 2
        istop = np.array([e in top for e in ee])
        if istop.sum() < 2 or (~istop).sum() < 2:
            continue

        def con(mask):
            return (np.sum(w[mask] * val[mask]) / np.sum(w[mask])) - (np.sum(w[~mask] * val[~mask]) / np.sum(w[~mask]))
        c = con(istop)
        gt = sum(w[i] * gd[e] for i, e in enumerate(ee) if istop[i]) / w[istop].sum()
        gr = sum(w[i] * gd[e] for i, e in enumerate(ee) if not istop[i]) / w[~istop].sum()
        gc = gt - gr; sc = float(np.sqrt(gc @ Cm @ gc))
        cp = np.array([con(rng.permutation(istop)) for _ in range(NPERM)])
        p = float((np.sum(np.abs(cp) >= abs(c)) + 1) / (NPERM + 1))
        out.append(dict(measure=meas, kind=kind, n=len(ee), n_top=int(istop.sum()), top=";".join(e for e in ee if e in top),
                        contrast=float(c), sig=sc, z=float(c / sc), p_scr=p, scr_c95=float(np.percentile(np.abs(cp), 95)),
                        top_mean=float(np.sum(w[istop] * val[istop]) / w[istop].sum()),
                        rest_mean=float(np.sum(w[~istop] * val[~istop]) / w[~istop].sum())))
    return out


def verdict(R, Cn):
    h = HEADLINE[C.REGION]
    r = R[(R.source == "headline") & (R.measure == h) & ~R.majors_only]
    if not len(r):
        return dict(verdict="NO DATA", measure=h)
    r = r.iloc[0]
    sc = R[(R.source.isin(["Scs2", "noSc"])) & (R.measure == h) & ~R.majors_only]
    sc_ok = bool(len(sc) and np.sign(sc.iloc[0].b) == np.sign(r.b) and abs(sc.iloc[0].z) >= 2)
    if abs(r.z) >= 3 and r.p_scr < 0.01 and sc_ok and r.loo_sign_kept:
        v = "SCALES"
    elif abs(r.z) >= 2 and r.p_scr < 0.05:
        v = "HINT"
    else:
        v = "NO SCALING SEEN"
    return dict(verdict=v, measure=h, unit=UNITS[h], b=r.b, sig_b=r.sig_b, z=r.z, p_scr=r.p_scr, recovery=r.get("recovery"),
                bound2=r.bound2, bound2_corr=r.get("bound2_corr"), b_det3=r.b_det3, dD_det3_q90_q10=r.dD_det3_q90_q10,
                dD_det3_corr=r.get("dD_det3_corr"), n=int(r.n), sc_column_same_sign_2sig=sc_ok, loo_sign_kept=bool(r.loo_sign_kept))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--test", action="store_true"); ap.add_argument("--amend", action="store_true")
    ap.add_argument("--form", default=None, help="post-hoc descriptive column (not pre-registered): ΔD source form, L-SEL + NOM")
    a = ap.parse_args()
    TB = M6.TB6 + ("/test" if a.test else ""); os.makedirs(TB, exist_ok=True)
    lg = M6.log("size_test" if a.test else ("size_posthoc" if a.form else ("size_amend" if a.amend else "size")))
    pre = "r6p" if a.form else ("r6a" if a.amend else "r6")
    t0 = time.time()
    adopted = json.load(open(f"{M6.R6}/adopted.json"))["task2_source"] if not a.test else "E4y"
    H, inj = R5.load_H(lg)
    S = sizes(H, lg)
    if a.test:                                             # code test: sizes shuffled first, so no real answer is seen
        rt = np.random.default_rng(123)
        for k in [c for c in S.columns if c not in ("major", "mid_PT")]:
            S[k] = rt.permutation(S[k].to_numpy())
    S.to_csv(f"{TB}/{pre}_sizes.csv")
    rows, crows, dd = [], [], []
    srcs = sources(adopted)
    if a.amend:                                            # §10 Amendment 1: the amendment form as a labelled column
        aa = a.form or json.load(open(f"{M6.R6}/adopted_amend.json"))["adopted"]
        srcs = [("amend", aa, "two3_20"), ("amend_nom", aa, "ell40")] if aa else []
        lg(f"amendment form: {aa}")
    for tag, name, setup in srcs:
        src = fit_source(H, inj, name, setup)
        for e in src["dD"]:
            dd.append(dict(source=tag, form=name, setup=setup, episode=e, dD=src["dD"][e], sig=src["sig"][e]))
        meas = [m for m in ("G", "M0", "G_WA", "T", "L", "Du", "Pc", "Fs", "Fc") if m in S.columns]
        for m in meas:
            rng = np.random.default_rng(SEED)
            rec, nrec = recovery(src, inj, S, m, setup)
            for mo in (False, True):
                r = regress(src, S, m, rng, lg, recov=rec, majors_only=mo)
                if r is None:
                    continue
                r.update(source=tag, form=name, setup=setup, n_recov_majors=nrec); rows.append(r)
            if tag in ("headline", "Scs2", "noSc", "amend"):
                for c in contrasts(src, S, m, np.random.default_rng(SEED)):
                    c.update(source=tag, form=name, setup=setup); crows.append(c)
        lg(f"source {tag} ({name}|{setup}) done")
    R = pd.DataFrame(rows); Cn = pd.DataFrame(crows)
    R.to_csv(f"{TB}/{pre}_size_regress.csv", index=False); Cn.to_csv(f"{TB}/{pre}_size_topthird.csv", index=False)
    pd.DataFrame(dd).to_csv(f"{TB}/{pre}_size_dD.csv", index=False)
    if a.amend:
        lg(f"done {time.time() - t0:.0f}s")
        for r in R[~R.majors_only].itertuples() if len(R) else []:
            lg(f"REG {r.source:9s} {r.measure:4s} n {r.n:2d}: b {r.b:+.4f} ± {r.sig_b:.4f} % per {r.unit} (z {r.z:+.2f}, scrambled p {r.p_scr:.3f}); "
               f"recovery {r.recovery:.2f}; LOO sign kept {r.loo_sign_kept}")
        for c in Cn.itertuples() if len(Cn) else []:
            lg(f"TOP3 {c.source:9s} {c.measure:4s} {c.kind:6s}: {c.contrast:+.4f} ± {c.sig:.4f} (z {c.z:+.2f}, scrambled p {c.p_scr:.3f})")
        open(f"{M6.R6}/SIZE_{'POSTHOC' if a.form else 'AMEND'}_{C.REGION}_DONE", "w").write(R5.pt_now() + "\n")
        return
    for r in R[~R.majors_only].itertuples():
        lg(f"REG {r.source:8s} {r.measure:4s} n {r.n:2d}: b {r.b:+.4f} ± {r.sig_b:.4f} % per {r.unit} (z {r.z:+.2f}, scrambled p {r.p_scr:.3f}, "
           f"99% |b_scr| {r.scr_b99:.4f}); GLS {r.b_gls:+.4f} ± {r.sig_b_gls:.4f}; per SD {r.b_std:+.4f} ± {r.sig_b_std:.4f} %; "
           f"recovery {r.recovery:.2f}; 3σ change q90-q10 {r.dD_det3_q90_q10:.3f} %; LOO sign kept {r.loo_sign_kept}")
    for c in Cn.itertuples():
        lg(f"TOP3 {c.source:8s} {c.measure:4s} {c.kind:6s}: top {c.top_mean:+.4f} vs rest {c.rest_mean:+.4f} -> {c.contrast:+.4f} ± {c.sig:.4f} "
           f"(z {c.z:+.2f}, scrambled p {c.p_scr:.3f}; n {c.n}, top {c.n_top})")
    v = verdict(R, Cn); v.update(region=C.REGION, adopted=adopted, written=R5.pt_now())
    R5.dump(v, f"{TB}/r6_size_verdict.json")
    lg(f"VERDICT {C.REGION}: {v}")
    lg(f"done {time.time() - t0:.0f}s")
    if not a.test:
        open(f"{M6.R6}/SIZE_{C.REGION}_DONE", "w").write(R5.pt_now() + "\n")


if __name__ == "__main__":
    main()
