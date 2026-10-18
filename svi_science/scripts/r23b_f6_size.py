#!/usr/bin/env python
"""r23b_f6_size (2026-10-11) -- Finding 6b: LFE size per detection (S6 dM, fixed CN station set) 2003-2026 within ETS.
Rules: svi_science/results/r23b/PLAN.md §3b (pre-registered). Inputs: S6 forward amplitudes (mf_station.parquet,
targets.parquet, window_check.csv) + r23b backward amplitudes (tables/f6/mf_cn_backward.parquet).

    python svi_science/scripts/r23b_f6_size.py
"""
from __future__ import annotations
import importlib.util, json, time
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
TB = f"{ROOT}/svi_science/results/r23b/tables"
OUT = f"{TB}/f6"
S6 = f"{ROOT}/svi_science/results/s6"
CN = ["CN.NLLB", "CN.PGC", "CN.SNB", "CN.VGZ"]
CHANGES = {"CN.PGC": ["2015-01-12", "2017-08-10"], "CN.SNB": ["2006-04-20", "2008-12-11", "2011-03-17", "2017-06-08"],
           "CN.VGZ": ["2007-02-28", "2009-07-31", "2017-03-03"], "CN.NLLB": ["2007-01-27", "2007-02-27", "2017-02-23"]}
NB = 2000; NP = 2000
rng = np.random.default_rng(66)
LOGF = open(f"{ROOT}/svi_science/results/r23b/logs/r23b_f6_size.log", "w")
spec = importlib.util.spec_from_file_location("s3b", f"{ROOT}/svi_science/scripts/s3_build.py")
s3b = importlib.util.module_from_spec(spec); spec.loader.exec_module(s3b)


def log(m):
    s = f"[{time.strftime('%H:%M:%S')}] {m}"; print(s, flush=True); LOGF.write(s + "\n"); LOGF.flush()


def epoch_of(sta, t):
    return int(np.sum([t >= pd.Timestamp(c) for c in CHANGES[sta]]))


def load():
    W = pd.read_csv(f"{TB}/windows.csv", parse_dates=["w0", "w1"])
    M = W[W.major]
    F = pd.read_csv(f"{ROOT}/svi_science/results/s0/families.csv")
    # forward
    T = pd.read_parquet(f"{S6}/targets.parquet")
    T = T[T.kind == "event"]
    T["head"] = (T.P_same >= 0.75) & T.keep_xfam8
    Sm = pd.read_parquet(f"{S6}/mf_station.parquet")
    Sm = Sm[Sm.sta.isin(CN)]
    Rf = Sm.merge(T[["det_id", "fam", "time_utc", "head"]], on="det_id").rename(columns={"det_id": "ev_id"})
    Rf["ep"] = s3b.assign(Rf.time_utc, M)
    Rf["seg"] = "FWD"; Rf["src"] = "fw"
    # FWD-S15 = forward family-days on the S15 station set in majors starting before 2019
    fr = pd.read_parquet(f"{ROOT}/lfe_extension/results/refine/rates/famday.parquet", columns=["fam", "date", "sset"])
    s15 = set(zip(fr[fr.sset == "S15"].fam, fr[fr.sset == "S15"].date))
    Rf["date"] = Rf.time_utc.dt.strftime("%Y-%m-%d")
    st = W.set_index("episode").start
    Rf["fwd_s15"] = np.array([(f, d) in s15 for f, d in zip(Rf.fam, Rf.date)]) & (Rf.ep.map(st).fillna("9999") < "2019-01-01").to_numpy()
    # backward
    X = pd.read_parquet(f"{OUT}/wave_targets.parquet")
    Mb = pd.read_parquet(f"{OUT}/mf_cn_backward.parquet")
    Mb = Mb.sort_values("reader").drop_duplicates(["ev_id", "sta"], keep="first")      # "cut" < "raw": cut store wins (identical anyway)
    Rb = Mb.merge(X[X.kind == "bw"][["ev_id", "fam", "time_utc", "seg", "ep"]], on="ev_id")
    Rb["head"] = False; Rb["src"] = "bw"; Rb["fwd_s15"] = False
    cols = ["ev_id", "sta", "a_mf", "cc_mf", "sig_noise", "fam", "time_utc", "head", "ep", "seg", "src", "fwd_s15"]
    R = pd.concat([Rf[cols], Rb[cols]], ignore_index=True)
    R["seg"] = R.seg.replace({"CN4f": "ERA"})
    R["neg"] = ~(R.a_mf > 0)
    R["logN"] = np.log10(R.sig_noise)
    R["logS"] = np.where(R.neg, R.logN, np.log10(R.a_mf.clip(lower=1e-30)))
    R["sn"] = R.a_mf / R.sig_noise
    R.loc[~np.isfinite(R.logN), ["logS", "logN"]] = np.nan
    R["epoch"] = [epoch_of(s, t) for s, t in zip(R.sta, R.time_utc)]
    wc = pd.read_csv(f"{S6}/window_check.csv")
    R = R.merge(wc[["fam", "sta", "window_ok"]], on=["fam", "sta"], how="left"); R["window_ok"] = R.window_ok.fillna(False)
    R = R.merge(F[["fid", "group", "dep"]].rename(columns={"fid": "fam"}), on="fam", how="left")
    return W, M, F, R


def best_stations(R):
    H = R[R["head"]]
    nh = H.groupby("fam").ev_id.nunique()
    FS = H.groupby(["fam", "sta"]).agg(med_sn=("sn", "median"), n=("logS", "count")).reset_index()
    FS["presence"] = FS.n / FS.fam.map(nh)
    wok = R.drop_duplicates(["fam", "sta"]).set_index(["fam", "sta"]).window_ok
    FS["eligible"] = [bool(wok.get((f, s), False)) for f, s in zip(FS.fam, FS.sta)]
    best, weak = {}, {}
    for f, g in FS.groupby("fam"):
        g = g[g.eligible].sort_values("med_sn", ascending=False)
        q = g[(g.med_sn >= 2) & (g.presence >= 0.5)]
        if len(q) >= 3:
            best[f] = q.sta.head(5).tolist(); weak[f] = False
        else:
            best[f] = g[g.presence >= 0.5].sta.head(3).tolist() or g.sta.head(3).tolist(); weak[f] = True
    FS["best"] = [s in best.get(f, []) for f, s in zip(FS.fam, FS.sta)]
    FS.to_csv(f"{OUT}/cn_family_stations.csv", index=False)
    return best, weak


def dm(B, col):
    ref = B[B["head"]].groupby(["fam", "sta"])[col].median().rename("ref")
    B = B.join(ref, on=["fam", "sta"]); B["r"] = B[col] - B["ref"]
    g = B.groupby("ev_id").r
    out = pd.DataFrame(dict(dM=g.median(), n_used=g.count()))
    out.loc[out.n_used < 2, "dM"] = np.nan
    return B, out


def pipeline(R, best):
    R = R.copy(); R["best"] = [s in best.get(f, []) for f, s in zip(R.fam, R.sta)]
    B = R[R.best & R.logS.notna()].copy()
    B1, D1 = dm(B, "logS")
    B1 = B1.sort_values("ev_id", kind="stable")
    res = []
    for ev, g in B1.groupby("ev_id", sort=True).r:
        v = g.values
        res.append(np.full(len(v), np.nan) if len(v) < 3 else np.array([v[i] - np.median(np.delete(v, i)) for i in range(len(v))]))
    B1["resid"] = np.concatenate(res)
    use = B1["head"] | (B1.src == "bw")
    Ep = B1[use].groupby(["sta", "epoch"]).resid.agg(["median", "count"]).reset_index().rename(columns={"median": "e_res"})
    Ep.loc[Ep["count"] < 20, "e_res"] = 0.0
    Ep.to_csv(f"{OUT}/cn_station_epoch_residuals.csv", index=False)
    B = B.merge(Ep[["sta", "epoch", "e_res"]], on=["sta", "epoch"], how="left"); B["e_res"] = B.e_res.fillna(0.0)
    B["logS_c"] = B.logS - B.e_res; B["logN_c"] = B.logN - B.e_res
    _, D2 = dm(B, "logS_c"); _, DN = dm(B, "logN_c")
    D = D2.join(D1.dM.rename("dM_senscorr_only")).join(DN.dM.rename("dN")).join(B.groupby("ev_id").sn.median().rename("sn_med"))
    return D, Ep


def maxc(z):
    z = np.asarray(z, float); z = z[np.isfinite(z)]
    if len(z) < 5:
        return np.nan
    h, ed = np.histogram(z, bins=np.arange(np.floor(z.min() * 10) / 10, z.max() + 0.2, 0.1))
    return ed[np.argmax(h)] + 0.05 + 0.1


def window_check_segments(R, best):
    """S-window check per segment on backward rows (median normalised envelope peak inside -1.5..+4.5 s, search -5..+4.5)."""
    X = pd.read_parquet(f"{OUT}/mf_cn_backward.parquet", columns=["ev_id", "sta"])
    env = np.load(f"{OUT}/env_cn_backward.npy")
    tg = -14.5 + (np.arange(env.shape[1]) + 0.5) * 0.1
    X["row"] = np.arange(len(X))
    X = X.merge(R[R.src == "bw"][["ev_id", "sta", "fam", "seg"]].drop_duplicates(["ev_id", "sta"]), on=["ev_id", "sta"])
    with np.errstate(invalid="ignore", divide="ignore"):
        En = env / np.nanmean(env[:, tg < -5], axis=1, keepdims=True)
    rows = []
    for (sg, f, s), g in X.groupby(["seg", "fam", "sta"]):
        if len(g) < 10:
            continue
        m = np.nanmedian(En[g.row.to_numpy()], axis=0)
        sel = tg >= -5
        pk = tg[sel][int(np.nanargmax(m[sel]))]
        rows.append(dict(seg=sg, fam=f, sta=s, n=len(g), peak_s=pk, ok=-1.5 <= pk <= 4.5, best=s in best.get(f, [])))
    Wc = pd.DataFrame(rows)
    Wc.to_csv(f"{OUT}/cn_window_check_segments.csv", index=False)
    return Wc.groupby("seg").agg(n=("ok", "size"), frac_ok=("ok", "mean"),
                                 frac_ok_best=("ok", lambda s: float(s[Wc.loc[s.index, "best"]].mean()))).reset_index()


def ts_slope_pairs(t, y, pairs):
    i, j = pairs
    dt = t[j] - t[i]; m = np.abs(dt) > 1e-9
    return float(np.median((y[j] - y[i])[m] / dt[m])) if m.any() else np.nan


def trend_test(E, col, families_min=20, span_min=10.0, nperm=NP, label=""):
    """per family Theil-Sen of col vs time (per decade); whole-episode time permutation within family."""
    fams, data = [], []
    for f, g in E.groupby("fam"):
        g = g[np.isfinite(g[col])]
        if len(g) < families_min:
            continue
        t = (g.time_utc - pd.Timestamp("2003-01-01")).dt.total_seconds().to_numpy() / (86400 * 3652.5)
        if t.max() - t.min() < span_min / 10:
            continue
        n = len(g)
        if n * (n - 1) / 2 > 20000:
            i = rng.integers(0, n, 20000); j = rng.integers(0, n, 20000); k = i != j; pairs = (i[k], j[k])
        else:
            pairs = np.triu_indices(n, 1)
        epi = g.ep.to_numpy(); ue = np.unique(epi)
        mid = {e: np.median(t[epi == e]) for e in ue}
        off = t - np.array([mid[e] for e in epi])
        fams.append(f); data.append((t, g[col].to_numpy(), pairs, epi, ue, mid, off))
    if not fams:
        return dict(label=label, n_families=0)
    sl = np.array([ts_slope_pairs(t, y, p) for (t, y, p, *_ ) in data])
    obs = float(np.nanmean(sl))
    nul = []
    for _ in range(nperm):
        s_ = []
        for (t, y, p, epi, ue, mid, off) in data:
            pe = rng.permutation(ue); mp = {a: mid[b] for a, b in zip(ue, pe)}
            tn = off + np.array([mp[e] for e in epi])
            s_.append(ts_slope_pairs(tn, y, p))
        nul.append(np.nanmean(s_))
    nul = np.array(nul)
    bs = [float(np.nanmean(rng.choice(sl, len(sl)))) for _ in range(NB)]
    return dict(label=label, col=col, n_families=len(fams), n_events=int(sum(len(d[0]) for d in data)), mean_slope_per_decade=obs,
                ci95=[float(np.quantile(bs, 0.025)), float(np.quantile(bs, 0.975))],
                p_perm=float((np.sum(np.abs(nul - np.nanmedian(nul)) >= abs(obs - np.nanmedian(nul))) + 1) / (nperm + 1)),
                null_mean=float(np.nanmean(nul)))


def pair_test(E, early_mask, late_mask, col, floor_by_fam, nperm=NP, label=""):
    """mean over families of [median late - median early] above the family floor; episode labels swapped between sides."""
    E = E.copy(); E["fl"] = E.fam.map(floor_by_fam) if floor_by_fam is not None else np.nan
    A = E[(early_mask | late_mask) & np.isfinite(E[col])].copy()
    A["late"] = late_mask[A.index]
    if floor_by_fam is not None:
        A = A[A.dM >= A.fl]
    def stat(lab):
        out = []
        for f, g in A.groupby("fam"):
            l = lab[g.index]
            a, b = g[col][~l], g[col][l]
            if len(a) >= 3 and len(b) >= 3:
                out.append(b.median() - a.median())
        return np.array(out)
    lab0 = A.late
    d = stat(lab0); obs = float(np.mean(d)) if len(d) else np.nan
    eps_late = A[A.late].ep.unique(); eps_early = A[~A.late].ep.unique(); alle = np.r_[eps_early, eps_late]
    nul = []
    for _ in range(nperm):
        p = rng.permutation(alle); lset = set(p[:len(eps_late)])
        nul.append(np.mean(stat(A.ep.isin(lset))) if True else np.nan)
    nul = np.array(nul, float)
    bs = [float(np.mean(rng.choice(d, len(d)))) for _ in range(NB)] if len(d) else [np.nan]
    return dict(label=label, col=col, n_families=int(len(d)), effect=obs, ci95=[float(np.nanquantile(bs, 0.025)), float(np.nanquantile(bs, 0.975))],
                p_perm=float((np.sum(np.abs(nul) >= abs(obs)) + 1) / (np.isfinite(nul).sum() + 1)), n_eps=[int(len(eps_early)), int(len(eps_late))])


def seen(t_m, t_n, t_all):
    ok = t_m.get("p_perm", 1) < 0.01 and not (t_m["ci95"][0] <= 0 <= t_m["ci95"][1])
    k = "mean_slope_per_decade" if "mean_slope_per_decade" in t_m else "effect"
    ok = ok and np.sign(t_m[k]) == np.sign(t_all[k]) and abs(t_n[k]) < abs(t_m[k])
    return bool(ok)


def main():
    W, M, F, R = load()
    log(f"rows: {R.groupby(['src', 'sta']).size().to_dict()}")
    best, weak = best_stations(R)
    log(f"best CN stations: families {len(best)}, weak {sum(weak.values())}; per-station use "
        f"{pd.Series([s for v in best.values() for s in v]).value_counts().to_dict()}")
    D, Ep = pipeline(R, best)
    log("station-epoch residuals (e_res, count):\n" + Ep.round(3).to_string())
    ev = R.drop_duplicates("ev_id").set_index("ev_id")[["fam", "time_utc", "head", "ep", "seg", "src", "fwd_s15", "group"]]
    E = ev.join(D, how="inner")
    E["z"] = E.dM - E.dN
    # completeness: Mc_z per family from forward headline (S6 rule), group fallback
    H = E[E["head"] & E.z.notna()]
    mcz = {f: maxc(g.z.values) for f, g in H.groupby("fam") if len(g) >= 30}
    mcg = {k: maxc(g.z.values) for k, g in H.groupby("group")}
    E["Mc_z"] = [mcz.get(f, mcg.get(gp, np.nan)) for f, gp in zip(E.fam, E.group)]
    # analysis events: inside majors only
    A = E[E.ep.notna()].copy()
    A["seg4"] = A.seg.replace({"B2f": "B2B3", "B3f": "B2B3"})
    A.loc[A.fwd_s15, "seg4b"] = "FWD_S15"
    A["seg4b"] = A.seg4b.fillna(A.seg4)
    Mc = A.groupby(["fam", "seg4"]).dN.median().rename("dNseg").reset_index()
    Mc["Mc_seg"] = Mc.fam.map(A.groupby("fam").Mc_z.first()) + Mc.dNseg
    nseg = A.groupby(["fam", "seg4"]).size().rename("n").reset_index()
    Mc = Mc.merge(nseg, on=["fam", "seg4"])
    Mcp = Mc.pivot(index="fam", columns="seg4", values="Mc_seg")
    Mc15 = A[A.fwd_s15].groupby("fam").dN.median() + A[A.fwd_s15].groupby("fam").Mc_z.first()
    A.reset_index().to_parquet(f"{OUT}/sizes_2003_2026.parquet", index=False)
    Mc.to_csv(f"{OUT}/completeness_by_segment.csv", index=False)
    summ = dict(n_events_majors=A.groupby("seg4").size().to_dict(), n_events_with_dM=A[A.dM.notna()].groupby("seg4").size().to_dict(),
                n_head_forward=int(E["head"].sum()), weak_families=int(sum(weak.values())))
    # S-window check per segment
    wseg = window_check_segments(R, best)
    summ["window_check_segments"] = wseg.round(3).to_dict("records")
    log(f"events: {summ}")
    # floors
    okseg = Mc[Mc.n >= 3].pivot(index="fam", columns="seg4", values="Mc_seg")
    floor_all = okseg.max(axis=1)
    A["floor_all"] = A.fam.map(floor_all)
    T = {}
    above = A[A.dM >= A.floor_all]
    T["T6a"] = trend_test(above, "dM", label="T6a dM above shared floor 2003-26")
    T["T6a_dN"] = trend_test(above, "dN", label="T6a dN (noise)")
    T["T6a_all"] = trend_test(A[A.dM.notna()], "dM", label="T6a all rows (no floor)")
    T["T6a_big"] = trend_test(above[above.ep.isin(W[W.big].episode)], "dM", label="T6a big ETS only")
    # T6b like-for-like: S15f vs FWD-S15
    fl_b = pd.concat([Mcp.get("S15f"), Mc15], axis=1).max(axis=1)
    early = A.seg4b == "S15f"; late = A.seg4b == "FWD_S15"
    T["T6b"] = pair_test(A, early, late, "dM", fl_b, label="T6b S15f vs FWD-S15")
    T["T6b_dN"] = pair_test(A, early, late, "dN", fl_b, label="T6b dN")
    T["T6b_all"] = pair_test(A.assign(dMx=A.dM), early, late, "dM", None, label="T6b all rows")
    # T6c PGC alone 2003-01-01 .. 2015-01-11 (one instrument, no correction)
    P = R[(R.sta == "CN.PGC") & (R.time_utc < "2015-01-12") & R.logS.notna()].copy()
    refP = R[(R.sta == "CN.PGC") & R["head"] & (R.time_utc < "2015-01-12")].groupby("fam").logS.median()
    refP_all = R[(R.sta == "CN.PGC") & R["head"]].groupby("fam").logS.median()
    P["ref"] = P.fam.map(refP).fillna(P.fam.map(refP_all))
    P["dM"] = P.logS - P.ref; P["dN"] = P.logN - P.ref
    P = P[P.ep.notna()]
    P["z"] = P.dM - P.dN
    mczP = {f: maxc(g.z.values) for f, g in P.groupby("fam") if len(g) >= 30}
    P["McP"] = P.fam.map(mczP)
    P["seg4"] = P.seg.replace({"B2f": "B2B3", "B3f": "B2B3"})
    flP = (P.groupby(["fam", "seg4"]).dN.median().unstack() .add(pd.Series(mczP), axis=0)).max(axis=1)
    P["floor"] = P.fam.map(flP)
    Pa = P[P.dM >= P.floor]
    T["T6c"] = trend_test(Pa, "dM", label="T6c PGC alone 2003-2015-01")
    T["T6c_dN"] = trend_test(Pa, "dN", label="T6c dN")
    T["T6c_all"] = trend_test(P, "dM", label="T6c all rows")
    # POST-HOC (PLAN Revision 2): T6c on the frozen-detector span only (no era detector), floor recomputed without ERA
    Pf = P[(P.time_utc >= "2006-07-01")].copy()
    flPf = (Pf.groupby(["fam", "seg4"]).dN.median().unstack().add(pd.Series(mczP), axis=0)).max(axis=1)
    Pf["floor"] = Pf.fam.map(flPf)
    T["T6c_frozen_posthoc"] = trend_test(Pf[Pf.dM >= Pf.floor], "dM", span_min=5.0, label="T6c PGC frozen detector 2006-07..2015-01 (post-hoc)")
    T["T6c_frozen_posthoc_dN"] = trend_test(Pf[Pf.dM >= Pf.floor], "dN", span_min=5.0, label="T6c frozen dN (post-hoc)")
    T["T6c_frozen_posthoc_all"] = trend_test(Pf, "dM", span_min=5.0, label="T6c frozen all rows (post-hoc)")
    # columns: segment medians relative to FWD above each pair's shared floor
    for sg in ["ERA", "B2B3", "S15f"]:
        if sg not in Mcp:
            continue
        fl = pd.concat([Mcp[sg], Mcp["FWD"]], axis=1).max(axis=1)
        T[f"seg_{sg}_vs_FWD"] = pair_test(A, A.seg4 == sg, A.seg4 == "FWD", "dM", fl, nperm=500, label=f"{sg} vs FWD")
        T[f"seg_{sg}_vs_FWD_dN"] = pair_test(A, A.seg4 == sg, A.seg4 == "FWD", "dN", fl, nperm=200, label=f"{sg} vs FWD dN")
    for k, v in T.items():
        log(f"{k}: {json.dumps(v, default=float)}")
    sa = seen(T["T6a"], T["T6a_dN"], T["T6a_all"]); sb = seen(T["T6b"], T["T6b_dN"], T["T6b_all"]); sc = seen(T["T6c"], T["T6c_dN"], T["T6c_all"])
    inside = T["T6a"]["ci95"][0] >= -0.05 and T["T6a"]["ci95"][1] <= 0.05
    if not sa and inside and not sb:
        v = "holds (no change)"
    elif sa and ((sb and np.sign(T["T6b"]["effect"]) == np.sign(T["T6a"]["mean_slope_per_decade"])) or
                 (sc and np.sign(T["T6c"]["mean_slope_per_decade"]) == np.sign(T["T6a"]["mean_slope_per_decade"]))):
        v = "changes"
    else:
        v = "undecided"
    T["seen"] = dict(T6a=sa, T6b=sb, T6c=sc, T6a_ci_inside_0p05=bool(inside)); T["verdict_3b"] = v; T["summary"] = summ
    log(f"seen {T['seen']} VERDICT 3b: {v}")
    # yearly medians per group (figure)
    above2 = above.copy(); above2["year"] = above2.time_utc.dt.year
    yr = above2.groupby(["group", "year"]).agg(dM=("dM", "median"), dN=("dN", "median"), n=("dM", "size")).reset_index()
    yr.to_csv(f"{OUT}/yearly_group_medians.csv", index=False)
    ally = A[A.dM.notna()].assign(year=lambda x: x.time_utc.dt.year).groupby(["group", "year"]).agg(
        dM=("dM", "median"), dN=("dN", "median"), n=("dM", "size")).reset_index()
    ally.to_csv(f"{OUT}/yearly_group_medians_allrows.csv", index=False)
    Pa.assign(year=Pa.time_utc.dt.year).groupby("year").agg(dM=("dM", "median"), dN=("dN", "median"), n=("dM", "size")).reset_index().to_csv(
        f"{OUT}/yearly_pgc_medians.csv", index=False)
    json.dump(T, open(f"{OUT}/f6_size.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


if __name__ == "__main__":
    main()
