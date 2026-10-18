#!/usr/bin/env python
"""S0b: episode-consistency re-check with cross-family de-duplication (svi_science/results/s0b/PLAN_S0b.md).

Reuses lfe_extension_wa/scripts/episode_rates.py as a module (EPR_OUT -> svi_science/results/s0b/run_dedup); that file is
not edited. Only the per (family, date) count files are rebuilt with de-duplication; everything else is the first run's.

  python s0_s0b_dedup.py verify     # keep_xfam8 definition checks on P_by_candidate (SVI)
  python s0_s0b_dedup.py nwdup      # WA-NW: duplicate rates at +-2/4/8 s, de-dup count file (+ raw rebuild check)
  python s0_s0b_dedup.py prep       # run folder: frozen inputs copied, SVI counts with keep_xfam8 == 1 (+ raw rebuild check)
  python s0_s0b_dedup.py run        # episode_rates table / checks / consistency (EPR_ROBUST=1) / figures in run folder
  python s0_s0b_dedup.py compare    # before / after table + verdict inputs
Writes only under svi_science/results/s0b/. Never deletes.
"""
from __future__ import annotations
import os, sys, json, glob, shutil, warnings
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
S0B = f"{ROOT}/svi_science/results/s0b"
RUN = f"{S0B}/run_dedup"
os.environ["EPR_OUT"] = RUN
os.environ["EPR_ROBUST"] = "1"
sys.path.insert(0, f"{ROOT}/lfe_extension_wa/scripts")
import episode_rates as ER          # noqa: E402  (reads EPR_OUT at import)
import wa_common as W               # noqa: E402

EPR = f"{ROOT}/lfe_extension_wa/results/episode_rates"
ROB = f"{EPR}/robustness"
PBC = f"{ROOT}/lfe_extension/results/refine/conf/P_by_candidate.parquet"
WAP = f"{ROOT}/lfe_extension_wa/results/wa_prod"
CUTS = ER.CUTS


def log(*a):
    os.makedirs(S0B, exist_ok=True)
    W.log(f"{S0B}/s0b.log", *a)


# ------------------------------------------------------------------------------------------------ cross-family keep (numba)
from numba import njit  # noqa: E402


@njit(cache=True)
def _keep_sorted(ts, vs, fs, share, w):
    """Same rule as lfe_extension/scripts/c_common.local_max_keep_cc on time-sorted input; also returns the family of the
    strongest suppressor (-1 if kept)."""
    n = len(ts)
    keep = np.ones(n, np.bool_)
    win = np.full(n, -1, np.int64)
    lo = 0
    for i in range(n):
        while ts[lo] < ts[i] - w:
            lo += 1
        best = -1.0e300
        j = lo
        while j < n and ts[j] <= ts[i] + w:
            if j != i and share[fs[i], fs[j]]:
                if vs[j] > vs[i] or (vs[j] == vs[i] and j < i):
                    keep[i] = False
                    if vs[j] > best:
                        best = vs[j]; win[i] = fs[j]
            j += 1
    return keep, win


def keep_cc(t, v, fam, share, w):
    o = np.argsort(t, kind="stable")
    k, wf = _keep_sorted(t[o].astype(np.float64), v[o].astype(np.float64), fam[o].astype(np.int64), share, float(w))
    keep = np.zeros(len(t), bool); win = np.full(len(t), -1, np.int64)
    keep[o] = k; win[o] = wf
    return keep, win


# ------------------------------------------------------------------------------------------------ verify (SVI keep_xfam8)
def verify(a):
    sys.path.insert(0, f"{ROOT}/lfe_extension/scripts")
    import c_common as C
    share = C.share_ok_matrix()
    P = pd.read_parquet(PBC, columns=["run", "date", "series", "fam", "t", "P_same", "keep_xfam8", "net_mad"])
    rep = {"definition_source": "lfe_extension/scripts/svi_ref_mf.py:474-477 -> c_common.local_max_keep_cc(t, net_cc, fam, "
                                "share_ok_matrix(), w=2.0), per run x day x series (data and c1 alike); '8' = Merlin E5r-audit "
                                "finding F8, not an 8-s window"}
    # empirical: does every keep==0 row have a station-sharing other-family candidate within +-2 s (same run/date/series)?
    rows = []
    for (run, ser), g in P.groupby(["run", "series"]):
        g = g.sort_values("t")
        t = g.t.to_numpy(); f = g.fam.to_numpy().astype(int); k = g.keep_xfam8.to_numpy() > 0.5
        lo = np.searchsorted(t, t - 2.0, "left"); hi = np.searchsorted(t, t + 2.0, "right")
        has2 = np.zeros(len(t), bool); has8 = np.zeros(len(t), bool)
        lo8 = np.searchsorted(t, t - 8.0, "left"); hi8 = np.searchsorted(t, t + 8.0, "right")
        for i in np.flatnonzero((hi - lo > 1) | (hi8 - lo8 > 1)):
            for j in range(lo8[i], hi8[i]):
                if j != i and f[j] != f[i] and share[f[i], f[j]]:
                    has8[i] = True
                    if abs(t[j] - t[i]) <= 2.0:
                        has2[i] = True
        hp = g.P_same.to_numpy() >= 0.75
        rows.append(dict(run=run, series=ser, n=len(g), n_keep0=int((~k).sum()),
                         keep0_with_nbr2=int((~k & has2).sum()), keep0_without_nbr2=int((~k & ~has2).sum()),
                         keep1_with_nbr2=int((k & has2).sum()),
                         p075_n=int(hp.sum()), p075_keep0=int((hp & ~k).sum()), p075_keep0_with_nbr2=int((hp & ~k & has2).sum()),
                         p075_keep1_with_nbr2_to_8_only=int((hp & k & has8 & ~has2).sum())))
    V = pd.DataFrame(rows)
    V["p075_keep0_frac"] = V.p075_keep0 / V.p075_n
    V.to_csv(f"{S0B}/keep_xfam8_check.csv", index=False)
    rep["table"] = V.to_dict("records")
    json.dump(rep, open(f"{S0B}/keep_xfam8_check.json", "w"), indent=1, default=str)
    log("[verify]", V.to_string())
    print(V.to_string())


# ------------------------------------------------------------------------------------------------ WA-NW duplicates
def _wa_share(F, T, core):
    """share[era] bool matrix (fid-indexed): |A & B| / min(|A|, |B|) >= 0.5 over tpl_library sets ∩ era core."""
    nb = int(F.fid.max()) + 2
    out = {}
    for era in ("E1", "E2", "E3"):
        st = T[T.netsta.isin(core[era])].groupby("fid").netsta.apply(set).to_dict()
        S = np.zeros((nb, nb), np.bool_)
        for fa, A in st.items():
            for fb, B in st.items():
                S[fa, fb] = len(A & B) / max(min(len(A), len(B)), 1) >= 0.5
        out[era] = S
    return out


_G = {}


def _nw_file(path):
    F, share, nwset = _G["F"], _G["share"], _G["nw"]
    D = pd.read_parquet(path, columns=["date", "net", "series", "fid", "t", "net_cc", "P_own_WAVE_PNET", "tpl_insample"])
    D = D[D.net == "core"].reset_index(drop=True)
    if not len(D):
        return None
    for w in (2.0, 4.0, 8.0):
        D[f"keep{int(w)}"] = True
    D["keepNW2"] = True; D["win2"] = -1
    for (dt, ser), g in D.groupby(["date", "series"]):
        era = ER.wa_era(dt)
        S = share[era]
        ix = g.index.to_numpy()
        t = g.t.to_numpy(); v = g.net_cc.to_numpy(float); f = g.fid.to_numpy().astype(int)
        for w in (2.0, 4.0, 8.0):
            k, wf = keep_cc(t, v, f, S, w)
            D.loc[ix, f"keep{int(w)}"] = k
            if w == 2.0:
                D.loc[ix, "win2"] = wf
        m = np.isin(f, list(nwset))
        if m.any():
            k, _ = keep_cc(t[m], v[m], f[m], S, 2.0)
            D.loc[ix[m], "keepNW2"] = k
    D = D[D.fid.isin(nwset)].copy()
    D["half"] = W.half_of(D.t.to_numpy()); D["ins"] = D.tpl_insample.fillna(False).astype(bool)
    p = D.P_own_WAVE_PNET.to_numpy()
    out = {}
    for tag, km in (("raw", None), ("dedup", "keep2")):
        X = D if km is None else D[D[km]]
        X = X.copy(); pp = X.P_own_WAVE_PNET.to_numpy()
        for k_, c in CUTS.items():
            X[f"n{k_}"] = (pp >= c).astype(np.int32)
        out[tag] = X.groupby(["fid", "date", "series", "half", "ins"]).agg(
            sumP=("P_own_WAVE_PNET", "sum"), **{f"n{k_}": (f"n{k_}", "sum") for k_ in CUTS}).reset_index()
    hp = D[p >= 0.75][["date", "series", "fid", "t", "net_cc", "P_own_WAVE_PNET", "keep2", "keep4", "keep8", "keepNW2", "win2"]]
    return out["raw"], out["dedup"], hp


def nwdup(a):
    from multiprocessing import Pool
    F = pd.read_csv(f"{WAP}/setup/families.csv")
    T = pd.read_csv(f"{WAP}/setup/tpl_library.csv")
    core = json.load(open(f"{WAP}/setup/core_sets.json"))
    _G.update(F=F, share=_wa_share(F, T, core), nw=set(F[F.reg == "NW"].fid.astype(int)))
    # numba rule == c_common rule on one sample month (exact)
    sys.path.insert(0, f"{ROOT}/lfe_extension/scripts")
    import c_common as C
    D = pd.read_parquet(f"{WAP}/conf/P/prod/P_2016-01.parquet", columns=["date", "net", "series", "fid", "t", "net_cc"])
    D = D[(D.net == "core") & (D.date == "2016-01-10") & (D.series == "data")]
    S = _G["share"]["E3"]
    k1 = C.local_max_keep_cc(D.t.to_numpy(), D.net_cc.to_numpy(float), D.fid.to_numpy().astype(int), S, 2.0)
    k2, _ = keep_cc(D.t.to_numpy(), D.net_cc.to_numpy(float), D.fid.to_numpy().astype(int), S, 2.0)
    log(f"[nwdup] numba vs c_common keep rule on 2016-01-10 data: {len(D)} rows, {int((k1 != k2).sum())} differ")
    assert (k1 == k2).all()
    fs = sorted(glob.glob(f"{WAP}/conf/P/prod/P_*.parquet"))
    with Pool(16) as pool:
        parts = [x for x in pool.imap(_nw_file, fs) if x is not None]
    RAW = pd.concat([p[0] for p in parts], ignore_index=True)
    DED = pd.concat([p[1] for p in parts], ignore_index=True)
    HP = pd.concat([p[2] for p in parts], ignore_index=True)
    os.makedirs(f"{S0B}/nw", exist_ok=True)
    DED.to_parquet(f"{S0B}/nw/counts_wa_nw_fid_date_half_dedup.parquet", index=False)
    HP.to_parquet(f"{S0B}/nw/nw_p075_candidates_keep.parquet", index=False)
    # raw rebuild == first run's count file?
    O = pd.read_parquet(f"{EPR}/counts_wa_nw_fid_date_half.parquet")
    key = ["fid", "date", "series", "half", "ins"]
    J = O.merge(RAW, on=key, how="outer", suffixes=("_o", "_m"), indicator=True)
    diff = {c: int((J[f"{c}_o"].fillna(0) != J[f"{c}_m"].fillna(0)).sum()) for c in [f"n{k}" for k in CUTS]}
    diff["sumP_absdiff_max"] = float((J.sumP_o.fillna(0) - J.sumP_m.fillna(0)).abs().max())
    diff["rows_only_one_side"] = int((J._merge != "both").sum())
    # duplicate rates (data and controls), NW families, P >= 0.75
    nwF = F[F.reg == "NW"].set_index("fid")
    HP["era"] = [ER.wa_era(d) for d in HP.date]
    rep = dict(raw_rebuild_vs_first_run=diff, rule="drop iff a station-sharing core-net family candidate within +-w s has "
               "higher raw net_cc (ties: earlier), per day x series; share = overlap >= 0.5 of tpl_library sets ∩ era core")
    tab = []
    for era in ("E1", "E2", "E3", "all"):
        for ser in ("data", "rev", "c1"):
            h = HP[(HP.series == ser) & ((HP.era == era) if era != "all" else True)]
            if not len(h):
                continue
            r = dict(era=era, series=ser, n_p075=len(h))
            for c in ("keep2", "keep4", "keep8", "keepNW2"):
                r[f"dup_frac_{c}"] = float(1 - h[c].mean())
                r[f"n_dup_{c}"] = int((~h[c]).sum())
            tab.append(r)
    TB = pd.DataFrame(tab)
    TB.to_csv(f"{S0B}/nw/nw_duplicate_rates.csv", index=False)
    # winner distance / region for the E3 data duplicates at +-2 s
    d3 = HP[(HP.series == "data") & (HP.era == "E3") & ~HP.keep2].copy()
    allF = F.set_index("fid")
    d3["win_reg"] = d3.win2.map(allF.reg)
    la = allF.lat_use; lo = allF.lon_use
    d3["win_km"] = W.km(d3.fid.map(la).to_numpy(float), d3.fid.map(lo).to_numpy(float),
                        d3.win2.map(la).to_numpy(float), d3.win2.map(lo).to_numpy(float))
    rep["E3_data_dups_by_winner_region"] = d3.win_reg.value_counts().to_dict()
    rep["E3_data_dups_winner_km_quantiles"] = {str(q): float(np.nanquantile(d3.win_km, q)) for q in (0.1, 0.25, 0.5, 0.75, 0.9)}
    rep["E3_data_dups_winner_km_bins"] = pd.cut(d3.win_km, [0, 3, 6, 10, 20, 40, 1000]).value_counts().sort_index().astype(int).rename(str).to_dict()
    e3 = TB[(TB.era == "E3") & (TB.series == "data")].iloc[0]
    rep["E3_data_dup_frac_2s"] = float(e3.dup_frac_keep2)
    rep["material_gt_5pct"] = bool(e3.dup_frac_keep2 > 0.05)
    rep["table"] = TB.to_dict("records")
    json.dump(rep, open(f"{S0B}/nw/nw_duplicates.json", "w"), indent=1, default=str)
    log("[nwdup]", json.dumps({k: v for k, v in rep.items() if k != "table"}, default=str))
    print(TB.to_string())


# ------------------------------------------------------------------------------------------------ prep run folder
def _count_svi(P):
    p = P.P_same.to_numpy()
    P = P.copy()
    for k, c in CUTS.items():
        P[f"n{k}"] = (p >= c).astype(np.int32)
    g = P.groupby(["run", "fam", "date", "series"]).agg(sumP=("P_same", "sum"), **{f"n{k}": (f"n{k}", "sum") for k in CUTS}).reset_index()
    return g.rename(columns={"fam": "fid"})


def prep(a):
    os.makedirs(RUN, exist_ok=True)
    for f in ("episodes_svi_final.csv", "episodes_nw_final.csv", "exposure_svi.parquet", "exposure_wa_nw.parquet",
              "svi_day_stations.csv", "wa_core_famday.parquet"):
        shutil.copy2(f"{EPR}/{f}", f"{RUN}/{f}")
    P = pd.read_parquet(PBC, columns=["run", "date", "series", "fam", "t", "P_same", "sset_run", "keep_xfam8"])
    raw = _count_svi(P)
    O = pd.read_parquet(f"{EPR}/counts_svi_fid_date.parquet")
    key = ["run", "fid", "date", "series"]
    J = O.merge(raw, on=key, how="outer", suffixes=("_o", "_m"), indicator=True)
    chk = {c: int((J[f"{c}_o"].fillna(-1) != J[f"{c}_m"].fillna(-1)).sum()) for c in [f"n{k}" for k in CUTS]}
    chk["rows_only_one_side"] = int((J._merge != "both").sum())
    ded = _count_svi(P[P.keep_xfam8 > 0.5])
    ded.to_parquet(f"{RUN}/counts_svi_fid_date.parquet", index=False)
    nwmode = a.nw if hasattr(a, "nw") and a.nw else "dedup"
    if nwmode == "dedup":
        shutil.copy2(f"{S0B}/nw/counts_wa_nw_fid_date_half_dedup.parquet", f"{RUN}/counts_wa_nw_fid_date_half.parquet")
    else:
        shutil.copy2(f"{EPR}/counts_wa_nw_fid_date_half.parquet", f"{RUN}/counts_wa_nw_fid_date_half.parquet")
    rep = dict(svi_raw_rebuild_vs_first_run=chk, nw_counts=nwmode,
               svi_n075=raw.groupby(["run", "series"]).n075.sum().to_dict(), svi_n075_dedup=ded.groupby(["run", "series"]).n075.sum().to_dict(),
               svi_n050=raw.groupby(["run", "series"]).n050.sum().to_dict(), svi_n050_dedup=ded.groupby(["run", "series"]).n050.sum().to_dict())
    rep = {k: ({"|".join(kk) if isinstance(kk, tuple) else str(kk): int(vv) for kk, vv in v.items()} if isinstance(v, dict) else v)
           for k, v in rep.items()}
    json.dump(rep, open(f"{S0B}/prep_checks.json", "w"), indent=1, default=str)
    log("[prep]", json.dumps(rep, default=str))
    print(json.dumps(rep, indent=1, default=str))


def run(a):
    class A_: pass
    aa = A_(); aa.rule = ""; aa.tag = ""; aa.nproc = ""
    ER.table(aa); ER.checks(aa); ER.consistency(aa); ER.figures(aa)
    log("[run] table/checks/consistency/figures done in", RUN)


# ------------------------------------------------------------------------------------------------ compare
BINS = [0, 3, 6, 10, 20, 40]


def _summ(folder, block, variant):
    S = json.load(open(f"{folder}/consistency_summary.json"))
    k = f"{block}_{variant}"
    s = S[k]; m = s.get("mantel", {})
    P = pd.read_csv(f"{folder}/residual_pair_correlations.csv")
    q = P[(P.block == block) & (P.variant == variant)].copy()
    q["b"] = pd.cut(q.dist_km, BINS)
    cb = q.groupby("b", observed=True)["corr"].agg(["mean", "size"])
    r = dict(n_episodes=s["n_episodes"], tested=s["n_tested"], over_dispersed=s["n_overdispersed_q05"], consistent=s["n_consistent"],
             median_chi2_df=round(s["median_chi2_df"], 3), n_pairs=m.get("n_pairs"),
             coef_logdist=m.get("coef_logdist"), p_logdist=m.get("p_logdist"),
             coef_carrier=m.get("coef_jaccard"), p_carrier=m.get("p_jaccard"))
    for iv, row in cb.iterrows():
        r[f"corr_{int(iv.left)}_{int(iv.right)}km"] = round(float(row["mean"]), 3)
        r[f"npairs_{int(iv.left)}_{int(iv.right)}km"] = int(row["size"])
    return r


def compare(a):
    rows = []
    blocks = [("SVI_S13fixed_S13fixed", "SVI fixed-S13"), ("NW_E3", "WA-NW E3"), ("SVI_eff_S13", "SVI eff S13 (secondary)")]
    for blk, lab in blocks:
        for variant, cut in (("main", 0.75), ("cut_050", 0.5)):
            for when, folder in (("before (raw)", ROB), ("after (de-dup)", RUN)):
                r = dict(block=lab, cut=cut, when=when); r.update(_summ(folder, blk, variant)); rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(f"{S0B}/before_after.csv", index=False)
    # verdict per pre-registered criterion (0.75 de-dup)
    ver = {}
    for blk, lab in blocks[:2]:
        r = T[(T.block == lab) & (T.cut == 0.75) & (T.when == "after (de-dup)")].iloc[0]
        i = (r.coef_logdist < 0) and (r.p_logdist < 0.01)
        ii = (r["corr_3_6km"] > 0) and (r["corr_6_10km"] > 0)
        ver[lab] = dict(distance_negative_p_lt_001=bool(i), bins_3_6_and_6_10_positive=bool(ii),
                        verdict="survives" if (i and ii) else ("only co-located (< 3 km) coupling survives" if i else "does not survive"))
    json.dump(ver, open(f"{S0B}/verdict.json", "w"), indent=1)
    # same-event check, SVI: share of de-dup-removed P>=0.75 repeats and of residual pair corr change by bin
    log("[compare]", json.dumps(ver))
    pd.set_option("display.width", 250)
    print(T.T.to_string()); print(json.dumps(ver, indent=1))


# ------------------------------------------------------------------------------------------------ coincidence diagnostic
@njit(cache=True)
def _nbr_count(tq, fq, tp, fp, share, w):
    """for each query (tq, fq): number of pool rows (tp sorted, fp) of a station-sharing OTHER family within +-w s."""
    out = np.zeros(len(tq), np.int64)
    for i in range(len(tq)):
        lo = np.searchsorted(tp, tq[i] - w); hi = np.searchsorted(tp, tq[i] + w, side="right")
        c = 0
        for j in range(lo, hi):
            if fp[j] != fq[i] and share[fq[i], fp[j]]:
                c += 1
        out[i] = c
    return out


def _coinc_rows(Q, Pool_, share, km_of):
    """Q: P >= 0.75 rows (t, fam); Pool_: all candidates of the same series (t, fam). Fraction with a station-sharing
    other-family candidate within +-2 s, at true times and with the pool shifted by +-30/60/120 s (rate-preserving
    random-coincidence level); plus the distance to the nearest-in-time such neighbour."""
    tp = Pool_.t.to_numpy(float); o = np.argsort(tp); tp = tp[o]; fp = Pool_.fam.to_numpy().astype(np.int64)[o]
    tq = Q.t.to_numpy(float); fq = Q.fam.to_numpy().astype(np.int64)
    r = {}
    for sh in (0, 30, -30, 60, -60, 120, -120):
        r[f"frac_nbr2_shift{sh:+d}s" if sh else "frac_nbr2_true"] = float((_nbr_count(tq + sh, fq, tp, fp, share, 2.0) > 0).mean())
    return r


def coinc(a):
    """Diagnostic (not pre-registered): is the +-2 s cross-family 'duplicate' the same event, or a random coincidence
    of two candidates during busy days? Compare the true-time neighbour rate with time-shifted pools."""
    sys.path.insert(0, f"{ROOT}/lfe_extension/scripts")
    import c_common as C
    out = []
    # SVI
    share = C.share_ok_matrix()
    P = pd.read_parquet(PBC, columns=["run", "series", "fam", "t", "P_same", "keep_xfam8"])
    E = pd.read_csv(ER.FAMS12).set_index("bkey")
    for (run, ser), g in P.groupby(["run", "series"]):
        q = g[g.P_same >= 0.75]
        r = dict(region="SVI", run=run, series=ser, n_p075=len(q), dup_frac_keep=float(1 - q.keep_xfam8.mean()))
        r.update(_coinc_rows(q, g, share, None)); out.append(r)
    # SVI: distance from each P>=0.75 keep==0 data repeat (run prod) to the station-sharing neighbours within +-2 s
    g = P[(P.run == "prod") & (P.series == "data")].sort_values("t")
    t = g.t.to_numpy(); f = g.fam.to_numpy().astype(int); hp = (g.P_same.to_numpy() >= 0.75) & (g.keep_xfam8.to_numpy() < 0.5)
    lo = np.searchsorted(t, t - 2.0); hi = np.searchsorted(t, t + 2.0, "right")
    dd = []
    for i in np.flatnonzero(hp):
        for j in range(lo[i], hi[i]):
            if j != i and f[j] != f[i] and share[f[i], f[j]]:
                dd.append(W.km(E.lat[f[i]], E.lon[f[i]], E.lat[f[j]], E.lon[f[j]]))
    dd = np.array(dd, float)
    svi_d = pd.cut(dd, [-0.01, 3, 6, 10, 20, 40, 1000]).value_counts().sort_index()
    # NW
    F = pd.read_csv(f"{WAP}/setup/families.csv"); T = pd.read_csv(f"{WAP}/setup/tpl_library.csv")
    core = json.load(open(f"{WAP}/setup/core_sets.json")); sh = _wa_share(F, T, core)["E3"]
    nw = set(F[F.reg == "NW"].fid.astype(int))
    parts = []
    for p in sorted(glob.glob(f"{WAP}/conf/P/prod/P_*.parquet")):
        if os.path.basename(p)[2:9] < "2015-01":
            continue
        D = pd.read_parquet(p, columns=["net", "series", "fid", "t", "P_own_WAVE_PNET"])
        parts.append(D[D.net == "core"].drop(columns=["net"]))
    D = pd.concat(parts, ignore_index=True).rename(columns={"fid": "fam"})
    for ser, g in D.groupby("series"):
        q = g[(g.P_own_WAVE_PNET >= 0.75) & g.fam.isin(nw)]
        r = dict(region="NW_E3", run="prod", series=ser, n_p075=len(q))
        r.update(_coinc_rows(q, g, sh, None)); out.append(r)
    O = pd.DataFrame(out)
    O.to_csv(f"{S0B}/coincidence_diagnostic.csv", index=False)
    J = dict(svi_prod_data_keep0_p075_neighbour_km_bins={str(k): int(v) for k, v in svi_d.items()},
             svi_neighbour_km_median=float(np.median(dd)) if len(dd) else None)
    json.dump(J, open(f"{S0B}/coincidence_diagnostic.json", "w"), indent=1)
    log("[coinc]", O.to_string(), json.dumps(J))
    print(O.to_string()); print(J)


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("what"); ap.add_argument("--nw", default="")
    a = ap.parse_args()
    globals()[a.what](a)
