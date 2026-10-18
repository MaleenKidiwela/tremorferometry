#!/usr/bin/env python
"""Per-family x per-ETS-episode LFE rate consistency (WA-NW + southern Vancouver Island).
Spec: results/episode_rates/MERLIN_DESIGN_2026-10-06.md; frozen rules: results/episode_rates/PLAN.md.
Writes only under results/episode_rates/. Never deletes.

  python scripts/episode_rates.py exposure      # step 2: exposure tables + checks
  python scripts/episode_rates.py episodes      # step 3a: episode lists + acceptance check
  python scripts/episode_rates.py count_wa      # per (fid, date) counts from conf/P/prod (NW core, own_pnet), halves
  python scripts/episode_rates.py count_svi     # per (fam, date, run) counts from P_by_candidate
  python scripts/episode_rates.py table         # per (family, episode) table
  python scripts/episode_rates.py checks        # sanity checks (a)-(c)
  python scripts/episode_rates.py consistency   # step 5
  python scripts/episode_rates.py figures
"""
from __future__ import annotations
import os, sys, json, glob, argparse, warnings
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import wa_common as W

ROOT = W.ROOT
OUT = os.environ.get("EPR_OUT", f"{W.RES}/episode_rates")   # robustness runs write elsewhere
WAP = f"{W.RES}/wa_prod"
SVI = f"{ROOT}/lfe_extension/results/refine"
PNSN = f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv"
FAMS12 = f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv"
WECC = [("2007-01-12", "2007-02-02"), ("2008-04-30", "2008-06-02"), ("2009-04-28", "2009-06-03")]  # wa_prod_post.WECC
COV_WA = ("2009-08-06", "2026-05-27")
COV_SVI = ("2010-01-06", "2026-05-27")
SVI_T = ("2014-01-01", "2026-09-15")
WA_T = ("2007-01-01", "2026-09-15")
LOGF = f"{OUT}/episode_rates.log"

# frozen episode rule (PLAN.md sec. 2); revisions are recorded in PLAN.md and set here
RULE = dict(frac=0.05, gap=5, minlen=5, major=0.15)   # PLAN.md Revision 1 (frozen 0.10/3/5/0.25)
RULE_FROZEN = dict(frac=0.10, gap=3, minlen=5, major=0.25)


def log(*a):
    os.makedirs(OUT, exist_ok=True)
    W.log(LOGF, *a)


# ------------------------------------------------------------------------------------------------ families
def wa_families():
    F = pd.read_csv(f"{WAP}/setup/families.csv")
    F = F[F.reg == "NW"].copy()
    return F.rename(columns={"lat_use": "lat", "lon_use": "lon"})[["fid", "lat", "lon", "family_uid"]].reset_index(drop=True)


def svi_families():
    E = pd.read_csv(FAMS12)
    fams = pd.read_parquet(f"{SVI}/conf/P_by_candidate.parquet", columns=["fam"]).fam.unique()
    E = E[E.bkey.isin(fams)].rename(columns={"bkey": "fid"})
    return E[["fid", "lat", "lon"]].reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ step 2 exposure
def exposure(a):
    os.makedirs(OUT, exist_ok=True)
    # WA: existing table, restricted to NW
    F = wa_families()
    X = pd.read_parquet(f"{WAP}/post/tremor_family_day.parquet")
    Xn = X[X.fid.isin(F.fid)].copy()
    rep = {"wa_table_families_nw": int(Xn.fid.nunique()), "wa_nw_families": len(F),
           "wa_table_date_min": Xn.date.min(), "wa_table_date_max": Xn.date.max(),
           "wa_table_rows_nw": len(Xn), "wa_table_rows_n30_pos": int((Xn.n30 > 0).sum())}
    # check against rates_daily_v2 n_trem30 (label own_pnet, NW)
    R = pd.read_parquet(f"{WAP}/post_v2/rates_daily_v2.parquet", columns=["unit_kind", "label", "fid", "date", "region", "n_trem30"])
    R = R[(R.unit_kind == "family") & (R.label == "own_pnet") & (R.region == "NW")]
    M = R.merge(Xn[["fid", "date", "n30"]], on=["fid", "date"], how="left")
    cov = (M.date >= COV_WA[0]) & (M.date <= COV_WA[1])
    Mc = M[cov]
    diff = (Mc.n_trem30.fillna(-1) != Mc.n30.fillna(0))
    rep.update(rates_rows_covered=int(cov.sum()), rates_rows_n30_mismatch=int(diff.sum()),
               rates_rows_offcov_with_ntrem=int(M[~cov].n_trem30.notna().sum()))
    Xn[["fid", "date", "n30"]].to_parquet(f"{OUT}/exposure_wa_nw.parquet", index=False)
    # SVI: build
    S = svi_families()
    T = pd.read_csv(PNSN, usecols=["time", "lat", "lon"])
    T["date"] = T.time.astype(str).str[:10]
    lat0, lat1, lon0, lon1 = S.lat.min() - 0.5, S.lat.max() + 0.5, S.lon.min() - 0.7, S.lon.max() + 0.7
    Tb = T[(T.lat > lat0) & (T.lat < lat1) & (T.lon > lon0) & (T.lon < lon1)]
    la, lo = Tb.lat.to_numpy(float), Tb.lon.to_numpy(float)
    rows = []
    for r in S.itertuples():
        d = W.km(la, lo, r.lat, r.lon)
        c = Tb[d <= 30].groupby("date").size()
        rows += [(int(r.fid), dd, int(v)) for dd, v in c.items()]
    XS = pd.DataFrame(rows, columns=["fid", "date", "n30"])
    XS.to_parquet(f"{OUT}/exposure_svi.parquet", index=False)
    # PNSN coverage over SVI by year: epicentres within 30 km of any SVI family, and north of 48.5 N
    dmin = np.full(len(Tb), np.inf)
    for r in S.itertuples():
        dmin = np.minimum(dmin, W.km(la, lo, r.lat, r.lon))
    Tb = Tb.assign(dmin=dmin, year=Tb.date.str[:4])
    cy = Tb[Tb.dmin <= 30].groupby("year").agg(epicentres=("date", "size"), days=("date", "nunique"),
                                                 north_of_48p5=("lat", lambda x: int((x > 48.5).sum())),
                                                 max_lat=("lat", "max"))
    fam_y = XS.assign(year=XS.date.str[:4]).groupby("year").agg(families_exposed=("fid", "nunique"), family_days=("fid", "size"))
    cy = cy.join(fam_y, how="outer")
    cy.to_csv(f"{OUT}/svi_pnsn_coverage_by_year.csv")
    # same for WA NW
    TW = T[(T.lat > 46.3) & (T.lat < 49.2) & (T.lon > -125.3) & (T.lon < -121.5)]
    lw, ow = TW.lat.to_numpy(float), TW.lon.to_numpy(float)
    dw = np.full(len(TW), np.inf)
    for r in F.itertuples():
        dw = np.minimum(dw, W.km(lw, ow, r.lat, r.lon))
    cw = TW[dw <= 30].assign(year=lambda x: x.date.str[:4]).groupby("year").agg(epicentres=("date", "size"), days=("date", "nunique"))
    cw.to_csv(f"{OUT}/wa_nw_pnsn_coverage_by_year.csv")
    rep.update(svi_families=len(S), svi_family_days_exposed=len(XS), svi_first=XS.date.min(), svi_last=XS.date.max())
    json.dump(rep, open(f"{OUT}/exposure_checks.json", "w"), indent=1, default=str)
    log("[exposure]", json.dumps(rep, default=str))
    print(cy.to_string()); print(cw.to_string())


# ------------------------------------------------------------------------------------------------ step 3a episodes
def find_episodes(Xe, nfam, t0, t1, rule):
    """Xe: fid, date, n30 (exposure > 0 rows). -> regional episode table (MAJOR flagged)."""
    days = pd.date_range(t0, t1, freq="D")
    ds = days.strftime("%Y-%m-%d")
    nexp = Xe[Xe.n30 > 0].groupby("date").fid.nunique().reindex(ds, fill_value=0).to_numpy()
    q = nexp >= rule["frac"] * nfam
    idx = np.flatnonzero(q)
    eps = []
    if len(idx):
        s = idx[0]; prev = idx[0]
        for i in idx[1:]:
            if i - prev - 1 <= rule["gap"]:
                prev = i; continue
            eps.append((s, prev)); s = i; prev = i
        eps.append((s, prev))
    rows = []
    for s, e in eps:
        L = e - s + 1
        if L < rule["minlen"]:
            continue
        fed = int(nexp[s:e + 1].sum())
        rows.append(dict(start=ds[s], end=ds[e], span_days=L, qual_days=int(q[s:e + 1].sum()), family_exposure_days=fed,
                         max_families_exposed=int(nexp[s:e + 1].max()), major_threshold=rule["major"] * nfam * 14,
                         major=fed >= rule["major"] * nfam * 14, kind="pnsn"))
    return pd.DataFrame(rows)


def reference_entries():
    T = pd.read_csv(f"{ROOT}/data/extension_stations/tremor_episodes.csv")
    nw = T[(T.region == "NW_WA") & ((T.source != "PNSN") | (T.n_strong_days >= 14))].copy()
    nw = nw[nw.start <= COV_WA[1]]
    S = pd.read_csv(f"{ROOT}/fault_tomography/inversion/res_catalog_sse/tremor_episodes.csv")
    sv = S[S.band == 48.5].rename(columns={"onset": "start", "peak": "end"})[["key", "start", "end"]].copy()
    sv = sv[(sv.end >= SVI_T[0]) & (sv.start <= COV_SVI[1])]
    return nw[["region", "start", "end", "n_strong_days", "n_windows", "source"]].reset_index(drop=True), sv.reset_index(drop=True)


def acceptance(E, ref):
    M = E[E.major]
    out = []
    for r in ref.itertuples():
        hit = M[(M.start <= r.end) & (M.end >= r.start)]
        out.append(dict(ref_start=r.start, ref_end=r.end, in_major=len(hit) > 0,
                        major_episode=";".join(f"{a}..{b}" for a, b in zip(hit.start, hit.end))))
    return pd.DataFrame(out)


def episodes(a):
    rule = dict(RULE)
    if a.rule:
        rule.update(json.loads(a.rule))
    Xw = pd.read_parquet(f"{OUT}/exposure_wa_nw.parquet")
    Xs = pd.read_parquet(f"{OUT}/exposure_svi.parquet")
    nw_n = len(wa_families()); sv_n = len(svi_families())
    Ew = find_episodes(Xw, nw_n, COV_WA[0], COV_WA[1], rule)
    wecc = pd.DataFrame([dict(start=s, end=e, span_days=(pd.Timestamp(e) - pd.Timestamp(s)).days + 1, qual_days=np.nan,
                              family_exposure_days=np.nan, max_families_exposed=np.nan, major_threshold=np.nan, major=True,
                              kind="wecc") for s, e in WECC])
    Ew = pd.concat([wecc, Ew], ignore_index=True)
    Es = find_episodes(Xs, sv_n, SVI_T[0], COV_SVI[1], rule)
    ref_nw, ref_sv = reference_entries()
    Aw = acceptance(Ew, ref_nw); As = acceptance(Es, ref_sv)
    tag = a.tag or "final"
    # frozen-rule MAJOR spans for the in_frozen_major column
    Fw = find_episodes(Xw, nw_n, COV_WA[0], COV_WA[1], RULE_FROZEN); Fw = Fw[Fw.major]
    Fs = find_episodes(Xs, sv_n, SVI_T[0], COV_SVI[1], RULE_FROZEN); Fs = Fs[Fs.major]
    for nm, E, Fz in (("nw", Ew, Fw), ("svi", Es, Fs)):
        E.insert(0, "region", nm.upper())
        E["episode"] = E.start.str[:7]
        E["in_frozen_major"] = [(E.kind.iloc[i] == "wecc") or bool(((Fz.start <= E.end.iloc[i]) & (Fz.end >= E.start.iloc[i])).any())
                                for i in range(len(E))]
        # unique episode id: start month, suffixed when two start in the same month
        E["episode"] = E.groupby("episode").cumcount().map(lambda k: "" if k == 0 else f"b{k}").radd(E.episode)
        E.to_csv(f"{OUT}/episodes_{nm}_{tag}.csv", index=False)
    Aw.to_csv(f"{OUT}/acceptance_nw_{tag}.csv", index=False); As.to_csv(f"{OUT}/acceptance_svi_{tag}.csv", index=False)
    log(f"[episodes] rule {rule} tag {tag}: NW {len(Ew)} episodes, {int(Ew.major.sum())} MAJOR; acceptance NW "
        f"{int(Aw.in_major.sum())}/{len(Aw)}; SVI {len(Es)} episodes, {int(Es.major.sum())} MAJOR; acceptance SVI "
        f"{int(As.in_major.sum())}/{len(As)}")
    print(Ew.to_string()); print(Aw.to_string()); print(Es.to_string()); print(As.to_string())


# ------------------------------------------------------------------------------------------------ counts
CUTS = {"050": 0.5, "075": 0.75, "085": 0.85, "090": 0.9}


def _count_wa_file(path):
    nw = set(wa_families().fid.astype(int))
    D = pd.read_parquet(path, columns=["date", "net", "series", "fid", "t", "P_own_WAVE_PNET", "tpl_insample"])
    D = D[(D.net == "core") & D.fid.isin(nw)]
    if not len(D):
        return None
    D = D.assign(half=W.half_of(D.t.to_numpy()), ins=D.tpl_insample.fillna(False).astype(bool))
    p = D.P_own_WAVE_PNET.to_numpy()
    for k, c in CUTS.items():
        D[f"n{k}"] = (p >= c).astype(np.int32)
    g = D.groupby(["fid", "date", "series", "half", "ins"]).agg(sumP=("P_own_WAVE_PNET", "sum"),
                                                                **{f"n{k}": (f"n{k}", "sum") for k in CUTS}).reset_index()
    return g


def count_wa(a):
    from multiprocessing import Pool
    fs = sorted(glob.glob(f"{WAP}/conf/P/prod/P_*.parquet"))
    with Pool(int(a.nproc or 12)) as pool:
        parts = [x for x in pool.imap(_count_wa_file, fs) if x is not None]
    C = pd.concat(parts, ignore_index=True)
    C.to_parquet(f"{OUT}/counts_wa_nw_fid_date_half.parquet", index=False)
    # check against rates_daily_v2 (own_pnet, NW): all halves summed
    S = C.groupby(["fid", "date", "series"])[[f"n{k}" for k in CUTS]].sum().unstack("series", fill_value=0)
    S.columns = [f"{a_}_{b_}" for a_, b_ in S.columns]
    R = pd.read_parquet(f"{WAP}/post_v2/rates_daily_v2.parquet", columns=["unit_kind", "label", "fid", "date", "region"] +
                        [f"n{k}_{s}" for k in CUTS for s in ("data", "rev", "c1")])
    R = R[(R.unit_kind == "family") & (R.label == "own_pnet") & (R.region == "NW")].set_index(["fid", "date"])
    J = R.join(S, rsuffix="_mine", how="outer").fillna(0)
    rep = {}
    for k in CUTS:
        for s in ("data", "rev", "c1"):
            col = f"n{k}_{s}"
            rep[col] = dict(rates_daily=float(J[col].sum()), mine=float(J.get(col + "_mine", J[col] * 0).sum()),
                            rows_differ=int((J[col] != J.get(col + "_mine", 0)).sum()))
    json.dump(rep, open(f"{OUT}/counts_wa_check_vs_rates_daily.json", "w"), indent=1)
    log("[count_wa]", len(C), "rows;", json.dumps(rep["n075_data"]), json.dumps(rep["n075_rev"]))



def svi_day_stations():
    out = {}
    for run in ("prod", "prod_S13"):
        for f in glob.glob(f"{SVI}/feat/{run}/*.json"):
            J = json.load(open(f))
            st = set(J.get("stations", []))
            out[(run, J["date"])] = (st, "S15" if {"LRIV", "SNB"} <= st else ("S14" if "SNB" in st else "S13"))
    return out


def count_svi(a):
    P = pd.read_parquet(f"{SVI}/conf/P_by_candidate.parquet", columns=["run", "date", "series", "fam", "t", "P_same", "sset_run"])
    p = P.P_same.to_numpy()
    for k, c in CUTS.items():
        P[f"n{k}"] = (p >= c).astype(np.int32)
    g = P.groupby(["run", "fam", "date", "series"]).agg(sumP=("P_same", "sum"), **{f"n{k}": (f"n{k}", "sum") for k in CUTS}).reset_index()
    g = g.rename(columns={"fam": "fid"})
    g.to_parquet(f"{OUT}/counts_svi_fid_date.parquet", index=False)
    DS = svi_day_stations()
    rows = [(r, d, ";".join(sorted(st)), ss, len(st)) for (r, d), (st, ss) in DS.items()]
    pd.DataFrame(rows, columns=["run", "date", "stations", "sset", "n_sta"]).to_csv(f"{OUT}/svi_day_stations.csv", index=False)
    DT = pd.read_csv(f"{SVI}/day_table.csv", dtype=str)
    ins = DT[(DT.date >= SVI_T[0]) & (DT.role != "prod")]
    rep = dict(rows=len(g), prod_days=int(g[g.run == "prod"].date.nunique()), s13_days=int(g[g.run == "prod_S13"].date.nunique()),
               s13_last=g[g.run == "prod_S13"].date.max(), day_table_nonprod_rows_from_2014=len(ins),
               day_table_roles_from_2014=DT[DT.date >= SVI_T[0]].role.value_counts().to_dict(),
               sset_prod_after_2021_04=pd.Series([v[1] for (r, d), v in DS.items() if r == "prod" and d >= "2021-04-01"]).value_counts().to_dict(),
               n075={f"{r}_{s}": int(x) for (r, s), x in g.groupby(["run", "series"]).n075.sum().items()})
    json.dump(rep, open(f"{OUT}/counts_svi_check.json", "w"), indent=1, default=str)
    log("[count_svi]", json.dumps(rep, default=str))



# ------------------------------------------------------------------------------------------------ per (family, episode) table
INS_WIN = {"V2008": ("2008-05-01", "2008-06-18"), "V2010": ("2010-07-14", "2010-09-18")}   # wa_prod_common.WINDOW_BUILD
SEALED = {"2007-01": "sealed", "2009-05": "semi"}
FIT_MONTHS = {"2007-10", "2008-05", "2009-02", "2009-04", "2010-03", "2010-07", "2010-08", "2010-09", "2011-08", "2011-12"}
EXCL_MONTHS = {"2011-04"}
WA_ERAS = {"E1": ("2007-01-01", "2011-12-31"), "E2": ("2012-01-01", "2014-12-31"), "E3": ("2015-01-01", "2026-09-15")}


def test_half_factor():
    return 86400.0 / (12 / 2 * (W.BLOCK_S - 2 * W.GUARD_S))


def wa_era(d):
    for k, (a0, b0) in WA_ERAS.items():
        if a0 <= d <= b0:
            return k
    return None


def _windows(E, X, fids, wecc_kind="wecc"):
    """-> dict (fid, episode) -> (window days list, exposure days list, X)."""
    Xi = X[X.n30 > 0]
    byf = {f: g.set_index("date").n30 for f, g in Xi.groupby("fid")}
    out = {}
    for e in E.itertuples():
        if e.kind == wecc_kind:
            days = [str(d.date()) for d in pd.date_range(e.start, e.end)]
            for f in fids:
                out[(f, e.episode)] = (days, days, np.nan)
            continue
        for f in fids:
            s = byf.get(f)
            if s is None:
                out[(f, e.episode)] = ([], [], 0.0); continue
            ex = s[(s.index >= e.start) & (s.index <= e.end)]
            win = set()
            for d in ex.index:
                t = pd.Timestamp(d)
                win.update(str((t + pd.Timedelta(days=k)).date()) for k in range(-2, 3))
            out[(f, e.episode)] = (sorted(win), sorted(ex.index), float(ex.sum()))
    return out


def _cube(C, fids, keycols):
    """C: fid, date + count columns -> dict col -> DataFrame (date x fid)."""
    out = {}
    for c in keycols:
        out[c] = C.pivot_table(index="date", columns="fid", values=c, aggfunc="sum", fill_value=0).reindex(columns=fids, fill_value=0)
    return out


def _sum(cube, days, f):
    if not days:
        return 0.0
    s = cube[f]
    return float(s.reindex(days).fillna(0).sum())


def wa_daily():
    """per (fid, date) WA counts with the in-sample rule applied: on V2008/V2010 days test half x F; plus raw half splits."""
    C = pd.read_parquet(f"{OUT}/counts_wa_nw_fid_date_half.parquet")
    F = test_half_factor()
    ins_day = np.zeros(len(C), bool)
    for a0, b0 in INS_WIN.values():
        ins_day |= (C.date >= a0).to_numpy() & (C.date <= b0).to_numpy()
    w = np.where(ins_day, np.where(C.half.to_numpy() == 1, F, 0.0), 1.0)
    cols = {}
    for k in list(CUTS) + ["sumP"]:
        c = f"n{k}" if k != "sumP" else "sumP"
        for sr in ("data", "rev", "c1"):
            m = (C.series == sr).to_numpy()
            cols[f"{c}_{sr}"] = np.where(m, C[c].to_numpy() * w, 0.0)
    m075 = C.n075.to_numpy()
    for h in (0, 1):
        for sr in ("data", "rev"):
            cols[f"h{h}_n075_{sr}"] = np.where((C.series == sr).to_numpy() & (C.half.to_numpy() == h), m075, 0.0)
    cols["ins_rows"] = C.ins.to_numpy().astype(float)
    cols["ins_n075_data"] = np.where(C.ins.to_numpy() & (C.series == "data").to_numpy(), m075, 0.0)
    D = pd.DataFrame(cols); D["fid"] = C.fid.to_numpy(); D["date"] = C.date.to_numpy()
    return D.groupby(["fid", "date"]).sum().reset_index()


def wa_core_days():
    """date -> set of loaded core stations (mf/prod/run_<date>.json)."""
    out = {}
    for p in glob.glob(f"{WAP}/mf/prod/run_*.json"):
        J = json.load(open(p))
        out[J["date"]] = set(J.get("members", {}).get("core", []))
    return out



def _famday_file(path):
    D = pd.read_parquet(path, columns=["net", "fid", "date", "n_sta", "n_sites", "observable", "stations", "nvalid_data"])
    return D[D.net == "core"].drop(columns=["net"])


def wa_famday(a=None):
    """per (fid, date) core-net scan record from mf/prod/famday_<date>.parquet (cached)."""
    out = f"{OUT}/wa_core_famday.parquet"
    if os.path.exists(out):
        return pd.read_parquet(out)
    from multiprocessing import Pool
    fs = sorted(glob.glob(f"{WAP}/mf/prod/famday_*.parquet"))
    with Pool(16) as pool:
        parts = list(pool.imap(_famday_file, fs, chunksize=20))
    D = pd.concat(parts, ignore_index=True)
    D = D[D.fid.isin(set(wa_families().fid))]
    D.to_parquet(out, index=False)
    log("[famday]", len(D), "rows", D.date.nunique(), "days")
    return D


def table(a):
    rows = []
    # ---------------- WA NW
    F = wa_families(); fids = F.fid.astype(int).tolist()
    E = pd.read_csv(f"{OUT}/episodes_nw_final.csv")
    X = pd.read_parquet(f"{OUT}/exposure_wa_nw.parquet")
    Xc = X.pivot_table(index="date", columns="fid", values="n30", fill_value=0)
    Dd = wa_daily()
    vals = [c for c in Dd.columns if c not in ("fid", "date")]
    cube = _cube(Dd, fids, vals)
    T = pd.read_csv(f"{WAP}/setup/tpl_library.csv")
    core = json.load(open(f"{WAP}/setup/core_sets.json"))
    tpl = T.groupby("fid").netsta.apply(set).to_dict()
    CD = wa_core_days()
    CH = pd.read_csv(f"{WAP}/post_v2/core_set_changes_v2.csv")
    CH = CH[CH.netsta != "*"]
    K = pd.read_csv(f"{WAP}/post_v2/rates_k_factors_v2.csv")
    K = K[(K.net == "E3") & (K.cut_tag == 75) & (K.label == "own_pnet") & (K.unit_kind == "family")]
    kmap = dict(zip(K.unit.astype(int), K.k))
    FD = wa_famday()
    FDi = {f: g.set_index("date") for f, g in FD.groupby("fid")}
    # principal carriers: smallest set of era-core stations holding >= 50 % of the family's template weight (w)
    principal = {}
    for era_, st in core.items():
        if era_ == "rule" or not isinstance(st, list):
            continue
        Tc = T[T.netsta.isin(st)].groupby(["fid", "netsta"]).w.max().reset_index()
        for f_, g in Tc.groupby("fid"):
            g = g.sort_values("w", ascending=False)
            cw = g.w.cumsum() / g.w.sum()
            principal[(int(f_), era_)] = set(g.netsta[: int(np.searchsorted(cw.to_numpy(), 0.5) + 1)])
    win = _windows(E, X, fids)
    for e in E.itertuples():
        era = wa_era(e.start) if wa_era(e.start) == wa_era(e.end) else "cross"
        months = {str(m) for m in pd.period_range(e.start, e.end, freq="M")}
        label = []
        for m in sorted(months):
            if m in SEALED: label.append(SEALED[m])
            if m in FIT_MONTHS: label.append(f"fit:{m}")
            if m in EXCL_MONTHS: label.append(f"excluded:{m}")
        ins_halves = any(not (e.end < a0 or e.start > b0) for a0, b0 in INS_WIN.values())
        for f in fids:
            days, exd, xs = win[(f, e.episode)]
            carriers = tpl.get(f, set()) & set(core.get(era, [])) if era in core else tpl.get(f, set())
            r = dict(region="NW", fid=f, episode=e.episode, ep_start=e.start, ep_end=e.end, kind=e.kind, major=bool(e.major),
                     in_frozen_major=bool(e.in_frozen_major), era=era, exposed=len(exd) > 0,
                     win_start=days[0] if days else "", win_end=days[-1] if days else "", win_days=len(days),
                     D=len(exd), X=xs)
            for c in vals:
                r[c] = _sum(cube[c], days, f)
            r["n_data"] = r["n075_data"]; r["ctl_main"] = r["n075_rev"]; r["ctl_scr"] = r["n075_c1"]
            r["N"] = r["n075_data"] - r["n075_rev"]; r["N_scr"] = r["n075_data"] - r["n075_c1"]
            r["N_expdays"] = _sum(cube["n075_data"], exd, f) - _sum(cube["n075_rev"], exd, f)
            for k in ("050", "085", "090"):
                r[f"N{k}"] = r[f"n{k}_data"] - r[f"n{k}_rev"]
            r["sumP_net"] = r["sumP_data"] - r["sumP_rev"]
            # active days, lags
            if days:
                s = cube["n075_data"][f].reindex(days).fillna(0)
                act = s[s > 0].index.tolist()
            else:
                act = []
            r["active_days"] = len(act)
            if act and exd and e.kind != "wecc":
                r["onset_lag_d"] = (pd.Timestamp(act[0]) - pd.Timestamp(exd[0])).days
                r["offset_lag_d"] = (pd.Timestamp(act[-1]) - pd.Timestamp(exd[-1])).days
            else:
                r["onset_lag_d"] = np.nan; r["offset_lag_d"] = np.nan
            r["n_carriers"] = len(carriers)
            r["window_days_no_run"] = int(sum(1 for d in days if d not in CD))
            fd = FDi.get(f)
            if fd is not None and days:
                fw = fd.reindex(days)
                r["live_carrier_days"] = float(fw.n_sta.fillna(0).sum())
                r["mean_sites"] = float(fw.n_sites.fillna(0).mean())
                r["obs_days"] = int(fw.observable.fillna(False).astype(bool).sum())
                ok = fd.reindex(exd).observable.fillna(False).astype(bool) if exd else pd.Series(dtype=bool)
                r["D_obs"] = int(ok.sum())
                r["X_obs"] = float(Xc[f].reindex(ok[ok].index).fillna(0).sum()) if (e.kind != "wecc" and f in Xc and ok.any()) else (np.nan if e.kind == "wecc" else 0.0)
            else:
                r["live_carrier_days"] = 0.0; r["mean_sites"] = 0.0; r["obs_days"] = 0; r["D_obs"] = 0
                r["X_obs"] = np.nan if e.kind == "wecc" else 0.0
            pc = principal.get((f, era), set())
            r["principal_carriers"] = ";".join(sorted(pc))
            if days:
                lo = str((pd.Timestamp(days[0]) - pd.Timedelta(days=45)).date()); hi = str((pd.Timestamp(days[-1]) + pd.Timedelta(days=45)).date())
                chw = CH[(CH.date >= lo) & (CH.date <= hi)]
                ch = chw[chw.netsta.isin(carriers)]
                chp = chw[chw.netsta.isin(pc) & ((chw.change == "leaves") | chw.change.str.startswith("node"))]
            else:
                ch = CH.iloc[:0]; chp = ch
            r["station_change"] = len(ch) > 0
            r["station_changes"] = ";".join(f"{x.date}:{x.netsta}:{x.change}" for x in ch.itertuples())[:2000]
            r["changed_carriers"] = ";".join(sorted(set(ch.netsta)))
            r["station_change_principal"] = len(chp) > 0
            r["changed_principal"] = ";".join(sorted(set(chp.netsta)))
            r["insample"] = "insample_halves_x%.3f" % test_half_factor() if ins_halves else ""
            r["t2012_templates"] = bool(r["ins_rows"] > 0)
            lab = list(label)
            if ins_halves: lab.append("insample_halves")
            if r["t2012_templates"]: lab.append("t2012 templates")
            r["sample_label"] = ";".join(lab)
            r["excluded"] = any(x.startswith("excluded") for x in label)
            r["k_E3"] = kmap.get(f, np.nan) if era == "E3" else np.nan
            r["kN_E3"] = r["k_E3"] * r["N"] if era == "E3" else np.nan
            rows.append(r)
    TW = pd.DataFrame(rows)
    TW.to_csv(f"{OUT}/fam_episode_nw.csv", index=False)
    log("[table] NW", len(TW), "rows")
    # ---------------- SVI
    S = svi_families(); sfids = S.fid.astype(int).tolist()
    E = pd.read_csv(f"{OUT}/episodes_svi_final.csv")
    X = pd.read_parquet(f"{OUT}/exposure_svi.parquet")
    C = pd.read_parquet(f"{OUT}/counts_svi_fid_date.parquet")
    DSt = pd.read_csv(f"{OUT}/svi_day_stations.csv")
    stmap = {(r.run, r.date): set(str(r.stations).split(";")) for r in DSt.itertuples()}
    ssmap = {(r.run, r.date): r.sset for r in DSt.itertuples()}
    q = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_template_qc_2012.csv")
    stpl = q[q.usable & q.template].groupby("bkey").sta.apply(set).to_dict()
    sys.path.insert(0, f"{ROOT}/lfe_extension/scripts")
    CHANGES = [("2014-06-17", "LZB", "LZB ends"), ("2014-09-03", "YOUB", "YOUB ends"), ("2015-01-12", "PGC", "PGC response sub-change"),
               ("2017-02-23", "NLLB", "NLLB BH->HH"), ("2017-03-03", "VGZ", "VGZ sensor swap"), ("2017-06-08", "SNB", "SNB sensor swap"),
               ("2017-08-10", "PGC", "PGC BH->HH"), ("2018-12-20", "PFB", "PFB ends"), ("2019-05-30", "LRIV", "LRIV BH->HH, dropped"),
               ("2021-03-28", "SNB", "SNB ends")]   # svi_ref_common.CHANGES with the station named
    # effective-set transitions of the production run
    pr = DSt[DSt.run == "prod"].sort_values("date")
    tr = pr[pr.sset != pr.sset.shift()].iloc[1:]
    SETCH = [(r.date, "SET", f"set -> {r.sset}") for r in tr.itertuples()]
    win = _windows(E, X, sfids)
    Xs_c = X.pivot_table(index="date", columns="fid", values="n30", fill_value=0)
    for series in ("S13fixed", "eff"):
        if series == "S13fixed":
            Cs = pd.concat([C[(C.run == "prod_S13") & (C.date < "2021-04-01")], C[(C.run == "prod") & (C.date >= "2021-04-01")]])
            runof = lambda d: "prod_S13" if d < "2021-04-01" else "prod"
        else:
            Cs = C[C.run == "prod"]
            runof = lambda d: "prod"
        Wd = Cs.pivot_table(index=["fid", "date"], columns="series", values=["n050", "n075", "n085", "n090", "sumP"], aggfunc="sum", fill_value=0)
        Wd.columns = [f"{a_}_{b_}" for a_, b_ in Wd.columns]
        Wd = Wd.reset_index()
        for c in [f"{k}_{sr}" for k in ("n050", "n075", "n085", "n090", "sumP") for sr in ("data", "c1")]:
            if c not in Wd: Wd[c] = 0.0
        vals = [c for c in Wd.columns if c not in ("fid", "date")]
        cube = _cube(Wd, sfids, vals)
        rows = []
        for e in E.itertuples():
            for f in sfids:
                days, exd, xs = win[(f, e.episode)]
                carriers_all = stpl.get(f, set())
                sets = sorted({ssmap.get((runof(d), d), "") for d in days} - {""})
                r = dict(region="SVI", series=series, fid=f, episode=e.episode, ep_start=e.start, ep_end=e.end, kind=e.kind,
                         major=bool(e.major), in_frozen_major=bool(e.in_frozen_major),
                         era="S13fixed" if series == "S13fixed" else ("+".join(sets) if sets else ""), exposed=len(exd) > 0,
                         win_start=days[0] if days else "", win_end=days[-1] if days else "", win_days=len(days), D=len(exd), X=xs)
                for c in vals:
                    r[c] = _sum(cube[c], days, f)
                r["n_data"] = r["n075_data"]; r["ctl_main"] = r["n075_c1"]; r["ctl_scr"] = r["n075_c1"]
                r["N"] = r["n075_data"] - r["n075_c1"]; r["N_scr"] = r["N"]
                r["N_expdays"] = _sum(cube["n075_data"], exd, f) - _sum(cube["n075_c1"], exd, f)
                for k in ("050", "085", "090"):
                    r[f"N{k}"] = r[f"n{k}_data"] - r[f"n{k}_c1"]
                r["sumP_net"] = r["sumP_data"] - r["sumP_c1"]
                if days:
                    s_ = cube["n075_data"][f].reindex(days).fillna(0)
                    act = s_[s_ > 0].index.tolist()
                else:
                    act = []
                r["active_days"] = len(act)
                if act and exd:
                    r["onset_lag_d"] = (pd.Timestamp(act[0]) - pd.Timestamp(exd[0])).days
                    r["offset_lag_d"] = (pd.Timestamp(act[-1]) - pd.Timestamp(exd[-1])).days
                else:
                    r["onset_lag_d"] = np.nan; r["offset_lag_d"] = np.nan
                lc = 0; nostat = 0
                for d in days:
                    st = stmap.get((runof(d), d))
                    if st is None:
                        nostat += 1; continue
                    lc += len(carriers_all & st)
                r["n_carriers"] = len(carriers_all & (set(["PGC", "VGZ", "NLLB", "B001", "B003", "B004", "B005", "B006", "B007", "B009", "B010", "B011", "B926"]) if series == "S13fixed" else carriers_all))
                r["live_carrier_days"] = lc; r["window_days_no_run"] = nostat
                okd = [d for d in exd if len(carriers_all & stmap.get((runof(d), d), set())) >= 2]
                r["obs_days"] = int(sum(1 for d in days if len(carriers_all & stmap.get((runof(d), d), set())) >= 2))
                r["D_obs"] = len(okd)
                r["X_obs"] = float(Xs_c[f].reindex(okd).fillna(0).sum()) if (f in Xs_c and okd) else 0.0
                r["mean_sites"] = lc / max(len(days), 1)
                chs = []
                if days:
                    lo = str((pd.Timestamp(days[0]) - pd.Timedelta(days=45)).date()); hi = str((pd.Timestamp(days[-1]) + pd.Timedelta(days=45)).date())
                    use_set = SETCH if series == "eff" else []
                    for dd, sta, why in CHANGES + use_set:
                        if lo <= dd <= hi and (sta == "SET" or sta in carriers_all):
                            if series == "S13fixed" and sta in ("SNB", "LRIV"):
                                continue
                            chs.append((dd, sta, why))
                r["station_change"] = len(chs) > 0
                r["station_changes"] = ";".join(f"{a_}:{b_}:{c_}" for a_, b_, c_ in chs)
                r["changed_carriers"] = ";".join(sorted({b_ for _, b_, _ in chs}))
                r["station_change_principal"] = r["station_change"]; r["changed_principal"] = r["changed_carriers"]
                r["principal_carriers"] = ""
                r["insample"] = ""; r["t2012_templates"] = False; r["sample_label"] = ""; r["excluded"] = False
                r["k_E3"] = np.nan; r["kN_E3"] = np.nan
                rows.append(r)
        TS = pd.DataFrame(rows)
        TS.to_csv(f"{OUT}/fam_episode_svi_{series}.csv", index=False)
        log("[table] SVI", series, len(TS), "rows")



# ------------------------------------------------------------------------------------------------ sanity checks
MERLIN_REF = {"2011-08": 0.72, "2012-09": 0.97, "2013-09": 0.74, "2014-12": 0.52, "2016-01": 0.83, "2017-03": 0.81, "2018-06": 0.55,
              "2020-10": 0.80, "2021-10": 0.59, "2023-04": 0.71, "2024-08": 0.40, "2025-09": 0.69}


def load_tables():
    TW = pd.read_csv(f"{OUT}/fam_episode_nw.csv", low_memory=False)
    TS = pd.read_csv(f"{OUT}/fam_episode_svi_S13fixed.csv", low_memory=False)
    TE = pd.read_csv(f"{OUT}/fam_episode_svi_eff.csv", low_memory=False)
    for T in (TW, TS, TE):
        for c in ("sample_label", "insample", "station_changes", "changed_carriers", "principal_carriers", "changed_principal", "win_start", "win_end"):
            if c in T:
                T[c] = T[c].fillna("").astype(str)
    return TW, TS, TE


def regional(T):
    g = T.groupby("episode").agg(ep_start=("ep_start", "first"), ep_end=("ep_end", "first"), kind=("kind", "first"),
                                 major=("major", "first"), in_frozen_major=("in_frozen_major", "first"), era=("era", "first"),
                                 sample_label=("sample_label", lambda x: ";".join(sorted({y for v in x for y in str(v).split(";") if y and y != "t2012 templates"}))),
                                 n_t2012=("t2012_templates", "sum"), families_exposed=("exposed", "sum"),
                                 N=("N", "sum"), N_expdays=("N_expdays", "sum"), N_scr=("N_scr", "sum"), data=("n_data", "sum"),
                                 ctl_main=("ctl_main", "sum"), ctl_scr=("ctl_scr", "sum"), X=("X", "sum"), D=("D", "sum"),
                                 D_obs=("D_obs", "sum"), N050=("N050", "sum"), N085=("N085", "sum"), N090=("N090", "sum")).reset_index()
    g["ctl_main_share"] = g.ctl_main / g.data.where(g.data > 0)
    g["ctl_scr_share"] = g.ctl_scr / g.data.where(g.data > 0)
    g["N_per_D"] = g.N / g.D.where(g.D > 0)
    g["Nexp_per_D"] = g.N_expdays / g.D.where(g.D > 0)
    g["N_per_X"] = g.N / g.X.where(g.X > 0)
    return g.sort_values("ep_start")


def checks(a):
    from scipy.stats import spearmanr, beta as Beta
    TW, TS, TE = load_tables()
    res = {}
    # (a) WA build / test
    M = TW[TW.major & (TW.kind != "wecc")]
    rows = []
    for ep, g in TW[TW.major].groupby("episode"):
        lab = ";".join(sorted({y for v in g.sample_label for y in v.split(";") if y}))
        h0 = g.h0_n075_data.sum() - g.h0_n075_rev.sum(); h1 = g.h1_n075_data.sum() - g.h1_n075_rev.sum()
        oos = not any(k in lab for k in ("insample_halves", "sealed", "semi", "excluded"))
        n = h0 + h1
        if n > 0:
            lo_, hi_ = Beta.ppf(0.025, max(h0, 0) + 1e-9, h1 + 1), Beta.ppf(0.975, h0 + 1, max(h1, 0) + 1e-9)
            rlo, rhi = lo_ / (1 - lo_), hi_ / (1 - hi_)
        else:
            rlo = rhi = np.nan
        rows.append(dict(episode=ep, era=g.era.iloc[0], label=lab, out_of_sample=oos, build_net=h0, test_net=h1,
                         ratio=h0 / h1 if h1 > 0 else np.nan, ratio95_lo=rlo, ratio95_hi=rhi, n_half=n,
                         all_halves_data=g.n075_data.sum(), test_share_data=(g.h1_n075_data.sum() / g.n075_data.sum()) if g.n075_data.sum() > 0 else np.nan))
    A = pd.DataFrame(rows)
    A.to_csv(f"{OUT}/check_a_build_test.csv", index=False)
    o = A[A.out_of_sample]
    pooled = o.build_net.sum() / o.test_net.sum()
    big = o[o.n_half >= 400]
    per_ok = ((big.ratio95_hi >= 0.85) & (big.ratio95_lo <= 1.15))
    share = o.test_n = None
    ts = TW[TW.major & ~TW.sample_label.str.contains("insample_halves")]
    test_share = float(ts.h1_n075_data.sum() / ts.n075_data.sum())
    res["a"] = dict(pooled_ratio_out_of_sample=pooled, n_episodes_out_of_sample=len(o), n_big=len(big), big_all_overlap=bool(per_ok.all()),
                    big_failing=big[~per_ok].episode.tolist(), pass_=bool(0.85 <= pooled <= 1.15 and per_ok.all()),
                    out_of_sample_ratio_range=[float(o[o.n_half >= 100].ratio.min()), float(o[o.n_half >= 100].ratio.max())],
                    insample=A[A.label.str.contains("insample_halves")][["episode", "ratio", "n_half"]].to_dict("records"),
                    test_half_share_out_of_sample=test_share, empirical_factor=1 / test_share, design_factor=test_half_factor())
    # (b) control share per MAJOR episode
    RW = regional(TW); RS = regional(TS); RE = regional(TE)
    for nm, R in (("nw", RW), ("svi_S13fixed", RS), ("svi_eff", RE)):
        R.to_csv(f"{OUT}/regional_{nm}.csv", index=False)
    b = {}
    for nm, R in (("nw", RW), ("svi_S13fixed", RS), ("svi_eff", RE)):
        m = R[R.major & (R.data > 0)]
        fail = m[m.ctl_main_share > 0.02]
        b[nm] = dict(n_major=int(R.major.sum()), max_ctl_main_share=float(m.ctl_main_share.max()),
                     max_ctl_scr_share=float(m.ctl_scr_share.max()), failing=fail[["episode", "data", "ctl_main", "ctl_main_share"]].to_dict("records"),
                     pass_=len(fail) == 0, pooled_ctl_main_share=float(m.ctl_main.sum() / m.data.sum()))
    res["b"] = b
    # (c) Spearman within era
    c = {}
    def sp(R, label):
        m = R[R.major & R.X.notna() & (R.X > 0)]
        if len(m) < 4:
            return dict(n=len(m), rho=np.nan, p=np.nan, pass_=None)
        r_, p_ = spearmanr(m.N, m.X)
        r2, _ = spearmanr(m.N, m.D)
        r3, _ = spearmanr(m.N_expdays, m.X)
        return dict(n=len(m), rho=float(r_), p=float(p_), pass_=bool(r_ >= 0.5), rho_N_vs_D=float(r2), rho_Nexpdays_vs_X=float(r3),
                    frozen_subset=(lambda mm: dict(n=len(mm), rho=float(spearmanr(mm.N, mm.X)[0]) if len(mm) >= 4 else np.nan))(m[m.in_frozen_major]))
    for era in ("E1", "E2", "E3"):
        c[f"nw_{era}"] = sp(RW[RW.era == era], era)
    c["svi_S13fixed"] = sp(RS, "S13")
    for ss in ("S15", "S14", "S13"):
        c[f"svi_eff_{ss}"] = sp(RE[RE.era == ss], ss)
    c["svi_eff_all"] = sp(RE, "all")
    res["c"] = c
    # named cases
    def case(R, ep, era_mask):
        m = R[R.major & era_mask]
        r = R[R.episode == ep].iloc[0]
        return dict(episode=ep, start=r.ep_start, end=r.ep_end, N=float(r.N), X=float(r.X), D=float(r.D), N_per_D=float(r.N_per_D),
                    N_per_X=float(r.N_per_X), era_median_N_per_D=float(m.N_per_D.median()), era_median_N_per_X=float(m.N_per_X.median()),
                    rank_N_per_X=f"{int((m.N_per_X < r.N_per_X).sum()) + 1} of {len(m)} (1 = lowest)")
    res["named"] = dict(svi_2017_03=case(RS, "2017-03", RS.era == RS.era),
                        svi_2025_10=case(RS, "2025-09b1", RS.era == RS.era),
                        svi_eff_2017_03=case(RE, "2017-03", RE.era == RE.era),
                        svi_eff_2025_10=case(RE, "2025-09b1", RE.era == RE.era),
                        wa_2024_08=case(RW, "2024-08", RW.era == "E3"))
    # Merlin reference comparison: map each reference month to the MAJOR episode containing its 15th
    mr = []
    for mth, v in MERLIN_REF.items():
        d = f"{mth}-15"
        h = RW[RW.major & (RW.ep_start <= d) & (RW.ep_end >= f"{mth}-01")]
        if len(h):
            r = h.iloc[0]
            mr.append(dict(ref_month=mth, merlin=v, episode=r.episode, N_per_D=r.N_per_D, Nexp_per_D=r.Nexp_per_D, ratio=r.N_per_D / v))
    MR = pd.DataFrame(mr); MR.to_csv(f"{OUT}/merlin_reference_compare.csv", index=False)
    res["merlin_ref"] = dict(median_ratio=float(MR.ratio.median()), rows=MR.to_dict("records"))
    json.dump(res, open(f"{OUT}/sanity_checks.json", "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    log("[checks] a pass", res["a"]["pass_"], "pooled", round(pooled, 3), "| b", {k: v["pass_"] for k, v in b.items()},
        "| c", {k: (round(v["rho"], 3) if v["rho"] == v["rho"] else None, v["n"]) for k, v in c.items()})



def check_a_boot(a):
    """diagnostic (not pre-registered): day-bootstrap 95 % interval of the build/test ratio per out-of-sample MAJOR episode
    (episode span +- 2 d, all NW families), allowing for the day-scale clustering of LFE bursts that a Poisson interval ignores."""
    C = pd.read_parquet(f"{OUT}/counts_wa_nw_fid_date_half.parquet", columns=["date", "series", "half", "n075"])
    C = C[C.half >= 0]
    sgn = np.where(C.series == "data", 1, np.where(C.series == "rev", -1, 0))
    C = C.assign(v=C.n075 * sgn).groupby(["date", "half"]).v.sum().unstack(fill_value=0)
    A = pd.read_csv(f"{OUT}/check_a_build_test.csv")
    E = pd.read_csv(f"{OUT}/episodes_nw_final.csv").set_index("episode")
    rng = np.random.default_rng(5)
    out = []; allb = []
    for r in A[A.out_of_sample].itertuples():
        e = E.loc[r.episode]
        days = [str(d.date()) for d in pd.date_range(pd.Timestamp(e.start) - pd.Timedelta(days=2), pd.Timestamp(e.end) + pd.Timedelta(days=2))]
        M = C.reindex(days).fillna(0)
        h0, h1 = M[0].to_numpy(), M[1].to_numpy()
        allb.append((h0, h1))
        bs = []
        for _ in range(2000):
            ii = rng.integers(0, len(days), len(days))
            t1 = h1[ii].sum()
            bs.append(h0[ii].sum() / t1 if t1 > 0 else np.nan)
        bs = np.array(bs, float)
        out.append(dict(episode=r.episode, ratio_span=h0.sum() / h1.sum() if h1.sum() > 0 else np.nan, n_half_span=h0.sum() + h1.sum(),
                        boot95_lo=np.nanpercentile(bs, 2.5), boot95_hi=np.nanpercentile(bs, 97.5)))
    B = pd.DataFrame(out)
    B["overlaps_band"] = (B.boot95_hi >= 0.85) & (B.boot95_lo <= 1.15)
    # pooled, bootstrapping episodes (cluster = episode)
    h0s = np.array([x[0].sum() for x in allb]); h1s = np.array([x[1].sum() for x in allb])
    pb = []
    for _ in range(2000):
        ii = rng.integers(0, len(allb), len(allb)); pb.append(h0s[ii].sum() / h1s[ii].sum())
    B.to_csv(f"{OUT}/check_a_build_test_dayboot.csv", index=False)
    J = dict(pooled=float(h0s.sum() / h1s.sum()), pooled_episode_boot95=[float(np.percentile(pb, 2.5)), float(np.percentile(pb, 97.5))],
             big_n400=int((B.n_half_span >= 400).sum()), big_overlap=int(B[B.n_half_span >= 400].overlaps_band.sum()))
    json.dump(J, open(f"{OUT}/check_a_dayboot.json", "w"), indent=1)
    log("[check_a_boot]", json.dumps(J))
    print(B[B.n_half_span >= 400].to_string())



# ------------------------------------------------------------------------------------------------ step 5 consistency
def bh(p):
    p = np.asarray(p, float); q = np.full(len(p), np.nan); ok = np.isfinite(p)
    if ok.sum():
        pv = p[ok]; o = np.argsort(pv); n = len(pv)
        r = pv[o] * n / np.arange(1, n + 1)
        r = np.minimum.accumulate(r[::-1])[::-1]
        qq = np.empty(n); qq[o] = np.minimum(r, 1.0); q[ok] = qq
    return q


def share_tests(T, fams, block, ncol="N", nboot=2000, seed=0):
    """T: per (fid, episode) rows of one block (MAJOR episodes of one era). Returns (family table, per-row residual table)."""
    rng = np.random.default_rng(seed)
    eps = sorted(T.episode.unique())
    fids = fams.fid.astype(int).tolist()
    pos = fams.set_index("fid")[["lat", "lon"]]
    Nm = T.pivot_table(index="fid", columns="episode", values=ncol, aggfunc="sum").reindex(index=fids, columns=eps).fillna(0)
    Nm = np.maximum(np.round(Nm.to_numpy()), 0)
    Dm = T.pivot_table(index="fid", columns="episode", values="D", aggfunc="sum").reindex(index=fids, columns=eps).fillna(0).to_numpy()
    Xm = T.pivot_table(index="fid", columns="episode", values="X", aggfunc="sum").reindex(index=fids, columns=eps).to_numpy()
    if "test_ok" in T:
        OKm = T.assign(ok=T.test_ok.astype(float)).pivot_table(index="fid", columns="episode", values="ok", aggfunc="min").reindex(
            index=fids, columns=eps).fillna(1.0).to_numpy() > 0.5
    else:
        OKm = np.ones((len(fids), len(eps)), bool)
    lat = pos.reindex(fids).lat.to_numpy(); lon = pos.reindex(fids).lon.to_numpy()
    dist = np.array([W.km(lat, lon, la, lo) for la, lo in zip(lat, lon)])
    frow, rrow = [], []
    for i, f in enumerate(fids):
        nb = np.flatnonzero(dist[i] <= 40.0)
        nb_kind = "40km"
        if len(nb) - 1 < 5:
            nb = np.arange(len(fids)); nb_kind = "region"
        Tt = Nm[nb].sum(0)
        use = (Dm[i] >= 1) & (Tt > 0) & OKm[i]
        n_e = int(use.sum()); sN = float(Nm[i, use].sum())
        d = dict(block=block, fid=f, n_neighbours=len(nb) - 1, neighbourhood=nb_kind, n_episodes_tested=n_e, sumN=sN,
                 sumT=float(Tt[use].sum()), n_episodes_exposed=int((Dm[i] >= 1).sum()))
        if n_e < 3 or sN < 20:
            d.update(tested=False, reason="sumN<20" if sN < 20 else "n_e<3", chi2=np.nan, df=np.nan, chi2_df=np.nan, p_boot=np.nan,
                     disp=np.nan, disp_df=np.nan, p_share=np.nan)
            frow.append(d)
            for j in np.flatnonzero(Dm[i] >= 1):
                rrow.append(dict(block=block, fid=f, episode=eps[j], N=Nm[i, j], T=Tt[j], E=np.nan, resid=np.nan, tested=False))
            continue
        n, t = Nm[i, use], Tt[use]
        pf = n.sum() / t.sum()
        E = pf * t
        v = t * pf * (1 - pf)
        chi2 = float(((n - E) ** 2 / np.where(v > 0, v, np.nan)).sum())
        # parametric bootstrap
        bs = np.empty(nboot)
        for b in range(nboot):
            nb_ = rng.binomial(t.astype(np.int64), pf)
            pb = nb_.sum() / t.sum()
            vb = t * pb * (1 - pb)
            bs[b] = ((nb_ - pb * t) ** 2 / np.where(vb > 0, vb, np.nan)).sum()
        pboot = float((1 + (bs >= chi2).sum()) / (nboot + 1))
        # secondary Poisson dispersion vs X
        x = Xm[i, use]
        ok = np.isfinite(x) & (x > 0)
        if ok.sum() >= 3 and n[ok].sum() > 0:
            rho = n[ok].sum() / x[ok].sum(); ex = rho * x[ok]
            disp = float(((n[ok] - ex) ** 2 / ex).sum() / (ok.sum() - 1))
        else:
            disp = np.nan
        phi = float((((n - E) ** 2 - v).sum()) / (E ** 2).sum()) if (E ** 2).sum() > 0 else np.nan
        d.update(tested=True, reason="", p_share=pf, chi2=chi2, df=n_e - 1, chi2_df=chi2 / (n_e - 1), p_boot=pboot, disp=disp,
                 disp_df=int(ok.sum()) - 1, cv_extra=float(np.sqrt(max(phi, 0.0))) if np.isfinite(phi) else np.nan,
                 mean_E=float(E.mean()))
        frow.append(d)
        for j, k in zip(np.flatnonzero(use), range(n_e)):
            rrow.append(dict(block=block, fid=f, episode=eps[j], N=n[k], T=t[k], E=E[k],
                             resid=(n[k] - E[k]) / np.sqrt(v[k]) if v[k] > 0 else np.nan, tested=True))
    FT = pd.DataFrame(frow)
    FT["q_bh"] = bh(FT.p_boot)
    FT["label"] = np.where(~FT.tested, "not tested", np.where(FT.q_bh < 0.05, "over-dispersed", "consistent"))
    FT["exposure_explained"] = FT.tested & (FT.disp < 2) & (FT.chi2_df > 3)
    return FT, pd.DataFrame(rrow), dist


def weighted_jaccard_wa(fids, era):
    T = pd.read_csv(f"{WAP}/setup/tpl_library.csv")
    core = json.load(open(f"{WAP}/setup/core_sets.json"))[era]
    Wt = T[T.netsta.isin(core)].groupby(["fid", "netsta"]).w.max().unstack(fill_value=0.0).reindex(index=fids, columns=core).fillna(0).to_numpy()
    n = len(fids); J = np.ones((n, n))
    for i in range(n):
        mn = np.minimum(Wt[i], Wt).sum(1); mx = np.maximum(Wt[i], Wt).sum(1)
        J[i] = np.where(mx > 0, mn / mx, np.nan)
    return J


def set_jaccard_svi(fids, stations=None):
    q = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_template_qc_2012.csv")
    st = q[q.usable & q.template].groupby("bkey").sta.apply(set).to_dict()
    S = [st.get(f, set()) & (set(stations) if stations else st.get(f, set())) for f in fids]
    n = len(fids); J = np.ones((n, n))
    for i in range(n):
        for j in range(n):
            u = len(S[i] | S[j]); J[i, j] = len(S[i] & S[j]) / u if u else np.nan
    return J


def mantel(RT, fids, dist, J, nperm=1000, seed=1, min_shared=5):
    from scipy.stats import spearmanr
    R = RT[RT.tested].pivot_table(index="fid", columns="episode", values="resid").reindex(index=fids)
    A = R.to_numpy()
    n = len(fids)
    C = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i + 1, n):
            m = np.isfinite(A[i]) & np.isfinite(A[j])
            if m.sum() >= min_shared and np.std(A[i, m]) > 0 and np.std(A[j, m]) > 0:
                C[i, j] = np.corrcoef(A[i, m], A[j, m])[0, 1]
    iu = np.argwhere(np.isfinite(C))
    if len(iu) < 10:
        return dict(n_pairs=len(iu)), None
    y = C[iu[:, 0], iu[:, 1]]
    def fit(perm):
        ld = np.log10(np.maximum(dist[perm][:, perm][iu[:, 0], iu[:, 1]], 0.5))
        jj = J[perm][:, perm][iu[:, 0], iu[:, 1]]
        ok = np.isfinite(jj)
        Xd = np.column_stack([np.ones(ok.sum()), ld[ok], jj[ok]])
        b = np.linalg.lstsq(Xd, y[ok], rcond=None)[0]
        return b[1], b[2], spearmanr(ld[ok], y[ok])[0], spearmanr(jj[ok], y[ok])[0] if np.std(jj[ok]) > 0 else np.nan
    obs = fit(np.arange(n))
    rng = np.random.default_rng(seed)
    P = np.array([fit(rng.permutation(n)) for _ in range(nperm)])
    pv = [float((1 + (np.abs(P[:, k]) >= abs(obs[k])).sum()) / (nperm + 1)) if np.isfinite(obs[k]) else np.nan for k in range(4)]
    ld = np.log10(np.maximum(dist[iu[:, 0], iu[:, 1]], 0.5)); jj = J[iu[:, 0], iu[:, 1]]
    pairs = pd.DataFrame(dict(fid_i=np.array(fids)[iu[:, 0]], fid_j=np.array(fids)[iu[:, 1]], dist_km=dist[iu[:, 0], iu[:, 1]],
                              jaccard=jj, corr=y))
    # binned summary
    bins = [0, 10, 20, 40, 80, 160, 400]
    pairs["dist_bin"] = pd.cut(pairs.dist_km, bins)
    return dict(n_pairs=len(y), mean_corr=float(y.mean()), coef_logdist=float(obs[0]), p_logdist=pv[0], coef_jaccard=float(obs[1]),
                p_jaccard=pv[1], spearman_logdist=float(obs[2]), p_spearman_logdist=pv[2], spearman_jaccard=float(obs[3]) if np.isfinite(obs[3]) else None,
                p_spearman_jaccard=pv[3], jaccard_range=[float(np.nanmin(jj)), float(np.nanmax(jj))],
                corr_by_dist=pairs.groupby("dist_bin", observed=True)["corr"].agg(["mean", "size"]).reset_index().astype(str).to_dict("records")), pairs


def label_rows(RT, T, fams, dist, block, carrier_col="changed_principal", flag_col="station_change_principal", carriers=None):
    """per (family, episode) labels for |resid| >= 2 (pre-registered)."""
    fids = fams.fid.astype(int).tolist(); ix = {f: i for i, f in enumerate(fids)}
    meta = T.set_index(["fid", "episode"])
    R = RT[RT.tested].copy()
    rr = R.set_index(["fid", "episode"]).resid
    labs, why = [], []
    for r in R.itertuples():
        if not (np.isfinite(r.resid) and abs(r.resid) >= 2):
            labs.append("within noise"); why.append(""); continue
        i = ix[r.fid]
        same_ep = R[(R.episode == r.episode) & (R.fid != r.fid)]
        # spatial: others within 20 km
        near = same_ep[[dist[i, ix[g]] <= 20 for g in same_ep.fid]]
        sp = len(near) >= 3 and np.sign(near.resid.mean()) == np.sign(r.resid) and abs(near.resid.mean()) >= 1
        # station-coincident
        m = meta.loc[(r.fid, r.episode)]
        st = False; w = ""
        if bool(m[flag_col]):
            ch = set(x for x in str(m[carrier_col]).split(";") if x)
            if carriers is not None and ch:
                sh = same_ep[[len(ch & carriers.get(g, set())) > 0 for g in same_ep.fid]]
                if len(sh) >= 1 and np.sign(sh.resid.mean()) == np.sign(r.resid) and abs(sh.resid.mean()) >= 1:
                    st = True; w = f"{len(sh)} sharing {';'.join(sorted(ch))} mean {sh.resid.mean():.2f}"
        labs.append("spatially coherent + station-coincident" if (sp and st) else "spatially coherent" if sp else
                    "station-coincident" if st else "isolated")
        why.append(w + (f" | near20 n={len(near)} mean {near.resid.mean():.2f}" if len(near) else ""))
    R["label"] = labs; R["label_note"] = why
    return R


def consistency(a):
    TW, TS, TE = load_tables()
    FW = wa_families(); FS = svi_families()
    summ = {}
    allF, allR, allP = [], [], []
    # carriers for station-coincident labels
    Tl = pd.read_csv(f"{WAP}/setup/tpl_library.csv"); core = json.load(open(f"{WAP}/setup/core_sets.json"))
    def principal_of(era):
        out = {}
        for r in TW[TW.era == era].drop_duplicates("fid").itertuples():
            out[int(r.fid)] = set(x for x in str(r.principal_carriers).split(";") if x)
        return out
    q = pd.read_csv(f"{ROOT}/lfe_catalog_repair/results/e5r_template_qc_2012.csv")
    svi_car = q[q.usable & q.template].groupby("bkey").sta.apply(set).to_dict()
    S13 = ["PGC", "VGZ", "NLLB", "B001", "B003", "B004", "B005", "B006", "B007", "B009", "B010", "B011", "B926"]
    jobs = []
    for era in ("E3", "E2", "E1"):
        sub = TW[TW.major & (TW.era == era) & ~TW.excluded]
        jobs.append(("NW", era, sub, FW, "wa", era))
    jobs.append(("SVI_S13fixed", "S13fixed", TS[TS.major], FS, "svi13", None))
    for ss in ("S15", "S14", "S13"):
        jobs.append(("SVI_eff", ss, TE[TE.major & (TE.era == ss)], FS, "svieff", ss))
    for reg, blk, sub, fams, jk, era in jobs:
        if not len(sub) or sub.episode.nunique() < 3:
            summ[f"{reg}_{blk}"] = dict(n_episodes=int(sub.episode.nunique()), note="< 3 episodes, not tested"); continue
        VARS = [("main", sub, "N"), ("frozen_major", sub[sub.in_frozen_major], "N"), ("N_expdays", sub, "N_expdays"),
                ("ctl_main", sub, "ctl_main"), ("ctl_scr", sub, "ctl_scr")]
        if os.environ.get("EPR_ROBUST") == "1":
            VARS += [("cut_050", sub, "N050"), ("cut_085", sub, "N085"), ("flag_excl_075", None, "N")]
        main_labels = None
        for variant, subv, ncol in VARS:
            if variant == "flag_excl_075":
                subv = sub[~sub.episode.isin(["2017-03", "2025-09b1"])].copy() if reg.startswith("SVI") else sub.copy()
                med = subv[subv.win_days > 0].groupby("fid").mean_sites.median()
                low = subv.mean_sites.to_numpy() < 0.5 * subv.fid.map(med).fillna(0).to_numpy()
                sc = np.zeros(len(subv), bool)
                if main_labels is not None:
                    bad = set(map(tuple, main_labels[main_labels.label.str.contains("station-coincident") |
                                                     main_labels.get("label_broadflag", pd.Series("", index=main_labels.index)).fillna("").str.contains("station-coincident")]
                                  [["fid", "episode"]].to_numpy()))
                    sc = np.array([(f_, e_) in bad for f_, e_ in zip(subv.fid, subv.episode)])
                subv["test_ok"] = ~(low & (subv.D > 0)) & ~sc
                summ[f"{reg}_{blk}_flag_excl_075_removed"] = dict(episodes_dropped=["2017-03", "2025-09b1"] if reg.startswith("SVI") else [],
                                                                  cells_low_carriers=int((low & (subv.D > 0)).sum()),
                                                                  cells_station_coincident=int(sc.sum()))
            if subv.episode.nunique() < 3:
                continue
            FT, RT, dist = share_tests(subv, fams, f"{reg}_{blk}", ncol=ncol, seed=__import__("zlib").crc32(f"{reg}|{blk}|{variant}".encode()))
            FT.insert(0, "variant", variant); RT.insert(0, "variant", variant)
            t = FT[FT.tested]
            key = f"{reg}_{blk}_{variant}"
            summ[key] = dict(n_episodes=int(subv.episode.nunique()), n_families=len(FT), n_tested=len(t),
                             median_chi2_df=float(t.chi2_df.median()) if len(t) else None,
                             chi2_df_quartiles=[float(x) for x in t.chi2_df.quantile([0.25, 0.75])] if len(t) else None,
                             n_overdispersed_q05=int((t.q_bh < 0.05).sum()), n_consistent=int((t.q_bh >= 0.05).sum()),
                             n_chi2df_gt3=int((t.chi2_df > 3).sum()), n_exposure_explained=int(t.exposure_explained.sum()),
                             median_disp=float(t.disp.median()) if len(t) else None, n_disp_lt2=int((t.disp < 2).sum()),
                             median_cv_extra=float(t.cv_extra.median()) if len(t) else None,
                             cv_extra_quartiles=[float(x) for x in t.cv_extra.quantile([0.25, 0.75])] if len(t) else None,
                             median_mean_E=float(t.mean_E.median()) if len(t) else None,
                             sumN_tested=float(t.sumN.sum()) if len(t) else 0.0, sumN_all=float(FT.sumN.sum()))
            if variant in ("main", "cut_050", "cut_085", "flag_excl_075"):
                if jk == "wa":
                    fids = fams.fid.astype(int).tolist()
                    J = weighted_jaccard_wa(fids, era); carriers = principal_of(era)
                    fc, flc = "changed_principal", "station_change_principal"
                elif jk == "svi13":
                    fids = fams.fid.astype(int).tolist(); J = set_jaccard_svi(fids, S13); carriers = svi_car
                    fc, flc = "changed_carriers", "station_change"
                else:
                    fids = fams.fid.astype(int).tolist(); J = set_jaccard_svi(fids); carriers = svi_car
                    fc, flc = "changed_carriers", "station_change"
                mt, pairs = mantel(RT, fids, dist, J)
                summ[key]["mantel"] = mt
                if pairs is not None:
                    pairs.insert(0, "block", f"{reg}_{blk}"); pairs.insert(1, "variant", variant); allP.append(pairs.drop(columns=["dist_bin"]))
                    pairs["b"] = pd.cut(pairs.dist_km, [0, 3, 6, 10, 20, 40, 80, 400])
                    summ[key]["corr_bins"] = {str(k_): [round(float(v_["mean"]), 4), int(v_["size"])] for k_, v_ in
                                              pairs.groupby("b", observed=True)["corr"].agg(["mean", "size"]).iterrows()}
                LR = label_rows(RT, subv, fams, dist, f"{reg}_{blk}", carrier_col=fc, flag_col=flc, carriers=carriers)
                summ[key]["row_labels"] = LR[LR.tested].label.value_counts().to_dict()
                if jk == "wa":   # broad (pre-registered) station-change flag as a second label column
                    LB = label_rows(RT, subv, fams, dist, f"{reg}_{blk}", carrier_col="changed_carriers", flag_col="station_change",
                                    carriers={f: set(core[era]) for f in fids})
                    LR["label_broadflag"] = LB.label.to_numpy()
                    summ[key]["row_labels_broadflag"] = LB[LB.tested].label.value_counts().to_dict()
                RT = LR
                if variant == "main":
                    main_labels = LR
            allF.append(FT); allR.append(RT)
            log(f"[consistency] {key}: " + json.dumps({k: v for k, v in summ[key].items() if k not in ("mantel",)}, default=str))
    F = pd.concat(allF, ignore_index=True); R = pd.concat(allR, ignore_index=True)
    F.to_csv(f"{OUT}/consistency_families.csv", index=False)
    R.to_csv(f"{OUT}/consistency_rows.csv", index=False)
    if allP:
        pd.concat(allP, ignore_index=True).to_csv(f"{OUT}/residual_pair_correlations.csv", index=False)
    # SVI fixed-S13 ratio per (family, episode)
    m = TS[TS.major][["fid", "episode", "ep_end", "N"]].merge(TE[TE.major][["fid", "episode", "N", "era"]], on=["fid", "episode"], suffixes=("_S13", "_eff"))
    m = m[m.ep_end < "2021-04-01"].copy()
    m["ratio_S13_eff"] = np.where(m.N_eff >= 5, m.N_S13 / m.N_eff.where(m.N_eff > 0), np.nan)
    out = []
    for f, g in m.groupby("fid"):
        r = g.ratio_S13_eff.dropna()
        med = r.median() if len(r) else np.nan
        sig = 1.4826 * (r - med).abs().median() if len(r) >= 3 else np.nan
        g = g.assign(fam_median=med, fam_sigma_mad=sig)
        if len(r) >= 3 and not sig > 0:
            sig_p = np.sqrt(g.N_S13.clip(lower=1) / g.N_eff.clip(lower=1) ** 2 + g.N_S13.clip(lower=1) ** 2 / g.N_eff.clip(lower=1) ** 3)
            g["sigma_used"] = sig_p
        else:
            g["sigma_used"] = sig
        g["set_affected"] = (g.ratio_S13_eff - med).abs() > 2 * g.sigma_used
        out.append(g)
    SR = pd.concat(out, ignore_index=True)
    SR.to_csv(f"{OUT}/svi_fixedS13_ratio.csv", index=False)
    ep = SR.groupby("episode").agg(era=("era", "first"), N_S13=("N_S13", "sum"), N_eff=("N_eff", "sum"),
                                   n_rows_ratio=("ratio_S13_eff", lambda x: int(x.notna().sum())),
                                   n_set_affected=("set_affected", "sum")).reset_index()
    ep["ratio"] = ep.N_S13 / ep.N_eff
    ep.to_csv(f"{OUT}/svi_fixedS13_ratio_by_episode.csv", index=False)
    summ["svi_set_affected"] = dict(rows_with_ratio=int(SR.ratio_S13_eff.notna().sum()), n_set_affected=int(SR.set_affected.sum()),
                                    by_episode=ep.to_dict("records"))
    json.dump(summ, open(f"{OUT}/consistency_summary.json", "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))
    log("[consistency] done")



# ------------------------------------------------------------------------------------------------ figures
def figures(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                         "ytick.color": "#52514e", "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfb",
                         "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb"})
    CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
    div = LinearSegmentedColormap.from_list("div", ["#e34948", "#f0efec", "#2a78d6"])
    div.set_bad("#ffffff")
    TW, TS, TE = load_tables()
    FW, FS = wa_families(), svi_families()
    # ---- fig 1: heatmaps
    fig, axs = plt.subplots(1, 2, figsize=(15, 9), gridspec_kw=dict(width_ratios=[47, 39]))
    for ax, (T, F, title) in zip(axs, ((TW[TW.major & (TW.era == "E3")], FW, "WA-NW, era E3 (2015-2026), MAJOR episodes"),
                                       (TS[TS.major], FS, "SVI, fixed-S13 series (2014-2026), MAJOR episodes"))):
        eps = sorted(T.episode.unique(), key=lambda e: T[T.episode == e].ep_start.iloc[0])
        tot = T.groupby("fid").agg(N=("N", "sum"), X=("X", "sum"))
        keep = tot[(tot.N >= 20) & (tot.X > 0)]
        order = F[F.fid.isin(keep.index)].sort_values("lat", ascending=False).fid.tolist()
        rho = (keep.N / keep.X)
        Nm = T.pivot_table(index="fid", columns="episode", values="N", aggfunc="sum").reindex(index=order, columns=eps)
        Xm = T.pivot_table(index="fid", columns="episode", values="X", aggfunc="sum").reindex(index=order, columns=eps)
        Dm = T.pivot_table(index="fid", columns="episode", values="D", aggfunc="sum").reindex(index=order, columns=eps)
        V = np.log2((Nm.clip(lower=0) + 0.5) / (rho.reindex(order).to_numpy()[:, None] * Xm + 0.5))
        V = V.where(Dm > 0)
        im = ax.imshow(V.to_numpy(), aspect="auto", cmap=div, vmin=-4, vmax=4, interpolation="nearest")
        ax.set_xticks(range(len(eps))); ax.set_xticklabels(eps, rotation=90, fontsize=7)
        lat = F.set_index("fid").lat.reindex(order).to_numpy()
        yt = np.linspace(0, len(order) - 1, 8).astype(int)
        ax.set_yticks(yt); ax.set_yticklabels([f"{lat[i]:.2f} N" for i in yt])
        ax.set_title(f"{title}\n{len(order)} families with sum N >= 20, sorted north to south", fontsize=9, color="#0b0b0b")
        ax.set_xlabel("episode (start month)")
    cb = fig.colorbar(im, ax=axs, shrink=0.6, pad=0.01)
    cb.set_label("log2 [ (N + 0.5) / (rho_f X + 0.5) ]   (rho_f = family's sum N / sum X; white = family not exposed)")
    fig.suptitle("Per-family LFE count per tremor epicentre, relative to the family's own all-episode rate (blue = more, red = fewer)",
                 fontsize=10, color="#0b0b0b")
    fig.savefig(f"{OUT}/fig1_family_episode_heatmap.png", dpi=150, bbox_inches="tight"); plt.close(fig)
    # ---- fig 2: regional N vs X
    RW = pd.read_csv(f"{OUT}/regional_nw.csv"); RS = pd.read_csv(f"{OUT}/regional_svi_S13fixed.csv")
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), sharey=False)
    ax = axs[0]
    for k, era in enumerate(("E1", "E2", "E3")):
        m = RW[RW.major & (RW.era == era) & (RW.X > 0)]
        ax.scatter(m.X, m.N.clip(lower=0) + 1, s=36, color=CAT[k], edgecolor="#fcfcfb", linewidth=1.5, label=f"{era} ({len(m)} episodes)", zorder=3)
    r = RW[RW.episode == "2024-08"].iloc[0]
    ax.annotate("2024-08", (r.X, r.N + 1), xytext=(8, -12), textcoords="offset points", fontsize=8, color="#52514e")
    ax.set_title("WA-NW: regional N vs X per MAJOR episode", fontsize=10)
    ax2 = axs[1]
    m = RS[RS.major & (RS.X > 0)]
    ax2.scatter(m.X, m.N.clip(lower=0) + 1, s=36, color=CAT[0], edgecolor="#fcfcfb", linewidth=1.5, label=f"fixed-S13 ({len(m)} episodes)", zorder=3)
    for ep, lab in (("2017-03", "2017-03"), ("2025-09b1", "2025-10")):
        r = RS[RS.episode == ep].iloc[0]
        ax2.annotate(lab, (r.X, r.N + 1), xytext=(8, -4), textcoords="offset points", fontsize=8, color="#52514e")
    ax2.set_title("SVI: regional N vs X per MAJOR episode", fontsize=10)
    for a_ in axs:
        a_.set_xscale("log"); a_.set_yscale("log"); a_.grid(True, color="#e6e5e1", linewidth=0.6, zorder=0)
        a_.set_xlabel("X = sum over families of PNSN epicentres <= 30 km in the family windows")
        a_.set_ylabel("N + 1 (repeats at P >= 0.75 minus control)")
        a_.legend(frameon=False, loc="upper left")
    fig.tight_layout(); fig.savefig(f"{OUT}/fig2_regional_N_vs_X.png", dpi=150); plt.close(fig)
    # ---- fig 3: residual correlation vs distance and carrier overlap
    P = pd.read_csv(f"{OUT}/residual_pair_correlations.csv")
    if "variant" in P:
        P = P[P.variant == "main"]
    fig, axs = plt.subplots(1, 2, figsize=(12, 4.5))
    bins = [0, 3, 6, 10, 20, 40, 80, 160]
    for k, (blk, lab) in enumerate((("NW_E3", "WA-NW E3"), ("SVI_S13fixed_S13fixed", "SVI fixed-S13"))):
        q = P[P.block == blk].copy()
        q["b"] = pd.cut(q.dist_km, bins)
        g = q.groupby("b", observed=True)["corr"].agg(["mean", "std", "size"])
        mid = [np.sqrt(max(iv.left, 1) * iv.right) for iv in g.index]
        axs[0].errorbar(mid, g["mean"], yerr=g["std"] / np.sqrt(g["size"]), color=CAT[k], marker="o", ms=6, lw=2, capsize=3, label=lab)
        jb = np.quantile(q.jaccard.dropna(), np.linspace(0, 1, 7))
        q["jb"] = pd.cut(q.jaccard, np.unique(jb), include_lowest=True)
        g2 = q.groupby("jb", observed=True)["corr"].agg(["mean", "std", "size"])
        axs[1].errorbar([iv.mid for iv in g2.index], g2["mean"], yerr=g2["std"] / np.sqrt(g2["size"]), color=CAT[k], marker="o", ms=6,
                        lw=2, capsize=3, label=lab)
    axs[0].set_xscale("log"); axs[0].set_xlabel("family separation (km, bin geometric centre)")
    axs[1].set_xlabel("carrier overlap (WA: weighted Jaccard of template weights;\nSVI: Jaccard of template station sets)")
    for a_ in axs:
        a_.axhline(0, color="#52514e", lw=0.8); a_.set_ylabel("mean residual correlation across episodes (+- s.e.)")
        a_.grid(True, color="#e6e5e1", linewidth=0.6); a_.legend(frameon=False)
    axs[0].set_title("Share-test residuals: correlation vs distance", fontsize=10)
    axs[1].set_title("Share-test residuals: correlation vs carrier overlap", fontsize=10)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig3_residual_corr_vs_distance.png", dpi=150); plt.close(fig)
    log("[figures] done")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    ap = argparse.ArgumentParser()
    ap.add_argument("what")
    ap.add_argument("--rule", default="")
    ap.add_argument("--tag", default="")
    ap.add_argument("--nproc", default="")
    a = ap.parse_args()
    globals()[a.what](a)
