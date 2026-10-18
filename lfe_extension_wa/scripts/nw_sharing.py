#!/usr/bin/env python
"""NW event sharing between families: are some NW families not distinct sources on the permanent network?
Frozen rules: results/nw_sharing/PLAN.md (written before any pair count). Writes only under results/nw_sharing/.
Never deletes; imports wa_post_dedup / wa_prod_post_v2 / episode_rates read-only (none is edited or re-run).

  python nw_sharing.py pairs      # pass over conf/P/prod: per-pair same-event counts + chance (shifted) + NW headline rows
  python nw_sharing.py table      # pair table with location / catalogue / station columns
  python nw_sharing.py graph      # links (pre-registered rule) -> source groups; sensitivity; what groups follow
  python nw_sharing.py array      # array-era waveform test + array-era co-detection for linked + control pairs
  python nw_sharing.py episodes   # per-episode E3 consistency at source-group level (+ singleton acceptance check)
  python nw_sharing.py figures
"""
from __future__ import annotations
import os, sys, json, glob, warnings
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
HERE = os.path.dirname(os.path.abspath(__file__))
OD = f"{ROOT}/lfe_extension_wa/results/nw_sharing"
os.environ["EPR_OUT"] = OD                      # episode_rates reads it at import: its log goes here, nothing else
sys.path.insert(0, HERE)
import episode_rates as ER                      # noqa: E402  read-only
import wa_common as W                           # noqa: E402
import wa_post_dedup as DD                      # noqa: E402  read-only (share_mats, keep_flags)
from numba import njit                          # noqa: E402

WAP = f"{ROOT}/lfe_extension_wa/results/wa_prod"
CAT = f"{ROOT}/data/catalogue/families_v1_2.csv"
TRUTH = f"{ROOT}/lfe_array_discovery/results/ab_royer/truth/royer2014_templates_digitized.csv"
W_S, ID_TOL, HCUT = 2.0, 0.5, 0.75
SHIFTS = (60, -60, 120, -120, 300, -300, 600, -600)
# pre-registered link rule (PLAN.md sec. 2)
LINK_S, LINK_X, LINK_N = 0.30, 5.0, 10
NPROC = int(os.environ.get("NWS_NPROC", "6"))
LOGF = f"{OD}/nw_sharing.log"


def log(*a):
    os.makedirs(OD, exist_ok=True)
    W.log(LOGF, *a)


def nw_fams():
    F = pd.read_csv(f"{WAP}/setup/families.csv")
    return F[F.reg == "NW"].reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ pass over P tables
@njit(cache=True)
def _pairs(tq, fq, tp, fp, hp, Doff, Sh, w, tol):
    """for each query i (time tq, family fq) and pool row j (tp sorted) within +-w of a different family: bits
    1 = in window, 2 = arrival-consistent (Doff finite, |dt - Doff| <= tol), 4 = station sharing >= 0.5, 8 = partner headline.
    Returns flat arrays (query index, partner family, bits)."""
    qi = []; pf = []; bt = []
    for i in range(len(tq)):
        lo = np.searchsorted(tp, tq[i] - w)
        hi = np.searchsorted(tp, tq[i] + w, side="right")
        for j in range(lo, hi):
            if fp[j] == fq[i]:
                continue
            b = 1
            d = Doff[fq[i], fp[j]]
            if d == d and abs(tq[i] - tp[j] - d) <= tol:
                b |= 2
            if Sh[fq[i], fp[j]]:
                b |= 4
            if hp[j]:
                b |= 8
            qi.append(i); pf.append(fp[j]); bt.append(b)
    return np.array(qi, np.int64), np.array(pf, np.int64), np.array(bt, np.int64)


_G = {}


def _init():
    _G["S"] = DD.share_mats()
    _G["nw"] = set(nw_fams().fid.astype(int))


def _agg(qi, pf, bt, fq):
    """-> per (fid, partner) counts of distinct queries for each condition."""
    if not len(qi):
        return pd.DataFrame(columns=["fid", "partner", "n_win", "n_arr", "n_D", "n_share_win", "n_D_hh"])
    X = pd.DataFrame(dict(q=qi, partner=pf, win=1, arr=(bt & 2) > 0, sh=(bt & 4) > 0, hh=(bt & 8) > 0))
    X["D"] = X.arr & X.sh
    X["share_win"] = X.sh
    X["D_hh"] = X.D & X.hh
    G = X.groupby(["q", "partner"])[["win", "arr", "D", "share_win", "D_hh"]].max().reset_index()
    G["fid"] = fq[G.q.to_numpy()]
    return G.groupby(["fid", "partner"])[["win", "arr", "D", "share_win", "D_hh"]].sum().rename(
        columns={"win": "n_win", "arr": "n_arr", "D": "n_D", "share_win": "n_share_win", "D_hh": "n_D_hh"}).reset_index()


def _month(path):
    if not _G:
        _init()
    nw = _G["nw"]
    cols = ["date", "net", "series", "fid", "k", "t", "net_cc", "P_own_WAVE_PNET", "tpl_insample"]
    D = pd.read_parquet(path, columns=cols)
    D = D[(D.net == "core") & D.fid.isin(nw) & D.series.isin(["data", "rev"])].drop(columns=["net"]).reset_index(drop=True)
    if not len(D):
        return None
    era = ER.wa_era(D.date.iloc[0])
    assert era == ER.wa_era(D.date.iloc[-1])
    Sm = _G["S"][era]
    Doff = np.ascontiguousarray(Sm["Doff"]); Sh = np.ascontiguousarray(Sm["A"])
    D["hl"] = D.P_own_WAVE_PNET.to_numpy() >= HCUT
    # NW headline rows (data + rev), for the group-level counts
    H = D[D.hl].copy()
    H["half"] = W.half_of(H.t.to_numpy()); H["ins"] = H.tpl_insample.fillna(False).astype(bool)
    H["era"] = era
    ymd = H.date.str.replace("-", "", regex=False).astype(np.int64).to_numpy()
    H["det_id"] = ymd * 10_000_000 + H.k.to_numpy().astype(np.int64)
    H = H.drop(columns=["tpl_insample", "hl"])
    d = D[D.series == "data"]
    o = np.argsort(d.t.to_numpy(float), kind="stable")
    tp = d.t.to_numpy(float)[o]; fp = d.fid.to_numpy().astype(np.int64)[o]; hp = d.hl.to_numpy()[o]
    q = d[d.hl]
    tq = q.t.to_numpy(float); fq = q.fid.to_numpy().astype(np.int64)
    out = []
    a = _agg(*_pairs(tq, fq, tp, fp, hp, Doff, Sh, W_S, ID_TOL), fq)
    a["shift"] = 0; out.append(a)
    for sh in SHIFTS:
        a = _agg(*_pairs(tq, fq, tp + sh, fp, hp, Doff, Sh, W_S, ID_TOL), fq)
        a["shift"] = sh; out.append(a)
    C = pd.concat(out, ignore_index=True)
    C["era"] = era
    Nq = q.groupby("fid").size().rename("N").reset_index(); Nq["era"] = era
    # directional variant-D flag winners (post_v2/dedup/flags), NW headline data rows
    fl = f"{WAP}/post_v2/dedup/flags/flags_{os.path.basename(path)[2:9]}.parquet"
    Fw = pd.DataFrame()
    if os.path.exists(fl):
        Fl = pd.read_parquet(fl, columns=["det_id", "series", "fid", "dup_D", "win_fid_D"])
        Fl = Fl[(Fl.series == "data") & Fl.dup_D]
        hd = H[H.series == "data"][["det_id"]]
        Fw = Fl.merge(hd, on="det_id").groupby(["fid", "win_fid_D"]).size().rename("n_Dflag").reset_index()
        Fw["era"] = era
    return C, Nq, H, Fw


def pairs(a):
    from multiprocessing import Pool
    fs = sorted(glob.glob(f"{WAP}/conf/P/prod/P_*.parquet"))
    with Pool(NPROC, initializer=_init) as pool:
        res = [r for r in pool.imap(_month, fs) if r is not None]
    C = pd.concat([r[0] for r in res], ignore_index=True)
    for c in ("n_win", "n_arr", "n_D", "n_share_win", "n_D_hh"):
        C[c] = pd.to_numeric(C[c]).astype(np.int64)
    C = C.groupby(["era", "shift", "fid", "partner"]).sum(numeric_only=True).reset_index()
    C.to_parquet(f"{OD}/pair_counts_raw.parquet", index=False)
    N = pd.concat([r[1] for r in res], ignore_index=True).groupby(["era", "fid"]).N.sum().reset_index()
    N.to_csv(f"{OD}/headline_N_by_era.csv", index=False)
    H = pd.concat([r[2] for r in res], ignore_index=True)
    H.to_parquet(f"{OD}/nw_headline_rows.parquet", index=False)
    Fw = pd.concat([r[3] for r in res if len(r[3])], ignore_index=True)
    Fw.groupby(["era", "fid", "win_fid_D"]).n_Dflag.sum().reset_index().to_csv(f"{OD}/dflag_winner_counts.csv", index=False)
    # check: headline data counts == episode_rates count file (n075 data, all halves)
    O = pd.read_parquet(f"{ER.W.RES}/episode_rates/counts_wa_nw_fid_date_half.parquet")
    o = O.groupby(["fid", "series"]).n075.sum()
    m = H.groupby(["fid", "series"]).size()
    J = pd.concat([o.rename("epr"), m.rename("mine")], axis=1).fillna(0)
    chk = dict(rows=len(J), rows_differ=int((J.epr != J.mine).sum()), total_epr=float(J.epr.sum()), total_mine=float(J.mine.sum()))
    json.dump(chk, open(f"{OD}/check_headline_vs_episode_rates.json", "w"), indent=1)
    log("[pairs]", len(fs), "months;", len(C), "pair rows;", json.dumps(chk))


# ------------------------------------------------------------------------------------------------ pair table
def _hav(la1, lo1, la2, lo2):
    p1, p2 = np.radians(la1), np.radians(la2)
    a = np.sin((p2 - p1) / 2) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(np.radians(lo2 - lo1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def royer_match(F):
    from scipy.optimize import linear_sum_assignment
    T = pd.read_csv(TRUTH)
    T = T[(T.subarray == "NW") & (T.locator == "Hyp2000")]
    ok = F.lat_use.notna().to_numpy()
    Dm = _hav(F.lat_use.to_numpy()[ok][:, None], F.lon_use.to_numpy()[ok][:, None], T.lat.to_numpy()[None], T.lon.to_numpy()[None])
    r, c = linear_sum_assignment(Dm)
    km = np.full(len(F), np.nan); idx = np.flatnonzero(ok)
    km[idx[r]] = Dm[r, c]
    nn = np.full(len(F), np.nan); nn[ok] = Dm.min(1)
    return km, nn


def table(a):
    F = nw_fams()
    C = pd.read_csv(CAT)
    cc = ["family_uid", "source_run", "family_id", "union_nn_km", "union_flag_lt2km", "dup_group", "location_quality", "erh",
          "erz", "n_picks", "n_clear_sta", "depth_km", "lat", "lon", "stack_file", "det_source", "perm_carriers"]
    F = F.merge(C[cc].rename(columns={"depth_km": "depth_cat", "perm_carriers": "perm_carriers_cat"}), on="family_uid", how="left", suffixes=("", "_cat"))
    F["run"] = F.source_run.str.extract(r"\((\w+);")[0]
    F["royer_match_km"], F["royer_nn_km"] = royer_match(F)
    F["royer_match"] = F.royer_match_km <= 5.0
    F.to_csv(f"{OD}/nw_families_meta.csv", index=False)
    N = pd.read_csv(f"{OD}/headline_N_by_era.csv")
    R = pd.read_parquet(f"{OD}/pair_counts_raw.parquet")
    R = R[R.partner.isin(F.fid)]
    Dw = pd.read_csv(f"{OD}/dflag_winner_counts.csv")
    S = DD.share_mats()
    Tl = pd.read_csv(f"{WAP}/setup/tpl_library.csv")
    core = json.load(open(f"{WAP}/setup/core_sets.json"))
    fi = F.set_index("fid")
    fids = F.fid.astype(int).to_numpy()
    rows = []
    for era in ("E1", "E2", "E3"):
        Ne = N[N.era == era].set_index("fid").N.reindex(fids).fillna(0)
        Re = R[R.era == era]
        R0 = Re[Re["shift"] == 0].set_index(["fid", "partner"])
        Rs = Re[Re["shift"] != 0].groupby(["fid", "partner"])[["n_win", "n_arr", "n_D", "n_share_win", "n_D_hh"]].sum() / len(SHIFTS)
        Dwe = Dw[Dw.era == era].set_index(["fid", "win_fid_D"]).n_Dflag
        st = Tl[Tl.netsta.isin(core[era])].groupby(["fid", "netsta"]).w.max().unstack(fill_value=0.0)
        st = st.reindex(index=fids, fill_value=0.0)
        Wt = st.to_numpy(); stn = np.array(st.columns)
        for ii, i in enumerate(fids):
            for jj in range(ii + 1, len(fids)):
                j = fids[jj]
                Ni, Nj = Ne[i], Ne[j]
                r = dict(era=era, fid_i=i, fid_j=j, N_i=Ni, N_j=Nj)
                for c in ("n_win", "n_arr", "n_D", "n_share_win", "n_D_hh"):
                    r[f"{c}_ij"] = float(R0[c].get((i, j), 0)); r[f"{c}_ji"] = float(R0[c].get((j, i), 0))
                    r[f"{c}_ij_chance"] = float(Rs[c].get((i, j), 0)); r[f"{c}_ji_chance"] = float(Rs[c].get((j, i), 0))
                r["n_Dflag_i_by_j"] = float(Dwe.get((i, j), 0)); r["n_Dflag_j_by_i"] = float(Dwe.get((j, i), 0))
                mn = np.minimum(Wt[ii], Wt[jj]); mx = np.maximum(Wt[ii], Wt[jj])
                r["wjaccard"] = float(mn.sum() / mx.sum()) if mx.sum() > 0 else np.nan
                sh = stn[(Wt[ii] > 0) & (Wt[jj] > 0)]
                r["shared_stations"] = ";".join(s.split(".")[1] for s in sh)
                r["n_shared"] = len(sh)
                r["n_sta_i"] = int((Wt[ii] > 0).sum()); r["n_sta_j"] = int((Wt[jj] > 0).sum())
                r["top_shared"] = stn[mn.argmax()].split(".")[1] if mn.max() > 0 else ""
                r["share_ok"] = bool(S[era]["A"][i, j])
                rows.append(r)
    P = pd.DataFrame(rows)
    # small / big orientation and sharing fractions
    small_i = P.N_i <= P.N_j
    P["fid_small"] = np.where(small_i, P.fid_i, P.fid_j); P["fid_big"] = np.where(small_i, P.fid_j, P.fid_i)
    P["N_small"] = np.where(small_i, P.N_i, P.N_j); P["N_big"] = np.where(small_i, P.N_j, P.N_i)
    for c in ("n_win", "n_arr", "n_D", "n_D_hh"):
        P[f"{c}_small"] = np.where(small_i, P[f"{c}_ij"], P[f"{c}_ji"])
        P[f"{c}_small_chance"] = np.where(small_i, P[f"{c}_ij_chance"], P[f"{c}_ji_chance"])
        P[f"s_{c[2:]}"] = P[f"{c}_small"] / P.N_small.where(P.N_small > 0)
        P[f"c_{c[2:]}"] = P[f"{c}_small_chance"] / P.N_small.where(P.N_small > 0)
    P["s_i"] = P.n_D_ij / P.N_i.where(P.N_i > 0); P["s_j"] = P.n_D_ji / P.N_j.where(P.N_j > 0)
    P["s"] = P.s_D; P["c"] = P.c_D
    P["ratio_to_chance"] = P.s / P.c.where(P.c > 0)
    # geometry + catalogue
    la, lo = fi.lat_use, fi.lon_use
    P["dist_km"] = _hav(P.fid_i.map(la).to_numpy(float), P.fid_i.map(lo).to_numpy(float), P.fid_j.map(la).to_numpy(float),
                        P.fid_j.map(lo).to_numpy(float))
    P["dz_km"] = P.fid_i.map(fi.depth_km).to_numpy(float) - P.fid_j.map(fi.depth_km).to_numpy(float)
    for c in ("run", "erh", "erz", "n_picks", "location_quality", "union_flag_lt2km", "dup_group", "royer_match", "royer_match_km",
              "family_uid"):
        P[f"{c}_i"] = P.fid_i.map(fi[c]).to_numpy(); P[f"{c}_j"] = P.fid_j.map(fi[c]).to_numpy()
    P["sigma_h_km"] = np.hypot(P.erh_i.astype(float), P.erh_j.astype(float))
    P["sigma_z_km"] = np.hypot(P.erz_i.astype(float), P.erz_j.astype(float))
    P["same_run"] = P.run_i == P.run_j
    P.to_csv(f"{OD}/pair_table.csv", index=False)
    log("[table]", len(P), "pair rows", P.groupby("era").size().to_dict())


# ------------------------------------------------------------------------------------------------ graph
def links(P, s=LINK_S, x=LINK_X, n=LINK_N):
    return (P.s >= s) & (P.s >= x * P.c.fillna(0)) & (P.n_D_small >= n)


def components(fids, E):
    par = {f: f for f in fids}
    def fd(a):
        while par[a] != a:
            par[a] = par[par[a]]; a = par[a]
        return a
    for a_, b_ in E:
        ra, rb = fd(a_), fd(b_)
        if ra != rb:
            par[ra] = rb
    roots = {f: fd(f) for f in fids}
    return roots


def make_groups(P, fids, era="E3", **kw):
    Pe = P[P.era == era]
    L = Pe[links(Pe, **kw)]
    roots = components(list(fids), list(zip(L.fid_i, L.fid_j)))
    G = pd.DataFrame(dict(fid=list(fids), root=[roots[f] for f in fids]))
    sizes = G.groupby("root").size()
    order = sizes.sort_values(ascending=False).index
    gid = {r: f"G{k + 1:03d}" for k, r in enumerate(order)}
    G["group"] = G.root.map(gid); G["group_size"] = G.root.map(sizes)
    return G.drop(columns="root"), L


def graph(a):
    P = pd.read_csv(f"{OD}/pair_table.csv", low_memory=False)
    F = pd.read_csv(f"{OD}/nw_families_meta.csv")
    fids = F.fid.astype(int).tolist()
    G, L = make_groups(P, fids, "E3")
    fm = F.set_index("fid")
    G = G.merge(F[["fid", "family_uid", "run", "lat_use", "lon_use", "depth_km", "erh", "erz", "n_picks", "location_quality",
                   "royer_match", "royer_match_km", "dup_group", "union_flag_lt2km", "perm_carriers"]], on="fid")
    N3 = pd.read_csv(f"{OD}/headline_N_by_era.csv").query("era == 'E3'").set_index("fid").N
    G["N_E3"] = G.fid.map(N3).fillna(0)
    G.to_csv(f"{OD}/source_groups_E3.csv", index=False)
    L.to_csv(f"{OD}/links_E3.csv", index=False)
    rep = dict(rule=f"s >= {LINK_S} and s >= {LINK_X} x chance and n >= {LINK_N} (E3)", n_families=len(fids), n_links=len(L),
               n_groups=int(G.group.nunique()), n_singletons=int((G.group_size == 1).sum()),
               n_families_in_multi=int((G.group_size > 1).sum()),
               group_sizes=G.drop_duplicates("group").group_size.value_counts().sort_index().to_dict())
    # E3 families with headline repeats at all
    act = set(N3[N3 > 0].index)
    rep["n_families_with_E3_repeats"] = len(act)
    rep["n_groups_among_E3_active"] = int(G[G.fid.isin(act)].group.nunique())
    rep["n_families_N_ge_20"] = int((N3 >= 20).sum())
    # per multi-member group
    Pe = P[P.era == "E3"].set_index(["fid_i", "fid_j"])
    gl = []
    for g, m in G[G.group_size > 1].groupby("group"):
        f = sorted(m.fid)
        pp = [(x, y) for k_, x in enumerate(f) for y in f[k_ + 1:]]
        sub = Pe.loc[pp]
        lk = links(sub.reset_index())
        dmax = float(sub.dist_km.max())
        gl.append(dict(group=g, size=len(f), members=";".join(map(str, f)), uids=";".join(m.sort_values("fid").family_uid),
                       runs=";".join(m.sort_values("fid").run.astype(str)), n_member_pairs=len(pp), n_linked=int(lk.sum()),
                       frac_pairs_linked=float(lk.mean()), max_sep_km=dmax, median_sep_km=float(sub.dist_km.median()),
                       all_within_10km=bool(dmax <= 10), median_s=float(sub.s.median()), max_s=float(sub.s.max()),
                       N_E3_total=float(m.N_E3.sum()), lat=float(m.lat_use.mean()), lon=float(m.lon_use.mean()),
                       n_royer_match=int(m.royer_match.fillna(False).sum()),
                       dup_groups_v12=";".join(sorted(set(m.dup_group.dropna().astype(str)))),
                       top_shared=";".join(f"{k}:{v}" for k, v in sub.top_shared.value_counts().head(3).items())))
    GL = pd.DataFrame(gl).sort_values("size", ascending=False)
    GL.to_csv(f"{OD}/source_groups_E3_summary.csv", index=False)
    # what links follow (E3, links vs all pairs)
    Pa = P[(P.era == "E3") & (P.N_small >= LINK_N)]
    lk = links(Pa)
    def frac(m):
        return dict(linked=float(m[lk].mean()), all_pairs=float(m.mean()), n_linked=int(lk.sum()), n_all=int(len(Pa)))
    rep["follow"] = dict(
        within_10km=frac(Pa.dist_km <= 10), within_20km=frac(Pa.dist_km <= 20),
        same_run=frac(Pa.same_run), top_shared_B013_or_B943=frac(Pa.top_shared.isin(["B013", "B943"])),
        union_flag_lt2km_any=frac(Pa.union_flag_lt2km_i.fillna(False).astype(bool) | Pa.union_flag_lt2km_j.fillna(False).astype(bool)),
        dup_group_any=frac(Pa.dup_group_i.notna() | Pa.dup_group_j.notna()),
        both_quality_A=frac((Pa.location_quality_i == "A") & (Pa.location_quality_j == "A")),
        linked_dist_km_quantiles=[float(x) for x in Pa[lk].dist_km.quantile([0.1, 0.25, 0.5, 0.75, 0.9])],
        linked_wjaccard_median=float(Pa[lk].wjaccard.median()), all_wjaccard_median=float(Pa.wjaccard.median()))
    # B013/B943 dominance among the families of multi groups vs singletons
    pc = G.perm_carriers.fillna("")
    G["carries_B013_or_B943"] = pc.str.contains("B013") | pc.str.contains("B943")
    rep["follow"]["carriers_B013_B943"] = dict(multi=float(G[G.group_size > 1].carries_B013_or_B943.mean()),
                                               singletons=float(G[G.group_size == 1].carries_B013_or_B943.mean()))
    # run composition: pairs linked across runs vs within run
    rep["follow"]["linked_pairs_by_run_pair"] = Pa[lk].apply(lambda r: "|".join(sorted([str(r.run_i), str(r.run_j)])), axis=1).value_counts().to_dict()
    # sensitivity
    sens = []
    for era in ("E1", "E2", "E3"):
        Ne_ = P[P.era == era]
        act_e = set(Ne_[Ne_.N_i > 0].fid_i) | set(Ne_[Ne_.N_j > 0].fid_j)
        for s_ in (0.2, 0.3, 0.5):
            Gs, Ls = make_groups(P, fids, era, s=s_)
            sens.append(dict(era=era, link_s=s_, n_links=len(Ls), n_groups=int(Gs.group.nunique()),
                             n_groups_active=int(Gs[Gs.fid.isin(act_e)].group.nunique()), n_active=len(act_e),
                             largest=int(Gs.group_size.max()), n_in_multi=int((Gs.group_size > 1).sum())))
    SN = pd.DataFrame(sens); SN.to_csv(f"{OD}/groups_sensitivity.csv", index=False)
    rep["sensitivity"] = SN.to_dict("records")
    # cross-era persistence of E3 links
    key = set(zip(L.fid_i, L.fid_j))
    pers = {}
    for era in ("E1", "E2"):
        Pe2 = P[P.era == era].set_index(["fid_i", "fid_j"])
        sub = Pe2.loc[list(key)]
        ok = sub.N_small >= LINK_N
        pers[era] = dict(testable=int(ok.sum()), linked=int(links(sub[ok].reset_index()).sum()),
                         s_ge_0p3=int((sub[ok].s >= 0.3).sum()), median_s=float(sub[ok].s.median()) if ok.any() else None,
                         median_s_unlinked_same_dist=None)
    rep["E3_links_in_other_eras"] = pers
    json.dump(rep, open(f"{OD}/groups_summary.json", "w"), indent=1, default=str)
    log("[graph]", json.dumps({k: v for k, v in rep.items() if k not in ("sensitivity",)}, default=str))
    print(json.dumps(rep, indent=1, default=str)); print(GL.to_string())


# ------------------------------------------------------------------------------------------------ array-era tests
def _load_stacks(path, cache={}):
    if path not in cache:
        z = np.load(f"{ROOT}/{path}"); S = {}
        for k in z.files:
            g, s = k.split("_", 1); S.setdefault(int(g), {})[s] = z[k].astype(float)
        cache[path] = S
    return cache[path]


def _wf_dup(A, B, FS=40.0):
    """NW union rule (lfe_array_discovery/scripts/nw_union.py dup()), unchanged."""
    from scipy.signal import correlate
    ML = int(3 * FS)
    lags, ccs, allcc = [], [], []
    for s in set(A) & set(B):
        best = (0, 0)
        for c in (0, 1):
            x, y = A[s][c], B[s][c]; n = np.linalg.norm(x) * np.linalg.norm(y)
            if n == 0:
                continue
            cc = correlate(x, y, "full") / n; mid = len(y) - 1; w = cc[mid - ML:mid + ML + 1]; k = int(np.argmax(w))
            if w[k] > best[0]:
                best = (w[k], (k - ML) / FS)
        allcc.append(best[0])
        if best[0] >= 0.5:
            ccs.append(best[0]); lags.append(best[1])
    spread = float(np.ptp(lags)) if len(lags) > 1 else np.nan
    verdict = "duplicate" if (len(ccs) >= 2 and spread <= 0.15) else ("neighbour" if (len(ccs) >= 2 and spread > 0.3) else "unclear")
    return dict(wf_n_common=len(set(A) & set(B)), wf_n_good=len(ccs), wf_lag_spread=spread, wf_cc_med=float(np.median(allcc)) if allcc else np.nan,
                wf_cc_max=float(np.max(allcc)) if allcc else np.nan, wf_verdict=verdict)


def _load_det(path, cache={}):
    if path not in cache:
        cache[path] = pd.read_parquet(f"{ROOT}/{path}", columns=["grp", "t"])
    return cache[path]


def _codet(ta, tb):
    """array-era co-detection: modal offset (0.1 s bins) of b - a within +-10 s, share of a's detections with a b
    detection at mode +-0.25 s; control = b shifted by +-60, +-120 s (same mode)."""
    tb = np.sort(tb)
    if len(ta) == 0 or len(tb) == 0:
        return dict(codet_n_a=len(ta), codet_frac=np.nan, codet_ctrl=np.nan, codet_excess=np.nan, codet_offset=np.nan)
    lo = np.searchsorted(tb, ta - 10); hi = np.searchsorted(tb, ta + 10)
    dd = np.concatenate([tb[l:h] - x for x, l, h in zip(ta, lo, hi)]) if (hi > lo).any() else np.array([])
    if len(dd) < 5:
        return dict(codet_n_a=len(ta), codet_frac=0.0, codet_ctrl=np.nan, codet_excess=np.nan, codet_offset=np.nan)
    h, e = np.histogram(dd, bins=np.arange(-10, 10.01, 0.1)); mode = e[h.argmax()] + 0.05
    def fr(tbb):
        l = np.searchsorted(tbb, ta + mode - 0.25); u = np.searchsorted(tbb, ta + mode + 0.25, side="right")
        return float((u > l).mean())
    f0 = fr(tb); ctrl = float(np.mean([fr(tb + s) for s in (60, -60, 120, -120)]))
    return dict(codet_n_a=len(ta), codet_frac=f0, codet_ctrl=ctrl, codet_excess=f0 - ctrl, codet_offset=float(mode))


def array(a):
    P = pd.read_csv(f"{OD}/pair_table.csv", low_memory=False)
    F = pd.read_csv(f"{OD}/nw_families_meta.csv").set_index("fid")
    P3 = P[P.era == "E3"].copy()
    L = P3[links(P3)].copy(); L["set"] = "linked"
    # control: unlinked pairs with N_small >= 10, matched to the linked distance distribution (3 per linked pair, nearest distance)
    U = P3[~links(P3) & (P3.N_small >= LINK_N) & (P3.s < 0.05)].copy()
    rng = np.random.default_rng(3); ctl = []
    used = set()
    for d in L.dist_km:
        cand = U[~U.index.isin(used)]
        k = (cand.dist_km - d).abs().nsmallest(12).index.to_numpy()
        pick = rng.choice(k, size=min(3, len(k)), replace=False)
        used.update(pick); ctl.extend(pick)
    Cn = U.loc[ctl].copy(); Cn["set"] = "control_unlinked"
    X = pd.concat([L, Cn], ignore_index=True)
    rows = []
    for r in X.itertuples():
        fa, fb = F.loc[r.fid_i], F.loc[r.fid_j]
        d = dict(fid_i=r.fid_i, fid_j=r.fid_j, set=r.set, dist_km=r.dist_km, s=r.s, c=r.c)
        try:
            A = _load_stacks(fa.stack_file).get(int(fa.family_id)); B = _load_stacks(fb.stack_file).get(int(fb.family_id))
            d.update(_wf_dup(A, B) if (A is not None and B is not None) else dict(wf_verdict="no stack"))
        except Exception as e:  # noqa
            d.update(wf_verdict=f"error {e}")
        da = _load_det(fa.det_source); db = _load_det(fb.det_source)
        ta = da[da.grp == int(fa.family_id)].t.to_numpy(float); tb = db[db.grp == int(fb.family_id)].t.to_numpy(float)
        # common days only (runs cover different days)
        dya = np.floor(ta / 86400); dyb = np.floor(tb / 86400); cd = np.intersect1d(dya, dyb)
        ta = ta[np.isin(dya, cd)]; tb = tb[np.isin(dyb, cd)]
        if len(ta) > len(tb):   # smaller family as 'a'
            ta, tb = tb, ta
        d.update(_codet(ta, tb)); d["codet_common_days"] = len(cd)
        rows.append(d)
    A_ = pd.DataFrame(rows)
    A_.to_csv(f"{OD}/array_era_tests.csv", index=False)
    s = A_.groupby("set").agg(n=("fid_i", "size"), duplicate=("wf_verdict", lambda x: int((x == "duplicate").sum())),
                              neighbour=("wf_verdict", lambda x: int((x == "neighbour").sum())),
                              unclear=("wf_verdict", lambda x: int((x == "unclear").sum())),
                              median_cc_med=("wf_cc_med", "median"), median_codet_excess=("codet_excess", "median"),
                              codet_excess_gt_0p2=("codet_excess", lambda x: int((x > 0.2).sum())),
                              median_dist=("dist_km", "median"))
    s.to_csv(f"{OD}/array_era_tests_summary.csv")
    log("[array]", s.to_string())
    print(s.to_string())



# ------------------------------------------------------------------------------------------------ array-confirmed + location
def confirmed_links():
    """PLAN Revision 1: linked pair is array-confirmed same source if co-detection excess >= 0.1 or waveform duplicate on
    >= 5 stations."""
    A = pd.read_csv(f"{OD}/array_era_tests.csv")
    l = A[A.set == "linked"].copy()
    l["same_source"] = (l.codet_excess >= 0.1) | ((l.wf_verdict == "duplicate") & (l.wf_n_good >= 5))
    return l


def loc(a):
    P = pd.read_csv(f"{OD}/pair_table.csv", low_memory=False)
    F = pd.read_csv(f"{OD}/nw_families_meta.csv")
    fids = F.fid.astype(int).tolist()
    P3 = P[P.era == "E3"].copy()
    lk = links(P3)
    l = confirmed_links()[["fid_i", "fid_j", "same_source", "wf_verdict", "wf_n_good", "wf_lag_spread", "wf_cc_med", "codet_excess",
                            "codet_offset"]]
    L = P3[lk].merge(l, on=["fid_i", "fid_j"], how="left")
    L["within_2sigma_h"] = L.dist_km <= 2 * L.sigma_h_km
    L["within_2sigma_z"] = L.dz_km.abs() <= 2 * L.sigma_z_km
    L["max_erh"] = np.fmax(L.erh_i, L.erh_j); L["min_npicks"] = np.fmin(L.n_picks_i, L.n_picks_j)
    L["worst_quality"] = np.where((L.location_quality_i == "C") | (L.location_quality_j == "C"), "C",
                                  np.where((L.location_quality_i == "B") | (L.location_quality_j == "B"), "B", "A"))
    L.to_csv(f"{OD}/links_E3_classified.csv", index=False)
    U = P3[~lk & (P3.N_small >= LINK_N)].copy()
    U["within_2sigma_h"] = U.dist_km <= 2 * U.sigma_h_km
    rep = {}
    for nm, X in (("linked_same_source", L[L.same_source]), ("linked_distinct_on_array", L[~L.same_source]), ("unlinked_N_ge_10", U)):
        rep[nm] = dict(n=len(X), median_dist_km=float(X.dist_km.median()), frac_within_2sigma_h=float(X.within_2sigma_h.mean()),
                       median_sigma_h=float(X.sigma_h_km.median()), median_max_erh=float(np.fmax(X.erh_i, X.erh_j).median()),
                       median_min_npicks=float(np.fmin(X.n_picks_i, X.n_picks_j).median()),
                       frac_any_quality_C=float(((X.location_quality_i == "C") | (X.location_quality_j == "C")).mean()),
                       frac_both_A=float(((X.location_quality_i == "A") & (X.location_quality_j == "A")).mean()),
                       frac_dist_gt10=float((X.dist_km > 10).mean()), frac_same_run=float(X.same_run.mean()),
                       frac_abs_dz_le_2sigz=float((X.dz_km.abs() <= 2 * X.sigma_z_km).mean()))
    # array-confirmed groups
    E = L[L.same_source]
    roots = components(fids, list(zip(E.fid_i, E.fid_j)))
    Gc = pd.DataFrame(dict(fid=fids, root=[roots[f] for f in fids]))
    sz = Gc.groupby("root").size(); order = sz.sort_values(ascending=False).index
    gid = {r: f"C{k + 1:03d}" for k, r in enumerate(order)}
    Gc["group"] = Gc.root.map(gid); Gc["group_size"] = Gc.root.map(sz)
    Gc.drop(columns="root").to_csv(f"{OD}/source_groups_E3_array_confirmed.csv", index=False)
    N3 = pd.read_csv(f"{OD}/headline_N_by_era.csv").query("era == 'E3'").set_index("fid").N
    act = set(N3[N3 > 0].index)
    rep["array_confirmed_groups"] = dict(n_groups=int(Gc.group.nunique()), n_groups_active=int(Gc[Gc.fid.isin(act)].group.nunique()),
                                         n_in_multi=int((Gc.group_size > 1).sum()),
                                         sizes=Gc.drop_duplicates("group").group_size.value_counts().sort_index().to_dict())
    # within each array-confirmed group: separations (location scatter of copies of one source)
    fi = F.set_index("fid")
    seps = []
    for g, m in Gc[Gc.group_size > 1].groupby("group"):
        f = m.fid.tolist()
        la = fi.lat_use.reindex(f).to_numpy(); lo = fi.lon_use.reindex(f).to_numpy()
        d = _hav(la[:, None], lo[:, None], la[None], lo[None])
        seps.append(dict(group=g, size=len(f), members=";".join(map(str, f)), runs=";".join(fi.run.reindex(f).astype(str)),
                         max_sep_km=float(np.nanmax(d)), erh=";".join(f"{x:.1f}" for x in fi.erh.reindex(f)),
                         quality=";".join(fi.location_quality.reindex(f).astype(str)),
                         n_picks=";".join(str(int(x)) if np.isfinite(x) else "nan" for x in fi.n_picks.reindex(f)),
                         royer_match=";".join(str(bool(x)) for x in fi.royer_match.reindex(f).fillna(False)),
                         depth=";".join(f"{x:.0f}" for x in fi.depth_km.reindex(f))))
    S = pd.DataFrame(seps).sort_values("size", ascending=False)
    S.to_csv(f"{OD}/array_confirmed_groups_summary.csv", index=False)
    rep["array_confirmed_group_max_sep_km_quantiles"] = [float(x) for x in S.max_sep_km.quantile([0.1, 0.5, 0.9])]
    json.dump(rep, open(f"{OD}/location_vs_sharing.json", "w"), indent=1, default=str)
    print(json.dumps(rep, indent=1, default=str)); print(S.to_string())

# ------------------------------------------------------------------------------------------------ episodes (group level)
def _group_rows(H, gmap, era_S):
    """H: NW headline rows (data/rev) with fid. Members of one group count once: a member candidate is dropped if a
    higher-net_cc candidate of ANOTHER member of the same group lies within +-2 s with the arrival check."""
    keep = np.ones(len(H), bool)
    H = H.reset_index(drop=True)
    gm = H.fid.map(gmap).to_numpy()
    for (dt, ser), ix in H.groupby(["date", "series"]).indices.items():
        Sm = era_S[ER.wa_era(dt)]
        g = gm[ix]
        if pd.Series(g).duplicated().sum() == 0:
            continue
        f = H.fid.to_numpy().astype(np.int64)[ix]
        nb = Sm["Doff"].shape[0]
        # same-group matrix (different family, same group)
        gi = np.full(nb, -1, np.int64)
        for fid, gg in gmap.items():
            gi[fid] = int(gg[1:])
        Sg = (gi[:, None] == gi[None, :]) & (gi[:, None] >= 0)
        np.fill_diagonal(Sg, False)
        U = np.arange(nb, dtype=np.int64)
        k, _ = DD.keep_flags(H.t.to_numpy(float)[ix], H.net_cc.to_numpy(float)[ix], f, U, Sg, Sm["Doff"], ID_TOL)
        keep[ix] = k
    return H.assign(keep_group=keep)


def episodes(a):
    from scipy.stats import spearmanr  # noqa
    F = pd.read_csv(f"{OD}/nw_families_meta.csv")
    fids = F.fid.astype(int).tolist()
    G = pd.read_csv(f"{OD}/source_groups_E3.csv")
    H = pd.read_parquet(f"{OD}/nw_headline_rows.parquet")
    H = H[H.era == "E3"].reset_index(drop=True)
    S = DD.share_mats()
    E = pd.read_csv(f"{ER.W.RES}/episode_rates/episodes_nw_final.csv")
    X = pd.read_parquet(f"{ER.W.RES}/episode_rates/exposure_wa_nw.parquet")
    TW = pd.read_csv(f"{ER.W.RES}/episode_rates/robustness/fam_episode_nw.csv", low_memory=False)
    sub = TW[TW.major & (TW.era == "E3") & ~TW.excluded]
    eps = sorted(sub.episode.unique())
    Ee = E[E.episode.isin(eps)]
    win = ER._windows(Ee, X, fids)
    out = {}
    variants = {"family (singleton groups, acceptance)": {f: f"G{f}" for f in fids},
                "source groups E3 (pre-registered)": dict(zip(G.fid, G.group))}
    for s_ in (0.2, 0.5):
        Gs, _ = make_groups(pd.read_csv(f"{OD}/pair_table.csv", low_memory=False), fids, "E3", s=s_)
        variants[f"source groups E3 (link s >= {s_})"] = dict(zip(Gs.fid, Gs.group))
    Gc = pd.read_csv(f"{OD}/source_groups_E3_array_confirmed.csv")
    variants["array-confirmed same-source groups (Revision 1)"] = dict(zip(Gc.fid, Gc.group.str.replace("C", "G")))
    pos = F.set_index("fid")[["lat_use", "lon_use"]]
    Tl = pd.read_csv(f"{WAP}/setup/tpl_library.csv")
    core = json.load(open(f"{WAP}/setup/core_sets.json"))["E3"]
    Wt = Tl[Tl.netsta.isin(core)].groupby(["fid", "netsta"]).w.max().unstack(fill_value=0.0).reindex(index=fids, columns=core).fillna(0)
    allF = []; allP = []
    for name, gmap in variants.items():
        gmap = {int(k): (v if str(v).startswith("G") else f"G{v}") for k, v in gmap.items()}
        # numeric group ids for the keep matrix
        # unit id = smallest member fid (singletons keep their own fid; units ordered by it, as wa_families orders fids)
        gmin = {}
        for f, g in gmap.items():
            gmin[g] = min(gmin.get(g, f), f)
        gm2 = {f: f"G{gmin[g]}" for f, g in gmap.items()}
        Hk = _group_rows(H, gm2, S) if name != "family (singleton groups, acceptance)" else H.assign(keep_group=True)
        Hk = Hk[Hk.keep_group]
        Hk["unit"] = Hk.fid.map(gm2)
        cnt = Hk.groupby(["unit", "date", "series"]).size().unstack("series", fill_value=0)
        units = sorted(set(gm2.values()), key=lambda u: int(u[1:]))
        mem = {u: [f for f, g in gm2.items() if g == u] for u in units}
        rows = []
        for e in Ee.itertuples():
            for u in units:
                days = sorted(set().union(*[set(win[(f, e.episode)][0]) for f in mem[u]]))
                exd = set().union(*[set(win[(f, e.episode)][1]) for f in mem[u]])
                xs = np.nanmean([win[(f, e.episode)][2] for f in mem[u]])
                if days:
                    idx = pd.MultiIndex.from_product([[u], days])
                    c = cnt.reindex(idx).fillna(0)
                    nd = float(c["data"].sum()) if "data" in c else 0.0; nr = float(c["rev"].sum()) if "rev" in c else 0.0
                else:
                    nd = nr = 0.0
                rows.append(dict(fid=int(u[1:]), episode=e.episode, N=nd - nr, n_data=nd, n_rev=nr, D=len(exd), X=xs))
        T = pd.DataFrame(rows)
        fams = pd.DataFrame(dict(fid=[int(u[1:]) for u in units],
                                 lat=[pos.lat_use.reindex(mem[u]).mean() for u in units],
                                 lon=[pos.lon_use.reindex(mem[u]).mean() for u in units]))
        FT, RT, dist = ER.share_tests(T, fams, "NW_E3", ncol="N", seed=__import__("zlib").crc32(b"NW|E3|main"))
        Wg = np.array([Wt.reindex(mem[u]).max(0).to_numpy() for u in units])
        n = len(units); J = np.ones((n, n))
        for i in range(n):
            mn = np.minimum(Wg[i], Wg).sum(1); mx = np.maximum(Wg[i], Wg).sum(1)
            J[i] = np.where(mx > 0, mn / mx, np.nan)
        mt, pairs_ = ER.mantel(RT, fams.fid.tolist(), dist, J)
        t = FT[FT.tested]
        pairs_["b"] = pd.cut(pairs_.dist_km, [0, 3, 6, 10, 20, 40, 80])
        cb = pairs_.groupby("b", observed=True)["corr"].agg(["mean", "size"])
        out[name] = dict(units=n, tested=len(t), overdispersed=int((t.q_bh < 0.05).sum()), consistent=int((t.q_bh >= 0.05).sum()),
                         frac_overdispersed=float((t.q_bh < 0.05).mean()), median_chi2_df=float(t.chi2_df.median()),
                         median_cv_extra=float(t.cv_extra.median()), sumN_tested=float(t.sumN.sum()), sumN_all=float(FT.sumN.sum()),
                         total_data=float(T.n_data.sum()), coef_logdist=mt.get("coef_logdist"), p_logdist=mt.get("p_logdist"),
                         coef_carrier=mt.get("coef_jaccard"), p_carrier=mt.get("p_jaccard"), n_pairs=mt.get("n_pairs"),
                         corr_bins={str(k): [round(float(v["mean"]), 3), int(v["size"])] for k, v in cb.iterrows()})
        FT.insert(0, "variant", name); pairs_.insert(0, "variant", name)
        allF.append(FT); allP.append(pairs_.drop(columns="b"))
        log(f"[episodes] {name}: " + json.dumps({k: v for k, v in out[name].items() if k != "corr_bins"}, default=str))
    pd.concat(allF).to_csv(f"{OD}/episode_consistency_groups.csv", index=False)
    pd.concat(allP).to_csv(f"{OD}/episode_residual_pairs_groups.csv", index=False)
    # reference rows: family level V0 and de-duplicated (S0b)
    R0 = pd.read_csv(f"{ER.W.RES}/episode_rates/robustness/robustness_table.csv")
    r0 = R0[(R0.variant == "V0_rerun_075") & (R0.block == "NW_E3")].iloc[0].to_dict()
    BA = pd.read_csv(f"{ROOT}/svi_science/results/s0b/before_after.csv")
    rd = BA[(BA.block == "WA-NW E3") & (BA.cut == 0.75) & (BA.when == "after (de-dup)")].iloc[0].to_dict()
    out["_reference_family_V0"] = r0; out["_reference_dedup_S0b"] = rd
    json.dump(out, open(f"{OD}/episode_consistency_groups.json", "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


# ------------------------------------------------------------------------------------------------ figures
def figures(a):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    P = pd.read_csv(f"{OD}/pair_table.csv", low_memory=False)
    F = pd.read_csv(f"{OD}/nw_families_meta.csv").set_index("fid")
    G = pd.read_csv(f"{OD}/source_groups_E3.csv").set_index("fid")
    P3 = P[P.era == "E3"]
    lk = links(P3)
    # fig 1: map
    fig, ax = plt.subplots(figsize=(8, 9))
    show = P3[(P3.s >= 0.1) & (P3.s >= LINK_X * P3.c.fillna(0)) & (P3.n_D_small >= LINK_N)].sort_values("s")
    cm = plt.get_cmap("viridis_r")
    for r in show.itertuples():
        a_, b_ = F.loc[r.fid_i], F.loc[r.fid_j]
        ax.plot([a_.lon_use, b_.lon_use], [a_.lat_use, b_.lat_use], color=cm(min(r.s, 1.0)), lw=0.6 + 2.5 * r.s,
                alpha=0.9 if r.s >= LINK_S else 0.35, zorder=2)
    multi = G.group_size > 1
    ax.scatter(F.lon_use[~multi.reindex(F.index).fillna(False)], F.lat_use[~multi.reindex(F.index).fillna(False)], s=14,
               c="0.55", zorder=3, label="family (singleton group)")
    gg = G[multi]
    groups = sorted(gg.group.unique())
    tab = plt.get_cmap("tab20")
    for k, g in enumerate(groups):
        m = gg[gg.group == g].index
        ax.scatter(F.lon_use[m], F.lat_use[m], s=30, color=tab(k % 20), edgecolor="k", lw=0.4, zorder=4)
        ax.annotate(g, (F.lon_use[m].mean(), F.lat_use[m].mean()), fontsize=7, zorder=5)
    sm = plt.cm.ScalarMappable(cmap=cm, norm=plt.Normalize(0, 1)); sm.set_array([])
    plt.colorbar(sm, ax=ax, shrink=0.6, label="sharing fraction s (E3, smaller family's repeats shared)")
    ax.set_xlabel("lon"); ax.set_ylabel("lat")
    ax.set_title(f"NW families, E3 2015-2026: sharing links (s >= 0.1 thin/faded; linked s >= {LINK_S} bold)\n"
                 f"coloured dots = multi-family source groups")
    ax.set_aspect(1 / np.cos(np.radians(47.6)))
    fig.tight_layout(); fig.savefig(f"{OD}/fig1_sharing_map.png", dpi=150); plt.close(fig)
    # fig 2: s vs distance with chance
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    q = P3[P3.N_small >= LINK_N]
    ax[0].scatter(q.dist_km, q.s, s=4, c="0.6", label="all pairs (N_small >= 10)")
    ax[0].scatter(q.dist_km, q.c, s=3, c="tab:orange", alpha=0.5, label="chance level (shifted)")
    ax[0].scatter(q[links(q)].dist_km, q[links(q)].s, s=14, c="tab:red", label="linked")
    ax[0].axhline(LINK_S, color="k", ls="--", lw=0.8)
    ax[0].set_xlabel("distance (km)"); ax[0].set_ylabel("sharing fraction s"); ax[0].set_xlim(0, 80); ax[0].legend(fontsize=8)
    ax[0].set_title("E3 sharing fraction vs distance")
    bins = [0, 3, 6, 10, 15, 20, 30, 40, 60, 80]
    q = q.assign(b=pd.cut(q.dist_km, bins))
    gb = q.groupby("b", observed=True).agg(s=("s", "mean"), c=("c", "mean"), frac_linked=("s", lambda x: float(links(q.loc[x.index]).mean())))
    xm = [(i.left + i.right) / 2 for i in gb.index]
    ax[1].plot(xm, gb.s, "o-", label="mean s"); ax[1].plot(xm, gb.c, "s--", label="mean chance level")
    ax[1].plot(xm, gb.frac_linked, "^-", label="fraction of pairs linked")
    ax[1].set_xlabel("distance (km)"); ax[1].set_yscale("log"); ax[1].legend(fontsize=8); ax[1].set_title("binned")
    fig.tight_layout(); fig.savefig(f"{OD}/fig2_sharing_vs_distance.png", dpi=150); plt.close(fig)
    # fig 3: linked pairs separation vs combined erh
    L = P3[lk]
    fig, ax = plt.subplots(figsize=(6, 5))
    U = P3[~lk & (P3.N_small >= LINK_N)]
    ax.scatter(U.sigma_h_km, U.dist_km, s=3, c="0.7", label="unlinked")
    ax.scatter(L.sigma_h_km, L.dist_km, s=16, c="tab:red", label="linked")
    xx = np.linspace(0, 60, 10); ax.plot(xx, 2 * xx, "k--", lw=0.8, label="d = 2 sigma")
    ax.set_xlim(0, 60); ax.set_ylim(0, 80); ax.set_xlabel("sqrt(erh_i^2 + erh_j^2) (km)"); ax.set_ylabel("separation (km)")
    ax.legend(fontsize=8); fig.tight_layout(); fig.savefig(f"{OD}/fig3_separation_vs_erh.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("what")
    a = ap.parse_args()
    globals()[a.what](a)
