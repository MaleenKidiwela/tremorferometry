#!/usr/bin/env python
"""S0: shared foundations table for the SVI science studies (svi_science/results/s0/PLAN_S0.md).

  python s0_build.py groups     # coordinate frame + pre-registered family groups (locations only)
  python s0_build.py famday     # per run x family x day detectability table
  python s0_build.py table      # candidate tables (prod, prod_S13) + acceptance vs episode_rates fam_episode tables
  python s0_build.py counts     # families.csv totals, group x MAJOR episode counts and floors
  python s0_build.py figures
Writes only under svi_science/results/s0/ and svi_science/figures/s0_*.png. Never deletes.
"""
from __future__ import annotations
import os, sys, json, glob, warnings
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/s0"
FIG = f"{ROOT}/svi_science/figures"
EPR = f"{ROOT}/lfe_extension_wa/results/episode_rates"
REF = f"{ROOT}/lfe_extension/results/refine"
PBC = f"{REF}/conf/P_by_candidate.parquet"
FAMS12 = f"{ROOT}/lfe_catalog_repair/results/e5r_families_2012.csv"
TQC = f"{ROOT}/lfe_catalog_repair/results/e5r_template_qc_2012.csv"
SLAB = f"{ROOT}/data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz"
MC = f"{ROOT}/lfe_extension/results/trackD/mc_family_day.parquet"
S13 = ["PGC", "VGZ", "NLLB", "B001", "B003", "B004", "B005", "B006", "B007", "B009", "B010", "B011", "B926"]
KM = 111.195
CUT_TRY = [10, 12.5, 15, 17.5, 20, 25, 30]
STEPS = [("2015-01-12", "PGC", "pgc_resp_2015_01_12"), ("2017-02-23", "NLLB", "nllb_bh_hh_2017_02_23"),
         ("2017-03-03", "VGZ", "vgz_swap_2017_03_03"), ("2017-08-10", "PGC", "pgc_bh_hh_2017_08_10")]


def log(*a):
    os.makedirs(OUT, exist_ok=True)
    with open(f"{OUT}/s0.log", "a") as f:
        f.write(pd.Timestamp.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ ") + " ".join(str(x) for x in a) + "\n")
    print(*a)


def families():
    E = pd.read_csv(FAMS12)
    fams = pd.read_parquet(PBC, columns=["fam"]).fam.unique()
    return E[E.bkey.isin(fams)].rename(columns={"bkey": "fid"}).reset_index(drop=True)


# ------------------------------------------------------------------------------------------------ frame + groups
def slab_grid():
    S = pd.read_csv(SLAB, header=None, names=["lon", "lat", "z"])
    S["lon"] = S.lon - 360.0
    lons = np.sort(S.lon.unique()); lats = np.sort(S.lat.unique())
    Z = S.pivot_table(index="lat", columns="lon", values="z", dropna=False).reindex(index=lats, columns=lons).to_numpy()
    return lons, lats, Z          # Z negative km (Slab2 convention)


def _bilin(lons, lats, Z, lo, la):
    i = np.searchsorted(lons, lo) - 1; j = np.searchsorted(lats, la) - 1
    fx = (lo - lons[i]) / (lons[i + 1] - lons[i]); fy = (la - lats[j]) / (lats[j + 1] - lats[j])
    return ((1 - fx) * (1 - fy) * Z[j, i] + fx * (1 - fy) * Z[j, i + 1] + (1 - fx) * fy * Z[j + 1, i] + fx * fy * Z[j + 1, i + 1])


def frame(F):
    lons, lats, Z = slab_grid()
    lat0, lon0 = F.lat.mean(), F.lon.mean()
    F = F.copy()
    F["slab2_depth_km"] = -np.array([_bilin(lons, lats, Z, lo, la) for lo, la in zip(F.lon, F.lat)])
    # local strike from the depth gradient (finite difference 0.05 deg)
    h = 0.05
    az = []
    for lo, la in zip(F.lon, F.lat):
        dzdx = (_bilin(lons, lats, Z, lo + h, la) - _bilin(lons, lats, Z, lo - h, la)) / (2 * h * KM * np.cos(np.radians(la)))
        dzdy = (_bilin(lons, lats, Z, lo, la + h) - _bilin(lons, lats, Z, lo, la - h)) / (2 * h * KM)
        # Z is negative-down: down-dip direction = direction of decreasing Z = -grad
        az.append(np.degrees(np.arctan2(-dzdx, -dzdy)) % 360)     # azimuth (from north, clockwise) of steepest descent
    F["slab2_downdip_az"] = az
    F["slab2_dip_deg"] = np.nan
    dd = np.radians(np.array(az))
    mdd = np.degrees(np.arctan2(np.sin(dd).mean(), np.cos(dd).mean())) % 360    # down-dip is a direction, not axial
    strike = (mdd - 90) % 360
    x = (F.lon - lon0) * KM * np.cos(np.radians(lat0)); y = (F.lat - lat0) * KM
    F["x_km"] = x; F["y_km"] = y
    s = np.radians(strike); d = np.radians(mdd)
    F["along_strike_km"] = x * np.sin(s) + y * np.cos(s)
    F["downdip_km"] = x * np.sin(d) + y * np.cos(d)
    w, v = np.linalg.eigh(np.cov(np.vstack([x, y])))
    pca_az = np.degrees(np.arctan2(v[0, 1], v[1, 1])) % 180
    meta = dict(origin_lat=float(lat0), origin_lon=float(lon0), mean_downdip_az=float(mdd), strike_az=float(strike),
                strike_rule="circular mean of Slab2 down-dip azimuth at the family epicentres minus 90 deg; along-strike positive "
                            "toward the strike azimuth", pca_major_axis_az_mod180=float(pca_az),
                downdip_az_spread_deg=float(np.degrees(np.sqrt(-2 * np.log(np.hypot(np.sin(dd).mean(), np.cos(dd).mean()))))))
    return F, meta


def _km_matrix(F):
    x = F.x_km.to_numpy(); y = F.y_km.to_numpy()
    return np.hypot(x[:, None] - x[None], y[:, None] - y[None])


def _merge_small(lab, x, y, minsize=5):
    lab = lab.copy()
    while True:
        u, c = np.unique(lab, return_counts=True)
        if c.min() >= minsize or len(u) == 1:
            break
        small = u[np.argmin(c)]
        cx = {k: (x[lab == k].mean(), y[lab == k].mean()) for k in u}
        others = [k for k in u if k != small]
        tgt = min(others, key=lambda k: np.hypot(cx[k][0] - cx[small][0], cx[k][1] - cx[small][1]))
        lab[lab == small] = tgt
    return lab


def groups(a):
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    F, meta = frame(families())
    D = _km_matrix(F)
    L = linkage(squareform(D, checks=False), method="complete")
    tries = []
    for cut in CUT_TRY:
        lab0 = fcluster(L, t=cut, criterion="distance")
        lab = _merge_small(lab0, F.x_km.to_numpy(), F.y_km.to_numpy())
        tries.append(dict(cut_km=cut, clusters_before_merge=int(len(np.unique(lab0))), groups_after_merge=int(len(np.unique(lab))),
                          min_size=int(np.unique(lab, return_counts=True)[1].min()), lab=lab))
    ok = [t for t in tries if 8 <= t["groups_after_merge"] <= 12]
    pick = ok[0] if ok else min(tries, key=lambda t: (abs(t["groups_after_merge"] - 10), t["cut_km"]))
    lab = pick["lab"]
    cen = pd.DataFrame(dict(lab=lab, s=F.along_strike_km)).groupby("lab").s.mean().sort_values()
    remap = {k: i + 1 for i, k in enumerate(cen.index)}
    F["group"] = [remap[k] for k in lab]
    # group diameter and centroid
    G = F.groupby("group").agg(n_families=("fid", "size"), lat=("lat", "mean"), lon=("lon", "mean"), depth_km=("dep", "mean"),
                               along_strike_km=("along_strike_km", "mean"), downdip_km=("downdip_km", "mean")).reset_index()
    G["max_pair_km"] = [float(D[np.ix_(np.flatnonzero(F.group == g), np.flatnonzero(F.group == g))].max()) for g in G.group]
    meta.update(cut_tries=[{k: v for k, v in t.items() if k != "lab"} for t in tries], cut_chosen_km=pick["cut_km"],
                n_groups=int(F.group.nunique()), linkage="complete", min_group_size=5)
    F[["fid", "lat", "lon", "dep", "mag_med", "slab2_depth_km", "slab2_downdip_az", "x_km", "y_km", "along_strike_km", "downdip_km",
       "group"]].to_csv(f"{OUT}/family_groups.csv", index=False)
    G.to_csv(f"{OUT}/groups.csv", index=False)
    json.dump(meta, open(f"{OUT}/frame_groups.json", "w"), indent=1)
    log("[groups]", json.dumps({k: v for k, v in meta.items()}, default=str))
    print(G.to_string())


# ------------------------------------------------------------------------------------------------ family-day detectability
def _fd_file(p):
    run = p.split("/")[-2]
    D = pd.read_parquet(p, columns=["series", "observable", "mad", "fam", "date", "n_sites", "n_sta"])
    D = D[D.series == "data"].drop(columns=["series"])
    D["run"] = run
    return D


def famday(a):
    from multiprocessing import Pool
    fs = sorted(glob.glob(f"{REF}/feat/prod/*_fd.parquet")) + sorted(glob.glob(f"{REF}/feat/prod_S13/*_fd.parquet"))
    with Pool(16) as pool:
        parts = list(pool.imap(_fd_file, fs, chunksize=50))
    D = pd.concat(parts, ignore_index=True).rename(columns={"fam": "fid", "n_sites": "fam_n_sites", "n_sta": "fam_n_sta"})
    DS = pd.read_csv(f"{EPR}/svi_day_stations.csv")
    D = D.merge(DS.rename(columns={"n_sta": "n_sta_day", "stations": "day_stations"}), on=["run", "date"], how="left")
    med = D[D.observable].groupby(["run", "sset", "fid"]).mad.median().rename("mad_med_set")
    D = D.join(med, on=["run", "sset", "fid"])
    D["u_mad"] = np.log10(D.mad / D.mad_med_set)
    M = pd.read_parquet(MC, columns=["family_id", "date", "station_set", "Mc_T_core", "Pdet_T_core_M1.8", "Mc_T_ext",
                                     "day_terms_in_calib_range"]).rename(columns={"family_id": "fid", "station_set": "mc_station_set"})
    M["run"] = "prod"
    D = D.merge(M, on=["run", "fid", "date"], how="left")
    D.to_parquet(f"{OUT}/famday_detectability.parquet", index=False)
    rep = dict(rows=len(D), by_run=D.run.value_counts().to_dict(), observable_frac=float(D.observable.mean()),
               mc_rows_matched_prod=int(D[D.run == "prod"].Mc_T_core.notna().sum()), prod_rows=int((D.run == "prod").sum()),
               spearman_u_mad_vs_Mc_T_core=float(D[D.run == "prod"][["u_mad", "Mc_T_core"]].corr("spearman").iloc[0, 1]))
    json.dump(rep, open(f"{OUT}/famday_checks.json", "w"), indent=1, default=str)
    log("[famday]", json.dumps(rep, default=str))


# ------------------------------------------------------------------------------------------------ candidate table
def windows():
    """(fid, episode) -> window day set; plus per (fid, date) -> episode (family-window assignment)."""
    E = pd.read_csv(f"{EPR}/episodes_svi_final.csv")
    X = pd.read_parquet(f"{EPR}/exposure_svi.parquet")
    X = X[X.n30 > 0]
    rows = []
    for e in E.itertuples():
        ex = X[(X.date >= e.start) & (X.date <= e.end)]
        for f, g in ex.groupby("fid"):
            days = set()
            for d in g.date:
                t = pd.Timestamp(d)
                days.update(str((t + pd.Timedelta(days=k)).date()) for k in range(-2, 3))
            first = g.date.min()
            for d in days:
                rows.append((int(f), d, e.episode, first, g.date.max()))
    Wd = pd.DataFrame(rows, columns=["fid", "date", "episode", "fam_first_exposure", "fam_last_exposure"])
    dup = Wd.duplicated(["fid", "date"]).sum()
    return E, Wd, int(dup)


def _flags(T):
    d = T.date
    T["flag_ep_2017_03"] = (T.episode == "2017-03")
    T["flag_ep_2025_10"] = (T.episode == "2025-09b1")
    for dd, sta, nm in STEPS:
        lo = str((pd.Timestamp(dd) - pd.Timedelta(days=45)).date()); hi = str((pd.Timestamp(dd) + pd.Timedelta(days=45)).date())
        T[f"near45_{nm}"] = (d >= lo) & (d <= hi)
    T["pgc_epoch"] = np.where(d < "2015-01-12", 0, np.where(d < "2017-08-10", 1, 2)).astype(np.int8)
    T["nllb_hh"] = (d >= "2017-02-23")
    T["vgz_post_swap"] = (d >= "2017-03-03")
    return T


def table(a):
    E, Wd, ndup = windows()
    if ndup:
        raise SystemExit(f"family windows overlap on {ndup} (fid, date) pairs")
    Em = E.set_index("episode")
    maj = E[E.major].sort_values("start").reset_index(drop=True)
    mstart = pd.to_datetime(maj.start).to_numpy("datetime64[ns]")
    mend = (pd.to_datetime(maj.end) + pd.Timedelta(days=1)).to_numpy("datetime64[ns]")   # end of last day
    # episode-level spacing to previous / next MAJOR episode (days, start-to-end gaps)
    sp_prev, sp_next = {}, {}
    for e in E.itertuples():
        s = pd.Timestamp(e.start); en = pd.Timestamp(e.end) + pd.Timedelta(days=1)
        pv = maj[pd.to_datetime(maj.end) < pd.Timestamp(e.start)]
        nx = maj[pd.to_datetime(maj.start) > pd.Timestamp(e.end)]
        sp_prev[e.episode] = (s - (pd.Timestamp(pv.end.iloc[-1]) + pd.Timedelta(days=1))).days if len(pv) else np.nan
        sp_next[e.episode] = (pd.Timestamp(nx.start.iloc[0]) - en).days if len(nx) else np.nan
    X = pd.read_parquet(f"{EPR}/exposure_svi.parquet").rename(columns={"n30": "n_trem30"})
    FD = pd.read_parquet(f"{OUT}/famday_detectability.parquet")
    cols = ["run", "date", "series", "k", "fam", "t", "P_same", "P_same_WAVE_foldsd", "net_mad", "sset_run", "n_sta_day", "has_LRIV",
            "has_SNB", "quiet", "n_tremor_day", "keep_xfam8", "hfr", "n_sta_clbad", "in_v11"]
    P = pd.read_parquet(PBC, columns=cols).rename(columns={"fam": "fid"})
    P["keep_xfam8"] = P.keep_xfam8 > 0.5
    P["time_utc"] = pd.to_datetime(P.t, unit="s").astype("datetime64[ns]")
    P = P.merge(Wd, on=["fid", "date"], how="left")
    P = P.merge(X, on=["fid", "date"], how="left"); P["n_trem30"] = P.n_trem30.fillna(0).astype(np.int32)
    fdc = ["run", "fid", "date", "observable", "mad", "u_mad", "fam_n_sta", "fam_n_sites", "sset", "day_stations", "Mc_T_core",
           "Pdet_T_core_M1.8", "day_terms_in_calib_range"]
    P = P.merge(FD[fdc].rename(columns={"mad": "fam_day_mad", "observable": "fam_day_observable", "sset": "day_sset"}),
                on=["run", "fid", "date"], how="left")
    P["in_episode"] = P.episode.notna()
    P["major"] = P.episode.map(Em.major).fillna(False).astype(bool)
    P["in_frozen_major"] = P.episode.map(Em.in_frozen_major).fillna(False).astype(bool)
    P["ep_start"] = P.episode.map(Em.start); P["ep_end"] = P.episode.map(Em.end)
    tt = P.time_utc.to_numpy("datetime64[ns]")
    day = np.timedelta64(86400, "s")
    P["t_since_ep_start_d"] = (tt - pd.to_datetime(P.ep_start).to_numpy("datetime64[ns]")) / day
    P["t_since_fam_first_exposure_d"] = (tt - pd.to_datetime(P.fam_first_exposure).to_numpy("datetime64[ns]")) / day
    P["ep_gap_prev_major_d"] = P.episode.map(sp_prev); P["ep_gap_next_major_d"] = P.episode.map(sp_next)
    # from the candidate time: since the end of the last MAJOR episode that ended before it, to the start of the next one
    i_prev = np.searchsorted(mend, tt, "right") - 1
    i_next = np.searchsorted(mstart, tt, "right")
    P["d_since_prev_major_end"] = np.where(i_prev >= 0, (tt - mend[np.clip(i_prev, 0, None)]) / day, np.nan)
    P["d_to_next_major_start"] = np.where(i_next < len(mstart), (mstart[np.clip(i_next, None, len(mstart) - 1)] - tt) / day, np.nan)
    P = _flags(P)
    P = P.drop(columns=["t"]).rename(columns={"fid": "fam"})
    for c in ("day_stations", "day_sset", "sset_run", "series", "run", "episode", "ep_start", "ep_end"):
        P[c] = P[c].astype("category")
    for run, nm in (("prod", "analysis_table.parquet"), ("prod_S13", "analysis_table_S13.parquet")):
        Q = P[P.run == run].sort_values(["time_utc", "fam"]).reset_index(drop=True)
        Q.to_parquet(f"{OUT}/{nm}", index=False)
        log(f"[table] {nm}: {len(Q)} rows, series {Q.series.value_counts().to_dict()}")
    accept(P.rename(columns={"fam": "fid"}))


def accept(P):
    """reproduce episode_rates fam_episode_svi_{S13fixed,eff}.csv raw counts."""
    rep = {}
    for series, fn in (("S13fixed", "fam_episode_svi_S13fixed.csv"), ("eff", "fam_episode_svi_eff.csv")):
        if series == "S13fixed":
            Q = P[((P.run == "prod_S13") & (P.date < "2021-04-01")) | ((P.run == "prod") & (P.date >= "2021-04-01"))]
        else:
            Q = P[P.run == "prod"]
        Q = Q[Q.episode.notna()]
        g = {}
        for k, c in (("050", 0.5), ("075", 0.75), ("085", 0.85), ("090", 0.9)):
            for sr in ("data", "c1"):
                g[f"n{k}_{sr}"] = Q[(Q.P_same >= c) & (Q.series == sr)].groupby(["fid", "episode"], observed=True).size()
        G = pd.DataFrame(g).fillna(0)
        O = pd.read_csv(f"{EPR}/{fn}").set_index(["fid", "episode"])
        J = O.join(G, rsuffix="_s0", how="left").fillna({c + "_s0": 0 for c in g})
        r = {c: int((J[c] != J[c + "_s0"]).sum()) for c in g}
        r["N_differ"] = int(((J.n075_data_s0 - J.n075_c1_s0) != J.N).sum())
        r["rows"] = len(J); r["s0_cells_not_in_episode_rates"] = int(len(G.index.difference(O.index)))
        r["total_n075_data"] = [float(J.n075_data.sum()), float(J.n075_data_s0.sum())]
        rep[series] = r
    json.dump(rep, open(f"{OUT}/acceptance_fam_episode.json", "w"), indent=1)
    log("[accept]", json.dumps(rep))


# ------------------------------------------------------------------------------------------------ counts
def _s13fixed(T, T13):
    return pd.concat([T13[T13.date < "2021-04-01"], T[T.date >= "2021-04-01"]], ignore_index=True)


def counts(a):
    cols = ["date", "series", "fam", "P_same", "keep_xfam8", "episode", "major"]
    T = pd.read_parquet(f"{OUT}/analysis_table.parquet", columns=cols)
    T13 = pd.read_parquet(f"{OUT}/analysis_table_S13.parquet", columns=cols)
    FG = pd.read_csv(f"{OUT}/family_groups.csv")
    grp = dict(zip(FG.fid, FG.group))
    # families.csv
    q = pd.read_csv(TQC)
    car = q[q.usable & q.template].groupby("bkey").sta.apply(lambda s: ";".join(sorted(set(s)))).to_dict()
    F = FG.copy()
    F["carriers"] = F.fid.map(car).fillna("")
    F["carriers_S13"] = F.carriers.map(lambda s: ";".join(x for x in s.split(";") if x in S13))
    F["n_carriers"] = F.carriers.map(lambda s: len([x for x in s.split(";") if x]))
    F["n_carriers_S13"] = F.carriers_S13.map(lambda s: len([x for x in s.split(";") if x]))
    for nm, D in (("prod", T), ("S13run", T13), ("S13fixed", _s13fixed(T, T13))):
        for sr in ("data", "c1"):
            m = (D.series == sr) & (D.P_same >= 0.75)
            F[f"n075_{sr}_raw_{nm}"] = F.fid.map(D[m].groupby("fam").size()).fillna(0).astype(int)
            F[f"n075_{sr}_dedup_{nm}"] = F.fid.map(D[m & D.keep_xfam8].groupby("fam").size()).fillna(0).astype(int)
    F.to_csv(f"{OUT}/families.csv", index=False)
    # group x MAJOR episode
    E = pd.read_csv(f"{EPR}/episodes_svi_final.csv"); maj = E[E.major].sort_values("start")
    out = []
    for series, D in (("S13fixed", _s13fixed(T, T13)), ("eff", T)):
        D = D[D.major & D.episode.notna()].copy()
        D["group"] = D.fam.map(grp)
        for var, m_ in (("dedup", D.keep_xfam8), ("raw", pd.Series(True, index=D.index))):
            for cut in (0.75, 0.5):
                Dm = D[m_ & (D.P_same >= cut)]
                g = Dm.groupby(["group", "episode", "series"], observed=True).size().unstack("series", fill_value=0)
                for sr in ("data", "c1"):
                    if sr not in g: g[sr] = 0
                g = g.reindex(pd.MultiIndex.from_product([sorted(FG.group.unique()), maj.episode], names=["group", "episode"]), fill_value=0)
                g = g.reset_index().assign(series_set=series, variant=var, cut=cut)
                g["N"] = g.data - g.c1
                fe = Dm.groupby(["fam", "episode", "series"], observed=True).size().unstack("series", fill_value=0)
                for sr in ("data", "c1"):
                    if sr not in fe: fe[sr] = 0
                fe["N"] = fe.data - fe.c1
                fe = fe.reset_index(); fe["group"] = fe.fam.map(grp)
                g = g.merge(fe[fe.N >= 30].groupby(["group", "episode"]).size().rename("n_fam_cells_ge30").reset_index(),
                            on=["group", "episode"], how="left").fillna({"n_fam_cells_ge30": 0})
                out.append(g)
    C = pd.concat(out, ignore_index=True)
    C["ep_total_N"] = C.groupby(["series_set", "variant", "cut", "episode"]).N.transform("sum")
    C["meets_group_ep_100"] = C.N >= 100
    C["meets_group_ep_30"] = C.N >= 30
    C["episode_meets_300"] = C.ep_total_N >= 300
    C["usable"] = C.meets_group_ep_100 & C.episode_meets_300
    C = C.merge(maj[["episode", "start", "end"]], on="episode", how="left")
    C.to_csv(f"{OUT}/counts_group_episode.csv", index=False)
    summ = {}
    for (ss, var, cut), g in C.groupby(["series_set", "variant", "cut"]):
        ep = g.drop_duplicates("episode")
        summ[f"{ss}|{var}|{cut}"] = dict(total_N=int(g.N.sum()), episodes_ge300=int((ep.ep_total_N >= 300).sum()),
                                          group_ep_cells=len(g), cells_ge100=int(g.meets_group_ep_100.sum()),
                                          cells_ge30=int(g.meets_group_ep_30.sum()), usable_cells=int(g.usable.sum()),
                                          fam_ep_cells_ge30=int(g.n_fam_cells_ge30.sum()),
                                          groups_with_ge5_usable=int((g.groupby("group").usable.sum() >= 5).sum()),
                                          episodes_ge300_list=ep[ep.ep_total_N >= 300].episode.tolist())
    json.dump(summ, open(f"{OUT}/counts_summary.json", "w"), indent=1)
    log("[counts]", json.dumps(summ))


# ------------------------------------------------------------------------------------------------ figures
def figures(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 9, "axes.edgecolor": "#52514e", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                         "ytick.color": "#52514e", "axes.spines.top": False, "axes.spines.right": False, "figure.facecolor": "#fcfcfb",
                         "axes.facecolor": "#fcfcfb", "savefig.facecolor": "#fcfcfb"})
    F = pd.read_csv(f"{OUT}/families.csv"); G = pd.read_csv(f"{OUT}/groups.csv")
    meta = json.load(open(f"{OUT}/frame_groups.json"))
    ng = int(F.group.nunique())
    pal = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#8f5bd6", "#d6457a", "#3fb8c4", "#7a8c1f", "#b5651d", "#52514e",
           "#e34948", "#5c7cfa"]
    lons, lats, Z = slab_grid()
    fig, ax = plt.subplots(figsize=(8, 7.5))
    m = (lons > -125.2) & (lons < -122.4); n = (lats > 47.7) & (lats < 49.4)
    cs = ax.contour(lons[m], lats[n], -Z[np.ix_(n, m)], levels=np.arange(20, 65, 5), colors="#9b9a95", linewidths=0.8)
    ax.clabel(cs, fmt="%d km", fontsize=7)
    for g in sorted(F.group.unique()):
        s = F[F.group == g]
        ax.scatter(s.lon, s.lat, s=28, color=pal[(g - 1) % len(pal)], edgecolor="#fcfcfb", linewidth=0.8, zorder=3,
                   label=f"group {g} ({len(s)} families)")
        c = G[G.group == g].iloc[0]
        ax.text(c.lon, c.lat + 0.035, str(g), fontsize=10, fontweight="bold", color=pal[(g - 1) % len(pal)], ha="center", zorder=4)
    # strike arrow
    o_lo, o_la = meta["origin_lon"], meta["origin_lat"]; s = np.radians(meta["strike_az"])
    L = 25.0
    ax.annotate("", xy=(o_lo + L * np.sin(s) / (KM * np.cos(np.radians(o_la))), o_la + L * np.cos(s) / KM), xytext=(o_lo, o_la),
                arrowprops=dict(arrowstyle="->", color="#0b0b0b", lw=1.4))
    ax.text(o_lo + 1.1 * L * np.sin(s) / (KM * np.cos(np.radians(o_la))), o_la + 1.1 * L * np.cos(s) / KM,
            f"along-strike +\n(az {meta['strike_az']:.0f}°)", fontsize=8, ha="center")
    ax.set_xlim(F.lon.min() - 0.25, F.lon.max() + 0.25); ax.set_ylim(F.lat.min() - 0.15, F.lat.max() + 0.15)
    ax.set_aspect(1 / np.cos(np.radians(o_la)))
    ax.set_xlabel("longitude"); ax.set_ylabel("latitude")
    ax.set_title(f"SVI Bostock families: {ng} pre-registered location groups (complete linkage, cut {meta['cut_chosen_km']} km, "
                 f"groups ≥ 5)\ngrey = Slab2 depth contours (km)", fontsize=9)
    ax.legend(frameon=False, fontsize=7, loc="lower left", ncol=2)
    ax.grid(True, color="#e6e5e1", linewidth=0.5)
    fig.tight_layout(); fig.savefig(f"{FIG}/s0_family_groups_map.png", dpi=150); plt.close(fig)
    # heatmap
    C = pd.read_csv(f"{OUT}/counts_group_episode.csv")
    C = C[(C.series_set == "S13fixed") & (C.variant == "dedup") & (C.cut == 0.75)]
    eps = C.drop_duplicates("episode").sort_values("start")
    M = C.pivot_table(index="group", columns="episode", values="N", aggfunc="sum").reindex(columns=eps.episode)
    import matplotlib.colors as mcolors
    cmap = mcolors.LinearSegmentedColormap.from_list("seq", ["#f0efec", "#9cc3ee", "#2a78d6", "#0d3b75"])
    fig, ax = plt.subplots(figsize=(14, 4.8))
    im = ax.imshow(np.log10(M.clip(lower=0).to_numpy() + 1), aspect="auto", cmap=cmap, vmin=0, vmax=3, interpolation="nearest")
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.iat[i, j]
            ok = v >= 100 and eps.ep_total_N.iat[j] >= 300
            ax.text(j, i, f"{int(v)}", ha="center", va="center", fontsize=6.5,
                    color="#fcfcfb" if v >= 100 else "#0b0b0b", fontweight="bold" if ok else "normal")
            if ok:
                ax.add_patch(plt.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="#eb6834", lw=1.3))
    ax.set_xticks(range(len(eps))); ax.set_xticklabels([f"{e}{' *' if t >= 300 else ''}" for e, t in zip(eps.episode, eps.ep_total_N)],
                                                       rotation=90, fontsize=7)
    ax.set_yticks(range(M.shape[0])); ax.set_yticklabels([f"group {g}" for g in M.index])
    ax.set_xlabel("MAJOR episode (start month; * = episode total ≥ 300)")
    cb = fig.colorbar(im, ax=ax, pad=0.01); cb.set_label("log10(N + 1)")
    ax.set_title("De-duplicated repeats per group × MAJOR episode (P ≥ 0.75, data minus scrambled control, fixed-S13 series). "
                 "Orange box = group-episode ≥ 100 inside an episode ≥ 300", fontsize=9)
    fig.tight_layout(); fig.savefig(f"{FIG}/s0_counts_group_episode.png", dpi=150); plt.close(fig)
    log("[figures] done")


if __name__ == "__main__":
    warnings.filterwarnings("ignore")
    import argparse
    ap = argparse.ArgumentParser(); ap.add_argument("what")
    a = ap.parse_args()
    globals()[a.what](a)
