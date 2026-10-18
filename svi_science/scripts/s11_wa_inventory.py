#!/usr/bin/env python
"""s11_wa_inventory -- S11 step 7 (PLAN §5): Washington inventory in the S10-multi table format, so the unchanged S11
measurement / kernel / inversion code runs on WA. Run with S11_REGION=WA_NW (or WA_SW for the SW column).
Writes C.MTB (= s11_depth/<wa_nw|wa_sw>/tables/inv): families.csv, stations.csv, station_epochs.csv, offsets_sp.csv,
rows.parquet, row_station.parquet, active_days.csv, episode_units.csv; and C.TB/sensor_depths.csv (EarthScope FDSN).
Read-only on every input.
    S11_REGION=WA_NW python s11_wa_inventory.py
"""
from __future__ import annotations
import os, sys, json, time, io, urllib.request
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C

assert C.REGION != "SVI", "set S11_REGION=WA_NW or WA_SW"
R = C.REGION_CFG
WA = f"{C.ROOT}/lfe_extension_wa"
PROD = f"{WA}/results/wa_prod"
HYPO_MAX = 80.0
SITE_KM = 3.0
MIN_EVENTS = 30
SEED = C.SEED + 7


def km(la1, lo1, la2, lo2):
    return np.hypot((np.asarray(lo1) - lo2) * 111.195 * np.cos(np.radians((np.asarray(la1) + la2) / 2)),
                    (np.asarray(la1) - la2) * 111.195)


def slab2_depth(lat, lon):
    from scipy.interpolate import RegularGridInterpolator
    X = pd.read_csv(f"{C.ROOT}/data/slab_models/slab2/cas_slab2_dep_02.24.18.xyz", header=None, names=["lon", "lat", "z"])
    lons, lats = np.sort(X.lon.unique()), np.sort(X.lat.unique())
    Z = X.pivot(index="lat", columns="lon", values="z").reindex(index=lats, columns=lons).to_numpy()
    f = RegularGridInterpolator((lats, lons), Z, bounds_error=False, fill_value=np.nan)
    return -f(np.column_stack([lat, np.mod(lon, 360.0)]))


def families(log):
    cat = pd.read_csv(f"{C.ROOT}/data/catalogue/families_v1_3.csv")
    reg = cat[cat.region == R["famreg"]]
    log(f"catalogue {R['famreg']}: {len(reg)} families, quality {reg.location_quality.value_counts().to_dict()}")
    ab = reg[reg.location_quality.isin(["A", "B"])]
    F = pd.read_csv(f"{PROD}/setup/families.csv")
    F = F[F.region == R["famreg"]][["fid", "family_uid"]]
    T = ab.merge(F, on="family_uid", how="inner")
    log(f"A/B families {len(ab)}, in WA production {len(T)} (missing: {sorted(set(ab.family_uid) - set(T.family_uid))})")
    T = T[["fid", "family_uid", "family_id", "location_quality", "lat", "lon", "depth_km"]].rename(columns={"fid": "fam"})
    T["depth_cat"] = T.depth_km
    out = (T.depth_km < 25) | (T.depth_km > 45)
    T["slab2_km"] = slab2_depth(T.lat.to_numpy(), T.lon.to_numpy())
    T["depth_replaced"] = out
    T.loc[out, "depth_km"] = T.loc[out, "slab2_km"]
    log(f"depth outside 25-45 km replaced by Slab2: {int(out.sum())} families (catalogue {T.depth_cat[out].round(1).tolist()} -> "
        f"Slab2 {T.depth_km[out].round(1).tolist()}); still outside 25-45 after replacement: "
        f"{int(((T.depth_km < 25) | (T.depth_km > 45)).sum())}")
    return T.sort_values("fam").reset_index(drop=True)


def templates():
    L = pd.read_csv(f"{PROD}/setup/tpl_library.csv")[["fid", "netsta", "node", "n_events", "s_off", "tag"]]
    Q = pd.read_csv(f"{PROD}/templates/tpl_qc_w4b.csv")
    Q = Q[Q.use.astype(bool)][["fid", "netsta", "node", "n_events", "s_off"]].assign(tag="w4b")
    T = pd.concat([L, Q], ignore_index=True)
    return T[(T.n_events >= MIN_EVENTS) & T.s_off.notna()]


def vbar(depth, prof, col):
    """vertical-average slowness 0..depth (s/km)"""
    t = 0.0
    for r in prof.itertuples():
        top, bot = float(r.ztop_km), float(r.zbot_km)
        if top >= depth:
            break
        t += (min(bot, depth) - top) / float(getattr(r, col))
    return t / depth


def fdsn_depths(stas, log):
    rows = []
    for ns in stas:
        net, sta = ns.split(".")
        url = f"https://service.earthscope.org/fdsnws/station/1/query?net={net}&sta={sta}&level=channel&format=text"
        try:
            txt = urllib.request.urlopen(url, timeout=60).read().decode()
            D = pd.read_csv(io.StringIO(txt), sep="|")
            D.columns = [c.strip().lstrip("#").strip() for c in D.columns]
            D = D[D.Channel.str[:2].isin(["HH", "EH", "BH"])]
            rows.append(dict(Network=net, Station=sta, lat=float(D.Latitude.median()), lon=float(D.Longitude.median()),
                             elev_m=float(D.Elevation.median()), sensor_depth_m=float(D.Depth.max()), depth_min_m=float(D.Depth.min()),
                             channels=";".join(sorted(D.Channel.unique())), sensor=" | ".join(sorted(D.SensorDescription.astype(str).unique()))[:300],
                             source=f"EarthScope FDSN station service, level=channel, queried {C.pt_now()}"))
        except Exception as e:
            log(f"FDSN {ns}: {e}")
            rows.append(dict(Network=net, Station=sta, sensor_depth_m=np.nan, source=f"FDSN failed: {e}"))
    return pd.DataFrame(rows)


def epochs_table(stas):
    roots = [f"{C.ROOT}/data/extension_stations/{r}_bridge_nodes_strict.csv" for r in ("NW_WA", "SW_WA", "OLY2010")]
    N = pd.concat([pd.read_csv(p) for p in roots] + [pd.read_csv(f"{WA}/results/step4b/setup/extra_nodes.csv")])
    N["netsta"] = N.net + "." + N.sta
    E = pd.read_csv(f"{WA}/results/wa_pre/w2/epoch_changes.csv")
    rows = []
    for ns in stas:
        b = sorted(set(N[N.netsta == ns].start.astype(str)) | set(E[E.netsta == ns].date.astype(str)))
        why = {d: [] for d in b}
        for r in N[N.netsta == ns].itertuples():
            why[str(r.start)].append(f"node {r.node}")
        for r in E[E.netsta == ns].itertuples():
            why[str(r.date)].append(str(r.change)[:80])
        if not b:
            b = ["1990-01-01"]; why = {b[0]: ["no node table entry"]}
        for i, d in enumerate(b):
            rows.append(dict(netsta=ns, epoch=f"E{i + 1}", start=d if i else "1990-01-01", end=b[i + 1] if i + 1 < len(b) else "2100-01-01",
                             boundary=d, why="; ".join(why[d])))
    return pd.DataFrame(rows)


def epoch_of(EP, ns, dates):
    g = EP[EP.netsta == ns]
    d = pd.Series(np.asarray(dates, dtype=str)); out = np.full(len(d), "", dtype=object)
    for r in g.itertuples():
        out[((d >= r.start) & (d < r.end)).to_numpy()] = r.epoch
    return out


def episode_units(dates, log):
    E = pd.read_csv(f"{WA}/results/episode_rates/episodes_nw_final.csv")
    E = E[E.major.astype(bool)].copy()
    E["cls"] = "major"; E["big"] = True; E["src"] = "episodes_nw_final"
    E["w0"] = (pd.to_datetime(E.start) - pd.Timedelta(days=2)).dt.strftime("%Y-%m-%d")
    E["w1"] = (pd.to_datetime(E.end) + pd.Timedelta(days=2)).dt.strftime("%Y-%m-%d")
    E = E[["episode", "start", "end", "cls", "big", "src", "w0", "w1"]].sort_values("start").reset_index(drop=True)
    d = pd.Series(sorted(set(dates)))
    lab = pd.Series("", index=d.values, dtype=object)
    for r in E.itertuples():
        m = (d >= r.w0) & (d <= r.w1)
        lab[d[m].values] = np.where(lab[d[m].values] == "", r.episode, lab[d[m].values])
    rest = sorted(lab.index[lab == ""]); bursts = []; prev = None
    for x in rest:
        tx = pd.Timestamp(x)
        if prev is None or (tx - prev).days > 7:
            bursts.append(dict(episode=f"B{x}", start=x, end=x, cls="burst", big=False, src="bursts", w0=x, w1=x))
        bursts[-1]["end"] = x; bursts[-1]["w1"] = x; lab[x] = bursts[-1]["episode"]; prev = tx
    U = pd.concat([E, pd.DataFrame(bursts)], ignore_index=True).sort_values("start").reset_index(drop=True)
    log(f"episode units: majors {len(E)} ({len(set(lab) & set(E.episode))} with active days), bursts {len(bursts)}")
    return lab, U


def main():
    os.makedirs(C.MTB, exist_ok=True); os.makedirs(C.LOGD, exist_ok=True)
    log = C.Log(f"{C.LOGD}/wa_inventory.log"); t0 = time.time()
    log(f"region {C.REGION} -> {C.MTB}")
    Fm = families(log)
    fams = set(Fm.fam)
    # ---------------- repeats
    filt = [("region", "==", R["cat"])]
    filt += [("selected_v2", "==", True)] if C.REGION == "WA_NW" else [("P_own_pnet", ">=", 0.75), ("day_class", "==", "tremor")]
    cols = ["det_id", "time", "date", "family", "family_uid", "P_own_pnet", "role", "day_class", "selected_v2"]
    D = pq.read_table(f"{PROD}/post_v2/catalogue_v2.parquet", columns=cols, filters=filt).to_pandas()
    log(f"catalogue rows ({filt}): {len(D)}; of A/B families {int(D.family.isin(fams).sum())}; roles "
        f"{D[D.family.isin(fams)].role.value_counts().to_dict()}")
    D = D[D.family.isin(fams)].copy()
    D["t"] = D.time.astype("datetime64[ns]").astype("int64") / 1e9
    D = D.rename(columns={"family": "fam", "P_own_pnet": "P"})
    sp = D.groupby(["fam", "date"]).P.sum()
    act = sp[sp >= 3].rename("sumP").reset_index()
    key = set(zip(act.fam, act.date))
    on = np.array([(f, d) in key for f, d in zip(D.fam, D.date)])
    log(f"family-days {len(sp)}, active (sum P >= 3) {len(act)} on {act.date.nunique()} dates; rows on active days "
        f"{int(on.sum())} of {len(D)} ({100 * (1 - on.mean()):.1f} % not used)")
    D = D[on].reset_index(drop=True)
    act["n_data_rows"] = D.groupby(["fam", "date"]).size().reindex(pd.MultiIndex.from_frame(act[["fam", "date"]])).to_numpy()
    # ---------------- scrambled-control c1 rows (production detector, core network)
    rng = np.random.default_rng(SEED)
    nd = act.set_index(["fam", "date"]).n_data_rows
    c1 = []; chk = []
    for d in sorted(act.date.unique()):
        p = f"{PROD}/mf/prod/det_{d}.parquet"
        if not os.path.exists(p):
            continue
        M = pq.read_table(p, columns=["net", "series", "fid", "t", "k"], filters=[("net", "==", "core")]).to_pandas()
        q = M[(M.series == "c1") & M.fid.isin(fams)]
        for f, g in q.groupby("fid"):
            if (f, d) not in nd.index:
                continue
            n = int(nd[(f, d)])
            g = g.iloc[np.sort(rng.choice(len(g), n, replace=False))] if len(g) > n else g
            c1.append(pd.DataFrame(dict(fam=int(f), date=d, t=g.t.to_numpy(), k=g.k.to_numpy())))
        if len(chk) < 3:
            dd = M[M.series == "data"].set_index("k").t
            x = D[D.date == d]
            kk = (x.det_id % 10 ** 7).to_numpy()
            ok = np.isin(kk, dd.index)
            if ok.any():
                chk.append(float(np.max(np.abs(dd.reindex(kk[ok]).to_numpy() - x.t.to_numpy()[ok]))))
    C1 = pd.concat(c1, ignore_index=True) if c1 else pd.DataFrame(columns=["fam", "date", "t", "k"])
    log(f"c1 control rows kept {len(C1)} (<= n data per family-day); catalogue-vs-mf time check max |dt| {chk} s")
    X = pd.concat([D[["fam", "date", "t", "P", "det_id", "role", "day_class"]].assign(series="data"),
                   C1.assign(P=1.0, det_id=-1, role="", day_class="", series="c1")[["fam", "date", "t", "P", "det_id", "role", "day_class", "series"]]],
                  ignore_index=True)
    X["run"] = "wa_prod"; X["k"] = np.arange(len(X)); X["detector"] = "wa_prod"; X["src"] = "wa"
    X["row_id"] = np.arange(len(X)); X["det_id_cat"] = X.det_id; X["det_id"] = -1   # no cut store: raw reads only
    X["pilot_row_id"] = -1
    # ---------------- stations
    T = templates()
    T = T[T.fid.isin(fams)]
    best = T.sort_values(["fid", "netsta", "n_events"], ascending=[True, True, False]).groupby(["fid", "netsta"])
    O = best.first().reset_index()
    O["s_off_spread"] = best.s_off.agg(lambda x: float(x.max() - x.min())).to_numpy()
    O["n_nodes"] = best.size().to_numpy()
    ST = pd.read_csv(f"{PROD}/setup/stations.csv")
    O = O.merge(ST[["netsta", "lat", "lon"]], on="netsta", how="inner").merge(Fm[["fam", "lat", "lon", "depth_km"]].rename(
        columns={"fam": "fid", "lat": "flat", "lon": "flon"}), on="fid")
    O["epi_km"] = km(O.lat, O.lon, O.flat, O.flon); O["hypo_km"] = np.hypot(O.epi_km, O.depth_km)
    log(f"template (family, station) pairs {len(O)}; within {HYPO_MAX:.0f} km hypocentral {int((O.hypo_km <= HYPO_MAX).sum())}; "
        f"node s_off spread > 2 s: {int((O.s_off_spread > 2).sum())}")
    O = O[O.hypo_km <= HYPO_MAX].copy()
    prof = pd.read_csv(C.VEL)
    O["sp"] = [h * (vbar(z, prof, "vs") - vbar(z, prof, "vp")) for h, z in zip(O.hypo_km, O.depth_km)]
    off = O.rename(columns={"fid": "fam", "s_off": "offset"})[["fam", "netsta", "offset", "sp", "node", "n_events", "tag",
                                                                 "s_off_spread", "n_nodes", "epi_km", "hypo_km"]]
    off["offset_src"] = "wa_template_max_events"
    off.to_csv(f"{C.MTB}/offsets_sp.csv", index=False)
    stas = sorted(off.netsta.unique())
    # sites (single linkage < 3 km)
    S = ST[ST.netsta.isin(stas)].drop_duplicates("netsta").set_index("netsta").loc[stas].reset_index()
    lab = list(range(len(S)))
    for i in range(len(S)):
        for j in range(i + 1, len(S)):
            if km(S.lat[i], S.lon[i], S.lat[j], S.lon[j]) < SITE_KM:
                a, b = lab[i], lab[j]; lab = [a if x == b else x for x in lab]
    S["grp"] = lab
    S["site"] = S.groupby("grp").netsta.transform(lambda x: "+".join(sorted(s.split(".")[1] for s in x)))
    S["n_fams"] = S.netsta.map(off.groupby("netsta").fam.nunique())
    dep = fdsn_depths(stas, log)
    dep.to_csv(f"{C.TB}/sensor_depths.csv", index=False)
    S = S.merge(dep[["Network", "Station", "sensor_depth_m"]].assign(netsta=dep.Network + "." + dep.Station)[["netsta", "sensor_depth_m"]],
                on="netsta", how="left")
    S["sensor_depth_m"] = S.sensor_depth_m.fillna(0.0)
    EP = epochs_table(stas); EP.to_csv(f"{C.MTB}/station_epochs.csv", index=False)
    S["epochs"] = S.netsta.map(EP.groupby("netsta").apply(lambda g: " | ".join(f"{r.epoch} {r.boundary}" for r in g.itertuples()),
                                                         include_groups=False))
    S["use"] = True
    S.drop(columns="grp").to_csv(f"{C.MTB}/stations.csv", index=False)
    log(f"stations {len(S)} in {S.site.nunique()} sites; boreholes (sensor depth > 50 m): "
        f"{S[S.sensor_depth_m > 50].netsta.tolist()}; multi-station sites {S[S.site.str.contains('[+]')].site.unique().tolist()}")
    # ---------------- episode units
    lab, U = episode_units(act.date, log)
    act["episode"] = act.date.map(lab)
    act = act.merge(U[["episode", "cls"]], on="episode", how="left")
    U.to_csv(f"{C.MTB}/episode_units.csv", index=False); act.to_csv(f"{C.MTB}/active_days.csv", index=False)
    Fm.to_csv(f"{C.MTB}/families.csv", index=False)
    # ---------------- row x station expected S
    bad = pd.read_csv(f"{WA}/results/wa_pre/w2/r2_bad_for_dvv.csv")
    bad = bad[bad.reason == "dead"]          # PLAN §5 amendment: "flat" (>= 1 % repeated samples) flags whole years of GNW/HDW
    badset = set()
    for r in bad.itertuples():
        for dd in pd.date_range(r.start, r.end, freq="D"):
            badset.add((r.netsta, dd.strftime("%Y-%m-%d")))
    offk = off.set_index(["fam", "netsta"]).offset
    rows = []; nbad = 0
    for ns in stas:
        fs_ = set(off.fam[off.netsta == ns])
        Xs = X[X.fam.isin(fs_)]
        if not len(Xs):
            continue
        ep = epoch_of(EP, ns, Xs.date)
        o = np.array([offk[(f, ns)] for f in Xs.fam])
        good = np.array([(ns, d) not in badset for d in Xs.date])
        nbad += int((~good).sum())
        k = (ep != "") & good
        rows.append(pd.DataFrame(dict(row_id=Xs.row_id.to_numpy()[k], netsta=ns, epoch=ep[k], tS=(Xs.t.to_numpy() + o)[k],
                                      date=Xs.date.to_numpy()[k])))
    RS = pd.concat(rows, ignore_index=True)
    X.to_parquet(f"{C.MTB}/rows.parquet", index=False); RS.to_parquet(f"{C.MTB}/row_station.parquet", index=False)
    inv = RS.merge(X[["row_id", "series"]], on="row_id").groupby(["netsta", "series"]).size().unstack(fill_value=0)
    sd = RS.groupby("netsta").date.nunique()
    log(f"row-stations {len(RS)} (dead/flat station-days dropped: {nbad}); station-days to read {int(sd.sum())}\n" + inv.to_string())
    json.dump(dict(region=C.REGION, families=len(Fm), rows_data=int((X.series == "data").sum()), rows_c1=int((X.series == "c1").sum()),
                   active_family_days=len(act), dates=int(act.date.nunique()), stations=len(S), sites=int(S.site.nunique()),
                   row_stations=len(RS), station_days=int(sd.sum()), written=C.pt_now()),
              open(f"{C.MTB}/inventory_summary.json", "w"), indent=1)
    open(f"{C.OUT}/WA_INVENTORY_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log(f"inventory done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
