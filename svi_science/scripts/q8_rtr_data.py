"""Q8 RTR: unified event loaders. Each returns (events DataFrame, episodes DataFrame).

events: ep, time (UTC ns), x (along-strike km, + NW), dd (down-dip km), unit, mag (Bostock) , fam
episodes: ep, w0, w1 (window, UTC), onset (episode start date), n
"""
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
R = f"{ROOT}/svi_science/results"


def _ns(s):
    return pd.to_datetime(s).astype("datetime64[ns]")


def bostock():
    B = pd.read_parquet(f"{R}/p3_bostock/bostock_events_mapped.parquet")
    E = pd.read_csv(f"{R}/p3_bostock/episodes_bostock_final.csv")
    E["w0"] = _ns(E.start) - pd.Timedelta(days=3)
    E["w1"] = _ns(E.end) + pd.Timedelta(days=4)
    E["onset"] = _ns(E.start)
    B["time"] = B.time_utc.astype("datetime64[ns]")
    B["ep"] = None
    for r in E.itertuples():
        m = (B.time >= r.w0) & (B.time < r.w1)
        B.loc[m, "ep"] = r.episode
    B = B[B.ep.notna()].copy()
    ev = pd.DataFrame(dict(ep=B.ep.values, time=B.time.values, x=B.along_km.values, dd=B.downdip_km.values,
                           unit=B.fid.values, fam=B.fid.values, mag=B.mag.values, dep=B.dep.values))
    ev = ev.sort_values(["ep", "time"]).reset_index(drop=True)
    E = E.rename(columns={"episode": "ep"})
    E["n"] = E.ep.map(ev.groupby("ep").size())
    return ev, E[["ep", "w0", "w1", "onset", "n"]]


def _eps(region):
    E = pd.read_csv(f"{R}/s3/episodes_all.csv")
    E = E[(E.region == region) & E.big].rename(columns={"episode": "ep"})
    E["w0"] = _ns(E.w0); E["w1"] = _ns(E.w1); E["onset"] = _ns(E.start)
    return E


def ours_svi(variant="head"):
    T = pd.read_parquet(f"{R}/s3/events_svi.parquet")
    if variant == "head":
        T = T[(T.series == "data") & (T.P_same >= 0.75) & T.counted_once]
    elif variant == "p050":
        T = T[(T.series == "data") & (T.P_same >= 0.5) & T.counted_once]
    elif variant == "raw075":
        T = T[(T.series == "data") & (T.P_same >= 0.75)]
    elif variant == "scr075":       # scrambled-control catalogue (stations shifted to different days)
        T = T[(T.series == "c1") & (T.P_same >= 0.75) & T.counted_once]
    E = _eps("SVI")
    T = T[T.ep.isin(E.ep)]
    ev = pd.DataFrame(dict(ep=T.ep.values, time=T.time.values.astype("datetime64[ns]"), x=T.as_km.values, dd=T.dd_km.values,
                           unit=T.unit.values, fam=T.unit.values, mag=np.nan, dep=T.depth_km.values, P=T.P_same.values))
    ev = ev.sort_values(["ep", "time"]).reset_index(drop=True)
    E["n"] = E.ep.map(ev.groupby("ep").size()).fillna(0).astype(int)
    return ev, E[["ep", "w0", "w1", "onset", "n"]]


def ours_nw(variant="head"):
    T = pd.read_parquet(f"{R}/s3/events_nw.parquet")
    if variant == "head":
        T = T[T.selected_v2 & T.counted_once]
    elif variant == "p050":
        T = T[(T.P_same >= 0.5) & T.counted_once]
    elif variant == "raw075":
        T = T[T.selected_v2]
    T = T[T.as_km.notna()]
    E = _eps("NW")
    T = T[T.ep.isin(E.ep)]
    ev = pd.DataFrame(dict(ep=T.ep.values, time=T.time.values.astype("datetime64[ns]"), x=T.as_km.values, dd=T.dd_km.values,
                           unit=T.unit.values, fam=T.unit_fid.values, mag=np.nan, dep=T.depth_km.values, P=T.P_same.values))
    ev = ev.sort_values(["ep", "time"]).reset_index(drop=True)
    E["n"] = E.ep.map(ev.groupby("ep").size()).fillna(0).astype(int)
    return ev, E[["ep", "w0", "w1", "onset", "n"]]


LOADERS = {"bostock": bostock, "svi": ours_svi, "nw": ours_nw}


# ------------------------------------------------------------------ Q8b: PNSN tremor (one event = one tremor window)
_PNSN = {}


def _pnsn_all():
    if "P" not in _PNSN:
        P = pd.read_csv(f"{ROOT}/catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
        P["time"] = pd.to_datetime(P.time).astype("datetime64[ns]")
        _PNSN["P"] = P
    return _PNSN["P"]


def pnsn(region):
    """region in SVI, NW, SW. Box = <= 10 km from the region's family/source positions; S3 frames."""
    import json
    import s3_build as s3
    P = _pnsn_all()
    if region == "SVI":
        F = pd.read_csv(f"{R}/s0/families.csv")
        frs = json.load(open(f"{R}/s0/frame_groups.json"))
        fr = dict(origin_lat=frs["origin_lat"], origin_lon=frs["origin_lon"], downdip_az=frs["mean_downdip_az"], strike_az=frs["strike_az"])
        plat, plon = F.lat.values, F.lon.values
    elif region == "NW":
        S = pd.read_csv(f"{R}/s3/sources_nw.csv"); S = S[S.quality.isin(["A", "B"])]
        fr = json.load(open(f"{R}/s3/build_info.json"))["nw_frame"]
        plat, plon = S.lat.values, S.lon.values
    else:
        Fm = pd.read_csv(f"{ROOT}/data/catalogue/families_v1_3.csv")
        Fm = Fm[(Fm.region == "Royer_SW_WA") & Fm.pos_lat.notna()]
        plat, plon = Fm.pos_lat.values, Fm.pos_lon.values
        fr = s3.frame_from(plat, plon)
    m = s3.near_mask(P.lat.values, P.lon.values, plat, plon, 10.0)
    Q = P[m].copy()
    a, d, _, _ = s3.project(fr, Q.lat.values, Q.lon.values)
    Q["x"] = a; Q["dd"] = d
    E = pd.read_csv(f"{R}/b_ets_record/episodes_b.csv")
    E = E[E.big & (E.region == ("SVI" if region == "SVI" else "NW"))].copy()
    E["w0"] = _ns(E.w0); E["w1"] = _ns(E.w1) + pd.Timedelta(days=1)
    E = E[E.w1 > Q.time.min()]
    E["onset"] = _ns(E.start)
    Q["ep"] = None
    for r in E.itertuples():
        mm = (Q.time >= r.w0) & (Q.time < r.w1)
        Q.loc[mm, "ep"] = r.episode
    Q = Q[Q.ep.notna()]
    # nearest SVI family for tides
    if region == "SVI":
        dx = (Q.lon.values[:, None] - plon[None, :]) * 111.195 * np.cos(np.radians(48.4))
        dy = (Q.lat.values[:, None] - plat[None, :]) * 111.195
        tf = F.fid.values[np.argmin(dx ** 2 + dy ** 2, 1)]
    else:
        tf = np.full(len(Q), -1)
    unit = np.floor(Q.x.values / 2.0).astype(int)
    cell = (np.floor(Q.x.values / 4.0).astype(int) * 1000 + np.floor(Q.dd.values / 4.0).astype(int))
    ev = pd.DataFrame(dict(ep=Q.ep.values, time=Q.time.values, x=Q.x.values, dd=Q.dd.values, unit=unit, fam=cell,
                           tide_fid=tf, mag=np.nan, lat=Q.lat.values, lon=Q.lon.values))
    ev = ev.sort_values(["ep", "time"]).reset_index(drop=True)
    E = E.rename(columns={"episode": "ep"})
    E["n"] = E.ep.map(ev.groupby("ep").size()).fillna(0).astype(int)
    E["s3_big"] = E.s3_big.astype(bool)
    return ev, E[["ep", "w0", "w1", "onset", "n", "s3_big"]]


# ------------------------------------------------------------------ Q8c: all candidates (no P cut)
_AC = {}


def _pack(T, E, x, dd, unit, fam, P):
    ev = pd.DataFrame(dict(ep=T.ep.values, time=T.time.values.astype("datetime64[ns]"), x=x, dd=dd, unit=unit, fam=fam,
                           mag=np.nan, P=P))
    ev = ev.sort_values(["ep", "time"]).reset_index(drop=True)
    E = E.copy(); E["n"] = E.ep.map(ev.groupby("ep").size()).fillna(0).astype(int)
    return ev, E[["ep", "w0", "w1", "onset", "n"]]


def _assign(t, E):
    tv = np.asarray(t, dtype="datetime64[ns]")
    ep = np.full(len(tv), None, dtype=object)
    for r in E.itertuples():
        m = (tv >= np.datetime64(r.w0)) & (tv < np.datetime64(r.w1))
        ep[m] = r.ep
    return ep


def allcand_svi(variant="raw", series="data"):
    key = ("svi", series)
    if key not in _AC:
        T = pd.read_parquet(f"{R}/s0/analysis_table.parquet", columns=["series", "fam", "P_same", "keep_xfam8", "time_utc"])
        _AC[key] = T[T.series == series].copy()
    T = _AC[key]
    if variant == "once":
        T = T[T.keep_xfam8.astype(bool)]
    E = _eps("SVI")
    T = T.assign(time=T.time_utc.astype("datetime64[ns]"))
    T = T.assign(ep=_assign(T.time, E))
    T = T[T.ep.notna()]
    F = pd.read_csv(f"{R}/s0/families.csv").set_index("fid")
    return _pack(T, E, T.fam.map(F.along_strike_km).values, T.fam.map(F.downdip_km).values, T.fam.values, T.fam.values,
                 T.P_same.values)


def _wa_flags():
    if "wa" not in _AC:
        import glob
        fl = []
        for f in sorted(glob.glob(f"{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/dedup/flags/flags_20*.parquet")):
            if f.split("_")[-1][:7] < "2015-01":
                continue
            q = pd.read_parquet(f, columns=["det_id", "series", "fid", "t", "dup_D"])
            fl.append(q[q.series.isin(["data", "c1"])])
        fl = pd.concat(fl, ignore_index=True)
        fl["time"] = pd.to_datetime(fl.t, unit="s").astype("datetime64[ns]")
        C = pd.read_parquet(f"{ROOT}/lfe_extension_wa/results/wa_prod/post_v2/catalogue_v2.parquet", columns=["det_id", "P_own_pnet"],
                            filters=[("region", "in", ["NW", "SW"])])
        fl = fl.merge(C, on="det_id", how="left")
        # Q8d: frozen-model P for the scrambled-control rows (new file, written by q8d_wa_scrambled_P.py)
        import os
        fp = f"{R}/q8_rtr/tables/wa_scrambled_P_2015_2026.parquet"
        if os.path.exists(fp):
            Pc = pd.read_parquet(fp, columns=["det_id", "P_own_WAVE_PNET"]).rename(columns={"P_own_WAVE_PNET": "P_c1"})
            fl = fl.merge(Pc, on="det_id", how="left")
            fl["P_own_pnet"] = np.where(fl.series == "c1", fl.P_c1, fl.P_own_pnet)
            fl = fl.drop(columns="P_c1")
        _AC["wa"] = fl
    return _AC["wa"]


def allcand_wa(region="NW", variant="raw", series="data"):
    import json
    import s3_build as s3
    fl = _wa_flags()
    fl = fl[fl.series == series]
    if variant == "once":
        fl = fl[~fl.dup_D.fillna(False).astype(bool)]
    E = _eps("NW")
    setup = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/wa_prod/setup/families.csv")
    if region == "NW":
        SG = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/nw_sharing/source_groups_E3.csv")
        AC = pd.read_csv(f"{ROOT}/lfe_extension_wa/results/nw_sharing/source_groups_E3_array_confirmed.csv").rename(columns={"group": "source"})
        f2s = SG.merge(AC, on="fid").set_index("fid").source
        Sx = pd.read_csv(f"{R}/s3/sources_nw.csv").set_index("source")
        fl = fl[fl.fid.isin(setup[setup.reg == "NW"].fid) & fl.fid.isin(f2s.index)]
        src = fl.fid.map(f2s)
        x = src.map(Sx.as_km).values; dd = src.map(Sx.dd_km).values; unit = src.values
    else:
        sw = setup[setup.reg == "SW"].set_index("fid")
        fr = s3.frame_from(sw.lat_use.values, sw.lon_use.values)
        a, d, _, _ = s3.project(fr, sw.lat_use.values, sw.lon_use.values)
        sw["x"] = a; sw["dd"] = d
        fl = fl[fl.fid.isin(sw.index)]
        x = fl.fid.map(sw.x).values; dd = fl.fid.map(sw.dd).values; unit = fl.fid.values
    T = fl.assign(x=x, dd=dd, unit=unit)
    T = T.assign(ep=_assign(T.time, E))
    T = T[T.ep.notna() & np.isfinite(T.x)]
    return _pack(T, E, T.x.values, T.dd.values, T.unit.values, T.fid.values, T.P_own_pnet.values)


def back_2010(variant="raw", series="data"):
    """Backward production forward-detector candidates in the Bostock windows 2010-08, 2011-08, 2012-08, 2013-09."""
    import glob
    key = ("back", series)
    _, EB = bostock()
    EB = EB[EB.ep.isin(["2010-08", "2011-08", "2012-08", "2013-09"])]
    if key not in _AC:
        out = []
        for r in EB.itertuples():
            for d in pd.date_range(r.w0, r.w1 - pd.Timedelta(days=1), freq="D"):
                f = f"{R}/backward/mf/feat/p5f_{d.year}/{d.date()}.parquet"
                try:
                    q = pd.read_parquet(f, columns=["series", "fam", "t", "keep_xfam8"])
                except Exception:
                    continue
                out.append(q)
        T = pd.concat(out, ignore_index=True)
        _AC[key] = T[T.series == series].copy()
    T = _AC[key]
    if variant == "once":
        T = T[T.keep_xfam8.astype(bool)]
    T = T.assign(time=pd.to_datetime(T.t, unit="s").astype("datetime64[ns]"))
    T = T.assign(ep=_assign(T.time, EB))
    T = T[T.ep.notna()]
    F = pd.read_csv(f"{R}/s0/families.csv").set_index("fid")
    T = T[T.fam.isin(F.index)]
    return _pack(T, EB, T.fam.map(F.along_strike_km).values, T.fam.map(F.downdip_km).values, T.fam.values, T.fam.values,
                 np.full(len(T), np.nan))


# ------------------------------------------------------------------ Q8e: WA-NW diagnostic variants
Q8E_STRIKE_A1 = 345.24
Q8E_SINGLE_FRONT = None


def q8e_nw(variant, series="data"):
    """variant in V0, Pa, Pb, Pc, F1, F2, F3, A1 (PLAN.md Q8e). Base = allcand_wa('NW', 'raw', series)."""
    import json
    ev, E = allcand_wa("NW", "once" if variant in ("Pd", "Pe") else "raw", series)
    Sx = pd.read_csv(f"{R}/s3/sources_nw.csv").set_index("source")
    if variant == "Pe":
        ev = ev[ev.unit.isin(Sx[Sx.quality == "A"].index)]
    if variant in ("Pa", "Pc"):
        ok = set(Sx.best_fid.astype(int)) if variant == "Pa" else set(Sx[Sx.quality == "A"].best_fid.astype(int))
        ev = ev[ev.fam.isin(ok)]
    if variant == "Pb":
        ev = ev[ev.unit.isin(Sx[Sx.quality == "A"].index)]
    if variant in ("F2", "F3"):
        K = pd.read_csv(f"{R}/s3/episodes_kinematics.csv")
        keep = K[(K.region == "NW") & K.front_class.str.startswith("one front")].episode
        ev = ev[ev.ep.isin(keep)]; E = E[E.ep.isin(keep)]
    if variant == "A1":
        fr = json.load(open(f"{R}/s3/build_info.json"))["nw_frame"]
        s = np.radians(Q8E_STRIKE_A1); d = np.radians(Q8E_STRIKE_A1 + 90)
        xx = (Sx.lon - fr["origin_lon"]) * 111.195 * np.cos(np.radians(fr["origin_lat"])); yy = (Sx.lat - fr["origin_lat"]) * 111.195
        ev = ev.assign(x=ev.unit.map(xx * np.sin(s) + yy * np.cos(s)).values, dd=ev.unit.map(xx * np.sin(d) + yy * np.cos(d)).values)
    ev = ev.sort_values(["ep", "time"]).reset_index(drop=True)
    E = E.copy(); E["n"] = E.ep.map(ev.groupby("ep").size()).fillna(0).astype(int)
    return ev, E


def q8e_tremor_overrides():
    """per NW episode: (tremor onset Series indexed by 4-km bin, hours rel. to the S3 window start), tremor main direction."""
    import q8_rtr_backzone as BZ
    Tr = pd.read_parquet(f"{R}/s3/tremor_nw.parquet")
    K = pd.read_csv(f"{R}/s3/episodes_kinematics.csv"); K = K[K.region == "NW"].set_index("episode")
    E = _eps("NW").set_index("ep")
    out = {}
    for ep, e in E.iterrows():
        g = Tr[Tr.ep == ep]
        t = ((pd.to_datetime(g.time) - e.w0).dt.total_seconds() / 3600).values
        out[ep] = (BZ.onset_map(t, g.as_km.values), float(np.sign(K.loc[ep, "trem_speed"])))
    return out
