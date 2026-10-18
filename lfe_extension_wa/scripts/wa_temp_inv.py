#!/usr/bin/env python
"""Inventory of TEMPORARY seismic networks 2012-01-01 .. 2026-09-30 within 75 km of a Washington family
(v1.1 representatives of Royer NW/SW and Olympic arm A/B in data/catalogue/families_v1_2.csv; position pos_lat/lon,
else approx_lat/lon).  Read-only on the providers: FDSN station service (EarthScope, NCEDC, NRCan) + one
10-minute waveform sample per station-year (the BRIDGE_PLAN.md convention; middle of the station's epoch inside
that year, vertical channel first).  Network Z5 is excluded by the user (restricted); it is listed only.

  python scripts/wa_temp_inv.py stations   -> results/wa_temp_networks/channels.csv, stations.csv
  python scripts/wa_temp_inv.py samples    -> results/wa_temp_networks/samples.csv (resumable)
  python scripts/wa_temp_inv.py report     -> results/wa_temp_networks/inventory.csv (+ by_network.csv)
"""
from __future__ import annotations
import os, sys, time, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/lfe_extension_wa/results/wa_temp_networks"
NRCAN = "https://earthquakescanada.nrcan.gc.ca"
PROVS = ["EARTHSCOPE", "NCEDC", NRCAN]
T0, T1 = "2012-01-01", "2026-09-30"
BOX = dict(minlatitude=45.7, maxlatitude=49.7, minlongitude=-125.9, maxlongitude=-120.9)
EXCLUDED = {"Z5": "excluded, restricted (user handles it separately)"}
NONSEIS = {"1H": "magnetotelluric (MOCHA), no seismic channels", "C2": "geomagnetic observatory"}
SEIS_BANDS = set("BHESCDML")          # band codes kept; instrument codes H (seismometer), N (accel), P/L/G noted
os.makedirs(OUT, exist_ok=True)


def is_temp(code):
    return code[0] in "XYZ0123456789" or code == "TA"


def families():
    F = pd.read_csv(f"{ROOT}/data/catalogue/families_v1_2.csv", low_memory=False)
    F = F[F.region.isin(["Royer_NW_WA", "Royer_SW_WA", "Olympic_2010", "Olympic_2010_armB"])
          & (F.representative_v1_1 == True)].copy()
    F["flat"] = F.pos_lat.fillna(F.approx_lat)
    F["flon"] = F.pos_lon.fillna(F.approx_lon)
    return F[F.flat.notna()].reset_index(drop=True)


def dist_km(F, la, lo):
    return np.hypot((F.flon.values - lo) * 111.195 * np.cos(np.radians((F.flat.values + la) / 2)),
                    (F.flat.values - la) * 111.195)


def stations():
    from obspy.clients.fdsn import Client
    from obspy import UTCDateTime as U
    F = families()
    rows = []
    for prov in PROVS:
        c = Client(prov, timeout=300)
        nets = c.get_stations(starttime=U(T0), endtime=U(T1), level="network", **BOX)
        for n in nets:
            if not is_temp(n.code) and not (prov == NRCAN and n.code not in ("CN", "PO", "C8")):
                continue
            try:
                inv = c.get_stations(network=n.code, starttime=U(T0), endtime=U(T1), level="channel", **BOX)
            except Exception as e:
                print(prov, n.code, "channel query failed:", e)
                rows.append(dict(provider=prov, net=n.code, sta="", note=f"channel query failed {e}"[:200],
                                 net_restricted=n.restricted_status))
                continue
            for nn in inv:
                for s in nn:
                    for ch in s:
                        rows.append(dict(
                            provider=prov, net=nn.code, net_start=str(nn.start_date)[:10] if nn.start_date else "",
                            net_end=str(nn.end_date)[:10] if nn.end_date else "",
                            net_restricted=nn.restricted_status, net_desc=(nn.description or "")[:80],
                            sta=s.code, sta_restricted=s.restricted_status, lat=s.latitude, lon=s.longitude,
                            elev=s.elevation, loc=ch.location_code, cha=ch.code, fs=ch.sample_rate,
                            start=str(ch.start_date)[:10], end=str(ch.end_date)[:10] if ch.end_date else "2599-12-31",
                            cha_restricted=ch.restricted_status, dip=ch.dip,
                            sensor=(ch.sensor.description if ch.sensor and ch.sensor.description else
                                    (ch.sensor.type if ch.sensor else "")) or "",
                            depth=ch.depth))
            print(prov, n.code, len(rows), flush=True)
    C = pd.DataFrame(rows)
    C = C[C.sta != ""].copy()
    C = C.drop_duplicates(["net", "sta", "loc", "cha", "start"])            # same epoch served by two providers
    d = {k: dist_km(F, la, lo) for k, la, lo in
         C[["net", "sta", "lat", "lon"]].drop_duplicates(["net", "sta"]).itertuples(index=False, name=None)
         and [((a, b), c, e) for a, b, c, e in C[["net", "sta", "lat", "lon"]].drop_duplicates(["net", "sta"])
              .itertuples(index=False, name=None)]}
    near = {}
    for k, dd in d.items():
        o = np.argsort(dd)
        near[k] = (float(dd[o[0]]), F.family_uid.values[o[0]], F.region.values[o[0]],
                   int((dd <= 10).sum()), int((dd <= 30).sum()), int((dd <= 75).sum()))
    C["dmin_km"] = [near[(a, b)][0] for a, b in zip(C.net, C.sta)]
    C = C[C.dmin_km <= 75].copy()
    C.to_csv(f"{OUT}/channels.csv", index=False)
    g = C.groupby(["net", "sta"])
    S = g.agg(provider=("provider", "first"), lat=("lat", "first"), lon=("lon", "first"), elev=("elev", "first"),
              start=("start", "min"), end=("end", "max"), dmin_km=("dmin_km", "first"),
              restricted=("cha_restricted", lambda v: ",".join(sorted({str(x) for x in v}))),
              net_desc=("net_desc", "first")).reset_index()
    S["start"] = S.start.clip(lower=T0)
    S["end"] = S.end.clip(upper=T1)
    S["nearest_family"] = [near[(a, b)][1] for a, b in zip(S.net, S.sta)]
    S["nearest_region"] = [near[(a, b)][2] for a, b in zip(S.net, S.sta)]
    S["n_fam_10km"] = [near[(a, b)][3] for a, b in zip(S.net, S.sta)]
    S["n_fam_30km"] = [near[(a, b)][4] for a, b in zip(S.net, S.sta)]
    S["n_fam_75km"] = [near[(a, b)][5] for a, b in zip(S.net, S.sta)]
    S["seafloor"] = S.elev < -10
    S.to_csv(f"{OUT}/stations.csv", index=False)
    print(S.groupby("net").size().to_string())


def _seis(C):
    b = C.cha.str[0]
    i = C.cha.str[1]
    return C[b.isin(list(SEIS_BANDS)) & i.isin(["H", "N", "P", "L", "G", "X", "D"])]


def best_channels(Cs):
    """per station epoch: the seismometer channel group to sample/download (highest-rate H instrument, else N)."""
    Cs = Cs.copy()
    Cs["inst"] = Cs.cha.str[1]
    Cs["pref"] = Cs.inst.map({"H": 0, "L": 1, "P": 1, "N": 2}).fillna(3)
    return Cs.sort_values(["pref", "fs"], ascending=[True, False])


_CL = {}


def client(prov):
    from obspy.clients.fdsn import Client
    k = (prov, threading.get_ident())
    if k not in _CL:
        _CL[k] = Client(prov, timeout=120)
    return _CL[k]


def sample_one(job):
    """one 10-min sample at the first of the times (';'-separated) that returns data"""
    prov, net, sta, loc, cha, ts = job
    r = None
    tl = ts.split(";")
    with _NLK:
        st = _NET.setdefault((net, cha[:2]), [0, 0])
        if st[0] == 0 and st[1] >= 10:          # network-group with >= 10 empty stations and none with data:
            tl = tl[:1]                        # one try only (the fallback times are for patchy recorders)
    for i, t in enumerate(tl):
        r = _sample_at(prov, net, sta, loc, cha, t)
        r["t_used"], r["n_tries"] = t, i + 1
        if r["status"] in ("data", "restricted"):
            break
    with _NLK:
        _NET[(net, cha[:2])][0 if r["status"] == "data" else 1] += 1
    return r


_NET, _NLK = {}, threading.Lock()


def _sample_at(prov, net, sta, loc, cha, t):
    from obspy import UTCDateTime as U
    last = ""
    for k in range(3):
        try:
            st = client(prov).get_waveforms(net, sta, loc if loc else "--", cha, U(t), U(t) + 600)
            npts = sum(tr.stats.npts for tr in st)
            fs = st[0].stats.sampling_rate if len(st) else 0
            import io
            bio = io.BytesIO()
            if len(st):
                st.write(bio, format="MSEED")
            return dict(status="data" if npts else "empty", npts=npts, fs_got=fs,
                        frac=round(npts / (600 * fs), 3) if fs else 0.0, bytes10=len(bio.getvalue()), note="")
        except Exception as e:
            last = f"{type(e).__name__}:{str(e)[:100]}".replace("\n", " ")
            if "No data" in last or "204" in last or "NoData" in last or "404" in last:
                return dict(status="nodata", npts=0, fs_got=0, frac=0.0, bytes10=0, note="")
            if "403" in last or "401" in last or "Unauthorized" in last or "Forbidden" in last:
                return dict(status="restricted", npts=0, fs_got=0, frac=0.0, bytes10=0, note=last)
            time.sleep(5 * (k + 1))
    return dict(status="error", npts=0, fs_got=0, frac=0.0, bytes10=0, note=last)


def samples():
    C = pd.read_csv(f"{OUT}/channels.csv", keep_default_na=False, dtype={"loc": str})
    S = pd.read_csv(f"{OUT}/stations.csv", keep_default_na=False)
    C = _seis(C)
    jobs = []
    for (net, sta), g in C.groupby(["net", "sta"]):
        if net in EXCLUDED or net in NONSEIS:
            continue
        g = best_channels(g)
        for y in range(2012, 2027):
            a, b = max(f"{y}-01-01", T0), min(f"{y}-12-31", T1)
            live = g[(g.start <= b) & (g.end > a)]
            if not len(live):
                continue
            r = live.iloc[0]
            # vertical (or first) channel of the preferred group
            same = live[(live.cha.str[:2] == r.cha[:2]) & (live["loc"] == r["loc"])]
            z = same[same.cha.str[2].isin(["Z", "3"])]
            ch = (z.iloc[0] if len(z) else same.iloc[0]).cha
            ea, eb = pd.Timestamp(max(a, r.start)), pd.Timestamp(min(b, r.end[:10]))
            fr = [0.5] + ([0.1, 0.3, 0.7, 0.9] if S.set_index(["net", "sta"]).dmin_km.get((net, sta), 99) <= 40
                          else [])
            tt = [(ea + (eb - ea) * f).normalize() + pd.Timedelta(hours=12) for f in fr]
            tt = [x if x <= eb + pd.Timedelta(hours=12) else ea + pd.Timedelta(hours=1) for x in tt]
            jobs.append((r.provider, net, sta, r["loc"], ch, y, ";".join(x.strftime("%Y-%m-%dT%H:%M:%S") for x in tt),
                         r.cha[:2], len(same.cha.unique()), r.fs))
    done = set()
    fp = f"{OUT}/samples.csv"
    if os.path.exists(fp):
        D = pd.read_csv(fp, keep_default_na=False)
        done = set(zip(D.net, D.sta, D.year.astype(int)))
    jobs = [j for j in jobs if (j[1], j[2], j[5]) not in done]
    print(len(jobs), "sample requests", flush=True)
    rows = []
    with ThreadPoolExecutor(4) as ex:
        fut = {ex.submit(sample_one, (j[0], j[1], j[2], j[3], j[4], j[6])): j for j in jobs}
        for i, f in enumerate(as_completed(fut)):
            j = fut[f]
            rows.append(dict(provider=j[0], net=j[1], sta=j[2], loc=j[3], cha=j[4], year=j[5], t=j[6], group=j[7],
                             n_cha_group=j[8], fs_meta=j[9], **f.result()))
            if len(rows) >= 50 or i == len(jobs) - 1:
                pd.DataFrame(rows).to_csv(fp, mode="a", header=not os.path.exists(fp), index=False)
                rows = []
                print(i + 1, "/", len(jobs), flush=True)


if __name__ == "__main__":
    if sys.argv[1] in ("stations", "samples"):
        {"stations": stations, "samples": samples}[sys.argv[1]]()


# ---------------------------------------------------------------------------------------------- report / selection
SEL_KM, SPARSE_N, MAX_KM = 30.0, 3, 75.0
CLOSED_NOTE = "excluded, restricted (closed at the provider; public data only)"


def _days(c):
    """set of UTC days in [T0, T1] covered by the channel epochs c"""
    out = set()
    for s, e in zip(c.start, c.end):
        a, b = max(pd.Timestamp(s[:10]), pd.Timestamp(T0)), min(pd.Timestamp(e[:10]), pd.Timestamp(T1))
        if b >= a:
            out |= set(pd.date_range(a, b).strftime("%Y-%m-%d"))
    return out


def report():
    C = pd.read_csv(f"{OUT}/channels.csv", keep_default_na=False, dtype={"loc": str})
    S = pd.read_csv(f"{OUT}/stations.csv", keep_default_na=False)
    D = pd.read_csv(f"{OUT}/samples.csv", keep_default_na=False, dtype={"loc": str})
    F = families()
    rows = []
    for r in S.itertuples():
        c = C[(C.net == r.net) & (C.sta == r.sta)]
        d = D[(D.net == r.net) & (D.sta == r.sta)]
        bands = ",".join(sorted(set(c.cha.str[:2])))
        restr = set(c.cha_restricted) | set(c.net_restricted)
        g = loc = sensor = ""
        ncomp, fs, mb, ndays = 0, 0.0, 0.0, 0
        if len(d):
            dd = d[d.status == "data"]
            pick = (dd if len(dd) else d).group.mode().iloc[0]
            g = pick
            loc = (dd if len(dd) else d)[lambda x: x.group == pick]["loc"].mode().iloc[0]
            cg = c[(c.cha.str[:2] == g) & (c["loc"] == loc)]
            ncomp = cg.cha.nunique()
            fs = float(cg.fs.max()) if len(cg) else 0.0
            sensor = cg.sensor.mode().iloc[0] if len(cg) else ""
            ndays = len(_days(cg))
            bz = dd[dd.group == g].bytes10
            mb = float(bz.median()) * 144 * ncomp / 1e6 if len(bz) else 0.0
        ys, yd = len(d), int((d.status == "data").sum())
        inst = g[1:2]
        if r.net in EXCLUDED:
            status = EXCLUDED[r.net]
        elif r.net in NONSEIS:
            status = "not seismic: " + NONSEIS[r.net]
        elif "closed" in restr or "partial" in restr:
            status = CLOSED_NOTE
        elif inst == "N":
            status = "accelerometer only (no seismometer channels)"
        elif inst in ("D", "X", "G") or (g == "" and not len(d)):
            status = "no seismometer channels (pressure/other only)" if g else "no seismometer channels"
        elif yd == 0:
            status = "metadata only: no data served (all samples empty)"
        else:
            status = "ok"
        dist = np.hypot((F.flon.values - r.lon) * 111.195 * np.cos(np.radians((F.flat.values + r.lat) / 2)),
                        (F.flat.values - r.lat) * 111.195)
        o = np.argsort(dist)[:3]
        rows.append(dict(net=r.net, sta=r.sta, provider=r.provider, lat=r.lat, lon=r.lon, elev=r.elev,
                         seafloor=r.elev < -10, start=max(r.start, T0), end=min(r.end[:10], T1), all_bands=bands,
                         group=g, loc=loc, n_comp=ncomp, three_c=ncomp >= 3, fs=fs, sensor=sensor,
                         live_days=ndays, sample_years=ys, sample_years_data=yd,
                         mb_day=round(mb, 2), est_gb=round(mb * ndays * (yd / ys if ys else 0) / 1e3, 2),
                         dmin_km=round(r.dmin_km, 2),
                         nearest3=";".join(f"{F.family_uid.values[i]}:{dist[i]:.1f}" for i in o),
                         nearest_region=F.region.values[o[0]], n_fam_10km=r.n_fam_10km,
                         n_fam_30km=r.n_fam_30km, n_fam_75km=r.n_fam_75km, status=status))
    I = pd.DataFrame(rows)
    # --- selection: ok AND <= 30 km; + nearest temporary station(s) of families with < 3 usable stations
    I["select"] = ""
    ok = I.status == "ok"
    I.loc[ok & (I.dmin_km <= SEL_KM), "select"] = "<=30 km"
    N = pd.concat([pd.read_csv(f"{ROOT}/data/extension_stations/{a}_bridge_nodes_{k}.csv", keep_default_na=False)
                   for a in ["NW_WA", "SW_WA", "OLY2010"] for k in ["strict", "lenient"]])
    N = N[(N.end > T0) & (N.start < T1)].drop_duplicates(["net", "sta"])
    N = N[~N.net.isin(list(EXCLUDED))]
    perm = N[~N.net.map(is_temp)]
    sparse = []
    for f in F.itertuples():
        def dk(la, lo):
            return np.hypot((lo - f.flon) * 111.195 * np.cos(np.radians((la + f.flat) / 2)), (la - f.flat) * 111.195)
        n_perm = int((dk(perm.lat.values, perm.lon.values) <= SEL_KM).sum())
        t = I[ok]
        dt = dk(t.lat.values, t.lon.values)
        n_tmp = int(((dt <= SEL_KM)).sum())
        if n_perm + n_tmp < SPARSE_N:
            need = SPARSE_N - n_perm - n_tmp
            cand = t.assign(d=dt)[(dt > SEL_KM) & (dt <= MAX_KM)].sort_values("d").head(need)
            for i in cand.index:
                if I.at[i, "select"] == "":
                    I.at[i, "select"] = f"nearest for sparse family {f.family_uid}"
            sparse.append(dict(family_uid=f.family_uid, region=f.region, n_perm_30km=n_perm, n_temp_30km=n_tmp,
                               added=";".join(I.loc[cand.index, "net"] + "." + I.loc[cand.index, "sta"])))
    pd.DataFrame(sparse).to_csv(f"{OUT}/sparse_families.csv", index=False)
    I.to_csv(f"{OUT}/inventory.csv", index=False)
    sel = I[I.select != ""].copy()
    # download rank under the 2.0 TB queue cap (lower first; within a rank nearest station first):
    # 1 XU (the discovery array itself, 2012-01..05)  2 TA (2012-16, all WA ETS 2012-2016)  3 X4 (BB 2016-20,
    # 4 ETS)  4 XD broadband (2014-16, ETS 2014/2016)  5 ZZ, YG (2012 summer, no ETS)  6 7D OBS (2012-14)
    # 7 X6, YR (10-day OBS)  8 XD nodes (1-C, 4-7 days in 2014-07/08, no ETS)
    rk = {"XU": 1, "TA": 2, "X4": 3, "ZZ": 5, "YG": 5, "7D": 6, "X6": 7, "YR": 7}
    sel["rank"] = [8 if (n == "XD" and g == "DP") else 4 if n == "XD" else rk.get(n, 9)
                   for n, g in zip(sel.net, sel.group)]
    sel = sel.sort_values(["rank", "dmin_km"]).reset_index(drop=True)
    sel["prio"] = np.arange(1, len(sel) + 1)
    sel["cum_gb"] = sel.est_gb.cumsum().round(1)
    sel.to_csv(f"{OUT}/selected_for_download.csv", index=False)
    # --- per network
    g = I.groupby("net")
    B = pd.DataFrame(dict(
        provider=g.provider.first(), stations=g.size(), start=g.start.min(), end=g.end.max(),
        bands=g.all_bands.apply(lambda v: ",".join(sorted(set(",".join(v).split(",")) - {""}))),
        three_c=g.three_c.sum(), seafloor=g.seafloor.sum(), fs=g.fs.apply(lambda v: ",".join(
            str(int(x)) for x in sorted(set(v)) if x)),
        dmin_km=g.dmin_km.min(), n_le10=g.dmin_km.apply(lambda v: int((v <= 10).sum())),
        n_le30=g.dmin_km.apply(lambda v: int((v <= 30).sum())),
        yrs_sampled=g.sample_years.sum(), yrs_data=g.sample_years_data.sum(),
        live_days=g.live_days.sum(), est_gb_all=g.est_gb.sum().round(1),
        selected=g.select.apply(lambda v: int((v != "").sum())),
        est_gb_selected=I[I.select != ""].groupby("net").est_gb.sum().round(1),
        status=g.status.apply(lambda v: "; ".join(f"{k} ({n})" for k, n in v.value_counts().items())),
        sensor=g.sensor.apply(lambda v: " | ".join(sorted(set(v) - {""}))[:90]))).reset_index()
    B["est_gb_selected"] = B.est_gb_selected.fillna(0)
    desc = S.groupby("net").net_desc.first()
    B["description"] = B.net.map(desc)
    B.to_csv(f"{OUT}/by_network.csv", index=False)
    print(B.drop(columns=["sensor", "description"]).to_string())
    print("selected", len(sel), "stations,", round(sel.est_gb.sum(), 1), "GB est;", len(sparse), "sparse families")


if __name__ == "__main__" and sys.argv[1] == "report":
    report()
