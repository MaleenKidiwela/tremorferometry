#!/usr/bin/env python
"""WA PRODUCTION W12: raw-count repeat windows for coda dv/v (-15 s .. +30 s around the expected S at EVERY recording
station), so the continuous raw can later be released (nothing is deleted here).

Uses lfe_south/scripts/cut_windows.py UNCHANGED (imported) with a run profile 'wa_ext' registered at run time (the SVI_REF
svi_ref_cut.py / duc_fin_cut.py pattern), local channel epochs (node tables; no FDSN query), and filtered root views
OUT/cut_roots/r<i>/ (symlinks to every station folder of each production raw root within 150 km of a Washington family,
except the excluded networks 7D, X6, YR, Z5 and the XD nodes; e.g. the California stations of lfe_south/traces are out).
  candidates  every production candidate (runs prod, rich; every network; data series) with P_max(own, nbr3, any; WAVE) >= the
              floor (0.2, or 0.1 when the refit risk >= 2 %: conf/refit_risk.json), + the top 10 per family-day (by P_max)
              for family-months with < 10 windows at that floor, + every candidate of the truth run (any P) + every array
              truth event of the truth months (open and sealed; the sealed look is already done when this runs).
              Deduplicated on (family, 20-Hz origin sample).
  offsets     S at a template station = origin + template offset + 0.5 s (template start + PRE; median over its nodes);
              any other station = straight-ray S time through the USGS CVM v1.7 1-D profile (royer_NW_WA / royer_SW_WA) from
              the family position (depth: catalogue, else 35 km) + the family's median (template - predicted) residual.
  pilot       the busiest month first; extrapolated store must be < 3.0 TB (pre-registered) or the cut stops for a decision.
  verify      20 random windows re-read vs the raw miniSEED (obspy slice): sample-identical.
  stackcheck  per region, stacks of cut windows (P_own_WAVE >= 0.75 core repeats, template stations, 1-8 Hz 20 Hz chain)
              vs the family's template and vs the same stack shifted 5-10 s (cut_windows rule: cc >= 0.5 and > control + 0.2).

  python scripts/wa_prod_cut.py pilot --workers 8 | cut --workers 8 | guard | verify --n 20 | stackcheck
"""
from __future__ import annotations
import os, sys, json, glob, argparse, warnings
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P
W = P.patch_wa()
sys.path.insert(0, f"{P.ROOT}/lfe_south/scripts")
import cut_windows as CW

RUN = "wa_ext"
MODE = {"kind": "full", "month": None}
LIMIT_TB = 3.0


DMAX_CUT_KM = 150.0


def station_coords():
    """netsta -> (lat, lon) from the roster, every bridge-region channel table and the temporary-network inventory."""
    co = {}
    for f in glob.glob(f"{W.EXT}/*_channels.csv"):
        C = pd.read_csv(f, usecols=["net", "sta", "lat", "lon"])
        for r in C.drop_duplicates(["net", "sta"]).itertuples():
            co.setdefault(f"{r.net}.{r.sta}", (float(r.lat), float(r.lon)))
    inv = pd.read_csv(f"{P.RES}/wa_temp_networks/inventory.csv", usecols=["net", "sta", "lat", "lon"])
    for r in inv.itertuples():
        co.setdefault(f"{r.net}.{r.sta}", (float(r.lat), float(r.lon)))
    for r in P.stations().itertuples():
        co[r.netsta] = (float(r.lat), float(r.lon))
    cx = f"{P.OUT}/cut_roots/coords_stationxml.csv"
    if not os.path.exists(cx):                      # local StationXML files (no query)
        from obspy import read_inventory
        rows = []
        for f in glob.glob(f"{P.ROOT}/lfe_array_discovery/results/stationxml/*.xml") + glob.glob(f"{P.ROOT}/lfe_picker/results/p4/stationxml/*.xml"):
            try:
                inv = read_inventory(f, level="station")
            except Exception:
                continue
            for n_ in inv:
                for st in n_:
                    rows.append((f"{n_.code}.{st.code}", st.latitude, st.longitude))
        P.ensure(os.path.dirname(cx))
        pd.DataFrame(rows, columns=["netsta", "lat", "lon"]).drop_duplicates("netsta").to_csv(cx, index=False)
    for r in pd.read_csv(cx).itertuples():
        co.setdefault(r.netsta, (float(r.lat), float(r.lon)))
    try:                                            # the SVI station table (CN.NLLB etc.)
        sys.path.insert(0, f"{P.ROOT}/lfe_extension/scripts")
        import c_common as CC
        for sta, r in CC.station_table().iterrows():
            co.setdefault(str(r.netsta), (float(r.lat), float(r.lon)))
    except Exception:
        pass
    return co


def roots_view():
    """symlink views of the production roots holding every recording station <= 150 km from a Washington family
    (excluded networks never; stations without known coordinates are left out and listed)."""
    out = []
    inv = pd.read_csv(f"{P.RES}/wa_temp_networks/inventory.csv")
    xdnodes = set((inv.net + "." + inv.sta.astype(str))[inv.all_bands.astype(str).str.startswith("DP")])
    co = station_coords()
    F = P.families().dropna(subset=["lat_use"])
    fla, flo = F.lat_use.to_numpy(float), F.lon_use.to_numpy(float)
    def near(ns):
        if ns not in co:
            return None
        return float(np.min(W.km(co[ns][0], co[ns][1], fla, flo))) <= DMAX_CUT_KM
    skipped = []
    for i, r in enumerate(P.prod_roots()):
        v = f"{P.OUT}/cut_roots/r{i}"
        P.ensure(v)
        if os.path.isdir(r):
            for ns in os.listdir(r):
                if "." not in ns or ns.split(".")[0] in P.BANNED_NETS or ns in xdnodes:
                    continue
                nr = near(ns)
                if not nr:
                    skipped.append((ns, r, "no coordinates" if nr is None else f"> {DMAX_CUT_KM:.0f} km from every family"))
                    continue
                dst = f"{v}/{ns}"
                if os.path.isdir(f"{r}/{ns}") and not os.path.lexists(dst):
                    os.symlink(f"{r}/{ns}", dst)
        out.append(v)
    json.dump({f"r{i}": r for i, r in enumerate(P.prod_roots())}, open(f"{P.OUT}/cut_roots/roots.json", "w"), indent=1)
    pd.DataFrame(skipped, columns=["netsta", "root", "why"]).drop_duplicates().to_csv(f"{P.OUT}/cut_roots/skipped_stations.csv", index=False)
    return out


def vs_harmonic(z, prof):
    d = 0.0; t = 0.0
    for r in prof.itertuples():
        top, bot = float(r.ztop_km), float(r.zbot_km)
        if top >= z:
            break
        h = min(bot, z) - top
        t += h / float(r.vs); d += h
    return d / t if t > 0 else 3.5


def offsets():
    F = P.families(); S = P.stations()
    V1 = pd.read_csv(f"{P.ROOT}/data/catalogue/families_v1_1.csv").set_index("family_uid")
    prof = {k: pd.read_csv(f"{P.ROOT}/data/velocity_models/usgs_cascadia_v1.7/profiles/royer_{k}.csv") for k in ("NW_WA", "SW_WA")}
    L = pd.read_csv(f"{P.OUT}/setup/tpl_library.csv", usecols=["fid", "netsta", "offset_s"])
    if os.path.exists(f"{P.OUT}/templates/tpl_qc_w4b.csv"):
        Q = pd.read_csv(f"{P.OUT}/templates/tpl_qc_w4b.csv")
        L = pd.concat([L, Q[Q.use.astype(bool)][["fid", "netsta", "offset_s"]]])
    tp = L.groupby(["fid", "netsta"]).offset_s.median()
    rows = []
    for r in F.itertuples():
        z = V1.depth_km.get(r.family_uid, np.nan)
        z = float(z) if np.isfinite(z) else 35.0
        pr = prof["SW_WA" if r.reg == "SW" else "NW_WA"]
        v = vs_harmonic(z, pr)
        x = W.km(S.lat.to_numpy(float), S.lon.to_numpy(float), r.lat_use, r.lon_use)
        pred = np.sqrt(x ** 2 + z ** 2) / v
        pm = dict(zip(S.netsta, pred))
        tps = {ns: float(o) + W.PRE for (f, ns), o in tp.items() if f == r.fid}
        res = [tps[ns] - pm[ns] for ns in tps if ns in pm]
        corr = float(np.median(res)) if res else 0.0
        for ns in S.netsta:
            if ns in tps:
                rows.append((str(r.fid), ns, tps[ns], "template"))
            else:
                rows.append((str(r.fid), ns, pm[ns] + corr, "moveout_cvm17"))
    return pd.DataFrame(rows, columns=["fam", "netsta", "offset_s", "offset_source"])


def floor():
    u = f"{P.OUT}/conf/floor_user.json"          # 2026-10-06: user decision overrides the automatic refit-risk floor
    if os.path.exists(u):
        return float(json.load(open(u))["floor"])
    p = f"{P.OUT}/conf/refit_risk.json"
    return float(json.load(open(p))["floor"]) if os.path.exists(p) else 0.2


def candidates():
    fl = floor()
    cols = ["run", "date", "net", "series", "fid", "t", "net_mad", "P_own_WAVE", "P_nbr3_WAVE", "P_any_WAVE", "P_max_WAVE"]
    sel, top = [], []
    months = [MODE["month"]] if MODE["kind"] == "pilot" else None
    for run in ("prod", "rich"):
        for p in sorted(glob.glob(f"{P.OUT}/conf/P/{run}/P_*.parquet")):
            if months and os.path.basename(p)[2:9] not in months:
                continue
            import pyarrow.parquet as _pq       # 2026-10-06: headline is the PNET model -> select on max(WAVE labels, PNET own)
            have_pn = "P_own_WAVE_PNET" in _pq.read_schema(p).names
            D = pd.read_parquet(p, columns=cols + (["P_own_WAVE_PNET"] if have_pn else []))
            D = D[D.series == "data"]
            if not have_pn:
                D["P_own_WAVE_PNET"] = np.nan
            D["P_sel"] = np.fmax(D.P_max_WAVE.to_numpy(float), D.P_own_WAVE_PNET.to_numpy(float))
            sel.append(D[D.P_sel >= fl])
            top.append(D.sort_values("P_sel", ascending=False).groupby(["fid", "date"]).head(10))
    S = pd.concat(sel, ignore_index=True); S["src_floor"] = True
    Tp = pd.concat(top, ignore_index=True)
    S["month"] = S.date.str[:7]; Tp["month"] = Tp.date.str[:7]
    n_fm = S.groupby(["fid", "month"]).size()
    few = set(n_fm[n_fm < 10].index) | (set(zip(Tp.fid, Tp.month)) - set(n_fm.index))
    Tp = Tp[[(f, m) in few for f, m in zip(Tp.fid, Tp.month)]].assign(src_top10=True)
    parts = [S, Tp]
    if MODE["kind"] != "pilot":
        tr = []
        for p in sorted(glob.glob(f"{P.OUT}/mf/truth/det_*.parquet")):
            D = pd.read_parquet(p, columns=["date", "net", "series", "fid", "t", "net_mad"])
            tr.append(D[D.series == "data"])
        if tr:
            parts.append(pd.concat(tr, ignore_index=True).assign(src_truth_run=True, run="truth"))
        T = pd.concat([pd.read_parquet(f"{P.OUT}/setup/truth.parquet"), pd.read_parquet(f"{P.OUT}/setup/sealed/truth_sealed.parquet")])
        parts.append(T[["fid", "t", "date"]].assign(src_truth_event=True, run="array_truth"))
    X = pd.concat(parts, ignore_index=True)
    X["key"] = np.round(X.t.astype(float) * W.FS).astype(np.int64)
    flags = ["src_floor", "src_top10", "src_truth_run", "src_truth_event"]
    for c in flags:
        X[c] = X[c].fillna(False).astype(bool) if c in X else False
    agg = {c: "max" for c in flags + ["net_mad", "P_own_WAVE", "P_nbr3_WAVE", "P_any_WAVE", "P_max_WAVE", "P_own_WAVE_PNET", "P_sel"] if c in X}
    agg["t"] = "first"; agg["run"] = lambda v: ",".join(sorted(set(map(str, v))))
    D = X.groupby(["fid", "key"]).agg(agg).reset_index().drop(columns=["key"])
    D["fam"] = D.fid.astype(int).astype(str)
    D["t"] = pd.to_datetime(D.t.astype(float), unit="s").astype("datetime64[ns]")
    D["runs"] = D.pop("run")
    P.log(f"[cut candidates] floor {fl}; {len(D)} candidates ({ {c: int(D[c].sum()) for c in flags} })")
    return D


def prof_wa():
    D = candidates()
    p = f"{P.OUT}/cut_offsets.csv"
    if not os.path.exists(p):
        offsets().to_csv(p, index=False)
    O = pd.read_csv(p, dtype={"fam": str})
    if MODE["kind"] == "pilot":
        m = MODE["month"]; t0 = f"{m}-01"; t1 = str((pd.Period(m, freq="M") + 1).start_time.date())
    else:
        t0, t1 = P.T0, str((pd.Timestamp(P.T1) + pd.Timedelta(days=1)).date())
    return dict(det=D, off=O, stations=sorted(P.stations().netsta), t0=t0, t1=t1, roots=roots_view(),
                flags=["src_floor", "src_top10", "src_truth_run", "src_truth_event"],
                stats=["net_mad", "P_own_WAVE", "P_nbr3_WAVE", "P_any_WAVE", "P_max_WAVE", "runs"], tpl=None,
                done="WA_PROD_DONE", failed="WA_PROD_FAILED", label="wa")


def local_channel_epochs(netstas, t0, t1, path):
    from wa_prod_setup import nodes_table
    N = nodes_table()
    rows = []
    for r in N[N.netsta.isin(netstas)].itertuples():
        for c in str(r.comps):
            rows.append(dict(netsta=r.netsta, loc="" if not isinstance(r.loc, str) else r.loc, cha=f"{r.band}{c}", start=r.start,
                             end=r.end, sr=r.fs, sensor=str(r.sensor), gain=np.nan, az=np.nan, dip=np.nan, prov="local:bridge_nodes"))
    E = pd.DataFrame(rows, columns=["netsta", "loc", "cha", "start", "end", "sr", "sensor", "gain", "az", "dip", "prov"])
    E.to_csv(path, index=False)
    E["loc"] = E["loc"].fillna("").astype(str)
    return E


def setup():
    CW.PROFILES[RUN] = prof_wa
    CW.channel_epochs = local_channel_epochs


def fix_offset_source(out):
    O = pd.read_csv(f"{P.OUT}/cut_offsets.csv", dtype={"fam": str})
    src = {(r.fam, r.netsta): r.offset_source for r in O.itertuples()}
    import pyarrow as pa, pyarrow.parquet as pq
    n = {}
    for p in sorted(glob.glob(f"{out}/index_*.parquet")):      # row group by row group: a month can be ~100 M rows
        pf = pq.ParquetFile(p)
        if not pf.metadata.num_rows or "offset_source" not in pf.schema_arrow.names:
            continue
        j = pf.schema_arrow.get_field_index("offset_source"); w = pq.ParquetWriter(p + ".tmp", pf.schema_arrow)
        for i in range(pf.num_row_groups):
            t = pf.read_row_group(i)
            s = [src.get((str(f), ns), "detection_time") if o != "detection_time" else "detection_time"
                 for f, ns, o in zip(t["fam"].to_pylist(), t["netsta"].to_pylist(), t["offset_source"].to_pylist())]
            w.write_table(t.set_column(j, pf.schema_arrow.field(j), pa.array(s, pf.schema_arrow.field(j).type)))
            for k, v in pd.Series(s).value_counts().items():
                n[k] = n.get(k, 0) + int(v)
        w.close(); os.replace(p + ".tmp", p)
    json.dump(n, open(f"{out}/offset_source_counts.json", "w"), indent=1)


def cutguard(out):
    SD = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{out}/stationdays_*.parquet"))], ignore_index=True)
    rows = [(r.day, r.netsta, "ok" if r.status == "ok" else ("nodata" if r.status == "no_band" else "error"), r.status)
            for r in SD.itertuples()]
    P.write_status("cut", "all", rows)
    P.guard("cut")


def to_chain(x, fs):
    """the frozen WA chain on a raw window: demean, 1-8 Hz Butterworth-4 zero-phase at the native rate, then 20 Hz."""
    from fractions import Fraction
    from scipy.signal import resample_poly
    y = W.bp(x.astype(np.float64) - x.mean(), fs, 1.0, 8.0)
    fr = Fraction(W.FS / fs).limit_denominator(100)
    return resample_poly(y, fr.numerator, fr.denominator) if fr != 1 else y


def stackcheck(out, nfam=8):
    from wa_prod_mf import load_library
    TPL = load_library()
    # filtered read (the whole production index does not fit in memory); same selection as before
    J = CW.load_index(out, columns=["det_id", "fam", "netsta", "cha", "h5_file", "dataset", "h5_start", "n_samples", "t0", "fs",
                                    "partial", "n_gap", "offset_source", "P_own_WAVE"],
                      filters=[("partial", "==", False), ("n_gap", "==", 0), ("offset_source", "==", "template"),
                               ("P_own_WAVE", ">=", P.P_CUT)])
    J["fam"] = J.fam.astype(str); J["month"] = J.h5_file.str[8:15]
    F = P.families(); reg = dict(zip(F.fid.astype(str), F.reg))
    g = J.groupby(["fam", "netsta", "month"]).det_id.nunique().reset_index(name="n").sort_values("n", ascending=False)
    g = g[g.n >= 50]
    g["reg"] = g.fam.map(reg)
    g = g.drop_duplicates("fam").groupby("reg").head(max(2, nfam // 4)).head(nfam)
    res = []
    i0 = int(round((CW.PRE_S - W.PRE) * W.FS)); n = W.NT
    slot = lambda c: 0 if c[-1] == "Z" else 1 if c[-1] in "N1" else 2
    for b in g.itertuples():
        keys = [k for k in TPL if k[0] == int(b.fam) and k[1] == b.netsta]
        if not keys:
            continue
        T = TPL[keys[0]]["tpl"]
        K = J[(J.fam == b.fam) & (J.netsta == b.netsta) & (J.month == b.month)]
        stk = np.zeros((3, n + 400)); cnt = 0
        for det, gg in K.groupby("det_id"):
            A = np.zeros((3, i0 + n + 220)); ok = True
            for r in gg.to_dict("records"):
                x, t0, fs = CW.read_window(r, out)
                y = to_chain(x, fs)
                if y.size < i0 + n + 220:
                    ok = False; break
                A[slot(r["cha"])] = y[:i0 + n + 220]
            if not ok:
                continue
            r_ = np.sqrt(np.mean(A[:, i0:i0 + n] ** 2)) + 1e-30
            seg = A[:, i0 - 200:i0 - 200 + n + 400]
            if seg.shape[1] != n + 400:
                continue
            stk += seg / r_; cnt += 1
        if cnt < 50:
            continue

        def ccs(lag):
            w = stk[:, 200 + lag:200 + lag + n]; o = []
            for c in range(3):
                u = w[c] - w[c].mean(); v = T[c] - T[c].mean()
                if np.any(u) and np.any(v):
                    o.append(float(u @ v / (np.sqrt(u @ u * (v @ v)) + 1e-30)))
            return float(np.mean(o)) if o else np.nan
        cc = max(ccs(l) for l in range(-10, 11))
        ctl = max(abs(ccs(l)) for l in list(range(-200, -99)) + list(range(100, 181)))
        res.append(dict(fam=b.fam, region=b.reg, netsta=b.netsta, month=b.month, n_windows=cnt, cc_stack_vs_template=round(cc, 3),
                        cc_control_shift_5_10s=round(ctl, 3), pass_=bool(cc >= 0.5 and cc > ctl + 0.2)))
    R = pd.DataFrame(res); R.to_csv(f"{out}/stackcheck_families.csv", index=False)
    ok = bool(len(R) >= 1 and R.pass_.mean() >= 0.5)
    json.dump(dict(run=RUN, n_families=len(R), n_pass=int(R.pass_.sum()) if len(R) else 0, families=res, **{"pass": ok},
                   rule="stack of cut raw windows (P_own_WAVE >= 0.75, template stations) through the WA 1-8 Hz 20 Hz chain vs the "
                        "family template (+-0.5 s); control = same stack shifted 5-10 s; pass = cc >= 0.5 and cc > control + 0.2 for "
                        ">= half of the families checked (per region, up to 8)"), open(f"{out}/STACKCHECK.json", "w"), indent=1)
    P.log(f"[stackcheck] {int(R.pass_.sum()) if len(R) else 0}/{len(R)} pass; overall {ok}")
    return ok


def busiest():
    best = (None, -1)
    for p in sorted(glob.glob(f"{P.OUT}/conf/P/prod/P_*.parquet")):
        D = pd.read_parquet(p, columns=["series", "P_max_WAVE"])
        n = int(((D.series == "data") & (D.P_max_WAVE >= floor())).sum())
        if n > best[1]:
            best = (os.path.basename(p)[2:9], n)
    return best[0]


def pilot(a):
    m = a.month or busiest()
    out = f"{P.OUT}/cut_pilot"; P.ensure(out)
    MODE.update(kind="pilot", month=m)
    CW.cut(RUN, a.workers, [m], out, a.max_days); fix_offset_source(out)
    ok_v = CW.verify(RUN, a.n, out)
    Sm = pd.read_csv(f"{out}/cut_summary.csv"); r = Sm[Sm.month == m].iloc[0]
    ratio = float(r.cut_GB) / max(float(r.raw_GB), 1e-9)
    tot = 0
    for v in roots_view():
        for ns in os.listdir(v):
            d = os.path.realpath(f"{v}/{ns}")
            for y in os.listdir(d):
                if y.isdigit() and P.T0[:4] <= y <= P.T1[:4]:
                    for f in os.scandir(f"{d}/{y}"):
                        if f.name.endswith(".mseed"):
                            try:
                                tot += f.stat().st_size
                            except OSError:
                                pass
    ex = dict(pilot_month=m, cut_GB=float(r.cut_GB), raw_GB=float(r.raw_GB), cut_over_raw=ratio, raw_total_TB=tot / 1e12,
              extrapolated_TB=ratio * tot / 1e12, n_windows=int(r.n_windows), minutes=float(r.minutes), limit_TB=LIMIT_TB,
              note="the pilot month is the busiest production month; full-store estimate = pilot cut/raw ratio x all raw in the roots")
    ex["pass"] = bool(ex["extrapolated_TB"] < LIMIT_TB)
    json.dump(dict(verify=ok_v, extrapolation=ex), open(f"{out}/PILOT.json", "w"), indent=1)
    P.log(f"[pilot] {ex}; verify {ok_v}")
    sys.exit(0 if (ok_v and ex["pass"]) else 6)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["pilot", "cut", "guard", "verify", "stackcheck"])
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--month", default="")
    ap.add_argument("--months", default=""); ap.add_argument("--max-days", type=int, default=0); ap.add_argument("--n", type=int, default=20)
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    setup()
    if a.what == "pilot":
        pilot(a)
    out = P.CUT; P.ensure(out)
    if a.what == "cut":
        CW.cut(RUN, a.workers, [m for m in a.months.split(",") if m] or None, out, a.max_days); fix_offset_source(out)
    elif a.what == "guard":
        cutguard(out)
    elif a.what == "verify":
        sys.exit(0 if CW.verify(RUN, a.n, out) else 5)
    elif a.what == "stackcheck":
        sys.exit(0 if stackcheck(out) else 6)
