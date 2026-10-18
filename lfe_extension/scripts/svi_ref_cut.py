#!/usr/bin/env python
"""SVI_REF R8 -- raw-count repeat windows for dv/v (-15 s .. +30 s around the expected S at every station) for the
SVI catalogue, so the continuous raw miniSEED can later be released (nothing is deleted here).

Uses lfe_south/scripts/cut_windows.py UNCHANGED (imported) with a run profile 'svi_ext' registered at run time
(the duc_fin_cut.py pattern):
  candidates  pilot: every >= 8 x MAD v1.1 candidate of the pilot month (results/trackC/detections_<M>.parquet)
              full : the UNION of v1.1 (all months) and the SVI_REF production run (OUT/feat/prod), deduplicated on
                     (family, origin sample); + every truth-year candidate (OUT/feat/t*); + Bostock's catalogued
                     events 2010-2012 (2013 is cut by `cut --bostock2013` only AFTER the one-time sealed check)
  offsets     S arrival at each station: template stations = frozen placement + 0.5 s (template start + PRE);
              every other station = Crosson 1-D S time (c_vel1d) from the family location + the family's median
              (template - predicted) correction -> offset_source per row (template / moveout_crosson1d)
  roots       lfe_extension/traces + a symlink root (OUT/cut_roots/rep) holding ONLY the 19 SVI stations of the
              read-only lfe_catalog_repair/traces (2010 / 2012 truth windows, 2017 'present_elsewhere' days)
  metadata    channel epochs from local tables (no FDSN query, nothing downloaded)
Outputs CUT (results/cut/svi_ext/): windows_<M>.h5, segments_<M>.parquet, index_<M>.parquet, detections.parquet,
stationdays*.parquet, cut_summary.csv, VERIFY.txt, STACKCHECK.json, pilot_extrapolation.json.

  python scripts/svi_ref_cut.py pilot --month 2016-09 --workers 6
  python scripts/svi_ref_cut.py cut --workers 8         python scripts/svi_ref_cut.py guard
  python scripts/svi_ref_cut.py verify --n 20            python scripts/svi_ref_cut.py stackcheck
"""
from __future__ import annotations
import os, sys, json, glob, argparse, warnings
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import svi_ref_common as R
import c_common as C
sys.path.insert(0, f"{R.SOUTH}/scripts")
import cut_windows as CW

RUN = "svi_ext"
MODE = {"kind": "full", "month": None, "out": R.CUT}
ROOT_REP = f"{R.OUT}/cut_roots/rep"


def rep_root():
    """symlink root with only the SVI stations of lfe_catalog_repair/traces (never written through)."""
    R.ensure(ROOT_REP)
    stas = sorted(set(C.epochs().netsta))
    for ns in stas:
        src = f"{R.REP}/traces/{ns}"; dst = f"{ROOT_REP}/{ns}"
        if os.path.isdir(src) and not os.path.lexists(dst):
            os.symlink(src, dst)
    return ROOT_REP


def offsets():
    import c_vel1d as v1d
    fz = C.frozen(); L = float(fz["template_length_s"])
    tpl_path, _ = C.resolve_templates(fz, "")
    tpl, sp, F, wgt, mv, ampref = C.load_frozen_tables(fz, tpl_path, L)
    off = C.placement(tpl, sp, F, mv, C.station_table(), L)
    ep = C.epochs(); co = ep.groupby("netsta")[["lat", "lon"]].first()
    nsof = R.netsta_of()
    Fi = F.set_index("bkey")
    rows = []
    for f in sorted(int(x) for x in Fi.index):
        m = Fi.loc[f]
        pred = {ns: float(v1d.tt(float(C.epi_km(float(m.lat), float(m.lon), float(co.loc[ns, "lat"]), float(co.loc[ns, "lon"]))),
                                 float(m.dep), "S")) for ns in co.index}
        tp = {nsof[s]: off[(float(f), s)] + C.PRE for (ff, s) in off if ff == float(f) and s in nsof}
        res = [tp[ns] - pred[ns] for ns in tp]
        corr = float(np.median(res)) if res else 0.0
        for ns in co.index:
            if ns in tp: rows.append((str(f), ns, tp[ns], "template"))
            else: rows.append((str(f), ns, pred[ns] + corr, "moveout_crosson1d"))
    return pd.DataFrame(rows, columns=["fam", "netsta", "offset_s", "offset_source"])


def v11(months=None):
    ps = sorted(glob.glob(f"{C.TRC}/detections_*.parquet"))
    if months: ps = [p for p in ps if os.path.basename(p)[11:18] in months]
    cols = ["fam", "t", "net_mad", "net_cc", "n_sites", "n_sta", "cat_core", "cat_ext", "tier_core", "tier_ext", "hf_pass", "tremor_flag"]
    D = pd.concat([pd.read_parquet(p, columns=cols) for p in ps], ignore_index=True)
    D["src_v11"] = True
    return D


def feat_rows(run_glob):
    ps = sorted(glob.glob(f"{R.OUT}/feat/{run_glob}/*.parquet"))
    ps = [p for p in ps if not p.endswith("_fd.parquet") and os.path.getsize(p) > 0]
    out = []
    for p in ps:
        try:
            A = pd.read_parquet(p, columns=["series", "fam", "t", "net_mad", "net_cc", "n_sites", "n_sta", "run"])
        except Exception:
            continue
        out.append(A[A.series == "data"].drop(columns=["series"]))
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=["fam", "t", "net_mad", "net_cc", "n_sites", "n_sta", "run"])


def candidates():
    if MODE["kind"] == "pilot":
        D = v11([MODE["month"]])
    else:
        A = v11()
        B = feat_rows("prod"); B["src_v2"] = True
        Tt = pd.concat([feat_rows("t*"), feat_rows("r3t*")], ignore_index=True); Tt["src_truth"] = True
        X = pd.concat([A, B, Tt], ignore_index=True)
        X["key"] = np.round(X.t.astype(float) * R.FS).astype(np.int64)
        agg = {c: "max" for c in ["net_mad", "net_cc", "n_sites", "n_sta"]}
        for c in ["src_v11", "src_v2", "src_truth", "cat_core", "cat_ext", "tier_core", "tier_ext", "hf_pass", "tremor_flag"]:
            X[c] = X[c].fillna(False).astype(bool) if c in X else False; agg[c] = "max"
        agg["t"] = "first"
        D = X.groupby(["fam", "key"]).agg(agg).reset_index().drop(columns=["key"])
        Bo = R.bostock([2010, 2011, 2012])
        if MODE.get("bostock2013"):
            Bo = R.bostock([2013], allow_sealed=True)
            D = D.iloc[0:0]
        Bo = Bo.rename(columns={"t_origin": "t"})[["fam", "t"]].assign(src_bostock=True)
        D = pd.concat([D, Bo], ignore_index=True)
        for c in D.columns:
            if c.startswith("src_") or c in ("cat_core", "cat_ext", "tier_core", "tier_ext", "hf_pass", "tremor_flag"):
                D[c] = D[c].fillna(False).astype(bool)
    D["fam"] = D.fam.astype(int).astype(str)
    D["t"] = pd.to_datetime(D.t.astype(float), unit="s").astype("datetime64[ns]")
    return D


def prof_svi():
    D = candidates()
    flags = [c for c in ["src_v11", "src_v2", "src_truth", "src_bostock", "cat_core", "cat_ext", "tier_core", "tier_ext", "hf_pass",
                         "tremor_flag"] if c in D.columns]
    p = f"{R.OUT}/cut_offsets.csv"
    if not os.path.exists(p):
        offsets().to_csv(p, index=False)
    O = pd.read_csv(p, dtype={"fam": str})
    if MODE["kind"] == "pilot":
        m = MODE["month"]; t0 = f"{m}-01"; t1 = str((pd.Period(m, freq="M") + 1).start_time.date())
    elif MODE.get("bostock2013"):
        t0, t1 = "2013-08-30", "2013-10-29"
    else:
        t0, t1 = "2010-08-01", R.PROD_T1
    stations = sorted(C.epochs().netsta.unique())
    return dict(det=D, off=O, stations=stations, t0=t0, t1=t1, roots=[C.TRACES, rep_root()], flags=flags,
                stats=["net_mad", "net_cc", "n_sites", "n_sta"], tpl=None, done="SVI_REF_DONE", failed="SVI_REF_FAILED", label="svi")


def local_channel_epochs(netstas, t0, t1, path):
    import svi_ref_r2 as R2
    G = R2.gain_table(); nsof = R.netsta_of()
    rows = []
    for r in G.itertuples():
        ns = nsof.get(r.sta)
        if ns is None: continue
        rows.append(dict(netsta=ns, loc="", cha=r.cha, start=r.start, end=r.end, sr=np.nan, sensor="", gain=r.gain, az=np.nan,
                         dip=np.nan, prov=f"local:{r.src}"))
    E = pd.DataFrame(rows)
    E.to_csv(path, index=False); E["loc"] = E["loc"].fillna("").astype(str)
    return E


def setup():
    CW.PROFILES[RUN] = prof_svi
    CW.channel_epochs = local_channel_epochs


def fix_offset_source(out):
    O = pd.read_csv(f"{R.OUT}/cut_offsets.csv", dtype={"fam": str})
    src = {(r.fam, r.netsta): r.offset_source for r in O.itertuples()}
    n = {}
    for p in sorted(glob.glob(f"{out}/index_*.parquet")):
        I = pd.read_parquet(p)
        if not len(I) or "offset_source" not in I: continue
        s = np.array([src.get((str(f), ns), "detection_time") if osrc != "detection_time" else "detection_time"
                      for f, ns, osrc in zip(I.fam.astype(str), I.netsta, I.offset_source)], dtype=object)
        I["offset_source"] = s; I.to_parquet(p, index=False)
        for k, v in pd.Series(s).value_counts().items(): n[k] = n.get(k, 0) + int(v)
    json.dump(n, open(f"{out}/offset_source_counts.json", "w"), indent=1)
    R.log(f"[cut] offset sources: {n}")


def cutguard(out):
    """every station-day the manifest (prod) or the raw roots (truth years) say has a raw file must have been read ok."""
    SD = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{out}/stationdays_*.parquet"))], ignore_index=True)
    months = {os.path.basename(p)[5:12] for p in glob.glob(f"{out}/done_*.json")}
    seen = {(r.netsta, r.day): r.status for r in SD.itertuples()}
    man = C.manifest()
    man = man[man.date.str[:7].isin(months) & man.status.isin(["data", "present_elsewhere"])]
    last = SD.day.max() if len(SD) else ""
    rows = []
    for r in man.itertuples():
        if r.date > last: continue
        st = seen.get((r.netsta, r.date))
        if st is None: rows.append((r.date, r.netsta, "missing", "manifest has a raw file, the cut did not see it"))
        elif st == "ok": rows.append((r.date, r.netsta, "ok", ""))
        else: rows.append((r.date, r.netsta, "error", st))
    have = {(d, n) for d, n, _, _ in rows}
    for r in SD.itertuples():
        if (r.day, r.netsta) in have: continue
        rows.append((r.day, r.netsta, "ok" if r.status == "ok" else "error", "" if r.status == "ok" else r.status))
    R.write_status("cut", "all", rows)
    R.guard("cut")


def to_chain(x, fs):
    """the frozen SVI chain on a raw window: 1-8 Hz Butterworth-4 zero-phase at the native rate, then 20 Hz."""
    from scipy.signal import butter, sosfiltfilt, resample_poly
    y = sosfiltfilt(butter(4, [1.0 / (fs / 2), 8.0 / (fs / 2)], btype="band", output="sos"), x.astype(float) - x.mean())
    q = int(round(fs / 20.0))
    return resample_poly(y, 1, q) if q > 1 else y


def stackcheck(out, nfam=8):
    """per-family check: stack of the cut windows (cat_core v1.1 repeats, template stations, full windows) through the
    frozen chain vs the frozen 2012 template (3-C, +-0.5 s); control = the same stack shifted 5-10 s."""
    import c_chanorder as CO
    fz = C.frozen(); L = float(fz["template_length_s"])
    tpl_path, _ = C.resolve_templates(fz, "")
    tpl = np.load(tpl_path, allow_pickle=True)["tpl"].item()
    nsof = R.netsta_of(); sof = {v: k for k, v in nsof.items()}
    I = CW.load_index(out); I["fam"] = I.fam.astype(str); I["month"] = I.h5_file.str[8:15]
    flag = "cat_core" if "cat_core" in I else "src_v11"
    J = I[(I[flag]) & (~I.partial) & (I.n_gap == 0) & (I.offset_source == "template")]
    g = J.groupby(["fam", "netsta", "month"]).det_id.nunique().reset_index(name="n").sort_values("n", ascending=False)
    g = g[g.n >= 50].drop_duplicates("fam").head(nfam)
    res = []
    n = int(round(L * 20)); i0 = int(round((CW.PRE_S - C.PRE) * 20))
    for b in g.itertuples():
        sta = sof.get(b.netsta); key = (float(b.fam), sta)
        if key not in tpl or L not in tpl[key]: continue
        T = tpl[key][L]
        K = J[(J.fam == b.fam) & (J.netsta == b.netsta) & (J.month == b.month)]
        stk = np.zeros((3, n + 400)); cnt = 0
        for det, gg in K.groupby("det_id"):
            chans = {r["cha"]: r for r in gg.to_dict("records")}
            try:
                order = CO.resolve_order(list(chans), station=b.netsta, strict=True)
            except Exception:
                continue
            per = []
            for c in order:
                x, t0, fs = CW.read_window(chans[c], out)
                y = to_chain(x, fs)
                if y.size < i0 + n + 220: break
                per.append(y[:i0 + n + 220])
            if len(per) != 3: continue
            A = np.vstack(per); r_ = np.sqrt(np.mean(A[:, i0:i0 + n] ** 2)) + 1e-30
            seg = A[:, max(i0 - 200, 0):max(i0 - 200, 0) + n + 400]
            if seg.shape[1] != n + 400: continue
            stk += seg / r_; cnt += 1
        if cnt < 50: continue                      # a stack of fewer windows cannot reach the 0.5 bar reliably
        def ccs(lag):
            w = stk[:, 200 + lag:200 + lag + n]
            out_ = []
            for c in range(3):
                u = w[c] - w[c].mean(); v = T[c] - T[c].mean()
                if not np.any(u) or not np.any(v): continue
                out_.append(float(u @ v / (np.sqrt(u @ u * (v @ v)) + 1e-30)))
            return float(np.mean(out_)) if out_ else np.nan
        cc = max(ccs(l) for l in range(-10, 11))
        ctl = max(abs(ccs(l)) for l in list(range(-200, -99)) + list(range(100, 181)))
        res.append(dict(fam=b.fam, netsta=b.netsta, month=b.month, n_windows=cnt, cc_stack_vs_template=round(cc, 3),
                        cc_control_shift_5_10s=round(ctl, 3), pass_=bool(cc >= 0.5 and cc > ctl + 0.2)))
    Rr = pd.DataFrame(res); Rr.to_csv(f"{out}/stackcheck_families.csv", index=False)
    ok = bool(len(Rr) >= 1 and Rr.pass_.mean() >= 0.5 and Rr.pass_.iloc[0])
    summ = dict(run=RUN, n_families=len(Rr), n_pass=int(Rr.pass_.sum()) if len(Rr) else 0, families=res, **{"pass": ok},
                note="stack of cut raw windows through the frozen 1-8 Hz / 20 Hz chain vs the frozen 2012 template "
                     "(+-0.5 s); control = the same stack shifted 5-10 s; pass = cc >= 0.5 and cc > control + 0.2 "
                     "(cut_windows.py rule) for the best-populated family and for >= half of the families with >= 50 windows")
    json.dump(summ, open(f"{out}/STACKCHECK.json", "w"), indent=1)
    R.log(f"[stackcheck] {summ['n_pass']}/{summ['n_families']} families pass; overall {ok}")
    return ok


def pilot_extrapolate(out, month):
    Sm = pd.read_csv(f"{out}/cut_summary.csv")
    r = Sm[Sm.month == month].iloc[0]
    tot = 0
    for root in (C.TRACES, ROOT_REP):
        for ns in os.listdir(root):
            for y in os.listdir(f"{root}/{ns}"):
                if not y.isdigit(): continue
                dd = f"{root}/{ns}/{y}"
                if root == ROOT_REP and y not in ("2010", "2012", "2017"): continue
                for fn in os.listdir(dd):
                    if fn.endswith(".mseed"): tot += os.path.getsize(f"{dd}/{fn}")
    ratio = float(r.cut_GB) / max(float(r.raw_GB), 1e-9)
    ex = dict(pilot_month=month, pilot_cut_GB=float(r.cut_GB), pilot_raw_GB=float(r.raw_GB), cut_over_raw=ratio,
              n_windows=int(r.n_windows), n_det=int(r.n_det), minutes=float(r.minutes), raw_total_TB=tot / 1e12,
              extrapolated_TB_upper=ratio * tot / 1e12,
              note="upper bound: the pilot month is the busiest v1.1 month, and the full cut also adds the v2 / truth "
                   "candidates, which mostly coincide with v1.1 ones")
    ex["pass_lt_0p5TB"] = bool(ex["extrapolated_TB_upper"] < 0.5)
    json.dump(ex, open(f"{out}/pilot_extrapolation.json", "w"), indent=1)
    R.log(f"[pilot] {ex}")
    return ex


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["pilot", "cut", "guard", "verify", "stackcheck", "fixsrc"])
    ap.add_argument("--workers", type=int, default=6); ap.add_argument("--month", default="")
    ap.add_argument("--months", default=""); ap.add_argument("--max-days", type=int, default=0); ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--bostock2013", action="store_true")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    setup()
    if a.what == "pilot":
        out = f"{R.OUT}/cut_pilot"; R.ensure(out); MODE.update(kind="pilot", month=a.month)
        CW.cut(RUN, a.workers, [a.month], out, a.max_days); fix_offset_source(out)
        ok_v = CW.verify(RUN, a.n, out); ok_s = stackcheck(out); ex = pilot_extrapolate(out, a.month)
        json.dump(dict(verify=ok_v, stackcheck=ok_s, extrapolation=ex), open(f"{out}/PILOT.json", "w"), indent=1)
        sys.exit(0 if ok_v else 5)
    out = R.CUT; R.ensure(out)
    if a.bostock2013:
        out = f"{R.CUT}_bostock2013"; R.ensure(out); MODE.update(bostock2013=True)
    if a.what == "cut":
        CW.cut(RUN, a.workers, [m for m in a.months.split(",") if m] or None, out, a.max_days); fix_offset_source(out)
    elif a.what == "fixsrc": fix_offset_source(out)
    elif a.what == "guard": cutguard(out)
    elif a.what == "verify": sys.exit(0 if CW.verify(RUN, a.n, out) else 5)
    elif a.what == "stackcheck": sys.exit(0 if stackcheck(out) else 6)
