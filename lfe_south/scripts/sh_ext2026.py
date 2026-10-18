#!/usr/bin/env python
"""SH_EXT2026 (2026-09-30): extend Shelly et al. (2025)'s 27 Mendocino LFE families past his catalogue end
(2024-11-30) to the latest available day, with the SH_V3B reproduction detector (it0 templates, Text S1
catalogue rules).  Nothing of v3b is edited; outputs go to results/sh_ext2026/.

  epochs   FDSN channel epochs of the 17 S1 stations (NCEDC for BK/NC, EarthScope for PB) ->
           station_epochs.csv + epoch_decisions.csv.  An extension-period epoch whose sensor / gain / azimuth /
           coordinates / sample rate / channel set differ from the epoch in use at 2024-11-30 => 'retemplate'
           (re-stack that station from our own >= 9xMAD S1 detections made WITHOUT it, Text S1 recipe);
           a pure metadata split (all attributes identical) => 'keep'.
  days     day-list JSONs (results/sh_ext2026/days/*.json, {"Shelly2025|<fam>": [dates]}).
  plan     download plan: 17 stations x 2023-12-01..LAST inside their FDSN epochs.
  prep     v2 cleaning (sh_v2_common.prep_day, unchanged) on the 17 S1 stations -> results/sh_ext2026/prep/.
  s1cat    Text S1 rules (9xMAD, mean cc >= 0.12, summed cc >= 2.75, 4 s cross-family declustering).
  tplset   compose the template sets (it0 without KMR; it0 with KMR re-stacked per new epoch).
  cont     continuity check on 2024 held-out test days vs his catalogue and vs v3b on the same days (gate).
  catalog  merge the MF segments -> S1 catalogue, raw 8x candidates, daily counts, monthly series, join table.
  tplcheck template stability: stacks of our own detections per window vs the template in use.
  report   figures + results/SH_EXT2026_RESULTS.md.
"""
import os, sys, json, argparse, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t_common import RES, WS, CACHE, FS, log, npz_path

EXT = os.environ.get("SH_EXT_DIR", f"{RES}/sh_ext2026"); DAYS = f"{EXT}/days"; PREPD = f"{EXT}/prep"; V3B = f"{RES}/sh_v3b"
STAS = ["BK.BLAS", "BK.BRIC", "BK.GHOP", "BK.HALS", "BK.JEWT", "BK.PRDS", "BK.SPRL", "BK.RVIT", "BK.HUNT",
        "NC.KMPB", "NC.KMR", "NC.KPR", "NC.KBU", "PB.B049", "PB.B932", "PB.B934", "PB.B933"]
SRC = "Shelly2025"; FAMS = [str(i) for i in range(1, 28)]
T_JOIN0 = pd.Timestamp("2023-12-01")      # start of our continuous run (12 months of overlap with his catalogue)
T_HIS1 = pd.Timestamp("2024-12-01")       # his catalogue ends 2024-11-30
LAST = pd.Timestamp(os.environ.get("SH_EXT_LAST", "2026-09-28"))
T_END = LAST + pd.Timedelta(days=1)
KMR_SW1 = pd.Timestamp("2025-09-17")      # NC.KMR STS-2 -> Trillium Horizon 120 (23:18Z)
KMR_SW2 = pd.Timestamp("2026-03-12")      # NC.KMR Horizon 120 -> Trillium Compact Posthole (20:43Z)
MIN_N_RETPL = 20                          # a re-stacked KMR template needs >= 20 events, else KMR is dropped for that family-epoch
TPL0 = f"{V3B}/tpl_it0"                   # the validated v3b it0 templates


def days_json(path, days):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    ds = sorted({str(pd.Timestamp(d).date()) for d in days})
    json.dump({f"{SRC}|{f}": ds for f in FAMS}, open(path, "w"))
    return ds


def rng(a, b):                           # inclusive a, exclusive b
    return list(pd.date_range(a, b, freq="D", inclusive="left"))


# ------------------------------------------------------------------------------------------------ epochs
def epochs():
    raw = f"{EXT}/station_epochs_raw.csv"
    if not os.path.exists(raw):
        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        rows = []
        for ns in STAS:
            n, s = ns.split(".")
            for prov in (["NCEDC", "EARTHSCOPE"] if n in ("BK", "NC") else ["EARTHSCOPE", "NCEDC"]):
                try:
                    inv = Client(prov, timeout=120).get_stations(network=n, station=s, channel="HH?,BH?,EH?,SH?",
                                                                 level="response", starttime=UTCDateTime("2017-06-01"))
                except Exception as e:
                    log(f"  {ns} {prov} ERR {str(e)[:80]}"); continue
                for net in inv:
                    for st in net:
                        for ch in st:
                            sens = ch.response.instrument_sensitivity
                            rows.append(dict(prov=prov, netsta=ns, loc=ch.location_code, cha=ch.code, sr=ch.sample_rate,
                                             start=str(ch.start_date)[:19], end=str(ch.end_date)[:19] if ch.end_date else "",
                                             sensor=(ch.sensor.description or ch.sensor.model or ch.sensor.type or "") if ch.sensor else "",
                                             gain=sens.value if sens else None, gfreq=sens.frequency if sens else None,
                                             az=ch.azimuth, dip=ch.dip, lat=ch.latitude, lon=ch.longitude,
                                             nstage=len(ch.response.response_stages)))
                break
        pd.DataFrame(rows).to_csv(raw, index=False)
    D = pd.read_csv(raw, dtype={"loc": str}); D["loc"] = D["loc"].fillna("")
    D["band"] = D.cha.str[:2]
    en = D.end.fillna("").astype(str).str[:19]
    D["t0"] = pd.to_datetime(D.start); D["t1"] = pd.to_datetime(en.where(~(en.eq("") | en.str.startswith("3000") | en.eq("nan")), "2100-01-01T00:00:00"))
    out = []; dec = []
    for ns in STAS:
        g = D[D.netsta == ns]
        band = next(b for b in ["HH", "BH", "EH", "SH"] if (g.band == b).any())     # t3_cache band preference
        g = g[g.band == band]
        E = []
        for (t0, lc), h in g.groupby(["t0", "loc"]):
            h = h.sort_values("cha")
            E.append(dict(netsta=ns, loc=lc, band=band, start=t0, end=h.t1.iloc[0], cha=",".join(h.cha), sr=float(h.sr.iloc[0]),
                          sensor=str(h.sensor.iloc[0]), gain=",".join(f"{x:.4g}" for x in h.gain), az=",".join(f"{x:g}" for x in h.az),
                          lat=round(float(h.lat.iloc[0]), 4), lon=round(float(h.lon.iloc[0]), 4)))
        E = pd.DataFrame(E).sort_values("start").reset_index(drop=True)
        ref = E[(E.start <= T_HIS1 - pd.Timedelta(seconds=1)) & (E.end > T_HIS1 - pd.Timedelta(seconds=1))]
        ref = ref.iloc[-1] if len(ref) else None
        for _, e in E.iterrows():
            row = dict(e)
            in_ext = (e.end > T_HIS1) and (e.start < T_END)
            in_build = (e.start < T_HIS1) and (e.end > pd.Timestamp("2018-01-01"))
            row.update(in_build_period=bool(in_build), in_extension=bool(in_ext))
            if ref is not None:
                gs = np.array([float(x) for x in e.gain.split(",")]); gr = np.array([float(x) for x in ref.gain.split(",")])
                row.update(sensor_changed=e.sensor != ref.sensor, cha_changed=e.cha != ref.cha, sr_changed=e.sr != ref.sr,
                           az_changed=e.az != ref.az, coord_changed=(e.lat != ref.lat) or (e.lon != ref.lon),
                           gain_changed=bool(len(gs) != len(gr) or np.any(np.abs(gs / gr - 1) > 0.01)), loc_changed=e["loc"] != ref["loc"])
            out.append(row)
        O = pd.DataFrame([r for r in out if r["netsta"] == ns])
        X = O[O.in_extension & (O.start >= T_HIS1)]
        flags = ["sensor_changed", "cha_changed", "sr_changed", "az_changed", "coord_changed", "gain_changed", "loc_changed"]
        for _, e in X.iterrows():
            ch = [f for f in flags if bool(e.get(f, False))]
            dec.append(dict(netsta=ns, epoch_start=e.start, epoch_end=e.end, sensor=e.sensor, changes=";".join(ch) or "none",
                            decision="retemplate" if ch else "keep (metadata-only epoch split: identical sensor, gain, azimuth, coordinates, rate)"))
        # epochs that began inside the template-build period but late (<= 12 months before its end): informational
        L = O[O.in_build_period & (O.start >= pd.Timestamp("2023-06-01")) & (O.start < T_HIS1)]
        for _, e in L.iterrows():
            prev = O[O.start < e.start]
            if not len(prev): continue
            p = prev.iloc[-1]
            ch = [k for k, v in dict(sensor=e.sensor != p.sensor, gain=e.gain != p.gain, az=e.az != p.az,
                                     coords=(e.lat != p.lat) or (e.lon != p.lon), cha=e.cha != p.cha).items() if v]
            dec.append(dict(netsta=ns, epoch_start=e.start, epoch_end=e.end, sensor=e.sensor,
                            changes=("within build period vs previous epoch: " + (";".join(ch) or "none")),
                            decision="keep (epoch already inside the template-build days; checked by the 2024 continuity days + tplcheck)"))
    pd.DataFrame(out).to_csv(f"{EXT}/station_epochs.csv", index=False)
    Dd = pd.DataFrame(dec); Dd.to_csv(f"{EXT}/epoch_decisions.csv", index=False)
    log("[epochs]\n" + Dd.to_string())
    rt = sorted(set(Dd[Dd.decision == "retemplate"].netsta))
    if rt != ["NC.KMR"]:
        log(f"[epochs] retemplate set {rt} != ['NC.KMR'] -- the queue's segment plan assumes KMR only; stop for a decision")
        sys.exit(3)
    k = Dd[(Dd.netsta == "NC.KMR") & (Dd.decision == "retemplate")].epoch_start.sort_values().dt.floor("D").tolist()
    if k != [KMR_SW1, KMR_SW2]:
        log(f"[epochs] KMR epoch starts {k} differ from the hard-coded switch days"); sys.exit(3)


# ------------------------------------------------------------------------------------------------ days
def his_days(path, y0="2024-01-01", y1="2024-12-01"):
    dd = json.load(open(path)); out = {}
    for k, v in dd.items():
        if k.split("|")[0] == SRC:
            out[k.split("|")[1]] = sorted(x for x in v if y0 <= x < y1)
    return out


def days():
    os.makedirs(DAYS, exist_ok=True)
    te = his_days(f"{RES}/aa_shelly/stage2_days_test.json"); bu = his_days(f"{RES}/aa_shelly/stage2_days_build.json")
    tu = set().union(*te.values()); bu_all = set().union(*his_days(f"{RES}/aa_shelly/stage2_days_build.json", "2018-01-01").values())
    cont = days_json(f"{DAYS}/cont.json", tu)
    strict = sorted(tu - bu_all)
    days_json(f"{DAYS}/cont_strict.json", strict)
    json.dump({f"{SRC}|{f}": v for f, v in te.items()}, open(f"{DAYS}/cont_fam.json", "w"))    # each family's own 2024 test days
    A = [d for d in rng(T_JOIN0, KMR_SW1) if str(d.date()) not in set(cont)]
    days_json(f"{DAYS}/A.json", A)
    days_json(f"{DAYS}/NK.json", rng(KMR_SW1, T_END))
    days_json(f"{DAYS}/B.json", rng(KMR_SW1 + pd.Timedelta(days=1), KMR_SW2))
    days_json(f"{DAYS}/C.json", rng(KMR_SW2 + pd.Timedelta(days=1), T_END))
    days_json(f"{DAYS}/all.json", rng(T_JOIN0, T_END))
    days_json(f"{DAYS}/W0.json", rng(T_JOIN0, T_HIS1))
    days_json(f"{DAYS}/W1.json", rng(T_HIS1, KMR_SW1))
    days_json(f"{DAYS}/W4.json", rng(pd.Timestamp("2026-07-02"), T_END))
    log(f"[days] cont {len(cont)} (strict held-out {len(strict)}), A {len(A)}, NK {len(rng(KMR_SW1, T_END))}, "
        f"B {len(rng(KMR_SW1 + pd.Timedelta(days=1), KMR_SW2))}, C {len(rng(KMR_SW2 + pd.Timedelta(days=1), T_END))}, LAST {LAST.date()}")


def read_days(path):
    dd = json.load(open(path)); s = set()
    for k, v in dd.items():
        if k.split("|")[0] == SRC: s |= set(v)
    return sorted(s)


# ------------------------------------------------------------------------------------------------ plan
def plan():
    D = pd.read_csv(f"{EXT}/station_epochs.csv", parse_dates=["start", "end"])
    rows = []
    for ns in STAS:
        g = D[D.netsta == ns]; net, sta = ns.split(".")
        for d in rng(T_JOIN0, T_END):
            if ((g.start < d + pd.Timedelta(days=1)) & (g.end > d)).any():
                rows.append(dict(net=net, sta=sta, day=str(d.date())))
    P = pd.DataFrame(rows); P.to_csv(f"{EXT}/download_plan.csv", index=False)
    have = sum(os.path.exists(npz_path(r.net, r.sta, r.day)) for r in P.itertuples())
    log(f"[plan] {len(P)} station-days ({have} already cached): {P.groupby('sta').size().to_dict()}")


# ------------------------------------------------------------------------------------------------ prep
_G = {}
def _prep_job(d):
    warnings.filterwarnings("ignore")
    import sh_v2_common as C2
    return C2.prep_day(pd.Timestamp(d), _G["stas"], _G["C"], _G["nref"])


def prep(days_file, workers, stas=None, outdir=None, summary=None):
    import sh_v2_common as C2
    from concurrent.futures import ProcessPoolExecutor
    stas = stas or STAS; outdir = outdir or PREPD
    os.makedirs(outdir, exist_ok=True)
    C2.PREP = outdir                                                 # prep_day writes here (workers inherit via fork)
    _G.update(stas=stas, C=C2.coords(), nref=int(round(6.0 * C2.FS)))
    miss = [s for s in stas if s not in _G["C"]]
    if miss: raise SystemExit(f"no coordinates for {miss}")
    dl = [d for d in read_days(days_file) if not os.path.exists(f"{outdir}/{d}.npz")]
    log(f"[prep] {len(dl)} days to prep ({days_file}) -> {outdir}, {len(stas)} stations {stas}, workers={workers}")
    rows = []
    with ProcessPoolExecutor(workers) as ex:
        for i, r in enumerate(ex.map(_prep_job, dl, chunksize=2)):
            rows += r
            if i % 100 == 0: log(f"  {i + 1}/{len(dl)}")
    p = summary or f"{EXT}/prep_summary.csv"; df = pd.DataFrame(rows)
    if os.path.exists(p): df = pd.concat([pd.read_csv(p), df]).drop_duplicates(["day", "netsta"], keep="last")
    df.to_csv(p, index=False)
    if len(df):
        log(f"[prep] station-days {df.groupby('netsta').size().to_dict()}; station-dead {df.groupby('netsta').station_dead.sum().to_dict()}")


# ------------------------------------------------------------------------------------------------ s1cat
def s1cat(inp, out, days_file="", k=9.0):
    from sh_v3_common import apply_rules
    D = pd.read_parquet(inp); D["fam"] = D.fam.astype(str)
    if days_file:
        ds = set(pd.to_datetime(read_days(days_file)))
        D = D[D.t.dt.floor("D").isin(ds)]
    E = apply_rules(D.reset_index(drop=True), k)
    E.to_parquet(out); log(f"[s1cat] {inp} -> {len(E)} S1 {k:g}x detections -> {out}; per family {E.groupby('fam').size().to_dict()}")


# ------------------------------------------------------------------------------------------------ tplset
def tcc(a, b, lag=5):
    """channel-mean of the max normalised cc over +-lag samples between two (C, L) templates."""
    out = []
    for c in range(min(a.shape[0], b.shape[0])):
        x, y = a[c].astype(float), b[c].astype(float)
        if not np.any(x) or not np.any(y): continue
        best = -1
        for s in range(-lag, lag + 1):
            u = x[max(0, s):len(x) + min(0, s)]; v = y[max(0, -s):len(y) + min(0, -s)]
            u = u - u.mean(); v = v - v.mean()
            best = max(best, float((u * v).sum() / (np.sqrt((u * u).sum() * (v * v).sum()) + 1e-30)))
        out.append(best)
    return float(np.mean(out)) if out else np.nan


def tplset(stage):
    M = pd.read_csv(TPL0 + ".csv"); Z = dict(np.load(TPL0 + ".npz"))
    if stage == "nokmr":
        m = M[M.netsta != "NC.KMR"]
        m.to_csv(f"{EXT}/tpl_noKMR.csv", index=False); np.savez_compressed(f"{EXT}/tpl_noKMR.npz", **{k: Z[k] for k in m.key})
        log(f"[tplset] noKMR {len(m)} templates"); return
    rows = []
    for ep in ("B", "C"):
        R = pd.read_csv(f"{EXT}/tpl_kmr{ep}_all.csv"); R["fam"] = R.fam.astype(str)
        ZR = dict(np.load(f"{EXT}/tpl_kmr{ep}.npz")) if os.path.exists(f"{EXT}/tpl_kmr{ep}.npz") else {}
        m = M[M.netsta != "NC.KMR"].copy(); T = {k: Z[k] for k in m.key}; add = []
        for r in R.itertuples():
            n = int(r.n_stack) if r.status == "ok" and np.isfinite(r.n_stack) else 0
            use = r.status == "ok" and n >= MIN_N_RETPL and r.key in ZR
            old = Z.get(r.key)
            rows.append(dict(epoch=ep, fam=r.fam, status=r.status, n_stack=n, cc50=getattr(r, "cc50", np.nan),
                             n_cand=getattr(r, "n_cand", np.nan), used=bool(use),
                             cc_new_vs_old=tcc(ZR[r.key], old) if (r.key in ZR and old is not None) else np.nan))
            if use:
                T[r.key] = ZR[r.key]; add.append(r.Index)
        R2 = R.loc[add, [c for c in M.columns if c in R.columns]]
        m = pd.concat([m, R2], ignore_index=True)
        m.to_csv(f"{EXT}/tpl_{ep}.csv", index=False); np.savez_compressed(f"{EXT}/tpl_{ep}.npz", **T)
        log(f"[tplset] {ep}: {len(m)} templates, KMR re-stacked for {len(R2)}/27 families")
    C = pd.DataFrame(rows); C.to_csv(f"{EXT}/tpl_kmr_compare.csv", index=False)
    log("[tplset] KMR re-template vs old STS-2 template:\n" + C.groupby("epoch").agg(
        used=("used", "sum"), n_med=("n_stack", "median"), cc_new_old_med=("cc_new_vs_old", "median"),
        cc_new_old_min=("cc_new_vs_old", "min")).to_string())


# ------------------------------------------------------------------------------------------------ continuity
def cont():
    from sh_v3_common import his_catalogue, apply_rules, score_all, family_days
    S = his_catalogue(); T = S.rename(columns={"meancc": "cc"})[["fam", "t", "cc"]]
    strict = set(pd.to_datetime(read_days(f"{DAYS}/cont_strict.json")))
    famd = family_days(f"{DAYS}/cont_fam.json")
    allc = set(pd.to_datetime(read_days(f"{DAYS}/cont.json")))
    src = {"ext (this run)": f"{EXT}/mf_cont.parquet", "v3b (same days)": f"{V3B}/mf_it0_test.parquet"}
    rows = []
    for name, p in src.items():
        O = pd.read_parquet(p); O["fam"] = O.fam.astype(str)
        O = O[O.t.dt.floor("D").isin(allc)].reset_index(drop=True)
        for k in (8, 9, 10):
            E = apply_rules(O, k)
            for mode in ("strict-source", "fam"):
                if mode == "strict-source":
                    Om = E[E.t.dt.floor("D").isin(strict)]; Tm = T[T.t.dt.floor("D").isin(strict)]
                else:
                    dO = E.t.dt.floor("D"); dT = T.t.dt.floor("D")
                    Om = E[[f in famd and d in famd[f] for f, d in zip(E.fam, dO)]]
                    Tm = T[[f in famd and d in famd[f] for f, d in zip(T.fam, dT)]]
                F, g = score_all(Om.reset_index(drop=True), Tm.reset_index(drop=True))
                rows.append(dict(run=name, xmad=k, mode=mode, days=len(strict) if mode == "strict-source" else sum(map(len, famd.values())),
                                 ours=g["ours"], his=g["his"], recall_fam=g["match_fam"] / max(g["his"], 1), prec_fam=g["match_fam"] / max(g["ours"], 1),
                                 recall_any=g["match_any"] / max(g["his"], 1), prec_any=g["match_any"] / max(g["ours"], 1),
                                 same_fam_share=g["same_fam"] / max(g["match_any"], 1)))
    R = pd.DataFrame(rows); R.to_csv(f"{EXT}/continuity.csv", index=False)
    log("[cont]\n" + R.round(3).to_string())
    a = R[(R.xmad == 9) & (R["mode"] == "strict-source")].set_index("run")
    e, v = a.loc["ext (this run)"], a.loc["v3b (same days)"]
    ok = (e.recall_any >= v.recall_any - 0.10) and (e.prec_any >= v.prec_any - 0.10)
    msg = (f"{'PASS' if ok else 'FAIL'}: strict held-out 2024 days, S1 9x any-family recall {e.recall_any:.3f} (v3b {v.recall_any:.3f}), "
           f"precision {e.prec_any:.3f} (v3b {v.prec_any:.3f}); rule = within 0.10 of v3b")
    open(f"{EXT}/CONT_GATE.txt", "w").write(msg + "\n"); log("[cont] " + msg)
    if not ok: sys.exit(4)


# ------------------------------------------------------------------------------------------------ catalog
SEGS = [("cont", "mf_cont", None), ("A", "mf_A", None), ("NK", "mf_NK", [KMR_SW1, KMR_SW2]), ("B", "mf_B", None), ("C", "mf_C", None)]


def load_segments():
    Ds, Gs = [], []
    for seg, f, only in SEGS:
        p = f"{EXT}/{f}.parquet"
        D = pd.read_parquet(p); D["fam"] = D.fam.astype(str); D["seg"] = seg
        G = pd.read_csv(p.replace(".parquet", "_diag.csv")); G["fam"] = G.fam.astype(str); G["day"] = pd.to_datetime(G.day); G["seg"] = seg
        if only is not None:
            D = D[D.t.dt.floor("D").isin(only)]; G = G[G.day.isin(only)]
        Ds.append(D); Gs.append(G)
    D = pd.concat(Ds, ignore_index=True); G = pd.concat(Gs, ignore_index=True)
    dup = G.groupby(["fam", "day"]).size()
    if (dup > 1).any(): raise SystemExit(f"overlapping segment days: {dup[dup > 1].head()}")
    return D.sort_values("t").reset_index(drop=True), G


def catalog():
    from sh_v3_common import apply_rules, his_catalogue, one_to_one_pairs
    D, G = load_segments()
    D = D[(D.t >= T_JOIN0) & (D.t < T_END)].reset_index(drop=True)
    cols = ["t", "fam", "ccsum", "meancc", "nlive", "net_mad", "nchan", "nsta", "seg"]
    R8 = D[D.net_mad >= 8.0][cols].reset_index(drop=True); R8.to_parquet(f"{EXT}/candidates_raw8x.parquet")
    CAT = apply_rules(D, 9)[cols].reset_index(drop=True)
    CAT.to_parquet(f"{EXT}/cat_S1_9x_internal.parquet")          # original column names (restack det-mode input)
    CAT = CAT.rename(columns={"nlive": "n_channels", "net_mad": "xMAD", "ccsum": "summed_cc", "meancc": "mean_cc"})
    CAT.to_parquet(f"{EXT}/catalogue_S1_9x.parquet"); CAT.to_csv(f"{EXT}/catalogue_S1_9x.csv", index=False, float_format="%.4f")
    os.makedirs(f"{EXT}/per_family", exist_ok=True)
    for f, g in CAT.groupby("fam"):
        g.to_csv(f"{EXT}/per_family/fam{int(f):02d}_detections.csv", index=False, float_format="%.4f")
    # days processed (a family-day exists in the MF diag) and live stations
    Gd = G[(G.day >= T_JOIN0) & (G.day < T_END)]
    live = Gd.groupby("day").agg(nsta=("nsta", "max"), nchan=("nchan", "max")).reset_index()
    live.to_csv(f"{EXT}/days_processed.csv", index=False)
    CAT["day"] = CAT.t.dt.floor("D")
    dc = CAT.groupby(["day", "fam"]).size().unstack(fill_value=0).reindex(columns=FAMS, fill_value=0)
    dc = dc.reindex(pd.to_datetime(live.day), fill_value=0); dc.index.name = "day"
    dc.insert(0, "nsta_live", live.set_index("day").nsta.values)
    dc.to_csv(f"{EXT}/daily_counts.csv")
    # monthly series
    S = his_catalogue(); S = S[S.t < T_HIS1]
    S["month"] = S.t.dt.to_period("M"); CAT["month"] = CAT.t.dt.to_period("M")
    months = pd.period_range("2018-01", LAST.to_period("M"), freq="M")
    ndays_ours = live.assign(month=pd.to_datetime(live.day).dt.to_period("M")).groupby("month").size()
    hisN = S.groupby(["month", "fam"]).size(); ourN = CAT.groupby(["month", "fam"]).size()
    rows = []
    for f in FAMS:
        for m in months:
            dim = m.days_in_month; nd = int(ndays_ours.get(m, 0))
            hn = int(hisN.get((m, f), 0)) if m < T_HIS1.to_period("M") else np.nan
            on = int(ourN.get((m, f), 0)) if nd else np.nan
            rows.append(dict(month=str(m), fam=f, his_n=hn, his_rate=hn / dim if np.isfinite(hn) else np.nan,
                             ours_n=on, ours_days=nd, ours_rate=on / nd if nd else np.nan))
    Mo = pd.DataFrame(rows)
    J = Mo[(pd.PeriodIndex(Mo.month, freq="M") >= T_JOIN0.to_period("M")) & (pd.PeriodIndex(Mo.month, freq="M") < T_HIS1.to_period("M"))]
    jt = J.groupby("fam").agg(his_n=("his_n", "sum"), ours_n=("ours_n", "sum"), ours_days=("ours_days", "sum")).reset_index()
    jt["his_rate"] = jt.his_n / 366.0; jt["ours_rate"] = jt.ours_n / jt.ours_days
    jt["ratio_ours_over_his"] = jt.ours_rate / jt.his_rate
    jt["r_monthly"] = [np.corrcoef(g.his_rate, g.ours_rate)[0, 1] if g.ours_rate.notna().sum() > 3 else np.nan
                       for _, g in J.groupby("fam")]
    jt = jt.sort_values("fam", key=lambda s: s.astype(int))
    sc = dict(zip(jt.fam, 1.0 / jt.ratio_ours_over_his.replace(0, np.nan)))
    Mo["ours_rate_his_scaled"] = Mo.ours_rate * Mo.fam.map(sc)
    Mo["combined_rate"] = np.where(pd.PeriodIndex(Mo.month, freq="M") < T_HIS1.to_period("M"), Mo.his_rate, Mo.ours_rate)
    Mo["combined_source"] = np.where(pd.PeriodIndex(Mo.month, freq="M") < T_HIS1.to_period("M"), "Shelly2025", "ours_S1_9x")
    Mo.to_csv(f"{EXT}/combined_monthly.csv", index=False)
    # event matching over the join year (all days; includes template-build days)
    Tj = S[(S.t >= T_JOIN0) & (S.t < T_HIS1)]
    Oj = CAT[(CAT.t < T_HIS1)]
    ot = Oj.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9; tt = Tj.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    P = one_to_one_pairs(ot, tt, 2.0)
    same = sum(Oj.fam.values[i] == Tj.fam.values[k] for i, k in P)
    tot = dict(fam="all", his_n=int(jt.his_n.sum()), ours_n=int(jt.ours_n.sum()), ours_days=int(ndays_ours[(ndays_ours.index >= T_JOIN0.to_period("M")) & (ndays_ours.index < T_HIS1.to_period("M"))].sum()))
    tot.update(his_rate=tot["his_n"] / 366.0, ours_rate=tot["ours_n"] / max(tot["ours_days"], 1))
    tot["ratio_ours_over_his"] = tot["ours_rate"] / tot["his_rate"]
    jt = pd.concat([jt, pd.DataFrame([tot])], ignore_index=True)
    jt.to_csv(f"{EXT}/join_comparison.csv", index=False)
    json.dump(dict(join_match_any=len(P), join_his=len(Tj), join_ours=len(Oj), join_same_family=int(same),
                   join_recall_any=len(P) / max(len(Tj), 1), join_precision_any=len(P) / max(len(Oj), 1),
                   n_cat=len(CAT), n_raw8=len(R8), days_processed=len(live), first_day=str(live.day.min())[:10], last_day=str(live.day.max())[:10]),
              open(f"{EXT}/catalog_summary.json", "w"), indent=1)
    log(f"[catalog] S1 9x {len(CAT)} detections, raw 8x {len(R8)}, {len(live)} days; join year any-fam recall {len(P) / max(len(Tj), 1):.3f} "
        f"precision {len(P) / max(len(Oj), 1):.3f}; per-family ratio ours/his median {np.nanmedian(jt.ratio_ours_over_his[:-1]):.2f}")


# ------------------------------------------------------------------------------------------------ tplcheck
def tplcheck():
    rows = []
    Z0 = dict(np.load(TPL0 + ".npz")); ZC = dict(np.load(f"{EXT}/tpl_C.npz")) if os.path.exists(f"{EXT}/tpl_C.npz") else {}
    for w, desc in (("W0", "2023-12..2024-11 (his period, our detections)"), ("W1", "2024-12..2025-09-16"), ("W4", "2026-07-02..end (KMPB split)")):
        p = f"{EXT}/tplcheck_{w}_all.csv"
        if not os.path.exists(p): continue
        R = pd.read_csv(p); Zn = dict(np.load(f"{EXT}/tplcheck_{w}.npz"))
        for r in R[R.status == "ok"].itertuples():
            ref = ZC.get(r.key) if (w == "W4" and r.netsta == "NC.KMR" and r.key in ZC) else Z0.get(r.key)
            rows.append(dict(window=w, desc=desc, fam=str(r.fam), netsta=r.netsta, n_stack=int(r.n_stack),
                             cc_vs_template=tcc(Zn[r.key], ref) if ref is not None else np.nan))
    T = pd.DataFrame(rows); T.to_csv(f"{EXT}/tplcheck.csv", index=False)
    if len(T):
        S = T.groupby(["netsta", "window"]).cc_vs_template.median().unstack()
        S.to_csv(f"{EXT}/tplcheck_station.csv"); log("[tplcheck] median cc (fresh stack vs template in use) per station x window\n" + S.round(3).to_string())


# ------------------------------------------------------------------------------------------------ report
def fig_monthly(Mo, live, out):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HIS, OURS = "#2a78d6", "#eb6834"
    Mo = Mo.copy(); Mo["tm"] = pd.PeriodIndex(Mo.month, freq="M").to_timestamp() + pd.Timedelta(days=14)
    fig, axs = plt.subplots(7, 4, figsize=(15, 17), sharex=True)
    panels = FAMS + ["all"]
    for ax, f in zip(axs.flat, panels):
        g = Mo[Mo.fam == f] if f != "all" else Mo.groupby("tm", as_index=False).agg(his_rate=("his_rate", lambda s: s.sum(min_count=1)),
                                                                                    ours_rate=("ours_rate", lambda s: s.sum(min_count=1)))
        ax.plot(g.tm, g.his_rate, color=HIS, lw=1.3, label="Shelly 2025 catalogue")
        ax.plot(g.tm, g.ours_rate, color=OURS, lw=1.3, label="ours (S1 rules, 9×MAD)")
        ax.axvline(T_HIS1, color="0.25", ls="--", lw=1.0)
        for t in (KMR_SW1, KMR_SW2): ax.axvline(t, color="0.6", ls=":", lw=0.8)
        ax.set_title(f"family {f}" if f != "all" else "all 27 families", fontsize=9, loc="left")
        ax.grid(alpha=0.25, lw=0.5); ax.tick_params(labelsize=7)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    axs.flat[0].legend(fontsize=7, frameon=False, loc="upper left")
    for ax in axs[:, 0]: ax.set_ylabel("events / day", fontsize=8)
    fig.suptitle("Mendocino LFE families: monthly detection rate 2018 to 2026. Dashed line = his catalogue end / our extension start "
                 "(2024-12-01); dotted = NC.KMR sensor changes", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.98)); fig.savefig(out, dpi=110); plt.close(fig)


def fig_total(Mo, live, out):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HIS, OURS = "#2a78d6", "#eb6834"
    Mo = Mo.copy(); Mo["tm"] = pd.PeriodIndex(Mo.month, freq="M").to_timestamp() + pd.Timedelta(days=14)
    g = Mo.groupby("tm", as_index=False).agg(his=("his_rate", lambda s: s.sum(min_count=1)), ours=("ours_rate", lambda s: s.sum(min_count=1)))
    g = g[g.tm >= pd.Timestamp("2022-01-01")]
    fig, (a, b) = plt.subplots(2, 1, figsize=(11, 6), sharex=True, gridspec_kw=dict(height_ratios=(2.2, 1)))
    a.plot(g.tm, g.his, color=HIS, lw=1.6, label="Shelly 2025 catalogue"); a.plot(g.tm, g.ours, color=OURS, lw=1.6, label="ours (S1 rules, 9×MAD)")
    a.axvline(T_HIS1, color="0.25", ls="--", lw=1); a.set_ylabel("events / day (27 families)"); a.legend(frameon=False, fontsize=8)
    L = live.copy(); L["day"] = pd.to_datetime(L.day)
    b.plot(L.day, L.nsta, color="0.35", lw=0.8); b.set_ylabel("live S1 stations"); b.set_ylim(0, 18)
    b.axvline(T_HIS1, color="0.25", ls="--", lw=1)
    for ax in (a, b):
        for t in (KMR_SW1, KMR_SW2): ax.axvline(t, color="0.6", ls=":", lw=0.8)
        ax.grid(alpha=0.25, lw=0.5)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(out, dpi=120); plt.close(fig)


def md(df, fmt="%.3f"):
    try:
        return df.to_markdown(index=False, floatfmt=fmt.replace("%", ""))
    except Exception:
        return "```\n" + df.to_string(index=False) + "\n```"


def report():
    Mo = pd.read_csv(f"{EXT}/combined_monthly.csv", dtype={"fam": str}); live = pd.read_csv(f"{EXT}/days_processed.csv")
    os.makedirs(f"{EXT}/figures", exist_ok=True)
    fig_monthly(Mo, live, f"{EXT}/figures/monthly_rate_per_family.png"); fig_total(Mo, live, f"{EXT}/figures/monthly_rate_total.png")
    C = pd.read_csv(f"{EXT}/continuity.csv"); J = pd.read_csv(f"{EXT}/join_comparison.csv", dtype={"fam": str})
    Sm = json.load(open(f"{EXT}/catalog_summary.json")); Dd = pd.read_csv(f"{EXT}/epoch_decisions.csv")
    K = pd.read_csv(f"{EXT}/tpl_kmr_compare.csv") if os.path.exists(f"{EXT}/tpl_kmr_compare.csv") else pd.DataFrame()
    TC = pd.read_csv(f"{EXT}/tplcheck_station.csv") if os.path.exists(f"{EXT}/tplcheck_station.csv") else pd.DataFrame()
    CAT = pd.read_parquet(f"{EXT}/catalogue_S1_9x.parquet")
    gate = open(f"{EXT}/CONT_GATE.txt").read().strip()
    c9 = C[(C.xmad == 9)].set_index(["run", "mode"])
    e = c9.loc[("ext (this run)", "strict-source")]
    # rates: extension period vs his full record, per family
    mo = Mo.copy(); mo["m"] = pd.PeriodIndex(mo.month, freq="M")
    his_full = mo[mo.m < T_HIS1.to_period("M")].groupby("fam").his_n.sum() / (T_HIS1 - pd.Timestamp("2018-01-01")).days
    ext = mo[mo.m >= T_HIS1.to_period("M")]
    ext_rate = ext.groupby("fam").ours_n.sum() / ext.groupby("fam").ours_days.sum()
    jr = J[J.fam != "all"].set_index("fam")
    PF = pd.DataFrame(dict(his_rate_2018_2024=his_full, ours_rate_join_year=jr.ours_rate, his_rate_join_year=jr.his_rate,
                           ratio_join=jr.ratio_ours_over_his, r_monthly_join=jr.r_monthly, ours_rate_extension=ext_rate))
    PF["ext_his_equiv"] = PF.ours_rate_extension / PF.ratio_join
    PF["ext_vs_his_2018_2024"] = PF.ext_his_equiv / PF.his_rate_2018_2024
    PF = PF.reset_index().rename(columns={"index": "fam"}).sort_values("fam", key=lambda s: s.astype(int))
    PF.to_csv(f"{EXT}/per_family_rates.csv", index=False)
    jall = J[J.fam == "all"].iloc[0]
    n_ext = int((CAT.t >= T_HIS1).sum()); ext_days = int((pd.to_datetime(live.day) >= T_HIS1).sum())
    kmr_used = K.groupby("epoch").used.sum().to_dict() if len(K) else {}
    kmr_cc = K.groupby("epoch").cc_new_vs_old.median().round(2).to_dict() if len(K) else {}
    up = int((PF.ext_vs_his_2018_2024 > 1.25).sum()); dn = int((PF.ext_vs_his_2018_2024 < 0.8).sum())
    L = ["# SH_EXT2026: Shelly (2025) Mendocino LFE families extended past 2024-11-30", "",
         f"(auto-generated by `scripts/sh_ext2026.py report`, {pd.Timestamp.now('UTC'):%Y-%m-%d %H:%MZ}; queue `scripts/sh_ext2026_queue.sh`)", "",
         "## Plain-language summary", "",
         f"- We ran our reproduction of Shelly's detector (the v3b 'it0' templates on his 17 stations, with his catalogue rules) on every day from "
         f"{Sm['first_day']} to {Sm['last_day']} ({Sm['days_processed']} days). The first 12 months overlap his catalogue, so the two can be "
         f"compared directly before we go past his end date.",
         f"- **Still reproduces him.** On 2024 days held out from template building, our run finds {e.recall_any:.0%} of his events (any family) "
         f"and {e.prec_any:.0%} of our detections are events he also has. The earlier v3b validation on all held-out 2018-2024 days gave 68 % and 78 %. ({gate})",
         f"- **Rates at the join.** Over Dec 2023 to Nov 2024 his catalogue has {int(jall.his_n)} events and ours has {int(jall.ours_n)} "
         f"(ratio {jall.ratio_ours_over_his:.2f}). Event by event over that year, we recover {Sm['join_recall_any']:.0%} of his events, and "
         f"{Sm['join_precision_any']:.0%} of ours are his. That year includes template-build days, so these numbers are a little optimistic. "
         f"Per family, the monthly rate curves correlate with a median r = {np.nanmedian(jr.r_monthly):.2f}.",
         f"- **Beyond his catalogue.** From 2024-12-01 to {Sm['last_day']}, we have {n_ext} catalogue-grade detections on {ext_days} days. "
         f"We rescaled each family's rate by its ours/his ratio at the join. On that basis, {up} families run more than 25 % above their 2018-2024 average "
         f"and {dn} run more than 20 % below it (table in section 4; see the monthly rate figure below).",
         f"- **Station changes.** Only one station changed hardware after his catalogue ends. NC.KMR swapped its STS-2 for a Trillium Horizon 120 "
         f"on 2025-09-17, then a Trillium Compact Posthole on 2026-03-12. For each new sensor, we built KMR's templates again from the 100 "
         f"strongest catalogue detections made without KMR, using the Text S1 recipe. Families with enough events: {kmr_used}. The median cc "
         f"between the new and old KMR templates is {kmr_cc}. KMR is left out on the two changeover days. NC.KMPB's 2026-07-01 epoch is a "
         f"metadata-only split (identical sensor, gain, orientation and coordinates), so its template is unchanged.",
         "- Precision caveat: validated precision is about 0.78 (any family) at 9×MAD. Roughly one detection in five is not an event he would list. "
         "Family labels agree with his for about 70 % of matched events (source mode). Treat per-family counts as our detector's rates, "
         "scaled to his only through the join-year ratio.", "",
         "![monthly rate per family](sh_ext2026/figures/monthly_rate_per_family.png)", "",
         "![network total + live stations](sh_ext2026/figures/monthly_rate_total.png)", "",
         "## 1. Station epochs (FDSN) and how they were handled", "", md(Dd.assign(epoch_start=Dd.epoch_start.astype(str).str[:16], epoch_end=Dd.epoch_end.astype(str).str[:16])), ""]
    if len(K):
        L += ["NC.KMR re-templates (restack det mode: the 100 strongest S1 9× detections from the run WITHOUT KMR in that epoch, "
              "±0.1 s shifts, best 50; S offset kept from it0). Used only if n_stack ≥ 20.", "",
              md(K.groupby("epoch").agg(families_used=("used", "sum"), n_stack_median=("n_stack", "median"), n_stack_min=("n_stack", "min"),
                                        cc_new_vs_old_median=("cc_new_vs_old", "median"), cc_new_vs_old_min=("cc_new_vs_old", "min")).reset_index()), ""]
    if len(TC):
        L += ["Template stability: median (over families) cc between a fresh stack of our own S1 detections in each window and the template "
              "in use (W0 = 2023-12..2024-11 inside his period; W1 = 2024-12..2025-09-16; W4 = 2026-07-02..end, after the KMPB epoch split).", "",
              md(TC), ""]
    L += ["## 2. Continuity check (2024, inside his catalogue)", "",
          "Detector re-run with this run's prep (17-station cleaning) on the 100 Shelly-2025 test days of 2024. strict-source = the 34 of them that "
          "are no family's build day, all his events; fam = each family on its own held-out days. v3b = the validated run on the same days.", "",
          md(C), "", gate, "",
          "## 3. The join: his catalogue vs ours, 2023-12..2024-11", "", md(J), "",
          f"Event matching over the join year (±2 s one-to-one, any family): recall {Sm['join_recall_any']:.3f}, precision {Sm['join_precision_any']:.3f}, "
          f"same-family share {Sm['join_same_family'] / max(Sm['join_match_any'], 1):.2f}.", "",
          "## 4. Per-family rates", "",
          "his_rate_2018_2024 = his events/day over 2018-01..2024-11; ours_rate_extension = our S1 9× events per processed day from 2024-12-01; "
          "ext_his_equiv = that divided by the join-year ratio; ext_vs_his_2018_2024 = ext_his_equiv / his_rate_2018_2024.", "", md(PF), "",
          "## 5. Files (results/sh_ext2026/)", "",
          "- `catalogue_S1_9x.{parquet,csv}`: t, fam, summed_cc, mean_cc, n_channels (live channels at the peak), xMAD, nchan, nsta, seg; "
          "`per_family/famNN_detections.csv`",
          "- `candidates_raw8x.parquet`: every ≥8×MAD 4-s-bin maximum before the S1 rules",
          "- `daily_counts.csv`: day × family counts (processed days only, with live-station count); `days_processed.csv`",
          "- `combined_monthly.csv`: month × family: his_n/his_rate (2018-01..2024-11), ours_n/ours_rate (2023-12..), ours_rate_his_scaled, combined_rate",
          "- `join_comparison.csv`, `continuity.csv`, `per_family_rates.csv`, `station_epochs.csv`, `epoch_decisions.csv`, `tpl_kmr_compare.csv`, `tplcheck*.csv`",
          "- MF segments `mf_{cont,A,NK,B,C}.parquet` (+ `_diag.csv`): cont/A = it0; NK = it0 without KMR (used for 2025-09-17 and 2026-03-12 and "
          "to re-template KMR); B = it0 + KMR re-stacked for the Horizon epoch; C = it0 + KMR re-stacked for the Compact Posthole epoch", "",
          "## 6. Deviations", "",
          "- The download and continuous detection start on 2023-12-01, not 2024-12-01, because the join comparison needs 12 months of continuous "
          "overlap with his catalogue (v3b only has build and test days cached).",
          "- Prep (burst muting neighbours) uses the 17 S1 stations only. v3b used the 23-station v2 set for its old stations. The continuity check "
          "measures the effect on the same days.",
          "- The KMR re-templates keep the it0 S offset (no new S pick), as in the v3 'it1' stage.",
          "- Everything else is the v3b it0 recipe unchanged (t2_fetch, t3_cache, sh_v2_common.prep_day, t11_mf_v2, sh_v3_common.apply_rules).", ""]
    bs = f"{RES}/sh_back/BACK_SECTION.md"
    L += ["<!-- BACK_BEGIN -->", open(bs).read().strip() if os.path.exists(bs) else "## Backward extension\n\n(pending: results/sh_back, queue scripts/sh_back_queue.sh)", "<!-- BACK_END -->", ""]
    out = os.environ.get("SH_EXT_REPORT", f"{RES}/SH_EXT2026_RESULTS.md")
    import fcntl
    with open(f"{WS}/locks/sh_results_md.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX); open(out, "w").write("\n".join(L))
    log(f"[report] wrote {out}")


def datasum():
    P = pd.read_csv(f"{EXT}/download_plan.csv")
    M = pd.read_csv(f"{EXT}/download_manifest.csv").drop_duplicates(["net", "sta", "day"], keep="last")
    C = pd.read_parquet(f"{EXT}/cache_ext.parquet")
    rows = []
    for (n, s), g in P.groupby(["net", "sta"]):
        ns = f"{n}.{s}"; m = M[(M.net == n) & (M.sta == s)]; c = C[C.netsta == ns]
        rows.append(dict(station=ns, planned=len(g), data=int((m.status == "data").sum()), nodata=int((m.status == "nodata").sum()),
                         errors=int((m.status == "error").sum()), GB=round(m.bytes.sum() / 1e9, 2),
                         cached=int(c.status.isin(["built", "cached"]).sum()), band=",".join(sorted(c.band.dropna().unique())) if "band" in c else ""))
    R = pd.DataFrame(rows); R.to_csv(f"{EXT}/data_summary.csv", index=False)
    new = M[M.note != "resumed"]
    log(f"[datasum] total miniSEED {M.bytes.sum() / 1e9:.1f} GB ({new.bytes.sum() / 1e9:.1f} GB newly fetched)\n" + R.to_string(index=False))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what")
    ap.add_argument("--days-file", default=""); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--inp", default=""); ap.add_argument("--out", default=""); ap.add_argument("--stage", default="")
    ap.add_argument("--stas", default="", help="prep: comma list (default the 17 S1 stations)")
    ap.add_argument("--prep-dir", default=""); ap.add_argument("--summary", default="")
    a = ap.parse_args(); os.makedirs(EXT, exist_ok=True); warnings.filterwarnings("ignore")
    dict(epochs=epochs, days=days, plan=plan, cont=cont, catalog=catalog, tplcheck=tplcheck, report=report, datasum=datasum,
         prep=lambda: prep(a.days_file, a.workers, a.stas.split(",") if a.stas else None, a.prep_dir or None, a.summary or None), s1cat=lambda: s1cat(a.inp, a.out, a.days_file),
         tplset=lambda: tplset(a.stage))[a.what]()
