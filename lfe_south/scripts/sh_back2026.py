#!/usr/bin/env python
"""SH_BACK2026 (2026-09-30): can Shelly (2025)'s 27 Mendocino families be followed back before 2018?

  epochs   full FDSN epoch history (1990..) of the 17 S1 stations + his location-stage stations (NC.KBS KBSB KKP KIP KCPB KCPC
           KBN KFP KCR KSM, PB.B045 B046 B047 B935) -> back_station_summary.csv.  A detection station is usable back to the date
           from which it has run the SAME SENSOR MODEL as during the template build (2018-2024); gain-only changes do not
           matter for a normalised cc.  S6 = S1 stations usable before 2018; S2 = the two usable in 2003-2008 (KMPB, KMR).
  tplsub   it0 templates restricted to a station subset.
  BK1      (queue) the v3b it0 detector + Text S1 rules on the SAME 385 held-out 2018-2024 test days, with ONLY the S6 (and S2)
           stations, scored exactly like v3b by sh_v3_rescore.py.
  gate     BK2 runs only if BK1 (S6) at 9xMAD, S1 rules, source mode gives any-family precision >= 0.5.
  plan2    BK2 download plan: S6, 2008-01-01..2017-12-31, inside each station's same-sensor period.
  catalog2 S1 9x catalogue 2008-2017, per-family daily counts, monthly rates (+ his-equivalent via the BK1 per-family ratio).
  report   results/sh_back/BACK_SECTION.md, inserted into results/SH_EXT2026_RESULTS.md between BACK markers.
"""
import os, sys, json, argparse, warnings, subprocess
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t_common import RES, WS, log, npz_path

BK = os.environ.get("SH_BACK_DIR", f"{RES}/sh_back"); EXT = f"{RES}/sh_ext2026"; SRC = "Shelly2025"; FAMS = [str(i) for i in range(1, 28)]
S1 = ["BK.BLAS", "BK.BRIC", "BK.GHOP", "BK.HALS", "BK.JEWT", "BK.PRDS", "BK.SPRL", "BK.RVIT", "BK.HUNT",
      "NC.KMPB", "NC.KMR", "NC.KPR", "NC.KBU", "PB.B049", "PB.B932", "PB.B934", "PB.B933"]
LOC = ["NC.KBS", "NC.KBSB", "NC.KKP", "NC.KIP", "NC.KCPB", "NC.KCPC", "NC.KBN", "NC.KFP", "NC.KCR", "NC.KSM",
       "PB.B045", "PB.B046", "PB.B047", "PB.B935"]
SETS = {"S6": ["NC.KMPB", "NC.KMR", "PB.B049", "PB.B932", "PB.B933", "PB.B934"], "S2": ["NC.KMPB", "NC.KMR"]}
TPL0 = f"{RES}/sh_v3b/tpl_it0"
B0, B1 = pd.Timestamp("2008-01-01"), pd.Timestamp("2018-01-01")
GATE_PREC = 0.5
MB_PER_STADAY = 26.0      # v3b fetch: 42.6 GB / 1639 station-days


def model(s):
    return str(s).split(",")[0].strip()


def epochs():
    raw = f"{EXT}/station_epochs_full_raw.csv"
    if not os.path.exists(raw):
        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        rows = []
        for ns in S1 + LOC:
            n, s = ns.split(".")
            for prov in (["NCEDC", "EARTHSCOPE"] if n in ("BK", "NC") else ["EARTHSCOPE", "NCEDC"]):
                try:
                    inv = Client(prov, timeout=120).get_stations(network=n, station=s, channel="HH?,BH?,EH?,SH?,HN?", level="response",
                                                                 starttime=UTCDateTime("1990-01-01"))
                except Exception as e:
                    log(f"  {ns} {prov} ERR {str(e)[:80]}"); continue
                for net in inv:
                    for st in net:
                        for ch in st:
                            sens = ch.response.instrument_sensitivity
                            rows.append(dict(prov=prov, netsta=ns, loc=ch.location_code, cha=ch.code, sr=ch.sample_rate,
                                             start=str(ch.start_date)[:19], end=str(ch.end_date)[:19] if ch.end_date else "",
                                             sensor=(ch.sensor.description or ch.sensor.model or ch.sensor.type or "") if ch.sensor else "",
                                             gain=sens.value if sens else None, az=ch.azimuth, dip=ch.dip, lat=ch.latitude, lon=ch.longitude))
                break
        pd.DataFrame(rows).to_csv(raw, index=False)
    D = pd.read_csv(raw, dtype={"loc": str}); D["loc"] = D["loc"].fillna("")
    D = D[~D.cha.str.startswith("HN")].copy(); D["band"] = D.cha.str[:2]
    en = D.end.fillna("").astype(str).str[:19]
    D["t0"] = pd.to_datetime(D.start, format="ISO8601")
    D["t1"] = pd.to_datetime(en.where(~(en.eq("") | en.str.startswith("3000") | en.eq("nan")), "2100-01-01T00:00:00"), format="ISO8601")
    rows, eps = [], []
    for ns in S1 + LOC:
        g = D[D.netsta == ns]
        if not len(g):
            rows.append(dict(netsta=ns, role="S1 detection" if ns in S1 else "location stage", note="no FDSN metadata")); continue
        E = (g.groupby(["band", "loc", "t0"]).agg(t1=("t1", "first"), cha=("cha", lambda s: ",".join(sorted(s))), sr=("sr", "first"),
                                                  sensor=("sensor", "first"), gain=("gain", "first"), lat=("lat", "first"), lon=("lon", "first"))
             .reset_index().sort_values("t0"))
        E["netsta"] = ns; eps.append(E)
        # the band t3_cache would use on each date: HH > BH > EH > SH (>= 20 sps)
        E = E[E.sr >= 20]
        mid = pd.Timestamp("2021-06-01")
        def at(t):
            h = E[(E.t0 <= t) & (E.t1 > t)]
            if not len(h): return None
            for b in ("HH", "BH", "EH", "SH"):
                if (h.band == b).any(): return h[h.band == b].iloc[0]
        # transitions of what t3_cache would read (preferred band, sensor model, rate), sampled every 10 days
        ch = []; prev = "start"
        for t in pd.date_range(E.t0.min().normalize(), pd.Timestamp("2026-09-30"), freq="10D"):
            e = at(t); key = "off" if e is None else f"{e.band} {model(e.sensor)} {e.sr:g}sps"
            if key != prev: ch.append(f"{t:%Y-%m} {key}")
            prev = key
        cur = at(mid)
        same_since = None
        if cur is not None:
            hb = E[(E.band == cur.band) & (E.sensor.map(model) == model(cur.sensor))].sort_values("t0")
            # walk back while contiguous (gap < 30 d) and same model
            s0 = cur.t0
            for _, e in hb[hb.t0 < cur.t0].iloc[::-1].iterrows():
                if e.t1 >= s0 - pd.Timedelta(days=30): s0 = min(s0, e.t0)
                else: break
            same_since = s0
        rows.append(dict(netsta=ns, role="S1 detection" if ns in S1 else "location stage",
                         first_start=str(E.t0.min())[:10], last_end=str(E.t1.max())[:10],
                         build_sensor=(f"{cur.band} {model(cur.sensor)}" if cur is not None else "started after 2021-06"),
                         same_sensor_since=str(same_since)[:10] if same_since is not None else "",
                         operating_before_2018=bool(E.t0.min() < B1),
                         usable_before_2018=bool(ns in S1 and same_since is not None and same_since < B1),
                         band_sensor_changes="; ".join(ch)))
    R = pd.DataFrame(rows); os.makedirs(BK, exist_ok=True); R.to_csv(f"{BK}/back_station_summary.csv", index=False)
    pd.concat(eps).to_csv(f"{BK}/back_station_epochs.csv", index=False)
    log("[back epochs]\n" + R.drop(columns=["band_sensor_changes"]).to_string())
    use = sorted(R[R.usable_before_2018].netsta)
    if use != sorted(SETS["S6"]):
        log(f"[back epochs] usable-before-2018 set {use} != S6 {SETS['S6']} -- stop for a decision"); sys.exit(3)


def tplsub(name):
    M = pd.read_csv(TPL0 + ".csv"); Z = dict(np.load(TPL0 + ".npz"))
    m = M[M.netsta.isin(SETS[name])]
    m.to_csv(f"{BK}/tpl_{name}.csv", index=False); np.savez_compressed(f"{BK}/tpl_{name}.npz", **{k: Z[k] for k in m.key})
    log(f"[tplsub] {name}: {len(m)} templates on {m.netsta.nunique()} stations")


def days2():
    ds = [str(d.date()) for d in pd.date_range(B0, B1, freq="D", inclusive="left")]
    json.dump({f"{SRC}|{f}": ds for f in FAMS}, open(f"{BK}/days_bk2.json", "w"))


def bk1_table():
    rows = []
    for name, p in (("17 stations (v3b)", f"{RES}/sh_v3b/rescore_it0.csv"), ("S6 (pre-2018 stations)", f"{BK}/rescore_bk1_S6.csv"),
                    ("S2 (KMPB+KMR)", f"{BK}/rescore_bk1_S2.csv")):
        if not os.path.exists(p): continue
        R = pd.read_csv(p); R = R[(R.variant == "S1") & R.xmad.isin([8, 9, 10])]
        for r in R.itertuples():
            rows.append(dict(stations=name, xMAD=r.xmad, mode=r.mode, ours=r.ours, his=r.his, recall_own=r.recall_fam, precision_own=r.prec_fam,
                             recall_any=r.recall_any, precision_any=r.prec_any, same_fam_share=r.same_fam / max(r.match_any, 1)))
    return pd.DataFrame(rows)


def gate():
    T = bk1_table()
    r = T[(T.stations.str.startswith("S6")) & (T.xMAD == 9) & (T["mode"] == "source")].iloc[0]
    ok = r.precision_any >= GATE_PREC
    msg = (f"{'PASS' if ok else 'FAIL'}: BK1 S6 at 9xMAD (S1 rules, source mode) any-family precision {r.precision_any:.3f}, recall {r.recall_any:.3f}; "
           f"rule: BK2 runs only if precision >= {GATE_PREC}")
    open(f"{BK}/BK1_GATE.txt", "w").write(msg + "\n"); log("[gate] " + msg)


def plan2():
    E = pd.read_csv(f"{BK}/back_station_summary.csv")
    rows = []; have = 0
    for ns in SETS["S6"]:
        s0 = pd.Timestamp(E[E.netsta == ns].same_sensor_since.iloc[0])
        first = max(B0, s0.normalize() + pd.Timedelta(days=1) if s0 != s0.normalize() else s0)
        net, sta = ns.split(".")
        for d in pd.date_range(first, B1, freq="D", inclusive="left"):
            rows.append(dict(net=net, sta=sta, day=str(d.date())))
            have += os.path.exists(npz_path(net, sta, d))
    P = pd.DataFrame(rows); P.to_csv(f"{BK}/download_plan_bk2.csv", index=False)
    new = len(P) - have
    df = subprocess.run(["df", "-h", "/home/jovyan"], capture_output=True, text=True).stdout.strip().splitlines()[-1]
    msg = (f"BK2 plan: {len(P)} station-days ({have} already cached, {new} new) -> est. {new * MB_PER_STADAY / 1e3:.0f} GB miniSEED + "
           f"~{new * 0.8 * MB_PER_STADAY / 1e3:.0f} GB 20-sps cache; df -h /home/jovyan: {df}")
    open(f"{BK}/BK2_VOLUME.txt", "w").write(msg + "\n" + P.groupby("sta").day.agg(["size", "min"]).to_string() + "\n")
    log("[plan2] " + msg)


def catalog2():
    from sh_v3_common import apply_rules
    D = pd.read_parquet(f"{BK}/mf_bk2.parquet"); D["fam"] = D.fam.astype(str)
    G = pd.read_csv(f"{BK}/mf_bk2_diag.csv"); G["day"] = pd.to_datetime(G.day)
    cols = ["t", "fam", "ccsum", "meancc", "nlive", "net_mad", "nchan", "nsta"]
    D[D.net_mad >= 8][cols].to_parquet(f"{BK}/candidates_bk2_raw8x.parquet")
    C = apply_rules(D, 9)[cols].rename(columns={"nlive": "n_channels", "net_mad": "xMAD", "ccsum": "summed_cc", "meancc": "mean_cc"})
    C.to_parquet(f"{BK}/catalogue_bk2_S1_9x.parquet"); C.to_csv(f"{BK}/catalogue_bk2_S1_9x.csv", index=False, float_format="%.4f")
    live = G.groupby("day").agg(nsta=("nsta", "max")).reset_index(); live.to_csv(f"{BK}/days_processed_bk2.csv", index=False)
    C["day"] = C.t.dt.floor("D")
    dc = C.groupby(["day", "fam"]).size().unstack(fill_value=0).reindex(columns=FAMS, fill_value=0).reindex(live.day, fill_value=0)
    dc.index.name = "day"; dc.insert(0, "nsta_live", live.nsta.values); dc.to_csv(f"{BK}/daily_counts_bk2.csv")
    # BK1 per-family calibration: S6 S1-9x detections / his events on the 2018-2024 test days (source mode)
    F = pd.read_csv(f"{BK}/rescore_bk1_S6_fam.csv"); F = F[(F.variant == "S1") & (F.xmad == 9) & (F["mode"] == "source")]
    F["fam"] = F.fam.astype(str); ratio = (F.set_index("fam").n_ours / F.set_index("fam").n_theirs.replace(0, np.nan)).to_dict()
    C["month"] = C.t.dt.to_period("M"); nd = live.assign(m=live.day.dt.to_period("M")).groupby("m").size()
    rows = []
    for f in FAMS:
        for m in pd.period_range("2008-01", "2017-12", freq="M"):
            n = int(((C.fam == f) & (C.month == m)).sum()); d = int(nd.get(m, 0))
            r = n / d if d else np.nan
            rows.append(dict(month=str(m), fam=f, n=n, days=d, rate=r, ratio_bk1=ratio.get(f, np.nan),
                             rate_his_equiv=r / ratio.get(f, np.nan) if d else np.nan))
    pd.DataFrame(rows).to_csv(f"{BK}/monthly_bk2.csv", index=False)
    log(f"[catalog2] {len(C)} S1 9x detections 2008-2017 on {len(live)} days")


def fig_back(out):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HIS, OURS, BACK = "#2a78d6", "#eb6834", "#1baf7a"
    B = pd.read_csv(f"{BK}/monthly_bk2.csv", dtype={"fam": str}); B["tm"] = pd.PeriodIndex(B.month, freq="M").to_timestamp() + pd.Timedelta(days=14)
    Fw = pd.read_csv(f"{EXT}/combined_monthly.csv", dtype={"fam": str}) if os.path.exists(f"{EXT}/combined_monthly.csv") else None
    from sh_v3_common import his_catalogue
    S = his_catalogue(); S = S[S.t < pd.Timestamp("2024-12-01")]; S["m"] = S.t.dt.to_period("M")
    fig, axs = plt.subplots(7, 4, figsize=(15, 17), sharex=True)
    for ax, f in zip(axs.flat, FAMS + ["all"]):
        if f == "all":
            b = B.groupby("tm", as_index=False).rate_his_equiv.sum(min_count=1); hs = S.groupby("m").size()
        else:
            b = B[B.fam == f]; hs = S[S.fam == f].groupby("m").size()
        hs = hs.reindex(pd.period_range("2018-01", "2024-11", freq="M"), fill_value=0)
        ax.plot(b.tm, b.rate_his_equiv, color=BACK, lw=1.1, label="ours, 6 pre-2018 stations (his-equivalent)")
        ax.plot(hs.index.to_timestamp() + pd.Timedelta(days=14), hs.values / hs.index.days_in_month, color=HIS, lw=1.1, label="Shelly 2025 catalogue")
        if Fw is not None:
            w = Fw[pd.PeriodIndex(Fw.month, freq="M") >= pd.Period("2024-12", "M")]
            w = w.groupby("month", as_index=False).ours_rate_his_scaled.sum(min_count=1) if f == "all" else w[w.fam == f]
            ax.plot(pd.PeriodIndex(w.month, freq="M").to_timestamp() + pd.Timedelta(days=14), w.ours_rate_his_scaled, color=OURS, lw=1.1,
                    label="ours, 17 stations (his-equivalent)")
        for t in ("2018-01-01", "2024-12-01"): ax.axvline(pd.Timestamp(t), color="0.3", ls="--", lw=0.9)
        ax.set_title(f"family {f}" if f != "all" else "all 27 families", fontsize=9, loc="left"); ax.grid(alpha=0.25, lw=0.5); ax.tick_params(labelsize=7)
        for s in ("top", "right"): ax.spines[s].set_visible(False)
    axs.flat[0].legend(fontsize=6.5, frameon=False, loc="upper left")
    for ax in axs[:, 0]: ax.set_ylabel("events / day", fontsize=8)
    fig.suptitle("Mendocino LFE families 2008 to 2026, monthly rate. Our series are divided by each family's ours/his count ratio "
                 "(BK1 test days for 2008-2017, join year for 2024-12 on)", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.98)); fig.savefig(out, dpi=110); plt.close(fig)


def md(df):
    try: return df.to_markdown(index=False, floatfmt=".3f")
    except Exception: return "```\n" + df.to_string(index=False) + "\n```"


def report():
    S = pd.read_csv(f"{BK}/back_station_summary.csv")
    T = bk1_table()
    g = open(f"{BK}/BK1_GATE.txt").read().strip() if os.path.exists(f"{BK}/BK1_GATE.txt") else "(BK1 not scored yet)"
    L = ["## Backward extension (before 2018)", "",
         "Question: can the 27 families be followed back before 2018 with the stations that existed then?",
         "A station counts for a period only if it ran the same sensor model then as during the 2018-2024 template build. "
         "Gain-only changes do not affect a normalised cc; other stations are dropped.", ""]
    if len(T):
        v = T[(T.xMAD == 9) & (T["mode"] == "source")].set_index("stations")
        L += ["**Plain summary.** Only 6 of the 17 S1 stations were running before 2018: NC.KMPB (CMG-3ESPC since 2008-04-29), NC.KMR "
              "(STS-2 since 2007-01-30) and PB.B049/B932/B933/B934 (borehole HS-1, installed Sep-Nov 2008). We re-ran the detector with only "
              "these stations on the same held-out 2018-2024 days used for the v3b score. "
              + "; ".join(f"{k}: any-family recall {r.recall_any:.2f}, precision {r.precision_any:.2f}" for k, r in v.iterrows())
              + f" (9×MAD, Text S1 rules). {g}", ""]
    L += ["### Stations: full FDSN epoch history", "",
          md(S[["netsta", "role", "first_start", "last_end", "build_sensor", "same_sensor_since", "operating_before_2018", "usable_before_2018",
                "band_sensor_changes"]]), "",
          "### BK1: known-answer test (385 held-out 2018-2024 test days, scored exactly as v3b)", "",
          md(T) if len(T) else "(pending)", "", g, ""]
    if os.path.exists(f"{BK}/BK2_VOLUME.txt"):
        L += ["### BK2: 2008-2017 scan", "", "```", open(f"{BK}/BK2_VOLUME.txt").read().strip(), "```", ""]
    if os.path.exists(f"{BK}/monthly_bk2.csv"):
        fig_back(f"{BK}/monthly_rate_2008_2026.png")
        M = pd.read_csv(f"{BK}/monthly_bk2.csv", dtype={"fam": str}); C = pd.read_parquet(f"{BK}/catalogue_bk2_S1_9x.parquet")
        live = pd.read_csv(f"{BK}/days_processed_bk2.csv")
        yr = C.assign(y=C.t.dt.year).groupby(["y"]).size()
        pf = M.groupby("fam").agg(n=("n", "sum"), days=("days", "sum"), ratio_bk1=("ratio_bk1", "first")).reset_index()
        pf["rate"] = pf.n / pf.days.replace(0, np.nan); pf["rate_his_equiv"] = pf.rate / pf.ratio_bk1
        pf = pf.sort_values("fam", key=lambda s: s.astype(int))
        L += [f"{len(C)} S1 9× detections 2008-2017 on {len(live)} processed days (live stations per day in `days_processed_bk2.csv`).",
              "Per year: " + ", ".join(f"{y}: {n}" for y, n in yr.items()), "",
              "![2008-2026](sh_back/monthly_rate_2008_2026.png)", "", md(pf), "",
              "Files (results/sh_back/): `catalogue_bk2_S1_9x.{parquet,csv}`, `candidates_bk2_raw8x.parquet`, `daily_counts_bk2.csv`, `monthly_bk2.csv`.", ""]
    elif "FAIL" in g:
        L += ["BK2 not run (gate failed).", ""]
    s = "\n".join(L); open(f"{BK}/BACK_SECTION.md", "w").write(s)
    out = f"{RES}/SH_EXT2026_RESULTS.md"
    import fcntl
    with open(f"{WS}/locks/sh_results_md.lock", "w") as lk:
        fcntl.flock(lk, fcntl.LOCK_EX)
        if os.path.exists(out):
            r = open(out).read(); a, b = "<!-- BACK_BEGIN -->", "<!-- BACK_END -->"
            r = (r[:r.index(a) + len(a)] + "\n" + s + "\n" + r[r.index(b):]) if (a in r and b in r) else r.rstrip() + f"\n\n{a}\n{s}\n{b}\n"
            open(out, "w").write(r)
    log(f"[back report] wrote {BK}/BACK_SECTION.md" + (" and updated SH_EXT2026_RESULTS.md" if os.path.exists(out) else ""))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what"); ap.add_argument("--set", default="S6")
    a = ap.parse_args(); os.makedirs(BK, exist_ok=True); warnings.filterwarnings("ignore")
    dict(epochs=epochs, tplsub=lambda: tplsub(a.set), days2=days2, gate=gate, plan2=plan2, catalog2=catalog2, report=report)[a.what]()
