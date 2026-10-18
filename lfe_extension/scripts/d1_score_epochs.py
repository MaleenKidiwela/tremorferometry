#!/usr/bin/env python
"""Track D, D1(a) -- EPOCH / STATION-SET ABLATION  (Merlin Step 6, component 1).

Per-family own recall of Bostock's events at T_core and T_ext for EVERY station set the 2014-2026
production run actually realises, measured on the referee year (2012, test half) and on the held-out
year (2010).  The station sets come from Track C's epoch policy
(`results/trackC/templates_by_epoch.csv`) and `results/station_epochs.csv`:

  E1  2014-01-01 .. 2019-05-29   15 stations  (the frozen 2012 network)
  E2  2019-05-30 .. 2021-03-27   14 stations  (LRIV dropped by policy: A-2 step 4 identity FAIL)
  E3  2021-03-28 .. 2026-09-15   13 stations  (LRIV dropped + SNB ends 2021-03-28)

CN.PFB (ends 2018-12-20) carries NO 2012 template and is not in the network, so the 2018-12 boundary of
the brief has no effect on this detector; it is recorded in the output as a no-op epoch note.

Matching: own family, +-2 s, one-to-one greedy by time (Track A step 0(b)).  The frozen HF veto is
applied exactly as in production (PGC 10-19 Hz / 2-8 Hz, cut 3.6384 in 2012/2010 = the BH epochs; an
undefined ratio is never vetoed)."""
from __future__ import annotations
import argparse, json, os, sys
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_tiers  # noqa: E402  (per-epoch T_core / T_ext, Track C's FA-only re-derivation)

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
OUT = f"{EXT}/results"
TD = f"{EXT}/results/trackD"

# (station-set id, tag, referee year, stations, note).  The axis is the STATION SET, which is what the
# ablation measures; Track C's network epochs are mapped onto it by `d_tiers` (two of them share S15).
RUNS = [
    ("S15", "a1_mergeB_L12", 2012, "PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926",
     "the frozen 15-station network (Track C epochs E0 and E1)"),
    ("S14", "d1_abl2012_E2", 2012, "PGC,VGZ,NLLB,SNB,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926",
     "LRIV dropped from 2019-05-30, A-2 epoch-identity FAIL (Track C epoch E2)"),
    ("S13", "d1_abl2012_E3", 2012, "PGC,VGZ,NLLB,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926",
     "LRIV dropped + CN.SNB ends 2021-03-28 (Track C epoch E3)"),
    ("P9", "a2_abl2012_7p", 2012, "PGC,VGZ,SNB,NLLB,B926,LRIV,B009,B010,B011",
     "the 9-station persistent set = the 2010 network (cross-year control, Track A step 2)"),
    ("S15_2010", "a2_xfer2010", 2010, "PGC,VGZ,SNB,NLLB,B926,LRIV,B009,B010,B011",
     "held-out year, full 2010 network (only these 9 stations are cached for 2010)"),
    ("S14_2010", "d1_abl2010_E2", 2010, "PGC,VGZ,SNB,NLLB,B926,B009,B010,B011", "held-out year, LRIV dropped"),
    ("S13_2010", "d1_abl2010_E3", 2010, "PGC,VGZ,NLLB,B926,B009,B010,B011", "held-out year, LRIV + SNB dropped"),
]
SET_OF_RUN = {"S15": "S15", "S14": "S14", "S13": "S13", "P9": "S15",
              "S15_2010": "S15", "S14_2010": "S14", "S13_2010": "S13"}
HF_CUT = 3.6384


def one2one(a, b_sorted, tol=2.0):
    a = np.asarray(a, float); out = np.full(len(a), np.inf)
    if len(a) == 0 or len(b_sorted) == 0:
        return out
    j = 0; nb = len(b_sorted)
    for i in np.argsort(a, kind="stable"):
        ti = a[i]
        while j < nb and b_sorted[j] < ti - tol:
            j += 1
        if j < nb and abs(b_sorted[j] - ti) <= tol:
            out[i] = abs(b_sorted[j] - ti); j += 1
    return out


QUIET = {2012: list(range(228, 243)), 2010: [223, 224, 226, 252, 254, 255, 256]}


def fa_curve(det, famdays, quiet, grid):
    """quiet-day detections per template-day vs threshold (the same FA definition Track A froze on)."""
    nd = max(int(famdays.doy.isin(quiet).sum()), 1)
    q = det[det.doy.isin(quiet)]
    return np.array([len(q[q.net_mad >= T]) / nd for T in grid])


def t_at_fa(grid, fa, target):
    """log-interpolate the threshold that gives `target` false alarms per template-day."""
    f = np.log(np.maximum(fa, 1e-6))
    t = np.log(max(target, 1e-6))
    if f[0] <= t:
        return float(grid[0])
    if f[-1] >= t:
        return float(grid[-1])
    i = int(np.searchsorted(-f, -t))
    i = min(max(i, 1), len(grid) - 1)
    w = (f[i - 1] - t) / max(f[i - 1] - f[i], 1e-9)
    return float(grid[i - 1] + w * (grid[i] - grid[i - 1]))


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def famdays(tag, L=12.0):
    b = pd.read_csv(f"{OUT}/e5_budget_{tag}.csv"); b = b[b.L == L]
    T0 = sorted(b["T"].unique())[0]
    return b[b["T"] == T0][["doy", "fam"]].drop_duplicates()


def load_det(tag, L=12.0, hf=True):
    det = pd.read_parquet(f"{OUT}/e5r_detections_{tag}.parquet")
    det = det[det.L == L].reset_index(drop=True)
    n0 = len(det); n_vetoed = 0
    hp = f"{OUT}/a1_hf_ratio_{tag}.parquet"
    if hf and os.path.exists(hp):
        H = pd.read_parquet(hp); H = H[H.L == L][["doy", "fam", "t", "hfr_PGC"]]
        det = det.merge(H, on=["doy", "fam", "t"], how="left")
        v = det.hfr_PGC.to_numpy(); drop = np.isfinite(v) & (v > HF_CUT)
        n_vetoed = int(drop.sum()); det = det[~drop].reset_index(drop=True)
    elif hf:
        raise SystemExit(f"[d1] missing HF parquet {hp} -- run scripts/d1_run_hf.sh first")
    return det, n0, n_vetoed


def gt_events(year, fams_present):
    if year == 2012:
        g = pd.read_parquet(f"{IN}/e5r_gt_bostock2012.parquet").reset_index(drop=True)
        ev = g[(g.half == "test") & (~g.guard)].copy().reset_index(drop=True)
        ev["bkey"] = ev.bkey.astype(int)
        return ev
    g = pd.read_parquet(f"{OUT}/a2_gt_bostock2010.parquet")
    off = pd.read_csv(f"{IN}/e3_bostock_family_offsets.csv"); off = off[off.accept.astype(bool)]
    omap = dict(zip(off.bkey.astype(int), off.offset_s.astype(float)))
    g = g[g.bkey.astype(int).isin(omap)].copy()
    g["time"] = g.time.astype("datetime64[ns]")
    g["t_origin"] = g.time.astype("int64") / 1e9 + g.bkey.astype(int).map(omap).to_numpy()
    g["bkey"] = g.bkey.astype(int)
    return g.reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--t-ext", type=float, default=np.nan)
    ap.add_argument("--t-core", type=float, default=np.nan)
    ap.add_argument("--no-hf", action="store_true")
    a = ap.parse_args()
    os.makedirs(TD, exist_ok=True)
    EPOCHS, TIERS, TSRC = d_tiers.load()
    # The STATION-SET ablation must use ONE threshold for every set, or "what does losing a station
    # cost" is confounded with "what threshold does that epoch run at".  It therefore uses the FROZEN
    # referee point throughout; the per-EPOCH table further down re-reads recall at each epoch's own
    # Track C-derived thresholds, and that is what the catalog merges on.
    TH_SET = {ss: dict(d_tiers.FROZEN) for ss in ("S15", "S14", "S13")}
    if np.isfinite(a.t_ext) or np.isfinite(a.t_core):      # explicit override, one T everywhere
        TH_SET = {ss: dict(T_ext=a.t_ext if np.isfinite(a.t_ext) else 8.25,
                           T_core=a.t_core if np.isfinite(a.t_core) else 8.75) for ss in TH_SET}
        TSRC = "command-line override"
    print(f"[tiers] {TSRC}\n        " + "; ".join(f"{k} T_core {v['T_core']} T_ext {v['T_ext']}"
                                                   for k, v in TH_SET.items()))
    # template coverage per family (which stations carry a 2012 template)
    qc = pd.read_csv(f"{IN}/e5r_template_qc_2012.csv")
    ucol = "usable" if "usable" in qc.columns else None
    if ucol:
        qc = qc[qc[ucol].astype(bool)]
    tpl = qc.groupby("bkey").sta.apply(lambda s: ",".join(sorted(set(s)))).to_dict()
    ntpl = qc.groupby("bkey").sta.nunique().to_dict()

    rows = []; meta = []; mrows = []
    grid = np.arange(7.0, 12.01, 0.25)
    fa_base = {}
    for eid, tag, year, stas, note in RUNS:
        p = f"{OUT}/e5r_detections_{tag}.parquet"
        if not os.path.exists(p):
            print(f"[skip] {eid}: {p} not on disk yet"); continue
        det, n0, nv = load_det(tag, hf=not a.no_hf)
        ev = gt_events(year, set(det.fam.unique()))
        stal = stas.split(",")
        FDy = famdays(tag)
        TH0 = TH_SET[SET_OF_RUN[eid]]
        fa = fa_curve(det, FDy, QUIET[year], grid)
        meta.append(dict(run_id=eid, station_set=SET_OF_RUN[eid], tag=tag, ref_year=year,
                         n_stations=len(stal), stations=stas,
                         n_detections_recorded=n0, n_hf_vetoed=nv, note=note,
                         T_ext=TH0["T_ext"], T_core=TH0["T_core"],
                         quiet_fa_at_T_ext=float(np.interp(TH0["T_ext"], grid, fa)),
                         quiet_fa_at_T_core=float(np.interp(TH0["T_core"], grid, fa))))
        if eid in ("S15", "S15_2010"):
            fa_base[year] = dict(T_ext=float(np.interp(TH0["T_ext"], grid, fa)),
                                 T_core=float(np.interp(TH0["T_core"], grid, fa)))
        TH = TH_SET[SET_OF_RUN[eid]]
        # ---- MATCHED-FALSE-ALARM leg (Merlin F3): at a FIXED MAD threshold a smaller network is not a
        # fairer detector, it is a LOWER threshold.  Re-read recall at the T that reproduces the E1
        # quiet-day false-alarm rate of the same year.
        for nm, T in (("T_ext", TH["T_ext"]), ("T_core", TH["T_core"])):
            tgt = fa_base[year][nm]
            Tm = t_at_fa(grid, fa, tgt)
            Dm = det[det.net_mad >= Tm]
            mrows.append(dict(run_id=eid, station_set=SET_OF_RUN[eid], ref_year=year, tier=nm,
                              T_frozen=T, T_matched=Tm,
                              fa_target=tgt, fa_at_frozen=float(np.interp(T, grid, fa)),
                              n_det_matched=int(len(Dm))))
        for nm, T in (("T_ext", TH["T_ext"]), ("T_core", TH["T_core"])):
            D = det[det.net_mad >= T]
            Tm = [m for m in mrows if m["run_id"] == eid and m["tier"] == nm][0]["T_matched"]
            Dm = det[det.net_mad >= Tm]
            nh_m = nt_m = 0
            for bk, g in ev.groupby("bkey"):
                td = np.sort(D[D.fam == float(bk)].t.to_numpy())
                hit = one2one(g.t_origin.to_numpy(), td) <= 2.0
                tdm = np.sort(Dm[Dm.fam == float(bk)].t.to_numpy())
                hitm = one2one(g.t_origin.to_numpy(), tdm) <= 2.0
                nh_m += int(hitm.sum()); nt_m += len(g)
                lo, hi = wilson(int(hit.sum()), len(g))
                rows.append(dict(run_id=eid, station_set=SET_OF_RUN[eid], ref_year=year, tier=nm,
                                 T=T, family_id=int(bk),
                                 recall_matched_fa=float(hitm.mean()), T_matched_fa=Tm,
                                 n_events_ref=int(len(g)), n_recalled=int(hit.sum()),
                                 recall=float(hit.mean()), recall_lo=lo, recall_hi=hi,
                                 median_mag=float(np.median(g.mag)) if "mag" in g else np.nan,
                                 n_stations=len(stal), n_template_sta=int(ntpl.get(bk, 0)),
                                 template_stations=tpl.get(bk, ""),
                                 template_less=bool(ntpl.get(bk, 0) == 0)))
        print(f"[{eid}] {tag} {year}: {len(det)} detections ({nv} HF-vetoed of {n0}), "
              f"{ev.bkey.nunique()} families, {len(ev)} reference events", flush=True)
    pd.DataFrame(mrows).to_csv(f"{TD}/completeness_matched_fa.csv", index=False)
    json.dump(dict(source=TSRC, epochs=EPOCHS, tiers=TIERS, thresholds_by_station_set=TH_SET),
              open(f"{TD}/tiers_used.json", "w"), indent=1, default=str)
    R = pd.DataFrame(rows)
    # ratio to the reference epoch of the same year
    base = {2012: "S15", 2010: "S15_2010"}
    key = ["ref_year", "tier", "family_id"]
    B = R[R.run_id.isin(base.values())][key + ["recall"]].rename(columns={"recall": "recall_base"})
    R = R.merge(B, on=key, how="left")
    R["recall_ratio_to_base"] = R.recall / R.recall_base.replace(0, np.nan)
    # ---- per-EPOCH scoring at that epoch's OWN re-derived thresholds.  E0 and E1 share the S15
    # station set but Track C now gives them different T (8.75/8.25 vs 8.25/8.00), so the station-set
    # table above (scored at the tightest T of the epochs that use the set) understates the looser
    # epoch.  This block re-reads recall for every epoch at its own numbers; it is the table the
    # catalog and the family table merge on.
    EX = []
    for e in EPOCHS:
        ss = e["station_set"]
        tag = next((t for i_, t, y_, st_, n_ in RUNS if i_ == ss and y_ == 2012), None)
        if tag is None or not os.path.exists(f"{OUT}/e5r_detections_{tag}.parquet"):
            continue
        det, _, _ = load_det(tag, hf=not a.no_hf)
        ev = gt_events(2012, set(det.fam.unique()))
        TH = TIERS.get(e["epoch_id"], d_tiers.FROZEN)
        for nm in ("T_ext", "T_core"):
            T = TH[nm]
            D = det[det.net_mad >= T]
            for bk, g in ev.groupby("bkey"):
                td = np.sort(D[D.fam == float(bk)].t.to_numpy())
                hit = one2one(g.t_origin.to_numpy(), td) <= 2.0
                lo, hi = wilson(int(hit.sum()), len(g))
                EX.append(dict(epoch_id=e["epoch_id"], epoch_t0=e["t0"], epoch_t1=e["t1"],
                               station_set=ss, n_stations_epoch=e["n_stations"], tier=nm, T=T,
                               ref_year=2012, family_id=int(bk), n_events_ref=int(len(g)),
                               n_recalled=int(hit.sum()), recall=float(hit.mean()),
                               recall_lo=lo, recall_hi=hi,
                               n_template_sta=int(ntpl.get(bk, 0)),
                               template_stations=tpl.get(bk, ""),
                               template_less=bool(ntpl.get(bk, 0) == 0)))
        print(f"[epoch {e['epoch_id']}] {ss} via {tag}: T_core {TH['T_core']} T_ext {TH['T_ext']}",
              flush=True)
    R2 = pd.DataFrame(EX) if EX else R.assign(epoch_id="")
    R2.to_csv(f"{TD}/completeness_family_epoch.csv", index=False)
    R.to_csv(f"{TD}/completeness_family_stationset.csv", index=False)
    pd.DataFrame(meta).to_csv(f"{TD}/completeness_epoch_runs.csv", index=False)

    # pooled summary + the cross-year validation of the STATION-SET term
    S = R.groupby(["run_id", "ref_year", "tier"]).apply(
        lambda g: pd.Series(dict(n_fam=len(g), n_events=int(g.n_events_ref.sum()),
                                 recall_pooled=g.n_recalled.sum() / max(g.n_events_ref.sum(), 1),
                                 recall_fam_median=float(g.recall.median()),
                                 recall_fam_p10=float(g.recall.quantile(.1)),
                                 recall_fam_p90=float(g.recall.quantile(.9)),
                                 recall_pooled_matched_fa=(g.recall_matched_fa * g.n_events_ref).sum()
                                 / max(g.n_events_ref.sum(), 1))), include_groups=False).reset_index()
    M = pd.read_csv(f"{TD}/completeness_matched_fa.csv")
    S = S.merge(M[["run_id", "tier", "T_matched", "fa_target", "fa_at_frozen"]],
                on=["run_id", "tier"], how="left")
    S.to_csv(f"{TD}/completeness_epoch_summary.csv", index=False)
    if "epoch_id" in R2.columns and len(R2):
        SE = R2.groupby(["epoch_id", "station_set", "tier"], as_index=False).agg(
            T=("T", "first"), n_fam=("family_id", "size"), n_events=("n_events_ref", "sum"),
            n_recalled=("n_recalled", "sum"), recall_fam_median=("recall", "median"),
            recall_fam_p10=("recall", lambda v: v.quantile(.1)),
            recall_fam_p90=("recall", lambda v: v.quantile(.9)))
        SE["recall_pooled"] = SE.n_recalled / SE.n_events.clip(lower=1)
        SE.to_csv(f"{TD}/completeness_by_epoch_summary.csv", index=False)
        print("\n[per-EPOCH recall, each at its own re-derived thresholds]")
        print(SE.round(4).to_string(index=False))
    pd.set_option("display.width", 220)
    print(S.round(4).to_string(index=False))
    # cross-year check: does dropping LRIV / LRIV+SNB cost the same fraction in 2012 and in 2010?
    val = []
    for tier in ("T_ext", "T_core"):
        s = S[S.tier == tier].set_index("run_id").recall_pooled
        for e12, e10, lab in (("S14", "S14_2010", "drop LRIV"), ("S13", "S13_2010", "drop LRIV+SNB")):
            if e12 in s and e10 in s:
                val.append(dict(tier=tier, ablation=lab,
                                ratio_2012=float(s[e12] / s["P9"]) if "P9" in s else np.nan,
                                ratio_2010=float(s[e10] / s["S15_2010"]),
                                ratio_2012_from15=float(s[e12] / s["S15"])))
    V = pd.DataFrame(val)
    if len(V):
        V["cross_year_delta"] = V.ratio_2010 - V.ratio_2012
        V.to_csv(f"{TD}/completeness_epoch_crossyear.csv", index=False)
        print("\n[cross-year station-set term] (2012 ratios are vs the 9-station P9 base, the 2010 network)")
        print(V.round(4).to_string(index=False))
    print(f"[write] {TD}/completeness_family_epoch.csv (+ _summary, _runs, _crossyear)")


if __name__ == "__main__":
    main()
