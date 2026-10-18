#!/usr/bin/env python
"""P0 -- the LABEL TABLE (results/labels_v0.parquet), Merlin's P0 specification.

POSITIVES     Bostock 2010-2013 events that are ON DISK (2008/2009 appended when their traces land),
              one row per (event, station with data that day and a 2012 template).
              S time  = DB time + FROZEN 2012 family offset + [ the a5 stack S pick where one exists,
                        else Crosson ttS + the a5 per-station correction ].
              provenance = picked | predicted; label sigma = 0.5 s for both (Merlin: far-borehole S
              labels are PREDICTED, not measured, and carry >= 0.5 s uncertainty).
              P rows ONLY where a stack P pick exists (a5_picks.t_P finite).
NEGATIVES     certain: PNSN-quiet station-days (no tremor in the family box that day AND >= 3 days
              from any tremor day), with +-120 s around EVERY Bostock event excluded by rule.
              hard:    catalogued earthquakes / blasts (Track C's eq_catalog.csv); HF-veto quiet-day
              transients; T_ext quiet-day detections of the frozen detector.
UNLABELLED    ETS-day windows that are not near a Bostock event -- recorded as a rule + a station-day
              inventory, NEVER as negatives (Merlin F6).
HOLDOUTS      the 30 hidden families and the 3 held-out stations are FLAGGED, never dropped.

  python scripts/o_labels.py
"""
from __future__ import annotations
import os, sys, json, glob, argparse
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import o_vel1d as V

HERE = os.path.dirname(os.path.abspath(__file__))
PICK = os.path.dirname(HERE)
RES = f"{PICK}/results"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
IN = f"{REP}/results"
BOST = f"{IN}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"
R = 111.195
SIGMA = 0.5
CACHES = [f"{PICK}/cache20hz", f"{REP}/cache20hz", f"{EXT}/cache20hz"]
WINDOWS = {2010: ("2010-08-10", "2010-09-14"), 2011: ("2011-08-05", "2011-09-26"),
           2012: ("2012-08-15", "2012-10-06"), 2013: ("2013-08-30", "2013-10-29"),
           2008: ("2008-05-01", "2008-06-11"), 2009: ("2009-05-01", "2009-06-16")}
DET_SRC = {2010: ("a2_xfer2010", EXT), 2012: ("a1_mergeB_L12", EXT)}
HF_CUT = 3.6384


def epi_km(la1, lo1, la2, lo2):
    return np.hypot((np.asarray(lo1) - lo2) * R * np.cos(np.radians(48.6)),
                    (np.asarray(la1) - la2) * R)


def cache_index():
    """{(sta, year, doy)} -> netsta, over every 20 Hz cache root on disk."""
    idx = {}
    for root in CACHES:
        for p in glob.glob(f"{root}/*/*/*.npz"):
            parts = p.split(os.sep)
            netsta, yr, fn = parts[-3], int(parts[-2]), parts[-1]
            idx.setdefault((netsta.split(".")[1], yr, int(fn[:3])), netsta)
    return idx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{RES}/labels_v0.parquet")
    a = ap.parse_args()

    HO = json.load(open(f"{RES}/holdouts.json"))
    HID = {int(h["bkey"]) for h in HO["hidden_families"]}
    HSTA = set(HO["held_out_stations"])
    FAM = pd.read_csv(f"{IN}/e5r_families_2012.csv").set_index("bkey")
    OFF = pd.read_csv(f"{IN}/e5r_family_offsets_2012.csv")
    omap = {int(r.bkey): float(r.offset_s) for r in OFF.itertuples() if bool(r.accept)}
    PK = pd.read_csv(f"{EXT}/results/a5_picks.csv")
    CORR = pd.read_csv(f"{EXT}/results/a5_station_corrections.csv").set_index("sta")
    STA = pd.read_csv(f"{IN}/e5r_stations.csv").set_index("sta")
    QC = pd.read_csv(f"{IN}/e5r_template_qc_2012.csv")
    HASTPL = {(int(r.bkey), str(r.sta)) for r in QC.itertuples() if r.usable and r.template}
    idx = cache_index()
    print(f"[cache] {len(idx)} station-days on disk")

    # ---------- per (family, station) S and P travel times + provenance ----------
    tt = {}
    for r in PK.itertuples():
        s = str(r.sta); b = int(r.bkey)
        cS = float(CORR.corr_S.get(s, np.nan)); cP = float(CORR.corr_P.get(s, np.nan))
        cS = 0.0 if not np.isfinite(cS) else cS
        cP = 0.0 if not np.isfinite(cP) else cP
        if np.isfinite(r.t_S):
            ttS, prov = float(r.t_S) + cS, "picked"          # a5 stores obs - corr; undo it
        else:
            ttS, prov = float(r.ttS) + cS, "predicted"
        ttP = float(r.t_P) + cP if np.isfinite(r.t_P) else np.nan
        tt[(b, s)] = (ttS, prov, ttP)
    # families x stations with a template but no a5 row -> pure Crosson + correction
    n_fill = 0
    for (b, s) in HASTPL:
        if (b, s) in tt or s not in STA.index or b not in FAM.index:
            continue
        d = epi_km(FAM.lat[b], FAM.lon[b], STA.lat[s], STA.lon[s])
        cS = float(CORR.corr_S.get(s, np.nan)); cS = 0.0 if not np.isfinite(cS) else cS
        tt[(b, s)] = (float(V.tt(d, float(FAM.dep[b]), "S")) + cS, "predicted", np.nan)
        n_fill += 1
    print(f"[tt] {len(tt)} (family, station) travel times "
          f"({sum(1 for v in tt.values() if v[1]=='picked')} picked, {n_fill} Crosson-filled)")

    # ---------- Bostock events, exact hypocentre -> bkey, frozen offsets ----------
    B = pd.read_csv(BOST)
    B["time"] = pd.to_datetime(dict(year=B.year, month=B.month, day=B.day, hour=B.hour,
                                    minute=B["min"], second=B.sec.astype(int)), errors="coerce") \
        + pd.to_timedelta(B.sec % 1.0, unit="s")
    key = {(round(FAM.lat[b], 4), round(FAM.lon[b], 4), round(FAM.dep[b], 4)): int(b) for b in FAM.index}
    B["bkey"] = [key.get((round(la, 4), round(lo, 4), round(dp, 4)), -1)
                 for la, lo, dp in zip(B.lat, B.lon, B.dep)]
    B = B[B.bkey >= 0].copy()
    B["t_origin"] = B.time.astype("datetime64[ns]").astype("int64") / 1e9 + B.bkey.map(omap)
    B = B.dropna(subset=["t_origin"])

    # ---------- PNSN-quiet days in the family box ----------
    lat0, lat1 = FAM.lat.min() - 0.2, FAM.lat.max() + 0.2
    lon0, lon1 = FAM.lon.min() - 0.2, FAM.lon.max() + 0.2
    TR = pd.read_csv(PNSN, usecols=["time", "lat", "lon"])
    TR["time"] = pd.to_datetime(TR.time, errors="coerce")
    TR = TR.dropna(subset=["time"])
    TR = TR[TR.lat.between(lat0, lat1) & TR.lon.between(lon0, lon1)]
    trem_days = set(TR.time.dt.normalize())
    # PNSN COVERAGE GUARD: the catalog starts 2010-01, so "no tremor detection that day" means
    # NOTHING for 2008/2009 -- without this guard every 2008/2009 day would be labelled a certain
    # negative, including the ETS days themselves (exactly Merlin's F6 failure).
    TR_T0, TR_T1 = TR.time.min().normalize(), TR.time.max().normalize()
    print(f"[pnsn] tremor catalog covers {TR_T0.date()} .. {TR_T1.date()} in the family box")

    rows = []
    inv = []
    for year, (w0, w1) in WINDOWS.items():
        days = pd.date_range(w0, pd.Timestamp(w1) - pd.Timedelta(days=1), freq="D")
        have = sorted({s for (s, y, d) in idx if y == year})
        if not have:
            inv.append(dict(year=year, n_days=len(days), n_stations=0, note="no cache on disk"))
            continue
        ev = B[(B.time >= w0) & (B.time < w1)]
        # quiet day = no tremor in the box that day AND >= 3 days from any tremor day
        ev_days = pd.DatetimeIndex(ev.time).normalize().value_counts()
        ets_days = set(ev_days[ev_days >= 100].index)
        covered = (TR_T0 <= pd.Timestamp(w0)) and (pd.Timestamp(w1) <= TR_T1 + pd.Timedelta(days=1))
        if covered:
            qd = {d for d in days
                  if not any((abs((d - t).days) < 3) for t in trem_days if abs((d - t).days) < 3)
                  and int(ev_days.get(d, 0)) == 0}
        else:
            qd = set()          # no tremor coverage -> NO certain negatives can be declared
        # ---- positives ----
        npos = 0
        for r in ev.itertuples():
            b = int(r.bkey); doy = r.time.dayofyear
            for s in have:
                if (s, year, doy) not in idx or (b, s) not in tt:
                    continue
                ttS, prov, ttP = tt[(b, s)]
                rows.append((("pos_S"), year, doy, r.time.strftime("%Y-%m-%d"),
                             idx[(s, year, doy)].split(".")[0], s, b,
                             float(r.t_origin + ttS), prov, SIGMA,
                             b in HID, s in HSTA, "bostock2015", float(r.mag)))
                npos += 1
                if np.isfinite(ttP):
                    rows.append((("pos_P"), year, doy, r.time.strftime("%Y-%m-%d"),
                                 idx[(s, year, doy)].split(".")[0], s, b,
                                 float(r.t_origin + ttP), "picked", SIGMA,
                                 b in HID, s in HSTA, "bostock2015+a5_stack_P", float(r.mag)))
        # ---- certain negatives: quiet station-days ----
        nneg = 0
        for d in days:
            if d not in qd:
                continue
            doy = d.dayofyear
            for s in have:
                if (s, year, doy) not in idx:
                    continue
                rows.append(("neg_certain_day", year, doy, d.strftime("%Y-%m-%d"),
                             idx[(s, year, doy)].split(".")[0], s, -1, np.nan, "pnsn_quiet",
                             np.nan, False, s in HSTA,
                             "quiet station-day; EXCLUDE +-120 s around every Bostock event",
                             np.nan))
                nneg += 1
        # ---- unlabelled: ETS station-days ----
        nunl = 0
        for d in days:
            if d not in ets_days:
                continue
            doy = d.dayofyear
            for s in have:
                if (s, year, doy) not in idx:
                    continue
                rows.append(("unlabelled_ets_day", year, doy, d.strftime("%Y-%m-%d"),
                             idx[(s, year, doy)].split(".")[0], s, -1, np.nan, "ets_day", np.nan,
                             False, s in HSTA,
                             "non-event windows on this day are UNLABELLED, never negatives", np.nan))
                nunl += 1
        inv.append(dict(year=year, n_days=len(days), n_stations=len(have), n_events=int(len(ev)),
                        pnsn_covered=bool(covered),
                        n_quiet_days=int(len(qd & set(days))), n_ets_days=int(len(ets_days)),
                        n_pos_S=npos, n_neg_certain_day=nneg, n_unlabelled_ets_day=nunl))
        print(f"[{year}] {len(have)} stations, {len(ev)} events -> {npos} pos_S, {nneg} quiet "
              f"station-days, {nunl} ETS station-days", flush=True)

    # ---------- hard negatives: catalogued earthquakes / blasts ----------
    EQ = pd.read_csv(f"{EXT}/results/trackC/eq_catalog.csv")
    EQ["ot"] = pd.to_datetime(EQ.ot, errors="coerce")
    EQ = EQ.dropna(subset=["ot"])
    neq = 0
    for year, (w0, w1) in WINDOWS.items():
        e = EQ[(EQ.ot >= w0) & (EQ.ot < w1)]
        have = sorted({s for (s, y, d) in idx if y == year})
        for r in e.itertuples():
            doy = r.ot.dayofyear
            for s in have:
                if (s, year, doy) not in idx:
                    continue
                rows.append(("hard_neg_eq", year, doy, r.ot.strftime("%Y-%m-%d"),
                             idx[(s, year, doy)].split(".")[0], s, -1,
                             float(r.ot.value) / 1e9, "catalog_eq", np.nan, False, s in HSTA,
                             "eq_catalog.csv (PNSN+USGS)",
                             float(getattr(r, "mag", np.nan)) if hasattr(r, "mag") else np.nan))
                neq += 1

    # ---------- hard negatives from the frozen detector on QUIET family-days ----------
    nhf = ntx = 0
    for year in (2010, 2011, 2012, 2013):
        if year in DET_SRC:
            tag, base = DET_SRC[year]
            D = pd.read_parquet(f"{base}/results/e5r_detections_{tag}.parquet")
            H = pd.read_parquet(f"{base}/results/a1_hf_ratio_{tag}.parquet")[["doy", "fam", "t", "hfr_PGC"]]
            D = D[D.L == 12.0].merge(H, on=["doy", "fam", "t"], how="left")
            D["hf_pass"] = ~(np.isfinite(D.hfr_PGC) & (D.hfr_PGC > HF_CUT))
            qdoy = set(range(228, 243)) if year == 2012 else {223, 224, 226, 252, 254, 255, 256}
            D = D[D.doy.isin(qdoy)]
            D["date"] = [(pd.Timestamp(f"{year}-01-01") + pd.Timedelta(days=int(d) - 1)).strftime("%Y-%m-%d")
                         for d in D.doy]
        else:
            months = {2011: ["2011-08", "2011-09"], 2013: ["2013-08", "2013-09", "2013-10"]}[year]
            D = pd.concat([pd.read_parquet(f"{RES}/trackO/detections_{m}.parquet") for m in months],
                          ignore_index=True)
            D = D[D.quiet_flag]
        for r in D.itertuples():
            b = int(r.fam)
            if not r.hf_pass:
                rows.append(("hard_neg_hfveto", year, int(r.doy), str(r.date), "", "", b,
                             float(r.t), "hf_veto_quietday", np.nan, b in HID, False,
                             "frozen detector, HF ratio > cut on a quiet family-day", float(r.net_mad)))
                nhf += 1
            elif r.net_mad >= 8.25:
                rows.append(("hard_neg_text", year, int(r.doy), str(r.date), "", "", b,
                             float(r.t), "T_ext_quietday", np.nan, b in HID, False,
                             "frozen detector, T_ext detection on a quiet family-day", float(r.net_mad)))
                ntx += 1

    cols = ["kind", "year", "doy", "date", "net", "sta", "fam", "t", "provenance", "sigma_s",
            "hidden_family", "held_out_station", "source", "aux"]
    L = pd.DataFrame(rows, columns=cols)
    L["time"] = pd.to_datetime(L.t, unit="s", errors="coerce")
    L.to_parquet(a.out, index=False)
    print(f"\n[labels] {len(L)} rows -> {a.out}")
    print(L.groupby("kind").size().to_string())
    S = L[L.kind == "pos_S"]
    print("\npositives by provenance:\n" + S.groupby("provenance").size().to_string())
    print("\npositives by station x provenance:\n" +
          S.pivot_table(index="sta", columns="provenance", values="t", aggfunc="count",
                        fill_value=0).to_string())
    print("\npositives by year:\n" + S.groupby("year").size().to_string())
    print(f"\nhidden-family rows: {int(L.hidden_family.sum())};  held-out-station rows: "
          f"{int(L.held_out_station.sum())}")
    summ = dict(n_rows=int(len(L)), by_kind=L.groupby("kind").size().to_dict(),
                pos_by_provenance=S.groupby("provenance").size().to_dict(),
                pos_by_station=S.groupby("sta").size().to_dict(),
                pos_by_year=S.groupby("year").size().to_dict(),
                pos_picked_by_station=S[S.provenance == "picked"].groupby("sta").size().to_dict(),
                n_hidden_family_rows=int(L.hidden_family.sum()),
                n_held_out_station_rows=int(L.held_out_station.sum()),
                n_hard_neg_eq=neq, n_hard_neg_hfveto=nhf, n_hard_neg_text=ntx,
                per_year=inv, sigma_s=SIGMA,
                rules=dict(S="DB time + frozen 2012 family offset + (a5 stack S pick | Crosson ttS) "
                             "+ a5 station correction",
                           P="only where a5_picks.t_P is finite",
                           certain_negative="PNSN-quiet station-day, +-120 s around every Bostock "
                                            "event excluded",
                           unlabelled="ETS-day non-event windows are UNLABELLED, never negatives"))
    json.dump(summ, open(f"{RES}/labels_v0_summary.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
