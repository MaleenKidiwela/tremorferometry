#!/usr/bin/env python
"""P0-O item 1 -- REFEREE FREEZE and the hashed holdouts.

Writes results/referee_protocol.json (which year may be looked at, how often, with which frozen
time base) and results/holdouts.json (30 hidden families in three nearest-neighbour strata + 3
held-out stations).  Both files are WRITE-ONCE: re-running refuses to overwrite unless --force.

The per-family time offsets for EVERY year are the frozen 2012 offsets
(../lfe_catalog_repair/results/e5r_family_offsets_2012.csv).  They are never refit on 2010, 2011
or 2013 -- refitting the time base on the referee year is the classic way to launder a held-out
year into a development year.
"""
from __future__ import annotations
import os, sys, json, hashlib, argparse
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
PICK = os.path.dirname(HERE)
RES = f"{PICK}/results"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
IN = f"{REP}/results"
BOST = f"{IN}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"
R = 111.195
SEED = 20260917


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    for p in (f"{RES}/referee_protocol.json", f"{RES}/holdouts.json"):
        if os.path.exists(p) and not a.force:
            raise SystemExit(f"[freeze] {p} exists -- refusing to overwrite (use --force deliberately)")

    FAM = pd.read_csv(f"{IN}/e5r_families_2012.csv")          # bkey 0..129, Bostock hypocentres
    OFF = pd.read_csv(f"{IN}/e5r_family_offsets_2012.csv")
    STA = pd.read_csv(f"{IN}/e5r_stations.csv")
    FZ = f"{EXT}/results/frozen_operating_point.json"
    TPL = f"{IN}/e5r_templates_2012.npz"

    # ---- Bostock event counts per family per referee year (exact hypocentre join, 130/130) ----
    B = pd.read_csv(BOST)
    B["time"] = pd.to_datetime(dict(year=B.year, month=B.month, day=B.day, hour=B.hour,
                                    minute=B["min"], second=B.sec.astype(int)), errors="coerce")
    key = {(round(r.lat, 4), round(r.lon, 4), round(r.dep, 4)): int(r.bkey) for r in FAM.itertuples()}
    B["bkey"] = [key.get((round(la, 4), round(lo, 4), round(dp, 4)), -1)
                 for la, lo, dp in zip(B.lat, B.lon, B.dep)]
    n_unmapped = int((B.bkey < 0).sum())
    Bm = B[B.bkey >= 0]
    yr_counts = {}
    for y in (2008, 2009, 2010, 2011, 2012, 2013):
        g = Bm[Bm.time.dt.year == y].groupby("bkey").size()
        yr_counts[y] = {int(k): int(v) for k, v in g.items()}

    # ---- nearest-neighbour family distance (epicentral, and hypocentral for the record) ----
    cs = np.cos(np.radians(48.6))
    X = np.c_[FAM.lat * R, FAM.lon * R * cs]
    Z = FAM.dep.to_numpy()
    n = len(FAM)
    DH = np.sqrt(((X[:, None, :] - X[None, :, :]) ** 2).sum(-1))
    D3 = np.sqrt(DH ** 2 + (Z[:, None] - Z[None, :]) ** 2)
    np.fill_diagonal(DH, np.inf); np.fill_diagonal(D3, np.inf)
    nn_h = DH.min(1); nn_3 = D3.min(1)
    nn_id = FAM.bkey.to_numpy()[DH.argmin(1)]

    # ---- eligibility: a hidden family must have something to rediscover ----
    # >= 50 Bostock events in the 2012 referee window (the window every stage is measured on) AND
    # >= 20 in 2011 (so the one-shot 2011 gate can see it).  Nothing about our detector enters here.
    n2012 = FAM.set_index("bkey").n_all.to_dict()
    n2011 = yr_counts[2011]
    elig = np.array([(n2012.get(int(b), 0) >= 50) and (n2011.get(int(b), 0) >= 20)
                     for b in FAM.bkey])
    strata = {"isolated_ge5km": (nn_h >= 5.0), "mid_2.5_5km": (nn_h >= 2.5) & (nn_h < 5.0),
              "close_lt2.5km": nn_h < 2.5}
    rng = np.random.default_rng(SEED)
    hidden, stratum_report = [], {}
    for name, m in strata.items():
        pool = FAM.bkey.to_numpy()[m & elig]
        stratum_report[name] = dict(n_eligible=int(len(pool)), n_requested=10)
        pick = rng.choice(np.sort(pool), size=min(10, len(pool)), replace=False)
        stratum_report[name]["n_picked"] = int(len(pick))
        for b in sorted(int(x) for x in pick):
            i = int(np.where(FAM.bkey.to_numpy() == b)[0][0])
            hidden.append(dict(bkey=b, stratum=name,
                               lat=float(FAM.lat.iloc[i]), lon=float(FAM.lon.iloc[i]),
                               dep=float(FAM.dep.iloc[i]),
                               nn_km_epicentral=float(nn_h[i]), nn_km_hypocentral=float(nn_3[i]),
                               nn_family=int(nn_id[i]),
                               n_bostock_2012=int(n2012.get(b, 0)),
                               n_bostock_2011=int(n2011.get(b, 0)),
                               n_bostock_2010=int(yr_counts[2010].get(b, 0)),
                               n_bostock_2013=int(yr_counts[2013].get(b, 0))))

    hold_sta = ["B926", "NLLB", "LRIV"]
    H = dict(
        created_utc=pd.Timestamp.utcnow().isoformat(),
        seed=SEED,
        rule=("30 hidden Bostock families: 10 per nearest-neighbour stratum (epicentral NN distance "
              ">=5 km / 2.5-5 km / <2.5 km), drawn with numpy default_rng(seed) from the eligible "
              "pool, sorted by bkey so the draw is reproducible; plus 3 held-out stations."),
        eligibility="n_bostock_2012 >= 50 AND n_bostock_2011 >= 20 (nothing from our detector)",
        distance="epicentral km at cos(48.6 deg); hypocentral recorded alongside",
        family_universe="the 130 Bostock families of e5r_families_2012.csv (bkey 0..129)",
        strata=stratum_report,
        hidden_families=hidden,
        held_out_stations=hold_sta,
        held_out_station_rule=("excluded from ALL training and from every label; evaluated as unseen "
                               "stations only (Merlin P4a)"),
        nn_distribution=dict(median_km=float(np.median(nn_h)), p10=float(np.percentile(nn_h, 10)),
                             p90=float(np.percentile(nn_h, 90)),
                             frac_lt_3km=float((nn_h < 3).mean())),
        usage=("hidden_families and held_out_stations are removed from the positive class of every "
               "P1-P3 training set and from every label table (flagged, rows kept); they are used "
               "ONLY to score rediscovery (P3) and unseen-station transfer (P4a)."))
    H["sha256_of_hidden_family_list"] = hashlib.sha256(
        json.dumps([h["bkey"] for h in hidden], sort_keys=True).encode()).hexdigest()
    json.dump(H, open(f"{RES}/holdouts.json", "w"), indent=1)

    P = dict(
        created_utc=pd.Timestamp.utcnow().isoformat(),
        authority=("lfe_picker/results/merlin_review_plan_2026-09-17.md, P0-O item 1. This file is "
                   "the protocol: any use of a year beyond its budget below invalidates the phase."),
        years=dict(
            y2013=dict(role="FINAL REFEREE", bostock_window="2013-09-13..2013-10-13",
                       data_window="2013-08-30..2013-10-28",
                       budget="ONE look, at the END of the picker phase",
                       spent_in_P0O=1,
                       spent_note=("P0-O item 2 spends the single allowed look to establish the "
                                   "P1a template baseline. It is RECORDED and must not be used to "
                                   "select, tune or compare anything afterwards."),
                       n_bostock_events_in_window=int(((Bm.time >= "2013-09-13") &
                                                       (Bm.time < "2013-10-14")).sum())),
            y2011=dict(role="STAGE GATE", bostock_window="2011-08-19..2011-09-09 (+-)",
                       data_window="2011-08-05..2011-09-25",
                       scoring_window="all 2011 Bostock events on disk",
                       budget="exactly ONE evaluation per stage gate (P1a, P1b, P2, P3, P4)",
                       spent_in_P0O=1,
                       n_bostock_events_2011=int((Bm.time.dt.year == 2011).sum())),
            y2010=dict(role="DEVELOPMENT", data_window="2010-08-10..2010-09-13", budget="unlimited"),
            y2012=dict(role="DEVELOPMENT", data_window="2012-08-15..2012-10-05",
                       budget="unlimited; the oracle (P0-O item 3) and all fitting live here"),
            y2008=dict(role="TRAINING BREADTH", data_window="2008-05-01..2008-06-10",
                       budget="training only, never evaluation", status="queued for download"),
            y2009=dict(role="TRAINING BREADTH", data_window="2009-05-01..2009-06-15",
                       budget="training only, never evaluation", status="queued for download")),
        time_base=dict(
            rule=("per-family time offsets for ALL years = the FROZEN 2012 offsets; NEVER refit on "
                  "2010, 2011 or 2013"),
            file=f"{IN}/e5r_family_offsets_2012.csv",
            n_families=int(len(OFF)), n_accepted=int(OFF.accept.astype(bool).sum()),
            t_origin="t_origin = Bostock DB time + offset_s(bkey)"),
        family_map=dict(rule="exact (lat, lon, dep) match of the Bostock row to e5r_families_2012.csv",
                        n_families=int(len(FAM)),
                        n_bostock_rows_unmapped_2003_2013=n_unmapped),
        matching=dict(rule="own family, +-2 s, one-to-one greedy by time (Track A step 0(b))"),
        detector=dict(configuration="the FROZEN operating point, unchanged",
                      variant="mergeB, 12-s 2012 templates, SNR weights, HF veto, T_core 8.75 / T_ext 8.25",
                      station_set_S15="PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926"),
        inputs={
            "bostock_catalog": dict(path=BOST, sha256=sha256(BOST)),
            "family_offsets_2012": dict(path=f"{IN}/e5r_family_offsets_2012.csv",
                                        sha256=sha256(f"{IN}/e5r_family_offsets_2012.csv")),
            "templates_2012": dict(path=TPL, sha256=sha256(TPL)),
            "frozen_operating_point": dict(path=FZ, sha256=sha256(FZ)),
            "stations": dict(path=f"{IN}/e5r_stations.csv", sha256=sha256(f"{IN}/e5r_stations.csv"),
                             n=int(len(STA)), codes=",".join(sorted(STA.sta.astype(str)))),
            "families_2012": dict(path=f"{IN}/e5r_families_2012.csv",
                                  sha256=sha256(f"{IN}/e5r_families_2012.csv")),
            "weights": dict(path=f"{IN}/e5r_weights_snr_2012.csv",
                            sha256=sha256(f"{IN}/e5r_weights_snr_2012.csv")),
            "moveout": dict(path=f"{IN}/e4_moveout_e5r_mv_build.csv",
                            sha256=sha256(f"{IN}/e4_moveout_e5r_mv_build.csv")),
        },
        holdouts=dict(path=f"{RES}/holdouts.json",
                      sha256_of_hidden_family_list=H["sha256_of_hidden_family_list"],
                      hidden_families=[h["bkey"] for h in hidden], held_out_stations=hold_sta),
    )
    json.dump(P, open(f"{RES}/referee_protocol.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in stratum_report.items()}, indent=1))
    print("hidden:", [h["bkey"] for h in hidden])
    print("nn median %.2f km, p10 %.2f, p90 %.2f, frac<3km %.2f"
          % (np.median(nn_h), np.percentile(nn_h, 10), np.percentile(nn_h, 90), (nn_h < 3).mean()))
    for k, v in P["inputs"].items():
        print(f"  {k}: {v['sha256'][:16]}  {os.path.basename(v['path'])}")


if __name__ == "__main__":
    main()
