#!/usr/bin/env python
"""bw_p7_catalogue -- catalogue v1.4-backward (rules: results/backward/p7/P7.md, "Catalogue v1.4-backward").
Every data and scrambled-control candidate of p5f_2003..2013 (frozen) and p5e_2003_kelb, p5e_2004..2006 (era), with the
forward per-candidate catalogue columns (refine/conf/P_by_candidate.parquet) plus era, net_set, tpl_epoch, detector, P, model,
burst_day, burst_ratio, n_z2, multi_site_z2, counted_once, primary, ctl_heavy_set, stations_day, row.
Writes data/catalogue/svi_catalogue_v1_4_backward.parquet (+ _meta.json) as NEW files; refuses to overwrite.
"""
import os, sys, glob, json
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

ROOT = "/home/jovyan/tremorferometry"
BW = f"{ROOT}/svi_science/results/backward"
P6 = f"{BW}/p6"
OUTF = f"{ROOT}/data/catalogue/svi_catalogue_v1_4_backward.parquet"
RUNS = [f"p5f_{y}" for y in range(2003, 2014)] + ["p5e_2003_kelb", "p5e_2004", "p5e_2005", "p5e_2006"]
FEATC = ["series", "fam", "date", "t", "k", "net_mad", "n_sta_day", "has_LRIV", "has_SNB", "conf_z", "conf_ze", "keep_xfam8",
         "hfr", "n_sta_clbad", "sset"]
FWD_P = ["P_same_WAVE", "P_any_WAVE", "P_same_FULL", "P_any_FULL", "P_same_NO_FAM", "P_any_NO_FAM", "P_same_NO_DENSITY",
         "P_any_NO_DENSITY", "P_same_MF_ONLY", "P_any_MF_ONLY", "P_same_NO_PNET", "P_same_WITH_PNET", "P_any_NO_PNET",
         "P_any_WITH_PNET", "P_same", "P_any", "P_same_WAVE_foldsd", "P_any_WAVE_foldsd"]
FWD_COLS = ["run", "date", "series", "k", "fam", "t", "net_mad", "sset_run", "n_sta_day", "has_LRIV", "has_SNB", "quiet",
            "n_tremor_day", "conf_z", "conf_ze", "keep_xfam8", "hfr", "n_sta_clbad"] + FWD_P + ["month", "det_row", "cat_core",
                                                                                                "cat_ext", "in_v11"]
ERA_BOUND = [("2003-01-01", "B0"), ("2005-07-01", "B1"), ("2005-11-02", "B2"), ("2007-11-13", "B3"), ("2009-08-27", "B4")]
ERA_DET = ("2003-02-01", "2006-06-30")


def era_of(d):
    e = "B0"
    for b, n in ERA_BOUND:
        if d >= b: e = n
    return e


def tpl_epoch(d, det):
    import bw_gate_common as G
    if det == "era":
        return "2012+era2004(KELB)"
    s = ["2012"]
    for sta, eps in G.INSTR.items():
        for a, b, lab in eps:
            if a <= d < b and "2012 build epoch" not in lab:
                s.append(f"{sta} {lab}")
    return "; ".join(s)


def one_day(args):
    run, ds = args
    p = f"{BW}/mf/feat/{run}/{ds}.parquet"
    if not os.path.exists(p) or os.path.getsize(p) == 0:
        return None
    x = pd.read_parquet(p)
    if not len(x):
        return None
    x = x[[c for c in FEATC if c in x.columns]].copy(); x["row"] = np.arange(len(x))
    js = f"{BW}/mf/feat/{run}/{ds}.json"
    x["stations_day"] = ",".join(json.load(open(js)).get("stations", [])) if os.path.exists(js) else ""
    return x


def main():
    if os.path.exists(OUTF) and "--force" not in sys.argv:
        raise SystemExit(f"[cat] {OUTF} exists; refusing to overwrite")
    Q = pd.read_csv(f"{BW}/p4/envflag/tables/quiet_flag_env.csv", dtype={"date": str})
    Q = Q.set_index("date")
    cal = json.load(open(f"{P6}/tables/calib_decision.json"))
    heavy = {k.split(":")[1] if k.startswith("frozen") else "ERA": bool(v["clearly_above"]) for k, v in cal["sets"].items()}
    parts = []
    for run in RUNS:
        ds = sorted(os.path.basename(p)[:10] for p in glob.glob(f"{BW}/mf/feat/{run}/*.json"))
        with ProcessPoolExecutor(8) as ex:
            X = pd.concat([r for r in ex.map(one_day, [(run, d) for d in ds], chunksize=8) if r is not None], ignore_index=True)
        X["run"] = run
        P = pd.read_parquet(f"{P6}/P/{run}_P.parquet")
        X = X.merge(P[["date", "row", "P", "P_same", "P_any", "P_same_foldsd", "P_any_foldsd", "model", "net_set"]], on=["date", "row"],
                    how="left", validate="one_to_one")
        S = pd.read_parquet(f"{BW}/p4b/sidecar/{run}_flags.parquet", columns=["date", "row", "series", "burst_day", "burst_ratio", "n_z2",
                                                                                 "multi_site_z2"])
        X = X.merge(S.drop(columns=["series"]), on=["date", "row"], how="left", validate="one_to_one")
        if X.P.isna().any() or X.burst_day.isna().any():
            raise SystemExit(f"[cat] {run}: {int(X.P.isna().sum())} rows without P, {int(X.burst_day.isna().sum())} without sidecar flags")
        parts.append(X); print(f"[cat] {run}: {len(X)} rows", flush=True)
    A = pd.concat(parts, ignore_index=True)
    A["detector"] = np.where(A.run.str.startswith("p5e_"), "era", "frozen")
    A["sset_run"] = A.sset.astype(str); A = A.drop(columns=["sset"])
    A["era"] = A.date.map({d: era_of(d) for d in A.date.unique()})
    te = {(d, det): tpl_epoch(d, det) for d, det in A[["date", "detector"]].drop_duplicates().itertuples(index=False)}
    A["tpl_epoch"] = [te[(d, det)] for d, det in zip(A.date, A.detector)]
    in_era = (A.date >= ERA_DET[0]) & (A.date <= ERA_DET[1])
    A["primary"] = np.where(A.detector == "era", in_era, ~in_era)
    A["keep_xfam8"] = A.keep_xfam8.astype(float)
    A["counted_once"] = (A.series == "data") & (A.keep_xfam8 > 0) & A.primary
    A["ctl_heavy_set"] = A.net_set.map(heavy).fillna(False).astype(bool)
    pq_ = Q.pnsn_quiet.reindex(A.date.unique()); eq = Q.quiet_env3_CN4.reindex(A.date.unique())
    post = {d: d >= "2010-01-10" for d in A.date.unique()}
    qmap = {d: (bool(pq_.get(d)) if post[d] and pd.notna(pq_.get(d)) else bool(eq.get(d)) if pd.notna(eq.get(d)) else False) for d in A.date.unique()}
    A["quiet"] = A.date.map(qmap).astype(bool)
    A["quiet_src"] = np.where(A.date >= "2010-01-10", "pnsn", "env3_CN4")
    A["quiet_env3_CN4"] = A.date.map(eq.fillna(False).astype(bool).to_dict()).astype(bool)
    ntd = Q.n_tremor_day.reindex(A.date.unique()).to_dict()
    A["n_tremor_day"] = [float(ntd.get(d)) if d >= "2010-01-06" else np.nan for d in A.date]
    A["P_same_WAVE"] = A.P_same; A["P_any_WAVE"] = A.P_any
    A["P_same_WAVE_foldsd"] = A.P_same_foldsd; A["P_any_WAVE_foldsd"] = A.P_any_foldsd
    for c in FWD_P:
        if c not in A: A[c] = np.nan
    A["month"] = A.date.str[:7]; A["det_row"] = np.nan; A["cat_core"] = None; A["cat_ext"] = None; A["in_v11"] = False
    A["has_LRIV"] = A.has_LRIV.astype(float); A["has_SNB"] = A.has_SNB.astype(float)
    extra = ["era", "net_set", "tpl_epoch", "detector", "P", "model", "burst_day", "burst_ratio", "n_z2", "multi_site_z2", "counted_once",
             "primary", "ctl_heavy_set", "quiet_src", "quiet_env3_CN4", "stations_day", "row"]
    A = A[FWD_COLS + extra].sort_values(["date", "detector", "t"]).reset_index(drop=True)
    A.to_parquet(OUTF, index=False)
    meta = dict(created_utc=pd.Timestamp.now("UTC").isoformat(), rows=len(A), by_series=A.series.value_counts().to_dict(),
                by_detector=A.detector.value_counts().to_dict(), rules="svi_science/results/backward/p7/P7.md",
                models="svi_science/results/backward/p6/models/MODELS_FROZEN_BW.json",
                counted_once="data & keep_xfam8 & primary detector of the date (era 2003-02-01..2006-06-30, frozen otherwise)",
                headline_P075_counted_once=int((A.counted_once & (A.P >= 0.75)).sum()),
                headline_P075_ctl_primary=int(((A.series == "c1") & A.primary & (A.keep_xfam8 > 0) & (A.P >= 0.75)).sum()),
                note="P = P_same of the row's model (B0/B1/B2/era fitted here, B3/B4 frozen forward WAVE). Forward model variants are NaN. "
                     "quiet = PNSN rule from 2010-01-10, env-quiet (quiet_env3_CN4) before. No v1.1 join before 2014.")
    json.dump(meta, open(OUTF.replace(".parquet", "_meta.json"), "w"), indent=1, default=str)
    print(json.dumps(meta, indent=1, default=str))


if __name__ == "__main__":
    main()
