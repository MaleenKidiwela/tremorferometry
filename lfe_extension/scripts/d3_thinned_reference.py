#!/usr/bin/env python
"""E6 step 4(iii) -- the COMPLETENESS-MATCHED reference catalog (Merlin E5r-audit F7): Bostock's 2012 catalog thinned per FAMILY by
Bernoulli sampling at OUR measured own-family recall of that family (results/e6_score_<tag>_recall_family.csv at the operating
threshold), plus a uniform-in-time FLOOR per templated family at our measured quiet-day rate (e5_budget_<tag>.csv, doy 228-242, mean
detections / template / day x 52 days, Poisson), 20 seeds; the unchanged e5_thomas_harness.py (--mc none) is run on every seed and on
OUR catalog at the same operating point.  Summary = mean, sd and 2.5-97.5 % range over seeds of the main-front first-burst speed,
the re-rupture interval quantiles and family fraction, and the secondary-front count (the count is reported but F6 says it is not
density-monotonic; step 4(iv) replaces it).  Nothing in results/e5* is touched; per-seed harness files go to results/e6_thinned/.
Writes results/e6_thinned_reference_<tag>_T<T>.csv (+ the per-seed catalogs and harness summaries under results/e6_thinned/)."""
from __future__ import annotations
import argparse, glob, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd
RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
TD = f"{EXT}/results/trackD"
OUT = f"{TD}/thinned"
BOS = f"{RES}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"
HARNESS = f"{EXT}/scripts/d3_thomas_harness.py"
KEYS = ["n_events", "front_firstburst_kmday", "front_firstburst_ci_lo", "front_firstburst_ci_hi", "front_firstburst_rho", "n_fam_firstburst_inwin",
        "front_lead_kmday", "n_secondary", "n_secondary_null_circ", "n_bursts", "n_reactivations", "n_fam_reactivated", "frac_fam_reactivated",
        "reactivation_interval_d_p10", "reactivation_interval_d_p50", "reactivation_interval_d_p90", "reactivation_frac_le5d"]


def run_harness(cat, tag, env):
    log = f"{OUT}/log_{tag}.txt"
    with open(log, "w") as fh:
        r = subprocess.run([sys.executable, HARNESS, "--catalog", cat, "--tag", tag, "--mc", "none", "--outdir", OUT,
                            "--t0", "2012-08-15", "--t1", "2012-10-16"], stdout=fh, stderr=subprocess.STDOUT, env=env)
    return tag, r.returncode


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="a1_mergeB_L12"); ap.add_argument("--T", type=float, required=True)
    ap.add_argument("--epoch", default="E1"); ap.add_argument("--tier", default="T_core")
    ap.add_argument("--seeds", type=int, default=20); ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--ours-catalog", default="", help="our catalog parquet at the operating point (e6_catalog_2012_sVI_<tag>_T<T>.parquet)")
    ap.add_argument("--no-floor", action="store_true", help="also build a variant WITHOUT the uniform floor (labelled)")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    env = dict(os.environ, PYTHONPATH="/opt/conda/lib/python3.13/site-packages", OMP_NUM_THREADS="2")
    rf = pd.read_csv(f"{TD}/completeness_family_epoch.csv")
    rf = rf[(rf.epoch_id == a.epoch) & (rf.tier == a.tier)]
    rec = dict(zip(rf.family_id.astype(int), rf.recall.fillna(0.0)))
    # quiet-day rate per family at the EXACT frozen threshold, with the frozen HF veto applied (the
    # budget grid does not carry 8.75, so this is read off the detections parquet itself)
    det = pd.read_parquet(f"{EXT}/results/e5r_detections_{a.tag}.parquet"); det = det[det.L == 12.0]
    H = pd.read_parquet(f"{EXT}/results/a1_hf_ratio_{a.tag}.parquet"); H = H[H.L == 12.0][["doy", "fam", "t", "hfr_PGC"]]
    det = det.merge(H, on=["doy", "fam", "t"], how="left")
    hv = det.hfr_PGC.to_numpy(); det = det[~(np.isfinite(hv) & (hv > 3.6384))]
    bud = pd.read_csv(f"{EXT}/results/e5_budget_{a.tag}.csv"); bud = bud[bud.L == 12.0]
    bud = bud[bud["T"] == sorted(bud["T"].unique())[0]]
    ndq = bud[bud.doy.isin(range(228, 243))].drop_duplicates(["fam", "doy"]).groupby("fam").size()
    q = det[(det.net_mad >= a.T) & det.doy.isin(range(228, 243))]
    rate = (q.groupby("fam").size().reindex(ndq.index).fillna(0) / ndq).to_dict()
    F = pd.read_csv(f"{RES}/e5r_families_2012.csv")
    c = pd.read_csv(BOS)
    t = pd.to_datetime(dict(year=c.year, month=c.month, day=c.day, hour=c.hour, minute=c["min"], second=c.sec.astype(int))) + pd.to_timedelta(c.sec % 1, unit="s")
    c = c[(t >= "2012-08-15") & (t < "2012-10-06")].copy().reset_index(drop=True)
    # bkey of every Bostock event = its (lat, lon, dep) family in e5r_families_2012.csv
    key = F.set_index([F.lat.round(3), F.lon.round(3), F.dep.round(2)]).bkey
    c["bkey"] = pd.Series(list(zip(c.lat.round(3), c.lon.round(3), c.dep.round(2)))).map(key.to_dict()).to_numpy()
    print(f"[ref] {len(c)} Bostock events in the window, {c.bkey.notna().sum()} keyed; per-family own recall p10/50/90 "
          f"{np.percentile(list(rec.values()), [10, 50, 90]).round(3)}; floor rate p50/90 {np.percentile(list(rate.values()), [50, 90]).round(2)} / template / day "
          f"({len(rate)} templated families) at T {a.T:g}", flush=True)
    t0 = pd.Timestamp("2012-08-15").value / 1e9; span = 52 * 86400.0
    jobs = []
    for seed in range(a.seeds):
        rng = np.random.default_rng(1000 + seed)
        p = c.bkey.map(rec).fillna(0.0).to_numpy()
        keep = rng.random(len(c)) < p
        thin = c[keep].copy()
        # uniform-in-time floor per templated family (Poisson count at the measured quiet-day rate over the 52 days)
        fl = []
        for fk, r_ in rate.items():
            nf = rng.poisson(r_ * 52.0)
            if nf == 0:
                continue
            tt = t0 + rng.uniform(0, span, nf); ts = pd.to_datetime(tt, unit="s")
            fr = F[F.bkey == int(fk)].iloc[0]
            fl.append(pd.DataFrame(dict(year=ts.year, month=ts.month, day=ts.day, hour=ts.hour, min=ts.minute, second=ts.second + ts.microsecond / 1e6,
                                        lat=fr.lat, lon=fr.lon, dep=fr.dep, mag=1.5, bkey=int(fk))).rename(columns={"second": "sec"}))
        floor = pd.concat(fl, ignore_index=True) if fl else thin.iloc[0:0]
        for lab, cat in (("floor", pd.concat([thin, floor], ignore_index=True)), ("nofloor", thin)):
            if lab == "nofloor" and not a.no_floor:
                continue
            fn = f"{OUT}/e6_ref_{a.tag}_T{a.T:g}_{lab}_s{seed}.csv"
            cat.drop(columns=["bkey"]).to_csv(fn, index=False)
            jobs.append((fn, f"e6ref_{a.tag}_T{a.T:g}_{lab}_s{seed}"))
        if seed == 0:
            print(f"[seed 0] thinned {len(thin)} events + floor {len(floor)} = {len(thin) + len(floor)} (Bostock 28,149; floor share {len(floor)/(len(thin)+len(floor)):.2f})", flush=True)
    if a.ours_catalog:
        jobs.append((a.ours_catalog, f"e6ours_{a.tag}_T{a.T:g}"))
    print(f"[harness] {len(jobs)} runs, {a.workers} parallel", flush=True)
    with ThreadPoolExecutor(a.workers) as ex:
        for tag, rc in ex.map(lambda j: run_harness(*j, env), jobs):
            if rc != 0:
                print(f"  [FAILED] {tag} (see {OUT}/log_{tag}.txt)", flush=True)
    rows = []
    for fn in sorted(glob.glob(f"{OUT}/e5_thomas_summary_e6ref_{a.tag}_T{a.T:g}_*.csv")) + sorted(glob.glob(f"{OUT}/e5_thomas_summary_e6ours_{a.tag}_T{a.T:g}.csv")):
        s = pd.read_csv(fn); s = s[s.mc.astype(str) == "none"]
        if len(s):
            r = {k: s.iloc[0].get(k, np.nan) for k in KEYS}; r["run"] = os.path.basename(fn)[len("e5_thomas_summary_"):-4]
            r["kind"] = "ours" if "e6ours" in fn else ("ref_nofloor" if "_nofloor_" in fn else "ref_floor"); rows.append(r)
    R = pd.DataFrame(rows)
    b = pd.read_csv(f"{RES}/e5_thomas_summary_bostock2012.csv"); b = b[b.mc.astype(str) == "none"].iloc[0]
    summ = []
    for kind, g in R.groupby("kind"):
        for k in KEYS:
            v = pd.to_numeric(g[k], errors="coerce").dropna()
            summ.append(dict(kind=kind, metric=k, n_runs=int(len(v)), mean=float(v.mean()) if len(v) else np.nan, sd=float(v.std()) if len(v) > 1 else np.nan,
                             p2_5=float(v.quantile(.025)) if len(v) else np.nan, p97_5=float(v.quantile(.975)) if len(v) else np.nan, bostock_full=float(b.get(k, np.nan)) if k in b else np.nan))
    Sm = pd.DataFrame(summ); Sm.to_csv(f"{TD}/thinned_reference_{a.tier}_T{a.T:g}.csv", index=False)
    R.to_csv(f"{TD}/thinned_reference_{a.tier}_T{a.T:g}_runs.csv", index=False)
    pd.set_option("display.width", 250)
    print(Sm.round(3).to_string(index=False))
    print(f"[write] {TD}/thinned_reference_{a.tier}_T{a.T:g}.csv (+ _runs.csv; per-seed files in {OUT}/)")


if __name__ == "__main__":
    main()
