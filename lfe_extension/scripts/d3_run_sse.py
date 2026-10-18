#!/usr/bin/env python
"""Track D, D3 -- SCIENCE QC: the Thomas-2025 harness per SSE (Merlin Step 8).

Runs the E5 harness (copied verbatim as `scripts/d3_thomas_harness.py`; the only edit is that the
window-dependent constants T_START and DIAG_DAY became command-line parameters) on

  * every SSE window of `results/trackD/sse_windows.csv`, at BOTH tiers (cat_core, cat_ext), and
  * the 2012 REFEREE as the control, with the same code: OUR frozen-operating-point 2012 catalog
    (a1_mergeB_L12 + HF veto + the frozen cross-family assignment = the production definition) and
    Bostock's own 2012 catalog.

Pre-registered bars (reported met / not met, NEVER tuned), evaluated on the 2012 referee only, because
they are comparisons against Bostock's catalog and Bostock stops in 2013:
  (1) main-front first-burst speed within 20 % of Bostock's 2012 value;
  (2) re-rupture (re-activation) interval median within 5 % of Bostock's 2012 value.
Everything else -- secondary-front counts against the circular and day-shuffle nulls, re-activation
fractions, per-SSE front speeds -- is REPORTED.

Family positions are Bostock's hypocentres for every family (not our relocations), so the along-strike
coordinate is one consistent frame; the 10 accepted relocations move families by a median 3.6 km, well
inside the front-tracking bin width."""
from __future__ import annotations
import argparse, json, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
import numpy as np, pandas as pd

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
TD = f"{EXT}/results/trackD"
TH = f"{TD}/thomas"
HARNESS = f"{EXT}/scripts/d3_thomas_harness.py"
BOSTOCK = f"{IN}/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv"
HF_CUT = 3.6384
TIERS = (("core", "cat_core"), ("ext", "cat_ext"), ("core_hp", "cat_core_hp"))
KEYS = ["n_events", "n_families", "win_how", "win_start", "win_end", "win_days", "n_bursts",
        "n_fam_firstburst", "n_fam_firstburst_inwin", "front_firstburst_kmday",
        "front_firstburst_ci_lo", "front_firstburst_ci_hi", "front_firstburst_rho",
        "front_lead_kmday", "front_lead_ci_lo", "front_lead_ci_hi", "front_centroid_kmday",
        "front_lead_travel_km", "front_lead_travel_days",
        "n_secondary", "n_secondary_null_circ", "n_secondary_null_circ_p95",
        "n_secondary_null_dayshuf", "n_secondary_null_dayshuf_p95",
        "sf_speed_kmh_median", "sf_span_km_median", "sf_resid_km_median", "sf_behind_km_median",
        "n_reactivations", "frac_bursts_reactivation", "n_fam_reactivated", "frac_fam_reactivated",
        "reactivation_interval_d_p10", "reactivation_interval_d_p50", "reactivation_interval_d_p90",
        "reactivation_frac_le5d"]


def write_cat(ev, fam, path):
    """harness input: time, lat, lon, dep, mag, family_id (Bostock hypocentres, magnitudes unknown)."""
    d = ev.merge(fam, on="family_id", how="left")
    out = pd.DataFrame(dict(time=pd.to_datetime(d.time), lat=d.lat_bostock, lon=d.lon_bostock,
                            dep=d.depth_bostock, mag=np.nan, family_id=d.family_id.astype(int)))
    out = out[out.lat.notna()].sort_values("time")
    out.to_csv(path, index=False)
    return len(out)


def run_one(job):
    cat, tag, t0, t1, nd = job
    log = f"{TH}/log_{tag}.txt"
    env = dict(os.environ, PYTHONPATH="/opt/conda/lib/python3.13/site-packages", OMP_NUM_THREADS="2")
    with open(log, "w") as fh:
        r = subprocess.run([sys.executable, HARNESS, "--catalog", cat, "--tag", tag, "--mc", "none",
                            "--t0", t0, "--t1", t1, "--outdir", TH, "--null-draws", str(nd)],
                           stdout=fh, stderr=subprocess.STDOUT, env=env)
    return tag, r.returncode


def ref2012_ours():
    """OUR 2012 catalog at the frozen operating point, built exactly as Track C builds a month."""
    sys.path.insert(0, f"{EXT}/scripts")
    import c_common as C
    det = pd.read_parquet(f"{EXT}/results/e5r_detections_a1_mergeB_L12.parquet")
    det = det[det.L == 12.0].reset_index(drop=True)
    H = pd.read_parquet(f"{EXT}/results/a1_hf_ratio_a1_mergeB_L12.parquet")
    H = H[H.L == 12.0][["doy", "fam", "t", "hfr_PGC"]]
    det = det.merge(H, on=["doy", "fam", "t"], how="left")
    v = det.hfr_PGC.to_numpy()
    det["hf_pass"] = ~(np.isfinite(v) & (v > HF_CUT))
    share_ok = C.share_ok_matrix()
    fam_i = det.fam.to_numpy().astype(int)
    out = {}
    for lbl, T in (("ext", 8.25), ("core", 8.75)):
        sel = (det.net_mad >= T).to_numpy()
        kc = np.zeros(len(det), bool)
        if sel.any():
            kc[sel] = C.local_max_keep_cc(det.t.to_numpy()[sel], det.net_cc.to_numpy()[sel].astype(float),
                                          fam_i[sel], share_ok, 2.0)
        out[lbl] = sel & kc & det.hf_pass.to_numpy()
    E = det.assign(cat_ext=out["ext"], cat_core=out["core"])
    E["time"] = pd.to_datetime(E.t, unit="s")
    return E.rename(columns={"fam": "family_id"})


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--null-draws", type=int, default=20)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--major-only", action="store_true")
    ap.add_argument("--catalog", default=f"{TD}/lfe_extension_catalog_v1.parquet")
    ap.add_argument("--suffix", default="", help="appended to thomas_per_sse / thomas_bars and to every "
                                                 "per-run tag, so a second catalog version does not "
                                                 "overwrite the first")
    # Two window SETS are run and both are reported.  "main" is the pre-registered definition
    # (PNSN >= 50/day, gaps <= 20 d, span >= 10 d).  Its 20-day gap rule merges sub-episodes that are
    # separated by a week of quiet, and the harness's single-main-front window finder then returns
    # nothing ("no departure") -- so "sub" repeats the SAME definition with a 5-day clustering gap as a
    # REPORTED sensitivity, never as a bar.
    ap.add_argument("--windows", default="sse_windows.csv:main,sse_windows_sub.csv:sub")
    a = ap.parse_args()
    os.makedirs(TH, exist_ok=True)
    F = pd.read_csv(f"{IN}/e5r_families_2012.csv").rename(
        columns={"bkey": "family_id", "lat": "lat_bostock", "lon": "lon_bostock", "dep": "depth_bostock"})
    fam = F[["family_id", "lat_bostock", "lon_bostock", "depth_bostock"]]
    E = pd.read_parquet(a.catalog)
    E["family_id"] = E.family_id.astype(int)
    # a third, higher-PURITY leg, reported as a sensitivity and never as a bar: cat_core restricted to
    # detections registered by >= 5 distinct SITES.  The harness's leading edge is the 10th percentile
    # of event position, so at ~50 % precision it can be pinned by the background rather than by the
    # front ("no departure"); this leg is the detector-internal way to test that without conditioning
    # on tremor (which would be circular for a propagation speed).
    E["cat_core_hp"] = E.cat_core & (E.n_sites >= 5)
    jobs, index = [], []
    for spec in a.windows.split(","):
        wf, wset = spec.split(":")
        if not os.path.exists(f"{TD}/{wf}"):
            print(f"[d3] window set {wset}: {wf} not on disk, skipped"); continue
        W = pd.read_csv(f"{TD}/{wf}")
        if a.major_only:
            W = W[W.major]
        for r in W.itertuples():
            for tier, col in TIERS:
                sub = E[(E.date >= r.t0) & (E.date <= r.t1) & E[col]]
                tag = f"sse_{wset}_{r.sse_id}_{tier}{a.suffix}"
                p = f"{TH}/cat_{tag}.csv"
                n = write_cat(sub[["time", "family_id"]], fam, p)
                index.append(dict(tag=tag, window_set=wset, sse_id=r.sse_id, tier=f"cat_{tier}",
                                  t0=r.t0, t1=r.t1, major=bool(r.major), n_in_catalog=n,
                                  n_tremor_in_box=r.n_tremor_in_box, n_active_days=r.n_active_days))
                if n >= 50:
                    jobs.append((p, tag, r.t0, str((pd.Timestamp(r.t1) + pd.Timedelta(days=1)).date()),
                                 a.null_draws))
    # ---- the 2012 referee control, same code
    O12 = ref2012_ours()
    O12["cat_core_hp"] = O12.cat_core & (O12.n_sites >= 5)
    for tier, col in TIERS:
        sub = O12[O12[col]]
        tag = f"ref2012_ours_{tier}{a.suffix}"
        p = f"{TH}/cat_{tag}.csv"
        n = write_cat(sub[["time", "family_id"]], fam, p)
        index.append(dict(tag=tag, window_set="referee", sse_id="REF2012", tier=f"cat_{tier}", t0="2012-08-15",
                          t1="2012-10-16", major=True, n_in_catalog=n, n_tremor_in_box=np.nan,
                          n_active_days=np.nan))
        jobs.append((p, tag, "2012-08-15", "2012-10-16", a.null_draws))
    jobs.append((BOSTOCK, f"ref2012_bostock{a.suffix}", "2012-08-15", "2012-10-16", a.null_draws))
    index.append(dict(tag=f"ref2012_bostock{a.suffix}", window_set="referee", sse_id="REF2012", tier="bostock", t0="2012-08-15",
                      t1="2012-10-16", major=True, n_in_catalog=np.nan, n_tremor_in_box=np.nan,
                      n_active_days=np.nan))
    print(f"[d3] {len(jobs)} harness runs ({len(index)} windows indexed), {a.workers} workers", flush=True)
    with ThreadPoolExecutor(a.workers) as ex:
        for tag, rc in ex.map(run_one, jobs):
            print(f"   {tag} rc={rc}", flush=True)
    rows = []
    for r in index:
        p = f"{TH}/e5_thomas_summary_{r['tag']}.csv"
        d = dict(r)
        if os.path.exists(p):
            s = pd.read_csv(p)
            s = s[s.mc == "none"]
            if len(s):
                for k in KEYS:
                    if k in s.columns:
                        d[k] = s.iloc[0][k]
        rows.append(d)
    S = pd.DataFrame(rows)
    S.to_csv(f"{TD}/thomas_per_sse{a.suffix}.csv", index=False)
    # ---- the two pre-registered bars, on the 2012 referee
    bars = []
    bos = S[S.tag == f"ref2012_bostock{a.suffix}"]
    if len(bos):
        b = bos.iloc[0]
        for tier in ("core", "ext"):
            o = S[S.tag == f"ref2012_ours_{tier}{a.suffix}"]
            if not len(o):
                continue
            o = o.iloc[0]
            for k, lab, tol in (("front_firstburst_kmday", "main-front first-burst speed", 0.20),
                                ("reactivation_interval_d_p50", "re-rupture interval median", 0.05)):
                if not (np.isfinite(b.get(k, np.nan)) and np.isfinite(o.get(k, np.nan))):
                    continue
                rel = abs(o[k] - b[k]) / abs(b[k])
                bars.append(dict(bar=lab, tier=f"cat_{tier}", ours=o[k], bostock=b[k],
                                 rel_diff=rel, tol=tol, met=bool(rel <= tol)))
    B = pd.DataFrame(bars)
    B.to_csv(f"{TD}/thomas_bars{a.suffix}.csv", index=False)
    pd.set_option("display.width", 220)
    print("\n[BARS -- 2012 referee, our catalog vs Bostock's, same code]")
    print(B.round(4).to_string(index=False) if len(B) else "  (not evaluated)")
    print(f"\n[write] {TD}/thomas_per_sse{a.suffix}.csv, thomas_bars{a.suffix}.csv; "
          f"per-run files under {TH}/")


if __name__ == "__main__":
    main()
