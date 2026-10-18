#!/usr/bin/env python
"""Track A step 1(c) -- define the TWO operating thresholds by FALSE ALARMS ONLY (Merlin plan-review F7/F1), then report
recall and precision at them.  Nothing here may be adjusted after seeing recall.

  T_ext  = smallest T on a 0.25-MAD grid with quiet-day FA <= 1.0 AND post-SSE (doy 270-279) FA <= 1.5 per template-day
  T_core = smallest T with quiet-day FA <= 0.1 per template-day
Quiet days = 2012 doy 228-242: no Bostock event and (verified here) ZERO PNSN tremor detections in the family box.
FA per template-day = declustered detections of that family's network / (families x days) over the window; the detections
parquet reproduces the budget counts exactly at every tabulated T (checked), so the 0.25-MAD grid is exact, not interpolated.
Recall = own-family, +-2 s, one-to-one greedy by time, Bostock 2012 TEST half.  Precision = own / all detections on ETS days
(>= 100 Bostock events) and over the whole run.  `--min-sites 2` (Merlin F8) is applied post hoc via the recorded n_sites /
n_site column: the gate only decides whether a family-day is scored, so filtering afterwards is exactly equivalent."""
from __future__ import annotations
import argparse, json, hashlib, os
import numpy as np, pandas as pd

REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
IN = f"{REP}/results"
OUT = "/home/jovyan/tremorferometry/lfe_extension/results"
QUIET = list(range(228, 243)); POST = list(range(270, 280))
PNSN = "/home/jovyan/tremorferometry/catalogs/pnsn_tremor_cascadia_full.csv"


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


def own_recall(det, test, T, extra=None):
    D = det[det.net_mad >= T] if extra is None else det[(det.net_mad >= T) & extra]
    hit = np.zeros(len(test), bool)
    for fk, g in test.groupby("bkey"):
        td = np.sort(D[D.fam == float(fk)].t.to_numpy())
        hit[test.index.get_indexer(g.index)] = one2one(g.t_origin.to_numpy(), td) <= 2.0
    return hit


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for c in iter(lambda: f.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()[:16]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--det", default=f"{IN}/e5r_detections_e6_mergeB.parquet")
    ap.add_argument("--bud", default=f"{IN}/e5_budget_e6_mergeB.csv")
    ap.add_argument("--tag", default="e6_mergeB")
    ap.add_argument("--L", type=float, default=8.0)
    ap.add_argument("--tmin", type=float, default=7.0); ap.add_argument("--tmax", type=float, default=12.0)
    ap.add_argument("--step", type=float, default=0.25)
    ap.add_argument("--min-sites", type=int, default=0, help="post-hoc F8 filter on the recorded n_sites")
    ap.add_argument("--out-prefix", default=f"{OUT}/a1_operating")
    # step 1(b): apply the adopted HF-ratio veto before the FA count (detections whose PGC-site 10-19 Hz / 2-8 Hz energy ratio
    # in the scored 8-s window exceeds the build-half own-matched 99th percentile are dropped; an undefined ratio is never vetoed)
    ap.add_argument("--hf-parquet", default="")
    ap.add_argument("--hf-col", default="hfr_PGC")
    ap.add_argument("--hf-cut", type=float, default=np.nan)
    a = ap.parse_args()
    det = pd.read_parquet(a.det); det = det[det.L == a.L].reset_index(drop=True)
    bud = pd.read_csv(a.bud); bud = bud[bud.L == a.L]
    nsc = "n_sites" if "n_sites" in det.columns else ("n_site" if "n_site" in det.columns else None)
    if a.min_sites and nsc:
        det = det[det[nsc] >= a.min_sites].reset_index(drop=True)
        bud = bud[bud[nsc.replace("n_sites", "n_site") if nsc == "n_sites" and "n_site" in bud.columns else nsc] >= a.min_sites] \
            if any(c in bud.columns for c in ("n_sites", "n_site")) else bud
    hf_info = None
    if a.hf_parquet:
        H = pd.read_parquet(a.hf_parquet)
        H = H[H.L == a.L][["doy", "fam", "t", a.hf_col]]
        det = det.merge(H, on=["doy", "fam", "t"], how="left")
        cut = a.hf_cut
        if not np.isfinite(cut):
            raise SystemExit("--hf-cut required with --hf-parquet")
        v = det[a.hf_col].to_numpy()
        drop = np.isfinite(v) & (v > cut)
        hf_info = dict(col=a.hf_col, cut=float(cut), n_dropped=int(drop.sum()), frac_dropped=float(drop.mean()),
                       frac_ratio_defined=float(np.isfinite(v).mean()))
        print(f"[HF veto] {a.hf_col} > {cut:g}: {drop.sum()} of {len(det)} detections dropped "
              f"({drop.mean():.1%}); ratio defined for {np.isfinite(v).mean():.1%}")
        det = det[~drop].reset_index(drop=True)
    Tgrid = sorted(bud["T"].unique())
    fam_days = bud[bud["T"] == Tgrid[0]][["doy", "fam"]].drop_duplicates()
    nq = int(fam_days.doy.isin(QUIET).sum()); npo = int(fam_days.doy.isin(POST).sum()); nall = len(fam_days)
    gt = pd.read_parquet(f"{IN}/e5r_gt_bostock2012.parquet").reset_index(drop=True)
    gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    per_day = gt.groupby("doy").size(); ETS = set(per_day[per_day >= 100].index)
    nets = int(fam_days.doy.isin(ETS).sum())
    test = gt[(gt.half == "test") & (~gt.guard)].copy().reset_index(drop=True)
    # PNSN verification of the quiet window
    tr = pd.read_csv(PNSN); tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    F = pd.read_csv(f"{IN}/e5r_families_2012.csv")
    box = tr[(tr.time.dt.year == 2012) & tr.lat.between(F.lat.min() - .2, F.lat.max() + .2) &
             tr.lon.between(F.lon.min() - .2, F.lon.max() + .2)]
    qtr = box[box.time.dt.dayofyear.isin(QUIET)]
    print(f"[quiet] doy {QUIET[0]}-{QUIET[-1]}: {len(qtr)} PNSN tremor detections in the family box "
          f"(lat {F.lat.min():.2f}-{F.lat.max():.2f}, lon {F.lon.min():.2f}-{F.lon.max():.2f}) -> "
          f"{'all 15 days qualify as quiet' if len(qtr) == 0 else 'NOT ALL QUIET'}")
    rows = []
    for T in np.arange(a.tmin, a.tmax + 1e-9, a.step):
        D = det[det.net_mad >= T]
        rows.append(dict(T=round(float(T), 4),
                         quiet_fa=len(D[D.doy.isin(QUIET)]) / nq, post_fa=len(D[D.doy.isin(POST)]) / npo,
                         ets_det_per_tpl_day=len(D[D.doy.isin(ETS)]) / max(nets, 1), n_catalog=len(D),
                         ets_precision_own=float(D[D.doy.isin(ETS)].match_own.mean()) if len(D[D.doy.isin(ETS)]) else np.nan,
                         ets_precision_any=float(D[D.doy.isin(ETS)].match_any.mean()) if len(D[D.doy.isin(ETS)]) else np.nan,
                         precision_own_all=float(D.match_own.mean()) if len(D) else np.nan))
    G = pd.DataFrame(rows)
    G.to_csv(f"{a.out_prefix}_fa_vs_T_{a.tag}.csv", index=False)
    ext = G[(G.quiet_fa <= 1.0) & (G.post_fa <= 1.5)]
    core = G[G.quiet_fa <= 0.1]
    T_ext = float(ext["T"].min()) if len(ext) else np.nan
    T_core = float(core["T"].min()) if len(core) else np.nan
    print(f"\n[thresholds, FA-ONLY]  T_ext = {T_ext}  (quiet FA {float(G[G['T']==T_ext].quiet_fa.iloc[0]):.3f}, "
          f"post-SSE FA {float(G[G['T']==T_ext].post_fa.iloc[0]):.3f})   "
          f"T_core = {T_core}  (quiet FA {float(G[G['T']==T_core].quiet_fa.iloc[0]):.4f})")
    out = dict(tag=a.tag, L=a.L, T_ext=T_ext, T_core=T_core, min_sites_applied=a.min_sites,
               n_quiet_family_days=nq, n_post_family_days=npo, n_ets_family_days=nets, n_family_days=nall,
               pnsn_tremor_in_quiet_window=int(len(qtr)),
               det_file=a.det, det_sha256_16=sha(a.det), hf_veto=hf_info)
    for nm, T in (("T_ext", T_ext), ("T_core", T_core)):
        if not np.isfinite(T):
            continue
        g = G[G["T"] == T].iloc[0]
        hit = own_recall(det, test, T)
        blk = {}
        for lo, hi, lbl in [(-9, 1.6, "M<1.6"), (1.6, 1.8, "M1.6-1.8"), (1.8, 2.0, "M1.8-2.0"), (2.0, 9, "M>=2.0")]:
            k = (test.mag >= lo) & (test.mag < hi); blk[lbl] = float(hit[k.to_numpy()].mean())
        out[nm + "_report"] = dict(T=T, quiet_fa=float(g.quiet_fa), post_fa=float(g.post_fa),
                                   ets_det_per_tpl_day=float(g.ets_det_per_tpl_day), n_catalog=int(g.n_catalog),
                                   recall_own_test=float(hit.mean()), recall_by_mag=blk,
                                   precision_own_ets=float(g.ets_precision_own), precision_any_ets=float(g.ets_precision_any),
                                   precision_own_all=float(g.precision_own_all))
        print(f"  {nm} = {T}: own recall (test half, one-to-one) {hit.mean():.4f}; ETS own precision {g.ets_precision_own:.3f}; "
              f"any {g.ets_precision_any:.3f}; catalog {int(g.n_catalog)}; by mag " +
              ", ".join(f"{k} {v:.3f}" for k, v in blk.items()))
    # site-rolled / data on the quiet window, at the budget grid points bracketing each threshold
    if "n_det_roll_site" in bud.columns:
        sr = []
        for T in Tgrid:
            b = bud[(bud["T"] == T) & bud.doy.isin(QUIET)]
            sr.append(dict(T=T, quiet_data=float(b.n_det.mean()), quiet_roll_site=float(b.n_det_roll_site.mean()),
                           site_over_data=float(b.n_det_roll_site.mean()) / max(float(b.n_det.mean()), 1e-9),
                           quiet_roll_sta=float(b.n_det_roll.mean())))
        S = pd.DataFrame(sr); S.to_csv(f"{a.out_prefix}_siterolled_{a.tag}.csv", index=False)
        out["site_rolled_grid"] = S.to_dict("records")
        print("\n[site-rolled control, quiet window, budget grid]"); print(S.round(3).to_string(index=False))
    json.dump(out, open(f"{a.out_prefix}_{a.tag}.json", "w"), indent=1, default=float)
    print(f"\n[write] {a.out_prefix}_{a.tag}.json, _fa_vs_T_{a.tag}.csv")


if __name__ == "__main__":
    main()
