#!/usr/bin/env python
"""Track C — RE-FLAG the two tiers on processed months from results/trackC/tiers_by_epoch.csv,
WITHOUT re-running the matched filter.

`net_mad` is stored, so a change of threshold is pure post-processing.  What is rewritten:
    tier_ext, tier_core          from the month's NETWORK-EPOCH T_ext / T_core
    keep_xfam_ext, keep_xfam_core, keep_mad_ext, keep_mad_core
                                 the cross-family assignment is evaluated ON the thresholded set, so a
                                 new threshold changes which detections can suppress which - it MUST
                                 be recomputed, with the same frozen rule (raw weighted cc among
                                 families sharing >= 50 % of their 2012 stations, +-2 s)
    cat_ext, cat_core            tier AND hf_pass AND keep_xfam_*
    T_core_used, T_ext_used, tier_epoch   provenance on every row
and in family_days_<YYYY-MM>.parquet, n_det_ext / n_det_core.
NOT touched: net_mad, net_cc, rel_amp, hfr, hf_cut, hf_pass, the templates, the stacks, the tremor /
EQ / quiet flags.

  python scripts/c_reflag_tiers.py --all --dry-run
  python scripts/c_reflag_tiers.py --all
"""
from __future__ import annotations
import os, sys, glob, json, argparse, subprocess
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C

PY = "/opt/conda/bin/python"
HERE = os.path.dirname(os.path.abspath(__file__))


def reflag(ym, tiers, share_ok, dry=False, out_dir=None):
    OUT = out_dir or C.TRC
    p = f"{OUT}/detections_{ym}.parquet"
    if not os.path.exists(p):
        return None
    ep = C.epoch_of_month(ym)
    if ep not in tiers:
        return dict(ym=ym, epoch=ep, status="no tier row", n=0)
    T_core, T_ext = tiers[ep]
    D = pd.read_parquet(p)
    if not len(D):
        return dict(ym=ym, epoch=ep, status="empty", n=0)
    before = dict(cat_ext=int(D.cat_ext.sum()), cat_core=int(D.cat_core.sum()),
                  tier_ext=int(D.tier_ext.sum()), tier_core=int(D.tier_core.sum()),
                  T_core=float(D.T_core_used.iloc[0]) if "T_core_used" in D else np.nan,
                  T_ext=float(D.T_ext_used.iloc[0]) if "T_ext_used" in D else np.nan)
    tier_ext = (D.net_mad >= T_ext).to_numpy()
    tier_core = (D.net_mad >= T_core).to_numpy()
    fam_i = D.fam.to_numpy().astype(int)
    t = D.t.to_numpy(); cc = D.net_cc.to_numpy().astype(float); mad = D.net_mad.to_numpy()
    keeps = {}
    for lbl, sel in (("ext", tier_ext), ("core", tier_core)):
        kc = np.zeros(len(D), bool); km = np.zeros(len(D), bool)
        if sel.any():
            kc[sel] = C.local_max_keep_cc(t[sel], cc[sel], fam_i[sel], share_ok, 2.0)
            km[sel] = C.local_max_keep(t[sel], mad[sel], 2.0)
        keeps[lbl] = (kc, km)
    cat_ext = tier_ext & D.hf_pass.to_numpy() & keeps["ext"][0]
    cat_core = tier_core & D.hf_pass.to_numpy() & keeps["core"][0]
    res = dict(ym=ym, epoch=ep, status="reflagged", n=len(D), T_core=T_core, T_ext=T_ext,
               T_core_before=before["T_core"], T_ext_before=before["T_ext"],
               tier_ext_before=before["tier_ext"], tier_ext_after=int(tier_ext.sum()),
               tier_core_before=before["tier_core"], tier_core_after=int(tier_core.sum()),
               cat_ext_before=before["cat_ext"], cat_ext_after=int(cat_ext.sum()),
               cat_core_before=before["cat_core"], cat_core_after=int(cat_core.sum()))
    if dry:
        return res
    D["tier_ext_prev"] = D.tier_ext; D["tier_core_prev"] = D.tier_core
    D["tier_ext"] = tier_ext; D["tier_core"] = tier_core
    D["keep_xfam_ext"], D["keep_mad_ext"] = keeps["ext"]
    D["keep_xfam_core"], D["keep_mad_core"] = keeps["core"]
    D["cat_ext"] = cat_ext; D["cat_core"] = cat_core
    D["T_core_used"] = T_core; D["T_ext_used"] = T_ext; D["tier_epoch"] = ep
    tmp = p + ".tmp"; D.to_parquet(tmp, index=False); os.replace(tmp, p)
    fp = f"{OUT}/family_days_{ym}.parquet"
    if os.path.exists(fp):
        F = pd.read_parquet(fp)
        ne = D[tier_ext].groupby(["fam", "date"]).size()
        nc = D[tier_core].groupby(["fam", "date"]).size()
        idx = pd.MultiIndex.from_arrays([F.fam.astype(int), F.date])
        F["n_det_ext"] = ne.reindex(idx).fillna(0).to_numpy().astype(int)
        F["n_det_core"] = nc.reindex(idx).fillna(0).to_numpy().astype(int)
        F["T_core_used"] = T_core; F["T_ext_used"] = T_ext; F["tier_epoch"] = ep
        tmp = fp + ".tmp"; F.to_parquet(tmp, index=False); os.replace(tmp, fp)
    mp = f"{OUT}/run_{ym}.json"
    if os.path.exists(mp):
        m = json.load(open(mp))
        m["tiers_reflagged"] = res
        m["T_core"], m["T_ext"] = T_core, T_ext
        m["n_cat_ext"], m["n_cat_core"] = res["cat_ext_after"], res["cat_core_after"]
        json.dump(m, open(mp, "w"), indent=1, default=str)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("months", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-qc", action="store_true")
    ap.add_argument("--out-dir", default="")
    a = ap.parse_args()
    OUT = a.out_dir or C.TRC
    tiers = C.tiers_by_epoch()
    if not tiers:
        raise SystemExit("[tiers] results/trackC/tiers_by_epoch.csv is missing — "
                         "run scripts/c_tiers_by_epoch.py first")
    print("[tiers] " + "; ".join(f"{k}: T_core {v[0]} / T_ext {v[1]}" for k, v in sorted(tiers.items())))
    share_ok = C.share_ok_matrix()
    months = a.months
    if a.all or not months:
        months = sorted(os.path.basename(f).rsplit("_", 1)[1][:-8]
                        for f in glob.glob(f"{OUT}/detections_*.parquet"))
    rows = []
    for ym in months:
        r = reflag(ym, tiers, share_ok, a.dry_run, OUT)
        if not r:
            continue
        rows.append(r)
        if r.get("status") == "reflagged":
            print(f"  {ym} [{r['epoch']}] T {r['T_core_before']}/{r['T_ext_before']} -> "
                  f"{r['T_core']}/{r['T_ext']}: cat_core {r['cat_core_before']} -> {r['cat_core_after']}, "
                  f"cat_ext {r['cat_ext_before']} -> {r['cat_ext_after']}", flush=True)
            if not a.dry_run and not a.no_qc:
                subprocess.call([PY, f"{HERE}/c_qc_month.py", ym], cwd=C.EXT, stdout=subprocess.DEVNULL)
    if rows and not a.dry_run:
        D = pd.DataFrame(rows)
        D.to_csv(f"{OUT}/tier_reflag_log.csv", index=False)
        g = D[D.status == "reflagged"].groupby("epoch").agg(
            months=("ym", "count"), cat_core_before=("cat_core_before", "sum"),
            cat_core_after=("cat_core_after", "sum"), cat_ext_before=("cat_ext_before", "sum"),
            cat_ext_after=("cat_ext_after", "sum"))
        print("\n" + g.to_string())
        print(f"[tiers] -> {OUT}/tier_reflag_log.csv")


if __name__ == "__main__":
    main()
