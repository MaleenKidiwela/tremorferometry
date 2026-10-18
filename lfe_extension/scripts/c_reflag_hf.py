#!/usr/bin/env python
"""Track C — RE-FLAG the HF veto on months that are already processed, WITHOUT re-running the
matched filter.

The HF ratio `hfr` and the HF station's native rate `hf_fs` are stored on every detection, so a new
per-epoch cut (results/a4_hf_veto_cuts.csv) only changes `hf_cut`, `hf_pass` and the two catalog
columns `cat_ext` / `cat_core`.  Nothing else in the file is touched: net_mad, net_cc, rel_amp, the
tiers, the cross-family assignment, the tremor / EQ / quiet flags and every stack stay as computed.
The previous flags are kept as `hf_pass_prev` / `hf_cut_prev` so the change is auditable, and the
month's QC is regenerated.

  python scripts/c_reflag_hf.py --all            # every processed month
  python scripts/c_reflag_hf.py 2018-06 2018-07  # named months
  python scripts/c_reflag_hf.py --all --dry-run  # report what would change
"""
from __future__ import annotations
import os, sys, glob, json, argparse, subprocess
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C

PY = "/opt/conda/bin/python"
HERE = os.path.dirname(os.path.abspath(__file__))


def reflag(ym, HFC, hf_sta, default_cut, dry=False, out_dir=None):
    OUT = out_dir or C.TRC
    p = f"{OUT}/detections_{ym}.parquet"
    if not os.path.exists(p):
        return None
    D = pd.read_parquet(p)
    if not len(D):
        return dict(ym=ym, n=0)
    if "hfr" not in D or "hf_fs" not in D:
        raise SystemExit(f"[reflag] {p} has no hfr/hf_fs column — it predates the stored ratio; re-run the month")
    old_pass = D.hf_pass.to_numpy().copy()
    old_cut = D.hf_cut.to_numpy().copy() if "hf_cut" in D else np.full(len(D), default_cut)
    hc = [C.hf_cut_for(HFC, hf_sta, d_, f_, default_cut)
          for d_, f_ in zip(D.date.to_numpy(), D.hf_fs.to_numpy())]
    cut = np.array([x[0] for x in hc], float)
    src = [x[1] for x in hc]
    new_pass = ~(np.isfinite(D.hfr.to_numpy()) & (D.hfr.to_numpy() > cut))
    res = dict(ym=ym, n=len(D), n_changed=int((new_pass != old_pass).sum()),
               cut_before=float(np.nanmedian(old_cut)), cut_after=float(np.nanmedian(cut)),
               cat_core_before=int(D.cat_core.sum()), cat_ext_before=int(D.cat_ext.sum()))
    if dry:
        res["cat_core_after"] = int((D.tier_core & new_pass & D.keep_xfam_core).sum())
        res["cat_ext_after"] = int((D.tier_ext & new_pass & D.keep_xfam_ext).sum())
        return res
    D["hf_pass_prev"] = old_pass
    D["hf_cut_prev"] = old_cut
    D["hf_cut"] = cut
    D["hf_cut_src"] = src
    D["hf_pass"] = new_pass
    D["cat_ext"] = D.tier_ext & D.hf_pass & D.keep_xfam_ext
    D["cat_core"] = D.tier_core & D.hf_pass & D.keep_xfam_core
    tmp = p + ".tmp"
    D.to_parquet(tmp, index=False)
    os.replace(tmp, p)
    res["cat_core_after"] = int(D.cat_core.sum()); res["cat_ext_after"] = int(D.cat_ext.sum())
    mp = f"{OUT}/run_{ym}.json"
    if os.path.exists(mp):
        m = json.load(open(mp))
        m["hf_reflagged"] = res
        m["n_cat_ext"] = res["cat_ext_after"]; m["n_cat_core"] = res["cat_core_after"]
        json.dump(m, open(mp, "w"), indent=1, default=str)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("months", nargs="*")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--hf-sta", default="PGC")
    ap.add_argument("--no-qc", action="store_true")
    ap.add_argument("--out-dir", default="")
    a = ap.parse_args()
    OUT = a.out_dir or C.TRC
    fz = C.frozen()
    default_cut = float(fz.get("hf_veto", {}).get("cut", np.nan))
    HFC = C.hf_cuts()
    if not HFC:
        print(f"[reflag] results/a4_hf_veto_cuts.csv is not there yet — nothing to re-flag "
              f"(the frozen cut {default_cut:.4f} stays in force).")
        return
    print(f"[reflag] {len(HFC)} per-epoch cuts; frozen fallback {default_cut:.4f}")
    months = a.months
    if a.all or not months:
        months = sorted(os.path.basename(f).rsplit("_", 1)[1][:-8] for f in glob.glob(f"{OUT}/detections_*.parquet"))
    rows = []
    for ym in months:
        r = reflag(ym, HFC, a.hf_sta, default_cut, a.dry_run, OUT)
        if r:
            rows.append(r)
            print(f"  {ym}: {r.get('n_changed', 0)} of {r['n']} flags change, cut "
                  f"{r['cut_before']:.4f} -> {r['cut_after']:.4f}, cat_core "
                  f"{r['cat_core_before']} -> {r.get('cat_core_after')}", flush=True)
            if not a.dry_run and not a.no_qc:
                subprocess.call([PY, f"{HERE}/c_qc_month.py", ym], cwd=C.EXT,
                                stdout=subprocess.DEVNULL)
    if rows and not a.dry_run:
        pd.DataFrame(rows).to_csv(f"{OUT}/hf_reflag_log.csv", index=False)
        print(f"[reflag] -> {OUT}/hf_reflag_log.csv (QC regenerated for each month)")


if __name__ == "__main__":
    main()
