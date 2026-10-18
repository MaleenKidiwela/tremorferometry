#!/usr/bin/env python
"""STAGE 2 -- the LITERATURE-GROUNDED acceptance procedure, applied to one chain.

Order, and only these steps (results/SCREEN_VALIDATION_NORTH.md sec.4 is the source for each):

  1. DEDUPE.  The one family-acceptance rule positively attested for the Cascadia catalogues:
     merge templates whose declustered hit sets overlap by >= 10 % of the smaller set OR whose
     stacks cross-correlate >= 0.3, keeping the one with more detections (tie-break: higher median
     cc).  Bostock 2012 culled 130 templates to 93 geographically independent families; Plourde
     2015 used "> 10 % of detections shared" or "cc > 0.3"; Chestler 2017 used R_com > 0.5 and
     kept the higher-median-cc family.
  2. LOCATION.  Located on >= 5 stations -- the locator's own bar, and the ONE clause measured to
     be justified (86 % of the 58 rediscovered Bostock families pass it; `sf_depthcheck.py`
     reproduces agent Q's number on this population).  Railing on a face of the search grid is
     REPORTED as a diagnostic, never a rejection: the interior clause rejects 62 % of Bostock's
     own rediscovered families because the rails are a defect of our grid.
     The DEPTH WINDOW (15-58 km) is likewise reported, NOT applied.  `sf_depthcheck.py` measured
     it: on the north bed the window is vacuous (that locator's grid is 18-56 km, so 100 % pass
     for a reason that has nothing to do with the families), while its actual effect on the
     Washington grid is to delete solutions sitting on the bottom depth rail (dep = 59 km) and
     nothing else -- 100 % of its rejections on every Washington chain are that rail, 0 are
     genuinely interior solutions outside the window.  The transferable form of the test, "reject
     a depth-rail solution", keeps only 66 % of the located true north families.  It is therefore
     a grid-rail clause in disguise and may not reject: `ACCEPT` is dedupe + location, and
     `accept_depth_window` is carried beside it so every count can be given both ways.
  3. CONFIRMATION on a sealed window -- `sf_confirm.py`, with the rule calibrated in stage 1.

No quiet-day false-alarm gate, no minimum-day rule, no Slab2 agreement, no interior clause and no
split-half stability test: none of them appears in any published Cascadia LFE catalogue, and each
was measured on the north bed to reject a large share of Bostock's own families.

  python scripts/sf_accept.py --det ... --tpl ... --loc ... --out results/sf_final/<tag> --tag ...
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
sys.path.insert(0, HERE)
import sf_common as SF                                                      # noqa: E402

DEP_MIN, DEP_MAX = 15.0, 58.0        # REPORTED, never applied -- see the module docstring
WA_RAILS = (15.0, 59.0)              # the wide grid's depth faces (arange(15, 60.1, 2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--det", required=True)
    ap.add_argument("--tpl", required=True)
    ap.add_argument("--loc", required=True)
    ap.add_argument("--out", required=True, help="output prefix (directory + tag)")
    ap.add_argument("--tag", required=True)
    ap.add_argument("--no-cc", action="store_true", help="skip the stack-cc leg (smoke only)")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)

    D = pd.read_parquet(a.det)
    TPL = SF.load_tpl(a.tpl)
    grps = sorted(TPL)
    if a.limit:
        grps = grps[:a.limit]
        TPL = {g: TPL[g] for g in grps}
        D = D[D.grp.isin(set(grps))]
    SF.log(f"[{a.tag}] {len(grps)} grown templates, {len(D)} detections")

    keep, DD, pairs = SF.dedupe(D, TPL, grps, do_cc=not a.no_cc)
    DD.to_csv(f"{a.out}_dedupe.csv", index=False)
    SF.log(f"[{a.tag}] dedupe: {pairs['n_templates']} -> {pairs['n_kept']} families "
           f"(pairs merged by shared detections {pairs['n_pairs_share']}, by stack cc "
           f"{pairs['n_pairs_cc']}; max stack cc {pairs['max_cc_overall']:.3f}, "
           f"max shared fraction {pairs['max_share_overall']:.3f})")

    L = pd.read_csv(a.loc).set_index("grp")
    NDET = DD.set_index("grp").n_det.to_dict()
    rows = []
    for g in grps:
        r = L.loc[g] if g in L.index else None
        rec = dict(grp=int(g), kept_dedupe=bool(g in set(keep)),
                   n_det=int(NDET.get(g, 0)),
                   n_days=int(r.n_days) if r is not None and "n_days" in L.columns else 0,
                   located=bool(r.located) if r is not None else False)
        if rec["located"]:
            rec.update(n_sta=int(r.n_sta), lat=float(r.lat), lon=float(r.lon), dep=float(r.dep),
                       rms=float(r.rms), ctrl_pass=bool(r.ctrl_pass),
                       edge=bool(r.edge) if "edge" in L.columns else False,
                       slab2=float(r.slab2) if "slab2" in L.columns else np.nan,
                       dz=float(r.dz) if "dz" in L.columns else np.nan)
            rec["depth_ok"] = bool(DEP_MIN <= rec["dep"] <= DEP_MAX)
            rec["dep_rail"] = bool(rec["dep"] in WA_RAILS)
        else:
            rec.update(n_sta=int(r.n_sta) if r is not None else 0, lat=np.nan, lon=np.nan,
                       dep=np.nan, rms=np.nan, ctrl_pass=False, edge=False, slab2=np.nan,
                       dz=np.nan, depth_ok=False, dep_rail=False)
        rec["ACCEPT"] = bool(rec["kept_dedupe"] and rec["located"])
        rec["accept_depth_window"] = bool(rec["ACCEPT"] and rec["depth_ok"])
        rows.append(rec)
    A = pd.DataFrame(rows)
    A.to_csv(f"{a.out}_accept.csv", index=False)
    A[A.ACCEPT][["grp"]].to_csv(f"{a.out}_accept_grps.csv", index=False)
    A[A.accept_depth_window][["grp"]].to_csv(f"{a.out}_accept_grps_depthwindow.csv", index=False)

    ded = A[A.kept_dedupe]
    summ = dict(tag=a.tag, n_grown=int(len(A)),
                n_dedupe=int(A.kept_dedupe.sum()),
                n_dedupe_located=int((A.kept_dedupe & A.located).sum()),
                n_accept=int(A.ACCEPT.sum()),
                n_accept_depth_window=int(A.accept_depth_window.sum()),
                n_accept_on_depth_rail=int((A.ACCEPT & A.dep_rail).sum()),
                n_accept_on_bottom_rail=int((A.ACCEPT & (A.dep == WA_RAILS[1])).sum()),
                n_accept_on_floor_rail=int((A.ACCEPT & (A.dep == WA_RAILS[0])).sum()),
                n_accept_interior=int((A.ACCEPT & ~A.edge).sum()),
                n_accept_edge=int((A.ACCEPT & A.edge).sum()),
                n_accept_ctrl_pass=int((A.ACCEPT & A.ctrl_pass).sum()),
                frac_edge_of_located=float((ded.located & ded.edge).sum()
                                           / max(int(ded.located.sum()), 1)),
                dedupe=pairs, dep_min=DEP_MIN, dep_max=DEP_MAX)
    json.dump(summ, open(f"{a.out}_accept_summary.json", "w"), indent=1, default=float)
    print(json.dumps(summ, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
