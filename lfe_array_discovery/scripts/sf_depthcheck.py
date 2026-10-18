#!/usr/bin/env python
"""Is the DEPTH WINDOW allowed to reject anything?  Measured on Bostock's own families.

The stage-2 acceptance procedure applied a "physically sensible depth, 15-58 km" bar.  Under the
project's binding rule (`SCREEN_VALIDATION_NORTH.md` sec.3), no criterion may REJECT a candidate
until it has been shown to pass >= 85 % of the 58 blind-rediscovered Bostock families on the 2005
bed.  This module measures it there, and measures the >= 5-station location bar on the same
population, and reports what the window actually removes on the Washington beds.

Two things have to be separated, because the two beds have different grids:
  * the LITERAL window, 15-58 km.  The north locator's grid is 18-56 km, so on that bed the window
    cannot reject anything at all -- it passes 100 % of the located true templates vacuously.  A
    100 % that is produced by the grid, not by the families, is not a validation.
  * what the window actually DOES on the Washington grid (15-60 km at 2-km spacing, i.e. nodes
    15, 17, ..., 59): its only effect there is to delete solutions sitting on the BOTTOM DEPTH
    RAIL, dep = 59.  The transferable form of the test is therefore "reject a solution that rails
    on a depth face of the search grid", and THAT is measured on the north bed against its own
    rails (18 km and 56 km).

`results/repro/classify15.csv` carries `cls` (kept_family / hidden_family / candidate_new /
unlocatable) beside the k_locate solution for every one of the 1,460 grown 2005 templates.
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
sys.path.insert(0, HERE)
import sf_common as SF                                                      # noqa: E402

CLS = f"{SF.RES}/repro/classify15.csv"
NORTH_RAILS = (18.0, 56.0)          # k_locate grid floor / bottom (arange(18, 56.1, 2))
WA_RAILS = (15.0, 59.0)             # n_locate grid floor / bottom (arange(15, 60.1, 2))
DEP_MIN, DEP_MAX = 15.0, 58.0
PASS_FRAC = 0.85


def north():
    C = pd.read_csv(CLS)
    t = C[C.cls == "kept_family"].copy()
    z = C[C.cls == "unlocatable"].copy()
    nfam = t.kept_bkey.nunique()

    def fam(mask):
        return int(t[mask].kept_bkey.nunique())

    loc = t.located == True                                                 # noqa: E712
    locz = z.located == True                                                # noqa: E712
    lit = loc & (t.dep >= DEP_MIN) & (t.dep <= DEP_MAX)
    norail = loc & (t.dep > NORTH_RAILS[0]) & (t.dep < NORTH_RAILS[1])
    litz = locz & (z.dep >= DEP_MIN) & (z.dep <= DEP_MAX)
    norailz = locz & (z.dep > NORTH_RAILS[0]) & (z.dep < NORTH_RAILS[1])
    out = dict(
        n_true_templates=int(len(t)), n_true_families=int(nfam),
        n_noise_templates=int(len(z)),
        grid_dep_min=float(C[C.located == True].dep.min()),                 # noqa: E712
        grid_dep_max=float(C[C.located == True].dep.max()),                 # noqa: E712
        located=dict(tpl=int(loc.sum()), fam=fam(loc), frac_fam=fam(loc) / nfam,
                     noise_tpl=int(locz.sum()), noise_frac=float(locz.mean())),
        literal_window=dict(lo=DEP_MIN, hi=DEP_MAX, tpl=int(lit.sum()), fam=fam(lit),
                            frac_fam=fam(lit) / nfam,
                            frac_fam_of_located=fam(lit) / max(fam(loc), 1),
                            noise_tpl=int(litz.sum()), noise_frac=float(litz.mean()),
                            vacuous=bool(int(lit.sum()) == int(loc.sum()))),
        no_depth_rail=dict(rails=list(NORTH_RAILS), tpl=int(norail.sum()), fam=fam(norail),
                           frac_fam=fam(norail) / nfam,
                           frac_fam_of_located=fam(norail) / max(fam(loc), 1),
                           noise_tpl=int(norailz.sum()), noise_frac=float(norailz.mean())),
        rail_counts=dict(floor=int((loc & (t.dep == NORTH_RAILS[0])).sum()),
                         bottom=int((loc & (t.dep == NORTH_RAILS[1])).sum())),
        depth_quantiles_located_true={str(q): float(np.quantile(t[loc].dep, q))
                                      for q in (0.05, 0.25, 0.5, 0.75, 0.95)},
    )
    # what window WOULD pass 85 % of the located true families?
    sweep = []
    for lo, hi in [(15, 58), (20, 55), (25, 50), (28, 45), (18, 56), (16, 58), (14, 60)]:
        m = loc & (t.dep >= lo) & (t.dep <= hi)
        sweep.append(dict(lo=lo, hi=hi, fam=fam(m), frac_fam_of_located=fam(m) / max(fam(loc), 1),
                          frac_fam_of_all=fam(m) / nfam))
    out["window_sweep"] = sweep
    out["verdict_literal"] = ("JUSTIFIED-BUT-VACUOUS"
                              if out["literal_window"]["vacuous"] else
                              ("JUSTIFIED" if out["literal_window"]["frac_fam_of_located"]
                               >= PASS_FRAC else "UNJUSTIFIED"))
    out["verdict_transferable"] = ("JUSTIFIED"
                                   if out["no_depth_rail"]["frac_fam_of_located"] >= PASS_FRAC
                                   else "UNJUSTIFIED")
    out["verdict_min5sta"] = ("JUSTIFIED" if out["located"]["frac_fam"] >= PASS_FRAC
                              else "UNJUSTIFIED")
    return out


def wa():
    """What the window removes on each Washington chain, from the saved acceptance tables."""
    out = {}
    for p in sorted(glob.glob(f"{SF.OUT}/*_accept.csv")):
        key = os.path.basename(p)[:-len("_accept.csv")]
        A = pd.read_csv(p)
        D = A[A.kept_dedupe.astype(bool) & A.located.astype(bool)]
        if not len(D):
            continue
        bad = D[~((D.dep >= DEP_MIN) & (D.dep <= DEP_MAX))]
        out[key] = dict(
            n_dedupe_located=int(len(D)),
            n_pass_window=int(((D.dep >= DEP_MIN) & (D.dep <= DEP_MAX)).sum()),
            n_rejected=int(len(bad)),
            rejected_on_bottom_rail=int((bad.dep == WA_RAILS[1]).sum()),
            rejected_on_floor_rail=int((bad.dep == WA_RAILS[0]).sum()),
            rejected_genuinely_interior_depth=int((~bad.dep.isin(WA_RAILS)).sum()),
            n_on_floor_rail_kept=int((D.dep == WA_RAILS[0]).sum()),
            frac_rejected=float(len(bad) / len(D)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{SF.OUT}/depth_check.json")
    a = ap.parse_args()
    os.makedirs(SF.OUT, exist_ok=True)
    R = dict(north=north(), washington=wa(), pass_frac_target=PASS_FRAC,
             dep_min=DEP_MIN, dep_max=DEP_MAX)
    json.dump(R, open(a.out, "w"), indent=1, default=float)
    print(json.dumps(R, indent=1, default=float), flush=True)


if __name__ == "__main__":
    main()
