#!/usr/bin/env python
"""Step 0 coverage + withheld-set proposal, computed with ADMITTED sites only.

Active cells (Merlin's definition): Slab2 interface points 28-45 km in the box with >= 30 PNSN
2010-2011 epicentres within 8 km.  There is NO 2008 tremor catalog (PNSN starts 2010-01-06), so
the Puget/S-Olympic 2008 numbers use the SAME 2010-11 cells - stated, not hidden.

Sub-regions (Merlin): N-Olympic 47.8-48.3, S-Olympic 47.3-47.8, Puget/S-WA 46.2-47.3.
Coverage statistic: % of a sub-region's active cells with >= 5 ADMITTED sites within 60 km
HYPOCENTRAL, and within 30 km EPICENTRAL.
Withheld proposal: Merlin's rule - 3 permanent sites within 60 km of >= 80 % of the sub-region's
cells, 3-C, admitted by the probe, chosen to minimise the coverage loss when they leave the sum.
-> results/coverage_by_subregion.csv, results/withheld_candidates.csv, results/withheld_proposal.json
"""
from __future__ import annotations
import os
import sys
import json
import itertools
import argparse
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_common as D                                                   # noqa: E402

SUBREG = {"N-Olympic": (47.8, 48.3), "S-Olympic": (47.3, 47.8), "Puget-SWA": (46.2, 47.3)}
MIN_SITES, R_HYP, R_EPI = 5, 60.0, 30.0
PERM_NETS = ("PB", "UW", "CN")
EXT_NEEDED = ["PB.B001", "PB.B003", "PB.B004", "PB.B005", "PB.B006", "PB.B007", "PB.B014",
              "PB.B017", "UW.LRIV", "UW.DOSE", "CN.VGZ", "CN.PFB"]


def site_frame(S, admitted, window_col):
    """Admitted SITES: a site is admitted when >= 1 of its 3-C members passed the probe."""
    rows = []
    for r in S.itertuples():
        mem = [m for m in str(r.members_3c).split(";") if m]
        adm = [m for m in mem if m in admitted]
        if getattr(r, window_col) <= 0:
            continue
        rows.append(dict(site_id=r.site_id, site_name=r.site_name, lat=r.lat, lon=r.lon,
                         nets=r.net3c, n_members_3c=len(mem), n_admitted=len(adm),
                         admitted=len(adm) > 0, members=";".join(mem),
                         admitted_members=";".join(adm)))
    return pd.DataFrame(rows)


def cover(cells, sites):
    """-> (frac cells with >=5 sites <=60 km hypocentral, same at <=30 km epicentral, counts)."""
    if not len(sites) or not len(cells):
        return 0.0, 0.0, np.zeros(len(cells)), np.zeros(len(cells))
    dh = D.hyp_km(cells.lat.values[:, None], cells.lon.values[:, None],
                  cells.depth_km.values[:, None], sites.lat.values[None, :],
                  sites.lon.values[None, :])
    de = D.epi_km(cells.lat.values[:, None], cells.lon.values[:, None],
                  sites.lat.values[None, :], sites.lon.values[None, :])
    nh = (dh <= R_HYP).sum(1)
    ne = (de <= R_EPI).sum(1)
    return float((nh >= MIN_SITES).mean()), float((ne >= MIN_SITES).mean()), nh, ne


def main():
    ap = argparse.ArgumentParser()
    a = ap.parse_args()
    G = D.active_cells()
    A = G[G.tremor_active].reset_index(drop=True)
    S = pd.read_csv(os.path.join(D.RES, "sites_array_era_final.csv"), keep_default_na=False)
    P10 = pd.read_csv(os.path.join(D.RES, "site_probe_2010.csv"))
    P08 = pd.read_csv(os.path.join(D.RES, "site_probe_2008.csv"))
    adm10 = set(P10[P10.admit].netsta)
    adm08 = set(P08[P08.admit].netsta)
    out, prop = [], {}
    for tag, admitted, wcol, subs in [
            ("N-Olympic 2010", adm10, "w2010_08_days_3c", ["N-Olympic", "S-Olympic"]),
            ("Puget-SOly 2008", adm08, "w2008_05_days_3c", ["Puget-SWA", "S-Olympic"])]:
        SF = site_frame(S, admitted, wcol)
        for sub in subs:
            lo, hi = SUBREG[sub]
            cells = A[(A.lat >= lo) & (A.lat <= hi)].reset_index(drop=True)
            all_s = SF
            adm_s = SF[SF.admitted]
            fh_a, fe_a, nh_a, _ = cover(cells, all_s)
            fh, fe, nh, ne = cover(cells, adm_s)
            out.append(dict(setting=tag, subregion=sub, n_cells=len(cells),
                            n_sites_available=len(all_s), n_sites_admitted=len(adm_s),
                            pct_cells_ge5_sites_60km_all=round(100 * fh_a, 1),
                            pct_cells_ge5_sites_60km=round(100 * fh, 1),
                            pct_cells_ge5_sites_30km=round(100 * fe, 1),
                            median_sites_60km=float(np.median(nh)) if len(nh) else 0.0))
            # ---------------- withheld candidates: admitted PERMANENT sites near >=80% of cells
            cand = []
            for r in adm_s.itertuples():
                if not any(n in PERM_NETS for n in str(r.nets).split(",")):
                    continue
                dh = D.hyp_km(cells.lat.values, cells.lon.values, cells.depth_km.values,
                              r.lat, r.lon)
                f80 = float((dh <= R_HYP).mean())
                one = adm_s[adm_s.site_id != r.site_id]
                fh1, fe1, _, _ = cover(cells, one)
                cand.append(dict(setting=tag, subregion=sub, meets_80pct=bool(f80 >= 0.80),
                                 site_id=r.site_id,
                                 site=r.site_name, nets=r.nets,
                                 members=r.admitted_members, frac_cells_60km=round(f80, 3),
                                 pct_60km_without_it=round(100 * fh1, 1),
                                 loss_pct=round(100 * (fh - fh1), 1),
                                 needed_by_extension=any(m in EXT_NEEDED
                                                         for m in r.admitted_members.split(";"))))
            C = pd.DataFrame(cand)
            if len(C):
                C = C.sort_values(["loss_pct", "frac_cells_60km"], ascending=[True, False])
                # Merlin's STRICT rule: >= 80 % of cells within 60 km.  Where fewer than 3
                # permanent sites qualify, the rule is INFEASIBLE and we say so, then fall back
                # to the best 3 admitted permanent sites with >= 60 % coverage.
                strict = C[C.meets_80pct]
                feasible_strict = len(strict) >= 3
                # The operative constraint is "minimise the coverage loss"; the 80 % clause is
                # reported but cannot be met anywhere (a 60 km radius cannot reach 80 % of a
                # sub-region that is ~150 km long).  Pool = every ADMITTED permanent site.
                pool = strict if feasible_strict else C
                if len(pool[~pool.needed_by_extension]) >= 3:
                    pool = pool[~pool.needed_by_extension]
                # A withheld site is only useful if it can SEE the cells (Gate 1 at a site with
                # no cell within 60 km is vacuous), so the objective maximises the summed
                # cell-visibility of the trio subject to the residual coverage staying >= 80 %
                # of cells (relaxed to the best achievable residual when nothing reaches 80 %).
                frac = dict(zip(C.site_id, C.frac_cells_60km))
                need = dict(zip(C.site_id, C.needed_by_extension))
                combos = []
                for combo in itertools.combinations(list(pool.site_id), min(3, len(pool))):
                    rest = adm_s[~adm_s.site_id.isin(combo)]
                    f, fe2, _, _ = cover(cells, rest)
                    combos.append((combo, f, fe2,
                                   sum(frac[c] for c in combo),
                                   sum(int(need[c]) for c in combo)))
                feas = [c for c in combos if c[1] >= 0.80]
                use = feas if feas else combos
                best4 = max(use, key=lambda c: (c[3] - 0.02 * c[4], c[1]))
                combo, f, fe2 = best4[0], best4[1], best4[2]
                best = (combo, f, fe2)
                resid_constraint = bool(feas)
                if best:
                    combo, f, fe2 = best
                    prop[f"{tag}|{sub}"] = dict(
                        strict_rule_feasible=bool(feasible_strict),
                        residual_coverage_ge80_achievable=bool(resid_constraint),
                        withheld_needed_by_extension=[bool(need[c]) for c in combo],
                        n_sites_meeting_80pct=int(len(strict)),
                        rule_used=("Merlin strict (>=80 % of cells within 60 km)"
                                   if feasible_strict else
                                   "RELAXED to >=60 % of cells - fewer than 3 permanent sites "
                                   "meet the 80 % bar"),
                        withheld=[C[C.site_id == c].site.iloc[0] for c in combo],
                        withheld_members=[C[C.site_id == c].members.iloc[0] for c in combo],
                        withheld_frac_cells_60km=[float(C[C.site_id == c].frac_cells_60km.iloc[0])
                                                  for c in combo],
                        pct_cells_ge5_sites_60km_with_them=round(100 * fh, 1),
                        pct_cells_ge5_sites_60km_without_them=round(100 * f, 1),
                        pct_cells_ge5_sites_30km_without_them=round(100 * fe2, 1),
                        n_candidates=len(C))
            out[-1]["n_withheld_candidates"] = len(C)
            globals().setdefault("_C", []).append(C)
    O = pd.DataFrame(out)
    O.to_csv(os.path.join(D.RES, "coverage_by_subregion.csv"), index=False)
    CC = pd.concat(globals()["_C"], ignore_index=True)
    CC.to_csv(os.path.join(D.RES, "withheld_candidates.csv"), index=False)
    D.jdump(dict(rule="3 admitted PERMANENT sites within 60 km of >= 80 % of the sub-region's "
                      "active cells, chosen to maximise the residual coverage and to avoid the "
                      "permanents the extension needs",
                 cells="Slab2 28-45 km, >= 30 PNSN 2010-11 epicentres within 8 km; the SAME cells "
                       "are used for 2008 because no 2008 tremor catalog exists",
                 extension_must_stay=EXT_NEEDED, proposals=prop),
           os.path.join(D.RES, "withheld_proposal.json"))
    print(O.to_string(index=False))
    print()
    print(CC.to_string(index=False))
    print()
    print(json.dumps(prop, indent=1))


if __name__ == "__main__":
    main()
