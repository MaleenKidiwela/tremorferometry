#!/usr/bin/env python
"""Assemble results/WA_LOCATION_REBUILD.md from whatever stages have finished.  Prefix `uu_`."""
from __future__ import annotations

import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import uu_common as U                                                       # noqa: E402

O = U.OUT


def jl(p, d=None):
    try:
        return json.load(open(p))
    except Exception:
        return d


def tl(p):
    try:
        return open(p).read().strip()
    except Exception:
        return ""


def md(df, cols=None, r=2):
    if df is None or not len(df):
        return "_(not produced)_\n"
    d = df[cols] if cols else df
    d = d.copy()
    for c in d.columns:
        if d[c].dtype.kind == "f":
            d[c] = d[c].round(r)
    head = "| " + " | ".join(str(c) for c in d.columns) + " |"
    sep = "|" + "|".join(["---"] * len(d.columns)) + "|"
    body = "\n".join("| " + " | ".join("" if pd.isna(v) else str(v) for v in row) + " |"
                     for row in d.itertuples(index=False))
    return f"{head}\n{sep}\n{body}\n"


def main():
    out = []
    A = out.append
    A("# Rebuilding the LFE family LOCATION step, calibrated on Bostock's published coordinates\n")
    A("Agent U, 2026-09-23.  Workspace `results/u_locate/`; queue `scripts/u_locate_queue.sh` "
      "(markers `results/U_LOCATE_DONE|FAILED`, interim `results/U_LOCATE_BASELINE_DONE`); code "
      "`scripts/uu_common.py`, `uu_stacks.py`, `uu_picks.py`, `uu_locate.py`, `uu_dd.py`, "
      "`uu_eval.py`, `uu_wa.py`.  **Detection and growth are untouched** -- every template, every "
      "detection and every stack member comes from the frozen chain on disk; this study only "
      "re-reads them and replaces what happens after.\n")
    A("No family is ever filtered out.  Every grown family keeps its row in every table and every "
      "criterion is a column (memory `rank-dont-discard`).\n")

    ts = jl(f"{O}/north_eval/truth_source.json")
    A("\n## 0. The calibration set, and the proof that it is the published answer\n")
    if ts:
        A(f"`lfe_catalog_repair/results/e5r_families_2012.csv` holds **{ts['n_families']} family "
          f"coordinates**, and **{ts['n_ours_in_catalog']} of them** are exactly among the "
          f"{ts['n_catalog_triples']} distinct (lat, lon, dep) triples carried by the event rows "
          f"of Bostock et al. (2015)'s own catalogue export "
          f"(`results/external/bostock/bostock2015_sloweq_FULL_2003-2013.csv`).  That file is "
          f"the published location of each family; `bkey` is this project's index into it.  "
          f"Depth convention: km positive down, {ts['dep_min']:.1f}-{ts['dep_max']:.1f} km "
          f"(median {ts['dep_median']:.1f}).\n")
    A("The 58 families the frozen recipe rediscovered blind on the 2005 POLARIS bed "
      "(`results/repro/classify15.csv`, column `kept_bkey`) are the DEVELOPMENT set.  The 18 "
      "located SEALED families (`hid_bkey`) are scored ONCE at the end and no choice is made on "
      "them.  One representative template per family is fixed before any variant runs (highest "
      "one-to-one match fraction), so every variant is scored on the same 58 objects.\n")
    A("The noise population is the 1,219 grown templates that match no catalogued family; it gets "
      "the identical table.\n")

    A("\n## 1. BASELINE -- what the current locator does, measured for the first time\n")
    L = None
    try:
        L = pd.read_csv(f"{O}/north_eval/ladder.csv")
    except Exception:
        pass
    if L is not None:
        base = L[L.variant.str.startswith("BASELINE")]
        A(md(base, ["variant", "n_families", "n_located", "epi_med", "epi_p90", "dep_med",
                    "dep_p90", "dep_bias", "d3_med", "rail", "within5", "within10", "within20",
                    "null_d3_med"]))
        A("\n`epi_*` = epicentral distance from our solution to Bostock's published coordinate "
          "for THAT family, in km; `dep_*` = |depth error|; `dep_bias` = median signed depth "
          "error; `rail` = fraction of located solutions sitting on a face of the search grid; "
          "`null_d3_med` = the same solutions scored against a RANDOMLY PERMUTED family label "
          "(the honest chance level).\n")

    pr = jl(f"{O}/north/pick_residuals.json")
    A("\n## 2. Where the error comes from -- per-station pick residuals against the published "
      "locations\n")
    if pr:
        A("Residual = observed pick minus the time predicted by Bostock's published coordinate "
          "for that family in the Crosson 1-D model, median-centred per family (the origin time "
          "is free).\n")
        A(md(pd.DataFrame(pr["by_threshold"])))
        A("\n")
        A(md(pd.DataFrame(pr["by_station"]).head(30)))
        A(f"\n**Reading.** {pr['headline']}\n")

    A("\n## 3. The ladder -- every change, measured on the same 58 families\n")
    if L is not None:
        lad = L[L.variant.str.startswith("L")]
        A(md(lad, ["variant", "n_located", "epi_med", "epi_p90", "dep_med", "dep_bias", "d3_med",
                   "rail", "within5", "within10", "within20", "noise_located", "noise_rail"]))

    A("\n## 4. Ablations -- what was tried and DISCARDED because it did not improve the error\n")
    if L is not None:
        ab = L[L.variant.str.startswith("A_")]
        A(md(ab, ["variant", "n_located", "epi_med", "epi_p90", "dep_med", "dep_bias", "d3_med",
                  "rail", "within10"]))

    ad = jl(f"{O}/north_eval/adopted.json")
    A("\n## 5. The adopted locator\n")
    if ad:
        A(f"**Adopted: `{ad['variant']}`** -- {ad.get('spec')}\n")
        A(md(pd.DataFrame(ad["rows"])))
        A(f"\nAcceptance (stated before the run): adopt only if the median distance to Bostock's "
          f"published coordinates improves materially AND the rail fraction drops, on the 58 "
          f"known families.  Measured: epicentral median "
          f"{ad['baseline']['epi_med']:.1f} -> {ad['best']['epi_med']:.1f} km, 3-D median "
          f"{ad['baseline']['d3_med']:.1f} -> {ad['best']['d3_med']:.1f} km, rail "
          f"{ad['baseline']['rail']:.2f} -> {ad['best']['rail']:.2f}, families located "
          f"{ad['baseline']['n_located']} -> {ad['best']['n_located']} of 58.  "
          f"Verdict: **{ad['verdict']}**.\n")

    A("\n## 6. Blind confirmation -- the 18 located SEALED families\n")
    if L is not None:
        bl = L[L.variant.str.contains("BLIND")]
        A(md(bl, ["variant", "n_families", "n_located", "epi_med", "epi_p90", "dep_med",
                  "d3_med", "within10"]))

    A("\n## 7. Per-family uncertainty\n")
    try:
        B = pd.read_csv(f"{O}/north_eval/per_family.csv")
        ad_name = (ad or {}).get("variant", "")
        b = B[B.variant == ad_name]
        if len(b) and "boot_km" in b.columns:
            A(f"Station bootstrap (resample the picking stations with replacement, 60 "
              f"realisations, 68th percentile of the epicentral scatter): median "
              f"{np.nanmedian(b.boot_km):.1f} km, p90 {np.nanquantile(b.boot_km, .9):.1f} km.  "
              f"Correlation with the TRUE error: "
              f"{np.corrcoef(b.dropna(subset=['boot_km', 'err_epi']).boot_km, b.dropna(subset=['boot_km', 'err_epi']).err_epi)[0, 1]:.2f} "
              f"-- this is the number the Washington catalogue carries per family.\n")
            A(md(b[["bkey", "grp", "n_sta", "lat", "lon", "dep", "boot_km", "boot_dz",
                    "err_epi", "err_dep"]].sort_values("err_epi").head(60)))
    except Exception as e:
        A(f"_(bootstrap table unavailable: {e})_\n")

    A("\n## 8. Double-difference pass\n")
    if L is not None and (L.variant == "ZZdd").any():
        A("`ZZdd` is the adopted solution after a Waldhauser-Ellsworth double-difference pass "
          "between neighbouring families, with the differential times measured by cross-"
          "correlating the two families' long stacks at the SAME station "
          "(`scripts/uu_dd.py`, observations in `results/u_locate/north/dd_ZZdd_obs.csv`).\n")
        ad_name = (ad or {}).get("variant", "")
        A(md(L[L.variant.isin(["ZZdd", ad_name])],
             ["variant", "n_located", "epi_med", "epi_p90", "dep_med", "d3_med", "rail",
              "within5", "within10", "within20"]))
        A("\nAcceptance, stated before the pass ran: adopt the double-difference solution only "
          "if it improves that table; otherwise it is reported and not used.\n")
        try:
            z = L[L.variant == "ZZdd"].iloc[0]
            b = L[L.variant == ad_name].iloc[0]
            better = (z.d3_med < b.d3_med) and (z.epi_med <= b.epi_med)
            A(f"\n**Measured: {'ADOPTED' if better else 'NOT ADOPTED'}.**  Median 3-D distance "
              f"{b.d3_med:.2f} km (`{ad_name}`) -> {z.d3_med:.2f} km (`ZZdd`), epicentral "
              f"{b.epi_med:.2f} -> {z.epi_med:.2f} km, within 10 km "
              f"{int(b.within10)} -> {int(z.within10)} of 58.  "
              + ("The pass improves the known answer and is carried forward.\n" if better else
                 "The pass does NOT improve the known answer, so the single-family solution "
                 "stands and Washington is located with it.  The reason is visible in "
                 "`logs/u_north_dd.log`: only 1,123 differential times over 179 linked pairs "
                 "survive the cc >= 0.4 / >= 4 shared-station requirement out of 3,363 "
                 "candidate pairs, so most families receive no differential constraint at all "
                 "(median |move| 0.00 km) and the few that move do not move toward the "
                 "published coordinates.\n"))
        except Exception:
            pass
    else:
        A("_(the double-difference stage produced no output; see `logs/u_north_dd.log`)_\n")

    A("\n## 9. Washington -- relocated, with the identical table for the null and for Arm B\n")
    try:
        S = pd.read_csv(f"{O}/wa/wa_chestler_summary.csv")
        A(md(S))
        A("\n`real` = the 2010 Olympic arrays-as-arrays discovery run "
          "(`results/o_oly2010/2010/`); `null` = its whole-chain site-shift null "
          "(`2010_null/`); `armB` = Chestler & Creager's OWN detection chain run literally on her "
          "own data (`chestler_arm/`, 739 families).  `chest_km_3d` = hypocentral distance to the "
          "nearest of the 69 published Olympic coordinates.  Ranked per-family tables with every "
          "criterion as a column: `results/u_locate/wa/wa_{real,null,armB}_ranked.csv`; "
          "real-vs-null separation curves in `wa_separation_curves.csv`.\n")
    except Exception as e:
        A(f"_(Washington stage not finished: {e})_\n")

    A("\n## 9b. THE PRIORITY QUESTION -- the 280 richest real 2010-Olympic templates, relocated\n")
    A("`results/sf_final/table_oly2010_{real,null}.csv` report **280 real 2010-Olympic templates "
      "with >= 20 dedupe partners and ZERO null templates that rich** (the null maximum is 16); "
      "AUC 0.758 in the correct direction, surviving a detection-richness control at matched "
      "`n_det`.  That is the strongest candidate family population Washington has produced, so "
      "it is the population the rebuilt locator has to be pointed at.  Below they are relocated "
      "with the adopted locator and carry their own bootstrap radius, and they are placed beside "
      "the **null's best-ranked rows** -- the equal-size head of the null ordered by the SAME "
      "criterion, which is the matched comparator.  This filters nothing: all 1,485 real "
      "families keep their rows in `wa_real_ranked.csv` and `dedupe_n_partners` / `sf_top` are "
      "simply two more columns (memory `rank-dont-discard`).\n")
    got = False
    for tag, d in (("adopted", "wa"), ("coverage companion", "wa_cov")):
        try:
            SS = pd.read_csv(f"{O}/{d}/wa_sf_top_summary.csv")
        except Exception:
            continue
        got = True
        A(f"\n**{tag} locator** (`{O}/{d}/`)\n")
        cols = [c for c in ("pop", "n_families", "n_located", "rail", "boot_km_med",
                            "boot_km_p90", "med_nearest_3d", "med_nearest_epi",
                            "p10_nearest_3d", "fam_within_5km", "fam_within_10km",
                            "fam_within_20km", "coords_hit_5km", "coords_hit_10km",
                            "coords_hit_20km", "n_coords") if c in SS.columns]
        A(md(SS, cols))
    if got:
        A("\n`pop` rows: `real[dedupe_n_partners>=20]` = the 280; "
          "`null[top 280 by dedupe_n_partners]` = the matched head of the null; "
          "`null[dedupe_n_partners>=20]` = empty by construction (no null template is that "
          "rich), kept as a row so the table is the identical table; `real[all]` / `null[all]` = "
          "the whole populations for scale.  `med_nearest_3d` = median hypocentral distance from "
          "a located family to the NEAREST of Chestler & Creager's 69 published Olympic "
          "coordinates; `coords_hit_k` = how many of those 69 have at least one family within "
          "k km, the reverse direction that a mere pile-up cannot fake.  `boot_km` = the "
          "per-family station-bootstrap epicentral radius, so a distance is only meaningful "
          "beside it.\n")
        A("\nReading the three outcomes this distinguishes: if the 280 sit on Chestler's "
          "coordinates within their bootstrap radii and cover many of the 69, they ARE her "
          "families rediscovered; if they sit in the same depth/along-strike band but away from "
          "her coordinates by more than their radii, they are new families on the same "
          "interface; if their distance distribution is indistinguishable from the null's "
          "best-ranked rows, the dedupe-partner excess is not locating anything.  Per-family "
          f"rows: `{O}/wa/wa_sf_top_real_located.csv`.\n")
    else:
        A("_(the Washington relocation of the 280 has not finished; see "
          "`logs/u_wa_chest_*.log`)_\n")

    A("\n## 10. Files\n")
    A("| what | path |\n|---|---|\n"
      f"| long stacks (44 s, 12 s pre-node) | `{O}/north_stacks.npz`, `{O}/wa_*_stacks.npz` |\n"
      f"| picks | `{O}/north/picks.parquet`, `{O}/wa_*/picks.parquet` |\n"
      f"| pick residuals vs published locations | `{O}/north/pick_residuals.csv` |\n"
      f"| every variant's locations | `{O}/north/loc_*.csv` |\n"
      f"| ladder + per-family errors | `{O}/north_eval/ladder.csv`, `per_family.csv` |\n"
      f"| adopted locator | `{O}/north_eval/adopted.json` |\n"
      f"| Washington ranked tables | `{O}/wa/wa_*_ranked.csv` |\n"
      f"| the 280 richest real templates, relocated | `{O}/wa/wa_sf_top_summary.csv`, "
      f"`wa_sf_top_real_located.csv` |\n"
      f"| double-difference observations | `{O}/north/dd_ZZdd_obs.csv` |\n")
    p = f"{U.RES}/WA_LOCATION_REBUILD.md"
    open(p, "w").write("\n".join(out))
    U.log(f"wrote {p}")


if __name__ == "__main__":
    main()
