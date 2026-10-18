#!/usr/bin/env python
"""Write `results/WA_2008_2007_ARRAYS.md`.  Prefix `y_`.  Robust to missing pieces: every section
prints what is on disk and says so when a stage has not produced its output yet."""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
sys.path.insert(0, HERE)
Y = f"{RES}/y_wa_arrays"
BEDS = ("2008", "2007")


def rd(p, **kw):
    try:
        return pd.read_csv(p, **kw) if os.path.exists(p) and os.path.getsize(p) else None
    except Exception:                                                       # noqa: BLE001
        return None


def rj(p):
    try:
        return json.load(open(p)) if os.path.exists(p) and os.path.getsize(p) else None
    except Exception:                                                       # noqa: BLE001
        return None


def _cell(v, rnd):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ""
    if isinstance(v, (float, np.floating)):
        return f"{round(float(v), rnd):g}"
    return str(v)


def tbl(T, cols=None, n=0, rnd=4):
    """A markdown table without the `tabulate` dependency (it is not installed on this pod)."""
    if T is None or not len(T):
        return "_(not produced)_\n"
    T = T[[c for c in cols if c in T.columns]] if cols else T
    if n:
        T = T.head(n)
    T = T.dropna(axis=1, how="all")
    if not len(T.columns):
        return "_(not produced)_\n"
    head = "| " + " | ".join(str(c) for c in T.columns) + " |"
    rule = "|" + "|".join("---" for _ in T.columns) + "|"
    body = ["| " + " | ".join(_cell(v, rnd) for v in r) + " |"
            for r in T.itertuples(index=False, name=None)]
    return "\n".join([head, rule] + body) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{RES}/WA_2008_2007_ARRAYS.md")
    a = ap.parse_args()
    L = []
    w = L.append
    w("# Washington finished: the Puget / S-Olympic 2008 and 2007 beds run AS ARRAYS\n")
    w("Agent Y, 2026-09-24.  Queue `scripts/y_wa_finish_queue.sh`; bed adapter "
      "`scripts/y_common.py`; markers `results/Y_WA_DONE|FAILED`, interim "
      "`results/Y_WA_MODE_DONE`, heartbeat `results/Y_WA_HEARTBEAT`.\n")
    w("**The gap this closes.**  The one change that made Washington work on the 2010 Olympic "
      "bed was treating each array AS AN ARRAY: agent O measured the tremor/quiet pair-count "
      "ratio at **1.39 member-wise against 15.81 site-averaged** "
      "(`results/WA_2010_ARRAY_RESULTS.md` sec.0d), and the bed built that way produced "
      "**61 of Chestler's 69 published Olympic coordinates with a real-chain family inside "
      "10 km, against 3 of 69 for the whole-chain site-shift null** "
      "(`results/WA_LOCATION_REBUILD.md` sec.9).  The 2008 Puget / S-Olympic bed had never been "
      "run that way -- agent N (`WA_2008_RESULTS.md`) and agent P (`WA_2008Q_RESULTS.md`) both "
      "ran all 51 stations MEMBER-WISE, exactly the configuration that failed on the Olympic "
      "bed, and both got a pair ratio of 0.86-0.89 and zero candidates.  The 2007 window (CAFE "
      "at its fullest) had never been run for discovery at all.  Both are run here, in the mode "
      "the data chooses, with the REBUILT locator (`L9_xcdiff` + `L12_comb`, "
      "`results/WA_LOCATION_REBUILD.md`) instead of the old `k_locate` whose median error "
      "against Bostock's published coordinates was 27.6 km.\n")
    w("**What is being tested, and why a zero is a bug report.**  Royer & Bostock report that "
      "LFE families EXIST in Washington on the plate interface at **29-46 km** depth.  That "
      "claim is falsifiable here WITHOUT their event list: run their data and either we produce "
      "locatable, repeating, interface-depth families above the whole-chain null in that band, "
      "or we have failed to reproduce them.  Their catalogue would only say WHICH families we "
      "recovered, never WHETHER any are there.  So this bed has a known positive answer, and it "
      "is treated the way the north was treated (memory `reproduce-known-answer-first`): the two "
      "previous Puget zeros (`results/WA_2008_RESULTS.md`, `results/WA_2008Q_RESULTS.md`) are "
      "**BUG REPORTS AGAINST OUR CODE, not evidence about Washington**, and if this run again "
      "produces nothing locatable above the null that is written as *our chain still fails on a "
      "bed where families are published to exist*, with every difference from the 2010 Olympic "
      "bed that DID work named and ranked as a candidate cause (sec.0b).  The 2010 bed's numbers "
      "are carried in the SAME tables as the positive control, so the reader sees directly "
      "whether 2008/2007 behave like the bed that worked or like its null.\n")
    w("**Rules.**  RANK, DO NOT DISCARD: every grown family keeps its row, every criterion is a "
      "column, the whole-chain site-shift null gets the identical table, and the evidence is the "
      "real-vs-null separation curve.  The recipe is agent K's, literally and unedited: 6-s "
      "windows, 1-8 Hz, per-component cc summed, 8 x MAD of the block, 4-s decluster, "
      "Bungum-Husebye realignment, transitive grouping >= 3, stack -> template -> 8 x MAD "
      "matched filter, pass 0, literature dedupe only.  `o_common.py` / `o_beam.py` (the array "
      "machinery) and `n_common.py` / `n_hours.py` (the 2008/2007 bed) are REUSED UNCHANGED; "
      "`y_common.py` only rebinds o_common's bed facts onto agent N's bed, the same device "
      "o_common uses on k_common.\n")

    # ------------------------------------------------------------------ 0. THE HEADLINE -------
    w("\n## 0. THE EXISTENCE TEST -- are there locatable interface-depth families, above the "
      "null, in the Royer & Bostock band?\n")
    w("Counts of LOCATED families that fall inside the bed's own sub-region AND inside "
      "**29-46 km**, real beside its whole-chain site-shift null, with bootstrap radii, and with "
      "the 2010 N-Olympic bed carried through the identical code as the positive control.  The "
      "band and the Slab2 interface are REFERENCE LINES drawn on the distribution, never gates "
      "(memory `rank-dont-discard`): the evidence is the real-vs-null curve of the next table, "
      "not any count past a threshold.\n")
    any_ex = False
    for b in BEDS:
        for v in ("L9_xcdiff", "L12_comb"):
            E = rd(f"{Y}/score_{b}_{v}/y_existence_{b}_{v}.csv")
            if E is None:
                continue
            any_ex = True
            w(f"\n**{b} -- {v}**\n\n")
            w(tbl(E, rnd=3))
    if not any_ex:
        w("\n_(the existence test has not been produced yet -- the queue has not reached "
          "stage 5)_\n")
    w("\n**How to read it.**  `n_in_RB_band` for `real` far above the same number for `null`, "
      "at bootstrap radii small enough for the depth to mean anything, is the reproduction.  "
      "`real` indistinguishable from `null` is a failure OF THIS CHAIN on a bed where families "
      "are published to exist -- not a statement about Washington -- and sec.0b is then the list "
      "of candidate causes.  The 2010 rows are the calibration of what success looks like in "
      "this table: on that bed the rebuilt locator put **788 of 1,479** real families on the "
      "map against **230 of 1,471** for the null, and **61 of Chestler's 69** published "
      "coordinates had a real family within 10 km against **3 of 69** for the null "
      "(`results/WA_LOCATION_REBUILD.md` sec.9).\n")
    w("\n**The curve** (`results/y_wa_arrays/score_<bed>_<variant>/y_existence_curve_<bed>_"
      "<variant>.csv`): in-band count for every value of bootstrap radius, cost and station "
      "count, per population.  No threshold is used anywhere; the whole curve is the evidence.\n")

    # ------------------------------------------------------------------ 0b. the differences ---
    w("\n## 0b. Every difference between these beds and the 2010 Olympic bed that WORKED\n")
    w("Ranked as candidate causes, each with whatever measurement this study permits.  This "
      "section is the one that matters if sec.0 comes out flat.\n\n")
    drows = []
    B10 = rj(f"{RES}/o_oly2010/bed_2010.json") or {}
    M10 = rj(f"{RES}/o_oly2010/mode_smoke.json") or {}
    r10 = {t["mode"]: t.get("tremor_quiet_ratio") for t in M10.get("table", [])}
    for b in BEDS:
        B = rj(f"{Y}/bed_{b}.json") or {}
        S = rd(f"{Y}/sites_{b}.csv")
        M = rj(f"{Y}/{b}/mode_smoke.json") or {}
        rr = {t["mode"]: t.get("tremor_quiet_ratio") for t in M.get("table", [])}
        ar = S[S.kind == "ARRAY"] if S is not None else None
        drows.append(dict(
            bed=b, n_member_stations=B.get("n_member_stations"),
            n_sites=B.get("n_sites"), n_array_sites=B.get("n_array_sites"),
            n_members_in_arrays=B.get("n_members_in_arrays"),
            frac_members_in_arrays=(round(B["n_members_in_arrays"] / B["n_member_stations"], 3)
                                    if B.get("n_member_stations") else None),
            median_array_aperture_km=(round(float(ar.aperture_km.median()), 3)
                                      if ar is not None and len(ar) else None),
            max_array_members=(int(ar.n_members.max()) if ar is not None and len(ar) else None),
            n_tremor_hours=B.get("n_tremor_hours"), n_tremor_days=B.get("n_tremor_days"),
            n_quiet_days=B.get("n_quiet_days"),
            hour_source="waveform audit (n_hours.py)",
            n_intraslab_corr=B.get("n_intraslab_corr"),
            n_active_cells=B.get("n_active_cells"),
            ratio_member=rr.get("member"), ratio_beam=rr.get("beam"),
            ratio_siteavg=rr.get("siteavg")))
    S10 = rd(f"{RES}/o_oly2010/sites.csv")
    a10 = S10[S10.kind == "ARRAY"] if S10 is not None else None
    drows.append(dict(
        bed="2010 (worked)", n_member_stations=B10.get("n_member_stations"),
        n_sites=B10.get("n_sites"), n_array_sites=B10.get("n_array_sites"),
        n_members_in_arrays=B10.get("n_members_in_arrays"),
        frac_members_in_arrays=(round(B10["n_members_in_arrays"] / B10["n_member_stations"], 3)
                                if B10.get("n_member_stations") else None),
        median_array_aperture_km=(round(float(a10.aperture_km.median()), 3)
                                  if a10 is not None and len(a10) else None),
        max_array_members=(int(a10.n_members.max()) if a10 is not None and len(a10) else None),
        n_tremor_hours=B10.get("n_tremor_hours"), n_tremor_days=B10.get("n_tremor_days"),
        n_quiet_days=B10.get("n_quiet_days"), hour_source="PNSN tremor catalogue",
        n_intraslab_corr=B10.get("n_intraslab_corr"),
        n_active_cells=B10.get("n_active_cells"),
        ratio_member=r10.get("member"), ratio_beam=r10.get("beam"),
        ratio_siteavg=r10.get("siteavg")))
    w(tbl(pd.DataFrame(drows), rnd=3))
    w("\nThe band (1-8 Hz), the window length (6 s), the threshold (8 x MAD of the block), the "
      "declustering (4 s), the grouping bar (>= 3) and the matched-filter bar are IDENTICAL "
      "across all three beds by construction -- they are agent K's frozen constants and no "
      "stage here can move them -- so they cannot be a difference.  The depth-grid geometry is "
      "also identical (`uu_locate.WA_BASE_GRID` / `WA_WIDE_GRID` covers 46.3-49.0 N, 6-76 km, "
      "for every Washington population including the 2010 control).  What remains, and is "
      "measured above, is: array fraction and aperture, station count and quality, the "
      "tremor-hour source, the number of hours, and the station-term coverage.\n")

    # ------------------------------------------------------------------ 1. site structure -----
    w("\n## 1. Site structure -- stations to sites, both beds\n")
    w("Single-link clustering at **1.5 km**, partitioned by NETWORK (a borehole and a surface "
      "sensor 100 m apart are not one array); every (cluster, network) sub-group with >= 2 "
      "members is ONE site, singletons stay single stations.  Identical code and identical "
      "distance to agent O's 2010 partition (`o_common._build_sites`).\n")
    for b in BEDS:
        S = rd(f"{Y}/sites_{b}.csv")
        B = rj(f"{Y}/bed_{b}.json")
        w(f"\n### {b}\n")
        if B:
            w(f"* **{B.get('n_member_stations')} member stations -> {B.get('n_sites')} sites** "
              f"({B.get('n_array_sites')} arrays covering {B.get('n_members_in_arrays')} "
              f"members + {B.get('n_single_sites')} single stations); window "
              f"{B.get('day0')}..{B.get('day1')} (tranche {B.get('tranche')}), sub-region "
              f"{B.get('subregion')}.\n")
            w(f"* audited tremor hours **{B.get('n_tremor_hours')}** on "
              f"{B.get('n_tremor_days')} days; literature-quiet days "
              f"{B.get('n_quiet_days')}; active cells {B.get('n_active_cells')}.\n")
            w(f"* members by network {B.get('members_by_net')}; arrays by network "
              f"{B.get('arrays_by_net')}.\n")
            w(f"* intraslab S station terms on {B.get('n_intraslab_corr')} sites; "
              f"{B.get('n_corr_flagged')} carry 0.0 and are FLAGGED, never imputed.\n")
        if S is not None:
            w("\n" + tbl(S[S.kind == "ARRAY"],
                         ["site", "net", "n_members", "aperture_km", "lat", "lon", "members"],
                         rnd=3))
            w("\nSingle stations: " + ", ".join(sorted(S[S.kind == "SINGLE"].site)) + "\n")
        else:
            w("_(sites not yet built)_\n")

    # ------------------------------------------------------------------ 2. mode ----------------
    w("\n## 2. The site treatment, chosen BY MEASUREMENT on each bed\n")
    w("The identical literal autocorrelation run three ways -- `member` (every station one term "
      "of the network mean, what agents N and P did), `beam` (one zero-slowness delay-and-sum "
      "trace per array, what Ghosh / Chestler / Sweet did), `siteavg` (one vote per site, member "
      "information kept) -- on the same tremor and quiet hours, restricted to the member set "
      "present in both populations.  The winner is the highest tremor/quiet pair-count ratio.  "
      "2010 N-Olympic for scale: member **1.385**, beam **0.991**, siteavg **15.81**.\n")
    for b in BEDS:
        M = rj(f"{Y}/{b}/mode_smoke.json") or rj(f"{Y}/mode_{b}/mode_smoke.json")
        w(f"\n### {b}\n")
        if M:
            w(tbl(pd.DataFrame(M["table"]),
                  ["mode", "n_blocks_tremor", "det_per_hour_tremor", "traces_tremor",
                   "thr_sd_tremor", "n_blocks_quiet", "det_per_hour_quiet", "traces_quiet",
                   "thr_sd_quiet", "tremor_quiet_ratio"], rnd=4))
            w(f"\n**Chosen: `{M['chosen']}`** -- {M['why']}\n")
        else:
            w("_(mode comparison not yet produced)_\n")

    # ------------------------------------------------------------------ 3. pair + grouping -----
    w("\n## 3. Pair stage, grouping and growth -- real beside the site-shift null\n")
    rows = []
    for b in BEDS:
        for tag, sub in (("real", b), ("null", f"{b}_null")):
            P = rj(f"{Y}/{sub}/pair_pairstats.json")
            g = rd(f"{Y}/{sub}/gb_summary.csv")
            d = f"{Y}/{sub}/grow/det_it0.parquet"
            r = dict(bed=b, chain=tag)
            if P:
                for k in ("n_blocks_tremor", "n_det_tremor", "det_per_hour_tremor",
                          "n_blocks_quiet", "n_det_quiet", "det_per_hour_quiet",
                          "tremor_quiet_ratio", "nmad_in_sd_median", "sta_per_block_median"):
                    if k in P:
                        r[k] = P[k]
            if g is not None:
                r["n_groups"] = len(g)
            if os.path.exists(d):
                D = pd.read_parquet(d, columns=["grp"])
                r["n_templates"] = int(D.grp.nunique())
                r["n_detections"] = int(len(D))
            rows.append(r)
    w(tbl(pd.DataFrame(rows), rnd=3))

    # ------------------------------------------------------------------ 4. screen --------------
    w("\n## 4. The five-clause screen -- identical table for real and for the null\n")
    w("The screen is agent N's, unedited (`n_screen.py`).  It is reported as a LEDGER, not as a "
      "gate: the clause counts are columns, every family keeps its row in the ranked tables of "
      "sec.5, and the evidence is real-minus-null.  Its clause (c) needs a coordinate, so the "
      "OLD `n_locate` is run once to feed this ledger and NOTHING ELSE -- every location "
      "reported in sec.0 and sec.5, and every distance, depth and bootstrap radius, comes from "
      "the REBUILT locator.\n")
    for b in BEDS:
        for tag, sub in (("real", b), ("null", f"{b}_null")):
            S = rd(f"{Y}/{sub}/screen.csv")
            w(f"\n**{b} {tag}** -- ")
            if S is None:
                w("_(not produced)_\n")
                continue
            cl = [c for c in S.columns if c.startswith("c_")]
            w(f"{len(S)} grown templates; clause counts (alone): "
              + ", ".join(f"`{c}` {int(S[c].astype(bool).sum())}" for c in cl) + "\n")
            if "pass_abcd" in S.columns:
                w(f"\npass (a)-(d): **{int(S.pass_abcd.astype(bool).sum())}**\n")

    # ------------------------------------------------------------------ 5. location ------------
    w("\n## 5. Location with the REBUILT locator, and the score\n")
    w("`L9_xcdiff` (adopted) and `L12_comb` (companion), `scripts/uu_locate.py`, with a "
      "40-draw station bootstrap radius and a 20-draw shuffled-station control as COLUMNS -- not "
      "the old `k_locate`, whose median epicentral error against Bostock's published "
      "coordinates is 27.6 km (`results/WA_LOCATION_REBUILD.md` sec.1).\n")
    for b in BEDS:
        for v in ("L9_xcdiff", "L12_comb"):
            d = f"{Y}/score_{b}_{v}"
            CH = rd(f"{d}/y_chestler_{b}_{v}.csv")
            DS = rd(f"{d}/y_depth_strike_{b}_{v}.csv")
            w(f"\n### {b} -- {v}\n")
            if CH is None and DS is None:
                w("_(not produced)_\n")
                continue
            w("\n**(b) distance to Chestler & Creager's 69 published Olympic coordinates** "
              "(`lat_min = 47.3` is the row restricted to where these beds overlap the Olympics; "
              "`lat_min = -90` is every located family).  Royer & Bostock released no "
              "Washington catalogue, so this is the only published coordinate set these beds can "
              "be scored against, and it only covers the northern edge of the Puget bed.\n\n")
            w(tbl(CH, rnd=3))
            w("\n (the band is a REFERENCE LINE, never a cut -- the fraction inside "
              "it is one more column, reported for the null too).\n\n")
            w(tbl(DS, rnd=3))
            SEP = rd(f"{d}/y_separation_{b}_{v}.csv")
            if SEP is not None and len(SEP):
                w("\n**(a) real-vs-null separation, best point of each criterion's curve** "
                  "(the full curves are in `" + os.path.relpath(f"{d}/y_separation_{b}_{v}.csv",
                                                                STUDY) + "`; no threshold is "
                  "used anywhere -- this row is only the reader's entry point into the curve).\n\n")
                best = (SEP[SEP.real >= 5].sort_values("excess", ascending=False)
                        .groupby("criterion", as_index=False).head(1)
                        .sort_values("excess", ascending=False))
                w(tbl(best, ["criterion", "sense", "threshold", "real", "null", "excess",
                             "ratio"], n=20, rnd=3))

    # ------------------------------------------------------------------ 6. files ---------------
    w("\n## 6. Files\n")
    w("* bed adapter `scripts/y_common.py`; dispatcher `scripts/y_run.py`; mode selection "
      "`scripts/y_mode.py`; scoring `scripts/y_score.py`; report `scripts/y_report.py`; queue "
      "`scripts/y_wa_finish_queue.sh`.\n")
    w("* REUSED UNCHANGED: `scripts/o_common.py`, `o_beam.py`, `o_auto.py` (agent O's array "
      "machinery), `n_common.py`, `n_hours.py`, `n_chunk.py`, `n_screen.py`, `n_locate.py`, "
      "`n_relocate.py`, `n_refire.py` (agent N's bed and stages), `k_auto.py`, `k_group.py`, "
      "`k_grow.py` (agent K's frozen recipe), `uu_*.py` (agent U's rebuilt locator).\n")
    for b in BEDS:
        w(f"* bed {b}: sites `results/y_wa_arrays/sites_{b}.csv`, chain "
          f"`results/y_wa_arrays/{b}/`, null `results/y_wa_arrays/{b}_null/`, scores "
          f"`results/y_wa_arrays/score_{b}_*/`.\n")
    w("* beams `beamcache_y/<SITE>/<year>/<doy>.npz` (zero-slowness delay-and-sum; the slowness "
      "search is not identifiable at 1 km aperture, `WA_2010_ARRAY_RESULTS.md` sec.0c).\n")
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w").write("".join(L))
    print(f"wrote {a.out} ({sum(len(x) for x in L)} bytes)", flush=True)


if __name__ == "__main__":
    main()
