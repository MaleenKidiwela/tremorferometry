#!/usr/bin/env python
"""Write results/WA_2008Q_RESULTS.md — the quiet-subset Puget/S-Olympic 2008 run.  Prefix `p_`."""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")

P = os.environ.get("P_OUTROOT", f"{RES}/p_2008")
NROOT = os.environ.get("N_OUTROOT", f"{RES}/n_puget2008")


def jload(p, d=None):
    try:
        return json.load(open(p))
    except Exception:                                                       # noqa: BLE001
        return {} if d is None else d


def cload(p):
    try:
        return pd.read_csv(p)
    except Exception:                                                       # noqa: BLE001
        return pd.DataFrame()


def md(df, cols=None, ff=2):
    if df is None or not len(df):
        return "_(empty)_\n"
    df = df[cols] if cols else df
    head = "| " + " | ".join(str(c) for c in df.columns) + " |"
    sep = "|" + "|".join("---" for _ in df.columns) + "|"
    rows = []
    for r in df.itertuples(index=False):
        cells = []
        for v in r:
            if isinstance(v, (float, np.floating)):
                cells.append(f"{v:.{ff}f}" if np.isfinite(v) else "--")
            elif v is None:
                cells.append("--")
            else:
                cells.append(str(v))
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join([head, sep] + rows) + "\n"


def g(d, k, dflt=np.nan):
    v = d.get(k, dflt)
    return v if v is not None else dflt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{RES}/WA_2008Q_RESULTS.md")
    a = ap.parse_args()
    R, NU, C7, BH = f"{P}/2008", f"{P}/2008_null", f"{P}/2007", f"{P}/bh"

    CEN = cload(f"{P}/station_noise_census.csv")
    BYN = cload(f"{P}/census_by_network.csv")
    SM = jload(f"{P}/subset_meta.json")
    bed = jload(f"{P}/bed_2008q.json")
    bed7 = jload(f"{P}/bed_2007q.json")
    ps51 = jload(f"{NROOT}/2008/pair_pairstats.json")
    ps = jload(f"{R}/pair_pairstats.json")
    psn = jload(f"{NU}/pair_pairstats.json")
    psbh = jload(f"{BH}/pair_pairstats.json")
    sr, sn = jload(f"{R}/screen_summary.json"), jload(f"{NU}/screen_summary.json")
    fr, fn = jload(f"{R}/screen_final_summary.json"), jload(f"{NU}/screen_final_summary.json")
    Vr = cload(f"{R}/screen_final.csv")
    Vn = cload(f"{NU}/screen_final.csv")
    if not len(Vr):
        Vr = cload(f"{R}/screen.csv")
    if not len(Vn):
        Vn = cload(f"{NU}/screen.csv")
    Cr, Cn = cload(f"{C7}/confirm_real.csv"), cload(f"{C7}/confirm_null.csv")
    cr, cn = jload(f"{C7}/confirm_real_summary.json"), jload(f"{C7}/confirm_null_summary.json")
    refire = jload(f"{C7}/refire_meta.json")
    spr, spn = jload(f"{R}/chunks/_split.json"), jload(f"{NU}/chunks/_split.json")

    L = []
    w = L.append
    w("# Washington Puget / S-Olympic 2008, QUIET-SUBSET bed — station quality as the lever\n")
    w("Agent P, 2026-09-21.  **The recipe is unchanged** (`results/REPRO_BOSTOCK.md` §0): 6-s "
      "windows, 1–8 Hz, per-component cc averaged over all admitted stations of the bed, "
      "median + 8 × MAD of the block's own distribution, |lag| ≥ 6 s, 4-s declustering, "
      "connected-component grouping at the 4-s node merge (≥3 members), Bungum–Husebye "
      "realignment, 12-s station-tapered templates, network matched filter at 8 × MAD, no "
      "waveform-similarity merge.  Agent K's `k_*.py` are unedited and agent N's `n_*.py` drive "
      "them; **the only thing that changed is the STATION SET** (`n_common.admitted()` honours an "
      "optional `N_STATION_SUBSET` whitelist — unset it and every agent-N number is reproduced "
      "byte-for-byte).  The audited tremor hours, the literature-quiet days, the health veto, the "
      "screen and the 2007 confirmation rule are agent N's, reused unmodified.\n")
    w("**Why.**  On the 51-station CAFE+PB bed the pair stage had a tremor/quiet ratio of **0.89** "
      "(`results/WA_2008_RESULTS.md` §3) — it fired as often on a literature-quiet day as inside "
      "the audited May-2008 ETS — and every grown template was a noise repeater.  Royer & Bostock "
      "(2014) built *composite arrays* from the stations that worked.  §1 below measures, station "
      "by station, who is manufacturing those pairs; §2–§5 re-run the unchanged recipe on the "
      "quiet subset the census defines.\n")

    # ---------------------------------------------------------------- 1. census
    w("\n## 1. The per-station noise-repeater census\n")
    w("`scripts/p_census.py` recomputes, at **exactly** the (i, j) pairs agent N's pair stage "
      "stored, each station's own 3-component mean cc, and compares it with that station's own "
      "median + 8 × MAD over 4,000 random pairs of the same block (the recipe's own threshold "
      "rule applied to one station instead of to the network).  Validation: the mean of the "
      "per-station values reproduces the stored network cc to 7 × 10⁻⁸.  Two numbers per station:\n")
    w("* **ratio** = contributions-above-its-own-8 × MAD per audited TREMOR hour ÷ the same per "
      "literature-QUIET hour (counts Laplace-smoothed by 0.5 on both sides; `ratio4`, the "
      "10×-denser 4 × MAD version, is the stability check and is in the CSV);\n")
    w("* **top3_share_quiet** = the fraction of the QUIET-hour network detections in which this "
      "station is one of the top 3 contributors (chance = 3/N ≈ "
      f"{g(g(SM, 'rule', {}), 'top3_chance', float('nan')):.3f}).\n")
    w("A **noise-repeater station** (the brief's definition) has ratio < 2 and top3_share_quiet "
      f"above twice chance: **{g(SM, 'n_noise_repeater', '--')} of {len(CEN)} stations** qualify.\n")
    w("\n### 1.1 By network\n")
    w(md(BYN, ff=3))
    w("\n### 1.2 The eight worst (lowest ratio) and the eight best\n")
    if len(CEN):
        cols = ["netsta", "ratio", "ratio4", "per_hour_tremor", "per_hour_quiet",
                "top3_share_quiet", "noise_db", "noise_repeater"]
        cols = [c for c in cols if c in CEN.columns]
        w("**Worst (the pair stage's own noise engine):**\n")
        w(md(CEN.sort_values("ratio").head(8)[cols], ff=3))
        w("\n**Best (tremor-locked contributors):**\n")
        w(md(CEN.sort_values("ratio", ascending=False).head(8)[cols], ff=3))
    w("\n`noise_db` is the P4-0 statistic — the station's median 2–8 Hz hourly power on the "
      "literature-quiet days, in dB above the QUIETEST PB borehole "
      f"(0 dB = {g(SM, 'noise_ref_netsta_db0', float('nan')):.3g} in the cache's units).  "
      "Full table: `results/p_2008/station_noise_census.csv`; per-station-block rows in "
      "`census_tremor.parquet` / `census_quiet.parquet`.\n")

    # ---------------------------------------------------------------- 2. subset
    w("\n## 2. The quiet subset\n")
    rule = g(SM, "rule", {})
    w(f"Rule, fixed in `scripts/p_subset.py` before any subset result existed: **the "
      f"{g(SM, 'n_pb', '--')} PB boreholes in the bed, plus the CAFE (XU)/YW stations with "
      f"ratio ≥ {g(rule, 'ratio_bar', '--')} and top3_share_quiet ≤ "
      f"{g(rule, 'top3_bar', float('nan')):.3f}** (= 2 × chance), ranked by (ratio desc, "
      f"top3 asc), floor {g(rule, 'min_n', '--')} and cap {g(rule, 'max_n', '--')} stations.\n")
    w(f"* CAFE/YW stations clearing both bars: **{g(SM, 'n_cafe_qualified', '--')}**; "
      f"taken: {g(SM, 'n_cafe_picked', '--')}.\n")
    if g(SM, "n_below_bar_added", 0):
        w(f"* ⚠ {g(SM, 'shortfall_note', '')}.\n")
    w(f"* **SUBSET ({g(SM, 'n_subset', '--')} stations): "
      f"{', '.join(g(SM, 'subset', []))}**\n")
    w(f"* boreholes-only control bed ({g(SM, 'n_pb', '--')}): {', '.join(g(SM, 'pb', []))}\n")
    cov = pd.DataFrame([dict(bed="51-station (agent N)", **g(SM, "coverage_51", {})),
                        dict(bed=f"quiet subset ({g(SM, 'n_subset', '--')})",
                             **g(SM, "coverage_subset", {})),
                        dict(bed=f"boreholes only ({g(SM, 'n_pb', '--')})",
                             **g(SM, "coverage_pb_only", {}))])
    w("\n**Coverage of the Puget/S-Olympic active cells** (Slab2 28–45 km with ≥30 PNSN "
      "**2010–2011** epicentres within 8 km — there is no 2008 catalogue, so these are 2010–11 "
      "cells used for a 2008 bed, stated not hidden), reported as `P_STEP2_SITES.md` §4 does:\n")
    w(md(cov, ff=1))
    pc = g(g(SM, "coverage_subset", {}), "pct_cells_ge5_60km_hyp", float("nan"))
    if np.isfinite(pc) and pc < 60:
        w(f"\n⚠ **Read the screen with this number.**  Only {pc:.1f} % of the active cells have "
          f"≥5 subset stations within 60 km, against "
          f"{g(g(SM, 'coverage_51', {}), 'pct_cells_ge5_60km_hyp', float('nan')):.1f} % on the "
          "51-station bed.  Clause (c) of the screen requires a location on ≥5 stations, so a "
          "template from a real family in a poorly covered cell can fail this run for a "
          "GEOMETRIC reason.  The pair ratio in §3 is the measurement this bed was built to "
          "make; the screen counts are a lower bound, not a census of what is there.\n")
    if bed:
        w(f"\nBed as actually run: **{g(bed, 'n_stations', '--')} stations** "
          f"{json.dumps(g(bed, 'by_net', {}))}, {g(bed, 'n_tremor_hours', '--')} audited tremor "
          f"hours on {g(bed, 'n_tremor_days', '--')} days, {g(bed, 'n_quiet_days', '--')} "
          f"literature-quiet days.  2007 confirmation bed: {g(bed7, 'n_stations', '--')} "
          f"stations {json.dumps(g(bed7, 'by_net', {}))} "
          f"(absent from 2007: {', '.join(g(bed7, 'absent', [])) or 'none'}).\n")

    # ---------------------------------------------------------------- 3. pair stage
    w("\n## 3. The pair stage — the number this whole run was launched to move\n")
    rows = []
    for nm, d in (("51-station bed (agent N)", ps51), ("QUIET SUBSET", ps),
                  ("quiet subset — site-shift null", psn), ("boreholes only", psbh)):
        if not d:
            continue
        rows.append(dict(bed=nm,
                         stations=g(d, "sta_per_block_median"),
                         blocks_tremor=g(d, "n_blocks_tremor_clean"),
                         det_tremor=g(d, "n_det_tremor"),
                         per_hour_tremor=g(d, "det_per_hour_tremor"),
                         blocks_quiet=g(d, "n_blocks_quiet_clean"),
                         det_quiet=g(d, "n_det_quiet"),
                         per_hour_quiet=g(d, "det_per_hour_quiet"),
                         tremor_quiet_ratio=g(d, "tremor_quiet_ratio"),
                         nmad_in_sd=g(d, "nmad_in_sd_median")))
    w(md(pd.DataFrame(rows), ff=2))
    w("\nReference points: Vancouver Island POLARIS 2005 (the bed on which this recipe "
      "reproduces 58/96 of Bostock's families) **11×**; N-Olympic XG/PB 2010 **4.2×**; "
      "Puget 51-station 2008 **0.89×**.\n")
    if ps:
        w(f"\nnmad of the subset's tremor-hour detections: p50 {g(g(ps, 'nmad', {}), 'p50', float('nan')):.1f}, "
          f"p90 {g(g(ps, 'nmad', {}), 'p90', float('nan')):.1f}, "
          f"p99 {g(g(ps, 'nmad', {}), 'p99', float('nan')):.1f}, "
          f"max {g(g(ps, 'nmad', {}), 'max', float('nan')):.1f}; 8 × MAD sits at "
          f"{g(ps, 'nmad_in_sd_median', float('nan')):.2f} sd of the network function.\n")
    w("\n⚠ **Caveat that must be read with the ratio:** the quiet days are calendar-remote from "
      "the ETS window and the CAFE stations are less completely cached on them, so the quiet "
      "blocks run on fewer stations than the tremor blocks (`stations` column above is the "
      "tremor-block median; the quiet-block median is in `pair_pairstats.json` as "
      "`sta_per_block_median_quiet`).  A bed with fewer stations has a noisier network function "
      "and a higher absolute 8 × MAD bar, which **suppresses** quiet-hour detections; part of any "
      "ratio above 1 can be this and not tremor.  The template-level clauses of §4 — where the "
      "SAME template is fired on both day sets through identical code — are not exposed to it.\n")
    gr = pd.DataFrame([dict(chain="real", components_ge3=g(spr, "n_groups_total"),
                            grown=g(spr, "n_groups_grown")),
                       dict(chain="site-shift null", components_ge3=g(spn, "n_groups_total"),
                            grown=g(spn, "n_groups_grown"))])
    w("\n" + md(gr, ff=0))

    # ---------------------------------------------------------------- 4. screen
    w("\n## 4. The five-clause screen (agent N's, unchanged), real vs whole-chain null\n")
    w("Every clause is a 2010→2011 failure mode turned into a test (`scripts/n_screen.py` "
      "docstring): (a) per-template quiet-day false alarms ≤ 0.1/template-day AND ≥3× below the "
      "tremor-day rate; (b) tremor co-timing above the same-day base rate; (c) a ≥5-station "
      "location that passes the shuffled-station control, is INTERIOR to the wide grid "
      "(46.5–48.8 N / −125.2…−121.8 / 15–60 km, used from the start) and sits within 8 km of "
      "Slab2; (d) ≥3 distinct tremor days; (e) a second-half re-stack that relocates within "
      "10 km.  The site-shift null is the same code on a network whose inter-station coherence "
      "has been destroyed by within-partition day shifts.\n")
    if sr:
        w(f"\nNetwork false-alarm rates measured directly: real chain "
          f"**{g(sr, 'network_fa_per_tpl_day', float('nan')):.3f}** detections per template-day on "
          f"the literature-quiet days against "
          f"**{g(sr, 'network_rate_tremor_per_tpl_day', float('nan')):.3f}** on the audited tremor "
          f"days; null chain {g(sn, 'network_fa_per_tpl_day', float('nan')):.3f} vs "
          f"{g(sn, 'network_rate_tremor_per_tpl_day', float('nan')):.3f}.  "
          "(The 51-station bed read 1.525 vs 1.541; the 2010 templates re-fired in 2011 at 1.89 "
          "quiet vs 0.66 tremor — the SIGN of this comparison is the whole test.)\n")
    tr = {r["clause"]: r for r in g(sr, "clause_table", [])}
    tn = {r["clause"]: r for r in g(sn, "clause_table", [])}
    rows = []
    for c in ["c_c1", "c_c2", "c_c3", "c_c4", "c_d", "c_a", "c_b"]:
        if c not in tr:
            continue
        rows.append(dict(clause=tr[c]["name"], real_alone=tr[c]["alone"],
                         real_cumulative=tr[c]["cumulative"],
                         null_alone=tn.get(c, {}).get("alone", 0),
                         null_cumulative=tn.get(c, {}).get("cumulative", 0)))
    if fr:
        rows.append(dict(clause="(e) stability: 2nd-half re-stack within 10 km",
                         real_alone=g(fr, "n_pass_e", 0), real_cumulative=g(fr, "n_candidate", 0),
                         null_alone=g(fn, "n_pass_e", 0), null_cumulative=g(fn, "n_candidate", 0)))
    w("\n" + md(pd.DataFrame(rows), ff=0))
    nc_r = int(g(fr, "n_candidate", 0) or 0)
    nc_n = int(g(fn, "n_candidate", 0) or 0)
    w(f"\n> **CANDIDATES after the whole screen: REAL {nc_r} vs NULL {nc_n}.**\n")
    cols = [c for c in ["grp", "n_sta", "n_det", "n_days", "lat", "lon", "dep", "dz",
                        "fa_quiet_per_tpl_day", "rate_tremor_per_tpl_day", "co_timing",
                        "co_base", "h2_km"] if len(Vr) and c in Vr.columns]
    w("\n### 4.1 real chain — templates passing the whole screen\n")
    w(md(Vr[Vr.get("candidate", pd.Series(False, index=Vr.index)).astype(bool)][cols], ff=3)
      if len(Vr) else "_(empty)_\n")
    w("\n### 4.2 site-shift null — templates passing the whole screen\n")
    w(md(Vn[Vn.get("candidate", pd.Series(False, index=Vn.index)).astype(bool)][cols], ff=3)
      if len(Vn) and cols else "_(empty)_\n")
    w(f"\nFull per-template tables: `{R}/screen_final.csv` (real), `{NU}/screen_final.csv` "
      "(null); every clause is a column.\n")

    # ---------------------------------------------------------------- 5. 2007
    w("\n## 5. Confirmation on 2007-01\n")
    w("The rule was fixed before any 2007 number existed (`scripts/n_confirm.py` docstring): "
      "**≥3 distinct 2007 tremor days with detections, quiet-day false alarms ≤ 0.1 per "
      "template-day on the 2007 literature-quiet days, and an INTERIOR 2007 location within "
      "10 km (3-D) of the 2008 one.**  Templates are read straight out of the 2008 run's "
      "`grow/tpl_it0.npz` — no re-stacking, no re-cutting, no re-thresholding.  The surviving "
      "null templates go through byte-identical code on the identical 2007 data: their "
      "confirmation count IS the false-confirmation rate.\n")
    if refire:
        w(f"\nRe-fired: {g(refire, 'n_real', 0)} real and {g(refire, 'n_null', 0)} null templates "
          f"over {g(refire, 'n_tremor_days', 0)} audited 2007 tremor days and "
          f"{g(refire, 'n_quiet_days', 0)} 2007 literature-quiet days."
          + (f"  ⚠ {g(refire, 'note', '')}" if g(refire, "note", "") else "") + "\n")
    if cr or cn:
        w("\n" + md(pd.DataFrame([
            dict(chain="real", refired=g(cr, "n_refired", 0),
                 ge3_days=g(cr, "n_ge3_days", 0), fa_ok=g(cr, "n_fa_ok", 0),
                 loc_ok=g(cr, "n_loc_ok", 0), CONFIRMED=g(cr, "n_CONFIRMED", 0),
                 rate_tremor=g(cr, "rate_tremor_per_tpl_day", float("nan")),
                 rate_quiet=g(cr, "rate_quiet_per_tpl_day", float("nan"))),
            dict(chain="null", refired=g(cn, "n_refired", 0),
                 ge3_days=g(cn, "n_ge3_days", 0), fa_ok=g(cn, "n_fa_ok", 0),
                 loc_ok=g(cn, "n_loc_ok", 0), CONFIRMED=g(cn, "n_CONFIRMED", 0),
                 rate_tremor=g(cn, "rate_tremor_per_tpl_day", float("nan")),
                 rate_quiet=g(cn, "rate_quiet_per_tpl_day", float("nan")))]),
            ff=0))
    if len(Cr):
        w("\n" + md(Cr.head(40), ff=3))
    else:
        w("\n_(nothing survived the 2008 screen, so there was nothing to re-fire)_\n")

    # ---------------------------------------------------------------- 6. plain language
    w("\n## 6. What the census says about Vancouver Island vs Puget\n")
    if len(CEN) and "ratio" in CEN.columns:
        med_pb = float(np.nanmedian(CEN[CEN.net == "PB"].ratio)) if (CEN.net == "PB").any() \
            else np.nan
        med_xu = float(np.nanmedian(CEN[CEN.net != "PB"].ratio)) if (CEN.net != "PB").any() \
            else np.nan
        db_pb = float(np.nanmedian(CEN[CEN.net == "PB"].noise_db)) if (CEN.net == "PB").any() \
            else np.nan
        db_xu = float(np.nanmedian(CEN[CEN.net != "PB"].noise_db)) if (CEN.net != "PB").any() \
            else np.nan
    else:
        med_pb = med_xu = db_pb = db_xu = np.nan
    w(f"The network autocorrelation asks one question of every 6-second window: *does this "
      f"waveform repeat somewhere else in the hour, coherently, at many stations at once?*  LFEs "
      f"answer yes because the same patch of the plate interface ruptures the same way.  So does "
      f"a pump, a mill, a road, a badly-clamped sensor and a wind-loaded tree root — and those "
      f"repeat all day, every day, tremor or no tremor.  The census measures exactly that: the "
      f"CAFE surface stations sit a median "
      f"{db_xu - db_pb if np.isfinite(db_xu) and np.isfinite(db_pb) else float('nan'):.0f} dB "
      f"above the quietest borehole in the 2–8 Hz band on quiet days, and their median "
      f"tremor/quiet contribution ratio is {med_xu:.2f} against {med_pb:.2f} for the boreholes "
      f"— i.e. on the surface sites the repeating-waveform population is very nearly the same "
      f"whether or not there is an ETS underway.  Because the recipe SUMS the per-component "
      f"cross-correlation over every station, adding more such stations does not average the "
      f"noise repeaters away: it adds their coherent, always-on pairs to the sum, so the 8 × MAD "
      f"bar rises to accommodate them and the LFE pairs — which are only present during tremor "
      f"and only at the stations that see that patch — sit below it.  Bostock's southern "
      f"Vancouver Island stations are POLARIS forest sites on rock with almost no cultural "
      f"neighbours; there the same sum is dominated by whatever genuinely repeats, which is the "
      f"LFEs, and the tremor/quiet ratio is 11×.  This is why 'more stations' made Washington "
      f"worse (0.89× on 51 CAFE+PB sites) and why the published Washington work built *composite "
      f"arrays* out of chosen stations rather than throwing the whole deployment into one sum.\n")
    w("\nWhat this run then tests is the direct consequence: if the obstruction is the noise "
      "repeaters, deleting them should restore the pair stage. §3 reports whether it did, and "
      "§4–§5 report whether a restored pair stage is enough to produce a template that survives "
      "a quiet-day false-alarm test, a location, a stability re-stack and an out-of-sample "
      "re-fire — which is the bar a family has to clear, not the pair ratio.\n")

    w("\n## 7. Files\n")
    w(md(pd.DataFrame([
        dict(what="per-station-block census rows",
             where=f"{P}/census_tremor.parquet, census_quiet.parquet"),
        dict(what="station census + by-network table",
             where=f"{P}/station_noise_census.csv, census_by_network.csv"),
        dict(what="subset definition, coverage, rule", where=f"{P}/subset_meta.json, "
                                                             "subset_stations.txt"),
        dict(what="pair stage (tremor / quiet)", where=f"{R}/auto/, auto_quiet/, "
                                                       "pair_all.parquet, pair_pairstats.json"),
        dict(what="groups", where=f"{R}/gb_nodes.csv, gb_summary.csv"),
        dict(what="templates + tremor-day detections", where=f"{R}/grow/tpl_it0.npz, "
                                                             "det_it0.parquet"),
        dict(what="quiet-day detections (clause a)", where=f"{R}/det_fa.parquet"),
        dict(what="locations (wide grid, edge-flagged)", where=f"{R}/loc_it0.csv"),
        dict(what="screen tables", where=f"{R}/screen.csv, screen_final.csv"),
        dict(what="whole-chain site-shift null", where=f"{NU}/ (offsets null_offsets.json)"),
        dict(what="boreholes-only pair stage", where=f"{BH}/"),
        dict(what="2007 re-fire + confirmation", where=f"{C7}/"),
    ])))

    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    open(a.out, "w").write("\n".join(L))
    print(f"wrote {a.out} ({sum(len(x) for x in L)} chars)")


if __name__ == "__main__":
    main()
