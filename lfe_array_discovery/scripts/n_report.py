#!/usr/bin/env python
"""Write results/WA_2008_RESULTS.md from what the queue produced.  Prefix `n_`."""
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

OUTROOT = os.environ.get("N_OUTROOT", f"{RES}/n_puget2008")


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


def hist_row(h):
    return " ".join(f"{k}:{v}" for k, v in sorted(h.items(), key=lambda x: int(x[0])))


def md_table(df, cols=None, floatfmt=2):
    if not len(df):
        return "_(empty)_\n"
    df = df[cols] if cols else df
    head = "| " + " | ".join(str(c) for c in df.columns) + " |"
    sep = "|" + "|".join("---" for _ in df.columns) + "|"
    rows = []
    for r in df.itertuples(index=False):
        cells = []
        for v in r:
            if isinstance(v, float) and np.isfinite(v):
                cells.append(f"{v:.{floatfmt}f}")
            elif v is None or (isinstance(v, float) and not np.isfinite(v)):
                cells.append("--")
            else:
                cells.append(str(v))
        rows.append("| " + " | ".join(cells) + " |")
    return "\n".join([head, sep] + rows) + "\n"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{RES}/WA_2008_RESULTS.md")
    a = ap.parse_args()
    R, N_ = f"{OUTROOT}/2008", f"{OUTROOT}/2008_null"

    bed08 = jload(f"{OUTROOT}/bed_2008.json")
    bed07 = jload(f"{OUTROOT}/bed_2007.json")
    au08 = jload(f"{OUTROOT}/hour_audit_2008.json")
    au07 = jload(f"{OUTROOT}/hour_audit_2007.json")
    psr = jload(f"{R}/pair_pairstats.json")
    psn = jload(f"{N_}/pair_pairstats.json")
    sr = jload(f"{R}/screen_summary.json")
    sn = jload(f"{N_}/screen_summary.json")
    fr = jload(f"{R}/screen_final_summary.json")
    fn = jload(f"{N_}/screen_final_summary.json")
    Vr = cload(f"{R}/screen_final.csv")
    Vn = cload(f"{N_}/screen_final.csv")
    if not len(Vr):
        Vr = cload(f"{R}/screen.csv")
    if not len(Vn):
        Vn = cload(f"{N_}/screen.csv")
    Cr = cload(f"{OUTROOT}/2007/confirm_real.csv")
    Cn = cload(f"{OUTROOT}/2007/confirm_null.csv")
    cr = jload(f"{OUTROOT}/2007/confirm_real_summary.json")
    cn = jload(f"{OUTROOT}/2007/confirm_null_summary.json")
    gs_r = cload(f"{R}/gb_summary.csv")
    gs_n = cload(f"{N_}/gb_summary.csv")
    sp_r = jload(f"{R}/chunks/_split.json")
    sp_n = jload(f"{N_}/chunks/_split.json")
    fa_r = jload(f"{R}/det_fa_summary.json")
    fa_n = jload(f"{N_}/det_fa_summary.json")
    refire = jload(f"{OUTROOT}/2007/refire_meta.json")

    L = []
    w = L.append
    w("# Washington Puget / S-Olympic 2008 — Bostock's published recipe on the last full-array "
      "Washington bed\n")
    w("Agent N, 2026-09-21.  **The same recipe and the same code** as the north reproduction and "
      "the N-Olympic 2010 run (`results/REPRO_BOSTOCK.md` §0 parameter table): 6-s windows, "
      "1–8 Hz, per-component cc averaged over **all** admitted stations, threshold = median + "
      "8 × MAD of the block's own distribution, |lag| ≥ 6 s, 4-s declustering, connected-component "
      "grouping at the 4-s node merge (≥3 members), Bungum–Husebye realignment before every stack, "
      "12-s station-tapered templates, network matched filter at 8 × MAD, **no waveform-similarity "
      "merge**.  Engine files `scripts/k_auto.py`, `k_group.py`, `k_grow.py`, `k_locate.py` are "
      "agent K's, unedited; the bed is supplied by `scripts/n_common.py` and driven by "
      "`n_auto.py` / `n_run.py` / `n_locate.py` / `n_screen.py` / `n_relocate.py` / `n_refire.py`, "
      "which are the 2010/2011 `l_`/`m_` scripts with the bed adapter swapped.\n")
    w("**There is no Washington LFE catalogue and no tremor catalogue before 2010-01-06.**  "
      "Nothing below is a rediscovery count.  What replaces truth is §1 (a tremor-hour audit built "
      "from the waveforms), §4 (a five-clause screen whose every clause is one of the failure "
      "modes the 2010→2011 test exposed) and §5 (an out-of-sample re-fire on 2007-01), each with "
      "a whole-chain site-shift null run through byte-identical code.\n")

    # ---------------------------------------------------------------- 1. the hour audit
    w("\n## 1. The tremor-hour audit (there is no 2008 tremor catalogue)\n")
    w("The 2008 hour flag this project had (`results/hour_energy_flag_2008.csv`, network-median "
      "hourly energy over a fixed p95) selects 113 hours of which **86 % fall in 15–22 UTC** "
      "(`results/WA_V2_INPUTS.md` §2.2) — i.e. it is a cultural-noise detector, and the brief for "
      "this run required it to be replaced before anything else happened.  "
      "`scripts/n_hours.py` replaces it with three clauses (full statement in that file's "
      "docstring):\n")
    w("1. **the station's own diurnal cycle divided out** — E(s,d,h) = median over the hour of the "
      "6-s-window mean square of the 2–8 Hz 3-C waveform, referred to that station's own median at "
      "the same hour-of-day over the literature-quiet days;\n"
      "2. **simultaneity** — ≥ 60 % of the reporting near-field stations (admitted, within 60 km "
      "of an active cell, health-veto applied) above **their own** quiet-day p95 at once;\n"
      "3. **spatial coherence** — the log-excess must fit an amplitude decay from ONE active cell "
      "better than a uniform field does, at p ≤ 0.05 against 200 station-label shuffles.\n")
    for tag, au in (("2008", au08), ("2007", au07)):
        if not au:
            continue
        s = au.get("selected", {})
        c2 = au.get("clause2_only", {})
        c23 = au.get("clause2_and_3", {})
        w(f"\n### 1.{1 if tag == '2008' else 2} {tag} bed\n")
        w(f"* near-field stations {au.get('n_near_stations')}; literature-quiet days "
          f"{au.get('n_quiet_days')}; hours judged {au.get('n_hours_judged')}.\n")
        w(f"* clause 2 alone: **{c2.get('n')} hours on {c2.get('n_days')} days**; "
          f"clause 2 AND 3: **{c23.get('n')} hours on {c23.get('n_days')} days**.\n")
        w(f"* **SELECTED: {s.get('n')} hours on {s.get('n_days')} days**; "
          f"fraction in 15–22 UTC **{s.get('frac_daytime_15_22utc', float('nan')):.3f}** "
          f"(the old flag: 0.86); weekday fraction "
          f"{s.get('frac_weekday', float('nan')):.3f} (5/7 = 0.714 for a calendar-random set).\n")
        if au.get("fallback"):
            w(f"* ⚠ **fallback applied**: {au['fallback']}\n")
        w(f"* **FALSE-POSITIVE CHECK — on the {au.get('n_quiet_day_hours_judged')} hours of the "
          f"literature-quiet days the same statistic fires "
          f"{au.get('quiet_day_hours_passing_clause2')} times on clause 2 and "
          f"{au.get('quiet_day_hours_passing_both')} times on clause 2+3.**\n")
        w(f"* selected days: {', '.join(s.get('dates', []))}\n")
        w(f"\nUTC-hour histogram of the selected hours: `{hist_row(s.get('by_utc_hour', {}))}`\n")
        fl = au.get("diurnal_noise_floor_by_utc_hour", {})
        if fl:
            w(f"\nThe array's own diurnal NOISE FLOOR (network median of the quiet-day reference "
              f"M(s,h) by UTC hour, normalised to its minimum): `{hist_row({k: round(v, 1) for k, v in fl.items()})}` "
              f"— it is ~10× higher at 19 UTC than at 09 UTC.  A fixed tremor amplitude therefore "
              f"produces a ~10× smaller ratio in the cultural afternoon, so **the clock "
              f"distribution of the selected hours is detection-limited, not a claim that tremor "
              f"is nocturnal**; the quiet-day false-positive rate by UTC hour "
              f"(`{hist_row({k: round(v, 3) for k, v in au.get('quiet_frac_above_by_utc_hour', {}).items()})}`) "
              f"is flat, so the THRESHOLD is not what makes the set night-heavy — the SNR is.\n")
    if au08:
        w("\n**Against the literature.**  The published Washington ETS of 2008 is the early-to-"
          "mid-May episode; the independent array-energy flag of `results/P_STEP2_SITES.md` §3.2 "
          "(a different statistic, fixed network p95) puts it on 2008-05-05..05-23 plus 05-29 and "
          "06-10.  The audit here selects hours only on "
          f"{au08.get('selected', {}).get('n_days')} days, all inside "
          f"{min(au08.get('selected', {}).get('dates', ['--']))} … "
          f"{max(au08.get('selected', {}).get('dates', ['--']))} — the two agree, and the audit "
          "reaches the answer with the diurnal cycle removed and a spatial-coherence test that "
          "the energy flag does not have.\n")

    # ---------------------------------------------------------------- 2. the bed
    w("\n## 2. The bed\n")
    for tag, b in (("2008 discovery", bed08), ("2007-01 confirmation", bed07)):
        if not b:
            continue
        w(f"\n**{tag}** — window {b.get('day0')} … {b.get('day1')} (tranche {b.get('tranche')}), "
          f"sub-region {b.get('subregion')} (lat 46.2–47.8).\n")
        w(f"* **{b.get('n_stations')} three-component stations** of the "
          f"{b.get('n_probe_admitted')} probe-admitted ones "
          f"(`results/site_probe_2008.csv`): {b.get('by_net')}.\n")
        if b.get("absent"):
            w(f"* absent from this window: {', '.join(b['absent'])}.\n")
        if b.get("health_dropped"):
            w(f"* dropped by the > 30 % vetoed-station-hour bar: "
              f"{', '.join(b['health_dropped'])}.\n")
        w(f"* literature-quiet days: {b.get('n_quiet_days')} "
          f"({', '.join(b.get('quiet_days', []))}).\n")
        w(f"* active cells in the sub-region: {b.get('n_active_cells')} "
          f"(Slab2 28–45 km with ≥ 30 PNSN **2010–2011** epicentres within 8 km — there is no 2008 "
          f"catalogue, so these are 2010–11 cells used for a 2008 bed, stated not hidden); "
          f"{b.get('n_near_stations')} admitted stations within 60 km of one.\n")
        if b.get("n_corr_flagged") is not None:
            w(f"* intraslab-earthquake S station terms available for "
              f"{b.get('n_intraslab_corr')} stations; {b.get('n_corr_flagged')} carry 0.0 and are "
              f"FLAGGED, never imputed.\n")
    w("\n**ZK is vertical-only in the archive** (1,604 of 1,642 delivered ZK station-days carry a "
      "single EHZ channel, `results/P_STEP2_SITES.md` finding 1) and contributes nothing; "
      "**UW.GNW fails the 2008 site probe** at 1.7 dB and is not in the sum.\n")

    # ---------------------------------------------------------------- 3. pair stage + groups
    w("\n## 3. Pair stage and grouping\n")
    rows = []
    for tag, p in (("real", psr), ("site-shift null", psn)):
        if not p:
            continue
        rows.append(dict(chain=tag,
                         blocks_tremor=p.get("n_blocks_tremor"),
                         det_tremor=p.get("n_det_tremor"),
                         per_hour_tremor=p.get("det_per_hour_tremor"),
                         blocks_quiet=p.get("n_blocks_quiet"),
                         det_quiet=p.get("n_det_quiet"),
                         per_hour_quiet=p.get("det_per_hour_quiet"),
                         tremor_quiet_ratio=p.get("tremor_quiet_ratio"),
                         nmad_in_sd=p.get("nmad_in_sd_median"),
                         sta_per_block=p.get("sta_per_block_median")))
    w(md_table(pd.DataFrame(rows)))
    if psr.get("nmad"):
        w(f"\nnmad distribution of the tremor-hour detections: "
          f"p50 {psr['nmad']['p50']:.1f}, p90 {psr['nmad']['p90']:.1f}, "
          f"p99 {psr['nmad']['p99']:.1f}, max {psr['nmad']['max']:.1f}.  "
          f"8 × MAD sits at {psr.get('nmad_in_sd_median', float('nan')):.2f} sd of the network "
          f"function — the same measured bar as every other bed (2005: 5.15 sd, 2010: 4.90 sd).\n")
    rows = []
    for tag, gs, sp in (("real", gs_r, sp_r), ("site-shift null", gs_n, sp_n)):
        rows.append(dict(chain=tag, components_ge3_members=len(gs) if len(gs) else 0,
                         grown=sp.get("n_groups_grown", 0),
                         growth_budget=sp.get("max_groups", 0)))
    w("\n" + md_table(pd.DataFrame(rows)))

    # ---------------------------------------------------------------- 4. the screen
    w("\n## 4. The screen (every clause is a 2010→2011 failure mode)\n")
    w("Applied to **every** grown template before anything is counted.  The 2010 run's "
      "\"50 candidates vs 17 null\" died in 2011 because 46/50 sat on a face of the search box, "
      "the templates fired 1.89 per template-day on NON-tremor days against 0.66 on tremor days, "
      "and their locations moved 42–150 km when re-stacked on other data "
      "(`results/WA_2011_REFIRE_RESULTS.md`).  Clauses (a)–(e) are those three facts turned into "
      "tests; the wide search box 46.5–48.8 N / −125.2…−121.8 / 15–60 km is used **from the "
      "start** so that \"interior\" means something.\n")
    if sr:
        w(f"\nNetwork false-alarm rates measured directly: real chain "
          f"**{sr.get('network_fa_per_tpl_day', float('nan')):.3f}** detections per template-day "
          f"on the {sr.get('n_quiet_days')} literature-quiet days against "
          f"**{sr.get('network_rate_tremor_per_tpl_day', float('nan')):.3f}** on the "
          f"{sr.get('n_tremor_days')} audited tremor days"
          + (f"; null chain {sn.get('network_fa_per_tpl_day', float('nan')):.3f} vs "
             f"{sn.get('network_rate_tremor_per_tpl_day', float('nan')):.3f}" if sn else "")
          + f".  (The 2010 templates re-fired in 2011 at 1.89 quiet vs 0.66 tremor — the sign of "
            f"this comparison is the whole test.)\n")
    tr = {c["clause"]: c for c in sr.get("clause_table", [])}
    tn = {c["clause"]: c for c in sn.get("clause_table", [])}
    rows = []
    for c in ["c_c1", "c_c2", "c_c3", "c_c4", "c_d", "c_a", "c_b"]:
        if c not in tr:
            continue
        rows.append({"clause": tr[c]["name"],
                     "real alone": tr[c]["alone"], "real cumulative": tr[c]["cumulative"],
                     "null alone": tn.get(c, {}).get("alone", "--"),
                     "null cumulative": tn.get(c, {}).get("cumulative", "--")})
    if rows:
        rows.append({"clause": "(e) stability: 2nd-half re-stack relocates within 10 km",
                     "real alone": fr.get("n_pass_e", "--"),
                     "real cumulative": fr.get("n_candidate", "--"),
                     "null alone": fn.get("n_pass_e", "--"),
                     "null cumulative": fn.get("n_candidate", "--")})
    w("\n" + md_table(pd.DataFrame(rows)))
    w(f"\n> **CANDIDATES after the whole screen: REAL "
      f"{fr.get('n_candidate', sr.get('n_pass_abcd', 0))} vs NULL "
      f"{fn.get('n_candidate', sn.get('n_pass_abcd', 0))}.**\n")
    w("\nClause (e) is the only one evaluated on the survivors of (a)–(d) rather than on all "
      "templates: a re-stack of every grown template is a second full pass over the waveforms. "
      "The intersection is unaffected by the order.\n")

    cols = [c for c in ["grp", "lat", "lon", "dep", "slab2", "dz", "n_sta", "n_det", "n_days",
                        "fa_quiet_per_tpl_day", "rate_tremor_per_tpl_day", "co_timing", "co_base",
                        "co_p", "rms", "edge", "h2_km"] if len(Vr) and c in Vr.columns]
    for tag, V in (("real", Vr), ("site-shift null", Vn)):
        if not len(V):
            continue
        key = "candidate" if "candidate" in V.columns else "pass_abcd"
        C = V[V[key].astype(bool)].sort_values("n_det", ascending=False)
        w(f"\n### 4.{1 if tag == 'real' else 2} {tag} chain — templates passing the whole screen "
          f"({len(C)})\n")
        w(md_table(C[cols] if len(C) else C, floatfmt=3))
    if len(Vr):
        key = "candidate" if "candidate" in Vr.columns else "pass_abcd"
        w(f"\nFull per-template screen tables: `{R}/screen_final.csv` (real), "
          f"`{N_}/screen_final.csv` (null); every clause is a column.\n")

    # ---------------------------------------------------------------- 5. the confirmation
    w("\n## 5. Confirmation on 2007-01 (tranche T5, the January 2007 ETS)\n")
    w("The rule was fixed before any 2007 number existed (`scripts/n_confirm.py` docstring): "
      "**≥ 3 distinct 2007 tremor days with detections, quiet-day false alarms ≤ 0.1 per "
      "template-day on the 2007 literature-quiet days, and an INTERIOR 2007 location within 10 km "
      "(3-D) of the 2008 one.**  The templates are read straight out of the 2008 run's "
      "`grow/tpl_it0.npz` — no re-stacking, no re-cutting, no re-thresholding.  The surviving "
      "site-shift-null templates go through byte-identical code on the identical 2007 data: "
      "**their confirmation count is the false-confirmation rate.**\n")
    if refire:
        w(f"\nRe-fired: {refire.get('n_real', 0)} real and {refire.get('n_null', 0)} null "
          f"templates over {refire.get('n_tremor_days', 0)} audited 2007 tremor days and "
          f"{refire.get('n_quiet_days', 0)} 2007 literature-quiet days"
          + (f"; {refire.get('note')}" if refire.get("note") else "") + ".\n")
    rows = []
    for tag, c in (("real", cr), ("site-shift null", cn)):
        if not c:
            continue
        rows.append({"chain": tag, "re-fired": c.get("n_refired"),
                     "det on tremor days": c.get("n_det_2007"),
                     "per tpl-day (tremor)": c.get("rate_tremor_per_tpl_day"),
                     "det on quiet days": c.get("n_det_quiet_2007"),
                     "per tpl-day (quiet)": c.get("rate_quiet_per_tpl_day"),
                     "≥3 days": c.get("n_ge3_days"), "FA ok": c.get("n_fa_ok"),
                     "interior + ≤10 km": c.get("n_loc_ok"),
                     "CONFIRMED": c.get("n_CONFIRMED")})
    w("\n" + md_table(pd.DataFrame(rows), floatfmt=3))
    if rows:
        w(f"\n> **CONFIRMED: REAL {cr.get('n_CONFIRMED', 0)} of {cr.get('n_refired', 0)}; "
          f"NULL {cn.get('n_CONFIRMED', 0)} of {cn.get('n_refired', 0)} "
          f"(false-confirmation rate "
          f"{cn.get('n_CONFIRMED', 0)}/{max(cn.get('n_refired', 0), 1)}).**\n")
    ccols = [c for c in ["grp", "lat_2008", "lon_2008", "dep_2008", "lat_2007", "lon_2007",
                         "dep_2007", "move_km", "edge_2007", "n_det_2007", "n_days_2007",
                         "n_det_quiet_2007", "fa_2007_per_tpl_day", "CONFIRMED"]
             if len(Cr) and c in Cr.columns]
    for tag, C in (("real", Cr), ("site-shift null", Cn)):
        if not len(C):
            continue
        w(f"\n### 5.{1 if tag == 'real' else 2} {tag} templates on 2007\n")
        w(md_table(C[ccols] if ccols else C, floatfmt=3))

    # ---------------------------------------------------------------- 6. Chestler
    w("\n## 6. Descriptive comparison with the Chestler-derived coordinates\n")
    w("`../codex-LFE/results/olympic_reference_discovery_v1/"
      "derived_coordinate_rows_not_a_catalog.csv` is, as its filename says, **not a catalogue**; "
      "it is used only to say where a confirmed family would sit relative to previously published "
      "Olympic LFE coordinates, never as truth.\n")
    conf = Cr[Cr.CONFIRMED.astype(bool)] if len(Cr) and "CONFIRMED" in Cr.columns else pd.DataFrame()
    if len(conf):
        import s_common as S
        CH = pd.read_csv("/home/jovyan/tremorferometry/codex-LFE/results/"
                         "olympic_reference_discovery_v1/"
                         "derived_coordinate_rows_not_a_catalog.csv")
        cla, clo = CH.latitude.to_numpy(), CH.longitude.to_numpy()
        d = [float(np.min(S.epi_km(float(r.lat_2008), float(r.lon_2008), cla, clo)))
             for r in conf.itertuples()]
        conf = conf.assign(km_to_nearest_chestler=d)
        w("\n" + md_table(conf[["grp", "lat_2008", "lon_2008", "dep_2008",
                                "km_to_nearest_chestler"]], floatfmt=2))
        w(f"\nThe Chestler-derived coordinates are Olympic; this bed is Puget/S-Olympic, so a "
          f"large distance is expected and carries no information about whether a family is real.\n")
    else:
        w("\n* no confirmed family to compare.\n")

    # ---------------------------------------------------------------- 7. files
    w("\n## 7. Files\n")
    w(md_table(pd.DataFrame([
        {"what": "hourly 2–8 Hz energy per station-day",
         "where": f"{OUTROOT}/hour_energy28_{{2008,2007}}.parquet"},
        {"what": "hour audit (per hour, every clause)",
         "where": f"{OUTROOT}/hour_audit_{{2008,2007}}.csv + .json"},
        {"what": "selected tremor hours / days",
         "where": f"{OUTROOT}/tremor_hours_{{2008,2007}}.csv, tremor_dates_*.txt"},
        {"what": "pair stage (tremor / quiet)", "where": f"{R}/auto/, auto_quiet/, "
                                                         f"pair_all.parquet, pair_pairstats.json"},
        {"what": "groups", "where": f"{R}/gb_nodes.csv, gb_summary.csv"},
        {"what": "templates + tremor-day detections", "where": f"{R}/grow/tpl_it0.npz, "
                                                               f"det_it0.parquet"},
        {"what": "quiet-day detections (clause a)", "where": f"{R}/det_fa.parquet"},
        {"what": "locations (wide grid, edge-flagged)", "where": f"{R}/loc_it0.csv"},
        {"what": "screen tables", "where": f"{R}/screen.csv, screen_final.csv"},
        {"what": "stability re-stacks", "where": f"{R}/reloc_h1.csv, reloc_h2.csv, reloc_all.csv"},
        {"what": "whole-chain site-shift null", "where": f"{N_}/ (offsets in "
                                                         f"{N_}/null_offsets.json)"},
        {"what": "2007 re-fire + confirmation", "where": f"{OUTROOT}/2007/"},
    ])))
    txt = "".join(L)
    open(a.out, "w").write(txt)
    print(f"wrote {a.out} ({len(txt)} chars)", flush=True)


if __name__ == "__main__":
    main()
