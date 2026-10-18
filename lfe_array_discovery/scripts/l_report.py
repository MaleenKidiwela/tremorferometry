#!/usr/bin/env python
"""Write results/WA_2010_RESULTS.md and append the run to NOTES.md.  Prefix `l_`."""
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
sys.path.insert(0, HERE)
import l_common as L                                                        # noqa: E402


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


def fmt(v, n=3):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "—"
    if isinstance(v, float):
        return f"{v:,.{n}f}"
    if isinstance(v, (int, np.integer)):
        return f"{int(v):,}"
    return str(v)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real", default=f"{STUDY}/results/l_wa2010")
    ap.add_argument("--null", default=f"{STUDY}/results/l_wa2010_null")
    ap.add_argument("--out", default=f"{STUDY}/results/WA_2010_RESULTS.md")
    ap.add_argument("--no-notes", action="store_true", help="dry run: do not append NOTES.md")
    a = ap.parse_args()
    R, N = a.real, a.null
    bed = jload(f"{R}/bed.json")
    ps = jload(f"{R}/pair_pairstats.json")
    gm = jload(f"{R}/gb_meta.json")
    sp = jload(f"{R}/chunks/_split.json")
    loc = cload(f"{R}/loc_it0.csv")
    ver = cload(f"{R}/verify_it0.csv")
    vs = jload(f"{R}/verify_it0_summary.json")
    nloc = cload(f"{N}/loc_it0.csv")
    nver = cload(f"{N}/verify_it0.csv")
    nvs = jload(f"{N}/verify_it0_summary.json")
    nsp = jload(f"{N}/chunks/_split.json")
    nps = jload(f"{N}/pair_pairstats.json")
    ngm = jload(f"{N}/gb_meta.json")
    par = jload(f"{R}/auto/_params.json")
    off = jload(f"{N}/null_offsets.json")
    hrs = cload(f"{R}/candidate_hours.csv")

    cand = ver[ver.candidate] if "candidate" in ver.columns else pd.DataFrame()
    ncand = nver[nver.candidate] if "candidate" in nver.columns else pd.DataFrame()

    o = []
    A = o.append
    A("# Washington N-Olympic 2010 — Bostock's published recipe, run blind (no catalogue truth)")
    A("")
    A("Agent L, 2026-09-20.  The **same recipe and the same code** as the north reproduction "
      "(`results/REPRO_BOSTOCK.md` §0 parameter table, `results/NORTH_DISCOVERY_SUMMARY.md`): "
      "6-s windows, 1–8 Hz, per-component cc averaged over **all** admitted stations, "
      "threshold = median + 8 × MAD of the block's own distribution, |lag| ≥ 6 s, 4-s "
      "declustering, connected-component grouping at the 4-s node merge (≥3 members), "
      "Bungum–Husebye realignment before every stack, 12-s station-tapered templates cut 2 s "
      "before the stack envelope peak, network matched filter at 8 × MAD, **no "
      "waveform-similarity merge**.  Engine files `scripts/k_auto.py`, `k_group.py`, `k_grow.py`, "
      "`k_locate.py` are agent K's, unedited; the WA bed is supplied by `scripts/l_common.py` "
      "(station table, borehole azimuths, net-qualified cache index, station-hour health mask, "
      "day-shift null) and driven by `l_auto.py` / `l_run.py` / `l_locate.py` / `l_verify.py`.")
    A("")
    A("**There is no Washington LFE catalogue.**  Bostock's ≥30 % one-to-one criterion cannot be "
      "applied, so nothing below is a rediscovery count and no family here is claimed to be real "
      "on the strength of a match.  What replaces truth is stated in §5: location with a "
      "shuffled-station control, repeat/co-timing statistics against their own base rates, and "
      "the identical chain run once on a site-shifted network.")
    A("")
    A("## 1. The bed")
    A("")
    A(f"* window **{bed.get('day0')} … {bed.get('day1')}**, sub-region {bed.get('subregion')} "
      f"({L.LAT0}–{L.LAT1} N, {L.LON0}…{L.LON1}).")
    A(f"* **{fmt(bed.get('n_stations'))} three-component stations**, every probe-admitted station "
      f"(`results/site_probe_2010.csv`) except the one that fails the station-level health bar "
      f"(> 30 % vetoed station-hours): {', '.join(bed.get('health_dropped', [])) or 'none'}.  "
      "No station is withheld — the north's one real lever is station count "
      "(11 → 15 stations raised true pairs 8.4 % → 12.9 %, ledger §8).")
    if bed.get("by_net"):
        A("")
        A("| network | stations |")
        A("|---|---|")
        for k, v in sorted(bed["by_net"].items()):
            A(f"| {k} | {v} |")
        A(f"| **total** | **{fmt(bed.get('n_stations'))}** |")
    A("")
    A(f"* PNSN candidate (tremor) hours: **{fmt(len(hrs))}** on "
      f"**{fmt(hrs.date.nunique() if len(hrs) else 0)}** days "
      f"(≥5 epicentres within 50 km of an active cell of the sub-region).")
    A(f"* hour-matched, earthquake-screened quiet hours: "
      f"**{fmt(ps.get('n_blocks_quiet', bed.get('n_quiet_hours')))}** "
      "(`results/wa_quiet_hours_2010.csv`, WA_V2_INPUTS §2.1).")
    A(f"* intraslab-earthquake S station terms available for "
      f"**{fmt(bed.get('n_intraslab_corr'))}** of the {fmt(bed.get('n_stations'))} stations; "
      f"{fmt(bed.get('n_corr_flagged'))} carry 0.0 and are FLAGGED, never imputed "
      "(`results/x_station_corrections.csv`).")
    A("* per-station-hour health veto applied as a mask in every autocorrelation block "
      "(`results/wa_station_hour_veto_2010.csv`).")
    A("")
    A("## 2. Pair stage (autocorrelation at 8 × MAD)")
    A("")
    A("| | tremor hours | quiet hours |")
    A("|---|---|---|")
    A(f"| blocks run | {fmt(ps.get('n_blocks_tremor'))} | {fmt(ps.get('n_blocks_quiet'))} |")
    A(f"| blocks rejected (8×MAD bar not > 3 sd) | {fmt(ps.get('n_blocks_tremor_rejected'))} | — |")
    A(f"| stations per block (median) | {fmt(ps.get('sta_per_block_median'), 0)} | "
      f"{fmt(ps.get('sta_per_block_median_quiet'), 0)} |")
    A(f"| 8 × MAD in sd of the network function (median) | {fmt(ps.get('nmad_in_sd_median'), 2)} | "
      f"{fmt(ps.get('nmad_in_sd_median_quiet'), 2)} |")
    A(f"| detections | {fmt(ps.get('n_det_tremor'))} | {fmt(ps.get('n_det_quiet'))} |")
    A(f"| detections per clean hour | {fmt(ps.get('det_per_hour_tremor'), 1)} | "
      f"{fmt(ps.get('det_per_hour_quiet'), 1)} |")
    A("")
    A(f"**Tremor-vs-quiet pair-count ratio at 8 × MAD = "
      f"{fmt(ps.get('tremor_quiet_ratio'), 2)}×** "
      f"(top-2,000 tier: {fmt(ps.get('tremor_quiet_ratio_above_nmad2000'), 2)}×).  "
      "The north's 11× is a *truth-based* ratio (12.85 % true pairs against a 1.17 % shifted-"
      "catalogue chance rate, ledger §8) and is **not the same quantity** — with no catalogue the "
      "only available contrast is tremor hours against catalogue-quiet hours, which is a weaker "
      "test because a quiet hour is quiet in tremor, not in noise.")
    A("")
    A("8 × MAD sits at "
      f"{fmt(ps.get('nmad_in_sd_median'), 2)} sd of the network function, against 5.13–5.24 sd on "
      "every north hour — the rule is a measured bar, not a high one, on this bed too.")
    if ps.get("nmad"):
        A("")
        A("nmad distribution of the tremor-hour detections: "
          + ", ".join(f"{k} {fmt(v, 1)}" for k, v in ps["nmad"].items()) + ".")
        A("")
        A("| tier by nmad | nmad floor | hours represented | days |")
        A("|---|---|---|---|")
        for n in (100, 500, 2000, 10000, 50000):
            if f"top{n}_nmad_min" in ps:
                A(f"| top {n:,} | {fmt(ps[f'top{n}_nmad_min'], 1)} | "
                  f"{fmt(ps[f'top{n}_hours'])} | {fmt(ps[f'top{n}_days'])} |")
    A("")
    A("Every autocorrelation block is **one hour**, so both windows of every pair lie inside the "
      "same hour by construction and the pair-stage \"fraction in a tremor-active hour\" is 1.0 "
      "with no information in it.  The informative version — the fraction of **grown-template** "
      "detections that land in a PNSN tremor-active hour, against that quantity's own base rate "
      "— is in §5.")
    A("")
    A("## 3. Grouping")
    A("")
    A(f"* nodes (detections merged at the 4-s declustering interval): {fmt(gm.get('n_nodes'))}")
    A(f"* connected components: {fmt(gm.get('n_components'))}")
    A(f"* components with ≥3 members: **{fmt(gm.get('n_kept'))}**")
    if sp:
        A(f"* grown: **{fmt(sp.get('n_groups_grown'))}** of {fmt(sp.get('n_groups_total'))} "
          f"(growth budget {fmt(sp.get('max_groups'))}; the GPU is shared, and the matched "
          "filter allocates one day-length accumulator row per template)")
    A("")
    A("## 4. Growth (stack → template → matched filter, pass 0)")
    A("")
    A("Pass 0 only: the north measured that iterating **loses** families (50 → 47 → 47 → 40, "
      "ledger §7) because an enriching template absorbs its neighbours' events, and the null must "
      "run the identical chain.")
    A("")
    A("| | real chain | site-shift null |")
    A("|---|---|---|")
    A(f"| groups grown | {fmt(sp.get('n_groups_grown'))} | {fmt(nsp.get('n_groups_grown'))} |")
    A(f"| templates | {fmt(vs.get('n_templates'))} | {fmt(nvs.get('n_templates'))} |")
    A(f"| matched-filter detections | {fmt(vs.get('n_det_total'))} | "
      f"{fmt(nvs.get('n_det_total'))} |")
    A(f"| located (≥5 stations) | {fmt(vs.get('n_located'))} | {fmt(nvs.get('n_located'))} |")
    A(f"| shuffled-station control passes | {fmt(vs.get('n_ctrl_pass'))} | "
      f"{fmt(nvs.get('n_ctrl_pass'))} |")
    A(f"| **CANDIDATE FAMILIES** | **{fmt(vs.get('n_candidate'))}** | "
      f"**{fmt(nvs.get('n_candidate'))}** |")
    A("")
    A("## 5. Verification without truth")
    A("")
    A("A grown template is a **CANDIDATE FAMILY** only if all of:")
    A("")
    A("1. it locates on **≥5 stations** (stack-and-pick S, ≥4 MAD stack envelope peak) on the "
      "Crosson 1-D S travel times with the **intraslab-earthquake S station terms** of "
      "`results/x_station_corrections.csv`; a station without a usable term gets 0.0 and is "
      "counted, never imputed;")
    A("2. the **shuffled-station control passes** — the observed demeaned-residual rms beats the "
      "5th percentile of 50 station-label permutations of its own picks;")
    A("3. **|Δz| ≤ 8 km** of the Slab2 interface under its epicentre;")
    A("4. **≥3 distinct days** of detections;")
    A("5. **tremor co-timing above the same-day base rate**: the fraction of its detections "
      "within ±30 min of a PNSN epicentre ≤60 km exceeds the fraction obtained by redrawing the "
      "same number of times uniformly on the same days (200 realisations).")
    A("")
    A("**(c) The whole-chain site-shift null**, run once: each station's day series is replaced by "
      "the same station's data from a station-specific whole-day offset, restamped with the "
      "nominal day's `t0`, and the **identical** autocorrelation → grouping → stacking → matched "
      "filter → locator → screen is applied.  Inter-station coherence is destroyed; everything "
      "else — noise, glitches, diurnal structure, station geometry, thresholds — is preserved.")
    A("")
    A(f"> **{fmt(vs.get('n_candidate'))} candidate families on the real network, "
      f"{fmt(nvs.get('n_candidate'))} on the null.**")
    A("")
    if vs:
        A(f"Fraction of grown detections inside a PNSN tremor-active hour: "
          f"**{fmt(vs.get('frac_det_in_tremor_hour'), 3)}** against a base rate of "
          f"{fmt(vs.get('tremor_hour_base_rate'), 3)} (null: "
          f"{fmt(nvs.get('frac_det_in_tremor_hour'), 3)} against "
          f"{fmt(nvs.get('tremor_hour_base_rate'), 3)}).  The matched filter is run over the "
          "tremor days only, so this base rate is high by construction and the test is weak.")
        A("")
    if len(cand):
        A("### Candidate families")
        A("")
        A("| grp | lat | lon | depth km | Slab2 km | Δz km | stations | n_det | days | "
          "co-timing | base | rms s |")
        A("|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in cand.sort_values("n_det", ascending=False).itertuples():
            A(f"| {r.grp} | {r.lat:.2f} | {r.lon:.2f} | {r.dep:.0f} | {r.slab2:.1f} | "
              f"{r.dz:+.1f} | {r.n_sta} | {r.n_det} | {r.n_days} | {r.co_timing:.3f} | "
              f"{r.co_base:.3f} | {r.rms:.2f} |")
        A("")
    else:
        A("No grown template passed all five clauses.  The per-template screen is in "
          f"`{os.path.relpath(f'{R}/verify_it0.csv', STUDY)}`.")
        A("")
    if len(ncand):
        A(f"The null produced {len(ncand)} \"family\" (families) passing the same screen; they are "
          f"listed in `{os.path.relpath(f'{N}/verify_it0.csv', STUDY)}`.")
        A("")
    A("### (d) Chestler comparison — descriptive only")
    A("")
    A("`../codex-LFE/results/olympic_reference_discovery_v1/"
      "derived_coordinate_rows_not_a_catalog.csv` is, as its filename says, **not a catalogue**; "
      "it is used here only to say where the candidates sit relative to previously published "
      "Olympic LFE coordinates, never as truth.")
    A("")
    ch = vs.get("chestler")
    if ch:
        A(f"* candidates: median **{fmt(ch.get('median_km'), 1)} km** to the nearest of "
          f"{fmt(ch.get('n_ref'))} derived coordinates (closest {fmt(ch.get('min_km'), 1)} km);")
        A(f"* 1,000 equal-sized random draws from the sub-region's active cells: median "
          f"{fmt(ch.get('random_cell_median_km'), 1)} km (5th percentile "
          f"{fmt(ch.get('random_cell_p05_km'), 1)} km);")
        A(f"* fraction of random draws at least as close: **{fmt(ch.get('p_better_than_random'), 3)}**.")
    else:
        A("Not computed (no candidate family to compare).")
    A("")
    A("## 6. Files")
    A("")
    A("| | |")
    A("|---|---|")
    A(f"| pair stage | `{os.path.relpath(R, STUDY)}/auto/`, `pair_all.parquet`, "
      "`pair_pairstats.json` |")
    A(f"| quiet control | `{os.path.relpath(R, STUDY)}/auto_quiet/`, `pair_quiet.parquet` |")
    A(f"| groups | `{os.path.relpath(R, STUDY)}/gb_nodes.csv`, `gb_summary.csv` |")
    A(f"| templates + detections | `{os.path.relpath(R, STUDY)}/grow/tpl_it0.npz`, "
      "`det_it0.parquet` |")
    A(f"| locations | `{os.path.relpath(R, STUDY)}/loc_it0.csv` |")
    A(f"| screen | `{os.path.relpath(R, STUDY)}/verify_it0.csv` |")
    A(f"| site-shift null (whole chain) | `{os.path.relpath(N, STUDY)}/` "
      f"(offsets in `null_offsets.json`) |")
    A(f"| engine parameters as run | `{os.path.relpath(R, STUDY)}/auto/_params.json` |")
    if par:
        A("")
        A("```json")
        A(json.dumps({k: v for k, v in par.items() if k != "stations"}, indent=1))
        A("```")
    if off:
        A("")
        A(f"Null day offsets (station → whole days): {len(off)} stations, "
          f"range {min(off.values())}–{max(off.values())} days, seed "
          f"{os.environ.get('L_NULL_SEED', '11')}.")
    A("")
    open(a.out, "w").write("\n".join(o) + "\n")
    print(f"wrote {a.out} ({len(o)} lines)")

    note = [
        "", "---", "",
        "## 2026-09-20 — agent L: Bostock's recipe on the Washington N-Olympic 2010 array bed",
        "",
        f"One detached queue (`scripts/l_wa2010_queue.sh`).  Bed: {bed.get('n_stations')} "
        f"probe-admitted 3-C stations ({', '.join(f'{k} {v}' for k, v in sorted(bed.get('by_net', {}).items()))}), "
        f"{bed.get('day0')}…{bed.get('day1')}, {len(hrs)} PNSN tremor hours, "
        f"{ps.get('n_blocks_quiet', 0)} hour-matched quiet hours.",
        "",
        f"* pair stage: {fmt(ps.get('n_det_tremor'))} detections above 8 × MAD "
        f"({fmt(ps.get('nmad_in_sd_median'), 2)} sd), tremor/quiet pair-count ratio "
        f"**{fmt(ps.get('tremor_quiet_ratio'), 2)}×**;",
        f"* grouping: {fmt(gm.get('n_kept'))} components with ≥3 members; grown "
        f"{fmt(sp.get('n_groups_grown'))};",
        f"* **{fmt(vs.get('n_candidate'))} candidate families** against "
        f"**{fmt(nvs.get('n_candidate'))}** from the whole-chain site-shift null;",
        f"* no catalogue truth exists here — nothing is a rediscovery count.  "
        f"Full ledger `results/WA_2010_RESULTS.md`.",
        "",
        "Engine reused unedited (`k_auto`/`k_group`/`k_grow`/`k_locate`); the WA bed is supplied "
        "by `l_common.py` (station table from `channels_2006_2012.csv`, PB borehole azimuths added "
        "to `s_common.AZ12` so the boreholes are rotated instead of dropped, net-qualified cache "
        "index because C8.PA02… shadows XG.PA02…, station-hour health mask, day-shift null).",
    ]
    if a.no_notes:
        print("[dry] NOTES.md not touched")
        return
    with open(f"{STUDY}/NOTES.md", "a") as f:
        f.write("\n".join(note) + "\n")
    print("appended NOTES.md")


if __name__ == "__main__":
    main()
