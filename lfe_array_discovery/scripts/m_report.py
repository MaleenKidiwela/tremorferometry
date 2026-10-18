#!/usr/bin/env python
"""The 2011 re-fire report: apply the PRE-STATED confirmation rule and write it up.  Prefix `m_`.

The rule was fixed BEFORE the seal was opened (agent M's brief, recorded in NOTES.md):

    a 2010 candidate is CONFIRMED if, in 2011,
      (i)   it has >= 3 distinct tremor days with detections,
      (ii)  its quiet-day detection rate is <= 0.1 per template-day,
      (iii) its 2011 location is INTERIOR to the wider search grid and within 10 km of its
            2010 location.

It is applied mechanically and identically to the 50 real and the 17 null candidate templates;
the null count IS the false-confirmation rate.

  python scripts/m_report.py --out results/WA_2011_REFIRE_RESULTS.md
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
STUDY = os.path.dirname(HERE)
RES = f"{STUDY}/results"
sys.path.insert(0, HERE)
sys.path.insert(0, "/home/jovyan/tremorferometry/lfe_picker/scripts")
import m_common as M                                                        # noqa: E402

MIN_DAYS_2011 = 3
MAX_QUIET_RATE = 0.1
MAX_MOVE_KM = 10.0
INTERIOR_2010 = [1656, 1193, 1152, 429]     # WA_2010_RESULTS.md addendum


def km3(la1, lo1, dp1, la2, lo2, dp2):
    import s_common as S
    h = S.epi_km(la1, lo1, la2, lo2)
    return float(np.sqrt(h ** 2 + (dp1 - dp2) ** 2)), float(h)


def load(tag, out):
    d = {}
    for k, p in [("loc2010", f"{out}/loc2010_{tag}.csv"),
                 ("loc2011", f"{out}/loc2011_{tag}.csv"),
                 ("tplsta", f"{out}/tpl_stations_{tag}_tremor.csv")]:
        d[k] = pd.read_csv(p) if os.path.exists(p) else pd.DataFrame()
    for k, p in [("tremor", f"{out}/det_{tag}_tremor.parquet"),
                 ("fa", f"{out}/det_{tag}_fa.parquet")]:
        d[k] = pd.read_parquet(p) if os.path.exists(p) else pd.DataFrame(
            columns=["grp", "t", "cc", "nmad", "date"])
    return d


def table(tag, grps, out, quiet_days, n_quiet, orig):
    d = load(tag, out)
    L10 = d["loc2010"].set_index("grp") if len(d["loc2010"]) else pd.DataFrame()
    L11 = d["loc2011"].set_index("grp") if len(d["loc2011"]) else pd.DataFrame()
    TS = d["tplsta"].set_index("grp") if len(d["tplsta"]) else pd.DataFrame()
    DT, DF = d["tremor"], d["fa"]
    O = orig.set_index("grp")
    nq = DF[DF.date.isin(quiet_days)] if len(DF) else DF
    rows = []
    for g in grps:
        dt = DT[DT.grp == g]
        nquiet = int((nq.grp == g).sum()) if len(nq) else 0
        nfa = int((DF.grp == g).sum()) if len(DF) else 0
        r = dict(grp=g,
                 n_sta_tpl_2010=float(TS.n_sta_tpl_2010.get(g, np.nan)) if len(TS) else np.nan,
                 n_sta_tpl_2011=float(TS.n_sta_tpl_2011.get(g, np.nan)) if len(TS) else np.nan,
                 orig_lat=float(O.lat.get(g, np.nan)), orig_lon=float(O.lon.get(g, np.nan)),
                 orig_dep=float(O.dep.get(g, np.nan)),
                 n_det_2010=int(O.n_det.get(g, 0)), n_days_2010=int(O.n_days.get(g, 0)),
                 n_det_2011=int(len(dt)), n_days_2011=int(dt.date.nunique()),
                 n_det_fa=nfa, n_det_quiet=nquiet,
                 quiet_rate=nquiet / max(n_quiet, 1))
        for yr, LL in (("2010", L10), ("2011", L11)):
            ok = len(LL) and g in LL.index and bool(LL.located.get(g, False))
            r[f"loc{yr}"] = bool(ok)
            for c in ("lat", "lon", "dep", "rms", "ctrl_pass", "edge", "n_sta", "dz"):
                r[f"{c}_{yr}"] = (LL[c].get(g, np.nan) if ok and c in LL.columns else np.nan)
        if r["loc2010"] and r["loc2011"]:
            r["move_km"], r["move_epi_km"] = km3(r["lat_2010"], r["lon_2010"], r["dep_2010"],
                                                 r["lat_2011"], r["lon_2011"], r["dep_2011"])
        else:
            r["move_km"] = r["move_epi_km"] = np.nan
        r["c_days"] = bool(r["n_days_2011"] >= MIN_DAYS_2011)
        r["c_quiet"] = bool(r["quiet_rate"] <= MAX_QUIET_RATE)
        r["c_loc"] = bool(r["loc2011"] and r["loc2010"]
                          and not bool(r["edge_2011"]) and np.isfinite(r["move_km"])
                          and r["move_km"] <= MAX_MOVE_KM)
        r["confirmed"] = bool(r["c_days"] and r["c_quiet"] and r["c_loc"])
        r["confirmed_epi"] = bool(r["c_days"] and r["c_quiet"] and r["loc2011"] and r["loc2010"]
                                  and not bool(r["edge_2011"])
                                  and np.isfinite(r["move_epi_km"])
                                  and r["move_epi_km"] <= MAX_MOVE_KM)
        rows.append(r)
    return pd.DataFrame(rows).sort_values("grp").reset_index(drop=True)


def md_table(T, cols, hdr):
    out = ["| " + " | ".join(hdr) + " |", "|" + "---|" * len(hdr)]
    for r in T.itertuples():
        cells = []
        for c in cols:
            v = getattr(r, c)
            if isinstance(v, (bool, np.bool_)):
                cells.append("**yes**" if v else "no")
            elif isinstance(v, (float, np.floating)):
                cells.append("--" if not np.isfinite(v) else f"{v:.2f}")
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def chestler(T, label):
    import s_common as S
    if not os.path.exists(M.L.CHESTLER) or not len(T):
        return None
    C = pd.read_csv(M.L.CHESTLER)
    cla, clo = C.latitude.to_numpy(), C.longitude.to_numpy()
    d = np.array([np.min(S.epi_km(float(x), float(y), cla, clo))
                  for x, y in zip(T.lat_2011, T.lon_2011)])
    AC = M.L.active_cells()
    rr = np.random.default_rng(5)
    draws = []
    for _ in range(1000):
        i = rr.integers(0, len(AC), len(T))
        draws.append(np.median([np.min(S.epi_km(float(AC.lat.iloc[j]), float(AC.lon.iloc[j]),
                                                cla, clo)) for j in i]))
    draws = np.array(draws)
    return dict(label=label, n=int(len(T)), n_ref=int(len(C)),
                median_km=float(np.median(d)), min_km=float(d.min()),
                per_family_km=[float(x) for x in d],
                random_cell_median_km=float(np.median(draws)),
                random_cell_p05_km=float(np.quantile(draws, .05)),
                p_random_at_least_as_close=float((draws <= np.median(d)).mean()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{RES}/WA_2011_REFIRE_RESULTS.md")
    a = ap.parse_args()
    out = M.OUT
    bed = json.load(open(f"{out}/bed.json"))
    quiet_days = bed["quiet_dates"]
    n_quiet = len(quiet_days)
    tremor_days = bed["tremor_dates"]

    Vr = pd.read_csv(f"{RES}/l_wa2010/verify_it0.csv")
    Vn = pd.read_csv(f"{RES}/l_wa2010_null/verify_it0.csv")
    gr = sorted(int(x) for x in Vr[Vr.candidate].grp)
    gn = sorted(int(x) for x in Vn[Vn.candidate].grp)

    Tr = table("real", gr, out, quiet_days, n_quiet, Vr)
    Tn = table("null", gn, out, quiet_days, n_quiet, Vn)
    Tr.to_csv(f"{out}/confirm_real.csv", index=False)
    Tn.to_csv(f"{out}/confirm_null.csv", index=False)

    seal = sorted(glob.glob(f"{RES}/SEALED_2011_OPENED_*"))
    seal_txt = open(seal[0]).read().strip() if seal else "(marker missing)"
    nconf_r, nconf_n = int(Tr.confirmed.sum()), int(Tn.confirmed.sum())
    Cr = Tr[Tr.confirmed]
    ch = chestler(Cr, "confirmed real families")

    interior = Tr[Tr.grp.isin(INTERIOR_2010)]
    edgy = Tr[~Tr.grp.isin(INTERIOR_2010)]

    S = []
    A = S.append
    A("# Washington N-Olympic — CONFIRMATION by re-fire on the SEALED 2011 window\n")
    A("Agent M, 2026-09-21.  The 2010 candidate templates are re-fired **unchanged** on data that "
      "no part of this project had read before: the sealed 2011-07-20…09-10 N-Olympic set.  "
      "Nothing was tuned, re-cut or re-thresholded — `k_grow.matched_filter` (network mean of the "
      "per-station normalised 3-C cc, threshold = median + 8 × MAD of that day's own function, "
      "4 s decluster) reads the template waveforms straight out of the 2010 run's "
      "`grow/tpl_it0.npz`.  The 17 site-shift-null candidate templates are re-fired on the same "
      "2011 data by the same code: **their confirmation count is the false-confirmation rate.**\n")
    A("## 0. The seal\n```\n" + seal_txt + "\n```\n")
    A("## 1. The 2011 bed\n")
    A(f"* window **{bed['window'][0]} … {bed['window'][1]}**, sub-region N-Olympic, tranche "
      f"{bed['tranche']}.\n"
      f"* **{bed['n_stations']} three-component stations** = the 2010 probe-admitted, "
      f"2010-health-kept set ({bed['n_kept_2010']}) minus {bed['n_absent_2011']} that do not "
      f"exist in 2011 ({', '.join(bed['absent_2011'])}) minus {bed['n_drop_2011_health']} that "
      f"fail the same > 30 % vetoed-station-hour bar on the 2011 census"
      + (f" ({', '.join(bed['drop_2011_health'])})" if bed['drop_2011_health'] else "") + ".\n"
      f"* by network: " + ", ".join(f"{k} {v}" for k, v in bed["by_net"].items()) + ".\n"
      f"* **{bed['n_cached_station_days']} station-days cached** with the frozen 20 Hz / 1–8 Hz "
      f"chain (rational twin for the 50 Hz XG/XU sites); median "
      f"{bed['cached_stations_on_tremor_days_median']} stations per tremor day.\n"
      f"* PNSN candidate (tremor) hours **{bed['n_candidate_hours']}** on "
      f"**{bed['n_tremor_days']}** days ({tremor_days[0]} … {tremor_days[-1]}); the same code on "
      f"the 2010 window reproduces that run's **{bed['selfcheck_2010_hours']} hours / "
      f"{bed['selfcheck_2010_days']} days** exactly.\n"
      f"* strictly-quiet days (no PNSN tremor within 100 km of the region centre, ≥ 3 days from "
      f"any such day): **{n_quiet}** ({', '.join(quiet_days)}).\n")
    A("> **Documented deviation (bed-forced).**  The brief asks for 20 quiet days of 2011.  The "
      "sealed set is only the 53-day window around the 2011 ETS and it contains exactly "
      f"{n_quiet} days that satisfy the frozen quiet rule; no more data may be downloaded.  The "
      f"false-alarm reference is therefore run over **all {bed['n_nontremor_days']} non-tremor "
      f"days** of the window and the pre-stated quiet-day rate is evaluated on the "
      f"{n_quiet} strictly-quiet ones.  Both are reported.  With {n_quiet} quiet days the "
      f"≤ 0.1 / template-day bar is met only by a template with **zero** quiet-day detections.\n")

    A("## 2. Station survival of the templates\n")
    for tag, T in (("real", Tr), ("null", Tn)):
        if len(T) and np.isfinite(T.n_sta_tpl_2011).any():
            A(f"* {tag}: template stations present in 2011, median "
              f"{T.n_sta_tpl_2011.median():.0f} of {T.n_sta_tpl_2010.median():.0f} "
              f"(min {T.n_sta_tpl_2011.min():.0f}, worst fraction "
              f"{(T.n_sta_tpl_2011 / T.n_sta_tpl_2010).min():.2f}).")
    A("")

    A("## 3. The re-fire\n")
    for tag, T, nd in (("real (50 candidates)", Tr, len(gr)), ("null (17 candidates)", Tn, len(gn))):
        A(f"* **{tag}**: {int(T.n_det_2011.sum())} detections on the "
          f"{bed['n_tremor_days']} tremor days; {int((T.n_det_2011 > 0).sum())}/{nd} templates "
          f"fire at all, {int((T.n_days_2011 >= 3).sum())} on ≥ 3 distinct days; "
          f"{int(T.n_det_fa.sum())} detections on the {bed['n_nontremor_days']} non-tremor days "
          f"({T.n_det_fa.sum() / max(nd * bed['n_nontremor_days'], 1):.3f} per template-day), "
          f"{int(T.n_det_quiet.sum())} on the {n_quiet} strictly-quiet days "
          f"({T.n_det_quiet.sum() / max(nd * n_quiet, 1):.3f} per template-day).")
    A("")

    A("## 4. Re-location on the wider grid\n")
    A("Grid **46.8–49.0 N / −125.2…−121.8 / 15–60 km** at the 2010 spacing (0.03°, 0.04°, 2 km) "
      "— the 2010 box was 47.50–48.80 / −124.80…−122.20 / 20–58 km, and 46 of the 50 real "
      "candidates railed on one of its faces.  Both years are re-stacked from their own "
      "detections with `k_grow.build_templates` (unchanged) and located with `k_locate`'s picks, "
      "demeaned-residual rms and 50-realisation shuffled-station control.\n")
    for tag, T in (("real", Tr), ("null", Tn)):
        n10, n11 = int(T.loc2010.sum()), int(T.loc2011.sum())
        e10 = int(T.edge_2010.fillna(False).astype(bool).sum())
        e11 = int(T.edge_2011.fillna(False).astype(bool).sum())
        A(f"* {tag}: located from 2010 detections {n10}/{len(T)} ({e10} still on an edge), "
          f"from 2011 detections {n11}/{len(T)} ({e11} on an edge); "
          f"median 3-D move {T.move_km.median():.1f} km "
          f"(epicentral {T.move_epi_km.median():.1f} km).")
    A("")

    A("## 5. CONFIRMATION (rule fixed before the seal was opened)\n")
    A(f"> ≥ {MIN_DAYS_2011} distinct 2011 tremor days with detections **and** quiet-day rate "
      f"≤ {MAX_QUIET_RATE}/template-day **and** a 2011 location that is interior to the wider "
      f"grid and within {MAX_MOVE_KM:.0f} km (3-D) of the 2010 location.\n")
    A(f"### **REAL {nconf_r} / {len(gr)} confirmed.  NULL {nconf_n} / {len(gn)} confirmed "
      f"(false-confirmation rate {nconf_n}/{len(gn)} = {nconf_n / max(len(gn), 1):.2f}).**\n")
    A(f"With the epicentral (horizontal-only) reading of the 10 km clause: real "
      f"{int(Tr.confirmed_epi.sum())}, null {int(Tn.confirmed_epi.sum())}.\n")
    A("Clause-by-clause (how many pass each on its own):\n")
    A("| clause | real | null |\n|---|---|---|")
    for c, lab in (("c_days", f"≥ {MIN_DAYS_2011} tremor days"),
                   ("c_quiet", f"quiet rate ≤ {MAX_QUIET_RATE}"),
                   ("c_loc", "interior + ≤ 10 km")):
        A(f"| {lab} | {int(Tr[c].sum())}/{len(gr)} | {int(Tn[c].sum())}/{len(gn)} |")
    A("")

    cols = ["grp", "orig_lat", "orig_lon", "orig_dep", "lat_2010", "lon_2010", "dep_2010",
            "lat_2011", "lon_2011", "dep_2011", "move_km", "edge_2011", "n_det_2011",
            "n_days_2011", "n_det_quiet", "quiet_rate", "confirmed"]
    hdr = ["grp", "2010 lat*", "lon*", "dep*", "2010 lat", "lon", "dep", "2011 lat", "lon",
           "dep", "move km", "2011 edge", "n_det 2011", "days 2011", "quiet det", "quiet FA/d",
           "CONFIRMED"]
    A("### 5a. The three interior 2010 candidates (the actual evidence of the 2010 run)\n")
    A("grp 1656, 1193 and the marginal 1152 are the interior, ≥5-day, co-timing-excess "
      "candidates of `WA_2010_RESULTS.md`'s addendum; 429 is the fourth interior one that failed "
      "the day/co-timing filter there.\n")
    A(md_table(interior, cols, hdr) + "\n")
    A("\\* `2010 lat/lon/dep` with a star = the ORIGINAL 2010 narrow-grid location from "
      "`loc_it0.csv`; the unstarred 2010 columns are the wide-grid re-location from the same "
      "2010 detections, which is the one the 10 km clause uses.\n")
    A("### 5b. The 46 edge-railing 2010 candidates\n")
    A(md_table(edgy, cols, hdr) + "\n")
    A("### 5c. The 17 site-shift-null candidates (false-confirmation control)\n")
    A(md_table(Tn, cols, hdr) + "\n")

    A("## 6. Descriptive — distance to the Chestler-derived coordinates\n")
    A("`../codex-LFE/results/olympic_reference_discovery_v1/"
      "derived_coordinate_rows_not_a_catalog.csv` is not a catalogue and is used only to say "
      "where the confirmed families sit relative to previously published Olympic LFE "
      "coordinates.\n")
    if ch:
        A(f"* confirmed families (n={ch['n']}): median **{ch['median_km']:.1f} km** to the "
          f"nearest of {ch['n_ref']} derived coordinates (closest {ch['min_km']:.1f} km); "
          f"per family " + ", ".join(f"{x:.1f}" for x in ch["per_family_km"]) + " km.\n"
          f"* 1,000 equal-sized random draws from the sub-region's active cells: median "
          f"{ch['random_cell_median_km']:.1f} km (5th pct "
          f"{ch['random_cell_p05_km']:.1f} km).\n"
          f"* fraction of random draws at least as close: "
          f"**{ch['p_random_at_least_as_close']:.3f}**.\n")
        json.dump(ch, open(f"{out}/chestler_confirmed.json", "w"), indent=1)
    else:
        A("* no confirmed family to compare.\n")

    A("## 7. Files\n")
    A("| | |\n|---|---|")
    for lab, p in [("bed / stations / hours / quiet days", f"{out}/bed.json, stations.csv, "
                    "candidate_hours.csv, quiet_days.csv"),
                   ("2011 re-fire detections (tremor days)", f"{out}/det_real_tremor.parquet, "
                    "det_null_tremor.parquet"),
                   ("2011 false-alarm detections (non-tremor days)", f"{out}/det_real_fa.parquet, "
                    "det_null_fa.parquet"),
                   ("template station survival", f"{out}/tpl_stations_*.csv"),
                   ("wide-grid re-locations", f"{out}/loc2010_*.csv, loc2011_*.csv "
                    "(+ _tpl.npz re-stacks)"),
                   ("confirmation tables", f"{out}/confirm_real.csv, confirm_null.csv"),
                   ("seal", "results/sealed_2011_manifest.sha256, SEALED_2011_OPENED_*")]:
        A(f"| {lab} | `{p}` |")
    A("")
    A("## 8. What this does and does not settle\n")
    A("* The re-fire is a genuine out-of-sample test: the 2011 data were hashed and sealed "
      "before the 2010 chain was written, and the templates crossing into 2011 were not "
      "re-fitted in any way.\n"
      "* The null templates are the controlled comparison for every number above — they were "
      "built by the identical chain on a site-shifted network and are fired by the identical "
      "code on the identical 2011 data.\n"
      "* There is still no Washington LFE catalogue, so nothing here is a rediscovery count and "
      "no family is claimed to be a catalogued LFE family.\n"
      f"* The quiet-day clause rests on {n_quiet} days (see the deviation in §1); the non-tremor "
      "false-alarm rate over the whole window is the more precise, weaker-definition version.\n")

    open(a.out, "w").write("\n".join(S) + "\n")
    print(f"[report] -> {a.out}")
    print(f"[report] REAL confirmed {nconf_r}/{len(gr)}; NULL confirmed {nconf_n}/{len(gn)}")
    print(Tr[Tr.confirmed][["grp", "lat_2011", "lon_2011", "dep_2011", "move_km", "n_det_2011",
                            "n_days_2011", "n_det_quiet"]].to_string(index=False))

    note = (f"\n## 2026-09-21 — agent M: sealed-2011 re-fire (confirmation of the WA 2010 "
            f"candidates)\n"
            f"Seal opened and recorded ({os.path.basename(seal[0]) if seal else 'MARKER MISSING'}"
            f"); {bed['n_cached_station_days']} station-days cached; {bed['n_stations']} stations "
            f"(of the 96 of 2010); {bed['n_tremor_days']} tremor days, {n_quiet} strictly-quiet "
            f"days (bed-forced deviation from the requested 20 — documented in "
            f"`results/WA_2011_REFIRE_RESULTS.md` §1).  Unchanged templates, unchanged 8×MAD "
            f"matched filter, wider relocation grid 46.8–49.0 / −125.2…−121.8 / 15–60 km.\n"
            f"**Result: REAL {nconf_r}/{len(gr)} confirmed, NULL {nconf_n}/{len(gn)} confirmed "
            f"(false-confirmation rate {nconf_n / max(len(gn), 1):.2f}).**  Interior 2010 "
            f"candidates: " + ", ".join(
                f"grp {int(r.grp)} {'CONFIRMED' if r.confirmed else 'not confirmed'}"
                for r in interior.itertuples()) + ".\n"
            f"Full write-up `results/WA_2011_REFIRE_RESULTS.md`; tables "
            f"`results/m_2011/confirm_{{real,null}}.csv`.\n")
    with open(f"{STUDY}/NOTES.md", "a") as f:
        f.write(note)
    print("[report] appended to NOTES.md")


if __name__ == "__main__":
    main()
