#!/usr/bin/env python
"""WA pre-production: write the W2 (station health) section and the W0 re-hash line into results/wa_pre/WA_PRE_RESULTS.md,
between the markers <!-- W2-AUTO-START --> / <!-- W2-AUTO-END --> (the rest of the file is hand-written and kept).
Also writes results/wa_pre/w2/W2_RESULTS.md (the same section, standalone).
"""
import os, sys, glob, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pre_common as P

O2 = f"{P.PRE}/w2"
MD = f"{P.PRE}/WA_PRE_RESULTS.md"


def fmt(x, n=3):
    try:
        return "nan" if x is None or not np.isfinite(float(x)) else f"{float(x):.{n}f}"
    except Exception:
        return str(x)


def w1_final():
    """W1 readability after the FULL reads: W1's quick check + W2 raw (full read) + SVI R2 raw status (full read)."""
    X = pd.read_parquet(f"{P.PRE}/w1/w1_station_days.parquet")
    X = X[X.expected]
    rs = pd.read_csv(f"{O2}/raw_status.csv")[["date", "netsta", "status"]].rename(columns={"status": "w2_full_read"})
    sv = []
    for f in glob.glob(f"{P.ROOT}/lfe_extension/results/refine/r2/jobs/*_raw.csv"):
        sv.append(pd.read_csv(f, dtype=str))
    sv = pd.concat(sv)[["date", "netsta", "status"]].rename(columns={"status": "svi_full_read"}) if sv else pd.DataFrame(columns=["date", "netsta", "svi_full_read"])
    X = X.merge(rs, on=["date", "netsta"], how="left").merge(sv.drop_duplicates(["date", "netsta"]), on=["date", "netsta"], how="left")
    X["full_read"] = X.w2_full_read.fillna(X.svi_full_read)
    X["full_ok"] = X.full_read.eq("ok")
    pres = X[X.raw]
    n_err = int((pres.full_read == "error").sum()); n_unk = int(pres.full_read.isna().sum())
    X[["netsta", "date", "raw", "readable", "full_read", "full_ok", "cache", "cache_ok"]].to_csv(f"{P.PRE}/w1/w1_readability_final.csv", index=False)
    return (f"- W1 final readability (full reads by W2 raw + SVI R2): of {len(pres)} present raw station-days, {int(pres.full_ok.sum())} read in full "
            f"without error, {n_err} errors, {n_unk} not fully read by either (`w1/w1_readability_final.csv`).")


def main():
    L = []
    M = pd.read_csv(f"{O2}/r2_health_monthly.csv")
    SW = pd.read_csv(f"{O2}/r2_swap_checks.csv")
    E = pd.read_csv(f"{O2}/r2_effective_sets.csv")
    X = pd.read_csv(f"{O2}/r2_exclusions.csv")
    B = pd.read_csv(f"{O2}/r2_bad_for_dvv.csv")
    pc = json.load(open(f"{O2}/r2_positive_control.json"))
    rs = pd.read_csv(f"{O2}/raw_status.csv")
    cm = pd.read_csv(f"{O2}/cache_manifest.csv") if os.path.exists(f"{O2}/cache_manifest.csv") else pd.DataFrame(dict(status=[]))
    T = M[M.tested]
    L.append("## Step 2 (W2) — station health of the 21 core stations, 2007 → 2026-09-15 (auto-written by scripts/pre_report.py)\n")
    L.append(f"- Raw metrics: {len(rs)} station-days read by this run ({(rs.status=='ok').sum()} ok, {(rs.status=='error').sum()} errors, "
             f"{(rs.status=='missing').sum()} expected but no file); the 13 SVI-shared stations 2014-01-01 → 2026-09-15 come from SVI R2's "
             f"daily table (same files, same code). Files: `w2/r2_health_daily.parquet`, `w2/r2_health_stationday.parquet`.")
    L.append(f"- Stack test: {len(pd.read_csv(f'{O2}/stack_days.csv'))} strongest-tremor days ({M.month.str[:4].min()}–{M.month.str[:4].max()}; "
             f"sealed 2007-01 / 2009-05 and 2011-04 never used), cache {cm.status.value_counts().to_dict()}; "
             f"{len(T)} station-months tested (≥ 3 families with ≥ 5 confident repeats, ≥ 30 repeats). "
             "Confident repeats for a station = detections whose network WITHOUT that station's site still reaches 8.75 × MAD "
             "with ≥ 2 other sites at site z ≥ 3. Scrambled control = the same station read on other selected days at the same time of day.")
    st = M.status_eff.value_counts().to_dict()
    L.append(f"- Station-month status (all {len(M)} metadata-alive station-months): {st}.")
    L.append(f"- Tested months at control level: {int((T.status=='control_level').sum())} of {len(T)}. Excluded station-months "
             f"(health reason, not plain no-data): {len(X)} (`w2/r2_exclusions.csv`).")
    per = T.groupby("netsta").agg(n_tested=("month", "size"), n_ctl=("status", lambda v: int((v == "control_level").sum())),
                                  cc_med=("cc_place", "median"), ctl95_med=("ctl_cc_place_p95", "median")).reset_index()
    L.append("\n| station | tested months | at control level | median stack cc at placement | median control p95 |\n|---|---|---|---|---|")
    for r in per.itertuples():
        L.append(f"| {r.netsta} | {r.n_tested} | {r.n_ctl} | {fmt(r.cc_med)} | {fmt(r.ctl95_med)} |")
    L.append("\nSwap / epoch-change checks (`w2/r2_swap_checks.csv`; carry-over = the old node's template used on the new node):\n")
    L.append("| station | date | change | cc before / ctl p95 | cc after / ctl p95 | best channel transform after | verdict |\n|---|---|---|---|---|---|---|")
    for r in SW.itertuples():
        L.append(f"| {r.netsta} | {r.date} | {str(r.change)[:70]} | {fmt(getattr(r,'cc_before',np.nan))} / {fmt(getattr(r,'ctl95_before',np.nan))} | "
                 f"{fmt(r.cc_after)} / {fmt(r.ctl95_after)} | {str(getattr(r,'perm_after',''))[:40]} | {r.verdict} |")
    L.append(f"\n- Positive control (LRIV after 2019-05-30, HH node, old BH template carried): tested months {pc['n_tested_post']}, "
             f"at control level {pc['n_control_level_post']}, median cc {fmt(pc['cc_place_post_median'])} vs control p95 "
             f"{fmt(pc['ctl_p95_post_median'])}; before 2019-05: {pc['n_tested_pre']} tested, {pc['n_control_level_pre']} at control, "
             f"median cc {fmt(pc['cc_place_pre_median'])}. Flagged: {pc['flags_lriv_post2019']}. ({pc['svi_reference']})")
    L.append(f"- Effective station set per month: `w2/r2_effective_sets.csv` (median {int(E.n_effective.median())} of 21, "
             f"min {int(E.n_effective.min())} in {E.loc[E.n_effective.idxmin(), 'month']}).")
    L.append(f"- Bad-for-dv/v windows (KMR style, `w2/r2_bad_for_dvv.csv`): {len(B)} windows; by reason {B.reason.str.split(' ').str[0].value_counts().to_dict()}.")
    L.append("- Read with care: the stack test has little power where a station's templates are weak for the families active "
             "on the test days (e.g. SW carriers B017/B018 on NW-dominated days); a 'control_level' month there means 'no "
             "evidence the station sees these families', not 'broken'. Months between test days inherit the nearest tested month "
             "of the same instrument epoch.")
    vf = sorted(glob.glob(f"{P.PRE}/w0/verify_*.csv"))
    if vf:
        V = pd.read_csv(vf[-1])
        L.append(f"- W0 re-hash at the end: {len(V)} frozen input files, changed {int((~V.same).sum())} (`{os.path.relpath(vf[-1], P.PRE)}`).")
    try:
        L.append(w1_final())
    except Exception as e:
        L.append(f"- (W1 final readability table not written: {type(e).__name__}: {e})")
    sec = "\n".join(L) + "\n"
    open(f"{O2}/W2_RESULTS.md", "w").write(sec)
    if os.path.exists(MD):
        s = open(MD).read()
        a, b = "<!-- W2-AUTO-START -->", "<!-- W2-AUTO-END -->"
        if a in s and b in s:
            s = s[:s.index(a) + len(a)] + "\n" + sec + s[s.index(b):]
        else:
            s += f"\n{a}\n{sec}{b}\n"
        open(MD, "w").write(s)
    print(sec)


if __name__ == "__main__":
    main()
