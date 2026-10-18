#!/usr/bin/env python
"""REPORT for the all-days detection pass (ab_alldays.py): real vs site-shift-null families.

Nothing is filtered: every number is a column, the null gets the same table, and thresholds on
nmad appear only as a separation CURVE (memory `rank-dont-discard`).  Counts are also given as
RATES per scanned day in the same timing bin, because the number of scanned days differs between
ETS, near-ETS and far-from-ETS periods.

  python scripts/ab_alldays_report.py
-> results/ab_royer/alldays/REPORT.md
-> results/ab_royer/alldays/<SUB>/family_month_<kind>.csv   detections per family per month
-> results/ab_royer/alldays/<SUB>/timing_<kind>.csv          timing bins, counts and rates
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import ab_alldays as AD                                                     # noqa: E402

OUT = AD.OUT
KINDS = ["real", "null1", "null2", "null3"]
BINS = [(-np.inf, -90, "pre >90 d"), (-90, -30, "pre 30-90 d"), (-30, -7, "pre 7-30 d"),
        (-7, 0, "pre 0-7 d"), (0, 0, "during"), (0, 7, "post 0-7 d"), (7, 30, "post 7-30 d"),
        (30, 90, "post 30-90 d"), (90, np.inf, "post >90 d")]
NMAD_CURVE = [8, 9, 10, 12, 15, 20, 30]


def tbin(phase, rel):
    out = np.full(len(rel), "unknown", dtype=object)
    rel = np.asarray(rel, float)
    phase = np.asarray(phase, object)
    out[phase == "during"] = "during"
    for lo, hi, name in BINS:
        if name == "during":
            continue
        if lo < 0:
            m = (phase == "pre") & (rel > lo) & (rel <= hi) if np.isfinite(lo) else \
                (phase == "pre") & (rel <= hi)
        else:
            m = (phase == "post") & (rel > lo) & (rel <= hi) if np.isfinite(hi) else \
                (phase == "post") & (rel > lo)
        out[m] = name
    return out


def fmt(D, floatfmt="{:.3g}"):
    D = D.copy()
    for c in D.columns:
        if D[c].dtype.kind == "f":
            D[c] = D[c].map(lambda v: "" if pd.isna(v) else floatfmt.format(v))
    cols = [str(c) for c in D.columns]
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for r in D.astype(str).itertuples(index=False):
        lines.append("| " + " | ".join(r) + " |")
    return "\n".join(lines)


def load(sub, kind):
    p = f"{OUT}/{sub}/detections_{kind}.parquet"
    dp = f"{OUT}/{sub}/{kind}/days.csv"
    if not os.path.exists(p):
        return None, None
    Dt = pd.read_parquet(p)
    Days = pd.read_csv(dp).drop_duplicates("date", keep="last") if os.path.exists(dp) else None
    return Dt, Days


def day_bins(sub, Days):
    L = pd.read_csv(f"{OUT}/{sub}/day_labels.csv")
    W = pd.read_csv(f"{OUT}/{sub}/ets_windows.csv")
    S = Days[Days.scanned.astype(bool)].copy()
    tmid = np.array([pd.Timestamp(d).timestamp() + 43200 for d in S.date])
    lab = AD.label_times(tmid, S.date.astype(str).to_numpy(), L, W)
    S["bin"] = tbin(lab["ets_phase"], lab["ets_rel_days"])
    S["in_ets"] = lab["in_ets_window"]
    return S


def section(sub):
    md = [f"## {sub}", ""]
    Wp = f"{OUT}/{sub}/ets_windows.csv"
    if os.path.exists(Wp):
        W = pd.read_csv(Wp)
        L = pd.read_csv(f"{OUT}/{sub}/day_labels.csv")
        md += [f"Day labels ({AD.LIFETIME[sub][0]}..{AD.LIFETIME[sub][1]}): "
               + ", ".join(f"{k} {v}" for k, v in L.label_source.value_counts().items())
               + f"; tremor days {int((L.tremor_day == 1).sum())}.  ETS windows ({len(W)}):", "",
               fmt(W), ""]
    summ, timing, curve, monthly = [], [], [], {}
    for kind in KINDS:
        Dt, Days = load(sub, kind)
        if Dt is None:
            continue
        nfam = Dt.family.nunique()
        known = Dt.in_ets_window.notna()
        out = (Dt.in_ets_window == 0)
        nscan = int(Days.scanned.astype(bool).sum()) if Days is not None else np.nan
        Sb = day_bins(sub, Days) if Days is not None else None
        d_in = int((Sb.in_ets == 1).sum()) if Sb is not None else np.nan
        d_out = int((Sb.in_ets == 0).sum()) if Sb is not None else np.nan
        summ.append(dict(kind=kind, families=nfam, detections=len(Dt), scanned_days=nscan,
                         det_known_label=int(known.sum()),
                         frac_outside_ets=float(out.sum() / max(known.sum(), 1)),
                         frac_unlabelled=float((~known).mean()) if len(Dt) else np.nan,
                         ets_days=d_in, non_ets_days=d_out,
                         rate_in_ets=float((Dt.in_ets_window == 1).sum() / max(d_in, 1)),
                         rate_outside_ets=float(out.sum() / max(d_out, 1)),
                         fam_with_outside_det=int(Dt[out].family.nunique())))
        for thr in NMAD_CURVE:
            m = Dt.nmad >= thr
            k = m & known
            curve.append(dict(kind=kind, nmad_ge=thr, detections=int(m.sum()),
                              frac_outside_ets=float((m & out).sum() / max(k.sum(), 1)),
                              rate_outside_per_day=float((m & out).sum() / max(d_out, 1)),
                              rate_in_ets_per_day=float((m & (Dt.in_ets_window == 1)).sum()
                                                        / max(d_in, 1))))
        b = tbin(Dt.ets_phase.to_numpy(object), Dt.ets_rel_days.to_numpy(float))
        cnt = pd.Series(b).value_counts()
        dcnt = Sb.bin.value_counts() if Sb is not None else pd.Series(dtype=int)
        for name in [x[2] for x in BINS] + ["unknown"]:
            timing.append(dict(kind=kind, bin=name, detections=int(cnt.get(name, 0)),
                               scanned_days=int(dcnt.get(name, 0)),
                               rate_per_day=float(cnt.get(name, 0) / dcnt.get(name, 0))
                               if dcnt.get(name, 0) else np.nan))
        mon = Dt.time.dt.strftime("%Y-%m")
        FM = pd.crosstab(Dt.family, mon)
        FM.to_csv(f"{OUT}/{sub}/family_month_{kind}.csv")
        monthly[kind] = mon.value_counts().sort_index()
        pd.DataFrame([t for t in timing if t["kind"] == kind]).to_csv(
            f"{OUT}/{sub}/timing_{kind}.csv", index=False)
    if not summ:
        return md + ["(no all-days detections yet)", ""]
    S = pd.DataFrame(summ)
    md += ["### Summary (every family, every scanned day; labels are columns, nothing filtered)", "",
           fmt(S), "",
           "`frac_outside_ets` = detections on LABELLED days outside every ETS window / detections "
           "with a known label.  `rate_*` = detections per scanned day of that kind.", ""]
    T = pd.DataFrame(timing)
    piv = T.pivot(index="bin", columns="kind", values="rate_per_day") \
        .reindex([x[2] for x in BINS] + ["unknown"])
    pc = T.pivot(index="bin", columns="kind", values="detections") \
        .reindex([x[2] for x in BINS] + ["unknown"])
    pdays = T[T.kind == "real"].set_index("bin").scanned_days if (T.kind == "real").any() else None
    tab = pd.concat({"det": pc, "det/day": piv}, axis=1)
    tab.columns = [f"{a}:{b}" for a, b in tab.columns]
    if pdays is not None:
        tab.insert(0, "scanned_days", pdays.reindex(tab.index))
    md += ["### Timing relative to ETS windows (pre / during / post), real vs null", "",
           fmt(tab.reset_index()), ""]
    C = pd.DataFrame(curve)
    cp = C.pivot(index="nmad_ge", columns="kind", values="frac_outside_ets")
    cn = C.pivot(index="nmad_ge", columns="kind", values="detections")
    cr = C.pivot(index="nmad_ge", columns="kind", values="rate_outside_per_day")
    ct = pd.concat({"n": cn, "frac_out": cp, "out/day": cr}, axis=1)
    ct.columns = [f"{a}:{b}" for a, b in ct.columns]
    md += ["### Separation curve over nmad (a curve, not a cut)", "", fmt(ct.reset_index()), ""]
    M = pd.DataFrame(monthly).fillna(0).astype(int)
    md += ["### Detections per month (all families; per-family counts in "
           f"`{sub}/family_month_<kind>.csv`)", "", fmt(M.reset_index().rename(
               columns={"index": "month", "time": "month"})), ""]
    Dt, _ = load(sub, "real")
    if Dt is not None and len(Dt):
        g = Dt.groupby("family").agg(n=("cc", "size"), n_out=("in_ets_window", lambda v: int((v == 0).sum())),
                                     n_during=("in_ets_window", lambda v: int((v == 1).sum())),
                                     nmad_med=("nmad", "median"), n_det_growth=(
                                         "n_det_growth", "first") if "n_det_growth" in Dt else ("cc", "size"))
        g["frac_out"] = g.n_out / (g.n_out + g.n_during).clip(lower=1)
        md += ["### Real families with the most detections outside ETS windows (top 15)", "",
               fmt(g.sort_values("n_out", ascending=False).head(15).reset_index()), ""]
        qs = []
        for kind in KINDS:
            Dk, _ = load(sub, kind)
            if Dk is None or not len(Dk):
                continue
            gk = Dk.groupby("family").in_ets_window.agg(
                lambda v: (v == 0).sum() / max(v.notna().sum(), 1))
            qs.append(dict(kind=kind, families=len(gk), **{f"q{int(q*100)}": float(gk.quantile(q))
                                                           for q in (0.1, 0.25, 0.5, 0.75, 0.9)}))
        md += ["### Per-family fraction of detections outside ETS windows (quantiles)", "",
               fmt(pd.DataFrame(qs)), ""]
    return md


def main():
    md = ["# All-days detection pass -- Royer & Bostock (2014) reproduction", "",
          f"Written {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())} by "
          "`scripts/ab_alldays_report.py`.", "",
          "Final growth templates (literal chain, built on tremor days) scanned over EVERY day with "
          ">= 4 template stations in the array lifetime; same 8 x MAD / 24 h network matched "
          "filter and +-15 s dedupe as growth.  Null = the chain's site-shift null final templates, "
          "same pass over the same real days.  Tremor-day and ETS labels are columns "
          "(rules in `scripts/ab_alldays.py` docstring).  Log: `ALLDAYS_LOG.md`.", ""]
    for sub in ("NW", "SW", "NVI"):
        md += section(sub)
    open(f"{OUT}/REPORT.md", "w").write("\n".join(md) + "\n")
    print(f"wrote {OUT}/REPORT.md")


if __name__ == "__main__":
    main()
