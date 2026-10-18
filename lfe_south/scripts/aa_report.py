#!/usr/bin/env python
"""Agent AA report: what the Shelly matched filter does before and after the two fixes.

Reads whatever rows of the run exist in results/aa_shelly/ and writes AA_SHELLY_REPORT.md.
RANK-DON'T-DISCARD: every family keeps its row at every threshold; the threshold is a column.
"""
import os,sys,glob,json
import numpy as np, pandas as pd
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from t_common import *
AA=f"{RES}/aa_shelly"
ROWS=[("A","mfA","12 s templates, cache band 1.5-9 Hz, moveout-aligned"),
      ("B","mfB","12 s templates, 2-9 Hz (Shelly's band), moveout-aligned"),
      ("C","mfC","6 s templates (Shelly's length), 2-9 Hz, moveout-aligned"),
      ("N","mfN","ABLATION: normalisation fixed but NOT moveout-aligned (t11_mf.py stacking)"),
      ("D","mfD","row A templates re-stacked on our own build-day detections (Shelly's iteration)")]
def qs(x,q):
    x=np.asarray(x,float); x=x[np.isfinite(x)]
    return float(np.percentile(x,q)) if x.size else np.nan
out=[]
out.append("# Shelly (2025) reproduction — the 12.8 M detections, diagnosed and fixed\n")
out.append("Agent AA, 2026-09-25.  Companion: `notes/PAPER_METHODS_TABLE.md` §S25 rows (his actual "
           "parameters, read from the paper) and `notes/CATALOGUE_SUMMARY.md` §10.\n")
out.append("## 0. The fault, measured\n")
out.append("`t11_mf.py` normalises the cc by `sd = sqrt(max(M*S2 - S1*S1, 1e-30))` from cumulative "
           "sums.  The 20 Hz cache zero-fills gaps (`t3_cache.py:31-32`, obspy `merge(fill_value=0)`"
           " + `trim(pad=True, fill_value=0)`) and stores no validity mask, so a gap is an exactly "
           "constant run: the true window variance is 0, the cumsum difference goes negative from "
           "round-off, the clamp gives sd = 1e-15, and the 'normalised' cc reaches 1e15.\n")
out.append("| measurement | value |\n|---|---|")
out.append("| BK.RVIT HHZ, family 23, 2022-05-11: constant-sample fraction | 0.333 |")
out.append("| same channel, max abs 'normalised' cc | 2.4e15 |")
out.append("| NC.KIP EHZ same day: constant fraction / max cc | 0.501 / 2.8e10 |")
out.append("| channel-summed cc series that day: MAD (cannot exceed nchan=20) | 25.7 |")
out.append("| fraction of samples above 8xMAD | 0.166 |")
out.append("| detections that ONE day, 27 families, t11_mf.py | 128,922 |")
out.append("| detections that ONE day, 27 families, fixed | 21 |")
out.append("| Shelly's own detections that day, 27 families | 71 |\n")
out.append("Second fault: `t11_mf.py` sums each station's cc at a common lag, but every station's "
           "template is cut at its own offset tt_s+tau_s (`t9_templates.py:93`), which for family 1 "
           "spans 1.49-5.85 s.  The stack was therefore incoherent and its time axis was not "
           "Shelly's.  After shifting each channel by (tt_s+tau_s-pre) the peak at his published "
           "times lands at dt median 0.00 s.\n")
out.append("Third fault, in the SCORING: `t8_score.py` used all 61,441 published events as the "
           "denominator while the filter ran only on the held-out days, where only 5,021 of his "
           "events live — the recall cap was 0.082, so the published '0.026' was 0.32 of the "
           "achievable maximum.  `t8_score_aa.py --days-mode source|fam` fixes the denominator.\n")
for tag,stem,desc in ROWS:
    dg=f"{AA}/{stem}_diag.csv"; sc=f"{AA}/{stem}_score.csv"; ts=f"{AA}/{stem}_truthstat.csv"
    if not os.path.exists(dg): continue
    D=pd.read_csv(dg)
    out.append(f"\n## Row {tag} — {desc}\n")
    out.append(f"* family-days processed **{len(D)}**, channels/family-day median "
               f"**{D.nchan.median():.0f}**, flat-masked fraction median **{D.flat_frac.median():.3f}**")
    out.append(f"* nmad of the daily statistic: p50 **{qs(D.nm_p50,50):.2f}**, p99 "
               f"**{qs(D.nm_p99,50):.2f}**, p99.9 **{qs(D.nm_p999,50):.2f}**, day max median "
               f"**{qs(D.nm_max,50):.1f}** (medians over family-days)")
    sweep=[c for c in D.columns if c.startswith("n") and c[1:].isdigit()]
    if sweep:
        out.append("* detections per family-day at each threshold (mean / median):")
        out.append("\n| threshold xMAD | mean det/family-day | median | max |\n|---|---|---|---|")
        for c in sorted(sweep,key=lambda s:int(s[1:])):
            out.append(f"| {c[1:]} | {D[c].mean():.1f} | {D[c].median():.0f} | {D[c].max():.0f} |")
    if "n_theirs" in D.columns:
        q=D[D.n_theirs==0]; nz=D[D.n_theirs>0]
        out.append(f"\n* **quiet-day false-alarm rate** (family-days where Shelly has NO detection, "
                   f"n={len(q)}): mean **{q.ndet.mean():.2f}** det/family-day at 9xMAD, "
                   f"vs **{nz.ndet.mean():.2f}** on his active family-days (n={len(nz)}); "
                   f"his own rate on active family-days is **{nz.n_theirs.mean():.2f}**")
        for c in sorted(sweep,key=lambda s:int(s[1:])):
            out.append(f"  * at {c[1:]}xMAD: quiet **{q[c].mean():.2f}** vs active **{nz[c].mean():.2f}** det/family-day")
    if os.path.exists(ts):
        T=pd.read_csv(ts)
        out.append(f"\n* **separation curve** — our statistic at his {len(T)} published times on "
                   f"these days: nmad p10 **{qs(T.nm_pm2s,10):.2f}**, p50 **{qs(T.nm_pm2s,50):.2f}**, "
                   f"p90 **{qs(T.nm_pm2s,90):.2f}** (best within +/-2 s), against a noise p99.9 of "
                   f"**{qs(D.nm_p999,50):.2f}**; timing offset dt median "
                   f"**{qs(T.dt_best,50):.2f}** s")
        for thr in [6,7,8,9,10,12]:
            out.append(f"  * fraction of his events above {thr}xMAD: **{float((T.nm_pm2s>=thr).mean()):.3f}**")
    if os.path.exists(sc):
        S=pd.read_csv(sc)
        out.append("\n| tmin | his events | ours | matched | pooled recall | precision | fam with >=1 match | fam >=30 % precision |\n|---|---|---|---|---|---|---|---|")
        for tm,g in S.groupby("tmin",dropna=False):
            out.append(f"| {tm} | {int(g.n_theirs.sum())} | {int(g.n_ours.sum())} | {int(g.n_match.sum())} | "
                       f"{g.n_match.sum()/max(g.n_theirs.sum(),1):.3f} | "
                       f"{g.n_match.sum()/max(g.n_ours.sum(),1):.4f} | {int((g.n_match>0).sum())}/{len(g)} | "
                       f"{int(g.bostock30.sum())} |")
        best=S[S.tmin==S.tmin.max()] if S.tmin.notna().any() else S
        out.append("\nPer-family rows at every threshold are in `"+os.path.basename(sc)+"` — no family is dropped.")
open(f"{AA}/AA_SHELLY_REPORT.md","w").write("\n".join(out)+"\n")
log(f"[report] {AA}/AA_SHELLY_REPORT.md ({len(out)} lines)")
