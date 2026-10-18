#!/usr/bin/env python
"""DUC_DIAG (2026-09-30) -- helper for the re-run after the cache GPU-OOM diagnosis
(results/duc_fix/DIAGNOSIS.md).  Nothing is deleted: superseded files are MOVED to *_pre_diag/.

  python scripts/duc_diag.py affected            print months (2008-01..2012-07) with any cache error
  python scripts/duc_diag.py setaside            move prep/<day>.npz of affected days, mf1/<m> of affected
                                                 months, retemplate outputs -> *_pre_diag/ (idempotent)
  python scripts/duc_diag.py cachecheck M        count residual cache errors of month M (plan csv)
  python scripts/duc_diag.py retpl_forced        re-templating forced at 2012-07 (diagnostic, scratch outputs
                                                 in results/duc_fix/diag_retpl/; production state untouched)
  python scripts/duc_diag.py summary             before/after table -> results/duc_fix/diag_after.csv and the
                                                 'Re-run results' section of DIAGNOSIS.md
"""
import os, sys, glob, shutil, json
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from duc_fix_common import *

M0, M1 = "2008-01", os.environ.get("DUC_DIAG_END", "2012-07")
PRE = f"{FIX}/pre_diag"


def plans():
    out = []
    for m in pd.period_range(M0, M1, freq="M"):
        for p in [f"{PRE}/plan_{m}.csv", f"{FIX}/data/plan_{m}.csv"]:   # the ORIGINAL plan if set aside
            if os.path.exists(p):
                P = pd.read_csv(p, dtype={"loc": str}).fillna({"loc": ""}); P["m"] = str(m); out.append(P); break
    return pd.concat(out, ignore_index=True)


def affected():
    P = plans(); P["err"] = P.cache.astype(str).str.startswith("error")
    days = sorted(P[P.err].day.unique()); months = sorted(P[P.err].m.unique())
    return days, months


def setaside():
    os.makedirs(PRE, exist_ok=True); os.makedirs(f"{FIX}/prep_pre_diag", exist_ok=True)
    os.makedirs(f"{FIX}/mf1_pre_diag", exist_ok=True)
    for m in pd.period_range(M0, M1, freq="M"):          # keep the original plans (they carry the errors)
        s, d = f"{FIX}/data/plan_{m}.csv", f"{PRE}/plan_{m}.csv"
        if os.path.exists(s) and not os.path.exists(d): shutil.copy2(s, d)
    days, months = affected(); nd = nm = 0
    for d in days:
        s = f"{FIX}/prep/{d}.npz"
        if os.path.exists(s): shutil.move(s, f"{FIX}/prep_pre_diag/{d}.npz"); nd += 1
    for m in months:
        for ext in [".parquet", "_diag.csv"]:
            s = f"{FIX}/mf1/{m}{ext}"
            if os.path.exists(s): shutil.move(s, f"{FIX}/mf1_pre_diag/{m}{ext}"); nm += 1
    for f in ["templates_new.csv", "templates_new.npz", "retemplate_state.csv", "rerun_months.txt"]:
        s = f"{FIX}/{f}"
        if os.path.exists(s): shutil.move(s, f"{PRE}/{f}")
    for s in glob.glob(f"{FIX}/data/retpl_*"):
        shutil.move(s, f"{PRE}/{os.path.basename(s)}")
    print(f"set aside: {nd} prep days, {nm} mf1 files; affected months {len(months)}")


def cachecheck(m):
    P = pd.read_csv(f"{FIX}/data/plan_{m}.csv", dtype={"loc": str})
    c = P.cache.astype(str).str.split(":").str[0].value_counts().to_dict()
    print(json.dumps(c))
    return c


def retpl_forced():
    """run duc_fix_retemplate --final on a scratch FIX (mf1 symlinked) -- what re-templating WOULD give if
    every new epoch were finalised with the data scanned so far (<= 2012-07, i.e. <= 9 months of window)."""
    import duc_fix_retemplate as R
    import duc_fix_common as C
    D = f"{FIX}/diag_retpl"; os.makedirs(D, exist_ok=True)
    if not os.path.exists(f"{D}/mf1"): os.symlink(f"{FIX}/mf1", f"{D}/mf1")
    if not os.path.exists(f"{D}/station_epochs.csv"): os.symlink(f"{FIX}/station_epochs.csv", f"{D}/station_epochs.csv")
    for f in ["templates_new.csv", "templates_new.npz", "retemplate_state.csv", "rerun_months.txt"]:
        if os.path.exists(f"{D}/{f}"): shutil.move(f"{D}/{f}", f"{D}/{f}.prev")
    R.FIX = D; R.STATE = f"{D}/retemplate_state.csv"; C.FIX = D
    sys.argv = ["x", "--upto", M1, "--final", "--workers", "4"]
    R.main()


def summary():
    rows = []
    for m in pd.period_range(M0, M1, freq="M"):
        m = str(m)
        a = f"{FIX}/mf1_pre_diag/{m}_diag.csv"; b = f"{FIX}/mf1/{m}_diag.csv"
        pa = f"{FIX}/mf1_pre_diag/{m}.parquet"; pb = f"{FIX}/mf1/{m}.parquet"
        if not os.path.exists(a): a, pa = b, pb          # unaffected month: before == after
        r = dict(month=m)
        for tag, dg, pq in [("before", a, pa), ("after", b, pb)]:
            if os.path.exists(dg) and os.path.exists(pq):
                g = pd.read_csv(dg); Dq = pd.read_parquet(pq, columns=["nchannel", "pass_his", "thr_perm"])
                r[f"ntpl_{tag}"] = round(float(g.n_tpl.mean()), 2)
                r[f"raw_{tag}"] = int(len(Dq)); r[f"rule_{tag}"] = int(Dq.pass_his.sum())
                r[f"nch_det_{tag}"] = round(float(Dq.nchannel.median()), 1) if len(Dq) else np.nan
        rows.append(r)
    S = pd.DataFrame(rows); S.to_csv(f"{FIX}/diag_after.csv", index=False)
    J = pd.read_csv(f"{FIX}/joins_2008_2011.csv") if os.path.exists(f"{FIX}/joins_2008_2011.csv") else None
    txt = ["", "## Re-run results (auto-written by `scripts/duc_diag.py summary`)", "",
           f"_{pd.Timestamp.utcnow():%Y-%m-%d %H:%M} UTC._ Per month: mean templates (channels) available per family-day, "
           "raw 8 x MAD detections, detections passing his cc x nchannel rule, median nchannel of the detections; "
           "before = the OOM-damaged scan (moved to `mf1_pre_diag/`), after = the re-run.", ""]
    try: txt.append(S.to_markdown(index=False))
    except Exception: txt.append("```\n" + S.to_string(index=False) + "\n```")
    Y = S.assign(y=S.month.str[:4]).groupby("y")[[c for c in S.columns if c.startswith(("raw_", "rule_"))]].sum()
    txt += ["", "Per year:", "", Y.to_markdown() if hasattr(Y, "to_markdown") else Y.to_string()]
    jy = f"{FIX}/joins_by_year_diag.csv"
    if os.path.exists(f"{FIX}/catalogue_all.parquet"):
        thr = thresholds()
        D = pd.read_parquet(f"{FIX}/catalogue_all.parquet", columns=["fam", "t", "pass_his"]); D["fam"] = D.fam.astype(str)
        H = his_catalogue(); H["thr"] = H.fam.map(thr).fillna(0.0); D["thr"] = D.fam.map(thr).fillna(0.0)
        H["hp"] = (H.thr > 0) & (H.cc * H.nch >= H.thr)
        jr = []
        for y in range(2008, 2012):
            d = D[pd.to_datetime(D.t).dt.year == y]; h = H[H.t.dt.year == y]
            jr.append(dict(year=y, ours_rule=int(d.pass_his.sum()), his_rule=int(h.hp.sum()),
                           ratio_rule=d.pass_his.sum() / max(h.hp.sum(), 1),
                           ours_raw=len(d), his_raw=len(h), ratio_raw=len(d) / max(len(h), 1),
                           ours_raw_60fam=int((d.thr > 0).sum()), his_raw_60fam=int((h.thr > 0).sum()),
                           ratio_raw_60fam=(d.thr > 0).sum() / max((h.thr > 0).sum(), 1)))
        pd.DataFrame(jr).round(3).to_csv(jy, index=False)
    if os.path.exists(jy):
        txt += ["", "**New joins with his catalogue (his rule both sides; from `duc_fix_report.py`):**", "",
                pd.read_csv(jy).to_markdown(index=False)]
    rs = f"{FIX}/diag_retpl/templates_new.csv"
    if os.path.exists(rs):
        T = pd.read_csv(rs, dtype={"fam": str})
        txt += ["", f"Forced re-templating diagnostic (all new epochs finalised at {M1}, scratch `diag_retpl/`): "
                f"{int(T.keep.sum())}/{len(T)} family-channels kept, {T[T.keep == 1].epoch.nunique()}/{T.epoch.nunique()} epochs; "
                f"median N windows {T.N.median():.0f}, median candidates {T.n_cand.median():.0f}."]
    st = f"{FIX}/retemplate_state.csv"
    if os.path.exists(st):
        s = pd.read_csv(st); txt += ["", f"Production re-templating state (no --final): {len(s)} epochs finalised."]
    else:
        txt += ["", "Production re-templating state (no --final): no epoch finalised yet (all windows still open at "
                f"{M1}; they finalise when the scan passes E.start + 24 months or 150 candidates per family)."]
    p = f"{FIX}/DIAGNOSIS.md"
    body = open(p).read() if os.path.exists(p) else ""
    k = body.find("\n## Re-run results")
    if k >= 0: body = body[:k]
    open(p, "w").write(body.rstrip("\n") + "\n" + "\n".join(txt) + "\n")
    print(S.to_string(index=False))


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "affected": print(" ".join(affected()[1]))
    elif c == "setaside": setaside()
    elif c == "cachecheck": cachecheck(sys.argv[2])
    elif c == "retpl_forced": retpl_forced()
    elif c == "summary": summary()
