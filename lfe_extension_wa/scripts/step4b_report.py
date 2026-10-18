#!/usr/bin/env python
"""Step 4b report: per-family before/after table + results/STEP4B_RESULTS.md.

Runs (mf tag -> window, pass): step4b_v2010b / step4b_v2010 (V2010 before / after, validation),
step4b_<e>b / step4b_<e>p1 / step4b_<e>p2 for e in t2012, t2013, t2014 (before = the 42-station templates
only; p1 = + templates from array-era detections and earlier episodes; p2 = + templates from this
episode's pass-1 detections of the separating families).

Per family and run, over the family's OBSERVABLE (>= 2 sites) days of the window that are tremor days of
its region (PNSN tremor-day table, a label only) -- and separately its quiet (non-episode) days:
  sites: median / max observable sites; clear-template sites
  counts at each tier (frozen T_core 8.75, T_ext 8.25; per-epoch core_ep): data (tier), data catalogue
  (cat: HF pass + cross-family kept), scrambled control, second scrambled control (own shifts)
  separation = (data + 1) / (control + 1); the SAME statistic on the control: (control2 + 1) / (control + 1)
  validation (V2010): recall / precision against the test-half array detections (wa_score columns)
No row is dropped; OLY:938 is flagged.
  -> results/step4b/step4b_by_family.csv (long), results/step4b/step4b_summary.csv, results/STEP4B_RESULTS.md
"""
from __future__ import annotations
import os, sys, glob, json
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

PFX = os.environ.get("STEP4B_PREFIX", "step4b_")   # smoke: step4b_smoke_
OUT = f"{W.RES}/step4b" + ("" if PFX == "step4b_" else "/smoke")
os.makedirs(OUT, exist_ok=True)
TREM = {"Royer_NW_WA": "NW_WA", "Royer_SW_WA": "SW_WA", "Olympic_2010": "OLY2010", "Olympic_2010_armB": "OLY2010"}
RUNS = [("V2010", "before", f"{PFX}v2010b"), ("V2010", "after", f"{PFX}v2010")]
for e in ["T2012", "T2013", "T2014"]:
    RUNS += [(e, "before", f"{PFX}{e.lower()}b"), (e, "after_p1", f"{PFX}{e.lower()}p1"),
             (e, "after", f"{PFX}{e.lower()}p2")]
TPL_BEFORE = ["build", "t2012new"]
TPL_AFTER = TPL_BEFORE + [f"{PFX}a", f"{PFX}t2012", f"{PFX}t2013", f"{PFX}t2014"]
SEPS = [1.5, 2.0, 3.0, 5.0]


def f2(x, n=2):
    return "–" if x is None or (isinstance(x, (float, np.floating)) and not np.isfinite(x)) else f"{x:.{n}f}"


def tremor_days():
    out = {}
    for rg, tr in TREM.items():
        t = pd.read_csv(f"{W.EXT}/{tr}_tremor_days.csv")
        out[rg] = set(pd.to_datetime(t[t.columns[0]]).dt.strftime("%Y-%m-%d")[t[t.columns[1]] > 0])
    return out


def clear_sites(tags):
    Q = []
    for t in tags:
        p = f"{W.RES}/templates/tpl_qc_{t}.csv"
        if os.path.exists(p):
            Q.append(pd.read_csv(p))
    if not Q:
        return pd.Series(dtype=int), pd.DataFrame()
    Q = pd.concat(Q)
    Q = Q[Q.clear.fillna(False).astype(bool)].drop_duplicates(["fid", "netsta", "node"])
    S = W.load_stations()
    Q["site"] = Q.netsta.map(dict(zip(S.netsta, S.site)))
    return Q.groupby("fid").site.nunique(), Q


def run_table(win, pas, mf, F, td):
    od = f"{W.RES}/mf_{mf}"
    if not glob.glob(f"{od}/famday_*.parquet"):
        return None
    FD = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{od}/famday_*.parquet"))], ignore_index=True)
    DET = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{od}/det_*.parquet"))], ignore_index=True)
    sd = f"{W.RES}/score_{mf}"
    TH = pd.read_csv(f"{sd}/tiers_by_epoch.csv") if os.path.exists(f"{sd}/tiers_by_epoch.csv") else pd.DataFrame()
    tce = float(TH.T_core_ep.iloc[0]) if len(TH) else W.T_CORE
    tee = float(TH.T_ext_ep.iloc[0]) if len(TH) else W.T_EXT
    DF = pd.read_parquet(f"{sd}/detections_flagged.parquet", columns=["fid", "t", "date", "cat_core", "cat_ext", "cat_core_ep"]) \
        if os.path.exists(f"{sd}/detections_flagged.parquet") else DET.iloc[:0].assign(cat_core=False, cat_ext=False, cat_core_ep=False)
    thr = {"core": W.T_CORE, "ext": W.T_EXT, "core_ep": tce, "ext_ep": tee}
    reg = dict(zip(F.fid, F.region))
    FD = FD[FD.fid.isin(reg)]
    FD["trem"] = [d in td[reg[f]] for f, d in zip(FD.fid, FD.date)]
    rows = []
    DG = dict(tuple(DET.groupby("fid"))); CG = dict(tuple(DF.groupby("fid"))); FG = dict(tuple(FD.groupby("fid")))
    for r in F.itertuples():
        fd = FG.get(r.fid, FD.iloc[:0])
        fo = fd[fd.observable]
        tdays = set(fo[fo.trem].date); qdays = set(fo[fo.kind != "episode"].date)
        row = dict(fid=r.fid, family_uid=r.family_uid, region=r.region, flag=r.flag if isinstance(r.flag, str) else "",
                   window=win, pass_=pas, mf=mf, n_days=len(fd), n_obs_days=len(fo), n_obs_tremor_days=len(tdays),
                   n_obs_quiet_days=len(qdays), med_sites=float(fo.n_sites.median()) if len(fo) else 0.0,
                   max_sites=int(fd.n_sites.max()) if len(fd) else 0,
                   T_core_ep=tce, T_ext_ep=tee)
        dd = DG.get(r.fid, DET.iloc[:0]); cc = CG.get(r.fid, DF.iloc[:0])
        for tier, T in thr.items():
            for lab, ser in [("data", "data"), ("ctl", "control"), ("ctl2", "control2")]:
                x = dd[(dd.series == ser) & (dd.net_mad >= T)]
                row[f"n_{lab}_{tier}_trem"] = int(x.date.isin(tdays).sum())
                row[f"n_{lab}_{tier}_quiet"] = int(x.date.isin(qdays).sum())
            if f"cat_{tier}" in cc:
                x = cc[cc[f"cat_{tier}"]]
                row[f"n_cat_{tier}_trem"] = int(x.date.isin(tdays).sum())
            nt = max(len(tdays), 1)
            row[f"rate_data_{tier}"] = row[f"n_data_{tier}_trem"] / nt if tdays else np.nan
            row[f"rate_ctl_{tier}"] = row[f"n_ctl_{tier}_trem"] / nt if tdays else np.nan
            row[f"rate_ctl2_{tier}"] = row[f"n_ctl2_{tier}_trem"] / nt if tdays else np.nan
            row[f"sep_{tier}"] = (row[f"n_data_{tier}_trem"] + 1) / (row[f"n_ctl_{tier}_trem"] + 1) if tdays else np.nan
            row[f"sep_ctl_{tier}"] = (row[f"n_ctl2_{tier}_trem"] + 1) / (row[f"n_ctl_{tier}_trem"] + 1) if tdays else np.nan
        rows.append(row)
    R = pd.DataFrame(rows)
    vp = f"{sd}/validation_by_family.csv"
    if os.path.exists(vp):
        V = pd.read_csv(vp)
        V = V[V.window == win][["fid", "n_truth_test", "recall_cat_core", "prec_own_cat_core", "prec_own_ctl_core",
                                "recall_cat_core_ep", "prec_own_cat_core_ep", "recall_cat_ext", "prec_own_cat_ext"]]
        R = R.merge(V, on="fid", how="left")
    return R


def main():
    F = pd.read_csv(f"{W.SETUP}/families.csv")
    td = tremor_days()
    sb, _ = clear_sites(TPL_BEFORE)
    sa, QA = clear_sites(TPL_AFTER)
    F["tpl_sites_before"] = F.fid.map(sb).fillna(0).astype(int)
    F["tpl_sites_after"] = F.fid.map(sa).fillna(0).astype(int)
    tabs = [t for t in (run_table(w, p, m, F, td) for w, p, m in RUNS) if t is not None]
    if not tabs:
        print("no runs yet"); return
    L = pd.concat(tabs, ignore_index=True).merge(
        F[["fid", "tpl_sites_before", "tpl_sites_after", "sites_before", "p_ge3_after", "projected_ge3", "lat_use", "lon_use"]],
        on="fid", how="left")
    L.to_csv(f"{OUT}/step4b_by_family.csv", index=False)
    # ---------------- summary per window x pass
    S = []
    for (w, p), g in L.groupby(["window", "pass_"], sort=False):
        o = g[g.n_obs_tremor_days > 0]
        row = dict(window=w, pass_=p, fam_observable=len(o), fam_ge3_sites=int((o.med_sites >= 3).sum()))
        for tier in ["core", "core_ep"]:
            for s in SEPS:
                row[f"sep{s:g}_{tier}_data"] = int((o[f"sep_{tier}"] >= s).sum())
                row[f"sep{s:g}_{tier}_ctl"] = int((o[f"sep_ctl_{tier}"] >= s).sum())
                o3 = o[o.med_sites >= 3]
                row[f"sep{s:g}_{tier}_data_ge3"] = int((o3[f"sep_{tier}"] >= s).sum())
                row[f"sep{s:g}_{tier}_ctl_ge3"] = int((o3[f"sep_ctl_{tier}"] >= s).sum())
            row[f"data_rate_{tier}"] = o[f"n_data_{tier}_trem"].sum() / max(o.n_obs_tremor_days.sum(), 1)
            row[f"ctl_rate_{tier}"] = o[f"n_ctl_{tier}_trem"].sum() / max(o.n_obs_tremor_days.sum(), 1)
            row[f"ctl2_rate_{tier}"] = o[f"n_ctl2_{tier}_trem"].sum() / max(o.n_obs_tremor_days.sum(), 1)
        if "recall_cat_core" in g and g.n_truth_test.notna().any():
            v = g[g.n_truth_test.fillna(0) >= 20]
            row["val_fam_ge20truth"] = len(v)
            row["val_recall_core_med"] = float(v.recall_cat_core.median()) if len(v) else np.nan
            row["val_prec_core_med"] = float(v.prec_own_cat_core.median()) if len(v) else np.nan
        S.append(row)
    S = pd.DataFrame(S)
    S.to_csv(f"{OUT}/step4b_summary.csv", index=False)
    # ---------------- markdown
    C = pd.read_csv(f"{W.RES}/step4b/candidates.csv")
    CS = C[C.selected]
    fm = pd.read_csv(f"{W.RES}/fetch_manifest.csv")
    SD = pd.read_csv(f"{W.SETUP}/station_days.csv", keep_default_na=False)
    keys = set(zip(SD.netsta, SD.date))
    fm4 = fm[[k in keys for k in zip(fm.netsta, fm.date)]]
    o = ["# Step 4b — added stations for the NW and Olympic families", "",
         f"Written {pd.Timestamp.now('UTC'):%Y-%m-%d %H:%MZ} by `scripts/step4b_report.py`. Plan: `PLAN.md` §Step 4b; log `NOTES.md`.", "",
         "## Stations added", "",
         f"- {len(CS)} stations ({CS.groupby('net').size().to_dict()}); tier 1 = permanent / served later arrays "
         f"({int((CS.tier == 1).sum())}), tier 2 = TA / 7D taken only where they are the only option ({int((CS.tier == 2).sum())}).",
         f"- Station-days in the step-4b windows: {int(SD.meta_live.astype(str).isin(['True']).sum())} live; downloaded "
         f"(this step, incl. the 42 trial stations in 2013/2014/A2011): {int((fm4.status == 'data').sum())} station-days, "
         f"{fm4.bytes.sum() / 1e9:.1f} GB; nodata {int((fm4.status == 'nodata').sum())}.",
         "- Candidate table (every candidate, every column, exclusion reason): `results/step4b/candidates.csv`.", "",
         "| station | tier | km to nearest family | families ≤30 km | comps | live windows | clear templates (families) |",
         "|---|---|---|---|---|---|---|"]
    ncl = QA.groupby("netsta").fid.nunique() if len(QA) else pd.Series(dtype=int)
    for r in CS.sort_values("dmin_km").itertuples():
        lw = ",".join(w for w in ["V2010", "A2011", "T2012", "T2013", "T2014"] if getattr(r, f"ver_{w}"))
        o.append(f"| {r.netsta} | {r.tier} | {r.dmin_km:.1f} | {r.n_fam_30km} | {r.n_comp} | {lw} | {int(ncl.get(r.netsta, 0))} |")
    o += ["", "## Clear-template sites per family (before = the 42 trial stations; after = + added stations / epochs)", "",
          "| region | families | ≥3 sites before | ≥3 sites after | projected ≥3 (P ≥ 0.5) | ≥2 before | ≥2 after |",
          "|---|---|---|---|---|---|---|"]
    for rg, g in F.groupby("region"):
        o.append(f"| {rg} | {len(g)} | {int((g.tpl_sites_before >= 3).sum())} | {int((g.tpl_sites_after >= 3).sum())} | "
                 f"{int(g.projected_ge3.sum())} | {int((g.tpl_sites_before >= 2).sum())} | {int((g.tpl_sites_after >= 2).sum())} |")
    o.append(f"| ALL | {len(F)} | {int((F.tpl_sites_before >= 3).sum())} | {int((F.tpl_sites_after >= 3).sum())} | "
             f"{int(F.projected_ge3.sum())} | {int((F.tpl_sites_before >= 2).sum())} | {int((F.tpl_sites_after >= 2).sum())} |")
    o += ["", "## Before / after per episode", "",
          "Rates per family-day on the family's observable tremor days (label only). Separation = (data + 1) / "
          "(scrambled control + 1) per family; the same statistic on the control = (second scrambled control + 1) / "
          "(scrambled control + 1). `sep ≥ 2` counts families; `ctl` = how many the control produces by the same rule.", "",
          "| window | pass | fam observable | fam ≥3 sites | sep≥2 data / ctl (core) | sep≥2 data / ctl, ≥3 sites | sep≥3 data / ctl | data / ctl / ctl2 rate (core) | sep≥2 data / ctl (core_ep) | val recall / prec (med, fam ≥20 truth) |",
          "|---|---|---|---|---|---|---|---|---|---|"]
    for r in S.itertuples():
        v = f"{f2(getattr(r, 'val_recall_core_med', np.nan))} / {f2(getattr(r, 'val_prec_core_med', np.nan))}" \
            if hasattr(r, "val_recall_core_med") and np.isfinite(getattr(r, "val_recall_core_med", np.nan)) else "–"
        o.append(f"| {r.window} | {r.pass_} | {r.fam_observable} | {r.fam_ge3_sites} | {S.loc[r.Index, 'sep2_core_data']} / {S.loc[r.Index, 'sep2_core_ctl']} | "
                 f"{S.loc[r.Index, 'sep2_core_data_ge3']} / {S.loc[r.Index, 'sep2_core_ctl_ge3']} | {S.loc[r.Index, 'sep3_core_data']} / {S.loc[r.Index, 'sep3_core_ctl']} | "
                 f"{f2(r.data_rate_core)} / {f2(r.ctl_rate_core)} / {f2(r.ctl2_rate_core)} | {S.loc[r.Index, 'sep2_core_ep_data']} / {S.loc[r.Index, 'sep2_core_ep_ctl']} | {v} |")
    o += ["", "Separation curve (families with sep ≥ x, core tier; data vs control):", "",
          "| window | pass | " + " | ".join(f"≥{s:g}" for s in SEPS) + " |", "|---|---|" + "---|" * len(SEPS)]
    for r in S.itertuples():
        o.append(f"| {r.window} | {r.pass_} | " + " | ".join(
            f"{S.loc[r.Index, f'sep{x:g}_core_data']} / {S.loc[r.Index, f'sep{x:g}_core_ctl']}" for x in SEPS) + " |")
    o += ["", "## Families ranked by separation (after; core tier)", "",
          "Top 30 per window. Every family is in `results/step4b/step4b_by_family.csv` (all runs, all columns).", ""]
    for w in ["V2010", "T2012", "T2013", "T2014"]:
        a = L[(L.window == w) & (L.pass_ == "after") & (L.n_obs_tremor_days > 0)]
        if not len(a):
            continue
        b = L[(L.window == w) & (L.pass_ == "before")].set_index("fid")
        o += [f"### {w}", "", "| family | region | sites before → after (med obs) | tpl sites b/a | data / ctl rate | sep (ctl sep) | sep before | recall / prec | flag |",
              "|---|---|---|---|---|---|---|---|---|"]
        for r in a.sort_values("sep_core", ascending=False).head(30).itertuples():
            bb = b.loc[r.fid] if r.fid in b.index else None
            rp = f"{f2(getattr(r, 'recall_cat_core', np.nan))} / {f2(getattr(r, 'prec_own_cat_core', np.nan))}" if "recall_cat_core" in a else "–"
            o.append(f"| {r.family_uid} | {r.region} | {f2(bb.med_sites, 0) if bb is not None else '–'} → {f2(r.med_sites, 0)} | "
                     f"{r.tpl_sites_before} / {r.tpl_sites_after} | {f2(r.rate_data_core)} / {f2(r.rate_ctl_core)} | "
                     f"{f2(r.sep_core)} ({f2(r.sep_ctl_core)}) | {f2(bb.sep_core) if bb is not None else '–'} | {rp} | {r.flag} |")
        o.append("")
    o += ["## Flag", "", "- OLY:938: runaway (~200 detections/day in data AND in the scrambled control; its sites are the "
          "B013/B943 borehole pair + neighbour). Kept in every table, flagged in `flag`.", ""]
    open(f"{W.RES}/STEP4B_RESULTS.md" if PFX == "step4b_" else f"{OUT}/STEP4B_RESULTS_smoke.md", "w").write("\n".join(o) + "\n")
    print("\n".join(o[:80]))


if __name__ == "__main__":
    main()
