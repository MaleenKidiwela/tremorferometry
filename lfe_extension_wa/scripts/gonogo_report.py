#!/usr/bin/env python
"""GO/NO-GO report tables (2026-10-02): reads results/wa_gonogo/{truth_array_frac*.csv, score_*/} and writes
results/wa_gonogo/per_family.csv, comparison.csv and the numeric tables (tables.md) pasted into GONOGO_RESULTS.md."""
import os, sys
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
GO = f"{os.path.dirname(HERE)}/results/wa_gonogo"
RG = [("Royer_NW_WA", "NW"), ("Royer_SW_WA", "SW"), ("Olympic_2010", "Olympic Arm A"), ("Olympic_2010_armB", "Arm B")]
out = []


def famtab(run):
    B = pd.read_csv(f"{GO}/score_{run}/by_family.csv").copy()
    B["pass_both"] = B.sepboth_core >= 2
    B["pass_c1"] = B.sep_core >= 2
    B["ctl_pass"] = B.sep_ctl_core >= 2
    return B


def agg(B, sub):
    g = B[sub]
    g2 = g[g.n_days_obs > 0]
    no = g2.n_ours_cat_core.sum(); fd = g2.n_days_obs.sum()
    nr = g2.n_ours_raw_core.sum(); nc = g2.n_ours_ctl_core.sum()
    return dict(fam_obs=g2.fid.nunique(), famwin_label=int((g.n_label_days > 0).sum()),
                fw_both=int(g.pass_both.sum()), fam_both=int(g[g.pass_both].fid.nunique()),
                fw_c1=int(g.pass_c1.sum()), fw_ctl=int(g.ctl_pass.sum()),
                prec_own=g2.n_match_own_cat_core.sum() / no if no else np.nan,
                prec_nbr3=g2.n_prec_nbr3_cat_core.sum() / no if no else np.nan,
                prec_own_raw=g2.n_match_own_raw_core.sum() / nr if nr else np.nan,
                prec_own_ctl=g2.n_match_own_ctl_core.sum() / nc if nc else np.nan,
                prec_nbr3_ctl=g2.n_prec_nbr3_ctl_core.sum() / nc if nc else np.nan,
                recall_own=g2.n_match_own_cat_core.sum() / g2.n_truth_test.sum() if g2.n_truth_test.sum() else np.nan,
                data_fd=(g2.data_per_famday_core * g2.n_days_obs).sum() / fd if fd else np.nan,
                ctl_fd=(g2.ctl_per_famday_core * g2.n_days_obs).sum() / fd if fd else np.nan,
                n_ours=int(no), n_truth=int(g2.n_truth_test.sum()))


rows = []
for sup, runs in [("main", [("orig42", "42-station T1 (original)"), ("core21", "21 core (20 usable)")]),
                  ("V2008 supplement", [("orig42_v2008", "42-station T1 (original)"), ("core21_v2008", "21 core (20 usable)")])]:
    for run, lab in runs:
        B = famtab(run)
        for rg, rn in RG:
            for w in (["V2010", "A2011", "both"] if sup == "main" else ["V2008"]):
                sub = (B.region == rg) & ((B.window == w) if w != "both" else B.window.isin(["V2010", "A2011"]))
                if not sub.any():
                    rows.append(dict(test=sup, run=run, network=lab, region=rn, window=w, note="no family-window"))
                    continue
                rows.append(dict(test=sup, run=run, network=lab, region=rn, window=w, **agg(B, sub)))
C = pd.DataFrame(rows)
C.to_csv(f"{GO}/comparison.csv", index=False)

# per family (main test)
c = famtab("core21"); o = famtab("orig42")
keep = ["fid", "family_uid", "region", "n_tpl30_core", "n_sites30_core", "n_clear_sites_core", "tpl30_core_stations"]
P = c.drop_duplicates("fid")[keep].copy()
for run, B in [("core21", c), ("orig42", o)]:
    for w in ["V2010", "A2011"]:
        b = B[B.window == w].set_index("fid")
        P[f"sepboth_{w}_{run}"] = P.fid.map(b.sepboth_core)
        P[f"lab_data_{w}_{run}"] = P.fid.map(b.lab_data_core)
        P[f"lab_ctl_{w}_{run}"] = P.fid.map(b.lab_ctl_core)
        P[f"lab_ctl2_{w}_{run}"] = P.fid.map(b.lab_ctl2_core)
    g = B.groupby("fid")
    P[f"pass_any_{run}"] = P.fid.map(g.pass_both.any())
    P[f"n_pass_windows_{run}"] = P.fid.map(g.pass_both.sum())
    P[f"prec_own_{run}"] = P.fid.map(g.n_match_own_cat_core.sum() / g.n_ours_cat_core.sum().replace(0, np.nan))
    P[f"n_ours_{run}"] = P.fid.map(g.n_ours_cat_core.sum())
F = pd.read_csv(f"{GO}/setup/families.csv")
P = P.merge(F[["fid", "located", "lat_use", "lon_use"]], on="fid", how="left")
P.to_csv(f"{GO}/per_family.csv", index=False)

# ---- markdown tables
L = []
def f3(x):
    return "–" if pd.isna(x) else f"{x:.3f}"
L.append("| test | region | window | network | obs fam | fam-windows ≥ 2× both controls (distinct families) | control-vs-control passes | own prec | ≤3-km prec | own prec, scrambled control | data / control per fam-day | own recall |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
for r in C.itertuples():
    if isinstance(getattr(r, "note", np.nan), str):
        L.append(f"| {r.test} | {r.region} | {r.window} | {r.network} | – | no family-window (no truth in this window) | | | | | | |")
        continue
    L.append(f"| {r.test} | {r.region} | {r.window} | {r.network} | {int(r.fam_obs)} | {int(r.fw_both)} ({int(r.fam_both)}) | {int(r.fw_ctl)} | "
             f"{f3(r.prec_own)} | {f3(r.prec_nbr3)} | {f3(r.prec_own_ctl)} | {r.data_fd:.2f} / {r.ctl_fd:.2f} | {f3(r.recall_own)} |")
L.append("")
L.append("Per region, template coverage on the 21 core (all 329 families, build-half templates `results/templates/tpl_qc_build.csv`; sites merged at 3 km):")
L.append("")
L.append("| region | families | any ≥ 30-event template on core | ≥ 2 core sites with a ≥ 30-event template | ≥ 1 clear core site | ≥ 2 clear core sites | pass (≥ 2× both controls, V2010 or A2011) on core | of these with ≥ 2 clear sites | pass on 42-station T1 | pass on both |")
L.append("|---|---|---|---|---|---|---|---|---|---|")
q = pd.read_csv(f"{os.path.dirname(GO)}/templates/tpl_qc_build.csv"); q["clear"] = q.clear.fillna(False).astype(bool)
S = pd.read_csv(f"{GO}/setup/stations.csv"); cs = dict(zip(S.netsta, S.site))
for rg, rn in RG:
    fids = F[F.region == rg].fid
    q30 = q[(q.n_events >= 30) & q.offset_s.notna() & q.w.notna() & q.netsta.isin(cs) & q.fid.isin(fids)]
    qc = q[q.clear & q.netsta.isin(cs) & q.fid.isin(fids)]
    s30 = q30.groupby("fid").netsta.apply(lambda v: len({cs[x] for x in v}))
    scl = qc.groupby("fid").netsta.apply(lambda v: len({cs[x] for x in v}))
    pr = P[P.region == rg]
    pc = pr.pass_any_core21.fillna(False).astype(bool); po = pr.pass_any_orig42.fillna(False).astype(bool)
    L.append(f"| {rn} | {len(fids)} | {s30.size} | {(s30 >= 2).sum()} | {scl.size} | {(scl >= 2).sum()} | {int(pc.sum())} | "
             f"{int((pc & (pr.n_clear_sites_core >= 2)).sum())} | {int(po.sum())} | {int((pc & po).sum())} |")
L.append("")
L.append("NW own precision by the family-month's truth_array_frac (Step A; V2010 -> 2010-08, A2011 -> 2011-08):")
L.append("")
L.append("| window | truth_array_frac | family-windows | own prec 42-station | own prec core | pass 42-station | pass core |")
L.append("|---|---|---|---|---|---|---|")
A = pd.read_csv(f"{GO}/truth_array_frac.csv")
mm = {"V2010": "2010-08", "A2011": "2011-08"}
tb = {}
for run in ["orig42", "core21"]:
    B = famtab(run); B = B[B.region == "Royer_NW_WA"].copy(); B["month"] = B.window.map(mm)
    B = B.merge(A[["fid", "month", "truth_array_frac"]], on=["fid", "month"], how="left")
    B["bin"] = pd.cut(B.truth_array_frac, [-0.01, 0.45, 0.55, 0.65, 1.0], labels=["<= 0.45", "0.45-0.55", "0.55-0.65", "> 0.65"])
    for (w, b), g in B.groupby(["window", "bin"], observed=True):
        no = g.n_ours_cat_core.sum()
        tb.setdefault((w, b), {})[run] = (len(g), g.n_match_own_cat_core.sum() / no if no else np.nan, int(g.pass_both.sum()))
for (w, b), d in sorted(tb.items()):
    L.append(f"| {w} | {b} | {d['core21'][0]} | {f3(d['orig42'][1])} | {f3(d['core21'][1])} | {d['orig42'][2]} | {d['core21'][2]} |")
open(f"{GO}/tables.md", "w").write("\n".join(L) + "\n")
print("\n".join(L))
