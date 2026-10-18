#!/usr/bin/env python
"""GO/NO-GO Step A (Merlin 2026-10-02): truth independence per array-truth month.

For every NW/SW truth month, per WA family (the 211 Royer NW/SW representatives of the T1 setup) that has >= 1 own
array detection in that month (its det_source run), count the family's growth-template stations
(tpl_it<N>_stations.csv next to its det_it<N>.parquet) that are LIVE that month (>= 10 raw station-days with a
non-empty file on disk, any of the 9 raw roots of results/wa_prod_plan/inventory_scan.py, symlinks resolved), split
into ARRAY (temporary networks, wa_common.TEMP_NETS) vs PERMANENT.  truth_array_frac = array / live.

  python scripts/gonogo_truth_frac.py  -> results/wa_gonogo/truth_array_frac.csv, truth_array_frac_months.csv
Read-only on all inputs.  The stage fails (exit 1) if > 3 % of the template station-months cannot be resolved to a
network (station code not found under any raw root) -- those are 'missing', not 'not live'.
"""
import os, sys, collections
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

OD = f"{W.RES}/wa_gonogo"
os.makedirs(OD, exist_ok=True)
MONTHS = ["2007-01", "2007-10", "2008-05", "2009-02", "2009-04", "2009-05", "2010-03", "2010-08", "2011-04",
          "2011-08", "2011-12"]
# lfe_array_discovery/traces_sealed is EXCLUDED: it is that project's sealed holdout (never read, never counted)
RAW_ROOTS = ["lfe_array_discovery/traces", "lfe_array_discovery/traces_eq",
             "lfe_extension/traces", "lfe_catalog_repair/traces", "lfe_extension_wa/traces", "lfe_picker/traces",
             "lfe_south/traces", "data/waveforms"]
MERLIN = {"2007-01": (6, 3), "2007-10": (4, 3), "2008-05": (6, 3), "2009-02": (3, 3), "2009-04": (3, 3),
          "2009-05": (4, 3), "2010-03": (3, 4), "2010-08": (5, 4), "2011-04": (0, 4), "2011-08": (4, 4),
          "2011-12": (3, 4)}
years = sorted({int(m[:4]) for m in MONTHS})

F = pd.read_csv(f"{W.RES}/merlin_tests/t1_setup/families.csv")
F = F[F.region.isin(["Royer_NW_WA", "Royer_SW_WA"])].copy()

# ---- station-code -> netsta, and raw days on disk per (netsta, month)
sta_net = collections.defaultdict(set)
days = collections.defaultdict(set)          # (netsta, 'YYYY-MM') -> set of doy
for r in RAW_ROOTS:
    base = f"{W.ROOT}/{r}"
    if not os.path.isdir(base):
        continue
    for ns in os.listdir(base):
        if "." not in ns or not os.path.isdir(f"{base}/{ns}"):
            continue
        net, sta = ns.split(".", 1)
        if not any(y.isdigit() for y in os.listdir(f"{base}/{ns}")):
            continue                       # event-file layouts (traces_eq) hold no station-days
        sta_net[sta].add(ns)
        for y in years:
            py = f"{base}/{ns}/{y}"
            if not os.path.isdir(py):
                continue
            for e in os.scandir(py):
                if not e.name.endswith(".mseed"):
                    continue
                try:
                    rp = os.path.realpath(e.path)
                    if not (os.path.exists(rp) and os.path.getsize(rp) > 0):
                        continue
                except OSError:
                    continue
                doy = int(e.name.split(".")[0])
                m = (pd.Timestamp(f"{y}-01-01") + pd.Timedelta(days=doy - 1)).strftime("%Y-%m")
                if m in MONTHS:
                    days[(ns, m)].add(doy)

rows, unresolved, total = [], 0, 0
amb = {s: v for s, v in sta_net.items() if len(v) > 1}
for src, g in F.groupby("det_source"):
    stf = f"{W.ROOT}/{src.replace('det_it', 'tpl_it').replace('.parquet', '_stations.csv')}"
    S = pd.read_csv(stf)
    tpl = {int(c): s.split(",") for c, s in zip(S.cid, S.stations)}
    D = pd.read_parquet(f"{W.ROOT}/{src}", columns=["grp", "date"])
    D["m"] = D.date.astype(str).str[:7]
    D = D[D.m.isin(MONTHS)]
    cnt = D.groupby(["grp", "m"]).size()
    for r in g.itertuples():
        stas = tpl.get(int(r.family_id))
        for m in MONTHS:
            nd = int(cnt.get((int(r.family_id), m), 0))
            if nd == 0:
                continue
            if stas is None:
                rows.append(dict(fid=r.fid, family_uid=r.family_uid, region=r.region, det_source=src, month=m,
                                 n_truth=nd, n_tpl_sta=0, note="no template station row"))
                continue
            la = lp = 0; la_l, lp_l, dead_l = [], [], []
            for s in stas:
                total += 1
                nss = sta_net.get(s, set())
                if not nss:
                    unresolved += 1; dead_l.append(f"{s}:unresolved"); continue
                ns = sorted(nss, key=lambda x: -len(days.get((x, m), ())))[0]
                nday = len(days.get((ns, m), ()))
                arr = ns.split(".")[0] in W.TEMP_NETS
                if nday >= 10:
                    if arr:
                        la += 1; la_l.append(f"{ns}:{nday}")
                    else:
                        lp += 1; lp_l.append(f"{ns}:{nday}")
                else:
                    dead_l.append(f"{ns}:{nday}")
            n_live = la + lp
            rows.append(dict(fid=r.fid, family_uid=r.family_uid, region=r.region, det_source=src, month=m,
                             n_truth=nd, n_tpl_sta=len(stas), n_live=n_live, n_array=la, n_perm=lp,
                             truth_array_frac=la / n_live if n_live else np.nan,
                             live_array=";".join(la_l), live_perm=";".join(lp_l), not_live=";".join(dead_l)))
A = pd.DataFrame(rows)
A.to_csv(f"{OD}/truth_array_frac.csv", index=False)
out = []
for (m, sub), g in [((m, "ALL"), A[A.month == m]) for m in MONTHS] + \
        [((m, rg), A[(A.month == m) & (A.region == rg)]) for m in MONTHS for rg in ["Royer_NW_WA", "Royer_SW_WA"]] + \
        [((m, "NW_lag15+SW_lag15"), A[(A.month == m) & A.det_source.str.contains("_lag15/")]) for m in MONTHS]:
    if not len(g):
        continue
    mer = MERLIN[m]
    out.append(dict(month=m, subset=sub, n_families=len(g), n_truth=int(g.n_truth.sum()),
                    med_array=float(g.n_array.median()), med_perm=float(g.n_perm.median()),
                    med_live=float(g.n_live.median()), med_array_frac=float(g.truth_array_frac.median()),
                    mean_array_frac=float(g.truth_array_frac.mean()),
                    frac_fam_array0=float((g.n_array == 0).mean()),
                    merlin_array=mer[0], merlin_perm=mer[1]))
# ---- Merlin's reference set: EVERY grown family (all cids) of NW_lag15 it14 + SW_lag15 it7 with >= 1 detection
# in the month, pooled (not only the WA representatives); reproduces his medians (see GONOGO_RESULTS.md)
def live_counts(stas, m):
    la = lp = 0
    for s in stas:
        nss = sta_net.get(s, set())
        if not nss:
            continue
        ns = sorted(nss, key=lambda x: -len(days.get((x, m), ())))[0]
        if len(days.get((ns, m), ())) >= 10:
            if ns.split(".")[0] in W.TEMP_NETS:
                la += 1
            else:
                lp += 1
    return la, lp
RA = "lfe_array_discovery/results/ab_royer"
for sub, srcs in [("grow_all_NW_lag15", [("NW_lag15", 14)]), ("grow_all_SW_lag15", [("SW_lag15", 7)]),
                  ("grow_all_NW+SW_lag15", [("NW_lag15", 14), ("SW_lag15", 7)])]:
    acc = {m: [] for m in MONTHS}
    for d_, it in srcs:
        S_ = pd.read_csv(f"{W.ROOT}/{RA}/{d_}/grow/tpl_it{it}_stations.csv")
        D_ = pd.read_parquet(f"{W.ROOT}/{RA}/{d_}/grow/det_it{it}.parquet", columns=["grp", "date"])
        D_["m"] = D_.date.astype(str).str[:7]
        tp_ = {int(c): s.split(",") for c, s in zip(S_.cid, S_.stations)}
        for m in MONTHS:
            for c in sorted(set(D_[D_.m == m].grp)):
                if int(c) in tp_:
                    acc[m].append(live_counts(tp_[int(c)], m))
    for m in MONTHS:
        if not acc[m]:
            continue
        x = np.array(acc[m], float); lv = x.sum(1)
        fr = np.where(lv > 0, x[:, 0] / np.maximum(lv, 1), np.nan)
        out.append(dict(month=m, subset=sub, n_families=len(x), n_truth=np.nan, med_array=float(np.median(x[:, 0])),
                        med_perm=float(np.median(x[:, 1])), med_live=float(np.median(lv)),
                        med_array_frac=float(np.nanmedian(fr)), mean_array_frac=float(np.nanmean(fr)),
                        frac_fam_array0=float((x[:, 0] == 0).mean()), merlin_array=MERLIN[m][0],
                        merlin_perm=MERLIN[m][1]))
M = pd.DataFrame(out)
M.to_csv(f"{OD}/truth_array_frac_months.csv", index=False)
print(M.to_string())
print(f"template station-months {total}, unresolved station codes {unresolved} ({unresolved/max(total,1):.2%}); "
      f"ambiguous codes {amb}")
if unresolved / max(total, 1) > 0.03:
    print("FAILED: > 3 % of template station-months unresolved"); sys.exit(1)
