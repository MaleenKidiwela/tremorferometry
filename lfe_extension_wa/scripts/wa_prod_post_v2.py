#!/usr/bin/env python
"""WA PRODUCTION post-processing v2 (user decisions 2026-10-06).  Reads the v1 inputs (conf/P/*, conf/truth_rows.parquet,
setup/*, post/tremor_family_day.parquet), writes ONLY under results/wa_prod/post_v2/ (log post_v2/post_v2.log).  Never
writes v1 files, never touches queue markers, CPU only.

Day classes v2 (per family and UTC day):
  tremor   >= 1 PNSN tremor epicentre within 30 km of the family position that day (catalogue days), or, before the
           catalogue starts (2009-08-06), a day inside the WECC ETS windows (2007-01-12..02-02, 2008-04-30..06-02,
           2009-04-28..06-03; lfe_array_discovery/results/ab_royer/POOL_RECONSTRUCTION.md sec. 2)
  quiet    the catalogue covers the day (2009-08-06 .. last catalogue day, 2026-05-27) and has no epicentre within 30 km
  unknown  every other day (before 2009-08-06 outside the WECC windows; after the catalogue's last day)
HEADLINE rule (user, final 2026-10-06; columns selected_v2 / cut_v2 / rule_v2, rates *_v2):
  NW and Arm B: P_own_WAVE_PNET >= 0.75 on every day class (tremor, quiet, unknown);
  SW: P_own_WAVE_PNET >= 0.5 on tremor days only (quiet / unknown: no headline single repeats, columns only);
  Olympic A: position-group rates on P_nbr3_WAVE >= 0.75 (no PNET version exists), every day class.
  'P_own_WAVE_PNET' = the saved PNET-adopted model (its out-of-fold truth-row column is P_own_WITH_PNET).
SECONDARY 'class rules' (the earlier decisions of the same day; columns *_classrules, rates *_v2cr), with either model:
  tremor days P >= 0.75 (SW 0.5); NW quiet days P >= 0.85; unknown days and SW / Arm B / Olympic A quiet days: no rule.
  'wave' = P_own_WAVE, 'pnet' = P_own_WAVE_PNET (Olympic A: P_nbr3_WAVE in both).
Counts at 0.5 / 0.75 / 0.85 / 0.9 are kept as columns everywhere.

  python scripts/wa_prod_post_v2.py all        (pass over the P tables -> cutoff + catalogue + unit-day counts; k; rates; steps)
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, warnings
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P
import wa_prod_post as V1          # import only (no side effects): reuse core_changes / _merged
W = P.W

O2 = os.environ.get("WA_POST_V2_DIR", f"{P.OUT}/post_v2")
LOG2 = f"{O2}/post_v2.log"
# fixed-10 re-run options (2026-10-07; every default = the production behaviour, unchanged):
#   WA_POST_NET       network whose rows make the catalogue / rates (default core)
#   WA_POST_HEAD_LAB  headline family label: own_pnet (default, P_own_WAVE_PNET) or own (P_own_WAVE, no PNET model)
#   WA_POST_KNET      '' (default): k per era from the truth rows of the era sets E1 / E2 / E3; else k from the truth rows of
#                     this network, applied to every era (one fixed network)
#   WA_POST_SCANNED   tpl (default): unit scanned when >= 1 template node is live in the core run; famday: unit scanned when
#                     the family is observable (>= 2 sites) in WA_POST_NET that day (mf famday tables)
NET = os.environ.get("WA_POST_NET", "core")
HEAD_LAB = os.environ.get("WA_POST_HEAD_LAB", "own_pnet")
KNET = os.environ.get("WA_POST_KNET", "")
SCANNED = os.environ.get("WA_POST_SCANNED", "tpl")
HEAD_NAME = "PNET" if HEAD_LAB == "own_pnet" else "WAVE"
WECC = V1.WECC
LAB_REG = {"NW": "own", "SW": "own", "ArmB": "own", "OlyA": "nbr3"}
CUTS = [0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9]
RCUTS = {"050": 0.5, "075": 0.75, "085": 0.85, "090": 0.9}            # rate columns
FULL = f"{P.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv"
P2009 = f"{P.ROOT}/lfe_array_discovery/results/ab_royer/tremor/pnsn_tremor_2009.csv"
LABCOL = {"own": "P_own_WAVE", "own_pnet": "P_own_WAVE_PNET", "nbr3": "P_nbr3_WAVE", "any": "P_any_WAVE"}
VARIANT_LAB = {"wave": "own", "pnet": "own_pnet"}
CLASSES = ["all_days", "tremor", "quiet", "unknown", "quiet_v1", "quiet_beyond60", "unknown_pre2009", "unknown_post_catalogue"]


def log(m):
    os.makedirs(O2, exist_ok=True)
    s = f"{pd.Timestamp.now('UTC').strftime('%Y-%m-%dT%H:%M:%SZ')} {m}"
    print(s, flush=True)
    with open(LOG2, "a") as f:
        f.write(s + "\n")


# ------------------------------------------------------------------------------------------------ catalogue span + classes
def catalogue_span():
    p = f"{O2}/catalogue_span.json"
    if os.path.exists(p):
        return json.load(open(p))
    T = pd.read_csv(FULL, usecols=["time"])
    d = T.time.astype(str).str[:10]
    T9 = pd.read_csv(P2009, skipinitialspace=True)
    T9.columns = [c.strip() for c in T9.columns]
    d9 = T9.starttime.astype(str).str.strip().str[:10]
    J = dict(full_first=d.min(), full_last=d.max(), pull2009_first=d9.min(), pull2009_last=d9.max(),
             covered_from=min(d9.min(), d.min()), covered_to=d.max(),
             note="covered = covered_from .. covered_to inclusive (the 2009 API pull starts 2009-08-06; the full file runs "
                  "2010-01-06 .. its last day). The 6 days 2009-12-31 .. 2010-01-05 lie between the two files and hold no "
                  "epicentre; they are treated as covered (catalogue-quiet), as the user's span says.")
    os.makedirs(O2, exist_ok=True)
    json.dump(J, open(p, "w"), indent=1)
    return J


_X = None


def _tfd():
    global _X
    if _X is None:
        X = pd.read_parquet(f"{P.OUT}/post/tremor_family_day.parquet")
        _X = X.set_index(["fid", "date"])
    return _X


def day_class_v2(D, span):
    """adds day_class (tremor / quiet / unknown), n_trem30, n_trem60 (NaN off catalogue), quiet_v1 (the v1 definition:
    no tremor <= 60 km, uncovered days counted as quiet), unknown_kind, cut_v2 (NaN = no validated rule)."""
    Xi = _tfd()
    k = pd.MultiIndex.from_arrays([D.fid.astype(int).to_numpy(), D.date.to_numpy()])
    Y = Xi.reindex(k)
    n30 = Y["n30"].fillna(0).to_numpy() if "n30" in Y else np.zeros(len(D))
    n60 = Y["n60"].fillna(0).to_numpy() if "n60" in Y else np.zeros(len(D))
    dt = D.date.to_numpy().astype(str)
    wecc = np.zeros(len(D), bool)
    for a0, b0 in WECC:
        wecc |= (dt >= a0) & (dt <= b0)
    pre = dt < span["covered_from"]
    post = dt > span["covered_to"]
    cov = ~pre & ~post
    cls = np.where(cov, np.where(n30 > 0, "tremor", "quiet"), np.where(pre & wecc, "tremor", "unknown"))
    D["day_class"] = cls
    D["n_trem30"] = np.where(cov, n30, np.nan)
    D["n_trem60"] = np.where(cov, n60, np.nan)
    D["quiet_v1"] = np.where(pre, ~wecc, n60 == 0)
    D["unknown_kind"] = np.where(cls != "unknown", "", np.where(pre, "pre2009", "post_catalogue"))
    return D


def cut_rule(region, cls):
    """vectorised v2 rule -> cut (NaN = unvalidated)."""
    region = np.asarray(region).astype(str); cls = np.asarray(cls).astype(str)
    c = np.full(len(cls), np.nan)
    c[cls == "tremor"] = 0.75
    c[(cls == "tremor") & (region == "SW")] = 0.5
    c[(cls == "quiet") & (region == "NW")] = 0.85
    return c


def cut_headline(region, cls):
    """vectorised HEADLINE rule -> cut (NaN = no headline single repeats)."""
    region = np.asarray(region).astype(str); cls = np.asarray(cls).astype(str)
    c = np.full(len(cls), 0.75)
    c[region == "SW"] = np.nan
    c[(region == "SW") & (cls == "tremor")] = 0.5
    return c


def class_masks(D):
    c = D.day_class.to_numpy()
    return {"all_days": np.ones(len(D), bool), "tremor": c == "tremor", "quiet": c == "quiet", "unknown": c == "unknown",
            "quiet_v1": D.quiet_v1.to_numpy().astype(bool),
            "quiet_beyond60": (c == "quiet") & (D.n_trem60.to_numpy() == 0),
            "unknown_pre2009": D.unknown_kind.to_numpy() == "pre2009",
            "unknown_post_catalogue": D.unknown_kind.to_numpy() == "post_catalogue"}


# ------------------------------------------------------------------------------------------------ unit-day counts
def unit_day_counts_v2(D):
    """per (unit_kind, unit, label, date, series): n050 / n075 / n085 / n090 (Olympic A groups: members merged at +-2 s) + Sum P.
    Family labels own (P_own_WAVE), own_pnet (P_own_WAVE_PNET, if present), nbr3; Olympic A groups nbr3."""
    F = P.families()
    grp = dict(zip(F.fid, F.pos_group.fillna("")))
    olyA = D.fid.map(lambda f: grp.get(f, "") != "").to_numpy()
    out = []
    Fm = D[~olyA]
    for lab in [x for x in ("own", "own_pnet", "nbr3") if LABCOL[x] in D]:
        pc = LABCOL[lab]
        cols = [Fm.groupby(["fid", "date", "series"])[pc].sum().rename("sumP")]
        for tg, c in RCUTS.items():
            cols.append(Fm[Fm[pc] >= c].groupby(["fid", "date", "series"]).size().rename(f"n{tg}"))
        X = pd.concat(cols, axis=1).fillna(0).reset_index()
        X["unit"] = X.fid.astype(int).astype(str); X["unit_kind"] = "family"; X["label"] = lab
        out.append(X.drop(columns=["fid"]))
    Om = D[olyA].copy()
    if len(Om):
        Om["unit"] = Om.fid.map(grp)
        cols = [Om.groupby(["unit", "date", "series"]).P_nbr3_WAVE.sum().rename("sumP")]
        for tg, c in RCUTS.items():
            hi = Om[Om.P_nbr3_WAVE >= c]
            cols.append(hi.groupby(["unit", "date", "series"]).t.apply(V1._merged).rename(f"n{tg}"))
        X = pd.concat(cols, axis=1).fillna(0).reset_index()
        X["unit_kind"] = "olyA_group"; X["label"] = "nbr3"
        out.append(X)
    return pd.concat(out, ignore_index=True)


# ------------------------------------------------------------------------------------------------ pass over P tables
CAT_SCHEMA = pa.schema([("det_id", pa.int64()), ("time", pa.timestamp("ns")), ("date", pa.string()), ("family", pa.int32()),
                        ("family_uid", pa.string()), ("region", pa.string()), ("pos_group", pa.string()),
                        ("P_own", pa.float32()), ("P_own_pnet", pa.float32()), ("P_nbr3", pa.float32()), ("P_any", pa.float32()),
                        ("P_max", pa.float32()), ("net_mad", pa.float32()), ("n_sta", pa.float32()), ("role", pa.string()),
                        ("day_class", pa.string()), ("unknown_kind", pa.string()), ("n_trem30", pa.float32()),
                        ("n_trem60", pa.float32()), ("label_v2", pa.string()), ("cut_v2", pa.float32()),
                        ("rule_v2", pa.string()), ("selected_v2", pa.bool_()),
                        ("cut_v2_classrules", pa.float32()), ("rule_v2_classrules", pa.string()),
                        ("selected_v2_wave_classrules", pa.bool_()), ("selected_v2_pnet_classrules", pa.bool_())])


def scan(a):
    span = catalogue_span()
    F = P.families()
    uid = dict(zip(F.fid.astype(int), F.family_uid.astype(str)))
    pg = dict(zip(F.fid.astype(int), F.pos_group.fillna("")))
    acc = {}; fams = {}; ndays = {}
    cat_path = f"{O2}/catalogue_v2.parquet"
    tmp = cat_path + ".partial"
    writer = pq.ParquetWriter(tmp, CAT_SCHEMA, compression="zstd")
    U_parts = []
    cols = ["date", "net", "series", "fid", "k", "t", "reg", "role", "net_mad", "n_sta", "P_own_WAVE", "P_nbr3_WAVE",
            "P_any_WAVE", "P_max_WAVE", "P_own_WAVE_PNET"]
    n_cat = 0
    for run in ("prod", "rich"):
        for p in sorted(glob.glob(f"{P.OUT}/conf/P/{run}/P_*.parquet")):
            have = set(pq.read_schema(p).names)
            D = pd.read_parquet(p, columns=[c for c in cols if c in have] if run == "prod" else [c for c in cols if c not in ("k", "net_mad", "n_sta", "P_max_WAVE", "P_own_WAVE_PNET") and c in have])
            if not len(D):
                continue
            D = day_class_v2(D, span)
            # ---- cut-off table (sealed months excluded, as v1)
            Dn = D[~D.role.astype(str).str.startswith("sealed")]
            M = class_masks(Dn)
            for cls in CLASSES:
                h = Dn[M[cls]]
                if not len(h):
                    continue
                for (net, reg), g in h.groupby(["net", "reg"]):
                    ndays.setdefault((run, net, reg, cls), set()).update(g.date.unique())
                    for lab in [x for x in ("own", "own_pnet", "nbr3", "any") if LABCOL[x] in g]:
                        pcv = g[LABCOL[lab]].to_numpy()
                        ser = g.series.to_numpy(); fid = g.fid.to_numpy()
                        for c in CUTS:
                            s = pcv >= c
                            key = (run, net, reg, lab, cls, c)
                            v = acc.setdefault(key, [0, 0, 0])
                            v[0] += int((s & (ser == "data")).sum()); v[1] += int((s & (ser == "rev")).sum()); v[2] += int((s & (ser == "c1")).sum())
                            fams.setdefault(key, set()).update(np.unique(fid[s & (ser == "data")]).tolist())
            if run != "prod":
                continue
            # ---- catalogue (core network, data series, every candidate incl. sealed months)
            C = D[(D.net == NET) & (D.series == "data")]
            reg = C.reg.astype(str).to_numpy()
            lab = np.where(reg == "OlyA", "nbr3", HEAD_LAB)
            ppn = C.P_own_WAVE_PNET.to_numpy() if "P_own_WAVE_PNET" in C else np.full(len(C), np.nan)
            hcut = cut_headline(reg, C.day_class.to_numpy())
            cut = cut_rule(reg, C.day_class.to_numpy())
            cc = np.nan_to_num(cut, nan=9.0)
            pl_w = np.where(reg != "OlyA", C.P_own_WAVE.to_numpy(), C.P_nbr3_WAVE.to_numpy())
            pl_p = np.where(reg != "OlyA", C[LABCOL[HEAD_LAB]].to_numpy(), C.P_nbr3_WAVE.to_numpy())
            sel_h = np.isfinite(hcut) & (pl_p >= np.nan_to_num(hcut, nan=9.0))
            hrule = np.where(np.isnan(hcut), "SW_no_headline_off_tremor_days",
                             np.where(reg == "OlyA", "OlyA_group_nbr3_WAVE>=0.75",
                                      np.where(reg == "SW", f"SW_tremor_{HEAD_NAME}>=0.5", f"{HEAD_NAME}>=0.75_all_days")))
            sel_w = np.isfinite(cut) & (pl_w >= cc)
            sel_p = np.isfinite(cut) & (pl_p >= cc)
            rule = np.where(cut == 0.5, "SW_tremor_P>=0.5", np.where(cut == 0.75, "tremor_P>=0.75",
                            np.where(cut == 0.85, "NW_quiet_P>=0.85", "unvalidated_no_rule")))
            ymd = C.date.str.replace("-", "", regex=False).astype(np.int64).to_numpy()
            tb = pa.table({
                "det_id": pa.array(ymd * 10_000_000 + C.k.to_numpy().astype(np.int64)),
                "time": pa.array((C.t.to_numpy() * 1e9).round().astype("int64").astype("datetime64[ns]")),
                "date": pa.array(C.date.astype(str).to_numpy()),
                "family": pa.array(C.fid.astype(np.int32).to_numpy()),
                "family_uid": pa.array([uid.get(int(f), "") for f in C.fid]),
                "region": pa.array(reg),
                "pos_group": pa.array([pg.get(int(f), "") for f in C.fid]),
                "P_own": pa.array(C.P_own_WAVE.to_numpy(np.float32)), "P_nbr3": pa.array(C.P_nbr3_WAVE.to_numpy(np.float32)),
                "P_own_pnet": pa.array(ppn.astype(np.float32)),
                "P_any": pa.array(C.P_any_WAVE.to_numpy(np.float32)), "P_max": pa.array(C.P_max_WAVE.to_numpy(np.float32)),
                "net_mad": pa.array(C.net_mad.to_numpy(np.float32)), "n_sta": pa.array(C.n_sta.to_numpy(np.float32)),
                "role": pa.array(C.role.astype(str).to_numpy()), "day_class": pa.array(C.day_class.to_numpy().astype(str)),
                "unknown_kind": pa.array(C.unknown_kind.to_numpy().astype(str)),
                "n_trem30": pa.array(C.n_trem30.to_numpy(np.float32)), "n_trem60": pa.array(C.n_trem60.to_numpy(np.float32)),
                "label_v2": pa.array(lab), "cut_v2": pa.array(hcut.astype(np.float32)), "rule_v2": pa.array(hrule),
                "selected_v2": pa.array(sel_h),
                "cut_v2_classrules": pa.array(cut.astype(np.float32)), "rule_v2_classrules": pa.array(rule),
                "selected_v2_wave_classrules": pa.array(sel_w), "selected_v2_pnet_classrules": pa.array(sel_p)}, schema=CAT_SCHEMA)
            writer.write_table(tb, row_group_size=2_000_000)
            n_cat += len(C)
            # ---- unit-day counts for rates (core, data / rev / c1)
            R = D[(D.net == NET) & D.series.isin(["data", "rev", "c1"])]
            if len(R):
                U_parts.append(unit_day_counts_v2(R))
            log(f"[scan] {run} {os.path.basename(p)} catalogue rows so far {n_cat:,}")
        log(f"[scan] {run} done")
    writer.close()
    os.replace(tmp, cat_path)
    U = pd.concat(U_parts, ignore_index=True)
    U.to_parquet(f"{O2}/unit_day_counts_v2.parquet", index=False)
    rows = []
    for (run, net, reg, lab, cls, c), (nd, nr, n1) in acc.items():
        rows.append(dict(run=run, net=net, region=reg, label=lab, day_class=cls, cut=c, n_days=len(ndays.get((run, net, reg, cls), ())),
                         data=nd, rev_control=nr, scrambled_control=n1, false_share_rev=nr / nd if nd else np.nan,
                         false_share_scrambled=n1 / nd if nd else np.nan, families_covered=len(fams[(run, net, reg, lab, cls, c)])))
    Rt = pd.DataFrame(rows).sort_values(["run", "net", "region", "label", "day_class", "cut"])
    oly = Rt.region == "OlyA"
    head_w = np.where(oly, Rt.label == "nbr3", Rt.label == "own")
    head_p = np.where(oly, Rt.label == "nbr3", Rt.label == "own_pnet")
    rule = np.where(Rt.day_class == "tremor", np.where(Rt.region == "SW", 0.5, 0.75),
                    np.where((Rt.day_class == "quiet") & (Rt.region == "NW"), 0.85, np.nan))
    hrule = np.where(Rt.region == "SW", np.where(Rt.day_class == "tremor", 0.5, np.nan), 0.75)
    hrule = np.where(Rt.day_class.isin(["tremor", "quiet", "unknown"]), hrule, np.nan)
    Rt["headline_cut"] = hrule
    Rt["is_headline"] = head_p & np.isclose(Rt.cut, hrule)
    Rt["classrules_cut"] = rule
    Rt["is_classrules_wave"] = head_w & np.isclose(Rt.cut, rule)
    Rt["is_classrules_pnet"] = head_p & np.isclose(Rt.cut, rule)
    Rt["variant"] = np.where(oly & (Rt.label == "nbr3"), "wave+pnet (Olympic A: nbr3 WAVE in both)",
                             np.where(Rt.label == "own", "wave", np.where(Rt.label == "own_pnet", "pnet", "")))
    Rt["class_kind"] = np.where(Rt.day_class.isin(["tremor", "quiet", "unknown"]), "v2 class", "diagnostic subset")
    Rt.to_csv(f"{O2}/false_share_by_cutoff_v2.csv", index=False)
    log(f"[scan] catalogue {n_cat:,} rows -> {cat_path}; cut-off table {len(Rt)} rows")



# ------------------------------------------------------------------------------------------------ k factors per cut and label
def k_factors_v2():
    """k per (era set, cut, label, unit): truth events / (data - time-reversed control) on the truth fit days.
    Labels: family own (P_own_WAVE), family own_pnet (truth column P_own_WITH_PNET = out-of-fold WAVE+PNET), Olympic A
    group nbr3."""
    tp = f"{P.OUT}/conf/truth_rows.parquet"
    kc = ["date", "net", "series", "fid", "t", "reg", "row_role", "P_own_WAVE", "P_nbr3_WAVE", "P_own_WITH_PNET"]
    A = pd.read_parquet(tp, columns=[c for c in kc if c in set(pq.read_schema(tp).names)])
    A = A[(A.row_role == "fit") & A.net.isin([KNET] if KNET else ["E1", "E2", "E3"])].rename(columns={"P_own_WITH_PNET": "P_own_WAVE_PNET"})
    T = pd.read_parquet(f"{P.OUT}/setup/truth.parquet")
    T = T[T.role == "fit"]
    # truth events only in the (day, build/test half) blocks whose detections are fit rows: on V2008 / V2010 days only the
    # test half (1) is 'fit' (build half 0 = insample, edges -1 = guard), but truth.parquet keeps every half -> k was ~2x
    # too high there (Merlin 2026-10-06; fixed same day)
    import wa_common as W
    fit_dh = set(zip(A.date.astype(str), W.half_of(A.t.to_numpy())))
    T = T[[(str(d), int(h)) in fit_dh for d, h in zip(T.date, W.half_of(T.t.to_numpy()))]]
    SC = pd.read_csv(f"{P.OUT}/setup/scan_days.csv")
    F = P.families()
    src = dict(zip(F.fid, F.det_source)); grp = dict(zip(F.fid, F.pos_group.fillna(""))); regm = dict(zip(F.fid, F.reg))
    scanned = set(zip(SC.det_source, SC.date))
    T["unit"] = [grp[f] if grp.get(f, "") else str(int(f)) for f in T.fid]
    T["unit_kind"] = ["olyA_group" if grp.get(f, "") else "family" for f in T.fid]
    TR = T.groupby(["unit_kind", "unit", "date"]).t.apply(lambda x: V1._merged(x)).rename("n_truth").reset_index()
    out = []
    for net, An in A.groupby("net"):
        U = unit_day_counts_v2(An)
        U = U[((U.unit_kind == "family") & U.label.isin(["own", "own_pnet"])) | (U.unit_kind == "olyA_group")]
        idx = ["unit_kind", "unit", "label", "date"]
        pv_s = U.pivot_table(index=idx, columns="series", values="sumP", fill_value=0.0)
        for c in ("data", "rev"):
            if c not in pv_s:
                pv_s[c] = 0.0
        net_s = (pv_s["data"] - pv_s["rev"])
        days_net = sorted(An.date.unique())
        for tg, cval in RCUTS.items():
            piv = U.pivot_table(index=idx, columns="series", values=f"n{tg}", fill_value=0)
            for c in ("data", "rev"):
                if c not in piv:
                    piv[c] = 0
            piv["net_c"] = piv["data"] - piv["rev"]
            piv["net_s"] = net_s.reindex(piv.index).fillna(0.0)
            piv = piv.reset_index()[idx + ["net_c", "net_s"]]
            for (uk, u), g in TR.groupby(["unit_kind", "unit"]):
                f0 = int(u) if uk == "family" else int(F[F.pos_group == u].fid.iloc[0])
                sd = [d for d in days_net if (src[f0], d) in scanned]
                if not sd:
                    continue
                tr_d = g.set_index("date").n_truth.reindex(sd).fillna(0).to_numpy()
                for lab in ((("own", "own_pnet") if "P_own_WAVE_PNET" in A else ("own",)) if uk == "family" else ("nbr3",)):
                    rs = np.random.default_rng(7)                 # same day draws for every label / cut
                    pu = piv[(piv.unit_kind == uk) & (piv.unit == u) & (piv.label == lab)].set_index("date")
                    pv = pu.net_c.reindex(sd).fillna(0).to_numpy()
                    ps_ = pu.net_s.reindex(sd).fillna(0).to_numpy()
                    n_tr, net_c, nets = float(tr_d.sum()), float(pv.sum()), float(ps_.sum())
                    bs = []
                    for _ in range(100):
                        ii = rs.integers(0, len(sd), len(sd))
                        bs.append((float(tr_d[ii].sum()), float(pv[ii].sum())))
                    out.append(dict(net=net, cut=cval, cut_tag=tg, label=lab, unit_kind=uk, unit=u, reg=regm[f0], n_truth=n_tr,
                                    net_count=net_c, n_days=len(sd), net_sumP=nets, k_sumP_raw=n_tr / nets if nets > 0 else np.nan,
                                    k_raw=n_tr / net_c if net_c > 0 else np.nan, boot_pairs=bs))
    K = pd.DataFrame(out)
    n0 = 20.0
    gk = ["net", "cut_tag", "label", "reg"]
    K["logk"] = np.log(K.k_raw.where(K.k_raw > 0))
    med = K.groupby(gk).logk.median()

    def shrunk(t_, c_, m_):
        lk = np.log(t_ / c_) if (c_ > 0 and t_ > 0) else m_
        return float(np.exp((max(c_, 0) * lk + n0 * m_) / (max(c_, 0) + n0)))
    K["k"] = [shrunk(r.n_truth, r.net_count, med.get((r.net, r.cut_tag, r.label, r.reg), np.nan)) for r in K.itertuples()]
    lo, hi = [], []
    for r in K.itertuples():
        v = [shrunk(t_, c_, med.get((r.net, r.cut_tag, r.label, r.reg), np.nan)) for t_, c_ in r.boot_pairs]
        v = [x for x in v if np.isfinite(x)]
        lo.append(float(np.percentile(v, 5)) if len(v) >= 20 else np.nan); hi.append(float(np.percentile(v, 95)) if len(v) >= 20 else np.nan)
    K["k_boot_lo"] = lo; K["k_boot_hi"] = hi
    K["logks"] = np.log(K.k_sumP_raw.where(K.k_sumP_raw > 0))
    meds = K.groupby(gk).logks.median()
    K["k_sumP"] = [np.exp((max(r.net_sumP, 0) * (r.logks if np.isfinite(r.logks) else meds.get((r.net, r.cut_tag, r.label, r.reg), np.nan))
                           + n0 * meds.get((r.net, r.cut_tag, r.label, r.reg), np.nan)) / (max(r.net_sumP, 0) + n0)) for r in K.itertuples()]
    K = K.drop(columns=["boot_pairs"])
    # the E1 / E2 truth run (truth_e12) has no PNET block: its P_own_WITH_PNET is the PNET model fed missing PNET features,
    # so a PNET k for E1 / E2 cannot be measured (0 % of those fit rows carry pnet_max; E3: 100 %)
    bad = (K.label == "own_pnet") & K.net.isin(["E1", "E2"])
    K["k_status"] = np.where(bad, "unavailable: no PNET block on the E1/E2 truth run", "ok")
    K.loc[bad, ["k", "k_boot_lo", "k_boot_hi", "k_sumP"]] = np.nan
    K["shrink_note"] = ("k = truth events / (data - time-reversed control) at the cut, truth fit days; log k shrunk to the "
                        "(era set, cut, label, region) median with pseudo-count 20 net detections; 90 % day bootstrap (100 draws)")
    return K




# ------------------------------------------------------------------------------------------------ rates
# views: (name, family label, prefix of the rule columns)
VIEWS = [("headline", "own_pnet", "v2"), ("wave_classrules", "own", "v2cr"), ("pnet_classrules", "own_pnet", "v2cr")]
if HEAD_LAB == "own":                     # WAVE-only run (no PNET model): the headline view is the WAVE model
    VIEWS = [("headline", "own", "v2"), ("wave_classrules", "own", "v2cr")]


def view_rows(R, name):
    lab = {n: l for n, l, _ in VIEWS}[name]
    return R[((R.unit_kind == "family") & (R.label == lab)) | (R.unit_kind == "olyA_group")].assign(view=name)


def scanned_unit_days(units):
    """(unit_kind, unit, date) for every production day on which the unit was scanned: at least one of its template nodes
    (setup/tpl_library.csv) was live in that day's core run (mf/prod/run_<date>.json members.core -> tkeys). Olympic A
    groups: any member family scanned. Used to add zero rows: rates_daily_v2 used to hold only unit-days with >= 1 candidate,
    so scanned days with no candidate were missing (NW E3 tremor family-days 10,691 of 57,831; all had 16-18 live template
    stations) and per-unit-day rates were ~1.2x too high (found 2026-10-06, episode check; fixed same day)."""
    F = P.families()
    if SCANNED == "famday":
        return scanned_unit_days_famday(units, F)
    T = pd.read_csv(f"{P.OUT}/setup/tpl_library.csv", usecols=["fid", "node"])
    fam_nodes = T.groupby("fid").node.apply(frozenset).to_dict()
    grp = dict(zip(F.fid, F.pos_group.fillna("")))
    members = {}
    for f, g in grp.items():
        if g:
            members.setdefault(g, []).append(f)
    out = []
    for path in sorted(glob.glob(f"{P.OUT}/mf/prod/run_*.json")):
        d = json.load(open(path))
        live = {d["tkeys"][m] for m in d["members"].get("core", []) if m in d["tkeys"]}
        if not live:
            continue
        ok = {f for f, nd in fam_nodes.items() if nd & live}
        for uk, u in units:
            if uk == "family":
                hit = int(u) in ok
            else:
                hit = any(f in ok for f in members.get(u, []))
            if hit:
                out.append((uk, u, d["date"]))
    return pd.DataFrame(out, columns=["unit_kind", "unit", "date"])


def scanned_unit_days_famday(units, F):
    """WA_POST_SCANNED=famday: a family is scanned on a day when it is observable (>= 2 sites) in network NET that day
    (mf/prod/famday_<date>.parquet); an Olympic A group when any member family is."""
    grp = dict(zip(F.fid, F.pos_group.fillna("")))
    members = {}
    for f, g in grp.items():
        if g:
            members.setdefault(g, []).append(f)
    out = []
    for path in sorted(glob.glob(f"{P.OUT}/mf/prod/famday_*.parquet")):
        X = pd.read_parquet(path, columns=["net", "fid", "date", "observable"])
        X = X[(X.net == NET) & X.observable.astype(bool)]
        if not len(X):
            continue
        ok = set(X.fid.astype(int)); d = str(X.date.iloc[0])
        for uk, u in units:
            hit = (int(u) in ok) if uk == "family" else any(f in ok for f in members.get(u, []))
            if hit:
                out.append((uk, u, d))
    return pd.DataFrame(out, columns=["unit_kind", "unit", "date"])


def rates(a):
    span = catalogue_span()
    K = k_factors_v2()
    K.to_csv(f"{O2}/rates_k_factors_v2.csv", index=False)
    log(f"[rates] k factors {len(K)} rows")
    U = pd.read_parquet(f"{O2}/unit_day_counts_v2.parquet")
    vals = ["sumP"] + [f"n{t}" for t in RCUTS]
    R = U.pivot_table(index=["unit_kind", "unit", "label", "date"], columns="series", values=vals, fill_value=0)
    R.columns = [f"{a_}_{b_}" for a_, b_ in R.columns]
    R = R.reset_index()
    for v in vals:
        for s in ("data", "rev", "c1"):
            if f"{v}_{s}" not in R:
                R[f"{v}_{s}"] = 0.0
    # zero rows for scanned unit-days without any candidate (see scanned_unit_days)
    labs = R.groupby("unit_kind").label.unique().to_dict()
    SD = scanned_unit_days(sorted(set(zip(R.unit_kind, R.unit))))
    Z = pd.concat([SD[SD.unit_kind == uk].assign(label=lab) for uk, ls in labs.items() for lab in ls], ignore_index=True)
    key = ["unit_kind", "unit", "label", "date"]
    Z = Z.merge(R[key], on=key, how="left", indicator=True)
    Z = Z[Z._merge == "left_only"].drop(columns="_merge")
    for c in R.columns:
        if c not in key:
            Z[c] = 0.0
    log(f"[rates] zero rows added for scanned unit-days without candidates: {len(Z)} (existing rows {len(R)})")
    R = pd.concat([R, Z[R.columns]], ignore_index=True)
    F = P.families()
    first = F[F.pos_group.fillna("") != ""].groupby("pos_group").fid.first().to_dict()
    R["fid"] = [int(u) if uk == "family" else int(first[u]) for uk, u in zip(R.unit_kind, R.unit)]
    R["region"] = R.fid.map(dict(zip(F.fid, F.reg)))
    R["era"] = R.date.map(P.era_of)
    R = day_class_v2(R, span)
    grpu = (R.unit_kind == "olyA_group").to_numpy()
    head = grpu | R.label.isin(["own", "own_pnet"]).to_numpy()
    R["model"] = np.where(grpu, "nbr3_WAVE (Olympic A group)", np.where(R.label == "own", "WAVE", np.where(R.label == "own_pnet", "WAVE+PNET", "nbr3_WAVE")))
    kidx = pd.MultiIndex.from_arrays([R.era if not KNET else pd.Series(KNET, index=R.index), R.label, R.unit_kind, R.unit])
    for tg in RCUTS:
        R[f"net{tg}"] = R[f"n{tg}_data"] - R[f"n{tg}_rev"]
        R[f"net{tg}_scr"] = R[f"n{tg}_data"] - R[f"n{tg}_c1"]
        Ki = K[K.cut_tag == tg].set_index(["net", "label", "unit_kind", "unit"])[["k", "k_boot_lo", "k_boot_hi"]]
        kk = Ki.reindex(kidx)
        R[f"k{tg}"] = np.where(head, kk.k.to_numpy(), np.nan)
        R[f"k{tg}_lo"] = np.where(head, kk.k_boot_lo.to_numpy(), np.nan)
        R[f"k{tg}_hi"] = np.where(head, kk.k_boot_hi.to_numpy(), np.nan)
        R[f"rate{tg}"] = R[f"k{tg}"] * R[f"net{tg}"]
    Ks = K[K.cut_tag == "075"].set_index(["net", "label", "unit_kind", "unit"]).k_sumP
    R["k_sumP"] = np.where(head, Ks.reindex(kidx).to_numpy(), np.nan)
    R["net_sumP"] = R.sumP_data - R.sumP_rev
    R["net_sumP_scr"] = R.sumP_data - R.sumP_c1
    R["sumP_scaled"] = R.k_sumP * R.net_sumP                      # second column (all classes; filter with day_class)
    reg_, cls_ = R.region.to_numpy(), R.day_class.to_numpy()
    ch = cut_headline(reg_, cls_); ch[~(grpu | (R.label == HEAD_LAB).to_numpy())] = np.nan
    cr = cut_rule(reg_, cls_); cr[~head] = np.nan
    for pre, cut in (("v2", ch), ("v2cr", cr)):
        R[f"cut_{pre}"] = cut
        tagv = np.select([cut == 0.5, cut == 0.75, cut == 0.85], ["050", "075", "085"], "")

        def pick(fmt):
            out = np.full(len(R), np.nan)
            for tg in ("050", "075", "085"):
                m = tagv == tg
                out[m] = R[fmt.format(tg)].to_numpy()[m]
            return out
        for s in ("data", "rev", "c1"):
            R[f"n_{pre}_{s}"] = pick("n{}_" + s)
        R[f"net_count_{pre}"] = R[f"n_{pre}_data"] - R[f"n_{pre}_rev"]        # main: minus the time-reversed control
        R[f"net_count_{pre}_scr"] = R[f"n_{pre}_data"] - R[f"n_{pre}_c1"]     # second: minus the scrambled control
        for x in ("", "_lo", "_hi"):
            R[f"k_{pre}{x}"] = pick("k{}" + x)
        R[f"rate_{pre}"] = R[f"k_{pre}"] * R[f"net_count_{pre}"]
        R[f"rate_{pre}_lo"] = R[f"k_{pre}_lo"] * R[f"net_count_{pre}"]
        R[f"rate_{pre}_hi"] = R[f"k_{pre}_hi"] * R[f"net_count_{pre}"]
        R[f"rate_{pre}_scr"] = R[f"k_{pre}"] * R[f"net_count_{pre}_scr"]
    R.to_parquet(f"{O2}/rates_daily_v2.parquet", index=False)
    sc = [c for c in R.columns if c.startswith(("n0", "net0", "rate", "n_v2", "net_count_v2", "sumP_", "net_sumP"))]
    M = R.assign(month=R.date.str[:7]).groupby(["unit_kind", "unit", "label", "model", "region", "month", "day_class"])[sc].sum(
        min_count=1).reset_index()
    M.to_csv(f"{O2}/rates_monthly_v2.csv", index=False)
    log(f"[rates] {R.unit.nunique()} units, {len(R)} unit-label-days")
    V = pd.concat([view_rows(R, n) for n, _, _ in VIEWS], ignore_index=True)
    summary(V)
    step_checks_v2(V[V.view.isin(["headline", "wave_classrules"])])


def summary(V):
    """view x region x era x day class: per unit-day mean of the net count (data minus time-reversed control) and of the
    k-scaled rate (NaN where k is unavailable), 90 % day bootstrap (1000 draws)."""
    rs = np.random.default_rng(11)
    rows = []
    pref = {n: p for n, _, p in VIEWS}

    def boot(s, n):
        if not len(s):
            return np.nan, np.nan, np.nan
        bs = []
        for _ in range(1000):
            ii = rs.integers(0, len(s), len(s))
            bs.append(s[ii].sum() / n[ii].sum())
        return s.sum() / n.sum(), float(np.percentile(bs, 5)), float(np.percentile(bs, 95))
    for (view, reg, era, cls), g in V.groupby(["view", "region", "era", "day_class"]):
        pre = pref[view]
        for rcol, tag, raw in [(f"rate_{pre}", "rule", (f"n_{pre}_data", f"n_{pre}_rev", f"n_{pre}_c1")),
                               (f"rate_{pre}_scr", "rule_scrambled_subtracted", (f"n_{pre}_data", f"n_{pre}_c1", f"n_{pre}_rev"))] + \
                              [(f"rate{t}", f"cut_{c:.2f}", (f"n{t}_data", f"n{t}_rev", f"n{t}_c1")) for t, c in RCUTS.items()]:
            gg = g[g[raw[0]].notna()]
            if not len(gg):
                continue
            netc = (gg[raw[0]] - gg[raw[1]]).to_numpy()
            d = pd.DataFrame(dict(date=gg.date.to_numpy(), net=netc, rate=gg[rcol].to_numpy(), one=1.0))
            dd = d.groupby("date").agg(net=("net", "sum"), rate=("rate", lambda x: x.sum(min_count=1)), n=("one", "sum"),
                                        k_ok=("rate", lambda x: x.notna().all()))
            n = dd.n.to_numpy()
            ne, nlo, nhi = boot(dd.net.to_numpy(), n)
            kok = bool(dd.k_ok.all())
            re_, rlo, rhi = boot(dd.rate.to_numpy(), n) if kok else (np.nan, np.nan, np.nan)
            rows.append(dict(view=view, region=reg, era=era, day_class=cls, column=tag, units=gg.unit.nunique(),
                             unit_days=int(n.sum()), days=len(dd), data=float(gg[raw[0]].sum()),
                             rev_control=float(gg[f"{raw[0][:-5]}_rev"].sum()), scrambled_control=float(gg[f"{raw[0][:-5]}_c1"].sum()),
                             net_per_unit_day=ne, net_boot90_lo=nlo, net_boot90_hi=nhi, k_available=kok,
                             rate_per_unit_day=re_, rate_boot90_lo=rlo, rate_boot90_hi=rhi,
                             subtracted="scrambled control" if tag == "rule_scrambled_subtracted" else "time-reversed control",
                             status=("headline rule" if view == "headline" else "class rule") if tag.startswith("rule") else "column only"))
        x = g[["date", "sumP_scaled"]].dropna()
        if len(x):
            dd = x.groupby("date").sumP_scaled.agg(["sum", "size"])
            e, lo, hi = boot(dd["sum"].to_numpy(), dd["size"].to_numpy())
            rows.append(dict(view=view, region=reg, era=era, day_class=cls, column="sumP", units=g.unit.nunique(),
                             unit_days=int(dd["size"].sum()), days=len(dd), rate_per_unit_day=e, rate_boot90_lo=lo,
                             rate_boot90_hi=hi, k_available=True, subtracted="time-reversed control", status="column only"))
    S = pd.DataFrame(rows)
    S.to_csv(f"{O2}/rates_summary_v2.csv", index=False)
    log(f"[rates] summary {len(S)} rows")


def step_checks_v2(V):
    C = V1.core_changes()
    C.to_csv(f"{O2}/core_set_changes_v2.csv", index=False)
    rows = []
    pref = {n: p for n, _, p in VIEWS}
    for (view, reg_kind), Rg in V.groupby(["view", "region"]):
        pre = pref[view]
        qcol = f"net_count_{pre}" if view == "headline" else "net085"
        Rg = Rg.assign(q=(Rg.day_class == "quiet"),
                       net_q=np.where(Rg.day_class == "quiet", Rg[qcol].fillna(0), 0.0),
                       net_t=np.where(Rg.day_class == "tremor", Rg[f"net_count_{pre}"].fillna(0), 0.0))
        tot = Rg.groupby("date").agg(net_q=("net_q", "sum"), net_t=("net_t", "sum"), quiet=("q", "mean"),
                                     n_trem=("n_trem60", "sum")).reset_index()
        tot["dt"] = pd.to_datetime(tot.date)
        tot = tot.set_index("dt")
        for r in C.itertuples():
            d = pd.Timestamp(r.date)
            for wdays in (21, 45):
                def stat(center):
                    b = tot[(tot.index >= center - pd.Timedelta(days=wdays)) & (tot.index < center)]
                    a_ = tot[(tot.index > center) & (tot.index <= center + pd.Timedelta(days=wdays))]
                    q = lambda x: x[x.quiet >= 0.5].net_q.mean() if (x.quiet >= 0.5).sum() >= 3 else np.nan
                    th = lambda x: x.net_t.sum() / (x.n_trem.sum() * 2.5 / 60) if x.n_trem.sum() > 24 else np.nan
                    return q(a_) - q(b), th(a_) - th(b)
                dq, dt_ = stat(d)
                ref_q, ref_t = [], []
                for y in range(2007, 2027):
                    if y == d.year:
                        continue
                    try:
                        c2 = d.replace(year=y)
                    except ValueError:
                        continue
                    x, z = stat(c2)
                    ref_q.append(x); ref_t.append(z)
                rq = np.array(ref_q, float); rt = np.array(ref_t, float)
                zq = (dq - np.nanmean(rq)) / np.nanstd(rq) if np.isfinite(rq).sum() >= 4 and np.nanstd(rq) > 0 else np.nan
                zt = (dt_ - np.nanmean(rt)) / np.nanstd(rt) if np.isfinite(rt).sum() >= 4 and np.nanstd(rt) > 0 else np.nan
                rows.append(dict(view=view, region=reg_kind, date=r.date, netsta=r.netsta, change=r.change, window_days=wdays,
                                 quiet_count_column=qcol, d_quiet_rate=dq, z_quiet=zq, d_rate_per_tremor_hour=dt_, z_tremor_hour=zt,
                                 flag=bool((np.isfinite(zq) and abs(zq) > 2) or (np.isfinite(zt) and abs(zt) > 2))))
    S = pd.DataFrame(rows)
    S.to_csv(f"{O2}/rates_step_checks_v2.csv", index=False)
    log(f"[steps] {len(C)} changes; flagged {int(S.flag.sum()) if len(S) else 0} of {len(S)} checks")


def tables(a):
    """summary tables for POST_V2_RESULTS.md -> post_v2/summary_*.csv"""
    cat = pq.read_table(f"{O2}/catalogue_v2.parquet", columns=[
        "date", "family", "family_uid", "region", "P_own", "P_own_pnet", "P_nbr3", "day_class", "unknown_kind",
        "selected_v2", "selected_v2_wave_classrules", "selected_v2_pnet_classrules"]).to_pandas()
    cat["year"] = cat.date.str[:4].astype(int)
    rows = []
    TV = (("headline", "selected_v2", "P_own_pnet"), ("wave_classrules", "selected_v2_wave_classrules", "P_own"),
          ("pnet_classrules", "selected_v2_pnet_classrules", "P_own_pnet"))
    if HEAD_LAB == "own":
        TV = (("headline", "selected_v2", "P_own"), ("wave_classrules", "selected_v2_wave_classrules", "P_own"))
    for view, selc, pcol in TV:
        pl = np.where(cat.region == "OlyA", cat.P_nbr3, cat[pcol])
        sel = cat[selc].to_numpy()
        for (reg, cls), ii in cat.groupby(["region", "day_class"]).indices.items():
            for per, m in (("all", np.ones(len(ii), bool)), ("2012-2026", cat.year.to_numpy()[ii] >= 2012)):
                jj = ii[m]
                d = dict(view=view, region=reg, day_class=cls, period=per, candidates=len(jj), days=cat.date.iloc[jj].nunique(),
                         kept=int(sel[jj].sum()), families_kept=cat.family.iloc[jj][sel[jj]].nunique())
                for c in (0.5, 0.75, 0.85, 0.9):
                    d[f"P>={c}"] = int((pl[jj] >= c).sum())
                rows.append(d)
    S1 = pd.DataFrame(rows)
    S1.to_csv(f"{O2}/summary_catalogue_region_class.csv", index=False)
    # headline: per family repeats 2012-2026 by day class
    h = cat[cat.selected_v2 & (cat.year >= 2012)]
    fam = h.pivot_table(index=["region", "family", "family_uid"], columns="day_class", values="date", aggfunc="size", fill_value=0)
    fam["total"] = fam.sum(axis=1)
    fam = fam.reset_index().sort_values(["region", "total"], ascending=[True, False])
    fam.to_csv(f"{O2}/summary_headline_families_2012_2026.csv", index=False)
    # headline false share 2012-2026 from the rate table (Olympic A: merged group counts)
    R = pd.read_parquet(f"{O2}/rates_daily_v2.parquet", columns=["unit_kind", "unit", "label", "date", "region", "day_class",
                                                                 "n_v2_data", "n_v2_rev", "n_v2_c1", "cut_v2"])
    R = R[R.cut_v2.notna() & (R.date >= "2012-01-01")]
    fs = R.groupby(["region", "day_class"]).agg(units=("unit", "nunique"), data=("n_v2_data", "sum"), rev_control=("n_v2_rev", "sum"),
                                                 scrambled_control=("n_v2_c1", "sum"), cut=("cut_v2", "first")).reset_index()
    fs["false_share_rev"] = fs.rev_control / fs.data
    fs["false_share_scrambled"] = fs.scrambled_control / fs.data
    fs.to_csv(f"{O2}/summary_headline_false_share_2012_2026.csv", index=False)
    # class-rule NW quiet (WAVE, 0.85) kept, for continuity with the quiet-day waveform test
    nq = cat[cat.selected_v2_wave_classrules & (cat.region == "NW") & (cat.day_class == "quiet")]
    nq.groupby("year").agg(repeats=("family", "size"), families=("family", "nunique")).reset_index().to_csv(
        f"{O2}/summary_nw_quiet_wave085_by_year.csv", index=False)
    for reg in ("NW", "SW", "ArmB", "OlyA"):
        f = fam[fam.region == reg]
        log(f"[tables] headline 2012-2026 {reg}: {int(f.total.sum())} repeats, {len(f)} families, {int((f.total >= 10).sum())} with >= 10")
    print(S1[S1.view == "headline"].to_string())
    print(fs.to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["scan", "rates", "all", "tables"])
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    os.makedirs(O2, exist_ok=True)
    if a.what in ("scan", "all"):
        scan(a)
    if a.what in ("rates", "all"):
        rates(a)
    if a.what in ("tables", "all"):
        tables(a)
