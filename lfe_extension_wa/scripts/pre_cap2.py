#!/usr/bin/env python
"""NW dominance-cap repeat on V2008 (test half) + A2011 (2026-10-03).  Rules: results/wa_pre/cap2/PREREGISTERED.md.

Helper stages for scripts/pre_cap2_queue.sh (read-only on everything outside results/wa_pre/cap2/; never reads
traces_dl*/traces_sealed):
  guard                      core station-day cache completeness on the V2008 + A2011 days (fails > 3 % on the 20)
  postguard <mfdir>          every cached core station-day loaded or dropped by the frozen > 5 %-guarded rule
  hfcheck <tag>              HF station-day jobs that returned ratios (fails > 3 % empty)
  equal <dirA> <dirB> <win>  matched-filter outputs identical on the days of window <win> (reproducibility)
  report                     cap2_comparison.csv, cap2_verdict.json, reversed-control columns
"""
import os, sys, glob, json, re
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pre_common as P
W = P.W

OD = f"{P.PRE}/cap2"
GOS = f"{P.STUDY}/results/wa_gonogo/setup"
EPIS = {"V2008": "v08", "A2011": "va"}           # episode -> run suffix (A2011 runs also hold V2010)
VARS = ["frozen", "cap30", "lodo"]
NW = "Royer_NW_WA"


def guard():
    core = P.core_list()
    usable = set(pd.read_csv(f"{GOS}/stations.csv").netsta)
    rows = []
    for w in ("V2008", "A2011"):
        for d, k in W.window_days(w):
            for ns in core:
                c = W.cache_path(ns, d)
                ok, err = False, ""
                if os.path.exists(c):
                    try:
                        z = np.load(c); _ = z["x"].shape; ok = True
                    except Exception as e:
                        err = repr(e)[:80]
                rows.append(dict(window=w, date=d, kind=k, netsta=ns, usable=ns in usable,
                                 raw_on_disk=bool(P.raw_files(ns, d)), cache=os.path.exists(c), cache_ok=ok, err=err))
    X = pd.DataFrame(rows)
    X.to_csv(f"{OD}/guard_station_days.csv", index=False)
    out = {}
    for lab, sub in [("all21", X), ("usable20", X[X.usable])]:
        for w in ("V2008", "A2011", "both"):
            s = sub if w == "both" else sub[sub.window == w]
            ond = s[s.raw_on_disk | s.cache]
            bad = ond[~ond.cache_ok]
            out[f"{lab}_{w}"] = dict(on_disk=len(ond), missing_or_error=len(bad), frac=len(bad) / max(len(ond), 1),
                                     by_station=bad.netsta.value_counts().to_dict())
    json.dump(out, open(f"{OD}/guard_summary.json", "w"), indent=1)
    print(json.dumps(out, indent=1))
    worst = max(out[f"usable20_{w}"]["frac"] for w in ("V2008", "A2011"))
    if worst > 0.03:
        P.marker("cap2/GUARD_FAILED", f"usable-20 missing/error fraction {worst:.2%}")
        sys.exit(1)


def postguard(mfdir):
    S = pd.read_csv(f"{GOS}/stations.csv")
    G = pd.read_csv(f"{OD}/guard_station_days.csv")
    G = G[G.netsta.isin(S.netsta) & G.cache_ok]
    rows = []
    for p in sorted(glob.glob(f"{mfdir}/run_*.json")):
        r = json.load(open(p))
        for ns in S.netsta:
            rows.append(dict(date=r["date"], netsta=ns, loaded=ns in r["nodes"], dead=ns in r.get("dead", {})))
    L = pd.DataFrame(rows).merge(G[["date", "netsta"]], on=["date", "netsta"], how="inner")
    other = int((~L.loaded & ~L.dead).sum())
    msg = (f"[postguard {os.path.basename(mfdir)}] cached core station-days {len(L)}; loaded {int(L.loaded.sum())}; "
           f"dropped (> 5 % guarded) {int(L.dead.sum())}; not loaded for any other reason {other}")
    print(msg)
    L.to_csv(f"{mfdir}/postguard.csv", index=False)
    if other / max(len(L), 1) > 0.03 or len(L) == 0:
        P.marker("cap2/GUARD_FAILED", msg)
        sys.exit(1)


def hfcheck(tag):
    log = open(f"{P.STUDY}/logs/hf.log").read().splitlines()
    jl = [l for l in log if f"[hf:{tag}]" in l and "station-day jobs" in l]
    nj = int(re.search(r"(\d+) station-day jobs", jl[-1]).group(1))
    H = pd.read_parquet(f"{OD}/mf_{tag}/hf.parquet")
    H["date"] = pd.to_datetime(np.floor(H.t / 86400) * 86400, unit="s").dt.strftime("%Y-%m-%d")
    got = H.groupby(["hf_sta", "date", "hf_node"]).ngroups
    msg = f"[hfcheck {tag}] {got} of {nj} station-day jobs returned ratios"
    print(msg)
    if (nj - got) / max(nj, 1) > 0.03:
        P.marker("cap2/GUARD_FAILED", msg)
        sys.exit(1)


def equal(da, db, win):
    days = [d for d, _ in W.window_days(win)]
    nd, maxd, nrow = 0, 0.0, 0
    for d in days:
        pa, pb = f"{da}/det_{d}.parquet", f"{db}/det_{d}.parquet"
        if not (os.path.exists(pa) and os.path.exists(pb)):
            print("missing", d, os.path.exists(pa), os.path.exists(pb)); nd += 1; continue
        A = pd.read_parquet(pa).sort_values(["series", "fid", "t"]).reset_index(drop=True)
        B = pd.read_parquet(pb).sort_values(["series", "fid", "t"]).reset_index(drop=True)
        if len(A) != len(B) or not (A[["series", "fid"]].equals(B[["series", "fid"]]) and np.allclose(A.t, B.t)):
            print("differ", d, len(A), len(B)); nd += 1; continue
        maxd = max(maxd, float(np.abs(A.net_mad - B.net_mad).max()) if len(A) else 0.0); nrow += len(A)
    msg = f"[equal {os.path.basename(da)} vs {os.path.basename(db)} {win}] {len(days)} days, {nd} differ, {nrow} candidates, max |d net_mad| {maxd:.2e}"
    print(msg)
    if nd or maxd > 1e-4:
        P.marker("cap2/REPRO_FAILED", msg)
        sys.exit(1)


# ---------------------------------------------------------------- report
def load_det(d, win):
    ps = [f"{d}/det_{x}.parquet" for x, _ in W.window_days(win)]
    return pd.concat([pd.read_parquet(p) for p in ps if os.path.exists(p)], ignore_index=True)


def load_fd(d, win):
    ps = [f"{d}/famday_{x}.parquet" for x, _ in W.window_days(win)]
    return pd.concat([pd.read_parquet(p) for p in ps if os.path.exists(p)], ignore_index=True)


def extra_cols(fwd, rev, win, F, T, scan, tdays):
    """controls 1+2 and the time-reversed control on the scorer's NW family-windows (same day sets)."""
    Fn = F[(F.region == NW) & F.val_windows.str.contains(win)]
    FD = load_fd(fwd, win); D = load_det(fwd, win); R = load_det(rev, win)
    D = D[D.net_mad >= W.T_CORE]; R = R[(R.series == "data") & (R.net_mad >= W.T_CORE)]
    scored = (lambda t: np.ones(len(t), bool)) if win == "A2011" else (lambda t: W.half_of(t) == 1)
    Dg = {k: g for k, g in D.groupby(["series", "fid"])}
    Rg = dict(tuple(R.groupby("fid")))
    Tg = dict(tuple(T[T.window == win].groupby("fid")))
    FDg = dict(tuple(FD.groupby("fid")))
    acc = dict(famdays=0, d=0, c1=0, c2=0, r=0, r_sc=0, r_match=0, sep_rev=0, n_lab=0)
    for r in Fn.itertuples():
        fd = FDg.get(r.fid)
        if fd is None:
            continue
        sd = set(scan[scan.det_source == r.det_source].date)
        days = set(fd[fd.observable & fd.date.isin(sd) & (fd.kind == "episode")].date)
        fo = fd[fd.observable]
        lab = set(fo[fo.kind == "episode"].date) if win == "V2008" else set(fo.date) & tdays
        if not days:
            continue
        acc["famdays"] += len(days)
        get = lambda s: Dg.get((s, r.fid), D.iloc[:0])
        acc["d"] += int(get("data").date.isin(days).sum())
        acc["c1"] += int(get("control").date.isin(days).sum())
        acc["c2"] += int(get("control2").date.isin(days).sum())
        rr = Rg.get(r.fid, R.iloc[:0])
        acc["r"] += int(rr.date.isin(days).sum())
        rs = rr[rr.date.isin(days)]
        rs = rs[scored(rs.t.to_numpy())] if len(rs) else rs
        tt = Tg.get(r.fid)
        for ds, g in rs.groupby("date"):
            own = np.sort(tt[(tt.date == ds)].t.to_numpy()) if tt is not None else np.zeros(0)
            own = own[scored(own)] if len(own) else own
            ia, _ = W.match_one_to_one(np.sort(g.t.to_numpy()), own)
            acc["r_match"] += len(ia)
        acc["r_sc"] += len(rs)
        if lab:
            acc["n_lab"] += 1
            ld = int(get("data").date.isin(lab).sum()); lr = int(rr.date.isin(lab).sum())
            acc["sep_rev"] += int((ld + 1) / (lr + 1) >= 2)
    fdn = max(acc["famdays"], 1)
    return dict(x_famdays=acc["famdays"], x_data_per_famday=acc["d"] / fdn, ctl1_per_famday=acc["c1"] / fdn,
                ctl2_per_famday=acc["c2"] / fdn, rev_per_famday=acc["r"] / fdn,
                rev_n_scored=acc["r_sc"], rev_prec_own=acc["r_match"] / acc["r_sc"] if acc["r_sc"] else np.nan,
                famwin_data_ge2x_rev=acc["sep_rev"], n_famwin_label=acc["n_lab"])


def report():
    F = pd.read_csv(f"{GOS}/families.csv")
    T = pd.read_parquet(f"{GOS}/truth.parquet")
    scan = pd.read_csv(f"{GOS}/scan_days.csv")
    t = pd.read_csv(f"{W.EXT}/NW_WA_tremor_days.csv")
    tdays = set(pd.to_datetime(t[t.columns[0]]).dt.strftime("%Y-%m-%d")[t[t.columns[1]] > 0])
    rows = []
    for win, suf in list(EPIS.items()) + [("V2010", "va")]:
        for v in VARS:
            sd = f"{OD}/score_{v}_{suf}/summary.csv"
            S = pd.read_csv(sd)
            s = S[(S.window == win) & (S.region == NW)].iloc[0]
            r = dict(episode=win, run={"frozen": "uncapped", "cap30": "capped_0.3", "lodo": "dominant_site_out"}[v],
                     n_fam_observable=int(s.n_fam_observable), n_truth_test=int(s.n_truth_test),
                     famdays_obs=int(s.famdays_obs), n_ours_cat_core=int(s.n_ours_cat_core),
                     prec_own=s.prec_own_cat_core, prec_nbr3=s.prec_nbr3_cat_core, recall_own=s.recall_own_cat_core,
                     prec_own_ctl=s.prec_own_ctl_core, data_per_famday=s.data_per_famday_core,
                     ctl_per_famday=s.ctl_per_famday_core, sepboth2_famwin=int(s.sepboth2_core_data),
                     sep2_ctl_vs_ctl=int(s.sep2_core_ctl), sep2_vs_ctl1=int(s.sep2_core_data))
            if win != "V2010":
                r.update(extra_cols(f"{OD}/mf_{v}_{suf}", f"{OD}/mf_rev{v}_{suf}", win, F, T, scan, tdays))
            fds = load_fd(f"{OD}/mf_{v}_{suf}", win).merge(F[["fid", "region"]], on="fid")
            g = fds[fds.observable & (fds.region == NW)]
            r["share_max_frozen_median"] = float(g.share_max_frozen.median()) if "share_max_frozen" in g else np.nan
            r["share_max_used_median"] = float(g.share_max_used.median()) if "share_max_used" in g and g.share_max_used.notna().any() else np.nan
            r["cap_infeasible_frac"] = float(g.cap_infeasible.fillna(False).astype(bool).mean()) if v == "cap30" else np.nan
            rows.append(r)
    C = pd.DataFrame(rows)
    C.to_csv(f"{OD}/cap2_comparison.csv", index=False)
    print(C.to_string())
    # reproducibility of V2010 vs results/wa_pre/cap
    old = pd.read_csv(f"{P.PRE}/cap/cap_comparison.csv")
    old = old[old.region == NW].set_index("run")
    rep = {}
    for v, lab in [("frozen", "uncapped"), ("cap30", "capped_0.3"), ("lodo", "dominant_site_out")]:
        n = C[(C.episode == "V2010") & (C.run == lab)].iloc[0]
        rep[v] = dict(prec_own_new=float(n.prec_own), prec_own_old=float(old.loc[v, "prec_own"]),
                      famwin_new=int(n.sepboth2_famwin), famwin_old=int(old.loc[v, "sepboth2_famwin"]),
                      recall_new=float(n.recall_own), recall_old=float(old.loc[v, "recall_own"]))
        rep[v]["same"] = bool(abs(rep[v]["prec_own_new"] - rep[v]["prec_own_old"]) < 1e-9 and
                              rep[v]["famwin_new"] == rep[v]["famwin_old"] and
                              abs(rep[v]["recall_new"] - rep[v]["recall_old"]) < 1e-9)
    # rule
    ver = dict(rule="adopt for NW iff on BOTH episodes: dprec >= 0.02 AND famwin_cap >= famwin_uncapped AND "
                    "recall_cap >= 0.8 recall_uncapped", v2010_reproduction=rep)
    allpass = True
    for win in EPIS:
        u = C[(C.episode == win) & (C.run == "uncapped")].iloc[0]
        c = C[(C.episode == win) & (C.run == "capped_0.3")].iloc[0]
        e = dict(prec_uncapped=float(u.prec_own), prec_capped=float(c.prec_own), dprec=float(c.prec_own - u.prec_own),
                 famwin_uncapped=int(u.sepboth2_famwin), famwin_capped=int(c.sepboth2_famwin),
                 recall_uncapped=float(u.recall_own), recall_capped=float(c.recall_own),
                 recall_rel_change=float(c.recall_own / u.recall_own - 1))
        e["pass_prec"] = e["dprec"] >= 0.02
        e["pass_famwin"] = e["famwin_capped"] >= e["famwin_uncapped"]
        e["pass_recall"] = e["recall_capped"] >= 0.8 * e["recall_uncapped"]
        e["pass"] = bool(e["pass_prec"] and e["pass_famwin"] and e["pass_recall"])
        allpass &= e["pass"]
        ver[win] = e
    ver["decision"] = "ADOPT cap 0.3 for NW production" if allpass else "NW stays uncapped"
    json.dump(ver, open(f"{OD}/cap2_verdict.json", "w"), indent=1, default=lambda x: bool(x) if isinstance(x, np.bool_) else x)
    print(json.dumps(ver, indent=1, default=lambda x: bool(x) if isinstance(x, np.bool_) else x))


if __name__ == "__main__":
    cmd = sys.argv[1]
    {"guard": guard, "postguard": postguard, "hfcheck": hfcheck, "equal": equal, "report": report}[cmd](*sys.argv[2:])
