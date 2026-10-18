#!/usr/bin/env python
"""WA PRODUCTION setup (W0 freeze, W1 completeness inputs, roster, day roles, truth, W4a template map, fixed rate sets).

  freeze     sha256 of every input the run reads (scripts, templates, setups, W2 tables, catalogue)  -> w0/frozen_inputs.csv
             `verify` re-hashes them at the end of the queue                                         -> w0/verify_<utc>.csv
  roster     stations: the 94 trial/step-4b stations + downloaded tiers A, B (traces_dl) + C, D, D2, T (traces_dl2);
             nets 7D, X6, YR, Z5 and XD nodes excluded (user 2026-10-03); coordinates from the bridge-node tables, the
             step-4b setup and the temporary-network inventory                                      -> setup/stations.csv
  families   the 329 v1.1 representatives (fids of every earlier WA run) + region, position, depth, v1.2 uid / dup group,
             the <= 3 km neighbour lists and the Olympic Arm A <= 3-km POSITION GROUPS (single linkage)  -> setup/families.csv
  days       one row per day 2007-01-01 .. 2026-09-15: era, month role (fit / sealed_headline 2007-01 / sealed_semi 2009-05
             / excluded 2011-04 / prod), truth-scan flag                                              -> setup/day_table.csv
  truth      every family's own array detections (det_source grp == family_id) in the truth months; scanned days per
             source run; the SEALED months go to setup/sealed/ (read once, by wa_prod_fit.py open_sealed)
  rawinv     per (station, day) the first raw root holding a file (directory listings; dl2 only when WA_DL2_DONE)
  templates  W4a: (1) D1 rebuild 'v10b' = the frozen wa_templates recipe on the V2010 BUILD-half array detections at the
             nodes whose only templates were tpl_step4b_a (built from A2011 events), so A2011 stays leak-free;
             (2) library index (tags build, v10b, t2012new, step4b_t2012/13/14; first tag wins; every >= 30-event
             template); (3) carry-over map = W2 swap checks 'above control (carry-over OK)' where the new node has no
             template of its own; (4) pending re-templates (W2 'AT CONTROL LEVEL' after a change: PGC 2015-01-12, VGZ
             2017-03-03, DOSE 2019-05-07, HDW 2020-12-04 -> W4b) ; B014 2023 check                     -> setup/tpl_*.csv
  core       fixed rate sets per era by a rule fixed here (before any run): permanent station, raw on >= 80 % of the era's
             days, a W4a template source (own node with >= 30 families, or approved carry-over) on >= 95 % of its raw days
             in the era, not XU.W060, LRIV out after 2019-05-30 (plan D10 / SVI policy)            -> setup/core_sets.json
  excl       W2 exclusions applied to the matched filter: only DEAD or broken station-months / days (W2 status dead /
             inherit_dead; day-level 'dead' windows) -- NOT 'control_level' months, NOT W2 'flat' (repeated-sample
             fraction >= 1 %, which flags years of normal PB data); flat / dead channels are also caught per window by the
             frozen energy-floor guard and the drop-dead (> 5 % guarded) rule      -> setup/exclusions.csv
"""
from __future__ import annotations
import os, sys, json, glob, argparse, time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P
W = P.W

SU = f"{P.OUT}/setup"
T1S = f"{P.RES}/merlin_tests/t1_setup"


def nodes_table():
    N = pd.concat([pd.read_csv(f"{W.EXT}/{r}_bridge_nodes_strict.csv") for r in ["NW_WA", "SW_WA", "OLY2010"]]
                  + [pd.read_csv(f"{P.RES}/step4b/setup/extra_nodes.csv")]).drop_duplicates("node")
    N["netsta"] = N.net + "." + N.sta
    return N


# ------------------------------------------------------------------------------------------------ W0
def freeze_list():
    fs = sorted(glob.glob(f"{P.STUDY}/scripts/wa_prod_*.py")) + [f"{P.STUDY}/scripts/{x}" for x in (
        "wa_common.py", "wa_cache.py", "wa_mf.py", "wa_hf.py", "wa_templates.py", "step1_features.py", "pnet_common.py",
        "wa_prod_queue.sh")]
    fs += sorted(glob.glob(f"{P.RES}/templates/tpl_*.npz")) + sorted(glob.glob(f"{P.RES}/templates/tpl_qc_*.csv"))
    fs += [f"{T1S}/{x}" for x in ("families.csv", "stations.csv", "truth.parquet", "scan_days.csv")]
    fs += [f"{P.RES}/step4b/setup/{x}" for x in ("stations.csv", "extra_nodes.csv", "families.csv")]
    fs += [f"{P.RES}/wa_pre/w2/{x}" for x in ("r2_swap_checks.csv", "r2_health_monthly.csv", "r2_bad_for_dvv.csv",
                                              "r2_exclusions.csv", "r2_effective_sets.csv")]
    fs += [f"{P.RES}/wa_pre/w1/{x}" for x in ("day_roles.csv", "family_month_roles.csv", "family_month_roles_olympic.csv")]
    fs += [f"{P.RES}/wa_gonogo/truth_array_frac.csv", f"{P.ROOT}/data/catalogue/families_v1_1.csv",
           f"{P.ROOT}/data/catalogue/families_v1_2.csv", f"{P.ROOT}/catalogs/pnsn_tremor_cascadia_full.csv",
           f"{W.EXT}/tremor_episodes.csv", f"{P.RES}/WA_PROD_PLAN.md", f"{P.RES}/pnet/models/ft1.pt"]
    return [f for f in fs if os.path.exists(f)]


def freeze(a):
    P.ensure(f"{P.OUT}/w0")
    p = f"{P.OUT}/w0/frozen_inputs.csv"
    if os.path.exists(p):
        P.log(f"[freeze] {p} exists; not rewritten"); return
    rows = [dict(path=f, bytes=os.path.getsize(f), mtime=os.path.getmtime(f), sha256=P.sha256(f)) for f in freeze_list()]
    pd.DataFrame(rows).to_csv(p, index=False)
    P.log(f"[freeze] {len(rows)} files, {sum(r['bytes'] for r in rows) / 1e9:.2f} GB hashed")


def verify(a):
    F = pd.read_csv(f"{P.OUT}/w0/frozen_inputs.csv")
    rows = []
    for r in F.itertuples():
        ok = os.path.exists(r.path)
        h = P.sha256(r.path) if ok else ""
        rows.append(dict(path=r.path, exists=ok, same=bool(ok and h == r.sha256)))
    V = pd.DataFrame(rows)
    V.to_csv(f"{P.OUT}/w0/verify_{time.strftime('%Y%m%dT%H%M', time.gmtime())}.csv", index=False)
    P.log(f"[verify] {len(V)} frozen inputs, changed {int((~V.same).sum())}: {V[~V.same].path.tolist()[:10]}")


# ------------------------------------------------------------------------------------------------ roster
def roster(a):
    P.ensure(SU)
    N = nodes_table()
    co = N.groupby("netsta")[["lat", "lon"]].first()
    S4 = pd.read_csv(f"{P.RES}/step4b/setup/stations.csv")
    rows = {}
    for r in S4.itertuples():
        rows[r.netsta] = dict(netsta=r.netsta, lat=r.lat, lon=r.lon, tier="trial" if r.netsta in set(
            pd.read_csv(f"{T1S}/stations.csv").netsta) else "step4b")
    dl = pd.read_csv(f"{P.RES}/wa_dl/completeness_by_station.csv")
    for r in dl.itertuples():
        rows.setdefault(r.netsta, dict(netsta=r.netsta, lat=np.nan, lon=np.nan))["tier"] = r.tier
    d2 = pd.read_csv(f"{P.RES}/wa_dl2/plan_summary.csv")
    plan2 = set(pd.read_csv(f"{P.RES}/wa_dl2/plan.csv", usecols=["netsta"]).netsta)
    for r in d2.itertuples():
        if r.netsta not in plan2:                       # removed from the plan (OBS / nodes, 2026-10-02)
            continue
        rows.setdefault(r.netsta, dict(netsta=r.netsta, lat=np.nan, lon=np.nan))["tier"] = r.tier
    inv = pd.read_csv(f"{P.RES}/wa_temp_networks/inventory.csv")
    inv["netsta"] = inv.net + "." + inv.sta.astype(str)
    td = pd.read_csv(f"{P.RES}/wa_dl2/tier_d_stations.csv")
    if "netsta" not in td:
        td["netsta"] = td.net + "." + td.sta.astype(str)
    ic = pd.concat([inv[["netsta", "lat", "lon"]], td[["netsta", "lat", "lon"]]]).groupby("netsta")[["lat", "lon"]].first()
    out = []
    for ns, r in rows.items():
        net = ns.split(".")[0]
        if net in P.BANNED_NETS:
            continue
        if net == "XD" and ns in set(inv[inv.all_bands.astype(str).str.startswith("DP")].netsta):
            continue
        la, lo = r.get("lat", np.nan), r.get("lon", np.nan)
        if not np.isfinite(la):
            if ns in co.index:
                la, lo = co.loc[ns, "lat"], co.loc[ns, "lon"]
            elif ns in ic.index:
                la, lo = ic.loc[ns, "lat"], ic.loc[ns, "lon"]
        out.append(dict(netsta=ns, net=net, sta=ns.split(".")[1], lat=la, lon=lo, tier=r.get("tier", ""),
                        temp=net in W.TEMP_NETS, no_template=ns in P.NO_TEMPLATE))
    S = pd.DataFrame(out).sort_values("netsta").reset_index(drop=True)
    miss = S[~np.isfinite(S.lat)]
    if len(miss):
        P.log(f"[roster] {len(miss)} stations without coordinates are dropped: {miss.netsta.tolist()[:20]}")
        S = S[np.isfinite(S.lat)].reset_index(drop=True)
    S.to_csv(f"{SU}/stations.csv", index=False)
    P.log(f"[roster] {len(S)} stations; tiers {S.tier.value_counts().to_dict()}; temp {int(S.temp.sum())}")


# ------------------------------------------------------------------------------------------------ families
def families(a):
    P.ensure(SU)
    F = pd.read_csv(f"{T1S}/families.csv")
    V1 = pd.read_csv(f"{P.ROOT}/data/catalogue/families_v1_1.csv").set_index("family_uid")
    V2 = pd.read_csv(f"{P.ROOT}/data/catalogue/families_v1_2.csv")
    v2 = V2.set_index("family_uid")
    F["reg"] = F.region.map(P.REGIONS)
    F["depth_km"] = F.family_uid.map(V1.depth_km)
    F["dup_group_v12"] = F.family_uid.map(v2.dup_group) if "dup_group" in v2 else np.nan
    F["region_v12"] = F.family_uid.map(v2.region)
    la, lo = F.lat_use.to_numpy(), F.lon_use.to_numpy()
    d = W.km(la[:, None], lo[:, None], la[None], lo[None])
    F["nbr3"] = [";".join(str(int(F.fid[j])) for j in np.nonzero(d[i] <= W.SITE_KM)[0]) for i in range(len(F))]
    # Olympic Arm A position groups (<= 3 km single linkage among Olympic_2010 families)
    o = np.flatnonzero((F.region == "Olympic_2010").to_numpy())
    lab = {i: i for i in o}
    for i in o:
        for j in o:
            if j > i and d[i, j] <= W.SITE_KM:
                a_, b_ = lab[i], lab[j]
                for k in lab:
                    if lab[k] == b_:
                        lab[k] = a_
    F["pos_group"] = ""
    for i, g in lab.items():
        F.loc[i, "pos_group"] = f"OLYA_g{g:03d}"
    F.to_csv(f"{SU}/families.csv", index=False)
    P.log(f"[families] {len(F)} {F.reg.value_counts().to_dict()}; Olympic A position groups "
          f"{F[F.pos_group != ''].pos_group.nunique()} (largest {F[F.pos_group != ''].pos_group.value_counts().max()})")


# ------------------------------------------------------------------------------------------------ days + truth
def truth(a):
    P.ensure(SU, f"{SU}/sealed")
    F = pd.read_csv(f"{SU}/families.csv")
    T, SC = [], []
    for src, g in F.groupby("det_source"):
        D = pd.read_parquet(f"{P.ROOT}/{src}", columns=["grp", "t", "date"])
        D["date"] = D.date.astype(str)
        SC.append(pd.DataFrame(dict(det_source=src, date=sorted(set(D.date)))))
        D = D[D.date.str[:7].isin(P.TRUTH_MONTHS)]
        m = dict(zip(g.family_id, g.fid))
        D = D[D.grp.isin(m)].copy()
        D["fid"] = D.grp.map(m).astype(int); D["det_source"] = src
        T.append(D[["fid", "t", "date", "det_source"]])
    T = pd.concat(T, ignore_index=True)
    T["month"] = T.date.str[:7]
    T["role"] = T.month.map(P.month_role)
    T["half"] = W.half_of(T.t.to_numpy())
    SCAN = pd.concat(SC, ignore_index=True)
    SCAN = SCAN[SCAN.date.str[:7].isin(P.TRUTH_MONTHS)]
    sealed = T.role.str.startswith("sealed")
    T[~sealed].to_parquet(f"{SU}/truth.parquet", index=False)
    T[sealed].to_parquet(f"{SU}/sealed/truth_sealed.parquet", index=False)
    open(f"{SU}/sealed/README.txt", "w").write("Sealed truth (2007-01 headline, 2009-05 semi-independent). Read ONCE by "
                                                 "wa_prod_fit.py open_sealed (SEALED_OPENED.json written first). Never fitted.\n")
    SCAN.to_csv(f"{SU}/scan_days.csv", index=False)
    # truth_array_frac per family-month (carried in every truth table)
    TA = pd.read_csv(f"{P.RES}/wa_gonogo/truth_array_frac.csv")
    TA.to_csv(f"{SU}/truth_array_frac.csv", index=False)
    P.log(f"[truth] open rows {int((~sealed).sum())} by month {T[~sealed].month.value_counts().sort_index().to_dict()}; "
          f"sealed rows {int(sealed.sum())} (not read here further); scan rows {len(SCAN)}")


def day_table(a):
    P.ensure(SU)
    SCAN = pd.read_csv(f"{SU}/scan_days.csv")
    sd = set(SCAN.date)
    rows = []
    for d in P.days():
        m = d[:7]
        rows.append(dict(date=d, month=m, era=P.era_of(d), role=P.month_role(m), truth_scanned=d in sd))
    D = pd.DataFrame(rows)
    p = f"{SU}/day_table.csv"
    if os.path.exists(p):
        P.log(f"[days] {p} exists (frozen); not rewritten"); return
    D.to_csv(p, index=False)
    P.log(f"[days] {len(D)} days; roles {D.role.value_counts().to_dict()}; truth-scanned {int(D.truth_scanned.sum())}")


# ------------------------------------------------------------------------------------------------ raw inventory
def _list_station(args):
    root, ns, years = args
    out = []
    for y in years:
        dd = f"{root}/{ns}/{y}"
        if not os.path.isdir(dd):
            continue
        try:
            with os.scandir(dd) as it:
                for e in it:
                    if not e.name.endswith(".mseed"):
                        continue
                    try:
                        if os.path.getsize(f"{dd}/{e.name}") <= 0:
                            continue
                    except OSError:
                        continue
                    doy = int(e.name[:3])
                    out.append((ns, str((pd.Timestamp(f"{y}-01-01") + pd.Timedelta(days=doy - 1)).date()), root))
        except OSError:
            continue
    return out


def rawinv(a):
    S = pd.read_csv(f"{SU}/stations.csv")
    roots = P.prod_roots()
    years = [str(y) for y in range(2007, 2027)]
    jobs = [(r, ns, years) for r in roots for ns in S.netsta if os.path.isdir(f"{r}/{ns}")]
    rows = []
    with ThreadPoolExecutor(16) as ex:
        for o in ex.map(_list_station, jobs):
            rows += o
    R = pd.DataFrame(rows, columns=["netsta", "date", "root"])
    R = R[(R.date >= P.T0) & (R.date <= P.T1)]
    R["rank"] = R.root.map({r: i for i, r in enumerate(roots)})
    R = R.sort_values("rank").drop_duplicates(["netsta", "date"]).drop(columns=["rank"])
    tag = "dl2" if P.dl2_ready() else "nodl2"
    R.to_parquet(f"{SU}/raw_inventory_{tag}.parquet", index=False)
    P.log(f"[rawinv] roots {len(roots)} (dl2 {'in' if P.dl2_ready() else 'NOT yet'}); {len(R)} station-days, "
          f"{R.netsta.nunique()} stations; by root {R.root.str.replace(P.ROOT + '/', '').value_counts().to_dict()}")


def raw_inventory():
    for tag in ("dl2", "nodl2"):
        p = f"{SU}/raw_inventory_{tag}.parquet"
        if os.path.exists(p):
            return pd.read_parquet(p), tag
    raise SystemExit("raw inventory missing")


# ------------------------------------------------------------------------------------------------ W4a templates
def build_v10b(a):
    """D1: rebuild the step-4b (A2011-event) templates from the V2010 BUILD-half array detections, frozen recipe."""
    od = f"{P.OUT}/templates"
    P.ensure(od)
    if os.path.exists(f"{od}/tpl_v10b.npz"):
        P.log("[v10b] exists"); return
    P.patch_wa()
    import wa_templates as WT
    Qa = pd.read_csv(f"{P.RES}/templates/tpl_qc_step4b_a.csv")
    Qb = pd.read_csv(f"{P.RES}/templates/tpl_qc_build.csv")
    only = set(Qa.node) - set(Qb.node)
    stas = sorted(set(Qa[Qa.node.isin(only)].netsta) - P.NO_TEMPLATE)
    T = pd.read_parquet(f"{T1S}/truth.parquet")
    EV = T[(T.window == "V2010") & (T.half == 0)][["fid", "t"]].copy()
    EV["date"] = pd.to_datetime(EV.t, unit="s").dt.strftime("%Y-%m-%d")
    P.log(f"[v10b] {len(EV)} V2010 build-half events on {EV.date.nunique()} days; {len(only)} nodes at {len(stas)} stations")
    from concurrent.futures import ProcessPoolExecutor
    acc = {}
    with ProcessPoolExecutor(a.workers) as ex:
        fut = [ex.submit(WT.stack_day, (d, g[["fid", "t"]], stas, only)) for d, g in EV.groupby("date")]
        for f in fut:
            d, out = f.result()
            for k, (s, n, live) in out.items():
                if k in acc:
                    s0, n0, l0 = acc[k]; acc[k] = (s0 + s, n0 + n, l0 | live)
                else:
                    acc[k] = (s.astype(np.float64), n, live)
    S4 = pd.read_csv(f"{P.RES}/step4b/setup/stations.csv")
    Q, TPL = WT.make_templates(acc, S4)
    for c_ in ["fid", "netsta", "node", "site", "n_events", "nch", "live", "clear", "imp", "dur", "why", "s_off", "offset_s", "snr_db"]:
        if c_ not in Q:
            Q[c_] = pd.Series(dtype=float if c_ not in ("netsta", "node", "site", "live", "why") else object)
    Q["clear"] = Q.clear.fillna(False).astype(bool)
    R = Qb[Qb.clear.fillna(False).astype(bool)]
    ref = float(np.median(10 ** (R.snr_db / 20.0)))
    Q["w"] = np.minimum(10 ** (Q.snr_db / 20.0), 3 * ref) / ref
    Q.to_csv(f"{od}/tpl_qc_v10b.csv", index=False)
    keys = list(TPL)
    np.savez_compressed(f"{od}/tpl_v10b.npz", keys=np.array([f"{f}|{ns}|{nd}" for f, ns, nd in keys]),
                        tpl=np.stack([TPL[k] for k in keys]) if keys else np.zeros((0, 3, W.NT), np.float32))
    n30 = Q[(Q.n_events >= W.MIN_EVENTS) & Q.offset_s.notna()]
    P.log(f"[v10b] {len(Q)} (family,node) stacks; >= 30 events {len(n30)} ({n30.fid.nunique()} families, {n30.node.nunique()} nodes); "
          f"step4b_a had {int((Qa.n_events >= 30).sum())}")


def tpl_dir(tag):
    return f"{P.OUT}/templates" if os.path.exists(f"{P.OUT}/templates/tpl_qc_{tag}.csv") else f"{P.RES}/templates"


def templates(a):
    """library index + carry-over map + pending re-template list."""
    build_v10b(a)
    rows = []
    for tag in P.TPL_TAGS_W4A:
        q = pd.read_csv(f"{tpl_dir(tag)}/tpl_qc_{tag}.csv", dtype={"live": str})
        q = q[(q.n_events >= W.MIN_EVENTS) & q.offset_s.notna() & q.w.notna()].copy()
        q["tag"] = tag
        rows.append(q)
    L = pd.concat(rows, ignore_index=True)
    L = L[~L.netsta.isin(P.NO_TEMPLATE) & ~L.netsta.str.split(".").str[0].isin(P.BANNED_NETS)]
    L = L.drop_duplicates(["fid", "netsta", "node"], keep="first")
    L.to_csv(f"{SU}/tpl_library.csv", index=False)
    nf = L.groupby("node").fid.nunique()
    # carry-over (W2 verdict) and pending re-templates
    N = nodes_table()
    SW = pd.read_csv(f"{P.RES}/wa_pre/w2/r2_swap_checks.csv", dtype=str)
    cm = []
    for r in SW.itertuples():
        if len(r.date) != 10:
            continue
        g = N[N.netsta == r.netsta].sort_values("start")
        new = g[g.start == r.date]
        if not len(new):
            continue
        for nd in new.node:
            if nf.get(nd, 0) > 0:
                continue                         # the new node has templates of its own
            prev = g[(g.start < r.date) & g.node.map(lambda x: nf.get(x, 0) > 0)]
            src = prev.node.iloc[-1] if len(prev) else ""
            if r.verdict.startswith("above control"):
                kind = "carry"
            elif r.verdict.startswith("AT CONTROL"):
                kind = "pending_retpl"
            else:
                kind = "untested"
            cm.append(dict(netsta=r.netsta, node=nd, start=r.date, kind=kind, src_node=src, w2_verdict=r.verdict,
                           cc_after=r.cc_after, ctl95_after=r.ctl95_after))
    C = pd.DataFrame(cm)
    C.to_csv(f"{SU}/tpl_nodemap.csv", index=False)
    # B014 2023 check (W2 monthly stack tests before / after its 2023-08 .. 2024-08 outage); no node change exists
    M = pd.read_csv(f"{P.RES}/wa_pre/w2/r2_health_monthly.csv")
    b = M[(M.netsta == "PB.B014") & M.status.isin(["ok", "control_level"])]
    pre = b[(b.month >= "2021-01") & (b.month <= "2023-08")]; post = b[b.month >= "2024-08"]
    chk = dict(pre_tested=len(pre), pre_ok=int((pre.status == "ok").sum()), post_tested=len(post),
               post_ok=int((post.status == "ok").sum()),
               rule="re-template only if the post-outage share of 'ok' tested months is below the pre-outage share minus 0.25")
    chk["pre_share"] = chk["pre_ok"] / max(chk["pre_tested"], 1); chk["post_share"] = chk["post_ok"] / max(chk["post_tested"], 1)
    chk["retemplate"] = bool(chk["post_share"] < chk["pre_share"] - 0.25)
    json.dump(chk, open(f"{SU}/b014_check.json", "w"), indent=1)
    P.log(f"[templates] library {len(L)} (family, node) templates on {L.node.nunique()} nodes {L.tag.value_counts().to_dict()}; "
          f"node map {C.kind.value_counts().to_dict() if len(C) else {}}; B014 check {chk}")


def tpl_nodes_ok():
    """{node: (kind, src_node)} usable in the W4a pass."""
    L = pd.read_csv(f"{SU}/tpl_library.csv", usecols=["fid", "node"])
    nf = L.groupby("node").fid.nunique()
    out = {nd: ("own", nd, int(n)) for nd, n in nf.items()}
    C = pd.read_csv(f"{SU}/tpl_nodemap.csv") if os.path.exists(f"{SU}/tpl_nodemap.csv") else pd.DataFrame()
    for r in C.itertuples():
        if r.kind == "carry" and isinstance(r.src_node, str) and r.src_node:
            out[r.node] = ("carry", r.src_node, int(nf.get(r.src_node, 0)))
    return out


# ------------------------------------------------------------------------------------------------ fixed rate sets
def node_at(N, ns, date):
    g = N[(N.netsta == ns) & (N.start <= date) & (N.end > date)]
    if not len(g):
        return None
    pref = {"HH": 0, "BH": 1, "EH": 2}
    g = g.assign(p=g.band.map(pref).fillna(9), z=g.comps.astype(str).str.contains("Z")).sort_values(["z", "p", "start"],
                                                                                                  ascending=[False, True, False])
    return g.node.iloc[0]


def core(a):
    S = pd.read_csv(f"{SU}/stations.csv")
    R, tag = raw_inventory()
    D2 = pd.DataFrame()
    if tag == "nodl2":       # dl2 not yet readable: its manifest says which days will exist (status data / linked)
        M = pd.read_csv(f"{P.RES}/wa_dl2/manifest.csv", usecols=["netsta", "date", "status"], dtype=str)
        D2 = M[M.status.isin(["data", "linked"])][["netsta", "date"]]
    have = pd.concat([R[["netsta", "date"]], D2]).drop_duplicates()
    N = nodes_table()
    ok = tpl_nodes_ok()
    X = excl_table()
    xs = set(zip(X.netsta, X.date))
    perm = S[~S.temp & ~S.no_template].netsta.tolist()
    res, rows = {}, []
    for era, (e0, e1) in P.ERAS.items():
        ed = pd.date_range(e0, e1, freq="D")
        nd_all = len(ed)
        members = []
        for ns in perm:
            dd = sorted(have[(have.netsta == ns) & (have.date >= e0) & (have.date <= e1)].date)
            dd = [d for d in dd if (ns, d) not in xs]
            if ns == "UW.LRIV" and era == "E3":
                rows.append(dict(era=era, netsta=ns, raw_frac=len(dd) / nd_all, tpl_frac=np.nan, member=False,
                                 why="LRIV after 2019-05-30 out of the rate network (plan D10, SVI policy)")); continue
            if len(dd) < 0.8 * nd_all:
                rows.append(dict(era=era, netsta=ns, raw_frac=len(dd) / nd_all, tpl_frac=np.nan, member=False,
                                 why="raw < 80 % of era days")); continue
            # node per month (first raw day of the month) -> template source
            mon = {}
            for d in dd:
                mon.setdefault(d[:7], d)
            good = 0; tot = 0
            for m, d in mon.items():
                nd = node_at(N, ns, d)
                n_m = sum(1 for x in dd if x[:7] == m)
                tot += n_m
                k = ok.get(nd)
                if k is not None and k[2] >= 30:
                    good += n_m
            tf = good / max(tot, 1)
            mem = tf >= 0.95
            rows.append(dict(era=era, netsta=ns, raw_frac=len(dd) / nd_all, tpl_frac=tf, member=mem,
                             why="" if mem else "W4a template source < 95 % of raw days"))
            if mem:
                members.append(ns)
        res[era] = sorted(members)
    T = pd.DataFrame(rows)
    T.to_csv(f"{SU}/core_set_table.csv", index=False)
    json.dump(dict(**res, rule="permanent, raw >= 80 % of era days (W2 dead days not counted), W4a template source "
                   "(own node >= 30 families or W2-approved carry-over) on >= 95 % of raw days, not XU.W060, LRIV out in E3",
                   raw_inventory=tag), open(f"{SU}/core_sets.json", "w"), indent=1)
    P.log(f"[core] {[(k, len(v)) for k, v in res.items()]}: " + "; ".join(f"{k}: {','.join(v)}" for k, v in res.items()))


# ------------------------------------------------------------------------------------------------ exclusions
def excl_table():
    p = f"{SU}/exclusions.csv"
    if os.path.exists(p):
        return pd.read_csv(p, dtype=str)
    M = pd.read_csv(f"{P.RES}/wa_pre/w2/r2_health_monthly.csv")
    dead = M[M.status.isin(["dead", "inherit_dead"]) | M.status_eff.astype(str).isin(["dead", "inherit_dead"])]
    rows = []
    for r in dead.itertuples():
        for d in pd.date_range(f"{r.month}-01", pd.Period(r.month, "M").end_time.normalize(), freq="D"):
            rows.append((r.netsta, str(d.date()), f"W2 month {r.status}"))
    B = pd.read_csv(f"{P.RES}/wa_pre/w2/r2_bad_for_dvv.csv")
    B = B[B.reason == "dead"]          # W2 "flat" (>= 1 % repeated samples) spans decades of PB borehole data: not a fault
    for r in B.itertuples():
        for d in pd.date_range(r.start, r.end, freq="D"):
            rows.append((r.netsta, str(d.date()), f"W2 day {r.reason}"))
    X = pd.DataFrame(rows, columns=["netsta", "date", "reason"]).drop_duplicates(["netsta", "date"])
    P.ensure(SU)
    X.to_csv(p, index=False)
    P.log(f"[excl] {len(X)} station-days excluded ({X.reason.value_counts().to_dict()})")
    return X.astype(str)


def excl(a):
    excl_table()


def daylists(a):
    """day lists for the queue; smoke runs restrict them with WA_PROD_TRUTH_DAYS / WA_PROD_PROD_DAYS (comma lists)."""
    env = lambda k: [x for x in os.environ.get(k, "").split(",") if x]
    D = pd.read_csv(f"{SU}/day_table.csv", dtype=str)
    truth = env("WA_PROD_TRUTH_DAYS") or D[D.month.isin(P.TRUTH_MONTHS)].date.tolist()
    prod = env("WA_PROD_PROD_DAYS") or D.date.tolist()
    json.dump(sorted(truth), open(f"{SU}/days_truth.json", "w"))
    json.dump(sorted(prod), open(f"{SU}/days_prod.json", "w"))
    P.log(f"[daylists] truth {len(truth)} days, prod {len(prod)} days")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["freeze", "verify", "roster", "families", "truth", "days", "rawinv", "templates",
                                     "core", "excl", "daylists", "all"])
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args()
    if a.what == "all":
        for f in (freeze, roster, families, truth, day_table, excl, rawinv, templates, core):
            f(a)
    else:
        dict(freeze=freeze, verify=verify, roster=roster, families=families, truth=truth, days=day_table, rawinv=rawinv,
             templates=templates, core=core, excl=excl, daylists=daylists)[a.what](a)
