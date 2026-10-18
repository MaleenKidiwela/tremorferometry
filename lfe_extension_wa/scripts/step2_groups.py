#!/usr/bin/env python
"""STEP 2 (Merlin, 2026-09-30) part 1 — <= 3-km GROUPS of families, member alignment, group setup, group templates.

Families: the 329 representatives used in T1 (results/merlin_tests/t1_setup/families.csv; NW, SW, Olympic Arm A, Arm B).

GROUP RULE (stated, fixed before any scoring): COMPLETE-LINKAGE agglomerative clustering of the epicentral
distances (lat_use / lon_use = location, or the approximate position for unlocated families), cut at 3.0 km
=> every pair of families inside a group is <= 3 km apart (no chaining across an area).  A family with no
position at all is its own group.  Groups may mix regions (the region of a group = its most common member
region, ties -> the reference member's).

MEMBER ALIGNMENT (deviation from "align on each member's own template offset", data-driven): the S picks of
the >= 30-event T1 stacks are noise-dominated at most stations (the within-group spread of own-offset differences
at shared stations is ~9 s), so members are aligned by cross-correlating their build-half LONG stacks
(results/templates/stacks_build.npz, 1-8 Hz, [t-20, t+32] s) with the group's reference member (most build events):
  delta_m   = common lag maximising the sum over shared station-nodes (both >= 10 events) of the 3-C normalised
              cross-correlation, lags +-10 s, window [t-10, t+22] s of the member stack
              (member event at t_m  <->  group-frame origin t_m + delta_m)
  tau_m,s   = per station-node residual: argmax of that station's cross-correlation within delta_m +- 1 s
              (member event cut at t_m + tau_m,s at station s); no stack >= 10 events -> tau = delta_m
The reference member has delta = tau = 0, so a singleton group is IDENTICAL to its family (check in the tpl stage).
  Co-detection prior: when >= 20 of the member's array detections have a reference-member detection within
  +-8 s and the 0.25-s histogram mode of those time differences holds >= max(5, 3x the uniform expectation n/64), the delta search is
  restricted to that mode +- 1 s (on the 31 members with >= 20 co-detections the unrestricted stack estimate
  agreed with the mode to <= 0.5 s for 65 %; where it did not, the stack curve was flat / multi-peaked).
Members without any >= 10-event stack keep delta = the co-detection mode if usable, else 0 (flag `aligned`).

GROUP TRUTH: the union of the members' array detections in the group frame (t_m + delta_m), merged within +-2 s
(start-anchored windows; the merged event lists every member that detected it).  half = wa_common.half_of(group t).
SCAN DAYS: a group-day of window w is scanned when EVERY member whose val_windows contains w was scanned that day.

GROUP TEMPLATES (stage tpl): per group and permanent station-epoch node, ALL members' build-half (V2008 + V2010,
half 0 in the member frame) array detections, merged co-detections counted once, cut at t_m + tau_m,s, long window
stacked with each event unit-max normalised (wa_templates rule); template = wa_templates.make_templates (S pick,
12-s window, clear flag, S-SNR); w = the e5r SNR weight normalised by the T1 build templates' clear median
(so group and family weights are on one scale).  Every >= 30-event template is kept (clear flag = a column).

  python scripts/step2_groups.py --stage groups
  python scripts/step2_groups.py --stage tpl --workers 12
  -> results/step2_groups/{groups.csv, members.csv, align_station.parquet, setup/*}, results/templates/*_step2grp.*
"""
from __future__ import annotations
import os, sys, argparse, time, json
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

G = f"{W.RES}/step2_groups"
SETUP_OUT = f"{G}/setup"
T1S = f"{W.RES}/merlin_tests/t1_setup"
LOG = f"{W.STUDY}/logs/step2.log"
TAG = "step2grp"
GID0 = 1000
RADIUS_KM = 3.0
MIN_ALIGN_EV = 10
LAG_MAX = int(10 * W.FS)          # +-10 s
RES_MAX = int(1.0 * W.FS)         # +-1 s per-station residual
WIN = (int(10 * W.FS), int(42 * W.FS))   # [t-10, t+22] s inside the [t-20, t+32] long window
NL = int(W.LONG_LEN * W.FS)
IPRE = int(W.LONG_PRE * W.FS)


def make_groups(F):
    from scipy.cluster.hierarchy import linkage, fcluster
    from scipy.spatial.distance import squareform
    ok = F.lat_use.notna().to_numpy() & F.lon_use.notna().to_numpy()
    la, lo = F.lat_use.to_numpy(float), F.lon_use.to_numpy(float)
    i = np.nonzero(ok)[0]
    D = W.km(la[i][:, None], lo[i][:, None], la[i][None], lo[i][None])
    np.fill_diagonal(D, 0.0)
    lab = fcluster(linkage(squareform(D, checks=False), "complete"), RADIUS_KM, "distance")
    g = np.full(len(F), -1)
    g[i] = lab
    nxt = lab.max() + 1 if len(lab) else 1
    for j in np.nonzero(~ok)[0]:
        g[j] = nxt; nxt += 1
    return g


def xcorr_curve(a, b):
    """3-C normalised cc of member stack a (3, NL) window WIN vs reference stack b shifted by lag l:
    c(l) = sum a(u) b(u - l) / norms, l = -LAG_MAX..LAG_MAX  (a(u) ~ b(u - delta))."""
    from numpy.lib.stride_tricks import sliding_window_view
    a0, a1 = WIN
    L = a1 - a0
    sa = a[:, a0:a1].astype(np.float64)
    na = np.sqrt((sa ** 2).sum())
    bb = b[:, a0 - LAG_MAX:a1 + LAG_MAX].astype(np.float64)       # starts 0..2*LAG_MAX  <-> l = LAG_MAX - start
    Vw = sliding_window_view(bb, L, axis=1)                        # (3, 2*LAG_MAX+1, L)
    num = np.einsum("ck,csk->s", sa, Vw)
    e = np.concatenate([np.zeros((bb.shape[0], 1)), np.cumsum(bb ** 2, 1)], 1)
    eb = (e[:, L:] - e[:, :-L]).sum(0)[:2 * LAG_MAX + 1]
    with np.errstate(invalid="ignore", divide="ignore"):
        c = num / (na * np.sqrt(eb))
    c = np.nan_to_num(c)[::-1]                                     # index k <-> l = k - LAG_MAX
    return c


def codet_mode(a, b):
    """mode (0.25-s histogram) of (nearest reference time - member time) within +-8 s; (mode, n, n_in_mode_bin)"""
    if len(a) == 0 or len(b) == 0:
        return np.nan, 0, 0
    j = np.clip(np.searchsorted(b, a), 1, max(len(b) - 1, 1))
    cand = np.stack([b[j - 1] - a, b[np.minimum(j, len(b) - 1)] - a], 1)
    d = cand[np.arange(len(a)), np.argmin(np.abs(cand), 1)]
    d = d[np.abs(d) <= 8]
    if len(d) < 5:
        return np.nan, len(d), 0
    h, e = np.histogram(d, bins=np.arange(-8, 8.01, 0.25))
    return float(e[np.argmax(h)] + 0.125), len(d), int(h.max())


def stage_groups():
    t0 = time.time()
    os.makedirs(SETUP_OUT, exist_ok=True)
    F = pd.read_csv(f"{T1S}/families.csv")
    F["g"] = make_groups(F)
    Qb = pd.read_csv(f"{W.RES}/templates/tpl_qc_build.csv")
    Z = np.load(f"{W.RES}/templates/stacks_build.npz", allow_pickle=False)
    keys = [str(k) for k in Z["keys"]]
    STK = Z["stacks"]
    kidx = {}
    for i, k in enumerate(keys):
        f, ns, nd = k.split("|")
        kidx[(int(f), ns, nd)] = i
    nev = {(int(f), ns, nd): int(n) for f, ns, nd, n in zip(Qb.fid, Qb.netsta, Qb.node, Qb.n_events)}
    T = pd.read_parquet(f"{T1S}/truth.parquet")
    Tg = {f: np.sort(g.t.to_numpy()) for f, g in T.groupby("fid")}
    nbuild = T[(T.half == 0) & T.window.isin(["V2008", "V2010"])].groupby("fid").size()
    F["n_build"] = F.fid.map(nbuild).fillna(0).astype(int)
    # group ids: ordered by region, size, latitude
    gg = F.groupby("g").agg(n=("fid", "size"), lat=("lat_use", "mean"))
    F["reg0"] = F.region
    order = sorted(gg.index, key=lambda x: (-gg.n[x], np.nan_to_num(gg.lat[x])))
    gid_of = {g: GID0 + k for k, g in enumerate(order)}
    F["gid"] = F.g.map(gid_of)
    mem_rows, st_rows, grp_rows = [], [], []
    for gid, gf in F.groupby("gid"):
        gf = gf.sort_values(["n_build", "n_detections"], ascending=False)
        ref = int(gf.fid.iloc[0])
        regs = gf.region.value_counts()
        region = regs.index[0] if (len(regs) == 1 or regs.iloc[0] > regs.iloc[1]) else gf.region.iloc[0]
        la, lo = gf.lat_use.mean(), gf.lon_use.mean()
        dmax = float(np.nanmax(W.km(gf.lat_use.to_numpy()[:, None], gf.lon_use.to_numpy()[:, None],
                                    gf.lat_use.to_numpy()[None], gf.lon_use.to_numpy()[None]))) if len(gf) > 1 else 0.0
        vws = sorted({w for v in gf.val_windows for w in v.split(",")}, key=["V2008", "V2010", "A2011"].index)
        grp_rows.append(dict(fid=gid, family_uid=f"G{gid}", region=region, regions=",".join(sorted(set(gf.region))),
                             det_source=f"G{gid}", val_windows=",".join(vws), lat_use=la, lon_use=lo,
                             n_members=len(gf), members=";".join(str(int(x)) for x in gf.fid),
                             member_uids=";".join(gf.family_uid), ref_fid=ref, max_pair_km=dmax,
                             n_unlocated=int((~gf.located.astype(bool)).sum()),
                             n_quality_C=int((gf.location_quality == "C").sum()),
                             n_build_events_sum=int(gf.n_build.sum())))
        for r in gf.itertuples():
            f = int(r.fid)
            row = dict(gid=gid, fid=f, family_uid=r.family_uid, region=r.region, is_ref=f == ref, n_build=int(r.n_build),
                       delta_s=0.0, aligned=f == ref, n_align_sta=0, align_peak=np.nan, align_peak2=np.nan,
                       codet_mode_s=np.nan, n_codet_8s=0, n_codet_mode_bin=0, codet_prior=False)
            if f != ref:
                cm_, cn_, cb_ = codet_mode(Tg.get(f, np.zeros(0)), Tg.get(ref, np.zeros(0)))
                use_prior = bool(cn_ >= 20 and cb_ >= max(5, 3.0 * cn_ / 64))
                row.update(codet_mode_s=cm_, n_codet_8s=cn_, n_codet_mode_bin=cb_, codet_prior=use_prior)
                if use_prior:
                    row.update(delta_s=cm_, aligned=True)
                common = [(ns, nd) for (ff, ns, nd) in kidx if ff == f and nev.get((f, ns, nd), 0) >= MIN_ALIGN_EV
                          and nev.get((ref, ns, nd), 0) >= MIN_ALIGN_EV and (ref, ns, nd) in kidx]
                if common:
                    curves = np.stack([xcorr_curve(STK[kidx[(f, ns, nd)]], STK[kidx[(ref, ns, nd)]]) for ns, nd in common])
                    C = curves.sum(0)
                    if use_prior:
                        kc = int(round(cm_ * W.FS)) + LAG_MAX
                        lo_, hi_ = max(0, kc - RES_MAX), min(len(C), kc + RES_MAX + 1)
                        k0 = lo_ + int(np.argmax(C[lo_:hi_]))
                    else:
                        k0 = int(np.argmax(C))
                    far = np.abs(np.arange(len(C)) - k0) > RES_MAX
                    row.update(delta_s=(k0 - LAG_MAX) / W.FS, aligned=True, n_align_sta=len(common),
                               align_peak=float(C[k0] / len(common)),
                               align_peak2=float(C[far].max() / len(common)) if far.any() else np.nan)
                    for (ns, nd), c in zip(common, curves):
                        lo_, hi_ = max(0, k0 - RES_MAX), min(len(c), k0 + RES_MAX + 1)
                        k1 = lo_ + int(np.argmax(c[lo_:hi_]))
                        st_rows.append(dict(gid=gid, fid=f, netsta=ns, node=nd, tau_s=(k1 - LAG_MAX) / W.FS,
                                            cc_at_delta=float(c[k0]), cc_at_tau=float(c[k1]),
                                            n_member=nev[(f, ns, nd)], n_ref=nev[(ref, ns, nd)]))
            mem_rows.append(row)
    GR = pd.DataFrame(grp_rows)
    M = pd.DataFrame(mem_rows)
    ST = pd.DataFrame(st_rows, columns=["gid", "fid", "netsta", "node", "tau_s", "cc_at_delta", "cc_at_tau",
                                        "n_member", "n_ref"])
    M["delta_minus_codet_s"] = M.delta_s - M.codet_mode_s
    M["align_ratio"] = M.align_peak / M.align_peak2        # stack cc peak over the best peak > 1 s away (1 = flat)
    # ---- group truth (union, group frame, merged +-2 s)
    dl = dict(zip(M.fid, M.delta_s)); g_of = dict(zip(M.fid, M.gid))
    pr = {f: k for k, f in enumerate(M.sort_values(["gid", "is_ref", "n_build"], ascending=[True, False, False]).fid)}
    T = T.copy()
    T["gid"] = T.fid.map(g_of); T["t_m"] = T.t; T["t_g"] = T.t + T.fid.map(dl)
    T["prio"] = T.fid.map(pr)
    rows = []
    for gid, g in T.groupby("gid"):
        g = g.sort_values(["t_g", "prio"])
        tg = g.t_g.to_numpy(); fm = g.fid.to_numpy(); pp = g.prio.to_numpy(); tm = g.t_m.to_numpy()
        wn = g.window.to_numpy()
        i = 0
        while i < len(g):
            j = i + 1
            while j < len(g) and tg[j] - tg[i] <= W.MATCH_S:
                j += 1
            k = i + int(np.argmin(pp[i:j]))
            mem = sorted(set(int(x) for x in fm[i:j]))
            rows.append((int(gid), float(tg[k]), int(fm[k]), float(tm[k]), ";".join(map(str, mem)), len(mem), j - i, wn[k]))
            i = j
    TG = pd.DataFrame(rows, columns=["fid", "t", "member", "t_member", "members", "n_members_det", "n_raw", "window"])
    TG["date"] = pd.to_datetime(TG.t, unit="s").dt.strftime("%Y-%m-%d")
    TG["half"] = W.half_of(TG.t.to_numpy())
    TG["det_source"] = "G" + TG.fid.astype(str)
    # ---- scan days
    S1 = pd.read_csv(f"{T1S}/scan_days.csv")
    sd = {s: set(g.date) for s, g in S1.groupby("det_source")}
    Fi = F.set_index("fid")
    scan = []
    for r in GR.itertuples():
        mem = [int(x) for x in r.members.split(";")]
        for w in r.val_windows.split(","):
            mw = [m for m in mem if w in Fi.val_windows[m].split(",")]
            days = set(d for d, _ in W.window_days(w))
            for m in mw:
                days &= sd.get(Fi.det_source[m], set())
            scan += [(r.det_source, d) for d in sorted(days)]
    SC = pd.DataFrame(scan, columns=["det_source", "date"]).drop_duplicates()
    GR.to_csv(f"{SETUP_OUT}/families.csv", index=False)
    GR.to_csv(f"{G}/groups.csv", index=False)
    M.to_csv(f"{G}/members.csv", index=False)
    ST.to_parquet(f"{G}/align_station.parquet", index=False)
    TG.to_parquet(f"{SETUP_OUT}/truth.parquet", index=False)
    SC.to_csv(f"{SETUP_OUT}/scan_days.csv", index=False)
    import shutil
    shutil.copy(f"{T1S}/stations.csv", f"{SETUP_OUT}/stations.csv")
    sz = GR.n_members.value_counts().sort_index().to_dict()
    W.log(LOG, f"[groups] {len(F)} families -> {len(GR)} groups, sizes {sz}; multi-member groups "
          f"{int((GR.n_members > 1).sum())} holding {int(GR.n_members[GR.n_members > 1].sum())} families; "
          f"aligned non-ref members {int((M.aligned & ~M.is_ref).sum())}/{int((~M.is_ref).sum())}; "
          f"truth raw {len(T)} -> merged {len(TG)}; scan rows {len(SC)}; {time.time()-t0:.0f}s")
    a = M[~M.is_ref]
    W.log(LOG, f"[groups] non-ref members {len(a)}: co-detection prior used {int(a.codet_prior.sum())}, stack-only "
          f"{int((a.aligned & ~a.codet_prior).sum())}, unaligned (delta 0) {int((~a.aligned).sum())}")


# ------------------------------------------------------------------ templates
def stack_day(args):
    date, ev, stations, tau = args
    D0 = pd.Timestamp(date).timestamp()
    out = {}
    fid = ev.fid.to_numpy(); gid = ev.gid.to_numpy(); tm = ev.t_m.to_numpy(); dlt = ev.delta.to_numpy()
    for ns in stations:
        p = W.cache_path(ns, date)
        if not os.path.exists(p):
            continue
        z = np.load(p)
        node = str(z["node"])
        X, m, live = z["x"], z["m"], z["live"]
        sh = np.array([tau.get((int(f), ns, node), d) for f, d in zip(fid, dlt)])
        i0 = np.round((tm + sh - D0) * W.FS).astype(int) - IPRE
        ok = (i0 >= 0) & (i0 + NL <= W.NDAY)
        if not ok.any():
            continue
        cm = np.concatenate([[0], np.cumsum(~m)])
        ii = i0[ok]
        ok2 = (cm[ii + NL] - cm[ii]) == 0
        gg = gid[ok][ok2]; ff = fid[ok][ok2]
        ii = ii[ok2]
        if len(ii) == 0:
            continue
        idx = ii[:, None] + np.arange(NL)[None]
        seg = X[:, idx].transpose(1, 0, 2).astype(np.float64)
        mx = np.abs(seg).max(axis=(1, 2))
        good = mx > 0
        seg = seg[good] / mx[good, None, None]
        gg = gg[good]; ff = ff[good]
        for g in np.unique(gg):
            s = seg[gg == g]
            out[(int(g), ns, node)] = (s.sum(0).astype(np.float32), len(s), live.copy(), set(int(x) for x in ff[gg == g]))
    return date, out


def stage_tpl(workers, days=None):
    import wa_templates as WT
    t0 = time.time()
    M = pd.read_csv(f"{G}/members.csv")
    ST = pd.read_parquet(f"{G}/align_station.parquet")
    tau = {(int(f), ns, nd): float(t) for f, ns, nd, t in zip(ST.fid, ST.netsta, ST.node, ST.tau_s)}
    TG = pd.read_parquet(f"{SETUP_OUT}/truth.parquet")
    # build-half union events (merged co-detections once): the representative member's own detection
    EV = TG[TG.window.isin(["V2008", "V2010"])].copy()
    T1 = pd.read_parquet(f"{T1S}/truth.parquet")
    hb = dict(zip(zip(T1.fid, T1.t), T1.half))
    EV["half_m"] = [hb.get((int(f), t), -9) for f, t in zip(EV.member, EV.t_member)]
    EV = EV[EV.half_m == 0]
    EV = pd.DataFrame(dict(gid=EV.fid.astype(int), fid=EV.member.astype(int), t_m=EV.t_member))
    EV["delta"] = EV.fid.map(dict(zip(M.fid, M.delta_s)))
    EV["date"] = pd.to_datetime(EV.t_m, unit="s").dt.strftime("%Y-%m-%d")
    if days:
        EV = EV[EV.date.isin(days)]
    S = pd.read_csv(f"{SETUP_OUT}/stations.csv")
    W.log(LOG, f"[tpl] {len(EV)} union build events on {EV.date.nunique()} days, {EV.gid.nunique()} groups")
    acc, contrib = {}, {}
    with ProcessPoolExecutor(workers) as ex:
        fut = [ex.submit(stack_day, (d, g, list(S.netsta), tau)) for d, g in EV.groupby("date")]
        for i, fu in enumerate(as_completed(fut)):
            d, out = fu.result()
            for k, (s, n, live, mem) in out.items():
                if k in acc:
                    s0, n0, l0 = acc[k]; acc[k] = (s0 + s, n0 + n, l0 | live); contrib[k] |= mem
                else:
                    acc[k] = (s.astype(np.float64), n, live); contrib[k] = set(mem)
            if (i + 1) % 10 == 0:
                W.log(LOG, f"[tpl] {i+1}/{len(fut)} days {time.time()-t0:.0f}s")
    Q, TPL = WT.make_templates(acc, S)
    Q["clear"] = Q.clear.fillna(False).astype(bool)
    Q["n_members_contrib"] = [len(contrib[(f, ns, nd)]) for f, ns, nd in zip(Q.fid, Q.netsta, Q.node)]
    R = pd.read_csv(f"{W.RES}/templates/tpl_qc_build.csv")
    R = R[R.clear.fillna(False).astype(bool)]
    ref = float(np.median(10 ** (R.snr_db / 20.0)))
    Q["w"] = np.minimum(10 ** (Q.snr_db / 20.0), 3 * ref) / ref
    # the most events any single member has at that node (T1 build stacks) -> "more events than any one family"
    Qb = pd.read_csv(f"{W.RES}/templates/tpl_qc_build.csv")
    nb = {(int(f), ns, nd): int(n) for f, ns, nd, n in zip(Qb.fid, Qb.netsta, Qb.node, Qb.n_events)}
    mem = {int(g): [int(x) for x in s.split(";")] for g, s in pd.read_csv(f"{SETUP_OUT}/families.csv")[["fid", "members"]].values}
    Q["n_events_best_member"] = [max([nb.get((m, ns, nd), 0) for m in mem[int(g)]]) for g, ns, nd in zip(Q.fid, Q.netsta, Q.node)]
    tag = TAG if not days else TAG + "_smoke"
    od = f"{W.RES}/templates"
    Q.to_csv(f"{od}/tpl_qc_{tag}.csv", index=False)
    keys = list(TPL)
    np.savez_compressed(f"{od}/tpl_{tag}.npz", keys=np.array([f"{f}|{ns}|{nd}" for f, ns, nd in keys]),
                        tpl=np.stack([TPL[k] for k in keys]) if keys else np.zeros((0, 3, W.NT), np.float32))
    np.savez_compressed(f"{od}/stacks_{tag}.npz", keys=np.array([f"{f}|{ns}|{nd}" for f, ns, nd in acc]),
                        stacks=np.stack([(v[0] / v[1]).astype(np.float32) for v in acc.values()]) if acc else np.zeros((0, 3, NL)))
    use = (Q.n_events >= W.MIN_EVENTS) & Q.offset_s.notna()
    W.log(LOG, f"[tpl:{tag}] (group,node) stacks {len(Q)}; >=30-event templates {int(use.sum())} "
          f"({Q[use].fid.nunique()} groups); clear {int(Q.clear.sum())}; {time.time()-t0:.0f}s")
    # singleton identity check vs the T1 build templates
    GR = pd.read_csv(f"{SETUP_OUT}/families.csv")
    single = dict(zip(GR[GR.n_members == 1].fid, GR[GR.n_members == 1].members.astype(int)))
    zb = np.load(f"{od}/tpl_build.npz", allow_pickle=False)
    zbk, zbt = zb["keys"], zb["tpl"]              # read once (NpzFile re-decompresses on every access)
    bk = {str(k): i for i, k in enumerate(zbk)}
    diffs, noff = [], 0
    ob = {(int(f), ns, nd): o for f, ns, nd, o in zip(Qb.fid, Qb.netsta, Qb.node, Qb.offset_s)}
    og = {(int(f), ns, nd): o for f, ns, nd, o in zip(Q.fid, Q.netsta, Q.node, Q.offset_s)}
    for i, (g, ns, nd) in enumerate(keys):
        if g in single:
            k = f"{single[g]}|{ns}|{nd}"
            if k in bk:
                diffs.append(float(np.abs(zbt[bk[k]] - TPL[(g, ns, nd)]).max()))
                noff += int(abs(og[(g, ns, nd)] - ob[(single[g], ns, nd)]) > 1e-9)
    W.log(LOG, f"[tpl:{tag}] singleton check: {len(diffs)} templates compared, max |diff| "
          f"{max(diffs) if diffs else float('nan'):.2e}, offset mismatches {noff}")
    json.dump(dict(n_compared=len(diffs), max_abs_diff=max(diffs) if diffs else None, offset_mismatch=noff),
              open(f"{G}/singleton_template_check{'_smoke' if days else ''}.json", "w"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["groups", "tpl"], required=True)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--days", default="", help="smoke: comma list of dates (tpl stage; writes *_smoke)")
    a = ap.parse_args()
    os.makedirs(G, exist_ok=True)
    if a.stage == "groups":
        stage_groups()
    else:
        stage_tpl(a.workers, a.days.split(",") if a.days else None)


if __name__ == "__main__":
    main()
