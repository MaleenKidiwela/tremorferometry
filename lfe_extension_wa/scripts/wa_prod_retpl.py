#!/usr/bin/env python
"""WA PRODUCTION W4b: templates for the station epochs the W4a library cannot serve, AFTER the confidence model exists.

Targets   (a) the epochs at control level after an equipment change (W2: PGC 2015-01-12, VGZ 2017-03-03, DOSE 2019-05-07,
              HDW 2020-12-04) -- re-templated, not carried over (user 2026-10-03);
          (b) every other (station, node) seen in the production caches with no W4a template (new 3-C UW stations, node
              changes of tier B / C stations, temporary networks after 2012), <= 50 km from a family; XU.W060 never.
Times     T2B rule (results/T2B_RESULTS.md): production candidates of the network WITHOUT the node (W4a 'all' network, in
          which the node has no template) with P_own_WAVE >= 0.3, families <= 50 km, on the node's build days = its first
          two WA ETS (tremor_episodes.csv, PNSN) with data, else its first 24 months.
Stack     the frozen wa_templates recipe ([t-20, t+32] s of the 20 Hz cache, unit max over 3 components, fully valid
          windows, >= 30 events, S pick from the 0.5-s envelope, 12-s template from S - 0.5 s, per-channel L2 norm, SNR weight
          normalised by the T1 build median).
Check     per channel (DUC_FIN validated rule, lfe_south/scripts/duc_fix_retemplate.py): coherence S = rms(stack) x sqrt(N)
          of the RMS-normalised 12-s segments; 20 scrambled controls = the same event times on days shifted by +-3..40 days
          inside the node's life; keep the channel only if N >= 20, >= 10 controls exist and S > every control; else the
          channel is zeroed (no channel kept -> the template is not used; its row stays).  T2B column: cc (+-2 s) between the
          template and the stack at the SCRAMBLED-CONTROL candidates with P >= 0.3 (control series c1).
Outputs   OUT/templates/tpl_w4b.npz, tpl_qc_w4b.csv, w4b_targets.csv, w4b_build_days.csv (in-sample days per node),
          rich_days.json (days with a used W4b node live -> the rich pass).

  python scripts/wa_prod_retpl.py targets | build --workers 12 | richdays
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, warnings, zlib
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("WA_EXTRA_NODES", "/home/jovyan/tremorferometry/lfe_extension_wa/results/step4b/setup/extra_nodes.csv")
import wa_prod_common as P
W = P.patch_wa()
import wa_templates as WT

OT = f"{P.OUT}/templates"
DMAX_KM = 50.0
NCTRL, NMIN_CH = 20, 20
SEG0 = WT.IPRE                     # index of t in the long window


def node_days():
    """(netsta, date, node) of every cached production station-day (manifests of every WA-chain cache)."""
    ms = glob.glob(f"{P.OUT}/cache/manifest_p*_*.csv") + [f"{P.RES}/cache_manifest.csv", f"{P.RES}/wa_pre/w2/cache_manifest.csv"]
    parts = []
    for m in ms:
        if os.path.exists(m):
            M = pd.read_csv(m, usecols=lambda c: c in ("netsta", "date", "status", "node"))
            M = M[M.status.isin(["ok", "cached"])]
            if "node" not in M.columns:      # 2026-10-06 fix: a month whose manifest has only 'empty' rows is written without 'node'
                if len(M): raise RuntimeError(f"{m}: ok rows but no node column")
                continue
            parts.append(M[["netsta", "date", "node"]])
    N = pd.concat(parts, ignore_index=True).drop_duplicates(["netsta", "date"])
    return N[(N.date >= P.T0) & (N.date <= P.T1)]


def targets(a):
    P.ensure(OT)
    S = P.stations(); F = P.families()
    L = pd.read_csv(f"{P.OUT}/setup/tpl_library.csv", usecols=["fid", "node"])
    own = set(L.node)
    C = pd.read_csv(f"{P.OUT}/setup/tpl_nodemap.csv")
    carry = set(C[C.kind == "carry"].node)
    pending = set(C[C.kind == "pending_retpl"].node)
    N = node_days()
    N = N[~(N.netsta.str.split(".").str[0].isin(["XU", "XG"]) & (N.date < "2012-01-01"))]
    g = N.groupby(["netsta", "node"]).date.agg(first="min", last="max", ndays="count").reset_index()
    co = S.set_index("netsta")
    rows = []
    for r in g.itertuples():
        if r.netsta in P.NO_TEMPLATE or r.netsta not in co.index or r.netsta.split(".")[0] in P.BANNED_NETS:
            continue
        if r.node in own or r.node in carry:
            continue
        la, lo = co.loc[r.netsta, "lat"], co.loc[r.netsta, "lon"]
        dkm = W.km(la, lo, F.lat_use.to_numpy(), F.lon_use.to_numpy())
        fams = F.fid[dkm <= DMAX_KM].astype(int).tolist()
        if not fams:
            continue
        rows.append(dict(netsta=r.netsta, node=r.node, tkey=r.node, kind="retemplate_control_level" if r.node in pending else "new_node_t2b",
                         first=r.first, last=r.last, n_days=int(r.ndays), n_fam_50km=len(fams), fams=";".join(map(str, fams))))
    T = pd.DataFrame(rows)
    # build days: first two WA ETS with data, else the first 24 months
    E = pd.read_csv(f"{W.EXT}/tremor_episodes.csv")
    E = E[E.region.isin(["NW_WA", "SW_WA", "OLY2010"])].sort_values("start")
    ep = [(r.start, r.end) for r in E.itertuples()]
    bd = []
    nd = N.groupby(["netsta", "node"]).date.apply(lambda x: sorted(set(x))).to_dict()
    for r in T.itertuples():
        dd = [d for d in nd[(r.netsta, r.node)] if d[:7] not in P.SEALED_MONTHS]     # the sealed months are never used to build
        sel = []
        nep = 0
        for s0, s1 in ep:
            if s1 < r.first:
                continue
            x = [d for d in dd if s0 <= d <= s1]
            if len(x) >= 3:
                sel += x; nep += 1
            if nep == 2:
                break
        if nep < 2:
            lim = str((pd.Timestamp(r.first) + pd.DateOffset(months=24)).date())
            sel = [d for d in dd if d <= lim]
        for d in sorted(set(sel)):
            bd.append(dict(netsta=r.netsta, tkey=r.tkey, date=d, rule="first two WA ETS" if nep == 2 else "first 24 months"))
    B = pd.DataFrame(bd)
    T.to_csv(f"{OT}/w4b_targets.csv", index=False)
    B.to_csv(f"{OT}/w4b_build_days.csv", index=False)
    P.log(f"[w4b targets] {len(T)} nodes ({T.kind.value_counts().to_dict() if len(T) else {}}); build station-days {len(B)}")


def _seg_rms(X, i0):
    s = X[:, i0:i0 + W.NT].astype(np.float64)
    r = np.sqrt((s ** 2).mean(1))
    return s, r


def build_one(args):
    ns, node, tkey, fams, bdays, ev, ev_ctl, ndays_all = args
    warnings.filterwarnings("ignore")
    t0 = time.time()
    fams = set(fams)
    ev = ev[ev.fid.isin(fams)]
    acc = {}
    # ---- the stack (frozen recipe): per build day
    for d, g in ev.groupby("date"):
        _, out = WT.stack_day((d, g[["fid", "t"]], [ns], {node}))
        for k, (s, n, live) in out.items():
            if k in acc:
                s0, n0, l0 = acc[k]; acc[k] = (s0 + s, n0 + n, l0 | live)
            else:
                acc[k] = (s.astype(np.float64), n, live)
    if not acc:
        return [], {}, ns, node, time.time() - t0
    Sdf = pd.DataFrame(dict(netsta=[ns], site=[ns]))
    Q, TPL = WT.make_templates(acc, Sdf)
    if not len(Q) or "offset_s" not in Q:
        return [], {}, ns, node, time.time() - t0
    # ---- per-channel scrambled-control check on the template segment
    rng = np.random.default_rng(zlib.crc32(f"{ns}|{node}".encode()))
    shifts = []
    while len(shifts) < NCTRL:
        k = int(rng.integers(3, 41)) * int(rng.choice([-1, 1]))
        if k not in shifts:
            shifts.append(k)
    ndset = set(ndays_all)
    cache = {}

    def day(d):
        if d not in cache:
            if len(cache) > 48:
                cache.pop(next(iter(cache)))
            p = P.find_cache(ns, d) if d in ndset else None
            if p is None:
                cache[d] = None
            else:
                z = np.load(p)
                cache[d] = (z["x"], z["m"], str(z["node"]))
        return cache[d]

    def segs(times, off):
        """RMS-normalised 12-s segments (3, n) at the template window, fully valid, at this node."""
        out = []
        for t in times:
            d = str(pd.Timestamp(t, unit="s").date())
            v = day(d)
            if v is None or v[2] != node:
                continue
            X, m, _ = v
            i0 = int(round((t - pd.Timestamp(d).timestamp() + off) * W.FS))
            if i0 < 0 or i0 + W.NT > W.NDAY or not m[i0:i0 + W.NT].all():
                continue
            s, r = _seg_rms(X, i0)
            out.append(np.where(r[:, None] > 0, s / np.where(r > 0, r, 1)[:, None], 0.0))
        return out

    rows, keep_tpl = [], {}
    for q in Q.itertuples():
        row = q._asdict(); row.pop("Index", None)
        row.update(tkey=tkey, kind="", S_ch="", ctl_max_ch="", keep_ch="", n_ctrl=0, cc_ctl_stack=np.nan, use=False)
        if (q.fid, ns, node) not in TPL or not np.isfinite(q.offset_s):
            rows.append(row); continue
        times = ev[ev.fid == q.fid].t.to_numpy(float)
        sg = segs(times, q.offset_s)
        N = len(sg)
        live = np.array([c == "1" for c in str(q.live).zfill(3)])
        keep = np.zeros(3, bool); Sv = np.full(3, np.nan); Cm = np.full(3, np.nan); nct = 0
        if N >= 2:
            A = np.stack(sg)                                        # (N, 3, n)
            Sv = np.sqrt((A.mean(0) ** 2).mean(1)) * np.sqrt(N)
            Cs = []
            for k in shifts:
                sc = segs(times + 86400.0 * k, q.offset_s)
                if len(sc) >= 2:
                    B = np.stack(sc)
                    Cs.append(np.sqrt((B.mean(0) ** 2).mean(1)) * np.sqrt(len(sc)))
            nct = len(Cs)
            if nct:
                Cm = np.max(np.stack(Cs), 0)
                keep = live & (N >= NMIN_CH) & (nct >= 10) & (Sv > Cm)
        tp = TPL[(q.fid, ns, node)].copy()
        tp[~keep] = 0.0
        # T2B column: cc vs the stack at scrambled-control candidate times
        ctl_t = ev_ctl[ev_ctl.fid == q.fid].t.to_numpy(float)
        cs = segs(ctl_t, q.offset_s)
        if len(cs) >= 3 and keep.any():
            C_ = np.stack(cs).mean(0)
            best = -1.0
            for lag in range(-40, 41):
                vals = []
                for c in np.flatnonzero(keep):
                    u = np.roll(C_[c], lag); v = tp[c]
                    nu = np.sqrt((u * u).sum()); nv = np.sqrt((v * v).sum())
                    if nu > 0 and nv > 0:
                        vals.append(float(u @ v / nu / nv))
                if vals:
                    best = max(best, float(np.mean(vals)))
            row["cc_ctl_stack"] = best
        row.update(S_ch=";".join(f"{x:.2f}" for x in Sv), ctl_max_ch=";".join(f"{x:.2f}" for x in Cm),
                   keep_ch="".join("1" if x else "0" for x in keep), n_ctrl=nct, n_seg=N,
                   live="".join("1" if x else "0" for x in keep), use=bool(keep.any()))
        rows.append(row)
        if keep.any():
            keep_tpl[(q.fid, ns, tkey)] = tp.astype(np.float32)
    return rows, keep_tpl, ns, node, time.time() - t0


def build(a):
    T = pd.read_csv(f"{OT}/w4b_targets.csv")
    B = pd.read_csv(f"{OT}/w4b_build_days.csv", dtype=str)
    N = node_days()
    nd = N.groupby(["netsta", "node"]).date.apply(lambda x: sorted(set(x))).to_dict()
    # candidate times: production 'all' network, P_own_WAVE >= 0.3, data (and c1 for the T2B column)
    ps = sorted(glob.glob(f"{P.OUT}/conf/P/prod/P_*.parquet"))
    need = set(B.date)
    parts = []
    for p in ps:
        D = pd.read_parquet(p, columns=["date", "net", "series", "fid", "t", "P_own_WAVE"])
        D = D[(D.net == "all") & D.series.isin(["data", "c1"]) & (D.P_own_WAVE >= 0.3) & D.date.isin(need)]
        parts.append(D)
    EV = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=["date", "net", "series", "fid", "t", "P_own_WAVE"])
    bd = B.groupby(["netsta", "tkey"]).date.apply(list).to_dict()
    jobs = []
    for r in T.itertuples():
        days = set(bd.get((r.netsta, r.tkey), []))
        fams = [int(x) for x in str(r.fams).split(";") if x]
        e = EV[EV.date.isin(days)]
        jobs.append((r.netsta, r.node, r.tkey, fams, sorted(days), e[e.series == "data"][["fid", "t", "date"]],
                     e[e.series == "c1"][["fid", "t", "date"]], nd.get((r.netsta, r.node), [])))
    P.log(f"[w4b build] {len(jobs)} nodes; candidate pool {len(EV)} rows ({EV.series.value_counts().to_dict() if len(EV) else {}})")
    allrows, alltpl = [], {}
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        fut = [ex.submit(build_one, j) for j in jobs]
        for i, f in enumerate(as_completed(fut)):
            rows, tp, ns, node, secs = f.result()
            allrows += rows; alltpl.update(tp)
            P.log(f"[w4b build] {i + 1}/{len(jobs)} {ns} {node}: {len(rows)} stacks, {len(tp)} used, {secs:.0f}s")
    Q = pd.DataFrame(allrows)
    if not len(Q):
        Q = pd.DataFrame(columns=["fid", "netsta", "node", "tkey", "kind", "n_events", "live", "clear", "offset_s", "snr_db", "w",
                                  "S_ch", "ctl_max_ch", "keep_ch", "n_ctrl", "cc_ctl_stack", "use"])
    if len(Q):
        tk = dict(zip(zip(T.netsta, T.node), T.kind))
        Q["kind"] = [tk.get((s, n), "") for s, n in zip(Q.netsta, Q.node)]
        Qb = pd.read_csv(f"{P.RES}/templates/tpl_qc_build.csv")
        ref = float(np.median(10 ** (Qb[Qb.clear.fillna(False).astype(bool)].snr_db / 20.0)))
        Q["w"] = np.minimum(10 ** (Q.snr_db.astype(float) / 20.0), 3 * ref) / ref
        Q.loc[~Q.use.astype(bool), "w"] = np.nan
    Q.to_csv(f"{OT}/tpl_qc_w4b.csv", index=False)
    keys = list(alltpl)
    np.savez_compressed(f"{OT}/tpl_w4b.npz", keys=np.array([f"{f}|{ns}|{tk}" for f, ns, tk in keys]),
                        tpl=np.stack([alltpl[k] for k in keys]) if keys else np.zeros((0, 3, W.NT), np.float32))
    u = Q[Q.use.astype(bool)] if len(Q) else Q
    P.log(f"[w4b build] {len(Q)} (family, node) stacks, {len(u)} used on {u.node.nunique() if len(u) else 0} nodes; "
          f"median cc vs control stack {u.cc_ctl_stack.median() if len(u) else np.nan:.2f}; {time.time() - t0:.0f}s")


def richdays(a):
    try:
        Q = pd.read_csv(f"{OT}/tpl_qc_w4b.csv")
    except pd.errors.EmptyDataError:
        Q = pd.DataFrame(columns=["netsta", "tkey", "use"])
    used = set(zip(Q[Q.use.astype(bool)].netsta, Q[Q.use.astype(bool)].tkey)) if len(Q) else set()
    N = node_days()
    N = N[[(s, n) in used for s, n in zip(N.netsta, N.node)]]
    N = N[~(N.netsta.str.split(".").str[0].isin(["XU", "XG"]) & (N.date < "2012-01-01"))]
    R = N.groupby("date").netsta.apply(lambda x: sorted(set(x))).to_dict()
    json.dump(R, open(f"{OT}/rich_days.json", "w"))
    P.log(f"[w4b richdays] {len(R)} days with a used W4b node live ({len(used)} nodes)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["targets", "build", "richdays"])
    ap.add_argument("--workers", type=int, default=12)
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    dict(targets=targets, build=build, richdays=richdays)[a.what](a)
