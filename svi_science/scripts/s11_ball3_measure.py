#!/usr/bin/env python
"""s11_ball3_measure -- S11 ballistic round 3, measurement pass (ballistic/round3/PLAN.md §0).
Round-1/2 machinery reused (pick_S, pick_P, boot, fshift; lag core re-implemented as cc curves + the identical parabolic
peak, checked against s11_ball_measure.lag_batch). New: region switch (SVI | WA_NW), ETS-phase bin schemes PH / PH5 / PHE,
200 episode-date permutation copies of PH (PHp, no bootstrap), estimators S, S3c, P3, P3s, Pz, anchor pairs for WA,
site-joint S/P (equal-weight mean of the members' cc curves on common rows).

    python s11_ball3_measure.py --region SVI|WA_NW --workers 8 [--fams 10] [--tag _test]
"""
from __future__ import annotations
import os, sys, argparse, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_multi_common as M
from s10_measure import pick_S
import s11_ball_measure as B1
from s11_ball_measure import boot, fshift, pick_P

ROOT3 = f"{S.ROOT}/svi_science/results/s11_depth/ballistic/round3"
NMIN, NBOOT, LMAX = 20, B1.NBOOT, B1.LMAX
LAGS = np.arange(-LMAX, LMAX + 1)
NPERM = 200
SEED = 20261013
SCHEMES = ("EP", "YR", "Y2", "PH", "PH5", "PHE")
CTRL_PAIR_SCHEMES = ("YR", "PH", "EP", "DH", "Q4")
if os.environ.get("S11B3_SCHEMES"):
    SCHEMES = tuple(os.environ["S11B3_SCHEMES"].split(","))
DO_PERM = "PH" in SCHEMES
ST_KEYS = ("S", "S3c", "P3", "P3s", "Pz")
PR_KEYS = ("S", "P3")
PHCODE = {"BTW": 0, "ARR": 1, "AFT": 2}


def region_cfg(region):
    if region == "SVI":
        st = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011", "PB.B926", "CN.VGZ", "CN.SNB"]
        quart = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]
        pairs = [(a, b) for i, a in enumerate(quart) for b in quart[i + 1:]] + \
                [(a, x) for x in ["PB.B926", "CN.SNB", "CN.VGZ"] for a in quart]
        groups = {"PGCsite": quart}
        return dict(tb=M.TB, win=M.WIN, stations=st, pairs=pairs, groups=groups)
    tb = f"{S.ROOT}/svi_science/results/s11_depth/wa_nw/tables/inv"
    Sx = pd.read_csv(f"{tb}/stations.csv"); Sx = Sx[Sx.use.astype(bool)]
    st = [s for s in Sx.netsta if os.path.exists(f"{tb}/win/{s}.npz")]
    anchors = ["PB.B005", "PB.B006", "PB.B007", "UW.LRIV", "PB.B013", "PB.B943"]
    pairs = []
    for i, a in enumerate(anchors):
        for b in st:
            if b == a or (b in anchors and anchors.index(b) < i):
                continue
            pairs.append((a, b))
    groups = {"G1": ["PB.B005", "PB.B006", "PB.B007", "UW.LRIV"], "G2": ["PB.B013", "PB.B943"]}
    return dict(tb=tb, win=f"{tb}/win", stations=st, pairs=pairs, groups=groups)


G_W, G_POS, G_PERM = {}, {}, {}


# ------------------------------------------------------------------------------------------ lag core
def cc_curves(X, g, seg, proj):
    """X (m, 3, NW), g (3, NW) or (m, 3, NW) references, seg, proj (k, 3) -> cc (m, len(LAGS))"""
    xs = np.einsum("kc,mcn->mkn", proj, X)
    if g.ndim == 2:
        gs = (proj @ g)[None]
    else:
        gs = np.einsum("kc,mcn->mkn", proj, g)
    gg = gs[:, :, seg]; gg = gg - gg.mean(2, keepdims=True); dg = (gg * gg).sum((1, 2))
    m = X.shape[0]; cc = np.zeros((m, len(LAGS)))
    for j, l in enumerate(LAGS):
        x = xs[:, :, seg + l]; x = x - x.mean(2, keepdims=True)
        cc[:, j] = (x * gg).sum((1, 2)) / np.sqrt(np.maximum((x * x).sum((1, 2)) * dg, 1e-30))
    return cc


def peak(cc):
    m = cc.shape[0]
    i = np.argmax(cc, 1); r = np.arange(m); ii = np.clip(i, 1, len(LAGS) - 2)
    y0, y1, y2 = cc[r, ii - 1], cc[r, ii], cc[r, ii + 1]; d2 = y0 - 2 * y1 + y2
    delta = np.where((d2 < 0) & (i == ii), 0.5 * (y0 - y2) / np.where(d2 == 0, -1, d2), 0.0)
    rail = (i == 0) | (i == len(LAGS) - 1)
    return (LAGS[i] + np.clip(delta, -1, 1)) / S.FS, cc[r, i], rail


def geom3(G, tS, tP, SP):
    segS = np.flatnonzero((S.TAU >= tS - 0.5) & (S.TAU <= tS + 1.5))
    segS3 = np.flatnonzero((S.TAU >= tS - 0.5) & (S.TAU <= tS + 2.5))
    segP = segPs = np.zeros(0, int)
    if np.isfinite(tP):
        segP = np.flatnonzero((S.TAU >= max(tP - 0.5, S.T0W + 0.5)) & (S.TAU <= min(tP + 2.5, tS - 1.0)))
        segPs = np.flatnonzero((S.TAU >= max(tP - 0.3, S.T0W + 0.5)) & (S.TAU <= min(tP + 1.2, tS - 1.0)))
    if len(segP) < 20:
        segP = segPs = np.zeros(0, int)
    v3 = np.array([1.0, 0, 0])
    if len(segP):
        H = G[:, segP]; H = H - H.mean(1, keepdims=True)
        v3 = np.linalg.eigh(H @ H.T)[1][:, -1]
    hz = np.array([[0, 1, 0], [0, 0, 1.0]])
    seg = {"S": segS, "S3c": segS3, "P3": segP, "P3s": segPs, "Pz": segP}
    proj = {"S": hz, "S3c": np.eye(3), "P3": v3[None, :], "P3s": v3[None, :], "Pz": np.array([[1.0, 0, 0]])}
    return dict(tS=tS, tP=tP, SP=SP, seg=seg, proj=proj, v3=v3,
                inc_deg=float(np.degrees(np.arccos(min(1.0, abs(v3[0]))))))


def meas(X, ref, gm, keys, curves=False):
    out = {}
    for k in keys:
        seg = gm["seg"][k]
        if len(seg) == 0:
            out[k] = (np.full(len(X), np.nan), np.full(len(X), np.nan), np.ones(len(X), bool)) + ((None,) if curves else ())
            continue
        cc = cc_curves(X, ref, seg, gm["proj"][k])
        out[k] = peak(cc) + ((cc,) if curves else ())
    return out


# ------------------------------------------------------------------------------------------ labels (main process)
def build_labels(X, Ad, U, rng_seed=SEED):
    """X rows (row_id, fam, date, t, series) -> frame of per-row labels for EP, YR, Y2, PH, PH5, PHE + PHp code matrix"""
    X = X.copy()
    A_map = {(int(f), d): e for f, d, e in zip(Ad.fam, Ad.date, Ad.episode)}
    du = Ad.groupby("date").episode.agg(lambda s: s.unique().tolist()); date_unit = {d: v[0] for d, v in du.items()}
    X["EP"] = [A_map.get((f, d)) if s == "data" else date_unit.get(d) for f, d, s in zip(X.fam, X.date, X.series)]
    X["EP"] = X.EP.where(X.EP.notna(), None)
    y = X.date.str[:4].astype(int)
    X["DH"] = np.where(y <= 2014, "E", "L")
    q4 = np.select([y <= 2006, y <= 2010, y <= 2014, y <= 2018, y <= 2022], ["Q0-2006", "Q2007-10", "Q2011-14", "Q2015-18", "Q2019-22"], "Q2023-26")
    X["Q4"] = q4
    X["YR"] = y.astype(str); X["Y2"] = [f"{2003 + 2 * ((v - 2003) // 2)}-{2004 + 2 * ((v - 2003) // 2)}" for v in y]
    Uc = U.set_index("episode")
    cls = X.EP.map(Uc.cls)
    X["cls"] = cls
    maj = U[U.cls == "major"].copy()
    for c in ("start", "end", "w0", "w1"):
        maj[c + "_t"] = pd.to_datetime(maj[c]).astype("datetime64[ns]").astype(np.int64) / 1e9
    maj = maj.set_index("episode")
    # family onsets: 10 % quantile of data-row times per (fam, major)
    D = X[(X.series == "data") & (X.cls == "major")]
    t0 = D.groupby(["fam", "EP"]).t.quantile(0.1)
    t0d = t0.to_dict()
    X["t0"] = [t0d.get((f, e), np.nan) if c == "major" else np.nan for f, e, c in zip(X.fam, X.EP, X.cls)]
    dt = (X.t - X.t0) / 86400.0
    ismaj = (X.cls == "major").to_numpy(); isb = (X.cls == "burst").to_numpy()
    ph = np.full(len(X), None, dtype=object)
    ph[isb] = "BTW"; ph[ismaj & (dt <= 1.0).to_numpy()] = "ARR"; ph[ismaj & (dt > 1.0).to_numpy()] = "AFT"
    ph[ismaj & ~np.isfinite(dt.to_numpy())] = None
    X["PH"] = ph
    # PH5
    bst = U[U.cls == "burst"].set_index("episode")
    bstart = pd.to_datetime(bst.start).astype("datetime64[ns]").astype(np.int64) / 1e9
    bend = pd.to_datetime(bst.end).astype("datetime64[ns]").astype(np.int64) / 1e9
    blab = {}
    for b in bst.index:
        before = maj.w0_t - bend[b]       # >0 if burst ends before major window start
        after = bstart[b] - maj.w1_t      # >0 if burst starts after major window end
        dpre = before[(before >= 0) & (before <= 60 * 86400)]; dpost = after[(after >= 0) & (after <= 60 * 86400)]
        if len(dpre) and (not len(dpost) or dpre.min() <= dpost.min()):
            blab[b] = "BPRE"
        elif len(dpost):
            blab[b] = "BPOST"
        else:
            blab[b] = "BMID"
    p5 = np.full(len(X), None, dtype=object)
    p5[isb] = [blab.get(e, "BMID") for e in X.EP[isb]]
    dtn = dt.to_numpy()
    p5[ismaj & (dtn <= 0)] = "EARLY"; p5[ismaj & (dtn > 0) & (dtn <= 1)] = "ARR0"
    p5[ismaj & (dtn > 1) & (dtn <= 3)] = "AFT1"; p5[ismaj & (dtn > 3)] = "AFT3"
    X["PH5"] = p5
    # PHE: thirds of [start, end + 1 day)
    pe = np.full(len(X), None, dtype=object); pe[isb] = "BTW"
    st_ = X.EP.map(maj.start_t).to_numpy(dtype=float); en_ = X.EP.map(maj.end_t).to_numpy(dtype=float) + 86400.0
    fr = (X.t.to_numpy() - st_) / (en_ - st_)
    pe[ismaj & (fr < 1 / 3)] = "E1"; pe[ismaj & (fr >= 1 / 3) & (fr < 2 / 3)] = "E2"; pe[ismaj & (fr >= 2 / 3)] = "E3"
    X["PHE"] = pe
    # PHp: episode-date permutation of PH (codes 0 BTW, 1 ARR, 2 AFT, -1 none)
    majs = sorted(set(D.EP))
    off = {(f, e): t0d[(f, e)] - maj.start_t[e] for (f, e) in t0d}
    offmed = {e: np.median([v for (f, ee), v in off.items() if ee == e]) for e in majs}
    P = np.full((len(X), NPERM), -1, dtype=np.int8)
    P[isb, :] = 0
    rng = np.random.default_rng(rng_seed)
    fam_a = X.fam.to_numpy(); ep_a = X.EP.to_numpy(); t_a = X.t.to_numpy()
    st_m = np.array([maj.start_t.get(e, np.nan) if c else np.nan for e, c in zip(ep_a, ismaj)])
    for k in range(NPERM):
        while True:
            perm = rng.permutation(len(majs))
            if not np.any(perm == np.arange(len(majs))):
                break
        pi = {e: majs[j] for e, j in zip(majs, perm)}
        fake = np.array([st_m[i] + off.get((fam_a[i], pi[ep_a[i]]), offmed[pi[ep_a[i]]]) if ismaj[i] and ep_a[i] in pi else np.nan
                         for i in range(len(X))])
        d_ = (t_a - fake) / 86400.0
        ok = ismaj & np.isfinite(d_)
        P[ok & (d_ <= 1.0), k] = 1; P[ok & (d_ > 1.0), k] = 2
    return X, P


# ------------------------------------------------------------------------------------------ one family
def run_family(fam, RF, sp_map, cfg):
    rng = np.random.default_rng(SEED + fam)
    STATIONS, PAIRS, GROUPS = cfg["stations"], cfg["pairs"], cfg["groups"]
    st_rows, pr_rows, geo_rows, pm_st, pm_pr, sj_rows = [], [], [], [], [], []
    per = {}
    for ns in RF.netsta.unique():
        if ns not in STATIONS:
            continue
        R = RF[RF.netsta == ns]
        for kind in ("data", "ctrl"):
            ser = "data" if kind == "data" else "c1"
            pos = G_POS.get((ns, kind))
            if pos is None:
                continue
            r = R[(R.series == ser) & R.row_id.isin(pos.index)].copy()
            if not len(r):
                continue
            r["ix"] = pos.reindex(r.row_id).to_numpy()
            W = G_W[(ns, kind)]
            ok = np.isfinite(W[r.ix.to_numpy(), :, 0]).all(1)
            r = r[ok]
            for ep, g in r.groupby("epoch"):
                per.setdefault((ns, ep), {})[kind] = g.reset_index(drop=True)
    for k in list(per):
        d = per[k]
        if "data" not in d or len(d["data"]) < NMIN + 6:
            per.pop(k); continue
        g = d["data"]; W = G_W[(k[0], "data")]; Wr = W[g.ix.to_numpy()].astype(np.float64); w = g.P.to_numpy()
        d.update(Gsum=np.einsum("n,ncs->cs", w, Wr), wtot=w.sum())
    pool = {}
    for ns in set(k[0] for k in per):
        ks = [k for k in per if k[0] == ns]
        Gp = sum(per[k]["Gsum"] for k in ks) / sum(per[k]["wtot"] for k in ks)
        sp = float(sp_map.get((fam, ns), np.nan)); tS = pick_S(Gp, 0.0)
        pool[ns] = (tS, pick_P(Gp[0], tS, sp), sp if np.isfinite(sp) else 6.0)
    for (ns, ep), d in per.items():
        G = d["Gsum"] / d["wtot"]; _, tP, SP = pool[ns]
        gm = geom3(G, pick_S(G, 0.0), tP, SP); d["gm"] = gm
        geo_rows.append(dict(fam=fam, netsta=ns, epoch=ep, n_rows=len(d["data"]), tS=gm["tS"], tP=gm["tP"], SP=SP,
                             SmP_pick=gm["tS"] - gm["tP"], nP=len(gm["seg"]["P3"]), inc_deg=gm["inc_deg"]))

    def leo(ns, ep, msk):
        d = per[(ns, ep)]; g = d["data"]; W = G_W[(ns, "data")]
        sel = g[msk]; Wr = W[sel.ix.to_numpy()].astype(np.float64); w = sel.P.to_numpy()
        return (d["Gsum"] - np.einsum("n,ncs->cs", w, Wr)) / (d["wtot"] - w.sum())

    refcache = {}
    # ---- station bins
    for (ns, ep), d in per.items():
        g = d["data"]; W = G_W[(ns, "data")]; gm = d["gm"]
        for bk in SCHEMES:
            for lab, ge in g[g[bk].notna()].groupby(bk):
                if len(ge) < NMIN or len(ge) > len(g) - 10:
                    continue
                ref = leo(ns, ep, (g[bk] == lab).to_numpy()); refcache[(ns, ep, bk, lab)] = ref
                Wr = W[ge.ix.to_numpy()].astype(np.float64); w = ge.P.to_numpy()
                St = np.einsum("n,ncs->cs", w, Wr) / w.sum()
                Bs, _ = boot(Wr, w, rng)
                Xs = np.concatenate([St[None], fshift(St, 0.005)[None], Bs])
                mm = meas(Xs, ref, gm, ST_KEYS)
                r = dict(fam=fam, netsta=ns, epoch=ep, bin=bk, lab=str(lab), kind="data", n=len(ge),
                         date=str(pd.to_datetime(ge.date).min().date()), tmed=float(ge.t.median()))
                for k in ST_KEYS:
                    l, c, rl = mm[k]
                    r[f"lag{k}"], r[f"cc{k}"], r[f"rail{k}"] = l[0] * 1e3, c[0], bool(rl[0])
                    r[f"sig{k}"] = np.nanstd(l[2:], ddof=1) * 1e3
                for a_, b_ in (("S", "P3"), ("S", "P3s"), ("S3c", "P3"), ("S", "Pz")):
                    r[f"sig{a_}m{b_}"] = np.nanstd(mm[a_][0][2:] - mm[b_][0][2:], ddof=1) * 1e3
                r["inj5S"] = (mm["S"][0][1] - mm["S"][0][0]) * 1e3; r["inj5P3"] = (mm["P3"][0][1] - mm["P3"][0][0]) * 1e3
                st_rows.append(r)
                if "ctrl" in d:
                    c = d["ctrl"]; mc = np.flatnonzero((c[bk] == lab).to_numpy())
                    if len(mc) > len(ge):
                        mc = np.sort(rng.choice(mc, len(ge), replace=False))
                    if len(mc) >= NMIN:
                        Wc = G_W[(ns, "ctrl")][c.ix.to_numpy()[mc]].astype(np.float64)
                        Sc = Wc.mean(0); Bc, _ = boot(Wc, np.ones(len(mc)), rng)
                        mc_ = meas(np.concatenate([Sc[None], Bc]), ref, gm, ST_KEYS)
                        rc = dict(fam=fam, netsta=ns, epoch=ep, bin=bk, lab=str(lab), kind="ctrl", n=len(mc), date=r["date"],
                                  tmed=r["tmed"])
                        for k in ST_KEYS:
                            l, cc_, rl = mc_[k]
                            rc[f"lag{k}"], rc[f"cc{k}"], rc[f"rail{k}"] = l[0] * 1e3, cc_[0], bool(rl[0])
                            rc[f"sig{k}"] = np.nanstd(l[1:], ddof=1) * 1e3
                        for a_, b_ in (("S", "P3"), ("S", "P3s"), ("S3c", "P3"), ("S", "Pz")):
                            rc[f"sig{a_}m{b_}"] = np.nanstd(mc_[a_][0][1:] - mc_[b_][0][1:], ddof=1) * 1e3
                        st_rows.append(rc)
        # ---- PHp: date-permutation copies (S, P3; no bootstrap)
        if not DO_PERM:
            continue
        Pm = G_PERM[g.row_id.to_numpy()]                    # (n, NPERM) codes
        Wr = W[g.ix.to_numpy()].astype(np.float64).reshape(len(g), -1); w = g.P.to_numpy()
        for code in (1, 2):
            Mk = (Pm == code)
            cnt = Mk.sum(0); ws = (w[:, None] * Mk)
            sums_all = (ws.T @ Wr).reshape(NPERM, 3, S.NW); wsum_all = ws.sum(0)
            with np.errstate(invalid="ignore", divide="ignore"):
                d.setdefault("permref", {})[code] = ((d["Gsum"][None] - sums_all) /
                                                     (d["wtot"] - wsum_all)[:, None, None]).astype(np.float32)
            okk = np.flatnonzero((cnt >= NMIN) & (cnt <= len(g) - 10))
            if not len(okk):
                continue
            sums = sums_all[okk]; wsum = wsum_all[okk]
            X_ = sums / wsum[:, None, None]; R_ = d["permref"][code][okk].astype(np.float64)
            mm = meas(X_, R_, gm, PR_KEYS)
            for j, kk in enumerate(okk):
                pm_st.append((fam, ns, ep, int(kk), code, int(cnt[kk]), mm["S"][0][j] * 1e3, mm["P3"][0][j] * 1e3,
                              bool(mm["S"][2][j]), bool(mm["P3"][2][j])))

    # ---- pairs on common rows
    ndat = {}
    for A_, B_ in PAIRS:
        ka = [(k, v) for k, v in per.items() if k[0] == A_]
        kb = [(k, v) for k, v in per.items() if k[0] == B_]
        for kind in ("data", "ctrl"):
            for (_, epA), dA in ka:
                for (_, epB), dB in kb:
                    if kind not in dA or kind not in dB:
                        continue
                    j = dA[kind].merge(dB[kind][["row_id", "ix"]], on="row_id", suffixes=("A", "B"))
                    if len(j) < NMIN:
                        continue
                    for bk in SCHEMES:
                        if kind == "ctrl" and bk not in CTRL_PAIR_SCHEMES:
                            continue
                        for lab, je in j[j[bk].notna()].groupby(bk):
                            key = (A_, B_, epA, epB, bk, lab)
                            if kind == "data" and len(je) < NMIN:
                                continue
                            if kind == "ctrl":
                                if key not in ndat:
                                    continue
                                if len(je) > ndat[key]:
                                    je = je.iloc[np.sort(rng.choice(len(je), ndat[key], replace=False))]
                                if len(je) < NMIN:
                                    continue
                            rA = refcache.get((A_, epA, bk, lab)); rB = refcache.get((B_, epB, bk, lab))
                            if rA is None or rB is None:
                                continue
                            if kind == "data":
                                ndat[key] = len(je)
                            w = je.P.to_numpy() if kind == "data" else np.ones(len(je))
                            WA = G_W[(A_, kind)][je.ixA.to_numpy()].astype(np.float64)
                            WB = G_W[(B_, kind)][je.ixB.to_numpy()].astype(np.float64)
                            SA = np.einsum("n,ncs->cs", w, WA) / w.sum(); SB = np.einsum("n,ncs->cs", w, WB) / w.sum()
                            BA, cnt = boot(WA, w, rng); BB, _ = boot(WB, w, rng, cnt)
                            mA = meas(np.concatenate([SA[None], BA]), rA, dA["gm"], PR_KEYS)
                            mB = meas(np.concatenate([SB[None], BB]), rB, dB["gm"], PR_KEYS)
                            r = dict(fam=fam, pair=f"{A_}-{B_}", A=A_, B=B_, epA=epA, epB=epB, bin=bk, lab=str(lab),
                                     kind=kind, n=len(je), date=str(pd.to_datetime(je.date).min().date()),
                                     tmed=float(je.t.median()))
                            for k in PR_KEYS:
                                lA, cA, raA = mA[k]; lB, cB, raB = mB[k]
                                r[f"lag{k}_A"], r[f"lag{k}_B"] = lA[0] * 1e3, lB[0] * 1e3
                                r[f"cc{k}_A"], r[f"cc{k}_B"] = cA[0], cB[0]
                                r[f"rail{k}"] = bool(raA[0] or raB[0])
                                r[f"d{k}"] = (lB[0] - lA[0]) * 1e3
                                r[f"sig_d{k}"] = np.nanstd(lB[1:] - lA[1:], ddof=1) * 1e3
                            pr_rows.append(r)
                    # PHp for pairs (data only)
                    if kind == "data" and DO_PERM:
                        Pm = G_PERM[j.row_id.to_numpy()]
                        WA = G_W[(A_, "data")][j.ixA.to_numpy()].astype(np.float64).reshape(len(j), -1)
                        WB = G_W[(B_, "data")][j.ixB.to_numpy()].astype(np.float64).reshape(len(j), -1)
                        w = j.P.to_numpy()
                        # station-level PHp references (leave-bin-out on the station's own series)
                        for code in (1, 2):
                            Mk = (Pm == code); cnt = Mk.sum(0)
                            okk = np.flatnonzero(cnt >= NMIN)
                            if not len(okk):
                                continue
                            ws = (w[:, None] * Mk)[:, okk]; wsum = ws.sum(0)
                            XA = (ws.T @ WA).reshape(len(okk), 3, S.NW) / wsum[:, None, None]
                            XB = (ws.T @ WB).reshape(len(okk), 3, S.NW) / wsum[:, None, None]
                            if code not in dA.get("permref", {}) or code not in dB.get("permref", {}):
                                continue
                            RA = dA["permref"][code][okk].astype(np.float64); RB = dB["permref"][code][okk].astype(np.float64)
                            fin = np.isfinite(RA[:, 0, 0]) & np.isfinite(RB[:, 0, 0])
                            if not fin.all():
                                okk, XA, XB, RA, RB = okk[fin], XA[fin], XB[fin], RA[fin], RB[fin]
                                if not len(okk):
                                    continue
                            mA = meas(XA, RA, dA["gm"], PR_KEYS); mB = meas(XB, RB, dB["gm"], PR_KEYS)
                            for jj, kk in enumerate(okk):
                                pm_pr.append((fam, f"{A_}-{B_}", epA, epB, int(kk), code, int(cnt[kk]),
                                              (mB["S"][0][jj] - mA["S"][0][jj]) * 1e3, (mB["P3"][0][jj] - mA["P3"][0][jj]) * 1e3,
                                              bool(mA["S"][2][jj] or mB["S"][2][jj]), bool(mA["P3"][2][jj] or mB["P3"][2][jj])))

    # ---- site-joint S / P3 (task 4)
    for gname, mem in GROUPS.items():
        nodes = [(k, v) for k, v in per.items() if k[0] in mem]
        if len(nodes) < 2:
            continue
        for bk in SCHEMES:
            labs = set()
            for _, v in nodes:
                labs |= set(v["data"][bk].dropna().unique())
            for lab in sorted(labs):
                sets = []
                for (ns, ep), v in nodes:
                    if (ns, ep, bk, lab) not in refcache:
                        continue
                    g = v["data"]; sets.append(((ns, ep), set(g.row_id[g[bk] == lab])))
                while len(sets) >= 2:
                    com = set.intersection(*[s for _, s in sets])
                    if len(com) >= NMIN:
                        break
                    sets.sort(key=lambda t: len(t[1])); sets = sets[1:]
                if len(sets) < 2:
                    continue
                com = np.array(sorted(com))
                for kind in ("data", "ctrl"):
                    if kind == "ctrl":
                        if bk not in ("YR", "EP", "PH", "DH", "Q4"):
                            continue
                        csets = []
                        for (ns, ep), _ in sets:
                            c = per[(ns, ep)].get("ctrl")
                            if c is None:
                                csets = []; break
                            csets.append(set(c.row_id[c[bk] == lab]))
                        if not csets:
                            continue
                        cc_ = np.array(sorted(set.intersection(*csets)))
                        if len(cc_) > len(com):
                            cc_ = np.sort(rng.choice(cc_, len(com), replace=False))
                        if len(cc_) < NMIN:
                            continue
                        rows_ = cc_
                    else:
                        rows_ = com
                    Xs, refs, gms, wts = [], [], [], None
                    cnt = None
                    curves = {"S": [], "P3": []}; ind = []
                    for (ns, ep), _ in sets:
                        v = per[(ns, ep)]; tab = v[kind].set_index("row_id").loc[rows_]
                        Wm = G_W[(ns, kind)][tab.ix.to_numpy()].astype(np.float64)
                        w = tab.P.to_numpy() if kind == "data" else np.ones(len(tab))
                        St = np.einsum("n,ncs->cs", w, Wm) / w.sum()
                        Bm, cnt = boot(Wm, w, rng, cnt)
                        mm = meas(np.concatenate([St[None], Bm]), refcache[(ns, ep, bk, lab)], v["gm"], ("S", "P3"), curves=True)
                        if mm["P3"][3] is None:
                            continue
                        curves["S"].append(mm["S"][3]); curves["P3"].append(mm["P3"][3])
                        ind.append(((ns, ep), mm))
                    if len(ind) < 2:
                        continue
                    jS = peak(np.mean(curves["S"], 0)); jP = peak(np.mean(curves["P3"], 0))
                    r = dict(fam=fam, group=gname, bin=bk, lab=str(lab), kind=kind, n=len(rows_), n_mem=len(ind),
                             members=";".join(f"{a}|{b}" for (a, b), _ in ind),
                             lagS=jS[0][0] * 1e3, ccS=jS[1][0], railS=bool(jS[2][0]), sigS=np.nanstd(jS[0][1:], ddof=1) * 1e3,
                             lagP3=jP[0][0] * 1e3, ccP3=jP[1][0], railP3=bool(jP[2][0]),
                             sigP3=np.nanstd(jP[0][1:], ddof=1) * 1e3,
                             sigSmP3=np.nanstd(jS[0][1:] - jP[0][1:], ddof=1) * 1e3)
                    best = None
                    for (ns, ep), mm in ind:
                        sp_ = np.nanstd(mm["P3"][0][1:], ddof=1) * 1e3
                        r[f"m_{ns}_lagSmP3"] = (mm["S"][0][0] - mm["P3"][0][0]) * 1e3
                        r[f"m_{ns}_sigP3"] = sp_; r[f"m_{ns}_ccP3"] = mm["P3"][1][0]
                        r[f"m_{ns}_sigSmP3"] = np.nanstd(mm["S"][0][1:] - mm["P3"][0][1:], ddof=1) * 1e3
                        if best is None or sp_ < best[0]:
                            best = (sp_, ns, mm["P3"][1][0], r[f"m_{ns}_sigSmP3"])
                    r.update(best_sigP3=best[0], best_member=best[1], best_ccP3=best[2], best_sigSmP3=best[3])
                    sj_rows.append(r)
    return (pd.DataFrame(st_rows), pd.DataFrame(pr_rows), pd.DataFrame(geo_rows), pd.DataFrame(pm_st, columns=PM_ST_COLS),
            pd.DataFrame(pm_pr, columns=PM_PR_COLS), pd.DataFrame(sj_rows))


def _run(args):
    try:
        return args[0], run_family(*args), None
    except Exception:
        import traceback
        return args[0], tuple(pd.DataFrame() for _ in range(6)), traceback.format_exc()[-1500:]


PM_ST_COLS = ["fam", "netsta", "epoch", "draw", "code", "n", "lagS", "lagP3", "railS", "railP3"]
PM_PR_COLS = ["fam", "pair", "epA", "epB", "draw", "code", "n", "dS", "dP3", "railS", "railP3"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--region", default="SVI"); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--fams", default=""); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    reg = a.region.upper(); cfg = region_cfg(reg)
    OUT = f"{ROOT3}/{reg.lower()}"; TB, LOGD = f"{OUT}/tables", f"{ROOT3}/logs"
    os.makedirs(TB, exist_ok=True); os.makedirs(LOGD, exist_ok=True)
    log = S.Log(f"{LOGD}/measure_{reg.lower()}{a.tag}.log")
    X = pd.read_parquet(f"{cfg['tb']}/rows.parquet", columns=["row_id", "fam", "P", "series", "date", "t"])
    X["date"] = X.date.astype(str); X["series"] = X.series.astype(str)
    Ad = pd.read_csv(f"{cfg['tb']}/active_days.csv"); U = pd.read_csv(f"{cfg['tb']}/episode_units.csv")
    XL, Pmat = build_labels(X, Ad, U)
    global G_PERM
    G_PERM_arr = np.full((int(XL.row_id.max()) + 1, NPERM), -1, dtype=np.int8); G_PERM_arr[XL.row_id.to_numpy()] = Pmat
    G_PERM = G_PERM_arr
    XL[["row_id", "fam", "series", "date", "t", "EP", "cls", "t0", "PH", "PH5", "PHE"]].to_parquet(f"{TB}/row_labels{a.tag}.parquet", index=False)
    log(f"{reg}: labels built; PH counts (data) {XL[XL.series == 'data'].PH.value_counts(dropna=False).to_dict()}")
    RS = pd.read_parquet(f"{cfg['tb']}/row_station.parquet", columns=["row_id", "netsta", "epoch", "date"])
    RS["netsta"] = RS.netsta.astype(str); RS["epoch"] = RS.epoch.astype(str); RS = RS.drop(columns="date")
    RS = RS[RS.netsta.isin(cfg["stations"])].merge(XL[["row_id", "fam", "P", "series", "date", "t", "EP", "YR", "Y2", "PH", "PH5", "PHE", "DH", "Q4"]], on="row_id")
    O = pd.read_csv(f"{cfg['tb']}/offsets_sp.csv"); sp_map = {(int(f), n): s for f, n, s in zip(O.fam, O.netsta, O.sp)}
    fams = [int(x) for x in a.fams.split(",")] if a.fams else sorted(pd.read_csv(f"{cfg['tb']}/families.csv").fam.astype(int))
    RS = RS[RS.fam.isin(fams)]
    t0 = time.time(); tot = 0
    for ns in cfg["stations"]:
        for kind, sfx in (("data", ""), ("ctrl", "__c2")):
            p = f"{cfg['win']}/{ns}{sfx}.npz"
            if not os.path.exists(p):
                continue
            z = np.load(p)
            rid = z["row_id"]; keep = np.isin(rid, RS.row_id.to_numpy())
            if not keep.any():
                continue
            G_W[(ns, kind)] = z["W"][keep]; G_POS[(ns, kind)] = pd.Series(np.arange(keep.sum()), index=rid[keep])
            tot += G_W[(ns, kind)].nbytes
    log(f"loaded windows: {len(G_W)} arrays, {tot / 1e9:.2f} GB, {time.time() - t0:.0f}s")
    import multiprocessing as mp
    res = {k: [] for k in ("st", "pr", "geo", "pms", "pmp", "sj")}
    args = [(f, RS[RS.fam == f], sp_map, cfg) for f in fams]
    args.sort(key=lambda t: -len(t[1]))
    with mp.get_context("fork").Pool(a.workers, maxtasksperchild=4) as p:
        for f, out, err in p.imap_unordered(_run, args):
            if err:
                log(f"ERROR fam {f}: {err}")
            for k, v in zip(res, out):
                if len(v):
                    res[k].append(v)
            log(f"  fam {f}: station {len(out[0])}, pair {len(out[1])}, perm {len(out[3])}/{len(out[4])}, joint {len(out[5])}  {time.time() - t0:.0f}s")
    cat = {k: (pd.concat(v, ignore_index=True) if v else pd.DataFrame()) for k, v in res.items()}
    cat["st"].to_parquet(f"{TB}/lags_station{a.tag}.parquet", index=False)
    cat["pr"].to_parquet(f"{TB}/lags_pair{a.tag}.parquet", index=False)
    cat["geo"].to_csv(f"{TB}/series_geom{a.tag}.csv", index=False)
    cat["pms"].to_parquet(f"{TB}/perm_station{a.tag}.parquet", index=False)
    cat["pmp"].to_parquet(f"{TB}/perm_pair{a.tag}.parquet", index=False)
    cat["sj"].to_parquet(f"{TB}/site_joint{a.tag}.parquet", index=False)
    log(f"done {time.time() - t0:.0f}s: " + ", ".join(f"{k} {len(v)}" for k, v in cat.items()))
    open(f"{ROOT3}/MEASURE_{reg}{a.tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
