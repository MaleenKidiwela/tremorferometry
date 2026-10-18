#!/usr/bin/env python
"""s11_ball2_measure -- S11 ballistic round 2, measurement pass (ballistic/round2/PLAN.md §0).
Round-1 machinery (s11_ball_measure: lag_batch, boot, fshift, pick_P) reused unchanged. New: bins EP / YR / Y2, P estimators
P3 (3C principal direction, headline), P3j (joint 3C), Pz, Pzr (round 1), split-half diagnostics, 18 common-row pairs
(6 quartet + 12 far-vs-quartet), scrambled control (c2) for station bins and pairs.

    python s11_ball2_measure.py --workers 8 [--fams 10] [--tag _test]
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
from s11_ball_measure import lag_batch, boot, fshift, pick_P

OUT = f"{S.ROOT}/svi_science/results/s11_depth/ballistic/round2"
TB, LOGD = f"{OUT}/tables", f"{OUT}/logs"
STATIONS = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011", "PB.B926", "CN.VGZ", "CN.SNB"]
QUART = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"]
FAR = ["PB.B926", "CN.SNB", "CN.VGZ"]
PAIRS = [(a, b) for i, a in enumerate(QUART) for b in QUART[i + 1:]] + [(a, x) for x in FAR for a in QUART]
NMIN, NBOOT = 20, B1.NBOOT
SEED = 20261013
PK = ("P3", "P3j", "Pz", "Pzr")
SVAR = os.environ.get("S11B2_SVAR", "")          # Amendment 4: S-estimator variants (quartet pairs, EP bins only)
if SVAR:
    PAIRS = [(a, b) for i, a in enumerate(QUART) for b in QUART[i + 1:]]
    STATIONS = QUART
BINS = ("EP",) if SVAR else ("EP", "YR", "Y2")
G_W, G_POS = B1.G_W, B1.G_POS          # shared globals (filled before fork)


def binlab(kind, date, A):
    if kind == "EP":
        return A
    y = int(str(date)[:4])
    return str(y) if kind == "YR" else f"{2003 + 2 * ((y - 2003) // 2)}-{2004 + 2 * ((y - 2003) // 2)}"


def geom2(G, tS, tP, SP):
    sw = {"s3c25": (0.5, 2.5), "sh07": (0.3, 0.7), "s3c40": (0.5, 4.0)}.get(SVAR, (0.5, 1.5))
    segS = np.flatnonzero((S.TAU >= tS - sw[0]) & (S.TAU <= tS + sw[1]))
    segP = np.zeros(0, int)
    if np.isfinite(tP):
        segP = np.flatnonzero((S.TAU >= max(tP - 0.5, S.T0W + 0.5)) & (S.TAU <= min(tP + 2.5, tS - 1.0)))
    if len(segP) < 20:
        segP = np.zeros(0, int)
    v3 = np.array([1.0, 0, 0]); v2 = np.array([1.0, 0])
    if len(segP):
        H = G[:, segP]; H = H - H.mean(1, keepdims=True)
        v3 = np.linalg.eigh(H @ H.T)[1][:, -1]
        H2 = H[1:]; v2 = np.linalg.eigh(H2 @ H2.T)[1][:, -1]
    proj = {"S": np.eye(3) if SVAR in ("s3c25", "s3c40") else np.array([[0, 1, 0], [0, 0, 1.0]]), "P3": v3[None, :], "P3j": np.eye(3),
            "Pz": np.array([[1.0, 0, 0]]), "Pzr": np.array([[1, 0, 0], [0, v2[0], v2[1]]])}
    return dict(tS=tS, tP=tP, SP=SP, segS=segS, segP=segP, proj=proj, v3=v3,
                inc_deg=float(np.degrees(np.arccos(min(1.0, abs(v3[0]))))))


def meas(X, ref, gm, keys):
    out = {}
    for k in keys:
        seg = gm["segS"] if k == "S" else gm["segP"]
        if len(seg) == 0:
            out[k] = (np.full(len(X), np.nan), np.full(len(X), np.nan), np.ones(len(X), bool)); continue
        out[k] = lag_batch(X, ref, seg, gm["proj"][k])
    return out


def zcc(a, b, seg, proj):
    xa = proj @ a[:, seg]; xb = proj @ b[:, seg]
    xa = xa - xa.mean(1, keepdims=True); xb = xb - xb.mean(1, keepdims=True)
    return float((xa * xb).sum() / np.sqrt((xa * xa).sum() * (xb * xb).sum() + 1e-30))


def run_family(fam, RF, sp_map, A_map, date_unit):
    rng = np.random.default_rng(SEED + fam)
    st_rows, pr_rows, geo_rows = [], [], []
    per = {}
    for ns in STATIONS:
        R = RF[RF.netsta == ns]
        for kind in ("data", "ctrl"):
            ser = "data" if kind == "data" else "c1"
            pos = G_POS.get((ns, kind))
            if pos is None:
                continue
            r = R[(R.series == ser) & R.row_id.isin(pos.index)].copy()
            r["ix"] = pos.reindex(r.row_id).to_numpy()
            W = G_W[(ns, kind)]
            ok = np.isfinite(W[r.ix.to_numpy(), :, 0]).all(1) if len(r) else np.zeros(0, bool)
            r = r[ok]
            r["EP"] = [A_map.get((fam, d)) for d in r.date] if kind == "data" else r.date.map(date_unit)
            r["YR"] = [binlab("YR", d, None) for d in r.date]
            r["Y2"] = [binlab("Y2", d, None) for d in r.date]
            for ep, g in r.groupby("epoch"):
                per.setdefault((ns, ep), {})[kind] = g.reset_index(drop=True)
    for k in list(per):
        d = per[k]
        if "data" not in d or len(d["data"]) < NMIN + 6:
            per.pop(k); continue
        g = d["data"]; W = G_W[(k[0], "data")]; Wr = W[g.ix.to_numpy()].astype(np.float64); w = g.P.to_numpy()
        d.update(Gsum=np.einsum("n,ncs->cs", w, Wr), wtot=w.sum())
        # split-half (odd/even rows in time order)
        o = np.argsort(g.date.to_numpy(), kind="stable"); h1, h2 = o[0::2], o[1::2]
        d["half"] = (np.einsum("n,ncs->cs", w[h1], Wr[h1]) / w[h1].sum(), np.einsum("n,ncs->cs", w[h2], Wr[h2]) / w[h2].sum())
    pool = {}
    for ns in STATIONS:
        ks = [k for k in per if k[0] == ns]
        if ks:
            Gp = sum(per[k]["Gsum"] for k in ks) / sum(per[k]["wtot"] for k in ks)
            sp = float(sp_map.get((fam, ns), np.nan)); tS = pick_S(Gp, 0.0)
            pool[ns] = (tS, pick_P(Gp[0], tS, sp), sp if np.isfinite(sp) else 6.0)
    for (ns, ep), d in per.items():
        G = d["Gsum"] / d["wtot"]; _, tP, SP = pool[ns]
        gm = geom2(G, pick_S(G, 0.0), tP, SP); d["gm"] = gm
        r = dict(fam=fam, netsta=ns, epoch=ep, n_rows=len(d["data"]), tS=gm["tS"], tP=gm["tP"], SP=SP,
                 SmP_pick=gm["tS"] - gm["tP"], nP=len(gm["segP"]), v3z=gm["v3"][0], v31=gm["v3"][1], v32=gm["v3"][2],
                 inc_deg=gm["inc_deg"])
        a, b = d["half"]
        r["half_ccS"] = zcc(a, b, gm["segS"], gm["proj"]["S"])
        for k in PK:
            r[f"half_cc{k}"] = zcc(a, b, gm["segP"], gm["proj"][k]) if len(gm["segP"]) else np.nan
        geo_rows.append(r)

    def leo(ns, ep, msk):
        d = per[(ns, ep)]; g = d["data"]; W = G_W[(ns, "data")]
        sel = g[msk]; Wr = W[sel.ix.to_numpy()].astype(np.float64); w = sel.P.to_numpy()
        return (d["Gsum"] - np.einsum("n,ncs->cs", w, Wr)) / (d["wtot"] - w.sum())

    keys_st = ("S",) + PK
    refcache = {}
    for (ns, ep), d in per.items():
        g = d["data"]; W = G_W[(ns, "data")]; gm = d["gm"]
        for bk in BINS:
            for lab, ge in g[g[bk].notna()].groupby(bk):
                if len(ge) < NMIN or len(ge) > len(g) - 10:
                    continue
                ref = leo(ns, ep, (g[bk] == lab).to_numpy()); refcache[(ns, ep, bk, lab)] = ref
                Wr = W[ge.ix.to_numpy()].astype(np.float64); w = ge.P.to_numpy()
                St = np.einsum("n,ncs->cs", w, Wr) / w.sum()
                Bs, _ = boot(Wr, w, rng)
                X = np.concatenate([St[None], fshift(St, 0.005)[None], Bs])
                mm = meas(X, ref, gm, keys_st)
                r = dict(fam=fam, netsta=ns, epoch=ep, bin=bk, lab=str(lab), kind="data", n=len(ge),
                         date=str(pd.to_datetime(ge.date).min().date()))
                for k in keys_st:
                    l, c, rl = mm[k]
                    r[f"lag{k}"], r[f"cc{k}"], r[f"rail{k}"] = l[0] * 1e3, c[0], bool(rl[0])
                    r[f"sig{k}"] = np.nanstd(l[2:], ddof=1) * 1e3
                    if k != "S":
                        r[f"sigSm{k}"] = np.nanstd(mm["S"][0][2:] - l[2:], ddof=1) * 1e3
                r["inj5S"] = (mm["S"][0][1] - mm["S"][0][0]) * 1e3
                r["inj5P3"] = (mm["P3"][0][1] - mm["P3"][0][0]) * 1e3
                st_rows.append(r)
                if "ctrl" in d:
                    c = d["ctrl"]; mc = np.flatnonzero((c[bk] == lab).to_numpy())
                    if len(mc) > len(ge):
                        mc = np.sort(rng.choice(mc, len(ge), replace=False))
                    if len(mc) >= NMIN:
                        Wc = G_W[(ns, "ctrl")][c.ix.to_numpy()[mc]].astype(np.float64)
                        Sc = Wc.mean(0); Bc, _ = boot(Wc, np.ones(len(mc)), rng)
                        mc_ = meas(np.concatenate([Sc[None], Bc]), ref, gm, keys_st)
                        rc = dict(fam=fam, netsta=ns, epoch=ep, bin=bk, lab=str(lab), kind="ctrl", n=len(mc), date=r["date"])
                        for k in keys_st:
                            l, cc_, rl = mc_[k]
                            rc[f"lag{k}"], rc[f"cc{k}"], rc[f"rail{k}"] = l[0] * 1e3, cc_[0], bool(rl[0])
                            rc[f"sig{k}"] = np.nanstd(l[1:], ddof=1) * 1e3
                        st_rows.append(rc)

    # ---- pairs on common rows
    ndat = {}
    for A_, B_ in PAIRS:
        for kind in ("data", "ctrl"):
            ka = [(k, v) for k, v in per.items() if k[0] == A_ and kind in v]
            kb = [(k, v) for k, v in per.items() if k[0] == B_ and kind in v]
            for (_, epA), dA in ka:
                for (_, epB), dB in kb:
                    j = dA[kind].merge(dB[kind][["row_id", "ix"]], on="row_id", suffixes=("A", "B"))
                    for bk in BINS:
                        if kind == "ctrl" and bk == "EP":
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
                            mA = meas(np.concatenate([SA[None], BA]), rA, dA["gm"], ("S", "P3"))
                            mB = meas(np.concatenate([SB[None], BB]), rB, dB["gm"], ("S", "P3"))
                            r = dict(fam=fam, pair=f"{A_}-{B_}", A=A_, B=B_, epA=epA, epB=epB, bin=bk, lab=str(lab),
                                     kind=kind, n=len(je), date=str(pd.to_datetime(je.date).min().date()))
                            for k in ("S", "P3"):
                                lA, cA, raA = mA[k]; lB, cB, raB = mB[k]
                                r[f"lag{k}_A"], r[f"lag{k}_B"] = lA[0] * 1e3, lB[0] * 1e3
                                r[f"cc{k}_A"], r[f"cc{k}_B"] = cA[0], cB[0]
                                r[f"rail{k}"] = bool(raA[0] or raB[0])
                                r[f"d{k}"] = (lB[0] - lA[0]) * 1e3
                                r[f"sig_d{k}"] = np.nanstd(lB[1:] - lA[1:], ddof=1) * 1e3
                            pr_rows.append(r)
    return st_rows, pr_rows, geo_rows


def _run(args):
    try:
        return args[0], run_family(*args), None
    except Exception:
        import traceback
        return args[0], ([], [], []), traceback.format_exc()[-1200:]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--fams", default=""); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    os.makedirs(TB, exist_ok=True); os.makedirs(LOGD, exist_ok=True)
    log = S.Log(f"{LOGD}/measure{a.tag}.log")
    X = pd.read_parquet(f"{M.TB}/rows.parquet", columns=["row_id", "fam", "P", "series"])
    RS = pd.read_parquet(f"{M.TB}/row_station.parquet", columns=["row_id", "netsta", "epoch", "date"])
    RS = RS[RS.netsta.isin(STATIONS)].merge(X, on="row_id")
    Ad = pd.read_csv(f"{M.TB}/active_days.csv")
    A_map = {(int(f), d): e for f, d, e in zip(Ad.fam, Ad.date, Ad.episode)}
    du = Ad.groupby("date").episode.agg(lambda s: s.unique().tolist())
    date_unit = {d: v[0] for d, v in du.items()}
    O = pd.read_csv(f"{M.TB}/offsets_sp.csv"); sp_map = {(int(f), n): s for f, n, s in zip(O.fam, O.netsta, O.sp)}
    fams = [int(x) for x in a.fams.split(",")] if a.fams else sorted(pd.read_csv(f"{M.TB}/families.csv").fam.astype(int))
    RS = RS[RS.fam.isin(fams)]
    t0 = time.time()
    for ns in STATIONS:
        for kind, sfx in (("data", ""), ("ctrl", "__c2")):
            p = f"{M.WIN}/{ns}{sfx}.npz"
            if not os.path.exists(p):
                log(f"{ns}{sfx}: no windows"); continue
            z = np.load(p)
            rid = z["row_id"]; keep = np.isin(rid, RS.row_id.to_numpy())
            G_W[(ns, kind)] = z["W"][keep]; G_POS[(ns, kind)] = pd.Series(np.arange(keep.sum()), index=rid[keep])
            log(f"loaded {ns}{sfx}: {keep.sum()} rows ({G_W[(ns, kind)].nbytes / 1e9:.2f} GB) {time.time() - t0:.0f}s")
    import multiprocessing as mp
    st, pr, geo = [], [], []
    args = [(f, RS[RS.fam == f], sp_map, A_map, date_unit) for f in fams]
    args.sort(key=lambda t: -len(t[1]))
    with mp.get_context("fork").Pool(a.workers, maxtasksperchild=4) as p:
        for f, (s_, p_, g_), err in p.imap_unordered(_run, args):
            if err:
                log(f"ERROR fam {f}: {err}")
            st += s_; pr += p_; geo += g_
            log(f"  fam {f}: station rows {len(s_)}, pair rows {len(p_)}  {time.time() - t0:.0f}s")
    pd.DataFrame(st).to_parquet(f"{TB}/lags_station{a.tag}.parquet", index=False)
    pd.DataFrame(pr).to_parquet(f"{TB}/lags_pair{a.tag}.parquet", index=False)
    pd.DataFrame(geo).to_csv(f"{TB}/series_geom{a.tag}.csv", index=False)
    log(f"done {time.time() - t0:.0f}s: station rows {len(st)}, pair rows {len(pr)}")
    open(f"{OUT}/MEASURE{a.tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
