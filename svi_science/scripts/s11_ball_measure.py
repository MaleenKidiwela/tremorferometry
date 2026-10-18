#!/usr/bin/env python
"""s11_ball_measure -- S11 ballistic test (svi_science/results/s11_depth/ballistic/PLAN.md §1-3, §7).
Per (family, station, epoch, episode): P-weighted episode stacks (P >= 0.5 rows, n >= 20), UNALIGNED windows, S lag
(N+E joint cc, [tS-0.5, tS+1.5]) and P lag (Z + polarised horizontal, full P wavetrain) against the station-epoch
leave-episode-out grand stack, +-100 ms search, parabolic refinement, 200 bootstrap draws. Inter-station dS / dP on rows
present at both stations (common bootstrap draws). Scrambled control (c2 windows) through the same path.

    python s11_ball_measure.py --workers 8 [--fams 10,87] [--tag _test]
"""
from __future__ import annotations
import os, sys, argparse, time, zlib
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_multi_common as M
from s10_measure import pick_S, env

OUT = f"{S.ROOT}/svi_science/results/s11_depth/ballistic"
TB, LOGD = f"{OUT}/tables", f"{OUT}/logs"
STATIONS = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011", "PB.B926", "CN.VGZ", "CN.SNB"]
PAIRS = [("CN.PGC", "PB.B926"), ("CN.PGC", "CN.SNB"), ("CN.PGC", "CN.VGZ"), ("CN.PGC", "PB.B009"),
         ("CN.PGC", "PB.B010"), ("CN.PGC", "PB.B011"), ("PB.B011", "PB.B009"), ("PB.B011", "PB.B010"),
         ("PB.B011", "PB.B926"), ("PB.B011", "CN.SNB"), ("PB.B011", "CN.VGZ")]          # B011 refs: Amendment 1
NMIN, NBOOT, LMAX = 20, 200, 4            # rows, draws, samples (4 / 40 Hz = 100 ms)
SEED = 20261012
INJ = 0.005                                # self-check shift (s)
FREQ = np.fft.rfftfreq(S.NW, 1 / S.FS)

G_W, G_POS = {}, {}                        # globals filled before fork: (netsta, kind) -> W, row_id -> index


# ------------------------------------------------------------------------------------------ lag core
def lag_batch(X, g, seg, proj):
    """X (m, 3, NW) stacks, g (3, NW) reference, seg sample indices, proj (k, 3) -> lag_s (m,), cc (m,), railed (m,)"""
    lags = np.arange(-LMAX, LMAX + 1)
    xs = np.einsum("kc,mcn->mkn", proj, X); gs = proj @ g
    gg = gs[:, seg]; gg = gg - gg.mean(1, keepdims=True); dg = (gg * gg).sum()
    m = X.shape[0]; cc = np.zeros((m, len(lags)))
    for j, l in enumerate(lags):
        x = xs[:, :, seg + l]; x = x - x.mean(2, keepdims=True)
        cc[:, j] = (x * gg[None]).sum((1, 2)) / np.sqrt(np.maximum((x * x).sum((1, 2)) * dg, 1e-30))
    i = np.argmax(cc, 1); r = np.arange(m); ii = np.clip(i, 1, len(lags) - 2)
    y0, y1, y2 = cc[r, ii - 1], cc[r, ii], cc[r, ii + 1]; d2 = y0 - 2 * y1 + y2
    delta = np.where((d2 < 0) & (i == ii), 0.5 * (y0 - y2) / np.where(d2 == 0, -1, d2), 0.0)
    rail = (i == 0) | (i == len(lags) - 1)
    return (lags[i] + np.clip(delta, -1, 1)) / S.FS, cc[r, i], rail


def fshift(x, dt):
    """x(tau - dt): content delayed by dt"""
    return np.fft.irfft(np.fft.rfft(x, axis=-1) * np.exp(-2j * np.pi * FREQ * dt), n=S.NW, axis=-1)


def boot(Wr, w, rng, cnt=None):
    """Wr (n, 3, NW), weights w -> (NBOOT, 3, NW) bootstrap weighted stacks; cnt reusable"""
    n = Wr.shape[0]
    if cnt is None:
        cnt = np.stack([np.bincount(rng.integers(0, n, n), minlength=n) for _ in range(NBOOT)]).astype(np.float64)
    ww = cnt * w[None]
    B = (ww @ Wr.reshape(n, -1)) / ww.sum(1, keepdims=True)
    return B.reshape(NBOOT, 3, S.NW), cnt


# ------------------------------------------------------------------------------------------ series setup
def pick_P(Gz_pool, tS, sp):
    """Amendment 1: Z-envelope max of the all-epoch pooled stack in [tS - 1.3 SP, tS - 1.0 SP]"""
    SP = sp if np.isfinite(sp) else 6.0
    e = env(Gz_pool)
    m = np.flatnonzero((S.TAU >= tS - 1.3 * SP) & (S.TAU <= tS - 1.0 * SP) & (S.TAU >= S.T0W + 1.0))
    return float(S.TAU[m[np.argmax(e[m])]]) if len(m) else np.nan


def series_geom(G, tS, tP, SP):
    segS = np.flatnonzero((S.TAU >= tS - 0.5) & (S.TAU <= tS + 1.5))
    segP = np.zeros(0, int)
    if np.isfinite(tP):
        segP = np.flatnonzero((S.TAU >= max(tP - 0.5, S.T0W + 0.5)) & (S.TAU <= min(tP + 2.5, tS - 1.0)))
    if len(segP) < 20:
        segP = np.zeros(0, int)
    if len(segP):
        H = G[1:, segP]; H = H - H.mean(1, keepdims=True); C = H @ H.T
        v = np.linalg.eigh(C)[1][:, -1]
    else:
        v = np.array([1.0, 0.0])
    return dict(tS=tS, tP=tP, SP=SP, segS=segS, segP=segP, r1=float(v[0]), r2=float(v[1]),
                projS=np.array([[0, 1, 0], [0, 0, 1.0]]), projP=np.array([[1, 0, 0], [0, v[0], v[1]]]),
                projPz=np.array([[1.0, 0, 0]]))


def measure_one(X, ref, gm, keys=("S", "P", "Pz")):
    """X (m, 3, NW) -> dict of arrays per key"""
    out = {}
    for k in keys:
        seg = gm["segS"] if k == "S" else gm["segP"]
        if len(seg) == 0:
            out[k] = (np.full(len(X), np.nan),) * 2 + (np.ones(len(X), bool),); continue
        out[k] = lag_batch(X, ref, seg, gm[{"S": "projS", "P": "projP", "Pz": "projPz"}[k]])
    return out


# ------------------------------------------------------------------------------------------ one family
def run_family(fam, RF, sp_map, A_map, date_unit, cls_map):
    """RF: row_station rows of this family (row_id, netsta, epoch, date, series, P)"""
    rng = np.random.default_rng(SEED + fam)
    st_rows, pr_rows, geo_rows = [], [], []
    per = {}                       # (ns, ep) -> dict(Gsum, wtot, gm, rows DataFrame (data), ctrl DataFrame)
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
            if kind == "data":
                r["episode"] = [A_map.get((fam, d), None) for d in r.date]
            else:
                r["episode"] = r.date.map(date_unit)
            for ep, g in r.groupby("epoch"):
                per.setdefault((ns, ep), {})[kind] = g.reset_index(drop=True)
    for (ns, ep), dct in list(per.items()):
        if "data" not in dct or len(dct["data"]) < NMIN + 6:
            per.pop((ns, ep)); continue
        g = dct["data"]; W = G_W[(ns, "data")]; Wr = W[g.ix.to_numpy()].astype(np.float64); w = g.P.to_numpy()
        dct.update(Gsum=np.einsum("n,ncs->cs", w, Wr), wtot=w.sum())
    pool = {}
    for ns in STATIONS:
        ks = [k for k in per if k[0] == ns]
        if not ks:
            continue
        Gp = sum(per[k]["Gsum"] for k in ks) / sum(per[k]["wtot"] for k in ks)
        sp = float(sp_map.get((fam, ns), np.nan)); tS = pick_S(Gp, 0.0)
        pool[ns] = (tS, pick_P(Gp[0], tS, sp), sp if np.isfinite(sp) else 6.0)
    for (ns, ep), dct in per.items():
        Gsum, wtot = dct["Gsum"], dct["wtot"]; tS_p, tP, SP = pool[ns]
        gm = series_geom(Gsum / wtot, pick_S(Gsum / wtot, 0.0), tP, SP)
        dct["gm"] = gm
        geo_rows.append(dict(fam=fam, netsta=ns, epoch=ep, n_rows=len(dct["data"]), tS=gm["tS"], tP=gm["tP"], SP=gm["SP"],
                             SmP_pick=gm["tS"] - gm["tP"], nP=len(gm["segP"]), r1=gm["r1"], r2=gm["r2"]))

    def leo(ns, ep, ids_mask):
        d = per[(ns, ep)]; g = d["data"]; W = G_W[(ns, "data")]
        sel = g[ids_mask]; Wr = W[sel.ix.to_numpy()].astype(np.float64); w = sel.P.to_numpy()
        return (d["Gsum"] - np.einsum("n,ncs->cs", w, Wr)) / (d["wtot"] - w.sum())

    # ---- per-station episode stacks
    for (ns, ep), d in per.items():
        g = d["data"]; W = G_W[(ns, "data")]; gm = d["gm"]
        for e_, ge in g[g.episode.notna()].groupby("episode"):
            if len(ge) < NMIN:
                continue
            msk = g.episode == e_
            ref = leo(ns, ep, msk)
            Wr = W[ge.ix.to_numpy()].astype(np.float64); w = ge.P.to_numpy()
            St = np.einsum("n,ncs->cs", w, Wr) / w.sum()
            Bs, _ = boot(Wr, w, rng)
            Xall = np.concatenate([St[None], fshift(St, INJ)[None], Bs])
            mm = measure_one(Xall, ref, gm)
            r = dict(fam=fam, netsta=ns, epoch=ep, episode=e_, cls=cls_map.get(e_, ""), kind="data", n=len(ge),
                     date=str(pd.to_datetime(ge.date).min().date()))
            for k in ("S", "P", "Pz"):
                l, c, rl = mm[k]
                r[f"lag{k}"], r[f"cc{k}"], r[f"rail{k}"] = l[0] * 1e3, c[0], bool(rl[0])
                r[f"sig{k}"] = np.nanstd(l[2:], ddof=1) * 1e3
                r[f"bootrail{k}"] = float(np.mean(rl[2:]))
            r["inj5S"] = (mm["S"][0][1] - mm["S"][0][0]) * 1e3
            r["sigSmP"] = np.nanstd((mm["S"][0][2:] - mm["P"][0][2:]), ddof=1) * 1e3
            st_rows.append(r)
            # scrambled control: same family/station/epoch, rows on the episode's dates, n <= n data, equal weights
            if "ctrl" in d:
                c = d["ctrl"]; mc = np.flatnonzero((c.episode == e_).to_numpy())
                if len(mc) > len(ge):
                    mc = np.sort(rng.choice(mc, len(ge), replace=False))
                if len(mc) >= NMIN:
                    Wc = G_W[(ns, "ctrl")][c.ix.to_numpy()[mc]].astype(np.float64); wc = np.ones(len(mc))
                    Sc = Wc.mean(0); Bc, _ = boot(Wc, wc, rng)
                    mc_ = measure_one(np.concatenate([Sc[None], Bc]), ref, gm)
                    rc = dict(fam=fam, netsta=ns, epoch=ep, episode=e_, cls=cls_map.get(e_, ""), kind="ctrl", n=len(mc),
                              date=r["date"])
                    for k in ("S", "P", "Pz"):
                        l, cc_, rl = mc_[k]
                        rc[f"lag{k}"], rc[f"cc{k}"], rc[f"rail{k}"] = l[0] * 1e3, cc_[0], bool(rl[0])
                        rc[f"sig{k}"] = np.nanstd(l[1:], ddof=1) * 1e3; rc[f"bootrail{k}"] = float(np.mean(rl[1:]))
                    rc["sigSmP"] = np.nanstd((mc_["S"][0][1:] - mc_["P"][0][1:]), ddof=1) * 1e3
                    st_rows.append(rc)

    # ---- pairs on common rows
    for A_, B_ in PAIRS:
        for kind in ("data", "ctrl"):
            ka = [(k, v) for k, v in per.items() if k[0] == A_ and kind in v]
            kb = [(k, v) for k, v in per.items() if k[0] == B_ and kind in v]
            for (_, epA), dA in ka:
                for (_, epB), dB in kb:
                    a = dA[kind]; b = dB[kind]
                    j = a.merge(b[["row_id", "ix"]], on="row_id", suffixes=("A", "B"))
                    j = j[j.episode.notna()]
                    for e_, je in j.groupby("episode"):
                        if kind == "data" and len(je) < NMIN:
                            continue
                        if kind == "ctrl":
                            ndat = [x for x in pr_rows if x["kind"] == "data" and x["pair"] == f"{A_}-{B_}" and
                                    x["epA"] == epA and x["epB"] == epB and x["episode"] == e_]
                            if not ndat:
                                continue
                            if len(je) > ndat[0]["n"]:
                                je = je.iloc[np.sort(rng.choice(len(je), ndat[0]["n"], replace=False))]
                            if len(je) < NMIN:
                                continue
                        refA = leo(A_, epA, dA["data"].episode == e_)
                        refB = leo(B_, epB, dB["data"].episode == e_)
                        w = je.P.to_numpy() if kind == "data" else np.ones(len(je))
                        WA = G_W[(A_, kind)][je.ixA.to_numpy()].astype(np.float64)
                        WB = G_W[(B_, kind)][je.ixB.to_numpy()].astype(np.float64)
                        SA = np.einsum("n,ncs->cs", w, WA) / w.sum(); SB = np.einsum("n,ncs->cs", w, WB) / w.sum()
                        BA, cnt = boot(WA, w, rng); BB, _ = boot(WB, w, rng, cnt)
                        mA = measure_one(np.concatenate([SA[None], BA]), refA, dA["gm"], ("S", "P"))
                        mB = measure_one(np.concatenate([SB[None], BB]), refB, dB["gm"], ("S", "P"))
                        r = dict(fam=fam, pair=f"{A_}-{B_}", A=A_, B=B_, epA=epA, epB=epB, episode=e_,
                                 cls=cls_map.get(e_, ""), kind=kind, n=len(je), date=str(pd.to_datetime(je.date).min().date()))
                        for k in ("S", "P"):
                            lA, cA, rA = mA[k]; lB, cB, rB = mB[k]
                            r[f"lag{k}_A"], r[f"lag{k}_B"] = lA[0] * 1e3, lB[0] * 1e3
                            r[f"cc{k}_A"], r[f"cc{k}_B"] = cA[0], cB[0]
                            r[f"rail{k}"] = bool(rA[0] or rB[0])
                            r[f"d{k}"] = (lB[0] - lA[0]) * 1e3
                            r[f"sig_d{k}"] = np.nanstd(lB[1:] - lA[1:], ddof=1) * 1e3
                            r[f"bootrail_d{k}"] = float(np.mean(rA[1:] | rB[1:]))
                        pr_rows.append(r)
    return st_rows, pr_rows, geo_rows


def _run(args):
    try:
        return args[0], run_family(*args), None
    except Exception:
        import traceback
        return args[0], ([], [], []), traceback.format_exc()[-800:]


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
    assert (du.str.len() == 1).all(), "date maps to >1 episode unit"
    date_unit = {d: v[0] for d, v in du.items()}
    U = pd.read_csv(f"{M.TB}/episode_units.csv"); cls_map = dict(zip(U.episode, U.cls))
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
    args = [(f, RS[RS.fam == f], sp_map, A_map, date_unit, cls_map) for f in fams]
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
