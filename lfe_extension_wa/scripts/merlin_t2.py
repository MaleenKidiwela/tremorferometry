#!/usr/bin/env python
"""Merlin test T2 — WITHHELD-STATION template recovery (can a new station's template be built from the
permanent network's candidate times alone?).  Array era, where the answer is known.

Pairs (family f, window w, withheld array station s): XU / XG stations with a CLEAR array-truth template for f
(stack of f's BUILD-half array detections at s, >= 30 events, catalogue clear rule), up to --cap stations per
(family, window) drawn at random (fixed seed; no selection on quality beyond "clear").  Windows: V2008 (XU;
SW + NW families), V2010 (XU + XG; NW, Olympic, Arm B).  No Oregon families are in this study (XA not used).

Pretend s is new:
  round 0   candidates = the family's PERMANENT-network matched filter (all >= 30-event build templates, frozen
            1-8 Hz chain, site merge B, valid-weighted network mean = stack_activity.py 'full' network), local
            maxima declustered +-4 s, top-K per day (K = 2 per scored hour), origin times in the TEST half only
            (the truth templates are from the build half -> disjoint events).  s is not part of it.
  stack     s's 1-8 Hz cache, long window [t-20, t+32] s at every candidate origin t, each window unit-max
            normalised over its 3 components (template-building rule), windows touching invalid data skipped.
  template  wa_templates.py rule on the stack: S = argmax of the 0.5-s smoothed 3-C envelope in [t-5, t+20],
            12 s from S-0.5 s, per-channel L2; SNR weight w (e5r rule, ref = build clear median); clear flag.
  round r>0 candidates re-picked from the permanent network PLUS s (its own site, weight n_ch x w, template of
            round r-1), restack -> round-r template.  --rounds 3 -> rounds 0, 1, 2.
  control   the same procedure on the SITE-ROLLED scrambled control: permanent sites rolled by their own shift
            in [300, 3000] s per day (candidates picked on the rolled network); s unrolled; its own template
            chain.  Expected: no resemblance to the truth template.
  ceiling   s stacked at the family's own TEST-half array detections (what a perfect candidate list gives).

Comparison (score stage), 2-8 Hz Butterworth-4 zero-phase designed at fs = 20 Hz (the cache rate):
  cc_al    mean-over-channels normalised cc between the own truth template (12 s) and the recovered stack,
           best lag within +-2 s of the truth template's own position (origin-aligned)
  cc_any   the same, any lag with the template start in [t-5.5, t+20] s
  rank     of the own truth template among ALL families' truth templates at s (>= 30 build events, any
           clarity) by cc_any; best family; best is a <= 3-km neighbour of f?
  purity   fraction of candidates within +-2 s of an own / any-family array detection (all halves)

  python scripts/merlin_t2.py --stage setup        # cache array station-days, truth stacks, pairs
  python scripts/merlin_t2.py --stage round --round 0|1|2
  python scripts/merlin_t2.py --stage score
  -> results/merlin_tests/t2/
"""
from __future__ import annotations
import os, sys, json, time, argparse
from concurrent.futures import ProcessPoolExecutor, as_completed, ThreadPoolExecutor
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

OD = os.environ.get("MT2_OD", f"{W.RES}/merlin_tests/t2")      # smoke: MT2_OD=results/merlin_tests/t2_smoke
LOG = f"{W.STUDY}/logs/merlin_t2.log"
NL = int(W.LONG_LEN * W.FS)          # 1040
IPRE = int(W.LONG_PRE * W.FS)        # 400
DECL = int(4.0 * W.FS)
WINS = ["V2008", "V2010"]
ARR_NETS = ("XU", "XG")
SERIES = ["data", "control"]


def arr_stations():
    S = pd.concat([pd.read_csv(f"{W.EXT}/{r}_stations.csv") for r in ["NW_WA", "SW_WA", "OLY2010"]])
    S = S[S.net.isin(ARR_NETS)]
    S["netsta"] = S.net + "." + S.sta
    return S.groupby("netsta").agg(lat=("lat", "first"), lon=("lon", "first")).reset_index()


# ------------------------------------------------------------------ setup
def cache_one_arr(job):
    """wa_cache.cache_one with the SH band accepted after HH > BH > EH (XG Array of Arrays records SH? 50 sps)"""
    import wa_cache
    wa_cache.BANDS = ["HH", "BH", "EH", "SH"]
    return wa_cache.cache_one(job)


def stack_station(job):
    """stack build / test half truth events of every family at one array station over its window days"""
    ns, w, days, ev = job
    out = {}
    node0 = None
    for ds in days:
        p = W.cache_path(ns, ds)
        if not os.path.exists(p):
            continue
        try:
            z = np.load(p)
        except Exception:
            continue
        X, m, live = z["x"], z["m"], z["live"]
        node0 = str(z["node"])
        e = ev[ev.date == ds]
        if not len(e):
            continue
        D0 = pd.Timestamp(ds).timestamp()
        i0 = np.round((e.t.to_numpy() - D0) * W.FS).astype(int) - IPRE
        ok = (i0 >= 0) & (i0 + NL <= W.NDAY)
        cm = np.concatenate([[0], np.cumsum(~m)])
        ok[ok] = (cm[i0[ok] + NL] - cm[i0[ok]]) == 0
        if not ok.any():
            continue
        ii = i0[ok]
        seg = X[:, ii[:, None] + np.arange(NL)[None]].transpose(1, 0, 2).astype(np.float64)
        mx = np.abs(seg).max(axis=(1, 2))
        g = mx > 0
        seg = seg[g] / mx[g, None, None]
        fid = e.fid.to_numpy()[ok][g]; hf = e.half.to_numpy()[ok][g]
        for h in (0, 1):
            for f in np.unique(fid[hf == h]):
                s = seg[(fid == f) & (hf == h)]
                k = (int(f), h)
                if k in out:
                    a, n, l = out[k]; out[k] = (a + s.sum(0), n + len(s), l | live)
                else:
                    out[k] = (s.sum(0), len(s), live.copy())
    return ns, w, node0, out


def tpl_from_stack(st, live, ref):
    """wa_templates.make_templates rule on a mean long stack (3, NL) -> dict"""
    from wa_templates import env_smooth, s_snr_db
    chl = np.array([bool(live[c]) and float((st[c] ** 2).sum()) > 0 for c in range(3)])
    if not chl.any():
        return None
    st28 = W.bp(st.astype(np.float64), W.FS, 2.0, 8.0)
    i_a = int((W.LONG_PRE - W.CLEAR_PRE) * W.FS)
    imp, dur = W.env_rule(st28[chl][:, i_a:i_a + int(W.CLEAR_LEN * W.FS)], W.FS)
    env = env_smooth(st[chl].astype(np.float64))
    a_, b_ = IPRE + int(-5 * W.FS), IPRE + int(20 * W.FS)
    k = a_ + int(np.argmax(env[a_:b_]))
    j0 = k - int(W.PRE * W.FS)
    tpl = st[:, j0:j0 + W.NT].astype(np.float64).copy()
    tpl[~chl] = 0.0
    nrm = np.sqrt((tpl ** 2).sum(1))
    chl2 = chl & (nrm > 0)
    tpl[chl2] /= nrm[chl2, None]
    snr = s_snr_db(st[chl].astype(np.float64), k)
    w = float(min(10 ** (snr / 20.0), 3 * ref) / ref) if np.isfinite(snr) else np.nan
    return dict(tpl=tpl.astype(np.float32), offset_s=(j0 - IPRE) / W.FS, live="".join("1" if c else "0" for c in chl2),
                imp=imp, dur=dur, clear=bool(np.isfinite(imp) and imp > 5 and dur < 3), snr_db=snr, w=w)


def wref():
    Q = pd.read_csv(f"{W.RES}/templates/tpl_qc_build.csv")
    Q = Q[Q.clear.fillna(False).astype(bool)]
    return float(np.median(10 ** (Q.snr_db / 20.0)))


def stage_setup(a):
    os.makedirs(OD, exist_ok=True)
    AS = arr_stations()
    F = W.load_families()
    T = pd.read_parquet(f"{W.SETUP}/truth.parquet")
    # 1. cache every XU / XG station-day of the two windows that has a raw trace
    jobs, sdays = [], []
    keep = set(a.days.split(",")) if a.days else None
    for w in WINS:
        for ds, kind in W.window_days(w):
            if keep is not None and ds not in keep:
                continue
            for ns in AS.netsta:
                p = W.raw_path(ns, ds)
                if p:
                    jobs.append((ns, ds, p)); sdays.append(dict(window=w, date=ds, kind=kind, netsta=ns))
    W.log(LOG, f"[t2 setup] {len(jobs)} array station-days with raw traces ({len(set(j[0] for j in jobs))} stations)")
    SD = pd.DataFrame(sdays)
    rows, t0 = [], time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        fut = [ex.submit(cache_one_arr, j) for j in jobs]
        for i, f in enumerate(as_completed(fut)):
            try:
                rows.append(f.result())
            except Exception as e:
                rows.append(dict(status=f"crash:{type(e).__name__}:{str(e)[:80]}"))
            if (i + 1) % 250 == 0:
                W.log(LOG, f"[t2 setup] cache {i+1}/{len(jobs)} {time.time()-t0:.0f}s")
    M = pd.DataFrame(rows)
    M.to_csv(f"{OD}/array_cache_manifest.csv", index=False)
    W.log(LOG, f"[t2 setup] cache {M.status.value_counts().to_dict()}")
    SD.to_csv(f"{OD}/array_station_days.csv", index=False)
    # 2. truth stacks (build + test half) of every family at every array station, per window
    ref = wref()
    jobs = []
    for w in WINS:
        wd = [d for d, _ in W.window_days(w) if keep is None or d in keep]
        ev = T[T.date.isin(wd) & T.half.isin([0, 1])][["fid", "t", "date", "half"]]
        for ns in sorted(SD[SD.window == w].netsta.unique()):
            jobs.append((ns, w, wd, ev))
    rows, stk = [], {}
    with ProcessPoolExecutor(a.workers) as ex:
        fut = [ex.submit(stack_station, j) for j in jobs]
        for i, f in enumerate(as_completed(fut)):
            ns, w, node, out = f.result()
            for (fid, h), (s, n, live) in out.items():
                row = dict(fid=fid, window=w, netsta=ns, node=node, half=h, n_events=n)
                if n >= a.min_events:
                    t = tpl_from_stack(s / n, live, ref)
                    if t is not None:
                        row.update({k: v for k, v in t.items() if k != "tpl"})
                        stk[f"{fid}|{w}|{ns}|{h}"] = (s / n).astype(np.float32)
                rows.append(row)
            if (i + 1) % 20 == 0:
                W.log(LOG, f"[t2 setup] truth stacks {i+1}/{len(jobs)} stations")
    TQ = pd.DataFrame(rows)
    TQ.to_csv(f"{OD}/truth_qc.csv", index=False)
    keys = sorted(stk)
    np.savez(f"{OD}/truth_stacks.npz", keys=np.array(keys), stacks=np.stack([stk[k] for k in keys]))
    # 3. pairs
    B = TQ[(TQ.half == 0) & (TQ.n_events >= a.min_events) & (TQ.clear == True)].copy()  # noqa: E712
    B = B.merge(AS, on="netsta")
    fam = F.set_index("fid")
    B["region"] = B.fid.map(fam.region); B["family_uid"] = B.fid.map(fam.family_uid)
    B["dist_km"] = W.km(B.fid.map(fam.lat_use).to_numpy(float), B.fid.map(fam.lon_use).to_numpy(float),
                        B.lat.to_numpy(), B.lon.to_numpy())
    rng = np.random.default_rng(a.seed)
    sel = []
    for (f, w), g in B.groupby(["fid", "window"]):
        idx = rng.permutation(len(g))[: a.cap]
        sel.append(g.iloc[np.sort(idx)].assign(n_clear_stations=len(g)))
    P = pd.concat(sel, ignore_index=True)
    P["pair"] = np.arange(len(P))
    P.to_csv(f"{OD}/pairs.csv", index=False)
    W.log(LOG, f"[t2 setup] clear truth (family, window, station): {len(B)}; pairs selected {len(P)} "
          f"(cap {a.cap}); by region {P.groupby(['window', 'region']).size().to_dict()}")


# ------------------------------------------------------------------ rounds (GPU)
def load_day(stations, ds):
    D = {}
    for ns in stations:
        p = W.cache_path(ns, ds)
        if os.path.exists(p):
            try:
                z = np.load(p)
                D[ns] = dict(x=z["x"], m=z["m"], live=z["live"], node=str(z["node"]))
            except Exception:
                pass
    return D


def stage_round(a):
    import torch
    import torch.nn.functional as Fnn
    from stack_activity import load_templates
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    r = a.round
    P = pd.read_csv(f"{OD}/pairs.csv")
    if a.max_pairs:
        P = P.head(a.max_pairs)
    S = W.load_stations()
    SITE = dict(zip(S.netsta, S.site))
    Q, T18, _ = load_templates(["build"])
    Q = Q[Q.netsta.isin(SITE)].reset_index(drop=True)
    TQ = {(x.fid, x.netsta, x.node): x for x in Q.itertuples()}
    T = pd.read_parquet(f"{W.SETUP}/truth.parquet")
    Tday = {d: g for d, g in T.groupby("date")}
    prev = {}
    if r > 0:
        z = np.load(f"{OD}/round{r-1}.npz")
        PQ = pd.read_csv(f"{OD}/round{r-1}_qc.csv", dtype={"live": str})
        for x in PQ.itertuples():
            if isinstance(x.live, str) and np.isfinite(x.offset_s) and np.isfinite(x.w):
                prev[(int(x.pair), x.series)] = dict(tpl=z["tpl"][int(x.row)], off=float(x.offset_s), w=float(x.w),
                                                     live=np.array([c == "1" for c in x.live.zfill(3)]))
    n = W.NT
    NV = W.NDAY - n + 1
    NF = int(2 ** np.ceil(np.log2(W.NDAY + n)))
    NB = (NV + DECL - 1) // DECL
    ar = torch.arange(NL, device=dev)
    ref = wref()
    shares = []
    acc = {}          # (pair, series) -> [sum (3, NL) float64, count, purity_own, purity_any, n_cand]
    ceil = {}
    days = [(w, d, k) for w in WINS for d, k in W.window_days(w)]
    if a.days:
        ks = set(a.days.split(","))
        days = [x for x in days if x[1] in ks]
    PW = {w: g for w, g in P.groupby("window")}
    W.log(LOG, f"[t2 r{r}] {len(P)} pairs, {len(days)} days, {len(prev)} previous templates, dev {dev}")

    def guard(X, live, m):
        xd = X.double()
        EN = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64), xd * xd], 1), 1)
        ew = (EN[:, n:] - EN[:, :-n])[:, :NV]
        del EN, xd
        lv = torch.from_numpy(live).to(dev)
        ewl = ew[lv]
        med_ch = ewl.median(dim=1).values
        e3 = ewl.sum(0)
        med3 = e3.median()
        ref_ch = torch.maximum(med_ch, 1e-2 * med3)
        bad = ((ewl < 1e-3 * ref_ch[:, None]).any(0) | (e3 < 1e-3 * med3) | (ewl.min(0).values <= 0))
        if float(med3) <= 0:
            bad = torch.ones_like(bad)
        mk = torch.from_numpy(~m).to(dev)
        cmk = torch.cat([torch.zeros(1, device=dev), torch.cumsum(mk.float(), 0)])
        bad |= (cmk[n:n + NV] - cmk[:NV]) > 0
        if bool(bad.any()):
            bad = Fnn.max_pool1d(bad.float()[None, None], kernel_size=2 * n + 1, stride=1, padding=n)[0, 0] > 0
        efl = 1e-10 * ewl.median()
        return bad, ew.clamp_min(efl).sqrt().float()

    def pick(x, ok, k):
        nr = x.shape[0]
        xp = torch.full((nr, NB * DECL), -float("inf"), device=dev)
        xp[:, :NV] = torch.where(ok, x, torch.full_like(x, -float("inf")))
        bm, ba = xp.view(nr, NB, DECL).max(2)
        del xp
        ninf = torch.full((nr, 1), -float("inf"), device=dev)
        left = torch.cat([ninf, bm[:, :-1]], 1); right = torch.cat([bm[:, 1:], ninf], 1)
        bm2 = torch.where((bm >= left) & (bm > right), bm, torch.full_like(bm, -float("inf")))
        v, ib = torch.topk(bm2, min(k, NB), dim=1)
        return ib * DECL + ba.gather(1, ib), torch.isfinite(v)

    def cc_series(P_, tpl, live_t, off):
        """origin-axis normalised cc of one station template (live channels) and its validity"""
        ch = live_t & P_["live"]
        if not ch.any():
            return None, None
        chi = [c for c in range(3) if ch[c]]
        WH = torch.fft.rfft(torch.from_numpy(tpl[chi].astype(np.float32)).to(dev), n=NF)
        num = torch.fft.irfft(WH.conj() * P_["DH"][chi], n=NF)[:, :NV]
        cc = (num / P_["E"][chi]).mean(0)
        del WH, num
        cc[P_["bad"]] = 0.0
        dd = int(round(off * W.FS))
        a0, b0 = max(0, -dd), min(NV, NV - dd)
        tmp = torch.zeros(NV, device=dev); tv = torch.zeros(NV, device=dev)
        tmp[a0:b0] = cc[a0 + dd:b0 + dd]
        tv[a0:b0] = (~P_["bad"][a0 + dd:b0 + dd]).float()
        return tmp, tv

    def stack_at(P_, idx_t, ok_t):
        """long windows [t-20, t+32] of the station's 1-8 Hz data at origin indices -> (sum, count)"""
        st = idx_t[ok_t] - IPRE
        st = st[(st >= 0) & (st + NL <= W.NDAY)]
        if st.numel() == 0:
            return None, 0
        cmk = P_["cm"]
        st = st[(cmk[st + NL] - cmk[st]) == 0]
        if st.numel() == 0:
            return None, 0
        seg = P_["X"][:, st[:, None] + ar[None]].permute(1, 0, 2)      # (k, 3, NL)
        mx = seg.abs().amax(dim=(1, 2))
        g = mx > 0
        seg = seg[g] / mx[g, None, None]
        return seg.sum(0).double().cpu().numpy(), int(g.sum())

    all_perm = sorted(Q.netsta.unique())
    days = [x for x in days if x[0] in PW]
    cand_rows = []
    ex = ThreadPoolExecutor(1)
    need = {w: all_perm + sorted(PW[w].netsta.unique()) for w in PW}
    fut = ex.submit(load_day, need[days[0][0]], days[0][1]) if days else None
    for it, (w, ds, kind) in enumerate(days):
        t_day = time.time()
        Pw = PW[w]
        arr_need = sorted(Pw.netsta.unique())
        Dn = fut.result()
        if it + 1 < len(days):
            fut = ex.submit(load_day, need[days[it + 1][0]], days[it + 1][1])
        if not any(ns in Dn for ns in arr_need):
            continue
        D0 = pd.Timestamp(ds).timestamp()
        tt = D0 + np.arange(NV) / W.FS
        scored_np = W.half_of(tt) == 1
        scored = torch.from_numpy(scored_np).to(dev)
        K = int(round(a.kph * scored_np.sum() / (3600 * W.FS)))
        PP = {}
        for ns, d in Dn.items():
            X = torch.from_numpy(d["x"]).to(dev)
            bad, E = guard(X, d["live"], d["m"])
            ent = dict(DH=torch.fft.rfft(X, n=NF), E=E, bad=bad, live=d["live"], node=d["node"])
            if ns in arr_need:
                ent["X"] = X
                ent["cm"] = torch.cat([torch.zeros(1, device=dev, dtype=torch.int32),
                                       torch.cumsum(torch.from_numpy(~d["m"]).to(dev).int(), 0)])
            else:
                del X
            PP[ns] = ent
        del Dn
        rng = np.random.default_rng(a.seed + 11 + int(pd.Timestamp(ds).strftime("%Y%m%d")))
        perm_sites = sorted({SITE[ns] for ns in PP if ns in SITE})
        shift = {q: int(rng.uniform(*W.ROLL_S) * W.FS) for q in perm_sites}
        tg = Tday.get(ds, T.iloc[:0])
        t_any = np.sort(tg.t.to_numpy())
        npairs = 0
        for f, gp in Pw.groupby("fid"):
            gp = gp[gp.netsta.isin(PP)]
            if not len(gp):
                continue
            ss = []
            for ns in PP:
                if ns not in SITE:
                    continue
                key = (f, ns, PP[ns]["node"])
                if key in TQ:
                    tl = np.array([c == "1" for c in TQ[key].live])
                    if (tl & PP[ns]["live"]).any():
                        ss.append((ns, key, tl))
            sites = sorted({SITE[ns] for ns, _, _ in ss})
            if len(sites) < 2:
                continue
            si = {q: i for i, q in enumerate(sites)}
            nq = len(sites)
            wm = {ns: float((tl & PP[ns]["live"]).sum()) * float(TQ[k].w) for ns, k, tl in ss}
            Wq = torch.zeros(nq, device=dev)
            for ns, _, _ in ss:
                Wq[si[SITE[ns]]] += wm[ns] ** 2
            Wq = Wq.sqrt()
            snum = torch.zeros((nq, NV), device=dev); sden = torch.zeros((nq, NV), device=dev)
            for ns, k, tl in ss:
                tmp, tv = cc_series(PP[ns], T18[k], tl, TQ[k].offset_s)
                if tmp is None:
                    continue
                q = si[SITE[ns]]
                snum[q] += wm[ns] * tv * tmp; sden[q] += wm[ns] * tv
                del tmp, tv
            vq = sden > 0
            mq = torch.where(vq, snum / sden.clamp_min(1e-20), torch.zeros_like(snum))
            del snum, sden
            Wtot = float(Wq.sum())
            Wv = Wq[:, None] * vq.float()
            base = {"data": ((Wv * mq).sum(0), Wv.sum(0))}
            shl = [shift[q] for q in sites]
            mr = torch.stack([torch.roll(mq[i], shl[i]) for i in range(nq)])
            vr = torch.stack([torch.roll(Wv[i], shl[i]) for i in range(nq)])
            base["control"] = ((vr * mr).sum(0), vr.sum(0))
            del mq, vq, Wv, mr, vr
            own_all = np.sort(tg[tg.fid == f].t.to_numpy())
            if r == 0:
                picks = {}
                for sr in SERIES:
                    nm, dn = base[sr]
                    ix, gd = pick((nm / dn.clamp_min(1e-20))[None], ((dn >= 0.5 * Wtot) & scored)[None], K)
                    picks[sr] = (ix[0], gd[0])
                own_test = tg[(tg.fid == f) & (tg.half == 1)].t.to_numpy()
            for p in gp.itertuples():
                Ps = PP[p.netsta]
                for sr in SERIES:
                    if r == 0:
                        ix, gd = picks[sr]
                    else:
                        pv = prev.get((p.pair, sr))
                        if pv is None:
                            continue
                        nm, dn = base[sr]
                        ws = float((pv["live"] & Ps["live"]).sum()) * pv["w"]
                        tmp, tv = cc_series(Ps, pv["tpl"], pv["live"], pv["off"])
                        if tmp is None:
                            continue
                        nm2 = nm + ws * tv * tmp; dn2 = dn + ws * tv
                        shares.append((int(p.pair), sr, ws / (Wtot + ws)))
                        ix, gd = pick((nm2 / dn2.clamp_min(1e-20))[None],
                                      ((dn2 >= 0.5 * (Wtot + ws)) & scored)[None], K)
                        ix, gd = ix[0], gd[0]
                        del tmp, tv, nm2, dn2
                    s_, c_ = stack_at(Ps, ix, gd)
                    tsel = D0 + ix[gd].cpu().numpy() / W.FS
                    cand_rows.append((int(p.pair), sr, tsel.astype(np.float64)))
                    po = int(np.sum(np.min(np.abs(own_all[None] - tsel[:, None]), 1) <= W.MATCH_S)) if len(own_all) and len(tsel) else 0
                    pa = int(np.sum(np.min(np.abs(t_any[None] - tsel[:, None]), 1) <= W.MATCH_S)) if len(t_any) and len(tsel) else 0
                    e = acc.setdefault((p.pair, sr), [np.zeros((3, NL)), 0, 0, 0, 0, np.zeros(3, bool)])
                    if s_ is not None:
                        e[0] += s_; e[1] += c_
                    e[2] += po; e[3] += pa; e[4] += len(tsel); e[5] |= Ps["live"]
                if r == 0 and len(own_test):
                    ti = torch.from_numpy(np.round((own_test - D0) * W.FS).astype(np.int64)).to(dev)
                    s_, c_ = stack_at(Ps, ti, torch.ones_like(ti, dtype=torch.bool))
                    e = ceil.setdefault(p.pair, [np.zeros((3, NL)), 0, np.zeros(3, bool)])
                    if s_ is not None:
                        e[0] += s_; e[1] += c_; e[2] |= Ps["live"]
                npairs += 1
            del base
        del PP
        if dev == "cuda":
            torch.cuda.empty_cache()
        W.log(LOG, f"[t2 r{r}] {ds} {kind}: {npairs} pair-days, K {K}, {time.time()-t_day:.0f}s")
    # ---- templates of this round
    rows, tpls, stks = [], [], []
    for (pr, sr), (s_, c_, po, pa, nc, lv) in sorted(acc.items()):
        row = dict(pair=pr, series=sr, round=r, n_stacked=c_, n_cand=nc, pur_own=po / nc if nc else np.nan,
                   pur_any=pa / nc if nc else np.nan)
        st = s_ / c_ if c_ else np.zeros((3, NL))
        t = tpl_from_stack(st, lv, ref) if c_ >= 1 else None
        if t is not None:
            row.update({k: v for k, v in t.items() if k != "tpl"}); tp = t["tpl"]
        else:
            tp = np.zeros((3, W.NT), np.float32)
        row["row"] = len(tpls)
        rows.append(row); tpls.append(tp); stks.append(st.astype(np.float32))
    ex.shutdown()
    R = pd.DataFrame(rows)
    if shares:
        SH = pd.DataFrame(shares, columns=["pair", "series", "share"]).groupby(["pair", "series"]).share.mean()
        R["s_weight_share"] = [SH.get((pr, sr), np.nan) for pr, sr in zip(R.pair, R.series)]
    if cand_rows:
        pd.DataFrame(dict(pair=np.concatenate([[p_] * len(t_) for p_, _, t_ in cand_rows]).astype(int),
                          series=np.concatenate([[s_] * len(t_) for _, s_, t_ in cand_rows]),
                          t=np.concatenate([t_ for _, _, t_ in cand_rows]))).to_parquet(f"{OD}/round{r}_cand.parquet", index=False)
    np.savez(f"{OD}/round{r}.npz", tpl=np.stack(tpls) if tpls else np.zeros((0, 3, W.NT), np.float32),
             stack=np.stack(stks) if stks else np.zeros((0, 3, NL), np.float32))
    R.to_csv(f"{OD}/round{r}_qc.csv", index=False)
    if r == 0:
        crow, cst = [], []
        for pr, (s_, c_, lv) in sorted(ceil.items()):
            crow.append(dict(pair=pr, n_stacked=c_, row=len(cst), live="".join("1" if v else "0" for v in lv)))
            cst.append((s_ / max(c_, 1)).astype(np.float32))
        pd.DataFrame(crow).to_csv(f"{OD}/ceiling_qc.csv", index=False)
        np.savez(f"{OD}/ceiling.npz", stack=np.stack(cst) if cst else np.zeros((0, 3, NL), np.float32))
    W.log(LOG, f"[t2 r{r}] wrote {len(R)} (pair, series) stacks; median n_stacked "
          f"{R.groupby('series').n_stacked.median().to_dict() if len(R) else {}}")


# ------------------------------------------------------------------ score (CPU)
def bp28(x):
    return W.bp(np.asarray(x, np.float64), W.FS, 2.0, 8.0)


def tpl28(stack_long, off, live):
    """12-s 2-8 Hz template from a long stack at the given offset, per-channel L2"""
    s28 = bp28(stack_long)
    j0 = IPRE + int(round(off * W.FS))
    t = s28[:, j0:j0 + W.NT].copy()
    t[~live] = 0
    nr = np.sqrt((t ** 2).sum(1))
    ok = live & (nr > 0)
    t[ok] /= nr[ok, None]
    t[~ok] = 0
    return t, ok


def slide_cc(tpl, ok_t, x28, live_x):
    """normalised cc (mean over common live channels) of a 12-s template against every start j of a long
    2-8 Hz stack -> array over j = 0..NL-NT"""
    ch = ok_t & live_x
    if not ch.any():
        return None
    nj = NL - W.NT + 1
    out = np.zeros(nj)
    for c in np.where(ch)[0]:
        num = np.correlate(x28[c], tpl[c], mode="valid")
        e = np.concatenate([[0], np.cumsum(x28[c] ** 2)])
        den = np.sqrt(np.maximum(e[W.NT:] - e[:-W.NT], 1e-30))
        out += num / den
    return out / ch.sum()


def stage_score(a):
    P = pd.read_csv(f"{OD}/pairs.csv")
    F = W.load_families()
    la, lo = F.lat_use.to_numpy(float), F.lon_use.to_numpy(float)
    dk = W.km(la[:, None], lo[:, None], la[None], lo[None])
    fi = {f: i for i, f in enumerate(F.fid)}
    TQ = pd.read_csv(f"{OD}/truth_qc.csv", dtype={"live": str})
    zt = np.load(f"{OD}/truth_stacks.npz")
    TK = {k: i for i, k in enumerate(zt["keys"])}
    TS = zt["stacks"]
    B = TQ[(TQ.half == 0) & (TQ.n_events >= a.min_events) & TQ.offset_s.notna()].copy()
    B["live"] = B.live.astype(str).str.zfill(3)
    # truth templates (2-8 Hz) per (window, station): all families
    TT = {}
    for (w, ns), g in B.groupby(["window", "netsta"]):
        lst = []
        for x in g.itertuples():
            k = f"{x.fid}|{w}|{ns}|0"
            if k not in TK:
                continue
            lv = np.array([c == "1" for c in x.live])
            t, ok = tpl28(TS[TK[k]], x.offset_s, lv)
            lst.append((int(x.fid), t, ok, int(round(x.offset_s * W.FS)) + IPRE))
        TT[(w, ns)] = lst
    # famday carrying info from the round-0 run? use stack_activity valA fam files (all >= 30-event templates)
    import glob
    FDA = pd.concat([pd.read_parquet(p) for p in glob.glob(f"{W.RES}/stack_activity/run_valA/fam_*.parquet")])
    FDA["window"] = np.where(FDA.date < "2009-01-01", "V2008", "V2010")
    carry = FDA.groupby(["fid", "window"]).agg(med_sites=("n_sites", "median"), med_sites_clear=("n_sites_clear", "median")).reset_index()
    P = P.merge(carry, on=["fid", "window"], how="left")
    CQ = pd.read_csv(f"{OD}/ceiling_qc.csv", dtype={"live": str}) if os.path.exists(f"{OD}/ceiling_qc.csv") else pd.DataFrame()
    CZ = np.load(f"{OD}/ceiling.npz")["stack"] if len(CQ) else None
    CQd = {int(x.pair): x for x in CQ.itertuples()} if len(CQ) else {}
    import re
    rounds = sorted(int(m.group(1)) for m in (re.search(r"/round(\d+)\.npz$", p) for p in glob.glob(f"{OD}/round*.npz")) if m)
    RQ = {r: (pd.read_csv(f"{OD}/round{r}_qc.csv", dtype={"live": str}), np.load(f"{OD}/round{r}.npz")["stack"]) for r in rounds}
    lag2 = int(W.MATCH_S * W.FS)
    j_lo, j_hi = IPRE + int(-5.5 * W.FS), IPRE + int(20 * W.FS)

    def compare(p, x_long, live_x):
        x28 = bp28(x_long)
        own = None
        res = []
        for fid, t, ok, j0 in TT.get((p.window, p.netsta), []):
            c = slide_cc(t, ok, x28, live_x)
            if c is None:
                continue
            cany = float(c[j_lo:j_hi + 1].max())
            res.append((fid, cany))
            if fid == p.fid:
                a0, b0 = max(0, j0 - lag2), min(len(c) - 1, j0 + lag2)
                own = (float(c[a0:b0 + 1].max()), float(c[min(max(j0, 0), len(c) - 1)]), cany,
                       int(np.argmax(c[a0:b0 + 1])) + a0 - j0)
        if own is None:
            return {}
        res.sort(key=lambda v: -v[1])
        rank = 1 + sum(1 for f_, c_ in res if f_ != p.fid and c_ > own[2])
        best = res[0][0]
        bd = dk[fi[p.fid], fi[best]] if best in fi else np.nan
        return dict(cc_al=own[0], cc_zero=own[1], cc_any=own[2], lag_s=own[3] / W.FS, rank=rank, n_competitors=len(res),
                    best_fid=best, best_is_own=best == p.fid, best_km=bd, best_is_nbr3=bool(best == p.fid or (np.isfinite(bd) and bd <= 3)))

    rows = []
    for p in P.itertuples():
        base = dict(pair=p.pair, fid=p.fid, family_uid=p.family_uid, region=p.region, window=p.window, netsta=p.netsta,
                    dist_km=p.dist_km, truth_n_build=p.n_events, truth_imp=p.imp, n_clear_stations=p.n_clear_stations,
                    med_sites=p.med_sites, med_sites_clear=p.med_sites_clear)
        c = CQd.get(int(p.pair))
        if c is not None and c.n_stacked > 0:
            rows.append(dict(base, series="ceiling", round=-1, n_stacked=int(c.n_stacked),
                             **compare(p, CZ[int(c.row)], np.array([ch == "1" for ch in c.live.zfill(3)]))))
        for r in rounds:
            R, Z = RQ[r]
            for x in R[R.pair == p.pair].itertuples():
                if x.n_stacked <= 0:
                    continue
                lv = np.array([ch == "1" for ch in x.live.zfill(3)]) if isinstance(x.live, str) else np.ones(3, bool)
                rows.append(dict(base, series=x.series, round=r, n_stacked=int(x.n_stacked), pur_own=x.pur_own,
                                 pur_any=x.pur_any, rec_clear=x.clear, rec_imp=x.imp, **compare(p, Z[int(x.row)], lv)))
    O = pd.DataFrame(rows)
    O.to_csv(f"{OD}/pair_results.csv", index=False)
    # summaries
    O["carry_bin"] = pd.cut(O.med_sites_clear.fillna(0), [-0.5, 0.5, 1.5, 2.5, 99], labels=["0", "1", "2", ">=3"])
    out = []
    for keys, g in [(("ALL",), O)] + [((rg,), g) for rg, g in O.groupby("region")]:
        for (sr, r), h in g.groupby(["series", "round"]):
            out.append(dict(group="region", value=keys[0], series=sr, round=r, n_pairs=len(h),
                            n_families=h.fid.nunique(), med_cc_al=h.cc_al.median(), frac_cc_gt07=(h.cc_al > 0.7).mean(),
                            frac_cc_gt05=(h.cc_al > 0.5).mean(), med_cc_any=h.cc_any.median(),
                            own_best=h.best_is_own.mean(), top3=(h["rank"] <= 3).mean(), best_nbr3=h.best_is_nbr3.mean(),
                            med_rank=h["rank"].median(), med_competitors=h.n_competitors.median(),
                            med_pur_own=h.pur_own.median() if "pur_own" in h else np.nan,
                            med_n_stacked=h.n_stacked.median()))
    for cb, g in O.groupby("carry_bin", observed=True):
        for (sr, r), h in g.groupby(["series", "round"]):
            out.append(dict(group="perm_clear_sites", value=str(cb), series=sr, round=r, n_pairs=len(h),
                            n_families=h.fid.nunique(), med_cc_al=h.cc_al.median(), frac_cc_gt07=(h.cc_al > 0.7).mean(),
                            frac_cc_gt05=(h.cc_al > 0.5).mean(), med_cc_any=h.cc_any.median(),
                            own_best=h.best_is_own.mean(), top3=(h["rank"] <= 3).mean(), best_nbr3=h.best_is_nbr3.mean(),
                            med_rank=h["rank"].median(), med_competitors=h.n_competitors.median(),
                            med_pur_own=h.pur_own.median() if "pur_own" in h else np.nan, med_n_stacked=h.n_stacked.median()))
    S = pd.DataFrame(out)
    S.to_csv(f"{OD}/summary.csv", index=False)
    W.log(LOG, "[t2 score]\n" + S[S.group == "region"].to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["setup", "round", "score"], required=True)
    ap.add_argument("--round", type=int, default=0)
    ap.add_argument("--workers", type=int, default=14)
    ap.add_argument("--cap", type=int, default=6, help="max withheld stations per (family, window)")
    ap.add_argument("--seed", type=int, default=20260930)
    ap.add_argument("--kph", type=float, default=2.0)
    ap.add_argument("--days", default="", help="smoke: comma list of dates")
    ap.add_argument("--max-pairs", type=int, default=0, help="smoke: first N pairs")
    ap.add_argument("--min-events", type=int, default=W.MIN_EVENTS, help="truth stack minimum (30; smoke only lower)")
    a = ap.parse_args()
    {"setup": stage_setup, "round": stage_round, "score": stage_score}[a.stage](a)


if __name__ == "__main__":
    main()
