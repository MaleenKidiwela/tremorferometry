#!/usr/bin/env python
"""Avenue A4 -- LOCATION-SEEDED POSITIVE CONTROL for Royer & Bostock (2014) NW / SW / NVI families.

Question: do repeating LFE families exist in OUR 40 sps 1-8 Hz data at THEIR template locations?
Independent of the literal detection chain (ab_*); it only reads that chain's day caches.

Per location (lon, lat, depth) of one condition:
  0. stations: the subarray's Fig.1D pool stations that have cache on the autocorrelation-episode
     days, ranked by epicentral distance; the nearest K_TPL (=10, the paper's "7 to 10") carry the
     template, the nearest K_LONG (=16) get diagnostic long stacks for location.
     Predicted P/S times: GSC onshore 1-D model (a4f_hypo.MODELS['GSC']), Vp/Vs 1.73, direct ray or
     head wave (first arrival), station elevation ignored.
  1. BEAM: per station-day, log envelope (0.5 s smoothing, 10 Hz) of the horizontals and of the
     vertical, minus its own 30 s running mean (a local log-STA/LTA, so tremor background and gain
     cancel).  Beam(tau) = mean over valid stations of [eH(tau + tS) + eZ(tau + tP)] / 2 (>= 4 valid
     stations).  Per day: robust z = (beam - median) / MAD.  ATTRIBUTION: the same beam is formed
     for every location of the subarray's reference set (all 122/54/47 locations of the condition),
     max-pooled +-2 s; a burst belongs to the location whose pooled z is the highest (argmax, a
     back-projection location).  Per location, the 20 largest own-attributed z peaks per day (15 s
     apart), then the M_CAND (=60) largest over all days (ranking; z and the margin over the
     runner-up location, `contrast`, are columns).
  2. SEED: every candidate is cut per template station, 3-C, window [tau + tP - 1 s, +L] with
     L = max(tS - tP) + 6 s.  Pairwise network similarity = max over a common lag (+-1 s) of the
     per-channel normalised cc averaged over the channels both have.  Seed = the candidate with the
     largest mean of its 5 best similarities; seed template = unit-max-3C stack of it and its 9 most
     similar candidates, each shifted by its best lag.  (Ranking only; the similarities are columns.)
  3. MATCHED FILTER (per-channel normalised cc averaged over valid channels, validity-masked: an
     invalid window gets an INFINITE denominator -> zero, never a clamped one), 8 x MAD over each
     24 h day, declustered by L, >= 4 valid stations; every detection restacked (unit max across
     3 C, Royer sec.3.3); iterate: seed -> MF0 -> stack -> MF1 -> stack -> MF2 -> final stack.
  4. DIAGNOSTICS (columns): detections/day, nmad/cc of detections, per-station stack SNR, observed
     vs predicted S (horizontal envelope peak) and P (vertical envelope peak before S), and a
     Hyp2000 location of the final long stacks (a4f_hypo.make_picks / write_hyp_inputs / run_hyp,
     GSC, Vp/Vs 1.73, R 1.5 -- the chain's headline) with distance from the seed location.

Conditions (identical code path):
  real   -- real locations, real data;
  nullT  -- real locations, TIME-SHIFTED station data: each station reads a different episode day
            (distinct seeded whole-day shift within the subarray's day list) AND is circularly rolled
            within the day by a distinct station offset (600 + 997*i s), so every station has its own
            time mapping; gaps, diurnal level and tremor level are kept, inter-station coherence is
            destroyed;
  nullR  -- rotated/shifted locations (the location file carries them), real data.

  python a4_seed.py --sub NW --locs results/a4_seed/locs_pilot.csv --cond real --out results/a4_seed/pilot
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
import zlib

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import a4f_common as A                                                      # noqa: E402
import a4f_hypo as H                                                        # noqa: E402

FS = A.FS
K_TPL, K_LONG = 10, 16
EFS = 10.0                 # envelope beam sampling
ENV_SMOOTH = 0.5
LTA_S = 30.0
PEAKS_PER_DAY = 20
DECL_BEAM = 15.0
M_CAND = 60
POOL_DT = 2.0
RANK_BY = "z"
TOPSIM = 5
N_SEED = 10
LAG_S = 1.0
PRE_P = 1.0
POST_S = 5.0
LONG_PRE, LONG_LEN = 2.0, 32.0      # diagnostic stack: origin - 2 s .. origin + 30 s
MIN_STA_VALID = 4
NMAD = 8.0
N_MF = 3


# ------------------------------------------------------------------ travel times ---------------
def tt_1d(x_km, z_km, vpvs=1.73, model="GSC", phase="P"):
    """First-arrival time (s) from a source at depth z to a surface receiver at epicentral x:
    min(direct ray by bisection on the ray parameter, head waves on interfaces below the source)."""
    zt, vp = H.MODELS[model]
    v = np.array(vp, float) / (vpvs if phase == "S" else 1.0)
    top = np.array(zt, float)
    bot = np.r_[top[1:], np.inf]
    # thicknesses above the source
    h = np.clip(np.minimum(bot, z_km) - top, 0, None)
    use = h > 0
    hh, vv = h[use], v[use]
    vmax = vv.max()
    x = float(x_km)
    if x < 1e-6:
        t_dir = float((hh / vv).sum())
    else:
        lo, hi = 0.0, (1.0 / vmax) * (1 - 1e-12)
        for _ in range(80):
            p = 0.5 * (lo + hi)
            X = (hh * p * vv / np.sqrt(1 - (p * vv) ** 2)).sum()
            if X < x:
                lo = p
            else:
                hi = p
        p = 0.5 * (lo + hi)
        t_dir = float((hh / (vv * np.sqrt(1 - (p * vv) ** 2))).sum())
    best = t_dir
    # head waves along interfaces below the source
    for k in range(len(top)):
        if top[k] <= z_km:
            continue
        vk = v[k]
        above = v[:k]
        if vk <= above.max():
            continue
        p = 1.0 / vk
        hdown = np.clip(np.minimum(bot[:k], top[k]) - np.maximum(top[:k], z_km), 0, None)
        hup = np.clip(np.minimum(bot[:k], top[k]) - top[:k], 0, None)
        q = np.sqrt(1.0 / above ** 2 - p ** 2)
        xc = ((hdown + hup) * p / q).sum()
        if x >= xc:
            best = min(best, float(x * p + ((hdown + hup) * q).sum()))
    return best


# ------------------------------------------------------------------ data ------------------------
def remap(sub, sta, date):
    """nullT: (source date, roll samples) for this station."""
    days = A.days(sub)
    n = len(days)
    pool = list(A.POOL[sub])
    i = pool.index(sta) if sta in pool else len(pool) + sorted(A.COORD).index(sta)
    rng = np.random.default_rng(zlib.crc32(f"A4nullT|{sub}".encode()))
    sh = rng.permutation(np.arange(1, n))
    h = int(sh[i % sh.size])
    src = days[(days.index(date) + h) % n]
    roll = int(((600 + 997 * i) % 86400) * FS)
    return src, roll


def load(sub, sta, date, cond):
    src, roll = (remap(sub, sta, date) if cond == "nullT" else (date, 0))
    p = A.cache_path(sta, src)
    if not os.path.exists(p):
        return None
    try:
        z = np.load(p)
        X = np.ascontiguousarray(z["x"], np.float32)
        m = z["m"].astype(bool) if "m" in z.files else np.ones(X.shape[1], bool)
        t0 = float(z["t0"])
    except Exception as e:          # a file being rewritten by the other chain: skip it
        A.log(f"[load] {sta} {src}: {e}")
        return None
    if X.shape[0] != 3:
        return None
    if roll:
        X = np.roll(X, roll, axis=1)
        m = np.roll(m, roll)
    t0 += (pd.Timestamp(date) - pd.Timestamp(src)).total_seconds()
    return t0, X, m


def load_dayset(sub, stas, date, cond):
    D = {}
    for s in stas:
        v = load(sub, s, date, cond)
        if v is not None:
            D[s] = v
    return D


def cut(D, sta, t_abs, n):
    v = D.get(sta)
    if v is None:
        return None
    t0, X, m = v
    i0 = int(round((t_abs - t0) * FS))
    if i0 < 0 or i0 + n > X.shape[1] or not m[i0:i0 + n].all():
        return None
    seg = X[:, i0:i0 + n].astype(np.float64)
    return seg if np.isfinite(seg).all() else None


def unit_max3(seg):
    mx = float(np.abs(seg).max())
    return seg / mx if mx > 0 else seg


# ------------------------------------------------------------------ envelopes ------------------
def log_env(X, m, dev):
    """(eH, eZ) at EFS Hz: log envelope minus its own 30 s running mean; NaN where invalid."""
    import torch
    N = X.shape[1]
    x = torch.from_numpy(X).to(dev, torch.float32)
    F = torch.fft.fft(x, dim=1)
    hwin = torch.zeros(N, device=dev)
    hwin[0] = 1
    hwin[1:(N + 1) // 2] = 2
    if N % 2 == 0:
        hwin[N // 2] = 1
    a = torch.fft.ifft(F * hwin[None], dim=1).abs()
    del F, x
    eH = torch.sqrt(a[1] ** 2 + a[2] ** 2)
    eZ = a[0]
    dec = int(round(FS / EFS))
    k = int(round(ENV_SMOOTH * FS))
    mm = torch.from_numpy(m).to(dev)
    out = []
    for e in (eH, eZ):
        e = torch.where(mm, e, torch.zeros_like(e))
        c = torch.nn.functional.avg_pool1d(e[None, None], k, stride=1, padding=k // 2)[0, 0][:N]
        c = c[::dec]
        valid = torch.nn.functional.avg_pool1d(mm.float()[None, None], k, stride=1,
                                               padding=k // 2)[0, 0][:N][::dec] > 0.999
        valid &= c > 0
        le = torch.where(valid, torch.log(c.clamp_min(1e-30)), torch.zeros_like(c))
        w = int(round(LTA_S * EFS))
        s1 = torch.nn.functional.avg_pool1d(le[None, None], w, stride=1, padding=w // 2)[0, 0][:le.numel()]
        s0 = torch.nn.functional.avg_pool1d(valid.float()[None, None], w, stride=1,
                                            padding=w // 2)[0, 0][:le.numel()]
        r = le - s1 / s0.clamp_min(1e-6)
        r = torch.where(valid & (s0 > 0.5), r, torch.full_like(r, float("nan")))
        out.append(r.cpu().numpy())
    return out[0], out[1]


# ------------------------------------------------------------------ matched filter -------------
def matched_filter(D, TPL, date, dev, nmad=NMAD, min_sta_valid=MIN_STA_VALID, batch=24):
    """Per-channel normalised cc averaged over valid channels; 8 x MAD per day; declustered by each
    template's own length.  TPL = {lid: {sta: (tpl (3,n), off_s)}} (all templates of one location
    share n).  Invalid windows: infinite denominator (-> cc 0, channel not counted)."""
    import torch
    stas = sorted({s for g in TPL for s in TPL[g]})
    have = [s for s in stas if s in D]
    if len(have) < min_sta_valid:
        return pd.DataFrame()
    N = min(D[s][1].shape[1] for s in have)
    t0s = {s: D[s][0] for s in have}
    tref = t0s[have[0]]
    out = []
    # group templates by length so each group shares NV
    bylen = {}
    for g, T in TPL.items():
        n = next(iter(T.values()))[0].shape[1]
        bylen.setdefault(n, []).append(g)
    Xf = {}
    for s in have:
        Xf[s] = torch.fft.rfft(torch.from_numpy(np.ascontiguousarray(D[s][1][:, :N], np.float64)).to(dev), n=N)
    for nt, gs_all in bylen.items():
        NV = N - nt + 1
        DEN, OK = {}, {}
        for s in have:
            t0, X, m = D[s]
            ew, vw = A.window_stats(X[:, :N], m[:N], nt, 1)
            e = torch.from_numpy(ew).to(dev)
            okc = torch.from_numpy(vw).to(dev)[None] & (e > 0)
            DEN[s] = torch.where(okc, e.clamp_min(0).sqrt(), torch.full_like(e, float("inf"))).float()
            OK[s] = okc
            del e
        for c0 in range(0, len(gs_all), 64):
            gs = gs_all[c0:c0 + 64]
            acc = torch.zeros((len(gs), NV), device=dev)
            cch = torch.zeros((len(gs), NV), device=dev)
            csta = torch.zeros((len(gs), NV), device=dev)
            for b0 in range(0, len(gs), batch):
                sub = gs[b0:b0 + batch]
                for s in have:
                    idx = [i for i, g in enumerate(sub) if s in TPL[g]]
                    if not idx:
                        continue
                    M = np.zeros((len(idx), 3, N), np.float64)
                    offs, nrm = [], []
                    for r, i in enumerate(idx):
                        tpl, off = TPL[sub[i]][s]
                        M[r, :, :nt] = tpl
                        # station clock offset (all caches start at midnight; kept general)
                        offs.append(off + (tref - t0s[s]))
                        nrm.append(np.sqrt((tpl ** 2).sum(1)))
                    Tf = torch.fft.rfft(torch.from_numpy(M).to(dev), n=N)
                    cc = torch.fft.irfft(Xf[s][None] * Tf.conj(), n=N)[..., :NV].float()
                    del Tf
                    nrmt = torch.tensor(np.array(nrm), device=dev, dtype=torch.float32)
                    tok = (nrmt > 0)[:, :, None] & OK[s][None]
                    cc = cc / (nrmt[:, :, None] * DEN[s][None])
                    cc = torch.where(tok, cc, torch.zeros_like(cc))
                    ccs = cc.sum(1)
                    nch = tok.float().sum(1)
                    del cc, tok
                    for r, i in enumerate(idx):
                        k = int(round(offs[r] * FS))
                        lo, hi = max(-k, 0), min(NV - k, NV)
                        if hi <= lo:
                            continue
                        acc[b0 + i, lo:hi] += ccs[r, lo + k:hi + k]
                        cch[b0 + i, lo:hi] += nch[r, lo + k:hi + k]
                        csta[b0 + i, lo:hi] += (nch[r, lo + k:hi + k] > 0).float()
                    del ccs, nch
            enough = csta >= float(min_sta_valid)
            fn = acc / cch.clamp_min(1.0)
            fn = torch.where(enough, fn, torch.full_like(fn, float("nan")))
            if not bool(enough.any()):
                continue
            med = torch.nanmedian(fn, dim=1, keepdim=True).values
            madv = torch.nanmedian((fn - med).abs(), dim=1, keepdim=True).values
            thr = med + nmad * madv.clamp_min(1e-9)
            dec = nt
            hit = torch.nan_to_num(fn, nan=-1e9) >= thr
            vals = torch.where(hit, fn, torch.full_like(fn, -1e9))
            pooled = torch.nn.functional.max_pool1d(vals[None], 2 * dec + 1, stride=1, padding=dec)[0]
            keep = hit & (vals >= pooled - 1e-9)
            ij = torch.nonzero(keep, as_tuple=False)
            if ij.numel():
                gi = ij[:, 0].cpu().numpy()
                ti = ij[:, 1].cpu().numpy()
                v = fn[ij[:, 0], ij[:, 1]].cpu().numpy()
                nm = ((fn[ij[:, 0], ij[:, 1]] - med[ij[:, 0], 0]) / madv[ij[:, 0], 0].clamp_min(1e-9)).cpu().numpy()
                ns = csta[ij[:, 0], ij[:, 1]].cpu().numpy().astype(np.int16)
                out.append(pd.DataFrame(dict(lid=[gs[i] for i in gi], t=tref + ti / FS, cc=v, nmad=nm,
                                             n_sta=ns, date=date)))
            # per-template day stats (for the null-comparable "how far is the tail" column)
            day_mad = madv[:, 0].cpu().numpy()
            day_med = med[:, 0].cpu().numpy()
            for i, g in enumerate(gs):
                out.append(pd.DataFrame(dict(lid=[g], t=[np.nan], cc=[np.nan], nmad=[np.nan], n_sta=[-1],
                                             date=[date], day_med=[day_med[i]], day_mad=[day_mad[i]],
                                             day_max_nmad=[float(((torch.nan_to_num(fn[i], nan=-1e9).max() - med[i, 0])
                                                                  / madv[i, 0].clamp_min(1e-9)).item())])))
            del acc, cch, csta, fn, hit, vals, pooled, keep, ij
        del DEN, OK
    del Xf
    if dev == "cuda":
        torch.cuda.empty_cache()
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame()


# ------------------------------------------------------------------ similarity -----------------
def pair_sim(W, Wp, L):
    """W: {sta: (M,3,n) cores or NaN}, Wp: {sta: (M,3,n+2L) padded}.  Returns (sim (M,M), lag (M,M)
    samples; lag = shift of j that best aligns it on i).  Channels valid in both are averaged."""
    M = next(iter(W.values())).shape[0]
    nl = 2 * L + 1
    num = np.zeros((M, M, nl))
    cnt = np.zeros((M, M, nl))
    for s in W:
        a, b = W[s], Wp[s]
        for c in range(3):
            x, y = a[:, c], b[:, c]                     # (M,n), (M,n+2L)
            okx = np.isfinite(x).all(1) & ((np.nan_to_num(x) ** 2).sum(1) > 0)
            oky = np.isfinite(y).all(1) & ((np.nan_to_num(y) ** 2).sum(1) > 0)
            if okx.sum() < 2:
                continue
            x = np.nan_to_num(x)
            y = np.nan_to_num(y)
            n = x.shape[1]
            nf = 1 << int(np.ceil(np.log2(y.shape[1] + n)))
            Fx = np.fft.rfft(x, nf)
            Fy = np.fft.rfft(y, nf)
            r = np.fft.irfft(np.conj(Fx)[:, None] * Fy[None], nf)[..., :nl]   # (M,M,nl)
            cs = np.concatenate([np.zeros((M, 1)), np.cumsum(y ** 2, 1)], 1)
            ey = cs[:, n:n + nl] - cs[:, :nl]                                  # (M,nl)
            ex = (x ** 2).sum(1)
            den = np.sqrt(ex[:, None, None] * ey[None])
            ok = okx[:, None, None] & oky[None, :, None] & (den > 0)
            num += np.where(ok, r / np.where(den > 0, den, 1), 0)
            cnt += ok
    cc = num / np.maximum(cnt, 1)
    cc[cnt < 3] = np.nan
    k = np.nanargmax(np.nan_to_num(cc, nan=-9), axis=2)
    sim = np.take_along_axis(cc, k[..., None], 2)[..., 0]
    return sim, k - L


# ------------------------------------------------------------------ main -----------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", required=True)
    ap.add_argument("--locs", required=True)
    ap.add_argument("--cond", required=True, choices=["real", "nullT", "nullR"])
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref", default=os.path.join(os.path.dirname(HERE), "results/a4_seed/locs_all.csv"),
                    help="reference location set for burst attribution (all locations of the subarray)")
    ap.add_argument("--days", default="", help="comma list; default = the sec.3.1 episode days")
    a = ap.parse_args()
    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    torch.backends.cuda.matmul.allow_tf32 = False
    out = f"{a.out}/{a.sub}_{a.cond}"
    os.makedirs(out, exist_ok=True)
    t_start = time.time()
    Lc = pd.read_csv(a.locs)
    Lc = Lc[(Lc.subarray == a.sub) & (Lc.cond == ("nullR" if a.cond == "nullR" else "real"))].copy()
    Lc = Lc.reset_index(drop=True)
    days = a.days.split(",") if a.days else A.days(a.sub)
    pool = [s for s in A.POOL[a.sub] if s in A.COORD]
    avail = {s: sum(os.path.exists(A.cache_path(s, d)) for d in days) for s in pool}
    pool = [s for s in pool if avail[s] > 0]
    A.log(f"[A4 {a.sub} {a.cond}] {len(Lc)} locations, {len(days)} days, pool {len(pool)} stations, dev {dev}")

    # --- stations + predicted times
    G = {}
    rows = []
    for r in Lc.itertuples():
        d = {s: float(A.epi_km(r.lat, r.lon, A.COORD[s][0], A.COORD[s][1])) for s in pool}
        near = sorted(pool, key=lambda s: d[s])
        tpl_st, long_st = near[:K_TPL], near[:K_LONG]
        tP = {s: tt_1d(d[s], r.depth, phase="P") for s in long_st}
        tS = {s: tt_1d(d[s], r.depth, phase="S") for s in long_st}
        G[r.lid] = dict(tpl=tpl_st, long=long_st, tP=tP, tS=tS, dist=d)
        rows.append(dict(lid=r.lid, tpl_stations=",".join(tpl_st),
                         dist_med_km=float(np.median([d[s] for s in tpl_st])),
                         dist_max_km=float(max(d[s] for s in tpl_st))))
    Sinfo = pd.DataFrame(rows)
    need = sorted({s for g in G.values() for s in g["long"]})
    # template length per location
    for lid, g in G.items():
        # template length, rounded UP to a 2 s multiple so locations share a few MF length groups
        Ls = max(g["tS"][s] - g["tP"][s] for s in g["tpl"]) + PRE_P + POST_S
        g["L"] = int(round(2.0 * np.ceil(Ls / 2.0) * FS))
        g["off"] = {s: g["tP"][s] - PRE_P for s in g["long"]}

    # --- 1. beam (GPU), each burst ATTRIBUTED to the reference location where its beam is highest.
    # Reference set = ALL locations of this subarray and condition in --ref (locs_all.csv: 122 NW /
    # 54 SW / 47 NVI; for nullR the rotated cloud), independent of which subset is run, so the pilot
    # and the full run use the same attribution.  Per reference point: beam on its own nearest 10
    # stations, per-day robust z, max-pooled over +-POOL_DT s (origin time trades off with position).
    # A location's candidate series = its z where its pooled z is the maximum over the reference set
    # (argmax attribution = a back-projection location of the burst; not a threshold), ranked by z.
    # Why: pilot v1 ranked the raw beam and 7/10 NW seeds were the same 2011-08-19 05:47 burst; a ring
    # contrast (smoke2) killed the bright events and left noise (seed similarity at the null level).
    bp = f"{out}/beam_cands.parquet"
    if os.path.exists(bp):
        C = pd.read_parquet(bp)
    else:
        import torch
        Rf = pd.read_csv(a.ref)
        Rf = Rf[(Rf.subarray == a.sub) & (Rf.cond == ("nullR" if a.cond == "nullR" else "real"))].reset_index(drop=True)
        ref_st, ref_k = [], []
        for r in Rf.itertuples():
            d = {s_: float(A.epi_km(r.lat, r.lon, A.COORD[s_][0], A.COORD[s_][1])) for s_ in pool}
            st_ = sorted(pool, key=lambda s_: d[s_])[:K_TPL]
            ref_st.append(st_)
            ref_k.append([(int(round(tt_1d(d[s_], r.depth, phase="S") * EFS)),
                           int(round(tt_1d(d[s_], r.depth, phase="P") * EFS))) for s_ in st_])
        rid = {int(l): i for i, l in enumerate(Rf.lid)}
        missing = [l for l in G if l not in rid]
        if missing:
            raise SystemExit(f"locations {missing} are not in the reference set {a.ref}")
        need_ref = sorted(set(need) | {s_ for st_ in ref_st for s_ in st_})
        cand = []
        w = int(round(DECL_BEAM * EFS))
        wp = int(round(POOL_DT * EFS))
        for date in days:
            D = load_dayset(a.sub, need_ref, date, a.cond)
            if len(D) < MIN_STA_VALID:
                continue
            t0 = D[next(iter(D))][0]
            E = {}
            for s_, (ts, X, m) in D.items():
                eH, eZ = log_env(X, m, dev)
                E[s_] = (eH, eZ, int(round((ts - t0) * EFS)))
            NE = min(v[0].size for v in E.values())
            EH = {s_: torch.from_numpy(np.ascontiguousarray(v[0][:NE])).to(dev) for s_, v in E.items()}
            EZ = {s_: torch.from_numpy(np.ascontiguousarray(v[1][:NE])).to(dev) for s_, v in E.items()}
            SH = {s_: v[2] for s_, v in E.items()}

            def shifted(e, k):
                o = torch.full_like(e, float("nan"))
                lo, hi = max(-k, 0), min(NE, NE - k)
                if hi > lo:
                    o[lo:hi] = e[lo + k:hi + k]
                return o

            ZP = torch.full((len(Rf), NE), -1e9, device=dev)
            Zs = {}
            for i, st_ in enumerate(ref_st):
                acc = torch.zeros(NE, device=dev)
                cnt = torch.zeros(NE, device=dev)
                for s_, (kS, kP) in zip(st_, ref_k[i]):
                    if s_ not in E:
                        continue
                    for e, k in ((EH[s_], kS), (EZ[s_], kP)):
                        seg = shifted(e, k - SH[s_])
                        ok = torch.isfinite(seg)
                        acc += torch.where(ok, seg, torch.zeros_like(seg))
                        cnt += ok.float() * 0.5
                b = torch.where(cnt >= MIN_STA_VALID, acc / (2 * cnt).clamp_min(1), torch.full_like(acc, float("nan")))
                if not bool(torch.isfinite(b).any()):
                    continue
                med = torch.nanmedian(b)
                mad = torch.nanmedian((b - med).abs()).clamp_min(1e-9)
                z = (b - med) / mad
                q = torch.nan_to_num(z, nan=-1e9)
                ZP[i] = torch.nn.functional.max_pool1d(q[None, None], 2 * wp + 1, stride=1, padding=wp)[0, 0]
                if int(Rf.lid[i]) in G:
                    Zs[int(Rf.lid[i])] = (z, cnt)
            top2 = torch.topk(ZP, 2, dim=0)
            for lid, (z, cnt) in Zs.items():
                i = rid[lid]
                other = torch.where(top2.indices[0] == i, top2.values[1], top2.values[0])
                marg = ZP[i] - other
                own = (marg >= 0) & torch.isfinite(z)
                for mode, sc in (("top", torch.nan_to_num(z, nan=-1e9).clone()),
                                 ("attr", torch.where(own, z, torch.full_like(z, -1e9)))):
                    for _ in range(PEAKS_PER_DAY):
                        j = int(torch.argmax(sc))
                        if float(sc[j]) <= -1e8:
                            break
                        cand.append(dict(lid=lid, mode=mode, date=date, t=t0 + j / EFS, z=float(z[j]),
                                         contrast=float(marg[j]), own=bool(own[j]), n_sta=float(cnt[j])))
                        sc[max(j - w, 0):j + w + 1] = -1e9
            del EH, EZ, ZP
            A.log(f"[beam] {date} {len(D)} sta; {len(cand)} cands")
        C = pd.DataFrame(cand)
        C.to_parquet(bp, index=False)

    # two SEED MODES share every later pass (one day load): key = 10 * lid + mode
    #   top  (0): the location's own beam, largest z (pilot v1 behaviour; bright bursts can be shared)
    #   attr (1): only bursts attributed to this location (argmax over the reference set)
    MODES = {"top": 0, "attr": 1}
    C = C.copy()
    C["lid0"] = C.lid
    C["lid"] = 10 * C.lid0 + C["mode"].map(MODES)
    G = {10 * l + v: g for l, g in G.items() for v in MODES.values()}

    # --- 2. seed
    sp = f"{out}/seed.npz"
    seed_rows = []
    TPL = {}
    if os.path.exists(sp):
        z = np.load(sp, allow_pickle=True)
        TPL = json.loads(str(z["tpl"]), object_hook=None)
        TPL = {int(l): {s: (np.array(v), G[int(l)]["off"][s]) for s, v in T.items()} for l, T in TPL.items()}
        seed_rows = pd.read_csv(f"{out}/seed_info.csv").to_dict("records")
    else:
        top = C.sort_values(RANK_BY, ascending=False).groupby("lid").head(M_CAND) if len(C) else C
        Lg = int(round(LAG_S * FS))
        W, Wp = {}, {}
        for lid, g in G.items():
            M = int((top.lid == lid).sum()) if len(top) else 0
            W[lid] = {s: np.full((M, 3, g["L"]), np.nan) for s in g["tpl"]}
            Wp[lid] = {s: np.full((M, 3, g["L"] + 2 * Lg), np.nan) for s in g["tpl"]}
        if len(top):
            top = top.copy()
            top["k"] = top.groupby("lid").cumcount()
            for date, dd in top.groupby("date"):
                D = load_dayset(a.sub, need, date, a.cond)
                for r in dd.itertuples():
                    g = G[r.lid]
                    for s in g["tpl"]:
                        c1 = cut(D, s, r.t + g["off"][s], g["L"])
                        c2 = cut(D, s, r.t + g["off"][s] - LAG_S, g["L"] + 2 * Lg)
                        if c1 is not None and c2 is not None:
                            W[r.lid][s][r.k] = c1
                            Wp[r.lid][s][r.k] = c2
        for lid, g in G.items():
            T = top[top.lid == lid] if len(top) else top
            M = len(T)
            row = dict(lid=lid, n_cand=M, beam_z_max=float(T.z.max()) if M else np.nan,
                       beam_z_med_top=float(T.z.median()) if M else np.nan,
                       beam_contrast_med_top=float(T.contrast.median()) if M else np.nan)
            if M < 3:
                seed_rows.append(row)
                continue
            sim, lag = pair_sim(W[lid], Wp[lid], Lg)
            np.fill_diagonal(sim, np.nan)
            simS = np.nanmean(np.stack([sim, sim.T]), 0)
            srt = -np.sort(-np.nan_to_num(simS, nan=-9), 1)
            sc = srt[:, :TOPSIM].mean(1)
            i0 = int(np.argmax(sc))
            order = np.argsort(-np.nan_to_num(sim[i0], nan=-9))
            mem = [i0] + [int(j) for j in order if j != i0 and np.isfinite(sim[i0, j])][:N_SEED - 1]
            T_ = {}
            for s in g["tpl"]:
                accs, n = np.zeros((3, g["L"])), 0
                for j in mem:
                    lj = 0 if j == i0 else int(lag[i0, j])
                    seg = Wp[lid][s][j][:, Lg + lj:Lg + lj + g["L"]]
                    if np.isfinite(seg).all() and (seg ** 2).sum() > 0:
                        accs += unit_max3(seg)
                        n += 1
                if n >= 1:
                    T_[s] = (accs / n, g["off"][s])
            if len(T_) >= MIN_STA_VALID:
                TPL[lid] = T_
            offd = simS[np.triu_indices(M, 1)]
            row.update(seed_score=float(sc[i0]), seed_sim_members=float(np.nanmean([sim[i0, j] for j in mem[1:]])),
                       cand_sim_median=float(np.nanmedian(offd)), cand_sim_p95=float(np.nanpercentile(offd, 95)),
                       seed_excess=float(sc[i0] - np.nanmedian(offd)),
                       seed_t=float(T.t.iloc[i0]), seed_z=float(T.z.iloc[i0]),
                       seed_contrast=float(T.contrast.iloc[i0]), n_seed_sta=len(T_))
            seed_rows.append(row)
        np.savez_compressed(sp, tpl=np.array(json.dumps({str(l): {s: v[0].tolist() for s, v in T.items()}
                                                          for l, T in TPL.items()})))
        pd.DataFrame(seed_rows).to_csv(f"{out}/seed_info.csv", index=False)
    A.log(f"[seed] {len(TPL)}/{len(G)} locations have a seed template")

    # --- 3. MF + restack iterations
    it_rows = []
    LONGN = int(round(LONG_LEN * FS))
    for it in range(N_MF):
        dp = f"{out}/det_it{it}.parquet"
        lp = f"{out}/long_it{it}.npz"
        if os.path.exists(dp) and os.path.exists(lp) and os.path.exists(f"{out}/tpl_it{it + 1}.json"):
            Dt = pd.read_parquet(dp)
            TPL = {int(l): {s: (np.array(v), G[int(l)]["off"][s]) for s, v in T.items()}
                   for l, T in json.load(open(f"{out}/tpl_it{it + 1}.json")).items()}
            continue
        if not TPL:
            break
        acc, cnt, lacc, lcnt = {}, {}, {}, {}
        dets = []
        for date in days:
            D = load_dayset(a.sub, need, date, a.cond)
            if len(D) < MIN_STA_VALID:
                continue
            R = matched_filter(D, TPL, date, dev)
            if not len(R):
                continue
            dets.append(R)
            Rd = R[R.t.notna()]
            for r in Rd.itertuples():
                g = G[int(r.lid)]
                for s in g["long"]:
                    if s in TPL[int(r.lid)] or s in g["tpl"]:
                        seg = cut(D, s, r.t + g["off"][s], g["L"])
                        if seg is not None and (seg ** 2).sum() > 0:
                            k = (int(r.lid), s)
                            acc[k] = acc.get(k, 0) + unit_max3(seg)
                            cnt[k] = cnt.get(k, 0) + 1
                    seg = cut(D, s, r.t - LONG_PRE, LONGN)
                    if seg is not None and (seg ** 2).sum() > 0:
                        k = (int(r.lid), s)
                        lacc[k] = lacc.get(k, 0) + unit_max3(seg)
                        lcnt[k] = lcnt.get(k, 0) + 1
            A.log(f"[mf it{it}] {date} {len(D)} sta -> {len(Rd)} det")
        Dt = pd.concat(dets, ignore_index=True) if dets else pd.DataFrame(columns=["lid", "t", "cc", "nmad", "n_sta", "date"])
        Dt.to_parquet(dp, index=False)
        np.savez_compressed(lp, **{f"{l}_{s}": (lacc[(l, s)] / lcnt[(l, s)]).astype(np.float32)
                                   for (l, s) in lacc if lcnt[(l, s)] >= 2},
                            counts=np.array(json.dumps({f"{l}_{s}": int(c) for (l, s), c in lcnt.items()})))
        newT = {}
        for (l, s), v in acc.items():
            if cnt[(l, s)] >= 2 and s in G[l]["tpl"]:
                newT.setdefault(l, {})[s] = (v / cnt[(l, s)], G[l]["off"][s])
        # a location whose restack lost stations keeps its previous template stations
        for l in TPL:
            if l not in newT or len(newT[l]) < MIN_STA_VALID:
                newT[l] = TPL[l]
        TPL = newT
        json.dump({str(l): {s: v[0].tolist() for s, v in T.items()} for l, T in TPL.items()},
                  open(f"{out}/tpl_it{it + 1}.json", "w"))
        nd = Dt[Dt.t.notna()].groupby("lid").size()
        it_rows.append(dict(it=it, n_det=int(nd.sum()), n_loc_with_det=int(nd.size)))
        A.log(f"[mf it{it}] total {int(nd.sum())} detections, {nd.size} locations with >=1")

    # --- 4. diagnostics
    last = N_MF - 1
    Dt = pd.read_parquet(f"{out}/det_it{last}.parquet") if os.path.exists(f"{out}/det_it{last}.parquet") else pd.DataFrame()
    D0 = pd.read_parquet(f"{out}/det_it0.parquet") if os.path.exists(f"{out}/det_it0.parquet") else pd.DataFrame()
    LS = {}
    if os.path.exists(f"{out}/long_it{last}.npz"):
        z = np.load(f"{out}/long_it{last}.npz")
        for k in z.files:
            if k == "counts":
                continue
            l, s = k.split("_", 1)
            LS.setdefault(int(l), {})[s] = z[k].astype(np.float64)
    ndays_scanned = {}
    if len(Dt):
        dd = Dt[Dt.t.isna()]
        ndays_scanned = dd.groupby("lid").date.nunique().to_dict()
    rows = []
    for lid, g in G.items():
        row = dict(lid=lid)
        for nm, DD in (("it0", D0), ("fin", Dt)):
            if len(DD):
                x = DD[(DD.lid == lid) & DD.t.notna()]
                dm = DD[(DD.lid == lid) & DD.t.isna()]
                nday = max(dm.date.nunique(), 1)
                row[f"n_det_{nm}"] = len(x)
                row[f"det_per_day_{nm}"] = len(x) / nday
                row[f"nmad_med_{nm}"] = float(x.nmad.median()) if len(x) else np.nan
                row[f"cc_med_{nm}"] = float(x.cc.median()) if len(x) else np.nan
                row[f"day_max_nmad_med_{nm}"] = float(dm.day_max_nmad.median()) if len(dm) else np.nan
                row[f"n_days_{nm}"] = int(dm.date.nunique())
                row[f"n_days_with_det_{nm}"] = int(x.date.nunique())
        T = LS.get(lid, {})
        snr, dS, dP, tSo, tSp = [], [], [], [], []
        for s, st in T.items():
            from scipy.signal import hilbert
            aH = np.sqrt(np.abs(hilbert(st[1])) ** 2 + np.abs(hilbert(st[2])) ** 2)
            aZ = np.abs(hilbert(st[0]))
            k = int(round(0.4 * FS))
            aH = np.convolve(aH, np.ones(k) / k, "same")
            aZ = np.convolve(aZ, np.ones(k) / k, "same")
            tot = aH + aZ
            if np.median(tot) > 0:
                snr.append(float(tot.max() / np.median(tot)))
            if (st[1] ** 2).sum() > 0:
                iS = int(np.argmax(aH))
                tS_obs = iS / FS - LONG_PRE
                dS.append(tS_obs - g["tS"][s])
                tSo.append(tS_obs)
                tSp.append(g["tS"][s])
                lo = 0
                hi = max(iS - int(1.5 * FS), 1)
                iP = int(np.argmax(aZ[lo:hi]))
                dP.append(iP / FS - LONG_PRE - g["tP"][s])
        dS, dP = np.array(dS), np.array(dP)
        row.update(n_long_sta=len(T), snr_med=float(np.median(snr)) if snr else np.nan,
                   dS_med_abs=float(np.median(np.abs(dS))) if dS.size else np.nan,
                   frac_S_1s=float(np.mean(np.abs(dS) <= 1.0)) if dS.size else np.nan,
                   dP_med_abs=float(np.median(np.abs(dP))) if dP.size else np.nan,
                   frac_P_1s=float(np.mean(np.abs(dP) <= 1.0)) if dP.size else np.nan,
                   moveout_corr=float(np.corrcoef(tSo, tSp)[0, 1]) if len(tSo) >= 4 and np.std(tSo) > 0 else np.nan)
        rows.append(row)
    Dg = pd.DataFrame(rows)
    # Hyp2000 of the final long stacks (common time base: sample 0 = origin - LONG_PRE)
    hyp = pd.DataFrame()
    if LS:
        try:
            P = H.make_picks({l: T for l, T in LS.items() if len(T) >= 3})
            P.to_csv(f"{out}/picks.csv", index=False)
            stas_all = sorted({s for T in LS.values() for s in T})
            H.write_hyp_inputs(f"{out}/hyp", P, stas_all, 1.73, 1.5, model="GSC")
            hyp, logtxt = H.run_hyp(f"{out}/hyp")
            if len(hyp):
                hyp = hyp.rename(columns={"grp": "lid", "lat": "hyp_lat", "lon": "hyp_lon", "dep": "hyp_dep",
                                          "num": "hyp_nph", "rms": "hyp_rms", "erh": "hyp_erh", "erz": "hyp_erz",
                                          "gap": "hyp_gap"})
                hyp = hyp.drop(columns=["ot_s"])
            np_ = P.groupby("grp").agg(n_picks_P=("rP", lambda x: int((x >= 1.5).sum())),
                                       n_picks_S=("rS", lambda x: int((x >= 1.5).sum()))).reset_index().rename(columns={"grp": "lid"}) \
                if "rS" in P and "rP" in P else pd.DataFrame(columns=["lid"])
            Dg = Dg.merge(np_, on="lid", how="left")
        except Exception as e:
            A.log(f"[hyp] failed: {e!r}")
    Mr = pd.DataFrame([dict(key=k, lid=k // 10, seed_mode={v: m for m, v in MODES.items()}[k % 10]) for k in G])
    Dg = Dg.rename(columns={"lid": "key"})
    SR = pd.DataFrame(seed_rows).rename(columns={"lid": "key"})
    if len(hyp):
        hyp = hyp.rename(columns={"lid": "key"})
    out_t = Lc.merge(Sinfo, on="lid").merge(Mr, on="lid").merge(SR, on="key", how="left").merge(Dg, on="key", how="left")
    if len(hyp):
        out_t = out_t.merge(hyp, on="key", how="left")
        out_t["hyp_dist_km"] = A.epi_km(out_t.hyp_lat, out_t.hyp_lon, out_t.lat, out_t.lon)
        out_t["hyp_dz_km"] = out_t.hyp_dep - out_t.depth
    out_t["cond"] = a.cond
    out_t.to_csv(f"{out}/loc_table.csv", index=False)
    json.dump(dict(sub=a.sub, cond=a.cond, n_loc=len(G), days=len(days), iters=it_rows,
                   runtime_s=time.time() - t_start), open(f"{out}/summary.json", "w"), indent=1)
    A.log(f"[A4 {a.sub} {a.cond}] done in {(time.time() - t_start) / 60:.1f} min -> {out}/loc_table.csv")


if __name__ == "__main__":
    main()
