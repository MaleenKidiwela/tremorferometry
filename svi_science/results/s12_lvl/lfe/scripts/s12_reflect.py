#!/usr/bin/env python
"""s12_reflect -- S12-LVL (LFE side), task 2: search the family grand stacks for a later phase reflected or converted at the
LVL base (or an oceanic Moho), using its predicted moveout across families and stations (PLAN.md §2).

Data: S11 grand stacks G1 (P-weighted, 3C, S-aligned, 40 Hz, 2-8 Hz, tau -14..+10 s) of every causal (family, station,
epoch) series (store + store_ext); epochs of one (family, station) are combined (n-weighted, aligned on the picked S).
Observable per pair: residual log-envelope r(t) = ln E(t) - smooth coda trend (quadratic in t over 0.6-6.5 s after S),
E = horizontal (H) or vertical (Z) energy envelope; and the S-wavelet matched filter |cc|(t) (H components).
Moveout stack: for each layer model (H km, Vs km/s) every pair's r(t) is shifted by its predicted lag (flat-layer ray
theory, s12_forward) and averaged; the peak = mean in +-0.1 s minus the flanks 0.3-0.6 s. Scrambled control: predicted lags
permuted among the pairs of the same station (keeps every station's lag distribution and its own reverberations), 1000
draws, max over the model grid (scan-corrected). Fixed-lag stacks per station show station reverberations.
    python s12_reflect.py [--region SVI]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
import sys, json, time, glob, argparse
import numpy as np
import pandas as pd
from scipy.signal import hilbert
ROOT = "/home/jovyan/tremorferometry"
sys.path.insert(0, f"{ROOT}/svi_science/scripts"); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s11_common as C
import s12_forward as FW

OUT = FW.OUT; TB, FIG = FW.TB, FW.FIG
FS = 40.0
TREL = np.arange(-12.0, 9.0 + 1e-9, 1 / FS)          # time after S
HG = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0]
VG = [2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8]
PHASES = {"SrS_b": ("H", "S"), "SrS_m": ("H", "S"), "SrP_b": ("Z", "S"), "PrP_b": ("Z", "P")}
NSCR = 1000
TAG = ""            # output tag ("_lb" for the 1-3 Hz restack)
SMN = 5              # envelope smoothing (samples)
CEN, FL0, FL1 = 0.1, 0.3, 0.6   # peak metric: centre half-width and flank range (s)
rng = np.random.default_rng(20261013)


def log(m):
    s = f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {m}"
    print(s, flush=True); open(f"{FW.LOG}/reflect.log", "a").write(s + "\n")


def load_pairs():
    ser = pd.concat([pd.read_csv(f"{C.TB}/series.csv"), pd.read_csv(f"{C.TB}/series_ext.csv")], ignore_index=True)
    ser["key"] = ser.fam.astype(str) + "|" + ser.netsta + "|" + ser.epoch
    ser = ser[(ser.n_data > 0) & (ser.causal25_3c > 1.5)].set_index("key")
    acc = {}
    for d in (f"{C.TB}/store", f"{C.TB}/store_ext"):
        for p in sorted(glob.glob(f"{d}/*.npz")):
            z = np.load(p); tau = z["tau"]
            for key, G, tS, sp in zip(z["keys"], z["G1"], z["tS"], z["sp"]):
                key = str(key)
                if key not in ser.index:
                    continue
                fam, ns, ep = key.split("|")
                X = np.array([np.interp(TREL + tS, tau, np.nan_to_num(G[c].astype(float)), left=np.nan, right=np.nan)
                              if np.isfinite(G[c]).all() else np.full(len(TREL), np.nan) for c in range(3)])
                n = float(ser.loc[key, "n_data"])
                k = (int(fam), ns)
                a = acc.setdefault(k, dict(X=np.zeros((3, len(TREL))), w=np.zeros(3), n=0.0, sp=[], eps=[]))
                for c in range(3):
                    if np.isfinite(X[c]).all():
                        x = X[c] / np.sqrt(np.mean(X[c] ** 2))
                        a["X"][c] += n * x; a["w"][c] += n
                a["n"] += n; a["sp"].append(float(sp)); a["eps"].append(ep)
    P = []
    for (fam, ns), a in acc.items():
        X = np.where(a["w"][:, None] > 0, a["X"] / np.maximum(a["w"][:, None], 1e-12), np.nan)
        P.append(dict(fam=fam, netsta=ns, X=X, n=a["n"], sp=float(np.median(a["sp"])), epochs=",".join(a["eps"])))
    return P


def envelopes(X):
    e = np.full((3, X.shape[1]), np.nan)
    for c in range(3):
        if np.isfinite(X[c]).all():
            e[c] = np.abs(hilbert(X[c])) ** 2
    EZ = e[0]
    EH = np.nansum(e[1:], 0) if np.isfinite(e[1:]).all(1).any() else np.full(X.shape[1], np.nan)
    return EZ, EH


def smooth(x, n=5):
    k = np.ones(n) / n
    return np.convolve(x, k, mode="same")


def resid_logenv(E, t0=0.6, t1=6.5):
    """ln E normalised by the S peak, minus a quadratic trend fitted over t0..t1 s after S -> residual on TREL"""
    if not np.isfinite(E).all():
        return np.full(len(TREL), np.nan)
    Es = smooth(E, SMN)
    pk = Es[(TREL >= -0.3) & (TREL <= 0.5)].max()
    L = np.log(np.maximum(Es / pk, 1e-8))
    m = (TREL >= t0) & (TREL <= t1)
    c = np.polyfit(TREL[m], L[m], 2)
    r = L - np.polyval(c, TREL)
    r[TREL < t0 - 0.2] = np.nan
    return r


def resid_logenv_P(E, tp, t0=0.25, t1=30.0):
    """same about P (tp = time of P after S, negative): trend over P+t0..min(P+t1, S-0.5)"""
    # 2026-10-14 fix: 10 of the 380 1-3 Hz pairs have no S-P (NaN in offsets_sp.csv) -> empty P window crashed the scan;
    # such pairs now get no P-side residual (the S-side phases are unaffected)
    if not np.isfinite(E).all() or not np.isfinite(tp) or not ((TREL >= tp - 0.3) & (TREL <= tp + 0.4)).any():
        return np.full(len(TREL), np.nan)
    Es = smooth(E, SMN)
    pk = Es[(TREL >= tp - 0.3) & (TREL <= tp + 0.4)].max()
    L = np.log(np.maximum(Es / pk, 1e-8))
    m = (TREL >= tp + t0) & (TREL <= min(tp + t1, -0.5))
    if m.sum() < 20:
        return np.full(len(TREL), np.nan)
    c = np.polyfit(TREL[m], L[m], 2)
    r = L - np.polyval(c, TREL)
    r[(TREL < tp + t0 - 0.1) | (TREL > min(tp + t1, -0.5))] = np.nan
    return r


def mf_cc(X):
    """S-wavelet matched filter on the horizontals: |cc|(t) of the S pulse (-0.2..+0.4 s) with the trace, mean of comps"""
    out = []
    w = (TREL >= -0.2) & (TREL <= 0.4); nw = w.sum()
    for c in (1, 2):
        if not np.isfinite(X[c]).all():
            continue
        s = X[c][w] - X[c][w].mean(); s /= np.linalg.norm(s)
        x = X[c]
        cc = np.full(len(TREL), np.nan)
        i0 = np.flatnonzero(w)[0]
        for i in range(len(TREL) - nw):
            seg = x[i:i + nw] - x[i:i + nw].mean(); n_ = np.linalg.norm(seg)
            if n_ > 0:
                cc[i - i0 + i0 + nw // 2] = abs(seg @ s) / n_
        out.append(cc)
    if not out:
        return np.full(len(TREL), np.nan)
    r = np.nanmean(out, 0)
    r[TREL < 0.5] = np.nan
    return r


def shift_stack(R, lags, w, half=1.2):
    """R (n, T) residual traces, lags (n,) -> weighted mean of R(t + lag) on a window +-half s"""
    tt = np.arange(-half, half + 1e-9, 1 / FS)
    S = np.zeros(len(tt)); W = np.zeros(len(tt))
    for r, L, wi in zip(R, lags, w):
        if not np.isfinite(L):
            continue
        v = np.interp(tt + L, TREL, np.nan_to_num(r, nan=np.nan), left=np.nan, right=np.nan)
        ok = np.isfinite(v)
        S[ok] += wi * v[ok]; W[ok] += wi
    return tt, np.where(W > 0, S / np.maximum(W, 1e-12), np.nan)


def peak_metric(tt, s):
    c = np.abs(tt) <= CEN; f = (np.abs(tt) >= FL0) & (np.abs(tt) <= FL1)
    return float(np.nanmean(s[c]) - np.nanmean(s[f]))


def pair_metric(R):
    """per pair: boxcar mean over +-0.1 s minus the mean over the flanks 0.3-0.6 s, as a function of the centre time"""
    n = int(round(CEN * FS)); M = np.full(R.shape, np.nan)
    cs = lambda x, a, b: x[a:b]
    for i, r in enumerate(R):
        x = r.copy()
        ok = np.isfinite(x)
        c0 = np.convolve(np.nan_to_num(x), np.ones(2 * n + 1), "same"); k0 = np.convolve(ok.astype(float), np.ones(2 * n + 1), "same")
        cen = np.where(k0 >= 2 * n + 1, c0 / np.maximum(k0, 1), np.nan)
        fl = np.zeros(len(x)); kf = np.zeros(len(x))
        for sgn in (-1, 1):
            lo, hi = int(round(FL0 * FS)), int(round(FL1 * FS))
            for d in range(lo, hi + 1):
                sh = np.roll(np.nan_to_num(x), -sgn * d); so = np.roll(ok.astype(float), -sgn * d)
                if sgn == -1:
                    sh[:d] = 0; so[:d] = 0
                else:
                    sh[-d:] = 0; so[-d:] = 0
                fl += sh; kf += so
        flm = np.where(kf >= 0.9 * 2 * (hi - lo + 1), fl / np.maximum(kf, 1), np.nan)
        M[i] = cen - flm
    return M


def predicted_lags(P):
    """flat-layer lags for every pair, phase and model -> dict phase -> (n, nH, nV); lag after S (or P for PrP)"""
    geo = json.load(open(f"{TB}/fwd_geometry.json"))
    L = {ph: np.full((len(P), len(HG), len(VG)), np.nan) for ph in PHASES}
    tl = {ph: np.full((len(P), len(HG), len(VG)), np.nan) for ph in PHASES}
    for i, p in enumerate(P):
        src = np.array(geo["src"][str(p["fam"])]); rec = np.array(geo["sta"][p["netsta"]])
        x = float(np.hypot(*(src[:2] - rec[:2]))); zs, zr = float(src[2]), float(rec[2])
        p["epi_km"] = x; p["zs"] = zs
        for a, H in enumerate(HG):
            for b, V in enumerate(VG):
                Ls = FW.model_layers(zs, H, V)
                base = {}
                for kind in ("S", "P"):
                    segs, _ = FW.phase_segs(Ls, zs, zr, kind)
                    base[kind] = FW.shoot_fast(segs, x)[0][1]
                for ph in PHASES:
                    segs, li = FW.phase_segs(Ls, zs, zr, ph)
                    (X_, T, Ti), _ = FW.shoot_fast(segs, x)
                    L[ph][i, a, b] = T - base[PHASES[ph][1]]
                    tl[ph][i, a, b] = sum(Ti[j] for j in li)
                if a == HG.index(3.5) and b == VG.index(3.0):
                    p["sp_pred"] = base["S"] - base["P"]
    return L, tl


LOADER = None


def main():
    global NSCR
    ap = argparse.ArgumentParser(); ap.add_argument("--region", default="SVI"); ap.add_argument("--nscr", type=int, default=1000)
    a = ap.parse_args(); NSCR = a.nscr
    t0 = time.time()
    P = LOADER()
    log(f"{len(P)} (family, station) pairs from causal series")
    L, TL = predicted_lags(P)
    log(f"lags predicted ({time.time() - t0:.0f}s)")
    st = np.array([p["netsta"] for p in P]); w = np.sqrt(np.array([p["n"] for p in P]))
    fam_arr = np.array([p["fam"] for p in P])
    RH, RZ, MF, RP = [], [], [], []
    for p in P:
        EZ, EH = envelopes(p["X"])
        RH.append(resid_logenv(EH)); RZ.append(resid_logenv(EZ)); MF.append(mf_cc(p["X"]))
        RP.append(resid_logenv_P(EZ, -p["sp"]))
    RH, RZ, MF, RP = map(np.array, (RH, RZ, MF, RP))
    np.savez_compressed(f"{TB}/refl_traces{TAG}.npz", RH=RH, RZ=RZ, MF=MF, RP=RP, TREL=TREL,
                        fam=np.array([p["fam"] for p in P]), netsta=st, n=np.array([p["n"] for p in P]),
                        epi=np.array([p["epi_km"] for p in P]), zs=np.array([p["zs"] for p in P]),
                        sp=np.array([p["sp"] for p in P]), sp_pred=np.array([p.get("sp_pred", np.nan) for p in P]),
                        X=np.array([p["X"] for p in P]).astype(np.float32),
                        **{f"lag_{ph}": L[ph] for ph in PHASES}, **{f"tlvl_{ph}": TL[ph] for ph in PHASES})
    obs = {"SrS_b": ("RH", RH), "SrS_m": ("RH", RH), "SrP_b": ("RP", RP), "PrP_b": ("RP", RP), "SrS_b_mf": ("MF", MF)}
    sp_arr = np.array([p["sp"] for p in P])
    rows = []; scr_rows = []
    groups = {s: np.flatnonzero(st == s) for s in np.unique(st)}
    nH, nV = len(HG), len(VG)
    for ph, (nm, R) in obs.items():
        phl = ph.replace("_mf", "")
        M = pair_metric(R)                                       # (n, T): central +-0.1 s minus flanks 0.3-0.6 s
        # lag relative to the anchor phase (S, or P for the P-side phases) and the anchor's place on the TREL axis
        sp_pred = np.array([p.get("sp_pred", np.nan) for p in P])
        if phl == "PrP_b":
            Lrel = L[phl].copy(); anc = -sp_arr
        elif phl == "SrP_b":
            Lrel = L[phl] + sp_pred[:, None, None]; anc = -sp_arr          # SrP time after P
        else:
            Lrel = L[phl].copy(); anc = np.zeros(len(P))
        if phl in ("PrP_b", "SrP_b"):
            Lrel = np.where(Lrel >= 0.6, Lrel, np.nan)                       # too close to the P pulse: not usable
        def grid_val(perm):
            Lp = Lrel[perm] + anc[:, None, None]
            ix = np.clip(np.round((np.nan_to_num(Lp, nan=-99) - TREL[0]) * FS).astype(int), 0, len(TREL) - 1)
            v = M[np.arange(len(P))[:, None, None], ix]
            v = np.where(np.isfinite(Lp), v, np.nan)
            ok = np.isfinite(v); ww = w[:, None, None] * ok
            return np.where(ok.sum(0) >= 10, np.nansum(ww * np.nan_to_num(v), 0) / np.maximum(ww.sum(0), 1e-30), np.nan)
        ix = np.clip(np.round((np.nan_to_num(Lrel + anc[:, None, None], nan=-99) - TREL[0]) * FS).astype(int), 0, len(TREL) - 1)
        grid = grid_val(np.arange(len(P)))
        smax = np.zeros(NSCR); snom = np.zeros(NSCR)
        for k in range(NSCR):
            perm = np.arange(len(P))
            for s_, idx in groups.items():
                perm[idx] = rng.permutation(idx)
            g2 = grid_val(perm)
            smax[k] = np.nanmax(g2); snom[k] = g2[HG.index(3.5), VG.index(3.0)]
        ib = np.unravel_index(np.nanargmax(grid), grid.shape)
        ino = (HG.index(3.5), VG.index(3.0))
        r = dict(phase=ph, obs=nm, best_H=HG[ib[0]], best_V=VG[ib[1]], best_peak=float(grid[ib]),
                 p_scan=float((np.sum(smax >= grid[ib]) + 1) / (NSCR + 1)), nom_peak=float(grid[ino]),
                 p_nom=float((np.sum(snom >= grid[ino]) + 1) / (NSCR + 1)),
                 scr_max_q95=float(np.quantile(smax, .95)), scr_nom_q95=float(np.quantile(snom, .95)))
        # odd / even families at the best model
        for par, msk in (("odd", fam_arr % 2 == 1), ("even", fam_arr % 2 == 0)):
            v = M[np.flatnonzero(msk), ix[np.flatnonzero(msk), ib[0], ib[1]]]
            v = np.where(np.isfinite(Lrel[np.flatnonzero(msk), ib[0], ib[1]]), v, np.nan)
            ok = np.isfinite(v); r[f"peak_{par}"] = float(np.sum(w[msk][ok] * v[ok]) / np.sum(w[msk][ok]))
        rows.append(r); log(f"{ph}: {r}")
        np.save(f"{TB}/refl_grid{TAG}_{ph}.npy", grid); np.save(f"{TB}/refl_scr{TAG}_{ph}.npy", np.c_[smax, snom])
        for s_, idx in groups.items():
            for lab, (a_, b_) in (("best", ib), ("nom", ino)):
                v = np.where(np.isfinite(Lrel[idx, a_, b_]), M[idx, ix[idx, a_, b_]], np.nan); ok = np.isfinite(v)
                # station's own scrambled control: lags permuted within the station only
                sc = []
                for k in range(200):
                    pj = rng.permutation(idx)
                    Lq = Lrel[pj, a_, b_] + anc[idx]
                    jx = np.clip(np.round((np.nan_to_num(Lq, nan=-99) - TREL[0]) * FS).astype(int), 0, len(TREL) - 1)
                    vv = np.where(np.isfinite(Lq), M[idx, jx], np.nan); o2 = np.isfinite(vv)
                    sc.append(np.sum(w[idx][o2] * vv[o2]) / max(np.sum(w[idx][o2]), 1e-30))
                pk = float(np.sum(w[idx][ok] * v[ok]) / max(np.sum(w[idx][ok]), 1e-30))
                scr_rows.append(dict(phase=ph, netsta=s_, site=C.SITE.get(s_, s_), model=lab, n_pairs=len(idx), peak=pk,
                                     scr_q95=float(np.quantile(sc, .95)), p=float((np.sum(np.array(sc) >= pk) + 1) / 201)))
    pd.DataFrame(rows).to_csv(f"{TB}/refl_scan{TAG}.csv", index=False)
    pd.DataFrame(scr_rows).to_csv(f"{TB}/refl_station{TAG}.csv", index=False)
    open(f"{OUT}/REFLECT{TAG.upper()}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))
    log(f"done {time.time() - t0:.0f}s")


if __name__ == "__main__":
    LOADER = load_pairs
    main()
