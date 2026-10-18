#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6 Amendment 4b): the Amendment 2b synthetics with the MEASURED LFE envelope broadening.
Envelope t^2 exp(-t/tau) on a 3 Hz carrier; tau set per row from the measured FWHM (P and S separately):
  W0 constant tau 0.25 s | W1 region trend w(d) (headline) | W2 own row's measured width | W3 3x trend slope |
  W2s scrambled control (W2 widths permuted across rows of the region).

  python s11_eq_synth_width.py [--workers 4] [--ndraw 3]
writes tomo/eq/tables/synthw_rows.csv, synthw_kappa.csv, synthw_closure.csv, synthw_tau_lookup.csv
"""
import argparse
import os
import sys
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.signal import hilbert  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import s11_eq_pickertest as PT  # noqa: E402
import s11_eq_synth as SY  # noqa: E402
import s11_eq_swidth as SW  # noqa: E402
import s11_eq_diff as D  # noqa: E402
import s11_tomo_pooled as PO  # noqa: E402
import s11_tomo_profile as PR  # noqa: E402

FS, PRE = SY.FS, SY.PRE
N = SY.N
WLEN = 6.0
TAUS = np.round(np.arange(0.04, 0.80001, 0.01), 3)


def wavelet(rng, tau, f0=3.0):
    t = np.arange(int(WLEN * FS)) / FS
    w = t ** 2 * np.exp(-t / tau) * np.cos(2 * np.pi * f0 * t + rng.uniform(0, 2 * np.pi))
    return w / np.abs(w).max()


def fwhm_of(w):
    x = np.zeros(int(20 * FS)); k0 = int(5 * FS)
    x[k0:k0 + len(w)] = w
    env = np.abs(hilbert(x))
    return SW.width(env, int(np.argmax(env)), search=0.05)[0]


def tau_lookup(rng):
    rows = []
    for tau in TAUS:
        rows.append(dict(tau=tau, fwhm=float(np.median([fwhm_of(wavelet(rng, tau)) for _ in range(40)]))))
    T = pd.DataFrame(rows)
    T["fwhm"] = np.maximum.accumulate(T.fwhm.values)     # monotone for inversion
    return T


def tau_for(T, w):
    w = np.clip(w, T.fwhm.min(), T.fwhm.max())
    return float(np.interp(w, T.fwhm.values, T.tau.values))


def synth(tp, ts, rng, tauP, tauS, leak=0.3):
    st = np.zeros((3, N))
    wp, ws = wavelet(rng, tauP), wavelet(rng, tauS)
    az = rng.uniform(0, 2 * np.pi); baz = rng.uniform(0, 2 * np.pi)
    kp, ks = (PRE + tp) * FS, (PRE + ts) * FS
    SY.place(st[0], kp, wp, 1.0)
    SY.place(st[1], kp, wp, 0.3 * np.cos(baz)); SY.place(st[2], kp, wp, 0.3 * np.sin(baz))
    SY.place(st[1], ks, ws, 4.0 * np.cos(az)); SY.place(st[2], ks, ws, 4.0 * np.sin(az))
    SY.place(st[0], ks, ws, 4.0 * leak)
    return st


def swidth_rec(st, tS):
    eh = np.sqrt(np.abs(hilbert(st[1])) ** 2 + np.abs(hilbert(st[2])) ** 2)
    return SW.width(eh, (tS + PRE) * FS)[0] if np.isfinite(tS) else np.nan


COMBOS = [("W0", 1.75), ("W1", 1.75), ("W2", 1.75), ("W3", 1.75), ("W2s", 1.75), ("W0", 1.90), ("W1", 1.90)]


def do_station(args):
    reg, sta, Q, ndraw, seed = args
    rng = np.random.default_rng(seed)
    bank = SY.noise_bank(sta, seed=seed)
    rows = []
    for r in Q.itertuples():
        for var, kap in COMBOS:
            tauP, tauS = getattr(r, f"tauP_{var}"), getattr(r, f"tauS_{var}")
            ts = kap * r.tp
            for d in range(ndraw + 1):
                st = synth(r.tp, ts, rng, tauP, tauS)
                if d > 0:
                    if not bank:
                        continue
                    st, _ = PT.noise_match(st, bank[int(rng.integers(len(bank)))], r.snr)
                pr = PT.pick(st)
                pr = {k: (v - PRE if k.startswith("t") else v) for k, v in pr.items()}
                rows.append(dict(reg=reg, fam=r.fam, sta=sta, epi=r.epi, variant=var, kappa_true=kap, draw=d,
                                 tauP=tauP, tauS=tauS, tp_true=r.tp, ts_true=ts,
                                 wS_meas=swidth_rec(st, pr.get("tS", np.nan)), **pr))
    return rows


def build(a):
    rng = np.random.default_rng(5)
    T = tau_lookup(rng)
    T.to_csv(f"{PT.EQ}/tables/synthw_tau_lookup.csv", index=False)
    print("tau lookup fwhm range", T.fwhm.min(), T.fwhm.max(), "tau0.25 ->", float(np.interp(0.25, T.tau, T.fwhm)), flush=True)
    W = pd.read_csv(f"{PT.EQ}/tables/swidth_rows.csv")
    W["fam"] = W.fam.astype(str)
    xy = PR.coords()
    jobs, trend = [], []
    for reg, (ps, prof) in SY.REG.items():
        w = W[W.reg == reg]
        ws = w[w.rS >= 1.2].dropna(subset=["wS", "epi_km"])
        bS = np.polyfit(ws.epi_km, ws.wS, 1)
        wp = w[(w.rP >= 1.2)].dropna(subset=["wP", "epi_km"])
        bP = np.polyfit(wp.epi_km, wp.wP, 1)
        trend.append(dict(reg=reg, S_slope_ms_km=1000 * bS[0], S_int=bS[1], P_slope_ms_km=1000 * bP[0], P_int=bP[1]))
        L = pd.read_csv(f"{PT.TT}/picks_{ps}.csv", low_memory=False)
        L = L[(L.depth_km >= 25) & (L.depth_km <= 50) & L.snr.notna() & L.sta.isin(list(xy))].copy()
        L["epi"] = [np.hypot((r.lon - xy[r.sta][1]) * 111.19 * np.cos(np.radians(r.lat)), (r.lat - xy[r.sta][0]) * 111.19)
                    for r in L.itertuples()]
        L = L[L.epi <= 100]
        pf = D.load_profile(prof)
        L["tp"] = np.nan
        for z, g in L.groupby(L.depth_km.round(0)):
            L.loc[g.index, "tp"] = D.ray_partials(pf, z, g.epi.values).sum(1)
        L["fam"] = L.fam.astype(str)
        L = L.drop(columns=[c for c in ("wS", "wP") if c in L]).merge(
            w[["fam", "sta", "wS", "wP", "rP"]].rename(columns={"rP": "rP_w"}).drop_duplicates(["fam", "sta"]),
            on=["fam", "sta"], how="left")
        L.loc[~(L.rP_w >= 1.2), "wP"] = np.nan          # own P width only where P was picked (rP >= 1.2)
        trS, trP = np.polyval(bS, L.epi), np.polyval(bP, L.epi)
        mS, mP = L.epi.mean(), L.epi.mean()
        tr3S = np.polyval(bS, mS) + 3 * bS[0] * (L.epi - mS)
        tr3P = np.polyval(bP, mP) + 3 * bP[0] * (L.epi - mP)
        own_S = np.where(np.isfinite(L.wS), L.wS, trS); own_P = np.where(np.isfinite(L.wP), L.wP, trP)
        perm = rng.permutation(len(L))
        for nm, vS, vP in (("W0", None, None), ("W1", trS, trP), ("W2", own_S, own_P), ("W3", tr3S, tr3P),
                           ("W2s", own_S[perm], own_P[perm])):
            if nm == "W0":
                L["tauS_W0"] = 0.25; L["tauP_W0"] = 0.25
            else:
                L[f"tauS_{nm}"] = [tau_for(T, v) for v in np.asarray(vS)]
                L[f"tauP_{nm}"] = [tau_for(T, v) for v in np.asarray(vP)]
        L["own_S_measured"] = np.isfinite(L.wS)
        cols = ["fam", "tp", "snr", "epi"] + [c for c in L.columns if c.startswith("tau")]
        for i, (sta, g) in enumerate(L.groupby("sta")):
            jobs.append((reg, sta, g[cols], a.ndraw, 300 + i + (0 if reg == "WA_NW" else 50)))
        print(reg, "rows", len(L), "own widths", int(L.own_S_measured.sum()), trend[-1], flush=True)
    pd.DataFrame(trend).to_csv(f"{PT.EQ}/tables/synthw_trend.csv", index=False)
    print("jobs", len(jobs), "rows", sum(len(j[2]) for j in jobs), flush=True)
    with Pool(a.workers) as pool:
        res = pool.map(do_station, jobs, chunksize=1)
    return pd.DataFrame([r for rr in res for r in rr])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--ndraw", type=int, default=3)
    a = ap.parse_args()
    f = f"{PT.EQ}/tables/synthw_rows.csv"
    if not os.path.exists(f):
        build(a).to_csv(f, index=False)
    R = pd.read_csv(f)
    R["fam"] = R.fam.astype(str)
    rng = np.random.default_rng(19)
    PO.NB, PO.NS, PO.EST, PO.NSTA_MIN = 300, 100, "inv", 4
    out, clo = [], []
    for (reg, var, kap, d), g in R.groupby(["reg", "variant", "kappa_true", "draw"]):
        for pk, tp, rp in (("V3", "tPv3", "rPv3"), ("frozen", "tP0", "rP0")):
            q = g[(g[rp] >= 1.2) & (g.rS >= 1.2)]
            q = pd.DataFrame(dict(fam=q.fam, sta=q.sta, tP=q[tp], tS=q.tS)).dropna()
            q = q[q.tS > q.tP].groupby("fam").filter(lambda x: len(x) >= 4)
            row = dict(reg=reg, variant=var, kappa_true=kap, draw=d, picker=pk, n_fam=q.fam.nunique())
            if q.fam.nunique() >= 5:
                o, _ = PO.zone(q, rng, "synthw", True)
                row.update({k: o[k] for k in ("n_rows", "kappa", "sig", "rms", "loso_max", "d_site", "pass")})
                row["bias"] = row["kappa"] - kap
            gg = g.dropna(subset=[tp])
            row.update(med_dtP=float((gg[tp] - gg.tp_true).median()), med_dtS=float((g.tS - g.ts_true).median()))
            # delay growth with distance (s/km) of the picks, and of S minus P
            for nm, dd, yv in (("dtS", g.dropna(subset=["tS"]), None), ("dtP", gg, None)):
                y = (dd.tS - dd.ts_true) if nm == "dtS" else (dd[tp] - dd.tp_true)
                ok = y.abs() < 1.0
                row[f"{nm}_slope_ms_km"] = 1000 * np.polyfit(dd.epi[ok], y[ok], 1)[0] if ok.sum() > 20 else np.nan
            out.append(row)
            print({k: (round(v, 3) if isinstance(v, float) else v) for k, v in row.items()}, flush=True)
        if d > 0:
            m = g.rS >= 1.2
            b = np.polyfit(g.epi[m & g.wS_meas.notna()], g.wS_meas[m & g.wS_meas.notna()], 1)
            clo.append(dict(reg=reg, variant=var, kappa_true=kap, draw=d, synth_wS_slope_ms_km=1000 * b[0],
                            synth_wS_med=float(g.wS_meas[m].median())))
    pd.DataFrame(out).to_csv(f"{PT.EQ}/tables/synthw_kappa.csv", index=False)
    pd.DataFrame(clo).to_csv(f"{PT.EQ}/tables/synthw_closure.csv", index=False)
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
