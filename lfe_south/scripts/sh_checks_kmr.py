#!/usr/bin/env python
"""SH_CHECKS (2026-10-01) CHECK 1: why do NC.KMR's re-stacked templates (Trillium Horizon 120 from 2025-09-17,
Trillium Compact Posthole from 2026-03-12) match the old STS-2 it0 templates at only cc ~0.13?

CPU only, reads the cut windows (results/cut/sh_ext) -- raw continuous data is gone.  Nothing is overwritten.
  python scripts/sh_checks_kmr.py extract   -> results/sh_checks/kmr_windows.npz (20 sps, 45 s, 3 ch, per event)
  python scripts/sh_checks_kmr.py analyze   -> results/sh_checks/kmr_*.csv + figures
"""
import os, sys, glob, json, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from scipy.signal import butter, sosfiltfilt, resample_poly
warnings.filterwarnings("ignore")
WS = "/home/jovyan/tremorferometry/lfe_south"; RES = f"{WS}/results"; CUT = f"{RES}/cut/sh_ext"; EXT = f"{RES}/sh_ext2026"
OUT = f"{RES}/sh_checks"; os.makedirs(OUT, exist_ok=True)
SW1 = pd.Timestamp("2025-09-17 23:18").timestamp(); SW2 = pd.Timestamp("2026-03-12 20:43").timestamp()
CH = ["HHE", "HHN", "HHZ"]; FS0 = 100.0; FS = 20.0
PRE, POST = 15.0, 30.0; NW = int((PRE + POST) * FS)        # 900 samples, arrival at index 300
IA = int(PRE * FS); S0 = IA - 20; L = 120                   # template slice = arrival-1 s .. +5 s
SOS_BP = butter(4, [1.0 / 50, 8.0 / 50], btype="band", output="sos"); SOS_HP = butter(4, 2.0 / 10.0, btype="highpass", output="sos")
NMAX = 300


def offsets():
    O = pd.read_csv(f"{CUT}/offsets_used.csv"); O["fam"] = O.fam.astype(str)
    return O[O.netsta == "NC.KMR"].set_index("fam").offset_s.to_dict()


def event_sets():
    """per set: DataFrame fam, t_epoch, rank_key (higher = stronger)."""
    CAT = pd.read_parquet(f"{EXT}/cat_S1_9x_internal.parquet"); CAT["fam"] = CAT.fam.astype(str)
    CAT["t_epoch"] = CAT.t.astype("datetime64[ns]").astype(np.int64) / 1e9
    old = CAT[CAT.seg.isin(["A", "cont"]) & (CAT.t_epoch < SW1 - 86400)].copy()
    mid = pd.Timestamp("2024-11-01").timestamp()
    S = {"OLD_all": old, "OLD_h1": old[old.t_epoch < mid], "OLD_h2": old[old.t_epoch >= mid]}
    for ep in ("B", "C"):
        N = pd.read_parquet(f"{EXT}/cat_NK_{ep}.parquet"); N["fam"] = N.fam.astype(str)
        S[f"NK_{ep}"] = N
    # confident detections (P_same) from the candidate table (B/C segments use the new KMR templates; A uses it0)
    D = pd.read_parquet(f"{CUT}/detections.parquet"); D["fam"] = D.fam.astype(str)
    C = pd.read_parquet(f"{CUT}/conf_P_by_det_id.parquet")[["det_id", "P_same", "P_same_NO_PNET", "in_sample"]]
    D = D.merge(C, on="det_id", how="left")
    S["P_OLD"] = D[D.seg.isin(["A", "cont"]) & (D.t_epoch < SW1 - 86400) & (~D.in_sample.fillna(True).astype(bool))].assign(ccsum=lambda x: x.P_same)
    S["P_B"] = D[(D.seg == "B")].assign(ccsum=lambda x: x.P_same)
    S["P_C"] = D[(D.seg == "C")].assign(ccsum=lambda x: x.P_same)
    out = {}
    for k, E in S.items():
        E = E[np.isfinite(E.ccsum)].sort_values("ccsum", ascending=False).groupby("fam").head(NMAX)
        out[k] = E[["fam", "t_epoch", "ccsum"]].reset_index(drop=True)
    return out


def load_month_kmr(m):
    import h5py
    p = f"{CUT}/segments_{m}.parquet"
    if not os.path.exists(p): return None
    S = pd.read_parquet(p); S = S[S.seed_id.str.startswith("NC.KMR..")]
    if not len(S): return None
    dat = {}
    with h5py.File(f"{CUT}/windows_{m}.h5", "r") as h:
        for c in CH:
            k = f"data/NC.KMR..{c}"
            dat[c] = h[k][:] if k in h else None
    segs = {c: S[S.seed_id == f"NC.KMR..{c}"].sort_values("t0").reset_index(drop=True) for c in CH}
    return segs, dat


def read_span(segs, dat, c, ta, n):
    """n raw samples at 100 sps starting at epoch ta from the union segments; None if not fully covered / gap."""
    G = segs[c]
    if G is None or not len(G) or dat[c] is None: return None
    j = np.searchsorted(G.t0.values, ta + 0.005) - 1
    if j < 0: return None
    g = G.iloc[j]
    i0 = int(round((ta - g.t0) * FS0))
    if i0 < 0 or i0 + n > g.n: return None
    for a, b in json.loads(g.gaps):
        if a < i0 + n and b > i0: return None
    return dat[c][int(g.h5_start) + i0:int(g.h5_start) + i0 + n]


def proc(x):
    y = sosfiltfilt(SOS_BP, x.astype(float) - x.mean())
    y = resample_poly(y, 1, 5)
    return sosfiltfilt(SOS_HP, y)


def extract():
    OFF = offsets(); SETS = event_sets()
    allev = pd.concat([E.assign(set=k) for k, E in SETS.items()], ignore_index=True)
    allev["ta"] = allev.t_epoch + allev.fam.map(OFF) - PRE
    allev["month"] = pd.to_datetime(allev.ta, unit="s").dt.strftime("%Y-%m")
    n = int((PRE + POST) * FS0)
    W = np.zeros((len(allev), 3, NW), np.float32); ok = np.zeros(len(allev), bool)
    for m, g in allev.groupby("month"):
        r = load_month_kmr(m)
        if r is None: continue
        segs, dat = r; cache = {}
        for i, ta in zip(g.index, g.ta.values):
            key = round(ta, 2)
            if key not in cache:
                xs = [read_span(segs, dat, c, ta, n) for c in CH]
                cache[key] = None if any(x is None for x in xs) or any(np.max(np.abs(x)) == 0 for x in xs) else np.stack([proc(x) for x in xs])
            if cache[key] is not None: W[i] = cache[key]; ok[i] = True
        print(f"[extract] {m}: {len(g)} events, ok {ok[g.index].sum()}", flush=True)
    allev["ok"] = ok
    np.savez_compressed(f"{OUT}/kmr_windows.npz", W=W)
    allev.to_parquet(f"{OUT}/kmr_windows_events.parquet")
    print(allev.groupby("set").ok.agg(["size", "sum"]))


# ================================================================================================ analysis
from sh_v3_restack import s1_stack
SENS = {"OLD": 629145000.0, "B": 504678000.0, "C": 316573000.0}
LAG = 40                                       # +-2 s at 20 sps


def ncc(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).sum() / (np.sqrt((a * a).sum() * (b * b).sum()) + 1e-30))


def tcc(a, b, lag=5):
    out = []
    for c in range(min(a.shape[0], b.shape[0])):
        x, y = a[c].astype(float), b[c].astype(float)
        if not np.any(x) or not np.any(y): continue
        best = -1
        for s in range(-lag, lag + 1):
            u = x[max(0, s):len(x) + min(0, s)]; v = y[max(0, -s):len(y) + min(0, -s)]
            best = max(best, ncc(u, v))
        out.append(best)
    return float(np.mean(out)) if out else np.nan


def lagcurve(T, X, lags=range(-LAG, LAG + 1), ti=None, xi=None):
    """cc(T[ti], X[xi, S0+l : S0+l+L]) for each lag (X = 900-sample stack)."""
    return np.array([ncc(T[ti], X[xi, S0 + l:S0 + l + L]) for l in lags])


def s1_template(Wset, slice_start=None):
    """Text S1 recipe on (n, 3, 900) windows ranked strongest first; slice at S0 or per-event random starts."""
    n = len(Wset); T = np.zeros((3, L)); info = []
    for c in range(3):
        if slice_start is None:
            X = Wset[:, c, :]; s0 = S0
        else:   # move each event's random slice to a common position
            X = np.zeros((n, L + 40), np.float32)
            for i, a in enumerate(slice_start): X[i] = Wset[i, c, a - 20:a + L + 20]
            s0 = 20
        ok = np.ones(n, bool)
        r = s1_stack(X, ok, s0, L, 2, 100, 50)
        if r is None: return None, None
        T[c], nn, cck, ccm = r; info.append(cck)
    return T, float(np.mean(info))


def direct_stack(Wset):
    X = Wset / (np.sqrt((Wset[:, :, S0:S0 + L] ** 2).mean(2, keepdims=True)) + 1e-30)
    return X.mean(0)


def rot(X, th):
    """rotate horizontals (rows 0=E, 1=N) of X by th degrees (component azimuths rotated clockwise)."""
    t = np.radians(th); E, N = X[0], X[1]
    Y = X.copy(); Y[0] = E * np.cos(t) - N * np.sin(t); Y[1] = E * np.sin(t) + N * np.cos(t)
    return Y


def analyze():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    W = np.load(f"{OUT}/kmr_windows.npz")["W"]; E = pd.read_parquet(f"{OUT}/kmr_windows_events.parquet")
    M = pd.read_csv(f"{RES}/sh_v3b/tpl_it0.csv"); M["fam"] = M.fam.astype(str); Z0 = dict(np.load(f"{RES}/sh_v3b/tpl_it0.npz"))
    ZB = dict(np.load(f"{EXT}/tpl_kmrB.npz")); ZC = dict(np.load(f"{EXT}/tpl_kmrC.npz"))
    rng = np.random.default_rng(1)
    rows = []; lagc = {}; ex = {}
    sets = ["OLD_all", "OLD_h1", "OLD_h2", "P_OLD", "NK_B", "P_B", "NK_C", "P_C"]
    epoch_of = {"OLD_all": "OLD", "OLD_h1": "OLD", "OLD_h2": "OLD", "P_OLD": "OLD", "NK_B": "B", "P_B": "B", "NK_C": "C", "P_C": "C"}
    for f in [str(i) for i in range(1, 28)]:
        key = f"Shelly2025|{f}|NC.KMR"; T0 = Z0[key]
        for sname in sets:
            g = E[(E.set == sname) & (E.fam == f) & E.ok].sort_values("ccsum", ascending=False)
            idx = g.index.values[:100]
            if len(idx) < 10: continue
            Ws = W[idx].astype(np.float64)
            ep = epoch_of[sname]
            # noise level (m/s, 2-8 Hz RMS) and single-event SNR
            noise = np.sqrt((Ws[:, :, 20:200] ** 2).mean(2)); sig = np.sqrt((Ws[:, :, S0:S0 + L] ** 2).mean(2))
            row = dict(fam=f, set=sname, epoch=ep, n=len(idx),
                       noise_nms=float(np.median(noise) / SENS[ep] * 1e9), snr_single_med=float(np.median(sig / noise)))
            Ds = direct_stack(Ws)
            row["dstack_snr"] = float(np.sqrt((Ds[:, S0:S0 + L] ** 2).mean()) / np.sqrt((Ds[:, 20:200] ** 2).mean()))
            row["dstack_cc_old_tcc"] = tcc(Ds[:, S0:S0 + L], T0)
            # lag scan, per channel and channel-mean
            lc = np.array([lagcurve(T0, Ds, ti=c, xi=c) for c in range(3)]); lagc[(f, sname)] = lc
            m = lc.mean(0); j = int(np.argmax(m)); row.update(lag_best_s=(j - LAG) / FS, cc_at_best_lag=float(m[j]), cc_at_0=float(m[LAG]),
                                                            cc_min_over_lags=float(m.min()))
            for c in range(3): row[f"cc0_{CH[c]}"] = float(lc[c, LAG])
            # channel permutation: max |cc| over lags between T0[i] and stack[j]
            P = np.array([[np.max(np.abs(lagcurve(T0, Ds, ti=i, xi=k))) for k in range(3)] for i in range(3)])
            row["perm_diag_mean"] = float(np.trace(P) / 3); row["perm_offdiag_max"] = float((P - np.diag(np.diag(P))).max())
            row["perm_best_is_identity"] = bool(all(P[i].argmax() == i for i in range(3)))
            # polarity: signed cc at the best |cc| lag per channel
            for c in range(3):
                l = lc[c]; jj = int(np.argmax(np.abs(l))); row[f"signed_best_{CH[c]}"] = float(l[jj])
            # rotation of horizontals (5 deg) -- tcc of horizontals, jointly
            ths = np.arange(0, 360, 5)
            rc = [np.mean([tcc(rot(Ds, th)[c:c + 1, S0:S0 + L], T0[c:c + 1]) for c in (0, 1)]) for th in ths]
            row["rot_best_deg"] = int(ths[int(np.argmax(rc))]); row["rot_best_cc_h"] = float(np.max(rc)); row["rot_0_cc_h"] = float(rc[0])
            # Text S1 recipe template from these windows
            Ts, cc50 = s1_template(Ws)
            row["s1_cc50"] = cc50; row["s1_cc_old"] = tcc(Ts, T0)
            if ep == "B": row["s1_cc_vs_tplB"] = tcc(Ts, ZB[key]) if key in ZB else np.nan
            if ep == "C": row["s1_cc_vs_tplC"] = tcc(Ts, ZC[key]) if key in ZC else np.nan
            row["dstack_cc_tplB"] = tcc(Ds[:, S0:S0 + L], ZB[key]) if key in ZB else np.nan
            row["dstack_cc_tplC"] = tcc(Ds[:, S0:S0 + L], ZC[key]) if key in ZC else np.nan
            # scrambled control: same windows, each event's slice at its own random time (5-14 s before / 8-23 s after S)
            for rep in range(3):
                cand = np.r_[np.arange(22, IA - 100), np.arange(IA + 160, NW - L - 22)]
                st = rng.choice(cand, size=len(idx))
                Tr, ccr = s1_template(Ws, st)
                row[f"scr_s1_cc50_{rep}"] = ccr; row[f"scr_s1_cc_old_{rep}"] = tcc(Tr, T0)
                Xr = np.stack([Ws[i, :, a - S0:a - S0 + NW] if False else Ws[i, :, a:a + L] for i, a in enumerate(st)])
                Dr = (Xr / (np.sqrt((Xr ** 2).mean(2, keepdims=True)) + 1e-30)).mean(0)
                row[f"scr_dstack_cc_old_{rep}"] = tcc(Dr, T0)
            if f in ("1", "8", "19", "23") : ex[(f, sname)] = (Ds, Ts)
            rows.append(row)
        # template-to-template
        rows.append(dict(fam=f, set="TPL", epoch="-", tplB_cc_old=tcc(ZB[key], T0) if key in ZB else np.nan,
                         tplC_cc_old=tcc(ZC[key], T0) if key in ZC else np.nan, tplB_cc_tplC=tcc(ZB[key], ZC[key]) if key in ZB and key in ZC else np.nan))
        print(f"[analyze] fam {f} done", flush=True)
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/kmr_per_family_set.csv", index=False)
    for c in ("scr_s1_cc50", "scr_s1_cc_old", "scr_dstack_cc_old"):
        R[c] = R[[f"{c}_{k}" for k in range(3)]].mean(1)
    cols = ["n", "noise_nms", "snr_single_med", "dstack_snr", "dstack_cc_old_tcc", "cc_at_0", "cc_at_best_lag", "lag_best_s", "cc0_HHE", "cc0_HHN", "cc0_HHZ",
            "perm_diag_mean", "perm_offdiag_max", "perm_best_is_identity", "rot_best_deg", "rot_best_cc_h", "rot_0_cc_h", "s1_cc50", "s1_cc_old", "s1_cc_vs_tplB", "s1_cc_vs_tplC",
            "dstack_cc_tplB", "dstack_cc_tplC", "scr_s1_cc50", "scr_s1_cc_old", "scr_dstack_cc_old"]
    S = R[R.set != "TPL"].groupby("set")[cols].median(numeric_only=False).reindex(sets)
    S.to_csv(f"{OUT}/kmr_summary_by_set.csv"); pd.set_option("display.width", 300); print(S.round(3).T.to_string())
    print(R[R.set == "TPL"][["tplB_cc_old", "tplC_cc_old", "tplB_cc_tplC"]].median())
    np.savez_compressed(f"{OUT}/kmr_lagcurves.npz", **{f"{k[0]}|{k[1]}": v for k, v in lagc.items()})
    # figures: lag curves (median over families) and example stacks
    f_, ax = plt.subplots(1, 2, figsize=(12, 4))
    lags = np.arange(-LAG, LAG + 1) / FS
    for sname, col in zip(sets, ["k", "0.4", "0.6", "#1f77b4", "#d62728", "#ff9896", "#2ca02c", "#98df8a"]):
        A = [lagc[(f, sname)].mean(0) for f in [str(i) for i in range(1, 28)] if (f, sname) in lagc]
        if A: ax[0].plot(lags, np.median(A, 0), color=col, label=f"{sname} (n fam {len(A)})")
    ax[0].set_xlabel("lag of the stack vs the old STS-2 template (s)"); ax[0].set_ylabel("cc (median over families, channel mean)"); ax[0].legend(fontsize=7); ax[0].axvline(0, color="0.8")
    for sname, col in zip(sets, ["k", "0.4", "0.6", "#1f77b4", "#d62728", "#ff9896", "#2ca02c", "#98df8a"]):
        if sname in S.index: ax[1].scatter([sname], [S.loc[sname, "dstack_cc_old_tcc"]], color=col, s=30)
    ax[1].scatter(sets, S.scr_dstack_cc_old.values, marker="x", color="0.5", label="scrambled control (random slices)")
    ax[1].set_ylabel("direct stack (top-100) cc vs old template"); ax[1].tick_params(axis="x", rotation=45); ax[1].legend(fontsize=7)
    f_.tight_layout(); f_.savefig(f"{OUT}/kmr_lag_and_cc.png", dpi=130)
    f_, ax = plt.subplots(len(ex) // 4 if len(ex) >= 4 else 1, 1, figsize=(1, 1)); plt.close(f_)
    fams_ex = sorted({k[0] for k in ex}, key=int)
    f_, ax = plt.subplots(len(fams_ex), 3, figsize=(13, 2.4 * len(fams_ex)), squeeze=False)
    t = (np.arange(L) - 20) / FS
    for i, fa in enumerate(fams_ex):
        key = f"Shelly2025|{fa}|NC.KMR"; T0 = Z0[key]
        for c in range(3):
            a = ax[i, c]; nz = lambda x: x / (np.abs(x).max() + 1e-30)
            a.plot(t, nz(T0[c]) + 3, color="k", lw=0.9)
            for k, (sn, colr) in enumerate([("OLD_all", "0.4"), ("NK_B", "#d62728"), ("NK_C", "#2ca02c")]):
                if (fa, sn) in ex: a.plot(t, nz(ex[(fa, sn)][0][c, S0:S0 + L]) + 1.5 - 1.5 * k, color=colr, lw=0.8)
            if key in ZB: a.plot(t, nz(ZB[key][c]) - 3, color="#d62728", lw=0.8, ls="--")
            if key in ZC: a.plot(t, nz(ZC[key][c]) - 4.5, color="#2ca02c", lw=0.8, ls="--")
            a.set_yticks([]); a.set_title(f"fam {fa} {CH[c]}", fontsize=8)
    ax[0, 0].text(-1, 3.5, "old it0 template", fontsize=7); ax[0, 0].text(-1, 2.0, "direct stacks: OLD, NK_B (red), NK_C (green)", fontsize=7)
    ax[0, 0].text(-1, -2.5, "new templates B / C (dashed)", fontsize=7); ax[-1, 1].set_xlabel("time after S (s)")
    f_.tight_layout(); f_.savefig(f"{OUT}/kmr_example_stacks.png", dpi=120)


if __name__ == "__main__":
    {"extract": extract, "analyze": analyze}[sys.argv[1]]()
