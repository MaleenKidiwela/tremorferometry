#!/usr/bin/env python
"""SH_KMRFIX task 1 (2026-10-01, CPU): NC.KMR during 2025-06 .. 2025-08-25 (old STS-2 failing), from the cut windows.
(a) daily per-channel raw health (std, distinct values, DC level) from every cut window;
(b) were KMR channels masked by the cleaning? (results/sh_ext2026/prep_summary.csv dead_ch / ref);
(c) per-channel signal: direct stacks of the strongest S1 9x events per family in each period vs the old it0 template,
    with a scrambled control (each event's slice taken at its own random time in the window);
(d) per-candidate per-channel KMR cc at the aligned position for every >= 8x cut candidate 2025-06-01..2025-08-25,
    validated against the stored channel-mean cc_NC.KMR (results/sh_conf/candidates_conf.parquet).
Outputs results/sh_kmrfix/t1_*.  Nothing is overwritten."""
import os, sys, json, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import sh_checks_kmr as K
RES = K.RES; CUT = K.CUT; EXT = K.EXT; OUT = f"{RES}/sh_kmrfix"; os.makedirs(OUT, exist_ok=True)
T_A = pd.Timestamp("2025-06-01"); T_END = pd.Timestamp("2025-08-26")
PER = [("healthy", "2025-06-01", "2025-07-25"), ("drift", "2025-07-25", "2025-08-03"), ("HHE_stuck", "2025-08-03", "2025-08-26")]


def per_of(t):
    for n, a, b in PER:
        if pd.Timestamp(a) <= t < pd.Timestamp(b): return n
    return None


def main():
    OFF = K.offsets()
    M = pd.read_csv(f"{RES}/sh_v3b/tpl_it0.csv"); M["fam"] = M.fam.astype(str); Z0 = dict(np.load(f"{RES}/sh_v3b/tpl_it0.npz"))
    T0 = {f: Z0[f"Shelly2025|{f}|NC.KMR"] for f in [str(i) for i in range(1, 28)]}
    # (b) cleaning flags
    P = pd.read_csv(f"{EXT}/prep_summary.csv")
    P = P[(P.netsta == "NC.KMR") & (P.day >= "2025-06-01") & (P.day <= "2025-09-17")][["day", "chans", "dead_plain", "dead_ch", "station_dead", "ref", "invalid_frac", "muted_frac"]]
    P.to_csv(f"{OUT}/t1_kmr_prep_flags.csv", index=False)
    # candidates: every >= 8x cut candidate in the window
    D = pd.read_parquet(f"{CUT}/detections.parquet"); D["fam"] = D.fam.astype(str)
    D = D[(D.t >= T_A) & (D.t < T_END)].reset_index(drop=True)
    D["ta"] = D.t_epoch + D.fam.map(OFF) - K.PRE
    D["month"] = pd.to_datetime(D.ta, unit="s").dt.strftime("%Y-%m")
    n = int((K.PRE + K.POST) * K.FS0)
    W = np.zeros((len(D), 3, K.NW), np.float32); ok = np.zeros(len(D), bool)
    hrows = []
    for m, g in D.groupby("month"):
        r = K.load_month_kmr(m)
        if r is None: continue
        segs, dat = r
        for i, ta in zip(g.index, g.ta.values):
            xs = [K.read_span(segs, dat, c, ta, n) for c in K.CH]
            if any(x is None for x in xs): continue
            hrows.append(dict(i=i, **{f"std_{c}": float(np.std(x)) for c, x in zip(K.CH, xs)},
                              **{f"nuniq_{c}": int(len(np.unique(x))) for c, x in zip(K.CH, xs)},
                              **{f"mean_{c}": float(np.mean(x)) for c, x in zip(K.CH, xs)}))
            W[i] = np.stack([K.proc(x) if np.ptp(x) > 0 else np.zeros(K.NW) for x in xs]); ok[i] = True
        print(f"[t1] {m}: {len(g)} candidates, windows {ok[g.index].sum()}", flush=True)
    H = pd.DataFrame(hrows).set_index("i")
    D = D.join(H)
    D["ok"] = ok
    # (a) daily health
    D["day"] = D.t.dt.strftime("%Y-%m-%d")
    agg = {}
    for c in K.CH:
        agg[f"std_med_{c}"] = (f"std_{c}", "median"); agg[f"frac_const_{c}"] = (f"nuniq_{c}", lambda v: float((v <= 2).mean()))
        agg[f"dc_med_{c}"] = (f"mean_{c}", "median")
    Hd = D[D.ok].groupby("day").agg(n=("ok", "size"), **agg).reset_index()
    Hd.to_csv(f"{OUT}/t1_kmr_daily_health.csv", index=False)
    # (d) per-candidate per-channel cc at lag 0 (template slice S-1 s .. S+5 s)
    S0, L = K.S0, K.L
    cc = np.full((len(D), 3), np.nan)
    for i in np.flatnonzero(ok):
        T = T0[D.fam.iat[i]]
        for c in range(3):
            x = W[i, c, S0:S0 + L].astype(float)
            if np.ptp(x) > 0: cc[i, c] = K.ncc(T[c], x)
    for c in range(3): D[f"cc_{K.CH[c]}"] = cc[:, c]
    # validation vs stored cc_NC.KMR (channel-mean over live channels) at the same candidates
    A = pd.read_parquet(f"{RES}/sh_conf/candidates_conf.parquet", columns=["run", "series", "fam", "t_epoch", "cc_NC.KMR", "nlive", "ccsum"])
    A = A[(A.run == "fwd") & (A.series == "data")]; A["fam"] = A.fam.astype(str)
    A["key"] = np.round(A.t_epoch * 20).astype(np.int64); D["key"] = np.round(D.t_epoch * 20).astype(np.int64)
    D = D.merge(A[["fam", "key", "cc_NC.KMR"]], on=["fam", "key"], how="left")
    D["period"] = [per_of(t) for t in D.t]
    D.drop(columns=["ta"]).to_parquet(f"{OUT}/t1_kmr_candidate_cc.parquet", index=False)
    vrows = []
    for p, g in D[D.ok].groupby("period"):
        for name, cols in [("mean3", ["cc_HHE", "cc_HHN", "cc_HHZ"]), ("mean_NZ", ["cc_HHN", "cc_HHZ"]), ("Z", ["cc_HHZ"]), ("N", ["cc_HHN"]), ("E", ["cc_HHE"])]:
            v = g[cols].mean(1); s = g["cc_NC.KMR"]; m = v.notna() & s.notna()
            vrows.append(dict(period=p, recon=name, n=int(m.sum()), r=float(np.corrcoef(v[m], s[m])[0, 1]) if m.sum() > 10 else np.nan,
                              med_abs_diff=float(np.median(np.abs(v[m] - s[m]))) if m.sum() else np.nan))
    V = pd.DataFrame(vrows); V.to_csv(f"{OUT}/t1_cc_reconstruction_check.csv", index=False); print(V.to_string())
    # (c) per-channel signal: direct stacks of the 100 strongest per family (by ccsum among S1-like candidates) per period
    rng = np.random.default_rng(7); srows = []
    for p, _, _ in PER:
        for f in [str(i) for i in range(1, 28)]:
            g = D[(D.period == p) & (D.fam == f) & D.ok].sort_values("ccsum", ascending=False).head(100)
            if len(g) < 10: continue
            Ws = W[g.index.values].astype(np.float64)
            r = dict(period=p, fam=f, n=len(g))
            for c in range(3):
                X = Ws[:, c, :]
                if not np.any(X): r[f"stack_cc_{K.CH[c]}"] = np.nan; r[f"ctl_cc_{K.CH[c]}"] = np.nan; continue
                Xn = X / (np.sqrt((X[:, S0:S0 + L] ** 2).mean(1, keepdims=True)) + 1e-30)
                st = Xn.mean(0)
                r[f"stack_cc_{K.CH[c]}"] = float(max(K.ncc(T0[f][c], st[S0 + l:S0 + l + L]) for l in range(-5, 6)))
                a = rng.integers(20, K.NW - L - 20, len(X))
                Xr = np.stack([X[k, a[k] - 5:a[k] + L + 5] for k in range(len(X))])
                Xr = Xr / (np.sqrt((Xr[:, 5:5 + L] ** 2).mean(1, keepdims=True)) + 1e-30)
                sr = Xr.mean(0)
                r[f"ctl_cc_{K.CH[c]}"] = float(max(K.ncc(T0[f][c], sr[5 + l:5 + l + L]) for l in range(-5, 6)))
                r[f"single_cc_med_{K.CH[c]}"] = float(np.nanmedian(g[f"cc_{K.CH[c]}"]))
            srows.append(r)
    S = pd.DataFrame(srows); S.to_csv(f"{OUT}/t1_kmr_channel_stacks.csv", index=False)
    Sm = S.groupby("period").median(numeric_only=True).reindex([p for p, _, _ in PER])
    Sm["n_fam"] = S.groupby("period").size()
    Sm.to_csv(f"{OUT}/t1_kmr_channel_stacks_summary.csv"); print(Sm.round(3).to_string())


if __name__ == "__main__":
    main()
