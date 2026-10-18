#!/usr/bin/env python
"""bw_p6_onelook -- P6: the ONE pre-registered truth look per year (rules: results/backward/p6/P6.md).
Writes p6/ONELOOK_OPENED.json FIRST and refuses to run again; refuses if the backward models changed after freezing.
Scored: frozen B0 2003, B1 2005-09, B2 2007, B3 2008 + 2009; era 2003 + 2005-09; reported as ALREADY SEEN: 2004 test blocks
(frozen B0, era).  Nothing is fitted or tuned on what is read here.
Out: p6/tables/onelook.csv, onelook_family.csv, onelook.json
"""
import os, sys, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_common as G
R, C = G.patch_forward()
import svi_ref_fit as F

BW = "/home/jovyan/tremorferometry/svi_science/results/backward"
P6 = f"{BW}/p6"
EP = pd.read_csv(f"{G.SVI}/results/backward_scope/truth_episodes.csv")
TCORE = {"CN4": [8.79, 9.82], "B1": [8.69, 9.56], "B2": [8.75], "B3": [8.75], "ERA": [8.25]}
LOOK = [(2003, "frozen", "CN4", "p5f_2003", False), (2003, "era", "ERA", "p5e_2003_kelb", False),
        (2005, "frozen", "B1", "p5f_2005", False), (2005, "era", "ERA", "p5e_2005", False),
        (2007, "frozen", "B2", "p5f_2007", False), (2008, "frozen", "B3", "p5f_2008", False),
        (2009, "frozen", "B3", "p5f_2009", False),
        (2004, "frozen", "CN4", "p5f_2004", True), (2004, "era", "ERA", "p5e_2004", True)]


def boot_recall(hit, ev_days, nb=1000, seed=4):
    g = pd.DataFrame(dict(h=hit.astype(float), d=ev_days)).groupby("d").h.agg(["sum", "count"])
    s, n = g["sum"].to_numpy(), g["count"].to_numpy(); rng = np.random.default_rng(seed); out = []
    for _ in range(nb):
        i = rng.integers(0, len(s), len(s)); out.append(s[i].sum() / n[i].sum())
    return np.percentile(out, [2.5, 97.5])


def main():
    dry = "--dry-seen" in sys.argv                                  # mechanics test on the already-seen 2004 rows only; no marker
    mk = f"{P6}/ONELOOK_OPENED.json"
    if not dry and os.path.exists(mk):
        raise SystemExit(f"[look] {mk} exists: the one look was already opened; refusing")
    info = json.load(open(f"{P6}/models/MODELS_FROZEN_BW.json"))
    if R.sha256(f"{P6}/models/bw_models.joblib") != info["sha256"]:
        raise SystemExit("[look] backward models changed after freezing; refusing")
    import subprocess
    pt = subprocess.run(["date"], env=dict(os.environ, TZ="America/Los_Angeles"), capture_output=True, text=True).stdout.strip()
    if not dry: json.dump(dict(opened_utc=pd.Timestamp.now("UTC").isoformat(), opened_pt=pt, models_sha256=info["sha256"],
                   forward_sha256=info["forward_sha256"], years=[2003, 2005, 2007, 2008, 2009], already_seen=[2004],
                   rule="one look per year; nothing fitted or tuned on it (P6.md)"), open(mk, "w"), indent=1)
    B = G.bostock_all()
    kfb = pd.read_csv(f"{P6}/tables/k_factors_bw.csv"); kff = pd.read_csv(f"{G.REFINE}/rates/k_factors.csv")
    setk = {s: float(g.set_mean_k.iloc[0]) for s, g in kfb.groupby("sset")}
    setk.update({{"S14": "B3", "S15": "S15"}[s]: float(g.set_mean_k.iloc[0]) for s, g in kff.groupby("sset") if s in ("S14", "S15")})
    rows, frows = [], []
    for year, det, ns, run, seen in LOOK:
        if dry and not seen: continue
        e = EP[EP.year == year].iloc[0]
        days = [str(x.date()) for x in pd.date_range(e.start, e.end)]
        A = pd.read_parquet(f"{P6}/P/{run}_P.parquet"); A = A[A.date.isin(days)]
        Bw = B[B.date.isin(days)]
        if seen:                                                  # 2004: TEST blocks only (build blocks = era fit rows)
            A = A[R.block_role(A.t.to_numpy()) == "fit"]; Bw = Bw[R.block_role(Bw.t_origin.to_numpy()) == "fit"]
        D = A[A.series == "data"]; Cc = A[A.series == "c1"]
        ms = (R.nearest_dt(D.t.to_numpy(), D.fam.to_numpy(), Bw.t_origin.to_numpy(), Bw.fam.to_numpy(), True) <= R.MATCH_S).astype(float)
        P = D.P.to_numpy()
        sl, ic, _ = F.calib(ms, P)
        r = dict(year=year, detector=det, net_set=ns, run=run, already_seen=seen, scope="test blocks" if seen else "all blocks",
                 n_days=len(days), n_events=len(Bw), n_data=len(D), n_ctl=len(Cc),
                 auc=F.auc(ms, P), auc_net_mad=F.auc(ms, D.net_mad), slope=sl, matched_over_sumP=float(ms.sum() / P.sum()),
                 ctl_frac_P_gt05=float((Cc.P > 0.5).mean()) if len(Cc) else np.nan,
                 n_data_P075=int((P >= 0.75).sum()), n_ctl_P075=int((Cc.P >= 0.75).sum()))
        sel = D[P >= 0.75]
        hit = R.one2one_recall(sel.t.to_numpy(), sel.fam.to_numpy(), Bw.t_origin.to_numpy(), Bw.fam.to_numpy())
        r["recall_P075"] = float(hit.mean()); lo, hi = boot_recall(hit, Bw.date.to_numpy()); r["recall_P075_lo95"], r["recall_P075_hi95"] = lo, hi
        r["precision_P075"] = float(ms[P >= 0.75].mean()) if (P >= 0.75).any() else np.nan
        fn = pd.Series(hit).groupby(Bw.fam.to_numpy()).sum(); r["fam_rediscovered_P075"] = int((fn >= 5).sum())
        r["fam_with_events"] = int(Bw.fam.nunique())
        for i, T in enumerate(TCORE[ns]):
            s2 = D[D.net_mad >= T]
            h2 = R.one2one_recall(s2.t.to_numpy(), s2.fam.to_numpy(), Bw.t_origin.to_numpy(), Bw.fam.to_numpy())
            sfx = "" if i == 0 else "_Thi"
            r[f"T_core{sfx}"] = T; r[f"recall_Tcore{sfx}"] = float(h2.mean())
            r[f"precision_Tcore{sfx}"] = float(ms[D.net_mad.to_numpy() >= T].mean()) if (D.net_mad >= T).any() else np.nan
        exc = P.sum() - Cc.P.sum()
        r["sumP_data"] = float(P.sum()); r["sumP_ctl"] = float(Cc.P.sum()); r["his_over_sumP_excess"] = float(len(Bw) / exc) if exc > 0 else np.nan
        r["his_over_sumP_data"] = float(len(Bw) / P.sum()); r["set_k_fitdays"] = setk.get(ns, np.nan)
        r["bar_slope"] = bool(np.isfinite(sl) and 0.8 <= sl <= 1.2); r["bar_matched_over_sumP"] = bool(0.8 <= r["matched_over_sumP"] <= 1.2)
        r["bar_ctl_le_1pct"] = bool(r["ctl_frac_P_gt05"] <= 0.01); r["all_bars"] = r["bar_slope"] and r["bar_matched_over_sumP"] and r["bar_ctl_le_1pct"]
        rows.append(r); print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()}, flush=True)
        ev = Bw.groupby("fam").size()
        for f_, n_ in ev.items():
            frows.append(dict(year=year, detector=det, net_set=ns, already_seen=seen, fam=int(f_), n_events=int(n_), matched_P075=int(fn.get(f_, 0))))
    T = pd.DataFrame(rows)
    if dry:
        print("[look] DRY (2004 already-seen rows only): nothing written"); return
    T.to_csv(f"{P6}/tables/onelook.csv", index=False)
    pd.DataFrame(frows).to_csv(f"{P6}/tables/onelook_family.csv", index=False)
    json.dump(T.to_dict("records"), open(f"{P6}/tables/onelook.json", "w"), indent=1, default=float)
    pd.set_option("display.width", 300)
    print(T[["year", "detector", "net_set", "already_seen", "n_events", "recall_P075", "precision_P075", "recall_Tcore", "auc", "slope",
             "matched_over_sumP", "ctl_frac_P_gt05", "his_over_sumP_excess", "set_k_fitdays", "all_bars"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
