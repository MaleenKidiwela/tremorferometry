#!/usr/bin/env python
"""r23_cycle (2026-10-11) -- R23 study 4: ETS cycle regularity and 'bigger episode -> longer wait' 2003-2026 (B Q1/Q2 rules).
Rules: svi_science/results/r23/PLAN.md §5.    python svi_science/scripts/r23_cycle.py
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r23_common as C

OUT = f"{C.TB}/study4"; os.makedirs(OUT, exist_ok=True)
log = C.Log(f"{OUT}/study4.log")
NP = 10000


def perm_spear(x, y, rng):
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y); x, y = x[m], y[m]
    if len(x) < 4:
        return dict(n=int(len(x)), rho=np.nan, p_two=np.nan)
    rho = stats.spearmanr(x, y).statistic
    rx = stats.rankdata(x); ry = stats.rankdata(y)
    rx = (rx - rx.mean()) / np.sqrt(((rx - rx.mean()) ** 2).sum()); ry = (ry - ry.mean()) / np.sqrt(((ry - ry.mean()) ** 2).sum())
    P = np.argsort(rng.random((NP, len(ry))), axis=1)
    pr = ry[P] @ rx
    return dict(n=int(len(x)), rho=float(rho), p_two=float((np.sum(np.abs(pr) >= abs(rho)) + 1) / (NP + 1)),
                p_pos=float((np.sum(pr >= rho) + 1) / (NP + 1)))


def detrended_rank(size, year):
    r = stats.rankdata(size); b = np.polyfit(year, r, 1); return r - np.polyval(b, year)


def table(E, flag, size_col="S_rate"):
    X = E[E[flag]].sort_values("start").reset_index(drop=True).copy()
    X["onset"] = X.start_ts
    X["year"] = X.onset.dt.year + X.onset.dt.dayofyear / 365.25
    X["int_after"] = (X.onset.shift(-1) - X.onset).dt.days
    X["int_before"] = (X.onset - X.onset.shift(1)).dt.days
    X["mid_after"] = X.year + X.int_after / 730.5
    X["size"] = X[size_col]
    return X


def analyse(X, period, label, rng, size_label):
    a, b = C.PERIODS[period]
    m = (X.start >= a) & (X.start <= b)
    Y = X[m]
    iv = Y.int_after.dropna().to_numpy()
    out = dict(set=label, period=period, size=size_label, n_episodes=int(len(Y)), n_intervals=int(len(iv)))
    if len(iv) >= 3:
        cv = float(iv.std(ddof=1) / iv.mean())
        on = Y.onset.to_numpy().astype("datetime64[D]").astype(float)
        span0, span1 = on.min(), on.min() + iv.sum()
        ctrl = []
        for _ in range(NP):
            u = np.sort(rng.uniform(span0, span1, len(iv) + 1)); d = np.diff(u); ctrl.append(d.std(ddof=1) / d.mean())
        ctrl = np.array(ctrl)
        out.update(mean_d=float(iv.mean()), median_d=float(np.median(iv)), sd_d=float(iv.std(ddof=1)), min_d=float(iv.min()), max_d=float(iv.max()),
                   CV=cv, ctrl_CV_median=float(np.median(ctrl)), ctrl_CV_p01=float(np.percentile(ctrl, 1)), p_regular=C.pct_lo(cv, ctrl),
                   quasi_periodic=bool(cv < np.percentile(ctrl, 1)))
        d = perm_spear(Y.int_after, Y.mid_after, rng); out.update(drift_rho=d["rho"], drift_p=d["p_two"])
    tp = perm_spear(Y["size"], Y.int_after, rng)
    sp = perm_spear(Y["size"], Y.int_before, rng)
    ok = Y.int_after.notna() & Y["size"].notna()
    if ok.sum() >= 4:
        dtp = perm_spear(detrended_rank(Y["size"][ok], Y.year[ok]), Y.int_after[ok], rng)
    else:
        dtp = dict(rho=np.nan, p_two=np.nan)
    out.update(time_pred_rho=tp["rho"], time_pred_p=tp["p_two"], time_pred_n=tp["n"], time_pred_detr_rho=dtp["rho"], time_pred_detr_p=dtp["p_two"],
               slip_pred_rho=sp["rho"], slip_pred_p=sp["p_two"], slip_pred_n=sp["n"])
    if np.isfinite(tp["rho"]):
        sup = tp["p_two"] < 0.01 and np.isfinite(dtp["rho"]) and np.sign(dtp["rho"]) == np.sign(tp["rho"]) and dtp["p_two"] < 0.01
        out["time_pred_verdict"] = "supported" if sup and tp["rho"] > 0 else ("suggestive" if tp["p_two"] < 0.05 and tp["rho"] > 0 else "not supported")
    return out


def main():
    rng = np.random.default_rng(0)
    E = C.episodes()
    E["gnss_D"] = E.gnss_D_mm
    E["big_B"] = E.big_B.fillna(False).astype(bool)
    rows = []
    X = table(E, "big")
    X.to_csv(f"{OUT}/big_ets_2003_2026.csv", index=False)
    log("big ETS 2003-2026:\n" + X[["episode", "start", "end", "seg", "S_rate", "n_once", "n_bostock", "gnss_D", "int_before", "int_after"]].round(1).to_string(index=False))
    for period in C.PERIODS:
        rows.append(analyse(X, period, "big ETS (rate split)", rng, "S_rate"))
        rows.append(analyse(table(E, "big", "S_rate_kc"), period, "big ETS (rate split)", rng, "S_rate_kc"))
        rows.append(analyse(table(E, "major"), period, "all majors", rng, "S_rate"))
    Xb = X.copy(); Xb["size"] = Xb.n_bostock.where(Xb.n_bostock > 0)
    rows.append(analyse(Xb, "2003-13", "big ETS (rate split)", rng, "Bostock events"))
    Xg = X.copy(); Xg["size"] = Xg.gnss_D
    rows.append(analyse(Xg, "2014-26", "big ETS (rate split)", rng, "GNSS D (mm)"))
    # reproduction of B: B big list 2010+, GNSS D (2014+)
    XB = table(E, "big_B"); XB["size"] = XB.gnss_D
    rows.append(analyse(XB, "2014-26", "B big list (reproduction)", rng, "GNSS D (mm)"))
    XB2 = table(E, "big_B"); rows.append(analyse(XB2, "2014-26", "B big list (reproduction)", rng, "S_rate"))
    # 2003-2026 big list restricted to 2010+ (PNSN-era clock) for comparison with B
    R = pd.DataFrame(rows)
    R.to_csv(f"{OUT}/cycle_tests.csv", index=False)
    log("tests:\n" + R.round(4).to_string(index=False))
    figure(X, E)
    log("done")


def figure(X, E):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt, matplotlib.ticker
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.3))
    M = E[E.major]
    ax[0].scatter(M.start_ts, M.S_rate, s=14, color="#bdbdbd", label="other majors")
    ax[0].scatter(X.onset, X.S_rate, s=30, color="#08519c", label="big ETS")
    ax[0].set_yscale("log"); ax[0].set_ylabel("episode size, Bostock-equivalent LFEs (Σ rate)"); ax[0].legend(fontsize=8)
    for s in ("2006-07-01", "2009-08-27", "2014-01-01"):
        ax[0].axvline(pd.Timestamp(s), color="k", ls=":", lw=0.7)
    ax[0].set_title("Episode sizes 2003–2026 (our record)")
    ax[1].plot(X.onset[:-1], X.int_after[:-1] / 30.44, "o-", color="#08519c")
    ax[1].set_ylabel("months to the next big ETS"); ax[1].axhline(14, color="k", ls=":", lw=0.7); ax[1].set_title("Intervals between big ETS")
    ok = X.int_after.notna()
    c = np.where(X.start[ok] < "2014-01-01", "#08519c", "#fd8d3c")
    ax[2].scatter(X.S_rate[ok], X.int_after[ok] / 30.44, c=c, s=30)
    for _, r in X[ok].iterrows():
        ax[2].annotate(r.episode, (r.S_rate, r.int_after / 30.44), fontsize=6)
    ax[2].set_xscale("log"); ax[2].set_xlabel("size (Σ rate)"); ax[2].set_ylabel("months to the next big ETS")
    ax[2].set_title("Size vs wait (blue 2003–13, orange 2014–26)")
    ax[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    fig.tight_layout(); fig.savefig(f"{C.FIG}/r23_cycle.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
