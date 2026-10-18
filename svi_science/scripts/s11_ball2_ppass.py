#!/usr/bin/env python
"""s11_ball2_ppass -- S11 ballistic round 2, §2 "better P" pass rule (round2/PLAN.md).
Per station x bin type x P estimator: median P cc of non-railed stacks (data, scrambled control), rail fraction,
median bootstrap sigma; split-half cc of the grand stacks; rule: P3 median cc >= 0.6 on EP, YR or Y2 for >= 3 stations,
control median < 0.3.

    python s11_ball2_ppass.py [--tag _test]
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = "/home/jovyan/tremorferometry/svi_science/results/s11_depth/ballistic/round2"
TB, FIG = f"{OUT}/tables", f"{OUT}/figures"
STATIONS = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011", "PB.B926", "CN.VGZ", "CN.SNB"]
PK = ["P3", "P3j", "Pz", "Pzr"]
BINS = ["EP", "YR", "Y2"]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default=""); a = ap.parse_args(); tag = a.tag
    st = pd.read_parquet(f"{TB}/lags_station{tag}.parquet"); geo = pd.read_csv(f"{TB}/series_geom{tag}.csv")
    rows = []
    for (ns, bk, kind), g in st.groupby(["netsta", "bin", "kind"]):
        for k in ["S"] + PK:
            ok = ~g[f"rail{k}"]
            rows.append(dict(netsta=ns, bin=bk, kind=kind, est=k, n_stacks=len(g), rail_frac=float(1 - ok.mean()),
                             cc_med=float(g.loc[ok, f"cc{k}"].median()) if ok.any() else np.nan,
                             cc_med_all=float(g[f"cc{k}"].median()),
                             sig_med=float(g.loc[ok, f"sig{k}"].median()) if ok.any() else np.nan,
                             sigSm_med=float(g.loc[ok, f"sigSm{k}"].median()) if (k != "S" and f"sigSm{k}" in g and ok.any()) else np.nan,
                             n_rows_med=float(g.n.median())))
    Q = pd.DataFrame(rows); Q.to_csv(f"{TB}/P_quality{tag}.csv", index=False)
    # rule 2
    verdict = {}
    for k in PK:
        d = Q[(Q.kind == "data") & (Q.est == k)].pivot(index="netsta", columns="bin", values="cc_med")
        c = Q[(Q.kind == "ctrl") & (Q.est == k)].pivot(index="netsta", columns="bin", values="cc_med")
        best = d.max(axis=1); passing = sorted(best[best >= 0.6].index)
        ctrl_max = float(np.nanmax(c.to_numpy())) if c.size else np.nan
        verdict[k] = dict(stations_ge_0p6=passing, n_pass=len(passing), ctrl_max_median=ctrl_max,
                          PASS=bool(len(passing) >= 3 and ctrl_max < 0.3), per_station_best=best.round(3).to_dict())
    H = geo.groupby("netsta")[["half_ccS"] + [f"half_cc{k}" for k in PK] + ["inc_deg", "SmP_pick"]].median()
    H.to_csv(f"{TB}/P_splithalf{tag}.csv")
    json.dump(verdict, open(f"{TB}/P_verdict{tag}.json", "w"), indent=1)
    with pd.option_context("display.width", 250, "display.max_rows", 300):
        print(Q[Q.kind == "data"].pivot_table(index=["netsta"], columns=["est", "bin"], values="cc_med").round(2).to_string())
        print(Q[Q.kind == "ctrl"].pivot_table(index=["netsta"], columns=["est", "bin"], values="cc_med").round(2).to_string())
        print(Q[Q.kind == "data"].pivot_table(index=["netsta"], columns=["est", "bin"], values="sig_med").round(1).to_string())
        print(H.round(2).to_string())
    print(json.dumps(verdict, indent=1))
    # figure: median cc vs bin type, line per station; panels per estimator; dashed = control
    fig, axs = plt.subplots(1, len(PK), figsize=(15, 4.2), sharey=True)
    cols = plt.get_cmap("tab10")
    for ax, k in zip(axs, PK):
        for i, ns in enumerate(STATIONS):
            d = Q[(Q.netsta == ns) & (Q.est == k) & (Q.kind == "data")].set_index("bin").reindex(BINS)
            c = Q[(Q.netsta == ns) & (Q.est == k) & (Q.kind == "ctrl")].set_index("bin").reindex(BINS)
            ax.plot(BINS, d.cc_med, marker="o", color=cols(i), label=ns.split(".")[1])
            ax.plot(BINS, c.cc_med, ls="--", lw=0.8, color=cols(i))
        ax.axhline(0.6, color="r", lw=0.8); ax.set_title(f"{k}" + (" (headline)" if k == "P3" else ""), fontsize=10)
        ax.set_xlabel("stack type (EP episode, YR year, Y2 two years)"); ax.grid(color="0.9")
    axs[0].set_ylabel("median P cc with reference"); axs[0].set_ylim(0, 1)
    axs[-1].legend(fontsize=8, loc="upper left", bbox_to_anchor=(1.01, 1))
    fig.suptitle("P waveform cc by stack type; solid = data, dashed = scrambled control; red = 0.6 rule", fontsize=10)
    fig.savefig(f"{FIG}/fig_P_cc{tag}.png", dpi=110, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    main()
