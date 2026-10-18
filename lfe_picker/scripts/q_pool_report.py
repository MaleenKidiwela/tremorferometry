#!/usr/bin/env python
"""Assemble the pooled (2010 + 2012) Step-1a repeat into the tables the coordinator asked for.

  * per isolated hidden family: same-family pairs above the pooled null's 1e-6 / 1e-5 / 1e-4
    quantile, in the protocol configuration (all 7 discovery sites);
  * the same-family enrichment over far-cross-family pairs, per configuration;
  * the arithmetic: pooled events per family, pairs available, pairs needed for the bar;
  * the SITE-SUBSET ablation rows (4 nearest discovery sites; the oracle 4 with B926), each family
    scored in its OWN subset;
  * the 2012-only control run through the same code, to show the cache-vs-mm-store path agrees with
    the original result;
  * the site-count extrapolation: the same-family shift in units of the null's spread, and what
    7 -> 20 sites would buy under the sqrt(N_sites) network scaling.
"""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
from scipy.stats import norm

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import p_common as P

P3 = f"{P.RES}/p3"
THRN = ["q_0.001", "q_0.0001", "q_1e-5", "q_1e-6"]
QV = dict(zip(THRN, [1e-3, 1e-4, 1e-5, 1e-6]))


def load(tag):
    return (json.load(open(f"{P3}/pairlink2_{tag}_null.json")),
            json.load(open(f"{P3}/pairlink2_{tag}_true.json")),
            pd.read_csv(f"{P3}/pairlink2_{tag}_true_by_family.csv"))


def main():
    ho = json.load(open(f"{P.RES}/holdouts.json"))
    strat = {int(h["bkey"]): h["stratum"] for h in ho["hidden_families"]}
    iso = [b for b, s in strat.items() if s == "isolated_ge5km"]
    N_, T_, F_ = load("pool")
    N12, T12, F12 = load("y2012")
    E = pd.read_parquet(f"{P3}/pairlink2_pool_true_events.parquet")
    thr = N_["thresholds"]
    all7 = F_[F_.config == "all7"].set_index("bkey")

    # ---------------- 1. isolated hidden families, protocol configuration
    rows = []
    ev = E.bkey.value_counts()
    ev10 = E[E.year == 2010].bkey.value_counts()
    ev12 = E[E.year == 2012].bkey.value_counts()
    for b in sorted(iso):
        r = all7.loc[b]
        rows.append(dict(bkey=b, ev_2010=int(ev10.get(b, 0)), ev_2012=int(ev12.get(b, 0)),
                         ev_pooled=int(ev.get(b, 0)), same_pairs=int(r.n_same_pairs),
                         **{k: int(r[k]) for k in THRN}))
    R = pd.DataFrame(rows)
    # 2012-only counts through the same code, for the same families
    a12 = F12[F12.config == "all7"].set_index("bkey")
    for k in THRN:
        R[f"{k}_2012only"] = [int(a12[k].get(b, 0)) for b in R.bkey]
    R["same_pairs_2012only"] = [int(a12.n_same_pairs.get(b, 0)) for b in R.bkey]
    R.to_csv(f"{P3}/pool_isolated_hidden.csv", index=False)

    # ---------------- 2. enrichment per configuration
    erows = []
    for c, meta in T_["configs"].items():
        sub = F_[F_.config == c]
        st = int(sub.n_same_pairs.sum())
        for k in THRN:
            sh = int(sub[k].sum())
            fh = int(T_["far_hit"][c][k]); ft = int(T_["far_tot"])
            ah = int(T_["alias_hit"][c][k]); at = int(T_["alias_tot"])
            erows.append(dict(config=c, kind=meta["kind"], n_sites=len(meta["sites"]),
                              n_families=meta["n_families"], threshold=thr[c][k], q=k,
                              same_pairs=st, same_above=sh,
                              same_rate=sh / max(1, st), far_rate=fh / max(1, ft),
                              alias_rate=ah / max(1, at),
                              enrichment=(sh / max(1, st)) / max(1e-15, fh / max(1, ft))))
    EN = pd.DataFrame(erows)
    EN.to_csv(f"{P3}/pool_enrichment_by_config.csv", index=False)

    # ---------------- 3. the ablation row: each family scored in its OWN 4-site subset
    arows = []
    for kind, pre in [("protocol, 7 sites", "all7"), ("4 nearest discovery sites", "n4_"),
                      ("oracle 4 incl. withheld B926", "o4_")]:
        sub = F_[F_.config.str.startswith(pre)] if pre != "all7" else F_[F_.config == "all7"]
        s_iso = sub[sub.bkey.isin(iso)]
        arows.append(dict(row=kind,
                          iso_same_pairs=int(s_iso.n_same_pairs.sum()),
                          **{f"iso_{k}": int(s_iso[k].sum()) for k in THRN},
                          iso_families_at_bar=int((s_iso["q_1e-6"] >= 30).sum()),
                          all_same_pairs=int(sub.n_same_pairs.sum()),
                          **{f"all_{k}": int(sub[k].sum()) for k in THRN}))
    AB = pd.DataFrame(arows)
    AB.to_csv(f"{P3}/pool_site_ablation.csv", index=False)

    # ---------------- 4. the arithmetic, and the site-count extrapolation
    # same-family shift delta, in units of the null's spread, from the tail rate at each quantile:
    #   rate_same(q) = Phi(-(z_q - delta)),  z_q = Phi^-1(1 - q)
    a7 = EN[(EN.config == "all7")]
    dlt = []
    for _, r in a7.iterrows():
        zq = norm.isf(QV[r.q])
        if r.same_rate > 0:
            dlt.append(zq - norm.isf(r.same_rate))
    delta = float(np.median(dlt))
    gains = {}
    for ns in (7, 10, 14, 20, 26):
        d = delta * np.sqrt(ns / 7.0)
        gains[ns] = dict(delta_sigma=d,
                         rate_at_1e6=float(norm.sf(norm.isf(1e-6) - d)),
                         pairs_for_30=float(30 / max(1e-30, norm.sf(norm.isf(1e-6) - d))))
    summary = dict(
        pooled_events=int(T_["n_events"]), by_year=T_["n_by_year"],
        pooled_pairs=int(T_["n_pairs"]), same_pairs=int(T_["same_tot"]),
        alias_pairs=int(T_["alias_tot"]), far_pairs=int(T_["far_tot"]),
        all7_thresholds=thr["all7"],
        isolated_at_bar=int((R["q_1e-6"] >= 30).sum()),
        isolated_best=int(R["q_1e-6"].max()),
        control_2012_events=int(T12["n_events"]), control_2012_pairs=int(T12["n_pairs"]),
        control_2012_thresholds=N12["thresholds"]["all7"],
        delta_sigma_same_family=delta, site_scaling=gains)
    json.dump(summary, open(f"{P3}/pool_summary.json", "w"), indent=1, default=str)

    pd.set_option("display.width", 200)
    print("== isolated hidden families, pooled 2010+2012, 7 discovery sites ==")
    print(R.to_string(index=False))
    print("\n== enrichment, protocol configuration ==")
    print(a7[["q", "threshold", "same_pairs", "same_above", "same_rate", "far_rate",
              "alias_rate", "enrichment"]].to_string(index=False))
    print("\n== site-subset ablation (each family in its own subset) ==")
    print(AB.to_string(index=False))
    print("\n== per-configuration enrichment at the 1e-6 quantile ==")
    print(EN[EN.q == "q_1e-6"][["config", "kind", "n_sites", "n_families", "threshold",
                                "same_rate", "far_rate", "enrichment"]].to_string(index=False))
    print("\n== site-count extrapolation (sqrt(N) network scaling) ==")
    print(f"same-family shift delta = {delta:.3f} sigma of the null at 7 sites")
    for ns, g in gains.items():
        print(f"  {ns:2d} sites: delta {g['delta_sigma']:.2f} sigma  "
              f"rate above 1e-6 {g['rate_at_1e6']:.3e}  "
              f"same-family pairs needed for 30 above threshold: {g['pairs_for_30']:.3g} "
              f"(events {np.sqrt(2*g['pairs_for_30']):.0f})")
    print(f"\n[done] -> {P3}/pool_*.csv/json")


if __name__ == "__main__":
    main()
