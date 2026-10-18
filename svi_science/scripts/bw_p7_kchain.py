#!/usr/bin/env python
"""bw_p7_kchain -- P7 item 5: completeness ladder k (bw_p4b_kladder.py, unchanged, with the P6 P) and the chain to the S15
(B4) level (rules: results/backward/p7/P7.md).
  primary k = env_quiet / rows all; fallback between, then all, if not finite or <= 0.
  K_B4 = 1; K_B3 = k(S15/B3); K_B2 = K_B3 k(B3/B2); K_B1 = K_B2 k(B2/B1); K_B0 = K_B1 k(B1/CN4);
  K_era = K_set / k(era/set) via CN4, B1, B2 overlaps (primary via B2).  sigma_ln per link = (ln hi68 - ln lo68)/2, quadrature.
Out: p7/tables/k_ladder.csv (copy), k_ladder_agreement.csv (copy), k_links.csv, k_chain.csv, k_chain.json
"""
import os, sys, json, shutil, subprocess
import numpy as np
import pandas as pd

S = "/home/jovyan/tremorferometry/svi_science/scripts"
BW = "/home/jovyan/tremorferometry/svi_science/results/backward"
P7 = f"{BW}/p7"


def pick(K, kind, pair, rows="all"):
    g = K[(K.kind == kind) & (K.pair == pair) & (K.rows == rows)].set_index("days")
    for cl in ("env_quiet", "between", "all"):
        if cl in g.index and np.isfinite(g.loc[cl, "k"]) and g.loc[cl, "k"] > 0:
            r = g.loc[cl]
            lo, hi = r.k_lo68, r.k_hi68
            sig = (np.log(hi) - np.log(lo)) / 2 if np.isfinite(lo) and np.isfinite(hi) and lo > 0 and hi > 0 else np.nan
            ks = g.k[(g.k > 0) & np.isfinite(g.k)]
            return dict(pair=pair, kind=kind, rows=rows, day_class=cl, k=float(r.k), sigma_ln=float(sig), n_days=int(r.n_days),
                        k_lo95=float(r.k_lo95), k_hi95=float(r.k_hi95), classes_min=float(ks.min()), classes_max=float(ks.max()),
                        classes_agree_x2=bool(ks.max() / ks.min() <= 2.0), k_env_quiet=float(g.k.get("env_quiet", np.nan)),
                        k_between=float(g.k.get("between", np.nan)), k_ets=float(g.k.get("ets", np.nan)), k_all=float(g.k.get("all", np.nan)))
    return dict(pair=pair, kind=kind, rows=rows, day_class=None, k=np.nan, sigma_ln=np.nan)


def main():
    r = subprocess.run([sys.executable, f"{S}/bw_p4b_kladder.py", "--p-dir", f"{BW}/p6/P"], capture_output=True, text=True)
    open(f"{P7}/logs/kladder.log", "w").write(r.stdout + r.stderr)
    if r.returncode:
        raise SystemExit(f"[kchain] kladder failed: {r.stderr[-2000:]}")
    for f in ("k_ladder.csv", "k_ladder_agreement.csv"):
        shutil.copy(f"{BW}/p4b/tables/{f}", f"{P7}/tables/{f}")
    K = pd.read_csv(f"{P7}/tables/k_ladder.csv")
    L = {}
    for rows in ("all", "ms"):
        for kind, pair in (("ladder", "S15/B3"), ("ladder", "B3/B2"), ("ladder", "B2/B1"), ("ladder", "B1/CN4"),
                           ("overlap", "era/CN4"), ("overlap", "era/B1"), ("overlap", "era/B2")):
            L[(pair, rows)] = pick(K, kind, pair, rows)
    pd.DataFrame(L.values()).to_csv(f"{P7}/tables/k_links.csv", index=False)
    out = []
    for rows in ("all", "ms"):
        def lk(p): return L[(p, rows)]
        chain = {"S15": (1.0, 0.0, [])}
        for new, old, pair in (("S15", "B3", "S15/B3"), ("B3", "B2", "B3/B2"), ("B2", "B1", "B2/B1"), ("B1", "CN4", "B1/CN4")):
            K0, s0, links = chain[new]; l = lk(pair)
            chain[old] = (K0 * l["k"], np.sqrt(s0 ** 2 + (l["sigma_ln"] if np.isfinite(l["sigma_ln"]) else np.nan) ** 2), links + [pair])
        eras = {}
        for via in ("B2", "B1", "CN4"):
            K0, s0, links = chain[via]; l = lk(f"era/{via}")
            eras[via] = (K0 / l["k"] if l["k"] > 0 else np.nan, np.sqrt(s0 ** 2 + l["sigma_ln"] ** 2), links + [f"era/{via}"])
        ev = np.array([v[0] for v in eras.values() if np.isfinite(v[0]) and v[0] > 0])
        chain["ERA"] = eras["B2"]
        for s, (Kv, sg, links) in chain.items():
            d = dict(rows=rows, net_set=s, K=Kv, sigma_ln=sg, K_lo68=Kv * np.exp(-sg), K_hi68=Kv * np.exp(sg),
                     K_lo95=Kv * np.exp(-1.96 * sg), K_hi95=Kv * np.exp(1.96 * sg), links=" x ".join(links))
            if s == "ERA":
                d.update(K_era_via_B2=eras["B2"][0], K_era_via_B1=eras["B1"][0], K_era_via_CN4=eras["CN4"][0],
                         era_estimates_agree_x2=bool(len(ev) >= 2 and ev.max() / ev.min() <= 2.0),
                         K_era_bracket_min=float(ev.min()) if len(ev) else np.nan, K_era_bracket_max=float(ev.max()) if len(ev) else np.nan)
            out.append(d)
    T = pd.DataFrame(out); T.to_csv(f"{P7}/tables/k_chain.csv", index=False)
    json.dump(T.to_dict("records"), open(f"{P7}/tables/k_chain.json", "w"), indent=1, default=float)
    pd.set_option("display.width", 250)
    print(K[K.rows == "all"].round(3).to_string(index=False)); print(pd.DataFrame(L.values()).round(3).to_string(index=False))
    print(T.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
