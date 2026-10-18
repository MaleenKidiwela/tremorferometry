#!/usr/bin/env python
"""E4 step 2 analysis — the q-control, the transfer, and the network gain, all held out.

Legs (all scored on the TEST half of the alternate-2-h-block split; the empirical moveout of each
leg was measured on the BUILD half only):
  A  the ring with the PUBLISHED SeisBench templates            -> reference registration
  B  the ring with OUR 2010 re-stacks built on the BUILD half    -> measured template quality q
  C  the ring (published) + the TRANSFERRED B011/B926/LRIV/STW   -> network gain

Outputs results/e4_step2_registration.csv, e4_step2_network.csv, e4_step2_q.csv
"""
from __future__ import annotations
import numpy as np, pandas as pd, os

RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
RING = ["PGC", "LZB", "VGZ", "SNB", "YOUB", "PFB"]
XTRA = ["B011", "B926", "LRIV", "STW"]
MB = [(-9, 1.5, "M<1.5"), (1.5, 1.8, "M1.5-1.8"), (1.8, 2.1, "M1.8-2.1"), (2.1, 9, "M>=2.1")]


def per_station(tag, stas):
    e = pd.read_parquet(f"{RES}/e4_tt_events_{tag}.parquet")
    st = pd.read_csv(f"{RES}/e4_tt_station_{tag}.csv")
    bg = st.groupby(["L", "sta"]).agg(rand=("rand3_median", "median"),
                                      randZ=("randZ_median", "median")
                                      if "randZ_median" in st.columns else ("rand3_median", "median"),
                                      mad_scale=("mad_scale", "median"),
                                      fa8=("fa8_per_day", "median")).reset_index()
    rows = []
    for s in stas:
        c = f"sta_{s}"
        if c not in e.columns:
            continue
        for L in sorted(e.L.unique()):
            r = bg[(bg.L == L) & (bg.sta == s)]
            if not len(r):
                continue
            rb, rz = float(r["rand"].iloc[0]), float(r["randZ"].iloc[0])
            g = e[(e.L == L) & e[c].notna()]
            if len(g) < 50:
                continue
            for lo, hi, lab in MB + [(-9, 9, "all")]:
                h = g[(g.mag > lo) & (g.mag <= hi)]
                if len(h) < 40:
                    continue
                rows.append(dict(leg=tag, sta=s, L=L, magbin=lab, n=len(h),
                                 mad3=float(h[c].median()), rand3=rb,
                                 excess3=float(h[c].median()) - rb,
                                 madZ=float(h[f"staZ_{s}"].median()), randZ=rz,
                                 excessZ=float(h[f"staZ_{s}"].median()) - rz,
                                 mad_scale=float(r.mad_scale.iloc[0]),
                                 fa8=float(r.fa8.iloc[0])))
    return pd.DataFrame(rows)


def network(tag):
    e = pd.read_parquet(f"{RES}/e4_tt_events_{tag}.parquet")
    nf = pd.read_csv(f"{RES}/e4_tt_netfam_{tag}.csv")
    nb = nf.groupby("L").agg(rand=("net_rand_med", "median"), fa6=("net_fa6", "median"),
                             fa8=("net_fa8", "median"), sfa8=("shift_fa8", "median"),
                             nch=("n_ch", "median")).reset_index()
    rows = []
    for L in sorted(e.L.unique()):
        r = nb[nb.L == L].iloc[0]
        g = e[e.L == L]
        for lo, hi, lab in MB + [(-9, 9, "all")]:
            h = g[(g.mag > lo) & (g.mag <= hi)]
            if len(h) < 40:
                continue
            rows.append(dict(leg=tag, L=L, magbin=lab, n=len(h), n_ch_median=float(r.nch),
                             net_mad=float(h.net_mad.median()), net_rand=float(r["rand"]),
                             excess=float(h.net_mad.median()) - float(r["rand"]),
                             p90=float(h.net_mad.quantile(.9)),
                             frac_ge6=float((h.net_mad >= 6).mean()),
                             frac_ge8=float((h.net_mad >= 8).mean()),
                             shift_ge8=float((h.net_mad_shift >= 8).mean()),
                             fa6=float(r.fa6), fa8=float(r.fa8), shift_fa8=float(r.sfa8)))
    return pd.DataFrame(rows)


def main():
    legs = {"A_test": RING, "B_test": RING, "C_test": RING + XTRA}
    S = pd.concat([per_station(t, s) for t, s in legs.items() if
                   os.path.exists(f"{RES}/e4_tt_events_{t}.parquet")], ignore_index=True)
    S.to_csv(f"{RES}/e4_step2_registration.csv", index=False)
    N = pd.concat([network(t) for t in legs if
                   os.path.exists(f"{RES}/e4_tt_events_{t}.parquet")], ignore_index=True)
    N.to_csv(f"{RES}/e4_step2_network.csv", index=False)

    # ---- q: the held-out ratio of our 2010 re-stack to the published stack, per station
    qc = pd.read_csv(f"{RES}/e4_transfer_stackqc_t2.csv")
    nev = qc.groupby("sta").n_events.median()
    r1 = qc.groupby("sta").r1.median()
    A = S[(S.leg == "A_test") & (S.magbin == "all")].set_index(["sta", "L"]).excess3
    B = S[(S.leg == "B_test") & (S.magbin == "all")].set_index(["sta", "L"]).excess3
    rows = []
    for (s, L), a in A.items():
        if (s, L) not in B.index:
            continue
        b = float(B.loc[(s, L)])
        n = float(nev.get(s, np.nan))
        rr = float(r1.get(s, np.nan))
        qpred = 1 / np.sqrt(1 + 1 / max(rr ** 2 * n, 1e-9)) if np.isfinite(rr) else np.nan
        rows.append(dict(sta=s, L=L, excess_seisbench=a, excess_2010=b,
                         ratio_2010_over_seisbench=b / a if a else np.nan,
                         n_events_stacked=n, r1=rr, q_pred_2010=qpred,
                         q_pred_ratio=qpred / 0.89 if np.isfinite(qpred) else np.nan,
                         cc_vs_seisbench=float(qc[qc.sta == s].cc_vs_seisbench.median())))
    Q = pd.DataFrame(rows)
    Q.to_csv(f"{RES}/e4_step2_q.csv", index=False)
    pd.set_option("display.width", 220)
    print("=== q-CONTROL: our 2010 re-stack vs the published SeisBench stack (held-out) ===")
    print(Q[Q.L == 8.0].round(3).to_string(index=False))
    print("\n=== per-station registration, M1.8-2.1, L = 8 s ===")
    print(S[(S.magbin == "M1.8-2.1") & (S.L == 8.0)].round(3).to_string(index=False))
    print("\n=== NETWORK (held-out test half) ===")
    print(N[N.magbin.isin(["M1.8-2.1", "all"])].round(4).to_string(index=False))

    # ---- network gain vs the per-station prediction
    gp = []
    for leg, stas in (("A_test", RING), ("C_test", RING + XTRA)):
        if not os.path.exists(f"{RES}/e4_tt_events_{leg}.parquet"):
            continue
        e = pd.read_parquet(f"{RES}/e4_tt_events_{leg}.parquet")
        st = pd.read_csv(f"{RES}/e4_tt_station_{leg}.csv")
        nf = pd.read_csv(f"{RES}/e4_tt_netfam_{leg}.csv")
        bg = st.groupby(["L", "sta"]).rand3_median.median()
        nb = nf.groupby("L").net_rand_med.median()
        for L in sorted(e.L.unique()):
            g = e[e.L == L]
            # per-FAMILY ideal from the stations that family actually has
            ideal_eq, ideal_opt, meas, w = [], [], [], []
            for fam, h in g.groupby("fam"):
                es = []
                for s in stas:
                    c = f"sta_{s}"
                    if c not in h.columns or h[c].notna().sum() < 10:
                        continue
                    es.append(float(h[c].median()) - float(bg.get((L, s), np.nan)))
                es = [x for x in es if np.isfinite(x)]
                if len(es) < 2 or h.net_mad.notna().sum() < 10:
                    continue
                ideal_eq.append(np.sum(es) / np.sqrt(len(es)))
                ideal_opt.append(np.sqrt(np.sum(np.square(es))))
                meas.append(float(h.net_mad.median()) - float(nb[L]))
                w.append(len(h))
            if not meas:
                continue
            gp.append(dict(leg=leg, L=L, n_fam=len(meas),
                           meas_excess=float(np.median(meas)),
                           ideal_equalweight=float(np.median(ideal_eq)),
                           ideal_optimal=float(np.median(ideal_opt)),
                           ratio_vs_equalweight=float(np.median(meas) / np.median(ideal_eq)),
                           ratio_vs_optimal=float(np.median(meas) / np.median(ideal_opt))))
    G = pd.DataFrame(gp)
    G.to_csv(f"{RES}/e4_step2_networkgain.csv", index=False)
    print("\n=== NETWORK GAIN vs the per-family prediction from the single-station values ===")
    print(G.round(3).to_string(index=False))
    print(f"\n[write] {RES}/e4_step2_{{registration,network,q,networkgain}}.csv")


if __name__ == "__main__":
    main()
