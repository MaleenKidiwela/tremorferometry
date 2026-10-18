#!/usr/bin/env python
"""bw_gate_g3_score -- G3 verdict: CN4 frozen templates (run a) vs CN4 + POLARIS/LZB era templates (run b) on the 2004-07
ETS, scored on the TEST (odd 2-h) blocks only (the era templates were built on the build blocks).

Recall = one-to-one, own family, +-2 s (svi_ref_common.one2one_recall) of Bostock's 2004 events in the test blocks of the
ETS days 2004-07-02..07-27; T_a = 8.25 / 8.75; T_b matched to (a)'s scrambled-control rate per observable family-day
(svi_ref_r3 rule, test-block rows).  Day bootstrap (1000) of recall(b) - recall(a).  Quiet days 2004-06-04..06-17:
data and control candidates per observable family-day at the same T.  Per family: carriers (frozen, era) and recall.

Out: OUT/tables/g3_recall.csv, g3_family_recall.csv, g3_quiet.csv, g3_decision.json; figures/g3_recall.png
"""
from __future__ import annotations
import os, sys, glob, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bw_gate_common as G
R, C = G.patch_forward()
TAB = f"{G.OUT}/tables"
ETS = ("2004-07-02", "2004-07-27"); QUIET = ("2004-06-04", "2004-06-17")
RUNS = dict(zip("ab", os.environ.get("BW_G3_RUNS", "g3a_CN4,g3b_CN4ERA").split(",")))
PFX = os.environ.get("BW_G3_TAG", "g3")      # output prefix (P3 SNB A/B re-uses this scorer with its own runs)


def load(run):
    fs = [p for p in sorted(glob.glob(f"{G.OUT}/mf/feat/{run}/*.parquet")) if not p.endswith("_fd.parquet") and os.path.getsize(p) > 0]
    A = pd.concat([pd.read_parquet(p, columns=["series", "fam", "date", "t", "net_mad", "n_sites"]) for p in fs], ignore_index=True)
    fd = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{G.OUT}/mf/feat/{run}/*_fd.parquet")) if os.path.getsize(p) > 0], ignore_index=True)
    return A, fd


def obs(fd, ser, days):
    return int(((fd.series == ser) & fd.observable.astype(bool) & fd.date.isin(days)).sum())


def main():
    ets_days = [str(x.date()) for x in pd.date_range(*ETS)]; q_days = [str(x.date()) for x in pd.date_range(*QUIET)]
    L = {k: load(v) for k, v in RUNS.items()}
    B = G.bostock_all(); B = B[B.date.isin(ets_days)]
    B = B[R.block_role(B.t_origin.to_numpy()) == "fit"]
    Sc = {}
    for k, (A, fd) in L.items():
        a = A[A.date.isin(ets_days)]; Sc[k] = a[R.block_role(a.t.to_numpy()) == "fit"]
    rows, frows, qrows, boots = [], [], [], {}
    for Ta in (8.25, 8.75):
        ca = Sc["a"]; cb = Sc["b"]
        target = ((ca.series == "c1") & (ca.net_mad >= Ta)).sum() / max(obs(L["a"][1], "c1", ets_days), 1)
        grid = np.round(np.arange(8.0, 14.0, 0.01), 2); mb = max(obs(L["b"][1], "c1", ets_days), 1)
        cvals = cb[cb.series == "c1"].net_mad.to_numpy()
        Tb = next((float(T) for T in grid if (cvals >= T).sum() / mb <= target + 1e-12), float(grid[-1]))
        res = {}
        for k, T in (("a", Ta), ("b", Tb)):
            D = Sc[k][(Sc[k].series == "data") & (Sc[k].net_mad >= T)]
            hit = R.one2one_recall(D.t.to_numpy(), D.fam.to_numpy(), B.t_origin.to_numpy(), B.fam.to_numpy())
            dt = R.nearest_dt(D.t.to_numpy(), D.fam.to_numpy(), B.t_origin.to_numpy(), B.fam.to_numpy(), True)
            fn = pd.Series(hit).groupby(B.fam.to_numpy()).sum()
            res[k] = dict(T=T, hit=hit, recall=hit.mean(), precision=float((dt <= 2.0).mean()) if len(D) else np.nan, n_data=len(D),
                          n_ctl=int(((Sc[k].series == "c1") & (Sc[k].net_mad >= T)).sum()), fam_red=int((fn >= 5).sum()), fn=fn,
                          ctl_rate=((Sc[k].series == "c1") & (Sc[k].net_mad >= T)).sum() / max(obs(L[k][1], "c1", ets_days), 1))
            Aq, fdq = L[k]; q = Aq[Aq.date.isin(q_days)]
            qrows.append(dict(T_ref=Ta, run=k, T=T, quiet_data_per_famday=((q.series == "data") & (q.net_mad >= T)).sum() / max(obs(fdq, "data", q_days), 1),
                              quiet_ctl_per_famday=((q.series == "c1") & (q.net_mad >= T)).sum() / max(obs(fdq, "c1", q_days), 1),
                              quiet_obs_famdays=obs(fdq, "data", q_days)))
        # day bootstrap of the recall difference
        dB = B.date.to_numpy(); ud = np.unique(dB); ix = {d: np.flatnonzero(dB == d) for d in ud}
        rs = np.random.default_rng(9); diffs = []
        for _ in range(1000):
            ii = np.concatenate([ix[d] for d in rs.choice(ud, len(ud), replace=True)])
            diffs.append(res["b"]["hit"][ii].mean() - res["a"]["hit"][ii].mean())
        ci = (float(np.percentile(diffs, 5)), float(np.percentile(diffs, 95)))
        rows.append(dict(T_ref=Ta, T_a=Ta, T_b=Tb, ctl_rate_a=res["a"]["ctl_rate"], ctl_rate_b=res["b"]["ctl_rate"], n_events=len(B),
                         recall_a=res["a"]["recall"], recall_b=res["b"]["recall"], ratio=res["b"]["recall"] / max(res["a"]["recall"], 1e-9),
                         diff_90ci_lo=ci[0], diff_90ci_hi=ci[1], precision_a=res["a"]["precision"], precision_b=res["b"]["precision"],
                         n_data_a=res["a"]["n_data"], n_data_b=res["b"]["n_data"], n_ctl_a=res["a"]["n_ctl"], n_ctl_b=res["b"]["n_ctl"],
                         fam_red_a=res["a"]["fam_red"], fam_red_b=res["b"]["fam_red"]))
        for f in range(130):
            frows.append(dict(T_ref=Ta, fam=f, n_events=int((B.fam == f).sum()), matched_a=int(res["a"]["fn"].get(f, 0)),
                              matched_b=int(res["b"]["fn"].get(f, 0))))
        G.log(f"[g3] T {Ta}: recall a {res['a']['recall']:.4f} -> b {res['b']['recall']:.4f} (x{rows[-1]['ratio']:.2f}, 90% CI of diff "
              f"{ci[0]:+.4f}..{ci[1]:+.4f}) at T_b {Tb}; precision {res['a']['precision']:.3f} -> {res['b']['precision']:.3f}; "
              f"families {res['a']['fam_red']} -> {res['b']['fam_red']}")
    Rr = pd.DataFrame(rows); Rr.to_csv(f"{TAB}/{PFX}_recall.csv", index=False)
    Fr = pd.DataFrame(frows)
    Ca = pd.read_csv(f"{TAB}/g3_carriers.csv")
    Fr = Fr.merge(Ca, on="fam", how="left"); Fr.to_csv(f"{TAB}/{PFX}_family_recall.csv", index=False)
    Q = pd.DataFrame(qrows); Q.to_csv(f"{TAB}/{PFX}_quiet.csv", index=False)
    ok = all((r.ratio >= 1.10) and (r.diff_90ci_lo > 0) and (r.fam_red_b >= r.fam_red_a) for r in Rr.itertuples())
    dec = dict(pass_=bool(ok), rows=Rr.to_dict("records"), quiet=Q.to_dict("records"),
               families_ge2_sites_cn4=int((Ca.n_frozen_cn4 >= 2).sum()), families_ge2_sites_cn4_era=int((Ca.n_total >= 2).sum()),
               median_carriers_cn4=float(Ca.n_frozen_cn4.median()), median_carriers_cn4_era=float(Ca.n_total.median()))
    json.dump(dec, open(f"{TAB}/{PFX}_decision.json", "w"), indent=1, default=float)
    G.log(f"[g3] PASS={ok}")
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    x = np.arange(len(Rr))
    ax[0].bar(x - .18, Rr.recall_a, .36, color="0.6", label="(a) CN4, frozen 2012 templates")
    ax[0].bar(x + .18, Rr.recall_b, .36, color="#2471a3", label="(b) CN4 + POLARIS/LZB era templates")
    ax[0].set_xticks(x); ax[0].set_xticklabels([f"control rate of (a) at {t}" for t in Rr.T_ref]); ax[0].legend(fontsize=8)
    ax[0].set_ylabel("held-out recall of Bostock's 2004-07 events"); ax[0].set_title("Test blocks, matched scrambled-control rate")
    h = Fr[Fr.T_ref == 8.75]
    ax[1].scatter(h.matched_a / h.n_events.clip(lower=1), h.matched_b / h.n_events.clip(lower=1), c=h.n_era, cmap="viridis", s=16)
    ax[1].plot([0, .6], [0, .6], "k-", lw=.6); ax[1].set_xlabel("per-family recall (a)"); ax[1].set_ylabel("per-family recall (b)")
    ax[1].set_title("per family (colour = number of era stations carrying it)")
    plt.tight_layout(); plt.savefig(f"{G.OUT}/figures/{PFX}_recall.png", dpi=130); plt.close()


if __name__ == "__main__":
    main()
