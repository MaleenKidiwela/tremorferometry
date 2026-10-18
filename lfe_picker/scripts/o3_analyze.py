#!/usr/bin/env python
"""P0-O item 3, part 2 -- read the oracle and apply the PRE-REGISTERED decision.

Inputs: results/oracle/o3_*_oracle2012.*  (written by scripts/o3_oracle.py)
Outputs: results/o3_oracle_decision.json, results/o3_station_roc.csv,
         results/o3_network_curve.csv, figures/o3_*.png
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PICK = os.path.dirname(HERE)
RES = f"{PICK}/results"
FIG = f"{PICK}/figures"
ORC = f"{RES}/oracle"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
HF_CUT = 3.6384
TARGETS = [12.0, 116.0]          # network-wide quiet-day FALSE EVENTS per day (pre-registered)
DIAG_RATES = [12.0, 116.0, 1000.0, 10000.0, 100000.0]
C_E, C_C, C_T = "#eb6834", "#2a78d6", "#52514e"


def auc(pos, neg):
    pos = np.asarray(pos, float); neg = np.asarray(neg, float)
    pos = pos[np.isfinite(pos)]; neg = neg[np.isfinite(neg)]
    if len(pos) < 10 or len(neg) < 10:
        return np.nan
    allv = np.concatenate([pos, neg])
    r = pd.Series(allv).rank().to_numpy()
    return float((r[:len(pos)].sum() - len(pos) * (len(pos) + 1) / 2) / (len(pos) * len(neg)))


def thr_at_rate(peaks, n_days, rate):
    """threshold whose declustered-peak count per day equals `rate`."""
    p = np.sort(np.asarray(peaks, float))[::-1]
    k = int(round(rate * n_days))
    if k <= 0:
        return np.inf
    if k >= len(p):
        return float(p[-1]) if len(p) else -np.inf
    return float(p[k - 1])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="oracle2012")
    a = ap.parse_args()
    os.makedirs(FIG, exist_ok=True)
    M = json.load(open(f"{ORC}/o3_meta_{a.tag}.json"))
    P = pd.read_parquet(f"{ORC}/o3_net_pos_{a.tag}.parquet")
    N = pd.read_parquet(f"{ORC}/o3_net_neg_{a.tag}.parquet")
    PS = pd.read_parquet(f"{ORC}/o3_sta_pos_{a.tag}.parquet")
    NS = pd.read_parquet(f"{ORC}/o3_sta_neg_{a.tag}.parquet")
    Z = np.load(f"{ORC}/o3_fa_peaks_{a.tag}.npz", allow_pickle=True)
    SP = np.load(f"{ORC}/o3_sta_peaks_{a.tag}.npz", allow_pickle=True)
    nqf = int(Z["n_quiet_famday"]); nfam = float(Z["nfam_median"])
    gt = pd.read_parquet(f"{IN}/e5r_gt_bostock2012.parquet")
    N_TEST = int(((gt.half == "test") & (~gt.guard)).sum())
    print(f"[oracle] {len(P)}/{N_TEST} test-half positives evaluated; {nqf} quiet family-days; "
          f"{nfam:.0f} observable families/day")

    # ---------------- the FROZEN TEMPLATE detector, same quiet days, same rate axis --------------
    D = pd.read_parquet(f"{EXT}/results/e5r_detections_a1_mergeB_L12.parquet")
    H = pd.read_parquet(f"{EXT}/results/a1_hf_ratio_a1_mergeB_L12.parquet")[["doy", "fam", "t", "hfr_PGC"]]
    D = D[D.L == 12.0].merge(H, on=["doy", "fam", "t"], how="left")
    D = D[~(np.isfinite(D.hfr_PGC) & (D.hfr_PGC > HF_CUT))]
    bud = pd.read_csv(f"{EXT}/results/e5_budget_a1_mergeB_L12.csv")
    bud = bud[bud.L == 12.0]
    T0b = sorted(bud["T"].unique())[0]
    fdq = bud[(bud["T"] == T0b) & bud.doy.isin(range(228, 243))][["doy", "fam"]].drop_duplicates()
    n_tplday_q = len(fdq)
    Dq = D[D.doy.isin(range(228, 243))]
    ev = gt[(gt.half == "test") & (~gt.guard)].reset_index(drop=True)

    def one2one(x, b_sorted, tol=2.0):
        x = np.asarray(x, float); out = np.full(len(x), np.inf)
        if len(x) == 0 or len(b_sorted) == 0:
            return out
        j = 0; nb = len(b_sorted)
        for i in np.argsort(x, kind="stable"):
            ti = x[i]
            while j < nb and b_sorted[j] < ti - tol:
                j += 1
            if j < nb and abs(b_sorted[j] - ti) <= tol:
                out[i] = abs(b_sorted[j] - ti); j += 1
        return out

    def tpl_point(T):
        fa_tpl = len(Dq[Dq.net_mad >= T]) / max(n_tplday_q, 1)
        Dt = D[D.net_mad >= T]
        nh = 0
        for bk, g in ev.groupby("bkey"):
            td = np.sort(Dt[Dt.fam == float(bk)].t.to_numpy())
            nh += int((one2one(g.t_origin.to_numpy(), td) <= 2.0).sum())
        return fa_tpl * nfam, nh / N_TEST

    tgrid = np.arange(6.0, 12.01, 0.25)
    tpl_curve = [dict(T=float(T), net_fa_day=len(Dq[Dq.net_mad >= T]) / max(n_tplday_q, 1) * nfam)
                 for T in tgrid]
    tpl_T = {}
    for tgt in TARGETS:
        f = np.array([c["net_fa_day"] for c in tpl_curve])
        lg = np.log(np.maximum(f, 1e-6))
        t = np.log(tgt)
        i = int(np.searchsorted(-lg, -t)); i = min(max(i, 1), len(tgrid) - 1)
        w = (lg[i - 1] - t) / max(lg[i - 1] - lg[i], 1e-9)
        tpl_T[tgt] = float(tgrid[i - 1] + w * (tgrid[i] - tgrid[i - 1]))

    # ---------------- the two ORACLE network legs ----------------
    legs = {}
    curves = []
    for nm, col in (("ENERGY", "E"), ("CC", "C")):
        pk = np.asarray(Z[col], float)
        grid = np.percentile(pk, np.linspace(0, 100, 400)) if len(pk) else np.zeros(0)
        grid = np.unique(np.concatenate([grid, np.linspace(pk.min(), pk.max(), 400)])) if len(pk) else grid
        v = P[col].to_numpy()
        fa = np.array([(pk >= T).sum() / nqf * nfam for T in grid])
        rc = np.array([(v >= T).sum() / N_TEST for T in grid])
        for T, f_, r_ in zip(grid, fa, rc):
            curves.append(dict(leg=nm, T=float(T), net_fa_day=float(f_), recall=float(r_)))
        pt = {}
        for tgt in TARGETS:
            T = thr_at_rate(pk, nqf / nfam, tgt)      # peaks/day over the whole network
            pt[tgt] = dict(T=float(T), recall=float((v >= T).mean() * len(v) / N_TEST),
                           net_fa_day=float((pk >= T).sum() / nqf * nfam))
        legs[nm] = pt
    pd.DataFrame(curves).to_csv(f"{RES}/o3_network_curve.csv", index=False)

    tpl_pt = {}
    for tgt in TARGETS:
        fa_, rc_ = tpl_point(tpl_T[tgt])
        tpl_pt[tgt] = dict(T=tpl_T[tgt], recall=rc_, net_fa_day=fa_)

    # ---------------- single-station ROC ----------------
    rocs = []
    for s in sorted(PS.sta.unique()):
        p = PS[PS.sta == s]; nq = NS[(NS.sta == s) & (NS.cls == "quiet")]
        ne = NS[(NS.sta == s) & (NS.cls == "ets")]
        row = dict(sta=s, n_pos=len(p), n_neg_quiet=len(nq), n_neg_ets=len(ne),
                   auc_E_quiet=auc(p.zE, nq.zE), auc_E_ets=auc(p.zE, ne.zE),
                   auc_C_quiet=auc(p.zC, nq.zC), auc_C_ets=auc(p.zC, ne.zC),
                   median_zE_pos=float(p.zE.median()), median_zC_pos=float(p.zC.median()))
        # recall at 100 declustered peaks per station-day (ENERGY) / per template-day (CC)
        nsd = json.loads(str(SP["n_quiet_staday"])); ntd = json.loads(str(SP["n_quiet_tpl_staday"]))
        for st_, key, nd, colp in (("E", f"E_{s}", nsd.get(s, 0), "zE"),
                                   ("C", f"C_{s}", ntd.get(s, 0), "zC")):
            pk = np.asarray(SP[key], float) if key in SP.files else np.zeros(0)
            T = thr_at_rate(pk, nd, 100.0) if nd else np.nan
            row[f"T100_{st_}"] = T
            row[f"recall100_{st_}"] = float((p[colp].to_numpy() >= T).mean()) if np.isfinite(T) else np.nan
            row[f"peaks_per_day_unit_{st_}"] = "station-day" if st_ == "E" else "template-day"
        rocs.append(row)
    RC = pd.DataFrame(rocs)
    RC.to_csv(f"{RES}/o3_station_roc.csv", index=False)

    # ---------------- decision ----------------
    dec = {}
    for tgt in TARGETS:
        rE = legs["ENERGY"][tgt]["recall"]; rC = legs["CC"][tgt]["recall"]
        rT = tpl_pt[tgt]["recall"]
        dec[str(tgt)] = dict(
            target_network_fa_per_day=tgt,
            energy=legs["ENERGY"][tgt], cc_same_combination=legs["CC"][tgt],
            frozen_template_detector=tpl_pt[tgt],
            ratio_energy_over_frozen_template=rE / rT if rT else np.nan,
            ratio_energy_over_cc_same_combination=rE / rC if rC else np.nan)
    # ---- diagnostics: the whole rate axis, and what it costs each leg to match the template ----
    diag = []
    for rate in DIAG_RATES:
        row = dict(net_fa_per_day=rate)
        for nm, col in (("ENERGY", "E"), ("CC", "C")):
            pk = np.asarray(Z[col], float)
            T = thr_at_rate(pk, nqf / nfam, rate)
            row[f"recall_{nm}"] = float((P[col].to_numpy() >= T).sum() / N_TEST)
            row[f"T_{nm}"] = float(T)
        f = np.array([c["net_fa_day"] for c in tpl_curve])
        if rate > f.max():
            # the Track A run only recorded detections down to ~7 MAD: a looser rate than the
            # detector ever produced cannot be read off, and is NOT extrapolated.
            row["recall_TEMPLATE"] = np.nan; row["T_TEMPLATE"] = np.nan
            row["note"] = f"beyond the recorded range (max {f.max():.0f}/day at T={tgrid[0]:g})"
        else:
            lg = np.log(np.maximum(f, 1e-6)); t = np.log(rate)
            i = int(np.searchsorted(-lg, -t)); i = min(max(i, 1), len(tgrid) - 1)
            w = (lg[i - 1] - t) / max(lg[i - 1] - lg[i], 1e-9)
            Tt = float(tgrid[i - 1] + w * (tgrid[i] - tgrid[i - 1]))
            _, rr_ = tpl_point(Tt)
            row["recall_TEMPLATE"] = rr_; row["T_TEMPLATE"] = Tt
            row["note"] = ""
        diag.append(row)
    DG = pd.DataFrame(diag)
    DG.to_csv(f"{RES}/o3_rate_axis.csv", index=False)
    # the FA rate at which each leg reaches the frozen T_core recall (0.1463)
    target_recall = tpl_pt[12.0]["recall"]
    fa_needed = {}
    for nm, col in (("ENERGY", "E"), ("CC", "C")):
        pk = np.sort(np.asarray(Z[col], float))[::-1]
        v = np.sort(P[col].to_numpy())[::-1]
        k = int(round(target_recall * N_TEST))
        T = float(v[k - 1]) if k <= len(v) else -np.inf
        fa_needed[nm] = dict(T=T, net_fa_per_day=float((pk >= T).sum() / nqf * nfam))
    nauc = {}
    for nm, col in (("ENERGY", "E"), ("CC", "C")):
        nauc[nm] = dict(vs_quiet=auc(P[col], N[N.cls == "quiet"][col]),
                        vs_ets=auc(P[col], N[N.cls == "ets"][col]))

    ratios = [dec[str(t)]["ratio_energy_over_frozen_template"] for t in TARGETS]
    rmax = float(np.nanmax(ratios))
    verdict = ("SENSITIVITY CLAIM ALLOWED" if rmax >= 1.5 else
               "DISCOVERY-ONLY" if rmax < 1.25 else "REPORT BOTH READINGS")
    OUT = dict(tag=a.tag, n_test_half_events=N_TEST, n_evaluated=int(len(P)),
               n_quiet_family_days=nqf, observable_families_per_day=nfam,
               statistic_definitions=dict(
                   ENERGY="log10 E(S-1..S+5 s) / E(same 6-s window 20 s earlier), 1-8 Hz, MAD-normalised "
                          "per station-day, site-merged by the mean, summed over the 9 sites",
                   CC="3-C normalised cc to the 2012 family template, MAD-normalised per "
                      "(family, station, day), same site merge and sum",
                   FROZEN_TEMPLATE="the production detector: SNR-weighted mergeB network cc in units "
                                   "of that family-day's own network MAD"),
               points=dec, pre_registered_rule=dict(discovery_only="<1.25x", sensitivity=">=1.5x",
                                                    both="1.25-1.5x"),
               max_ratio_energy_over_template=rmax, VERDICT=verdict,
               network_auc=nauc, rate_axis=DG.to_dict("records"),
               fa_needed_to_match_template_Tcore_recall=dict(target_recall=target_recall, **fa_needed),
               single_station=dict(
                   median_auc_E_quiet=float(RC.auc_E_quiet.median()),
                   median_auc_C_quiet=float(RC.auc_C_quiet.median()),
                   median_auc_E_ets=float(RC.auc_E_ets.median()),
                   median_auc_C_ets=float(RC.auc_C_ets.median()),
                   best_station_E=str(RC.loc[RC.auc_E_quiet.idxmax(), "sta"]),
                   best_station_C=str(RC.loc[RC.auc_C_quiet.idxmax(), "sta"]),
                   recall_at_100_peaks_per_day=RC.set_index("sta")[["recall100_E", "recall100_C"]]
                   .to_dict("index")))
    json.dump(OUT, open(f"{RES}/o3_oracle_decision.json", "w"), indent=1, default=str)

    # ---------------- figures ----------------
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
    w = 0.38; x = np.arange(len(RC))
    ax[0].bar(x - w / 2, RC.auc_E_quiet, w, color=C_E, label="energy excess")
    ax[0].bar(x + w / 2, RC.auc_C_quiet, w, color=C_C, label="template cc")
    ax[0].set_xticks(x); ax[0].set_xticklabels(RC.sta, rotation=90, fontsize=8)
    ax[0].axhline(0.5, color=C_T, lw=0.8, ls=":")
    ax[0].set_ylabel("AUC vs quiet-day windows"); ax[0].set_title("single-station, quiet negatives")
    ax[0].legend(fontsize=8, frameon=False); ax[0].set_ylim(0.4, 1.0)
    ax[1].bar(x - w / 2, RC.auc_E_ets, w, color=C_E)
    ax[1].bar(x + w / 2, RC.auc_C_ets, w, color=C_C)
    ax[1].set_xticks(x); ax[1].set_xticklabels(RC.sta, rotation=90, fontsize=8)
    ax[1].axhline(0.5, color=C_T, lw=0.8, ls=":")
    ax[1].set_title("single-station, ETS-day negatives (>= 120 s from any Bostock event)")
    ax[1].set_ylim(0.4, 1.0)
    fig.tight_layout(); fig.savefig(f"{FIG}/o3_station_roc.png", dpi=150); plt.close(fig)

    CU = pd.DataFrame(curves)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    for nm, c in (("ENERGY", C_E), ("CC", C_C)):
        g = CU[CU.leg == nm].sort_values("net_fa_day")
        g = g[g.net_fa_day > 0]
        ax.plot(g.net_fa_day, g.recall, color=c, lw=1.8, label=nm.lower() + " network")
    tg = pd.DataFrame(tpl_curve)
    rr = []
    for T in tg["T"]:
        _, r_ = tpl_point(float(T))
        rr.append(r_)
    tg["recall"] = rr
    tg = tg[tg.net_fa_day > 0]
    ax.plot(tg.net_fa_day, tg.recall, color=C_T, lw=1.8, ls="--", label="frozen template detector")
    for tgt in TARGETS:
        ax.axvline(tgt, color="#c3c2b7", lw=0.8)
        ax.text(tgt, 0.01, f" {tgt:g}/day", fontsize=7, color="#77746c", rotation=90, va="bottom")
    ax.set_xscale("log"); ax.set_xlabel("network-wide quiet-day false events / day")
    ax.set_ylabel("own recall of Bostock 2012 test-half events")
    ax.set_title("P0-O headroom oracle (2012)"); ax.legend(fontsize=8, frameon=False)
    ax.grid(alpha=0.25, lw=0.5)
    fig.tight_layout(); fig.savefig(f"{FIG}/o3_network_roc.png", dpi=150); plt.close(fig)
    tg.to_csv(f"{RES}/o3_template_curve.csv", index=False)

    print(json.dumps(dec, indent=1, default=str))
    print(DG.to_string(index=False))
    print("network AUC:", json.dumps(nauc))
    print("FA/day needed to match the template T_core recall:", json.dumps(fa_needed))
    print(f"\nVERDICT: {verdict}  (max ratio energy/template = {rmax:.3f})")
    print(RC[["sta", "auc_E_quiet", "auc_C_quiet", "auc_E_ets", "auc_C_ets",
              "recall100_E", "recall100_C"]].to_string(index=False))


if __name__ == "__main__":
    main()
