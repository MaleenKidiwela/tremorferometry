#!/usr/bin/env python
"""Track A-2 STEP 4 — combine the withheld-station stacks into the epoch-identity verdict.

For each (station, boundary) it joins the NEW-epoch leg with the same-procedure 2012 CONTROL leg and
applies the pre-registered bars per family-station:
    |cc| >= 0.6   (cc <= -0.6 -> POLARITY FLIP, the template is inverted before it is saved)
    |dS| <= 0.5 s
    registration >= 0.8 x the previous epoch  (the previous-epoch value is the control leg's
                                               registration for the SAME family and station)
The control leg also gives the FLOOR of the procedure: what |cc| and |dS| come out when the epoch has
NOT changed.  Both are reported, and a family that fails on the new leg while its control leg also
fails is flagged `floor_fail` rather than counted against the instrument change.

Writes results/a4_epoch_identity.csv (the deliverable), results/a4_templates_<STA>_<epoch>.npz for the
accepted families (--extra-tpl format, ready for Track C) and figures/a4_epoch_identity_<STA>.png.
"""
from __future__ import annotations
import argparse, glob, json, os, sys
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from a4_epoch_stack import best_cc                             # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)
REP = "/home/jovyan/tremorferometry/lfe_catalog_repair"
FS, NCMP = 20.0, 160


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", default="PGC,VGZ,NLLB,SNB,LRIV,VGZ17")
    ap.add_argument("--min-cc", type=float, default=0.6)
    ap.add_argument("--max-ds", type=float, default=0.5)
    ap.add_argument("--min-reg-ratio", type=float, default=0.8)
    a = ap.parse_args()
    EPOCH = {"PGC": ("2018_HH", "2017-08-10"), "VGZ": ("2018_HH", "2017-03-03"),
             "NLLB": ("2018_HH", "2017-02-23"), "SNB": ("2018_HH", "2017-06-08"),
             "LRIV": ("2021_HH", "2019-05-30"), "VGZ17": ("2017_HH", "2017-03-03"),
             "B926": ("2018_EH", "none (negative control)")}
    Z12 = np.load(f"{REP}/results/e5r_templates_2012.npz", allow_pickle=True)
    allrows, summary = [], []
    for key in a.stations.split(","):
        sta = "VGZ" if key == "VGZ17" else key
        epoch, bnd = EPOCH[key]
        ctl_leg = "oldm" if os.path.exists(f"{EXT}/results/a4_epoch_identity_{key}_oldm.csv") else "old"
        fn = {"new": f"{EXT}/results/a4_epoch_identity_{key}_new.csv",
              "old": f"{EXT}/results/a4_epoch_identity_{key}_{ctl_leg}.csv"}
        if not all(os.path.exists(v) for v in fn.values()):
            print(f"[{key}] missing legs {[k for k, v in fn.items() if not os.path.exists(v)]}, skipped")
            continue
        N = pd.read_csv(fn["new"]); O = pd.read_csv(fn["old"])
        N = N[N.status == "ok"]; O = O[O.status == "ok"]
        M = N.merge(O[["family", "cc", "dS", "registration", "n_events", "snr_db", "snr_null_db"]],
                    on="family", suffixes=("", "_ctl"), how="left")
        M["reg_ratio"] = M.registration / M.registration_ctl
        M["polarity_flip"] = M.cc <= -a.min_cc
        M["pass_cc"] = M.cc.abs() >= a.min_cc
        M["pass_dS"] = M.dS.abs() <= a.max_ds
        M["pass_reg"] = M.reg_ratio >= a.min_reg_ratio
        M["ctl_ok"] = (M.cc_ctl.abs() >= a.min_cc) & (M.dS_ctl.abs() <= a.max_ds)
        M["verdict"] = np.where(M.pass_cc & M.pass_dS & M.pass_reg, "ACCEPT",
                        np.where(~M.ctl_ok, "floor_fail", "REJECT"))
        M["station"] = sta; M["boundary"] = bnd; M["epoch"] = epoch; M["key"] = key
        allrows.append(M)
        tested = M[M.ctl_ok]
        summary.append(dict(
            key=key, station=sta, boundary=bnd, epoch=epoch, control_leg=ctl_leg,
            n_families=len(M), n_ctl_ok=int(len(tested)),
            n_accept=int((M.verdict == "ACCEPT").sum()),
            frac_accept_of_ctl_ok=float((tested.verdict == "ACCEPT").mean()) if len(tested) else np.nan,
            frac_fail_of_ctl_ok=float((tested.verdict == "REJECT").mean()) if len(tested) else np.nan,
            n_polarity_flip=int(M.polarity_flip.sum()),
            cc_median=float(M.cc.median()), cc_ctl_median=float(M.cc_ctl.median()),
            dS_median=float(M.dS.median()), dS_ctl_median=float(M.dS_ctl.median()),
            dS_mad=float((M.dS - M.dS.median()).abs().median()),
            reg_median=float(M.registration.median()), reg_ctl_median=float(M.registration_ctl.median()),
            reg_ratio_median=float(M.reg_ratio.median()),
            n_events_median=float(M.n_events.median()),
            audit_stationxml=bool(len(tested) and (tested.verdict == "REJECT").mean() > 0.5)))
        print(f"[{key}] {summary[-1]}", flush=True)

        # ---- CROSS-FAMILY chance level of the |cc| >= 0.6 bar (Merlin step-3 criterion (d)) -------
        src0 = f"{EXT}/results/a4_stacks_{key}_{epoch}_new.npz"
        xf = {}
        if os.path.exists(src0):
            Zx = np.load(src0, allow_pickle=True)
            qcall = pd.read_csv(f"{REP}/results/e5r_template_qc_2012.csv")
            qcall = qcall[(qcall.sta == sta) & qcall.usable & qcall.template].set_index("bkey")
            fams = [int(f) for f in M.family if f"{sta}|{int(f)}" in Zx.files]
            rng = np.random.default_rng(3)
            vals = []
            for _ in range(min(400, max(0, len(fams) * (len(fams) - 1)))):
                f1, f2 = rng.choice(fams, 2, replace=False)
                if f2 not in qcall.index:
                    continue
                ref = np.asarray(Z12[f"{sta}|{int(f2)}"], float)
                i_r = 400 + int(round(float(qcall.loc[f2, "s_pick_s"]) * FS)) - int(round(1.0 * FS))
                c, _ = best_cc(ref, np.asarray(Zx[f"{sta}|{int(f1)}"], float), i_r)
                if np.isfinite(c):
                    vals.append(abs(c))
            if vals:
                xf = dict(xfam_n=len(vals), xfam_p95=float(np.percentile(vals, 95)),
                          xfam_median=float(np.median(vals)),
                          xfam_frac_ge_mincc=float(np.mean(np.array(vals) >= a.min_cc)))
        summary[-1].update(xf)

        # ---- per-epoch template set for Track C (accepted families only, polarity applied) -------
        src = f"{EXT}/results/a4_stacks_{key}_{epoch}_new.npz"
        if os.path.exists(src):
            Zs = np.load(src, allow_pickle=True)
            tpl, sp, stk = Zs["tpl"].item(), Zs["sp"].item(), {}
            keep_f = set(M.loc[M.verdict == "ACCEPT", "family"].astype(int))
            flip_f = set(M.loc[M.polarity_flip, "family"].astype(int))
            T2, S2 = {}, {}
            for (f, s), d in tpl.items():
                if int(f) not in keep_f:
                    continue
                sgn = -1.0 if int(f) in flip_f else 1.0
                T2[(float(f), sta)] = {L: (sgn * v).astype(np.float32) for L, v in d.items()}
                S2[(float(f), sta)] = sp[(f, s)]
                k = f"{s}|{int(f)}"
                if k in Zs.files:
                    stk[f"{sta}|{int(f)}"] = np.asarray(Zs[k])
            np.savez_compressed(f"{EXT}/results/a4_templates_{key}_{epoch}.npz",
                                tpl=np.array(T2, dtype=object), sp=np.array(S2, dtype=object), **stk)
            print(f"  -> results/a4_templates_{key}_{epoch}.npz: {len(T2)} accepted templates "
                  f"({len(flip_f & keep_f)} polarity-inverted)", flush=True)

        # ---- figure ------------------------------------------------------------------------------
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        pick = M.dropna(subset=["cc"]).sort_values("n_events", ascending=False).head(6)
        fig, axs = plt.subplots(3, 3, figsize=(15, 8.5))
        for i, r in enumerate(pick.itertuples()):
            ax = axs.ravel()[i]
            k12 = f"{sta}|{int(r.family)}"
            if k12 not in Z12.files or not os.path.exists(src):
                continue
            ref = np.asarray(Z12[k12], float)
            qc = pd.read_csv(f"{REP}/results/e5r_template_qc_2012.csv")
            qc = qc[(qc.sta == sta) & (qc.bkey == int(r.family))]
            i_ref = 400 + int(round(float(qc.s_pick_s.iloc[0]) * FS)) - int(round(1.0 * FS))
            new = np.asarray(Zs[f"{sta}|{int(r.family)}"], float) if f"{sta}|{int(r.family)}" in Zs.files else None
            t = np.arange(NCMP) / FS - 1.0
            for c in range(3):
                A_ = ref[c, i_ref:i_ref + NCMP]
                A_ = A_ / (np.abs(A_).max() or 1)
                ax.plot(t, A_ + 2.2 * (2 - c), "k", lw=1.0, label="2012 (previous epoch)" if c == 0 else None)
                if new is not None:
                    j = i_ref + int(round(r.dS * FS))
                    if 0 <= j and j + NCMP <= new.shape[1]:
                        B_ = new[c, j:j + NCMP] * (-1 if r.polarity_flip else 1)
                        B_ = B_ / (np.abs(B_).max() or 1)
                        ax.plot(t, B_ + 2.2 * (2 - c), "r", lw=1.0,
                                label=f"{epoch} (new epoch)" if c == 0 else None)
            ax.set_title(f"fam {int(r.family)}  cc {r.cc:+.2f} (ctl {r.cc_ctl:+.2f})  "
                         f"dS {r.dS:+.2f} s  N {int(r.n_events)}", fontsize=8)
            ax.set_yticks([]); ax.set_xlabel("s from S")
            if i == 0:
                ax.legend(fontsize=7)
        ax = axs.ravel()[6]
        ax.hist(M.cc.dropna(), bins=np.arange(-1, 1.05, .1), color="tab:red", alpha=.7, label="new")
        ax.hist(M.cc_ctl.dropna(), bins=np.arange(-1, 1.05, .1), histtype="step", color="k", label="2012 control")
        ax.axvline(a.min_cc, ls=":", c="k"); ax.axvline(-a.min_cc, ls=":", c="k")
        ax.set_xlabel("signed cc vs the 2012 stack"); ax.legend(fontsize=7)
        ax = axs.ravel()[7]
        ax.hist(M.dS.dropna(), bins=np.arange(-2, 2.05, .1), color="tab:red", alpha=.7, label="new")
        ax.hist(M.dS_ctl.dropna(), bins=np.arange(-2, 2.05, .1), histtype="step", color="k", label="control")
        ax.axvline(a.max_ds, ls=":", c="k"); ax.axvline(-a.max_ds, ls=":", c="k")
        ax.set_xlabel("S shift dS (s)"); ax.legend(fontsize=7)
        ax = axs.ravel()[8]
        ax.scatter(M.registration_ctl, M.registration, s=14, c="tab:blue")
        lim = [0, max(0.3, float(np.nanmax(M[["registration", "registration_ctl"]].values)) * 1.1)]
        ax.plot(lim, lim, "k-", lw=.8); ax.plot(lim, [a.min_reg_ratio * x for x in lim], "k:", lw=.8)
        ax.set_xlabel("registration, 2012 control"); ax.set_ylabel(f"registration, {epoch}")
        fig.suptitle(f"{sta} epoch identity across {bnd}: {summary[-1]['n_accept']}/{summary[-1]['n_ctl_ok']} "
                     f"of control-OK families ACCEPT", fontsize=11)
        plt.tight_layout(rect=[0, 0, 1, .96])
        plt.savefig(f"{EXT}/figures/a4_epoch_identity_{key}.png", dpi=125)
        plt.close()

    if not allrows:
        raise SystemExit("no legs found")
    R = pd.concat(allrows, ignore_index=True)
    cols = ["station", "boundary", "epoch", "family", "n_events", "cc", "dS", "reg_ratio", "verdict",
            "polarity_flip", "registration", "registration_ctl", "cc_ctl", "dS_ctl", "n_events_ctl",
            "snr_db", "snr_null_db", "epi_km", "key"]
    R[cols].to_csv(f"{EXT}/results/a4_epoch_identity.csv", index=False)
    S = pd.DataFrame(summary)
    S.to_csv(f"{EXT}/results/a4_epoch_summary.csv", index=False)
    print("\n" + S.to_string(index=False))
    print("[write] results/a4_epoch_identity.csv, results/a4_epoch_summary.csv")


if __name__ == "__main__":
    main()
