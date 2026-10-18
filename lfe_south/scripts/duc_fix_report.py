#!/usr/bin/env python
"""DUC_FIX report -- assemble the catalogue and write results/DUC_EXT_FIXED_RESULTS.md.
  final detections per month = pass-2 (mf2/<m>.parquet) where it exists, else pass-1 (mf1/<m>)
  catalogue_all.parquet        every detection (8 x hourly MAD) with its statistics
  catalogue_his_rules.parquet  his catalogue rule: cc x nchannel >= threshold_perm (60 families; his six
                               eliminated families, threshold 0, keep their raw rows in catalogue_all)
Figures: duc_fix_monthly_counts.png (per-family monthly counts 2008-2026 + network total vs his),
         duc_fix_joins.png (our vs his monthly counts 2008-2011), duc_fix_old_vs_fixed.png.
Can be run at any time (partial results are labelled as such).
"""
import os, sys, json, glob, time
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from duc_fix_common import *
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUTMD = f"{RES}/DUC_EXT_FIXED_RESULTS.md"
BLUE, ORANGE, INK, INK2 = "#2a78d6", "#eb6834", "#0b0b0b", "#52514e"


def md(df, **kw):
    try: return df.to_markdown(**kw)
    except Exception: return "```\n" + df.to_string(**kw) + "\n```"


def final_months():
    ms = pd.period_range(T_START[:7], "2026-09", freq="M"); out = []
    for m in ms:
        p2, p1 = f"{FIX}/mf2/{m}.parquet", f"{FIX}/mf1/{m}.parquet"
        if os.path.exists(p2): out.append((m, p2, "pass2"))
        elif os.path.exists(p1): out.append((m, p1, "pass1"))
    return out


def main():
    thr = thresholds(); fp = families()
    lat = {}
    for ln in open(f"{DUC}/data/Plourde_2015/templates_list.txt").read().split("\n")[1:]:
        p = ln.split()
        if len(p) >= 5: lat[p[0]] = (float(p[2]), float(p[3]), float(p[4]))
    fams = sorted(fp.family, key=lambda f: -lat.get(f, (0, 0, 0))[0])
    L = ["# Ducellier & Creager (2022) southern-Cascadia families extended 2008 -> 2026 with the FIXED detector",
         "", f"_Auto-written by `scripts/duc_fix_report.py`, {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}. "
         "Queue `scripts/duc_ext_fixed_queue.sh`; all intermediate files in `results/duc_fix/`._", ""]
    # ---------------------------------------------------------------- stage 1
    S1 = pd.read_csv(f"{FIX}/repro_summary.csv") if os.path.exists(f"{FIX}/repro_summary.csv") else None
    ch = json.load(open(f"{FIX}/detector_choice.json")) if os.path.exists(f"{FIX}/detector_choice.json") else None
    M = final_months()
    det = []
    for m, p, src in M:
        D = pd.read_parquet(p); D["fam"] = D.fam.astype(str); D["src"] = src; det.append(D)
    D = pd.concat(det, ignore_index=True) if det else pd.DataFrame(columns=["fam", "t", "cc", "nchannel", "pass_his"])
    complete = len(M) == len(pd.period_range(T_START[:7], "2026-09", freq="M")) and os.path.exists(f"{RES}/DUC_FIX_DONE_PENDING")
    if len(D):
        D.to_parquet(f"{FIX}/catalogue_all.parquet")
        Dh = D[D.pass_his == 1]; Dh.to_parquet(f"{FIX}/catalogue_his_rules.parquet")
    else:
        Dh = D
    D["month"] = pd.to_datetime(D.t).dt.to_period("M"); Dh = D[D.pass_his == 1]
    H = his_catalogue(); H["thr"] = H.fam.map(thr).fillna(0.0)
    Hh = H[(H.thr > 0) & (H.cc * H.nch >= H.thr)].copy(); Hh["month"] = Hh.t.dt.to_period("M")
    # ---------------------------------------------------------------- epochs / re-templating
    E = epochs(); used = set(",".join(fp.stations).split(","))
    E = E[E.sta.isin(used)]
    TN = pd.read_csv(f"{FIX}/templates_new.csv", dtype={"fam": str}) if os.path.exists(f"{FIX}/templates_new.csv") else pd.DataFrame()
    dec = []
    for r in E.itertuples():
        if r.status == "his":
            chg = r.change if isinstance(r.change, str) else ""
            d_ = "HIS template kept (his epoch" + (f"; {chg}: his templates ran through this change in his catalogue, response removal only" if chg else "") + ")"
            k = n = np.nan
        else:
            t = TN[TN.epoch == r.epoch] if len(TN) else TN
            if not len(t): d_ = "new epoch: re-template PENDING"; k = n = np.nan
            else:
                k, n = int(t.keep.sum()), len(t)
                fk = t[t.keep == 1].fam.nunique(); fn = t.fam.nunique()
                d_ = (f"RE-TEMPLATED from our detections without it ({k}/{n} family-channels kept, {fk}/{fn} families)"
                      if k else f"DROPPED for this epoch (0/{n} family-channels pass the scrambled control)")
        dec.append(dict(epoch=r.epoch, start=r.start.date(), end=r.end.date(), loc=r.loc or "--", chans=r.chans,
                        sr=r.sr, sensor=r.sensor, status=r.status, change=r.change if isinstance(r.change, str) else "",
                        decision=d_))
    DEC = pd.DataFrame(dec); DEC.to_csv(f"{FIX}/station_epoch_decisions.csv", index=False)
    # ---------------------------------------------------------------- monthly counts
    months = pd.period_range("2008-01", "2026-09", freq="M")
    C = Dh.groupby(["fam", "month"]).size().unstack(fill_value=0).reindex(index=fams, columns=months, fill_value=0)
    Craw = D.groupby(["fam", "month"]).size().unstack(fill_value=0).reindex(index=fams, columns=months, fill_value=0)
    C.to_csv(f"{FIX}/monthly_counts_his_rules.csv"); Craw.to_csv(f"{FIX}/monthly_counts_raw.csv")
    Hm = Hh.groupby(["fam", "month"]).size().unstack(fill_value=0).reindex(index=fams, columns=months, fill_value=0)
    done_m = {m for m, _, _ in M}
    fig, ax = plt.subplots(2, 1, figsize=(15, 13), gridspec_kw=dict(height_ratios=[1, 4]), sharex=True)
    x = np.arange(len(months))
    tot = C.sum(0).values.astype(float); tot[[m not in done_m for m in months]] = np.nan
    htot = Hm.sum(0).values.astype(float); htot[months > pd.Period("2011-12", "M")] = np.nan
    ax[0].plot(x, tot, color=BLUE, lw=2, label="ours (fixed detector, his catalogue rule)")
    ax[0].plot(x, htot, color=ORANGE, lw=2, label="Ducellier & Creager 2022 (his rule)")
    ax[0].set_yscale("symlog", linthresh=10); ax[0].set_ylabel("detections / month\n(all families)", color=INK2)
    ax[0].legend(frameon=False, loc="upper right"); ax[0].grid(axis="y", color="#e6e6e3", lw=0.6)
    for s in ["top", "right"]: ax[0].spines[s].set_visible(False)
    Z = np.log10(1 + C.values.astype(float)); Z[:, [m not in done_m for m in months]] = np.nan
    im = ax[1].imshow(Z, aspect="auto", cmap="Blues", interpolation="nearest", vmin=0)
    ax[1].set_yticks(range(len(fams))); ax[1].set_yticklabels([f"{f} {lat.get(f, (0,))[0]:.2f}N" + (" *" if thr.get(f, 0) <= 0 else "") for f in fams], fontsize=6)
    yt = [i for i, m in enumerate(months) if m.month == 1]
    ax[1].set_xticks(yt); ax[1].set_xticklabels([str(months[i].year) for i in yt])
    ax[1].axvline([i for i, m in enumerate(months) if m == pd.Period("2011-11", "M")][0], color=INK2, lw=0.8, ls="--")
    ax[1].set_title("per-family monthly counts, his catalogue rule (cc x nchannel >= his threshold); families N->S; "
                    "* = one of his six eliminated families (threshold 0 -> empty); dashed = 2011-11 NC analogue->digital switch",
                    fontsize=8, color=INK2)
    cb = fig.colorbar(im, ax=ax[1], fraction=0.02, pad=0.01); cb.set_label("log10(1 + detections / month)")
    fig.tight_layout(); fig.savefig(f"{FIX}/duc_fix_monthly_counts.png", dpi=130); plt.close(fig)
    # ---------------------------------------------------------------- joins 2008-2011
    jm = pd.period_range("2008-01", "2011-12", freq="M")
    jm = [m for m in jm if m in done_m]
    J = []
    for f in fams:
        a = C.loc[f, jm].values if jm else np.array([]); b = Hm.loc[f, jm].values if jm else np.array([])
        ar = Craw.loc[f, jm].values if jm else np.array([])
        Hr = H[(H.fam == f)]; br = Hr.groupby(Hr.t.dt.to_period("M")).size().reindex(jm, fill_value=0).values if jm else np.array([])
        from scipy.stats import spearmanr
        rho = spearmanr(a, b).correlation if len(a) > 3 and a.std() > 0 and b.std() > 0 else np.nan
        rhor = spearmanr(ar, br).correlation if len(ar) > 3 and ar.std() > 0 and br.std() > 0 else np.nan
        J.append(dict(fam=f, lat=lat.get(f, (np.nan,))[0], thr_perm=thr.get(f, 0), ours_rule=int(a.sum()), his_rule=int(b.sum()),
                      ratio_rule=a.sum() / b.sum() if b.sum() else np.nan, spearman_rule=rho,
                      ours_raw=int(ar.sum()), his_raw=int(br.sum()), ratio_raw=ar.sum() / br.sum() if br.sum() else np.nan,
                      spearman_raw=rhor))
    JJ = pd.DataFrame(J); JJ.to_csv(f"{FIX}/joins_2008_2011.csv", index=False)
    JY = []
    for y in range(2008, 2012):
        mm = [m for m in jm if m.year == y]
        if not mm: continue
        JY.append(dict(year=y, months=len(mm), ours_rule=int(C[mm].values.sum()), his_rule=int(Hm[mm].values.sum()),
                       ours_raw=int(Craw[mm].values.sum()),
                       his_raw=int(H[H.t.dt.to_period("M").isin(mm)].shape[0])))
    JY = pd.DataFrame(JY)
    if len(JY):
        JY["ratio_rule"] = JY.ours_rule / JY.his_rule.replace(0, np.nan); JY["ratio_raw"] = JY.ours_raw / JY.his_raw.replace(0, np.nan)
    if jm:
        fig, ax = plt.subplots(figsize=(6.5, 6))
        a = C[jm].values.ravel(); b = Hm[jm].values.ravel(); k = (a + b) > 0
        ax.scatter(1 + b[k], 1 + a[k], s=9, color=BLUE, alpha=0.5, edgecolor="none")
        lim = [1, max(10, (1 + np.r_[a, b]).max()) * 1.5]
        ax.plot(lim, lim, color=INK2, lw=0.8, ls="--"); ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel("his detections per family-month + 1 (his rule)"); ax.set_ylabel("ours per family-month + 1 (his rule)")
        ax.set_title(f"joins 2008-2011: {len(jm)} months x {len(fams)} families", fontsize=9, color=INK2)
        for s in ["top", "right"]: ax.spines[s].set_visible(False)
        fig.tight_layout(); fig.savefig(f"{FIX}/duc_fix_joins.png", dpi=130); plt.close(fig)
    # ---------------------------------------------------------------- old vs fixed
    O = pd.read_csv(f"{RES}/catalog_Ducellier2022_per_family_year.csv", dtype={"fam": str})
    Oy = O.groupby("year")[["n_all", "n_core"]].sum()
    Fy = pd.DataFrame(dict(fixed_raw=D.groupby(pd.to_datetime(D.t).dt.year).size(),
                           fixed_his_rule=Dh.groupby(pd.to_datetime(Dh.t).dt.year).size()))
    Hy = pd.DataFrame(dict(his_raw=H.groupby(H.t.dt.year).size(), his_rule=Hh.groupby(Hh.t.dt.year).size()))
    OV = Oy.rename(columns={"n_all": "old_raw", "n_core": "old_core"}).join(Fy, how="outer").join(Hy, how="left")
    OV = OV.loc[2008:2026].fillna(0).astype(int)
    OV["inflation_raw"] = (OV.old_raw / OV.fixed_raw.replace(0, np.nan)).round(1)
    OV["inflation_core_vs_rule"] = (OV.old_core / OV.fixed_his_rule.replace(0, np.nan)).round(1)
    OV.to_csv(f"{FIX}/old_vs_fixed_per_year.csv")
    fig, ax = plt.subplots(figsize=(11, 4.5)); w = 0.2; xx = np.arange(len(OV))
    for i, (c, col, lab) in enumerate([("old_raw", "#86b6ef", "old extension, raw (unfixed t11_mf.py, our own templates)"),
                                        ("old_core", ORANGE, "old extension, 'core' tier"),
                                        ("fixed_raw", "#1c5cab", "fixed, raw 8 x hourly MAD"),
                                        ("fixed_his_rule", "#1baf7a", "fixed, his catalogue rule")]):
        ax.bar(xx + (i - 1.5) * w, OV[c].values + 1, w * 0.9, color=col, label=lab)
    ax.plot(xx, OV.his_rule.replace(0, np.nan) + 1, "o", color=INK, ms=5, label="his catalogue (his rule), to 2011")
    ax.set_yscale("log"); ax.set_xticks(xx); ax.set_xticklabels(OV.index, rotation=0, fontsize=8)
    ax.set_ylabel("detections / year + 1 (66 families)"); ax.legend(frameon=False, fontsize=7, ncol=2)
    for s in ["top", "right"]: ax.spines[s].set_visible(False)
    fig.tight_layout(); fig.savefig(f"{FIX}/duc_fix_old_vs_fixed.png", dpi=130); plt.close(fig)
    # ---------------------------------------------------------------- volume
    vol = {}
    for nm, pat in [("plans", f"{FIX}/data/plan_*.csv")]:
        P = [pd.read_csv(p) for p in glob.glob(pat)]
        if P:
            P = pd.concat(P); vol["station_days_planned"] = len(P)
            vol["downloaded_GB"] = round(P.bytes.sum() / 1e9, 1)
            vol["reused_from_traces"] = int((P.fetch == "have").sum())
            vol["nodata"] = int((P.fetch == "nodata").sum())
            if "cache" in P: vol["cache"] = P.cache.str.split(":").str[0].value_counts().to_dict()
    # ---------------------------------------------------------------- markdown
    nf_done = len(M); status = "COMPLETE" if nf_done == len(months) else f"PARTIAL ({nf_done}/{len(months)} months scanned)"
    kept_ep = DEC[DEC.decision.str.startswith("RE-TEMPLATED")].shape[0]; drop_ep = DEC[DEC.decision.str.startswith("DROPPED")].shape[0]
    L += ["## Plain-language summary", ""]
    if S1 is not None:
        s = S1.set_index("variant")
        L.append(f"* **Known answer first.** On Ducellier's own templates, stations and April 2008 window, the fixed detector "
                 f"(his hourly 8 x MAD, complete-hour channels) finds his events with mean per-family recall "
                 f"**{s.loc['full', 'mean_recall_all']:.3f}** with every fix, **{s.loc['nomute', 'mean_recall_all']:.3f}** with every fix but burst "
                 f"muting (old unfixed 12-h run: {s.loc['old_12h_unfixed', 'mean_recall_all']:.3f}); pooled precision "
                 f"{s.loc['full', 'pooled_precision']:.3f} / {s.loc['nomute', 'pooled_precision']:.3f} (old {s.loc['old_12h_unfixed', 'pooled_precision']:.3f}). "
                 + (f"Carried forward: **{ch['variant']}** ({ch['rule']})." if ch else "Stage 1 did not pass."))
    L.append(f"* **Scan status: {status}.** {len(Dh):,} detections pass his catalogue rule ({len(D):,} raw 8 x MAD detections kept with statistics).")
    L.append(f"* **Stations whose instruments changed** after his period (43 new epochs, mostly the 2011-10..12 NC analogue SHZ@50 Hz -> "
             f"digital EHZ@100 Hz switch) were re-templated from our own detections made without them: {kept_ep} epochs "
             f"re-templated, {drop_ep} dropped, the rest pending.")
    if len(JY):
        L.append(f"* **Joins with his catalogue 2008-2011** (his rule both sides): ours/his = " +
                 ", ".join(f"{int(r.year)} {r.ratio_rule:.2f}" for r in JY.itertuples()) + "; see the joins table.")
    if len(OV):
        L.append("* **The old extension was inflated**: old raw / fixed raw per year = " +
                 ", ".join(f"{y} x{v}" for y, v in OV.inflation_raw.items() if np.isfinite(v)) + ".")
    L += ["", "![monthly counts](duc_fix/duc_fix_monthly_counts.png)", ""]
    L += ["## 1. Known answer (2008-04-01..05-10, his templates + his stations, one-to-one +-2 s)", ""]
    if S1 is not None: L += [md(S1.round(3), index=False), ""]
    if os.path.exists(f"{FIX}/loo_retemplate_check.csv"):
        R = pd.read_csv(f"{FIX}/loo_retemplate_check.csv"); Rk = R[R.N >= 20]
        L += ["**Re-templating procedure checked on the known answer** (6 of his stations left out of the detector, re-templated "
              "from the detections made without them, compared with HIS template on the same channel):",
              f"{len(R)} family-channels, {len(Rk)} with >= 20 windows, {int(R.keep.sum())} kept by the scrambled control; "
              f"cc with his template: median {Rk.cc_his.median():.3f} (kept {R[R.keep == 1].cc_his.median():.3f}).", ""]
    L += ["## 2. Station epochs 2008 -> 2026-09 and what was done with each", "",
          "Epoch = constant (location, channels, sample rate, sensor model, horizontal azimuths) from FDSN StationXML. "
          "'his' = a key he used inside his catalogue period through one of his station rows; 'new' = anything else.", "",
          md(DEC, index=False), ""]
    L += ["## 3. Scan", "", f"Months scanned: {nf_done} of {len(months)}; pass-2 re-scans: {sum(1 for _, _, s in M if s == 'pass2')}.",
          f"Data: {json.dumps(vol)}", ""]
    L += ["## 4. Joins with his catalogue, 2008-2011", "", md(JY.round(3), index=False) if len(JY) else "(not yet scanned)", "",
          "Per family (his rule; raw = 8 x MAD without his cc x nchannel filter):", "", md(JJ.round(3), index=False) if len(JJ) else "", "",
          "![joins](duc_fix/duc_fix_joins.png)", ""]
    L += ["## 5. Old extension vs fixed", "",
          "Old = `results/ext/Ducellier2022_*.parquet` + `catalog_Ducellier2022*` (unfixed `t11_mf.py`, our own stacked templates on "
          "the 12 nearest modern stations, gap-normalisation + common-lag bugs) - SUPERSEDED, kept on disk.", "",
          md(OV.reset_index().rename(columns={"index": "year"}), index=False), "", "![old vs fixed](duc_fix/duc_fix_old_vs_fixed.png)", ""]
    open(OUTMD, "w").write("\n".join(L) + "\n")
    log(f"[report] {status}; wrote {OUTMD}")


if __name__ == "__main__":
    main()
