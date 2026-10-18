#!/usr/bin/env python
"""s11_wa_compare -- PLAN §5 extras for a region run (S11_REGION=WA_NW) after s11_v2 / s11_q3:
  * extra gates/columns: Q3 family split (odd/even), Q3 co-timed composite (WA majors within +-30 d of an SVI major),
    Q3 / Q1 / Q2 with the realign rule, Q1 / Q2 at n >= 50 and cc >= 0.8, single-site dominance (share of sum 1/sigma^2),
    per-station kernel-free annual amplitude outliers (|robust z| > 3);
  * SVI comparison (read-only on SVI tables; SVI Q1/Q2 recomputed at n >= 20 in a subprocess with S11_REGION=SVI):
    Q2 and Q3 composites side by side, calendar-year inverse-variance means of DeltaD, annual-cycle match;
  * figures fig_wa_*.png (region figures dir, copied to s11_depth/figures).
    S11_REGION=WA_NW python s11_wa_compare.py
    (internal) S11_REGION=SVI python s11_wa_compare.py --svi-part OUT.json
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, json, time, shutil, subprocess
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_invert as I
import s11_v2 as V
import s11_q3 as Q3

SVI_OUT = f"{C.ROOT}/svi_science/results/s11_depth"
MAINFIG = f"{SVI_OUT}/figures"
CC2 = 0.35


def demeaned_annual(cells):
    c = cells.copy(); c["dm"] = V.demeaned(c)
    return C.annual_fit(c.episode.map(V.U.mid), c.dm.to_numpy(), c.sig.to_numpy())


def svi_part(path):
    """runs with S11_REGION=SVI: read-only recomputation of SVI Q1/Q2 at n >= 20 (115 families, tag _all) + kernel-free
    pooled annual fit of the headline v2 and n >= 20 cells; nothing is written in the SVI folders."""
    out = {}
    for lab, tag, nmin in (("v2_n50_76fam", "", 50), ("n20_115fam", "_all", 20)):
        cells = I.prepare("raw", C.HEAD, tag, nmin=nmin, ccmin=CC2)
        q = V.Q(V.fit(cells, "smooth")); fa = demeaned_annual(cells)
        out[lab] = dict(Q=q, kf_amp=fa["amp"], kf_amp_sig=fa["amp_sig"], kf_peak=fa["peak_doy"], n_cells=len(cells))
        x, v = V.site_amp_curve(V.fit(cells, "smooth"), "")
        out[lab]["D_annual_curve"] = v.tolist()
    json.dump(out, open(path, "w"), indent=1, default=float)


def yearly(T, col="dD"):
    T = T.copy(); T["yr"] = pd.to_datetime(T.mid).dt.year; T["w"] = 1 / T.sig ** 2
    g = T.groupby("yr")
    return pd.DataFrame(dict(mean=g.apply(lambda q: np.average(q[col], weights=q.w), include_groups=False),
                             sig=g.w.sum() ** -0.5, n_eps=g.size(), wsum=g.w.sum())).reset_index()


def q3_composite(F, eps):
    gs = [Q3.dvec(F, e) for e in eps if ("D", "", e) in F["idx"]]
    gs = [g for g in gs if g is not None]
    if not gs:
        return np.nan, np.nan, 0
    gc = np.mean(gs, 0)
    return float(gc @ F["th"]), float(np.sqrt(gc @ F["C"] @ gc)), len(gs)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=14); ap.add_argument("--svi-part", default="")
    a = ap.parse_args()
    if a.svi_part:
        assert C.REGION == "SVI"; svi_part(a.svi_part); return
    assert C.REGION != "SVI"
    log = C.Log(f"{C.LOGD}/compare.log"); t0 = time.time()
    tag = C.REGION.lower().replace("_", "")[:4]                 # "wanw"
    pre = "fig_wa_" if C.REGION == "WA_NW" else f"fig_{C.REGION.lower()}_"
    rows = []
    # ---------------- headline objects
    cells = pd.read_csv(f"{C.TB}/q3_cells_b_n20.csv")
    F = Q3.fit(cells); T, comp = Q3.q3(F)
    MJ = pd.read_csv(f"{C.TB}/maj_list.csv")
    maj = list(MJ.episode.astype(str))
    rows.append(dict(item="Q3_composite", value=comp["composite"], sig=comp["composite_sig"], n=comp["n_majors"]))
    cot = list(MJ.episode[MJ.cotimed_svi.fillna("").astype(str) != ""].astype(str))
    v, s, n = q3_composite(F, cot)
    rows.append(dict(item="Q3_composite_cotimed_svi", value=v, sig=s, n=n))
    log(f"Q3 composite {comp['composite']:+.4f} ± {comp['composite_sig']:.4f} % ({comp['n_majors']} majors); co-timed {v:+.4f} ± {s:.4f} ({n})")
    # ---------------- Q3 family split
    To = Q3.q3(Q3.fit(cells[cells.fam % 2 == 1].reset_index(drop=True)))
    Te = Q3.q3(Q3.fit(cells[cells.fam % 2 == 0].reset_index(drop=True)))
    jo, je = To[0].set_index("episode"), Te[0].set_index("episode")
    com = [e for e in maj if e in jo.index and e in je.index]
    ok = [abs(jo.dD[e] - je.dD[e]) < 2 * np.hypot(jo.sig[e], je.sig[e]) for e in com]
    co, ce = To[1], Te[1]
    cpass = abs(co["composite"] - ce["composite"]) < 2 * np.hypot(co["composite_sig"], ce["composite_sig"]) if co and ce else None
    frac = float(np.mean(ok)) if ok else np.nan
    rows.append(dict(item="Q3_family_split", value=frac, n=len(com), odd=co["composite"] if co else np.nan,
                     odd_sig=co["composite_sig"] if co else np.nan, even=ce["composite"] if ce else np.nan,
                     even_sig=ce["composite_sig"] if ce else np.nan, passed=bool(cpass and frac >= 0.9) if ok else None))
    log(f"Q3 family split: {frac:.2f} of {len(com)} majors within 2σ; composites odd {rows[-1]['odd']:+.4f} ± {rows[-1]['odd_sig']:.4f}, "
        f"even {rows[-1]['even']:+.4f} ± {rows[-1]['even_sig']:.4f} → passed {rows[-1]['passed']}")
    pd.concat([jo.assign(half="odd"), je.assign(half="even")]).to_csv(f"{C.TB}/q3_family_split.csv")
    # ---------------- columns: realign rule, n >= 50, cc >= 0.8
    for lab, kw in (("realign", dict(nmin=20, realign_rule=True, ccmin=CC2)), ("n50", dict(nmin=50, ccmin=CC2)),
                    ("cc08_n20", dict(nmin=20, ccmin=0.8))):
        try:
            c = I.prepare("raw", C.HEAD, "", **kw)
            q = V.Q(V.fit(c, "smooth"))
            Fc = Q3.fit(c); vv, ss, nn = q3_composite(Fc, maj)
            rows.append(dict(item=f"column_{lab}", n_cells=len(c), Q1=q["Q1"], Q1_sig=q["Q1_sig"], Q1_peak=q["Q1_peak"], Q2=q["Q2"],
                             Q2_sig=q["Q2_sig"], value=vv, sig=ss, n=nn))
            log(f"column {lab}: cells {len(c)}, Q1 {q['Q1']:.4f} ± {q['Q1_sig']:.4f} (day {q['Q1_peak']:.0f}), Q2 {q['Q2']:+.4f} ± "
                f"{q['Q2_sig']:.4f}, Q3 composite {vv:+.4f} ± {ss:.4f} ({nn} majors)")
        except Exception as e:
            log(f"column {lab}: ERROR {e}")
    # ---------------- single-site dominance
    v2c = pd.read_csv(f"{C.TB}/v2_cells.csv")
    dom = []
    for nm, cc_ in (("q3_n20", cells), ("v2", v2c)):
        w = (1 / cc_.sig ** 2).groupby(cc_.site).sum(); sh = (w / w.sum()).sort_values(ascending=False)
        ncell = cc_.groupby("site").size()
        for s_, x in sh.items():
            dom.append(dict(cells=nm, site=s_, weight_share=float(x), n_cells=int(ncell[s_]), flag=bool(x > 0.5)))
    DOM = pd.DataFrame(dom); DOM.to_csv(f"{C.TB}/site_dominance.csv", index=False)
    top = DOM[DOM.cells == "q3_n20"].iloc[0]
    rows.append(dict(item="site_dominance_max", value=top.weight_share, site=top.site, passed=not bool(top.flag)))
    log(f"site dominance (Q3 cells): top {top.site} {100 * top.weight_share:.1f} % of Σ1/σ²; flagged > 50 %: "
        f"{DOM[DOM.flag].site.tolist()}")
    STA = pd.read_csv(f"{C.TB}/v2_station_amp.csv")
    med = STA.A.median(); mad = 1.4826 * np.median(np.abs(STA.A - med))
    STA["robust_z"] = (STA.A - med) / mad if mad > 0 else np.nan; STA["flag"] = STA.robust_z.abs() > 3
    STA.to_csv(f"{C.TB}/station_amp_outliers.csv", index=False)
    rows.append(dict(item="station_amp_outliers", value=int(STA.flag.sum()), site=";".join(STA.netsta[STA.flag])))
    log(f"station annual-amplitude outliers (|robust z| > 3): {STA.netsta[STA.flag].tolist()}")
    # ---------------- SVI comparison
    sj = f"{C.TB}/svi_recomputed.json"
    if not os.path.exists(sj):
        env = dict(os.environ, S11_REGION="SVI")
        subprocess.run([sys.executable, os.path.abspath(__file__), "--svi-part", sj], env=env, check=True)
    SV = json.load(open(sj))
    sm = pd.read_csv(f"{SVI_OUT}/tables/v2_models.csv"); sh_ = sm[sm.label == "head"].iloc[0]
    sq3 = json.load(open(f"{SVI_OUT}/tables/q3_headline.json"))["composite"]
    sT = pd.read_csv(f"{SVI_OUT}/tables/q3_dD_headline.csv")
    wm = pd.read_csv(f"{C.TB}/v2_models.csv"); wh = wm[wm.label == "head"].iloc[0]
    wv2 = V.Q(V.fit(v2c, "smooth")); wkf = demeaned_annual(v2c)
    side = [dict(quantity="Q1 deep annual (smooth, headline)", SVI=f"{sh_.Q1:.3f} ± {sh_.Q1_sig:.3f} % (n≥50), day {sh_.Q1_peak:.0f}",
                 SVI_n20=f"{SV['n20_115fam']['Q']['Q1']:.3f} ± {SV['n20_115fam']['Q']['Q1_sig']:.3f} %, day {SV['n20_115fam']['Q']['Q1_peak']:.0f}",
                 WA=f"{wh.Q1:.3f} ± {wh.Q1_sig:.3f} % (n≥20), day {wh.Q1_peak:.0f}"),
            dict(quantity="Q2 decade contrast of D", SVI=f"{sh_.Q2:+.3f} ± {sh_.Q2_sig:.3f} % (2015–25 − 2004–14)",
                 SVI_n20=f"{SV['n20_115fam']['Q']['Q2']:+.3f} ± {SV['n20_115fam']['Q']['Q2_sig']:.3f} %",
                 WA=f"{wh.Q2:+.3f} ± {wh.Q2_sig:.3f} % (2015–25 − 2007–14)"),
            dict(quantity="Q3 composite of majors", SVI=f"{sq3['composite']:+.3f} ± {sq3['composite_sig']:.3f} % (15 majors)", SVI_n20="same",
                 WA=f"{comp['composite']:+.3f} ± {comp['composite_sig']:.3f} % ({comp['n_majors']} majors); co-timed {v:+.3f} ± {s:.3f} ({n})"),
            dict(quantity="kernel-free pooled annual (series-demeaned cells)", SVI=f"{SV['v2_n50_76fam']['kf_amp']:.3f} ± {SV['v2_n50_76fam']['kf_amp_sig']:.3f} %, day {SV['v2_n50_76fam']['kf_peak']:.0f}",
                 SVI_n20=f"{SV['n20_115fam']['kf_amp']:.3f} ± {SV['n20_115fam']['kf_amp_sig']:.3f} %, day {SV['n20_115fam']['kf_peak']:.0f}",
                 WA=f"{wkf['amp']:.3f} ± {wkf['amp_sig']:.3f} %, day {wkf['peak_doy']:.0f}")]
    pd.DataFrame(side).to_csv(f"{C.TB}/compare_svi_side_by_side.csv", index=False)

    def match(A1, s1, p1, A2, s2, p2):
        dp = abs((p1 - p2 + 182.6) % 365.25 - 182.6)
        return dict(peak_diff_d=float(dp), amp_diff=float(A1 - A2), amp_diff_sig=float(np.hypot(s1, s2)),
                    match=bool(dp <= 60 and abs(A1 - A2) < 2 * np.hypot(s1, s2)))
    mt = [dict(what="Q1 WA vs SVI n20", **match(wh.Q1, wh.Q1_sig, wh.Q1_peak, SV['n20_115fam']['Q']['Q1'], SV['n20_115fam']['Q']['Q1_sig'],
                                                SV['n20_115fam']['Q']['Q1_peak'])),
          dict(what="Q1 WA vs SVI headline", **match(wh.Q1, wh.Q1_sig, wh.Q1_peak, sh_.Q1, sh_.Q1_sig, sh_.Q1_peak)),
          dict(what="kernel-free WA vs SVI n20", **match(wkf["amp"], wkf["amp_sig"], wkf["peak_doy"], SV['n20_115fam']['kf_amp'],
                                                         SV['n20_115fam']['kf_amp_sig'], SV['n20_115fam']['kf_peak'])),
          dict(what="kernel-free WA vs SVI headline", **match(wkf["amp"], wkf["amp_sig"], wkf["peak_doy"], SV['v2_n50_76fam']['kf_amp'],
                                                             SV['v2_n50_76fam']['kf_amp_sig'], SV['v2_n50_76fam']['kf_peak']))]
    MT = pd.DataFrame(mt); MT.to_csv(f"{C.TB}/compare_annual_match.csv", index=False)
    for r in mt:
        log("annual match " + ", ".join(f"{k} {x:.4f}" if isinstance(x, float) else f"{k} {x}" for k, x in r.items()))
    # calendar-year means of DeltaD (all episodes, and majors only)
    YW = yearly(T); YS = yearly(sT)
    YWm = yearly(T[T.episode.isin(maj)]); YSm = yearly(sT[sT.major15])
    Y = pd.concat([YW.assign(region="WA_NW", set="all"), YS.assign(region="SVI", set="all"),
                   YWm.assign(region="WA_NW", set="majors"), YSm.assign(region="SVI", set="majors")])
    Y.to_csv(f"{C.TB}/compare_yearly_dD.csv", index=False)
    j = YW.merge(YS, on="yr", suffixes=("_wa", "_svi"))
    r_y = float(np.corrcoef(j.mean_wa, j.mean_svi)[0, 1]) if len(j) > 3 else np.nan
    chi_w = float(np.mean((YW["mean"] / YW.sig) ** 2)); chi_s = float(np.mean((YS["mean"] / YS.sig) ** 2))
    rows.append(dict(item="yearly_dD_corr_WA_SVI", value=r_y, n=len(j), chi2_wa=chi_w, chi2_svi=chi_s))
    log(f"calendar-year ΔD: {len(YW)} WA years, {len(YS)} SVI years; corr over {len(j)} common years {r_y:.2f}; "
        f"mean (ΔD/σ)² WA {chi_w:.2f}, SVI {chi_s:.2f}")
    RW = pd.DataFrame(rows); RW.to_csv(f"{C.TB}/compare_extras.csv", index=False)
    figures(T, sT, maj, Y, SV, wv2, v2c, DOM, pre)
    # copies of the region figures
    os.makedirs(MAINFIG, exist_ok=True)
    for src, dst in (("fig_q3_dD.png", "q3_dD.png"), ("fig_v2_annual.png", "annual.png")):
        if os.path.exists(f"{C.FIG}/{src}"):
            shutil.copy(f"{C.FIG}/{src}", f"{MAINFIG}/{pre}{dst}")
    for f in ("compare.png", "sites.png"):
        shutil.copy(f"{C.FIG}/{pre}{f}", f"{MAINFIG}/{pre}{f}")
    log(f"compare done {time.time() - t0:.0f}s")
    open(f"{C.OUT}/WA_COMPARE_DONE", "w").write(C.pt_now() + "\n")


def figures(T, sT, maj, Y, SV, wv2, v2c, DOM, pre):
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    plt.rcParams.update({"axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#e6e6e6",
                         "grid.linewidth": 0.5, "font.size": 9})
    CW, CS = "#eb6834", "#2a78d6"
    fig, ax = plt.subplots(3, 1, figsize=(15, 12), gridspec_kw=dict(height_ratios=[1.3, 1, 1]))
    # (a) DeltaD of majors, both regions; points coloured by weight
    allw = np.r_[1 / T.sig ** 2, 1 / sT.sig ** 2]; vmax = np.quantile(allw, 0.95)
    for TT, mask, col, mk, lab in ((sT, sT.major15, CS, "s", "SVI (15 majors)"), (T, T.episode.isin(maj), CW, "o", f"WA-NW ({len(maj)} majors)")):
        q = TT[mask].copy(); q["mid"] = pd.to_datetime(q.mid); q = q.sort_values("mid")
        ax[0].plot(q.mid, q.dD, "-", color=col, lw=1.0, alpha=0.8, label=f"{lab}: mean {q.dD.mean():+.3f} %")
        sc = ax[0].scatter(q.mid, q.dD, c=1 / q.sig ** 2, cmap="viridis", vmin=0, vmax=vmax, marker=mk, s=40, edgecolor=col, lw=1.0, zorder=3)
    ax[0].axhline(0, color="grey", lw=0.6, ls=":")
    ax[0].set(ylabel="ΔD (%): deep-weighted dv/v of the episode\nminus neighbouring episodes (±3 yr)",
              title="Q3 per major slow-slip episode, SVI (squares, blue line) and WA-NW (circles, orange line); fill = weight 1/σ² (σ in tables, not drawn)")
    ax[0].legend(fontsize=8, loc="upper left"); fig.colorbar(sc, ax=ax[0], pad=0.01).set_label("weight 1/σ² (%⁻²)")
    # (b) calendar-year means
    for reg, col, mk in (("SVI", CS, "s"), ("WA_NW", CW, "o")):
        q = Y[(Y.region == reg) & (Y.set == "all")]
        ax[1].plot(q.yr, q["mean"], "-", color=col, lw=1.2, label=f"{reg}: all episodes")
        sc2 = ax[1].scatter(q.yr, q["mean"], c=q.wsum, cmap="viridis", vmin=0, vmax=np.quantile(Y.wsum, 0.95), marker=mk, s=45,
                            edgecolor=col, lw=1.0, zorder=3)
    ax[1].axhline(0, color="grey", lw=0.6, ls=":")
    ax[1].set(ylabel="calendar-year mean ΔD (%)", title="calendar-year inverse-variance means of ΔD (fill = Σ1/σ² of the year)")
    ax[1].legend(fontsize=8); fig.colorbar(sc2, ax=ax[1], pad=0.01).set_label("Σ weight")
    # (c) annual terms
    x = np.arange(1, 366)
    Fw = V.fit(v2c, "smooth"); _, vw = V.site_amp_curve(Fw, "")
    ax[2].plot(x, vw, "-", color=CW, lw=2.2, label=f"WA-NW deep D annual (smooth, n≥20): {wv2['Q1']:.3f} %, day {wv2['Q1_peak']:.0f}")
    q = SV["n20_115fam"]
    ax[2].plot(x, q["D_annual_curve"], "-", color=CS, lw=2.2, label=f"SVI deep D annual (smooth, n≥20): {q['Q']['Q1']:.3f} %, day {q['Q']['Q1_peak']:.0f}")
    q = SV["v2_n50_76fam"]
    ax[2].plot(x, q["D_annual_curve"], "--", color=CS, lw=1.2, label=f"SVI headline (v2 n≥50): {q['Q']['Q1']:.3f} %, day {q['Q']['Q1_peak']:.0f}")
    ax[2].axhline(0, color="grey", lw=0.6, ls=":")
    ax[2].set(xlabel="day of year", ylabel="dv/v (%)", title="annual term of the deep layer D(t) in each region")
    ax[2].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{C.FIG}/{pre}compare.png", dpi=130); plt.close(fig)
    # sites
    fig, ax = plt.subplots(1, 2, figsize=(15, 5.5))
    q = DOM[DOM.cells == "q3_n20"].sort_values("weight_share", ascending=False).head(25)
    ax[0].barh(q.site[::-1], 100 * q.weight_share[::-1], color=CW, height=0.7)
    ax[0].axvline(50, color="#52514e", lw=0.8, ls="--")
    for i, r in enumerate(q[::-1].itertuples()):
        ax[0].annotate(f"{r.n_cells} cells", (100 * r.weight_share, i), xytext=(3, -3), textcoords="offset points", fontsize=7)
    ax[0].set(xlabel="share of total weight Σ1/σ² (%)", title="WA-NW Q3 cells: weight share per site (top 25; dashed = 50 % flag)")
    ST = pd.read_csv(f"{C.TB}/station_amp_outliers.csv")
    w = 1 / ST.A_sig ** 2
    sc = ax[1].scatter(ST.a_d, ST.A, c=w / w.max(), cmap="viridis", vmin=0, vmax=1, s=50, edgecolor="k", lw=0.4, zorder=3)
    for r in ST.itertuples():
        if r.flag or r.n >= ST.n.quantile(0.75):
            ax[1].annotate(r.netsta.split(".")[1] + (" *" if r.flag else ""), (r.a_d, r.A), xytext=(4, 3), textcoords="offset points", fontsize=7)
    o = np.argsort(ST.a_d.to_numpy())
    ax[1].plot(ST.a_d.to_numpy()[o], ST.pred_a_d.to_numpy()[o], "-", color="k", lw=1.2, label="best fit ∝ a_d (deep)")
    ax[1].plot(ST.a_d, ST.pred_a_s, "x", color=CW, ms=7, mew=1.4, label="best fit ∝ a_s (shallow), per station")
    ax[1].set(xlabel="mean deep coefficient a_d of the station's cells", ylabel="in-phase annual dv/v amplitude (%)",
              title="WA-NW per-station kernel-free annual amplitude (* = |robust z| > 3)")
    ax[1].legend(fontsize=7); fig.colorbar(sc, ax=ax[1], label="relative weight 1/σ²")
    fig.tight_layout(); fig.savefig(f"{C.FIG}/{pre}sites.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
