#!/usr/bin/env python
"""r23_bursts (2026-10-11) -- R23 study 2: between-episode deep bursts 2003-2026 (finding 3; S7/Q5 rules).
Rules: svi_science/results/r23/PLAN.md §3.   python svi_science/scripts/r23_bursts.py
"""
from __future__ import annotations
import os, sys
import numpy as np
import pandas as pd
from scipy import stats
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import r23_common as C

OUT = f"{C.TB}/study2"; os.makedirs(OUT, exist_ok=True)
log = C.Log(f"{OUT}/study2.log")
RNG = np.random.default_rng(20261011)
NSH = 1000
SUBSEG = {"ERA": ("CN4f", "ERA"), "B2f+B3f": ("B2f", "B3f"), "S15f": ("S15f",), "FWD": ("FWD",)}
BURSTS2006 = [("2006-01-15", "2006-01-16"), ("2006-03-27", "2006-03-30"), ("2006-04-24", "2006-04-24"),
              ("2006-05-01", "2006-05-01"), ("2006-05-03", "2006-05-03"), ("2006-05-11", "2006-05-11"),
              ("2006-05-25", "2006-05-26"), ("2006-06-18", "2006-06-18")]


def in_period(df, p):
    a, b = C.PERIODS[p]; return (df.date >= a) & (df.date <= b)


def main():
    F = C.families(); E = C.episodes()
    V = C.load_events(F, E); D = C.load_famdays(F, E)
    dup = int(D.duplicated(["fam", "date"]).sum()); log(f"famday duplicates (fam,date): {dup}")
    maj = E[E.major]
    V["in_major"] = C.pad_mask(V.date, maj, 0)
    R = {}
    # ---------------- acceptance: S7 forward numbers
    s7 = V[in_period(V, "2014-26")]
    acc = dict(BE_data=int((s7.BE & (s7.series == "data")).sum()), BE_c1=int((s7.BE & (s7.series == "c1")).sum()),
               broad_data=int((s7.BE_broad & (s7.series == "data")).sum()))
    acc["pass"] = acc["BE_data"] == 76 and acc["BE_c1"] == 3 and acc["broad_data"] == 688
    R["acceptance_S7"] = acc; log(f"S7 acceptance (expect 76 / 3 / 688): {acc}")
    # ---------------- inventory per period / sub-segment
    inv = []
    units = [(p, None) for p in C.PERIODS] + [("2003-26", s) for s in SUBSEG]
    for p, s in units:
        m = in_period(V, p) if s is None else in_period(V, p) & V.seg.isin(SUBSEG[s])
        md = in_period(D, p) if s is None else in_period(D, p) & D.seg.isin(SUBSEG[s])
        x = V[m]
        row = dict(period=p, subseg=s or "all", elig_famdays_BE=int((md & D.BE).sum()), elig_famdays_broad=int((md & D.BE_broad).sum()))
        for nm in ("BE", "BE_broad", "BE_envsrc"):
            for ser in ("data", "c1"):
                xx = x[x[nm] & (x.series == ser)]
                row[f"{nm}_{ser}"] = int(len(xx))
                if nm == "BE":
                    row[f"BE_{ser}_no_burstday"] = int((~xx.burst_day_f).sum())
                    row[f"BE_{ser}_nz2ge2"] = int(((xx.n_z2.fillna(99)) >= 2).sum())
                    row[f"BE_{ser}_deep"] = int(xx.deep.sum())
        row["BE_per_1000_famdays"] = 1000 * row["BE_data"] / max(row["elig_famdays_BE"], 1)
        row["BE_excess_per_1000_famdays"] = 1000 * (row["BE_data"] - row["BE_c1"]) / max(row["elig_famdays_BE"], 1)
        row["exist"] = bool(row["BE_data"] >= 20 and row["BE_data"] >= 5 * max(row["BE_c1"], 0) and row["BE_data"] > 0
                            and (row["BE_c1"] == 0 or row["BE_data"] / row["BE_c1"] >= 5))
        inv.append(row)
    INV = pd.DataFrame(inv); INV.to_csv(f"{OUT}/inventory.csv", index=False)
    log("inventory (counts per eligible family-day are within-segment numbers, never trended):\n" + INV.to_string(index=False))
    R["inventory"] = INV.to_dict("records")

    # ---------------- per period: bursty, deep, bursts
    poolBE = C.pools(D, D.BE.to_numpy())
    BT = []
    for p, s in units:
        m = in_period(V, p) if s is None else in_period(V, p) & V.seg.isin(SUBSEG[s])
        tag = p if s is None else f"{p}:{s}"
        ev = V[m & V.BE & (V.series == "data")].copy()
        ets = V[m & V.in_major & (V.series == "data")]
        r = dict(n=int(len(ev)))
        if len(ev) < 5:
            R[f"tests_{tag}"] = r; continue
        # bursty
        keys = list(zip(ev.fam, ev.seg))
        sh = C.shuffle_dates(keys, poolBE, RNG, NSH)
        shd = np.array([len(set(x)) for x in sh])
        obs = ev.date.nunique()
        r.update(distinct_days=int(obs), shuffle_days_median=float(np.median(shd)), shuffle_days_p01=float(np.percentile(shd, 1)),
                 p_fewer_days=C.pct_lo(obs, shd), bursty=bool(obs < np.percentile(shd, 1)))
        # deep vs ETS shares (ETS events of the same period / sub-segment)
        sh_ets = ets.fam.value_counts(normalize=True)
        fams = sh_ets.index.to_numpy(); pr = sh_ets.to_numpy(); dep = F.dep.reindex(fams).to_numpy()
        draws = RNG.choice(len(fams), size=(10000, len(ev)), p=pr)
        dmean = dep[draws].mean(1)
        obsd = float(ev.dep.mean())
        r.update(mean_dep=obsd, expected_dep=float(dmean.mean()), expected_p99=float(np.percentile(dmean, 99)),
                 p_deeper=C.pct_hi(obsd, dmean), deeper=bool(obsd > np.percentile(dmean, 99)), n_ets=int(len(ets)),
                 share_dep_ge40=float((ev.dep >= 40).mean()), ets_share_dep_ge40=float((ets.dep >= 40).mean()),
                 share_deep_groups=float(ev.deep.mean()), ets_share_deep_groups=float(ets.deep.mean()))
        ev["burst"] = C.chain_bursts(ev.time_utc)
        bm = ev.groupby("burst").dep.mean()
        draws_b = RNG.choice(len(fams), size=(10000, len(bm)), p=pr)
        r.update(n_bursts=int(len(bm)), mean_dep_per_burst=float(bm.mean()), p_deeper_per_burst=C.pct_hi(float(bm.mean()), dep[draws_b].mean(1)))
        # scrambled control depth (listed)
        evc = V[m & V.BE & (V.series == "c1")]
        r.update(n_c1=int(len(evc)), c1_mean_dep=float(evc.dep.mean()) if len(evc) else np.nan)
        R[f"tests_{tag}"] = r
        log(f"[{tag}] {r}")
        # bursts table
        for b, g in ev.groupby("burst"):
            g = g.sort_values("time_utc")
            if len(g.fam.unique()) > 1:
                ff = F.loc[g.fam.unique()]
                dla = ff.lat.to_numpy(); dlo = ff.lon.to_numpy()
                sep = float(np.max(C_hav(dla[:, None], dlo[:, None], dla[None, :], dlo[None, :])))
            else:
                sep = 0.0
            BT.append(dict(unit=tag, period=p, subseg=s or "all", burst=int(b), start_PT=g.time_PT.iloc[0], end_PT=g.time_PT.iloc[-1],
                           n=len(g), size_k=float(g.k.sum()), dur_h=float((g.time_utc.iloc[-1] - g.time_utc.iloc[0]) / pd.Timedelta(hours=1)),
                           n_fam=int(g.fam.nunique()), max_sep_km=sep, dep_med=float(g.dep.median()),
                           groups=";".join(map(str, sorted(set(g.group)))), seg=";".join(sorted(set(g.seg))),
                           fams=";".join(map(str, sorted(set(g.fam)))), burst_day_any=bool((g.burst_day == True).any()),
                           n_z2_min=int(g.n_z2.min()) if g.n_z2.notna().any() else -1))
    BT = pd.DataFrame(BT); BT.to_csv(f"{OUT}/bursts.csv", index=False)

    # ---------------- burst-size comparisons (sub-segments vs FWD)
    cmp = []
    bf = BT[BT.unit == "2003-26:FWD"]
    for s in ("ERA", "B2f+B3f", "S15f"):
        bs = BT[BT.unit == f"2003-26:{s}"]
        for col in ("n", "size_k", "dur_h", "n_fam", "max_sep_km", "dep_med"):
            if len(bs) >= 3 and len(bf) >= 3:
                u = stats.mannwhitneyu(bs[col], bf[col], alternative="two-sided")
                cmp.append(dict(subseg=s, col=col, n_bursts=len(bs), n_bursts_fwd=len(bf), median=float(bs[col].median()),
                                median_fwd=float(bf[col].median()), mean=float(bs[col].mean()), mean_fwd=float(bf[col].mean()),
                                share_ge3=float((bs.n >= 3).mean()) if col == "n" else np.nan,
                                share_ge3_fwd=float((bf.n >= 3).mean()) if col == "n" else np.nan, p_MW=float(u.pvalue),
                                verdict="different" if u.pvalue < 0.01 else ("no detectable difference" if u.pvalue >= 0.05 else "hint")))
    CMP = pd.DataFrame(cmp); CMP.to_csv(f"{OUT}/burst_size_compare.csv", index=False)
    log("burst size comparison vs FWD:\n" + CMP.round(3).to_string(index=False))
    # 2003-13 pooled vs 2014-26
    b03 = BT[BT.unit == "2003-13"]; b14 = BT[BT.unit == "2014-26"]
    R["burst_pooled"] = {col: dict(med03=float(b03[col].median()), med14=float(b14[col].median()),
                                   p=float(stats.mannwhitneyu(b03[col], b14[col]).pvalue)) for col in ("n", "size_k", "dur_h", "n_fam", "max_sep_km", "dep_med")}
    log(f"burst pooled 2003-13 vs 2014-26: {R['burst_pooled']}")

    # ---------------- same patches: group shares
    G = []
    for p in ("2003-13", "2014-26"):
        ev = V[in_period(V, p) & V.BE & (V.series == "data")]
        evb = V[in_period(V, p) & V.BE_broad & (V.series == "data")]
        ets = V[in_period(V, p) & V.in_major & (V.series == "data")]
        for g in sorted(F.group.unique()):
            G.append(dict(period=p, group=g, dep=float(F[F.group == g].dep.mean()), BE=int((ev.group == g).sum()),
                          BE_share=float((ev.group == g).mean()), broad_share=float((evb.group == g).mean()),
                          ETS_share=float((ets.group == g).mean())))
    G = pd.DataFrame(G); G.to_csv(f"{OUT}/group_shares.csv", index=False)
    a = G[G.period == "2003-13"].set_index("group"); b = G[G.period == "2014-26"].set_index("group")
    out = {}
    for col in ("BE_share", "broad_share", "ETS_share"):
        rho = stats.spearmanr(a[col], b[col]).statistic
        perm = np.array([stats.spearmanr(a[col].to_numpy()[RNG.permutation(len(a))], b[col]).statistic for _ in range(10000)])
        out[col] = dict(rho=float(rho), p=C.pct_hi(rho, perm))
    # excess over ETS share agrees?
    ra = a.BE_share - a.ETS_share; rb = b.BE_share - b.ETS_share
    rho = stats.spearmanr(ra, rb).statistic
    perm = np.array([stats.spearmanr(ra.to_numpy()[RNG.permutation(len(ra))], rb).statistic for _ in range(10000)])
    out["BE_minus_ETS_share"] = dict(rho=float(rho), p=C.pct_hi(rho, perm))
    out["same_patches"] = bool(out["BE_share"]["rho"] > 0 and out["BE_share"]["p"] < 0.01)
    R["same_patches"] = out
    log("group shares:\n" + G.pivot(index="group", columns="period", values=["BE_share", "ETS_share", "broad_share"]).round(3).to_string())
    log(f"same patches: {out}")
    # top families per period
    top = {}
    for p in ("2003-13", "2014-26"):
        ev = V[in_period(V, p) & V.BE & (V.series == "data")]
        vc = ev.fam.value_counts().head(10)
        top[p] = [dict(fam=int(f), n=int(n), group=int(F.group[f]), dep=float(F.dep[f])) for f, n in vc.items()]
    R["top_families"] = top; log(f"top families: {top}")
    # per-family (BE+0.5)/(ETS+0.5) vs depth, per period
    for p in ("2003-13", "2014-26"):
        ev = V[in_period(V, p) & V.BE & (V.series == "data")]; ets = V[in_period(V, p) & V.in_major & (V.series == "data")]
        fe = pd.DataFrame(dict(BE=ev.fam.value_counts(), ETS=ets.fam.value_counts())).fillna(0)
        fe = fe[fe.ETS >= 10]
        rr = stats.spearmanr((fe.BE + .5) / (fe.ETS + .5), F.dep.reindex(fe.index))
        R[f"fam_ratio_vs_depth_{p}"] = dict(rho=float(rr.statistic), p=float(rr.pvalue), n_fam=int(len(fe)))

    # ---------------- early-2006 bursts
    rows = []
    q = pd.read_csv(f"{C.ROOT}/svi_science/results/backward/p4/envflag/tables/quiet_flag_env.csv", usecols=["date", "N_CN4", "N_all"])
    q = q.set_index("date")
    for a0, a1 in BURSTS2006:
        x = V[(V.date >= a0) & (V.date <= a1)]
        xd = x[x.series == "data"]; xc = x[x.series == "c1"]
        dd = pd.date_range(a0, a1).strftime("%Y-%m-%d")
        ep = E[(E.start <= a1) & (E.end >= a0)]
        rows.append(dict(dates=f"{a0}..{a1}", data=int(len(xd)), control=int(len(xc)), BE=int(xd.BE.sum()), BE_broad=int(xd.BE_broad.sum()),
                         in_listed_episode=";".join(ep.episode) if len(ep) else "", families=";".join(f"{f}({n})" for f, n in xd.fam.value_counts().head(8).items()),
                         groups=";".join(map(str, sorted(set(xd.group)))), dep_med=float(xd.dep.median()) if len(xd) else np.nan,
                         along_med=float(xd.along_km.median()) if len(xd) else np.nan,
                         env_N_CN4=";".join(str(q.N_CN4.get(d, np.nan)) for d in dd), env_N_all=";".join(str(q.N_all.get(d, np.nan)) for d in dd),
                         burst_day=bool((xd.burst_day == True).any()), first_PT=xd.time_PT.min() if len(xd) else "", last_PT=xd.time_PT.max() if len(xd) else ""))
    B6 = pd.DataFrame(rows); B6.to_csv(f"{OUT}/bursts_2006.csv", index=False)
    log("early-2006 bursts:\n" + B6.to_string(index=False))
    # all between-episode events list 2003-2013
    cols = ["time_PT", "date", "fam", "group", "dep", "along_km", "downdip_km", "P", "detector", "seg", "burst_day", "n_z2", "n30", "n30_env", "BE", "BE_broad"]
    V[(V.date <= "2013-12-31") & V.BE & (V.series == "data")][cols].to_csv(f"{OUT}/between_events_2003_2013.csv", index=False)
    V[(V.date <= "2013-12-31") & V.BE & (V.series == "c1")][cols].to_csv(f"{OUT}/between_events_2003_2013_control.csv", index=False)
    # ---------------- POST-HOC column (Revision 2): bursts of BE_broad events (between episodes, local tremor allowed)
    PB = []
    for s_ in SUBSEG:
        ev = V[V.seg.isin(SUBSEG[s_]) & V.BE_broad & (V.series == "data") & (V.date <= C.P_END)].copy()
        if len(ev) == 0:
            continue
        ev["burst"] = C.chain_bursts(ev.time_utc)
        g = ev.groupby("burst").agg(n=("fam", "size"), size_k=("k", "sum"), n_fam=("fam", "nunique"), dep=("dep", "median"),
                                    deep_share=("deep", "mean"), t0=("time_PT", "first"), groups=("group", lambda x: ";".join(map(str, sorted(set(x))))),
                                    n30=("n30", "max"))
        g["subseg"] = s_; PB.append(g.reset_index())
    PB = pd.concat(PB, ignore_index=True); PB.to_csv(f"{OUT}/posthoc_broad_bursts.csv", index=False)
    summ = PB.groupby("subseg").agg(n_bursts=("n", "size"), n_events=("n", "sum"), median_n=("n", "median"), max_n=("n", "max"),
                                    median_size_k=("size_k", "median"), max_size_k=("size_k", "max"), share_deep=("deep_share", "mean"),
                                    share_with_tremor30=("n30", lambda x: float((x > 0).mean())))
    log("POST-HOC broad bursts by sub-segment:\n" + summ.round(2).to_string())
    log("largest broad bursts:\n" + PB.sort_values("size_k", ascending=False).head(15).round(1).to_string(index=False))
    R["posthoc_broad_bursts"] = summ.reset_index().to_dict("records")
    for s_ in ("ERA", "B2f+B3f", "S15f"):
        a_, b_ = PB[PB.subseg == s_], PB[PB.subseg == "FWD"]
        R[f"posthoc_broad_MW_{s_}_vs_FWD"] = {c: float(stats.mannwhitneyu(a_[c], b_[c]).pvalue) for c in ("n", "size_k", "dep")}
    log(f"post-hoc MW: { {k: v for k, v in R.items() if k.startswith('posthoc_broad_MW')} }")
    C.dump(R, f"{OUT}/study2_summary.json")
    figures(V, BT, G, F, E)
    log("done")


def C_hav(lat1, lon1, lat2, lon2):
    p = np.deg2rad
    a = np.sin(p(lat2 - lat1) / 2) ** 2 + np.cos(p(lat1)) * np.cos(p(lat2)) * np.sin(p(lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(np.clip(a, 0, 1)))


def figures(V, BT, G, F, E):
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 1, figsize=(13, 7.5), gridspec_kw=dict(height_ratios=[2.2, 1]))
    ev = V[V.BE & (V.series == "data")]; evc = V[V.BE & (V.series == "c1")]
    t = pd.to_datetime(ev.date); tc = pd.to_datetime(evc.date)
    for _, r in E.iterrows():
        ax[0].axvspan(r.start_ts, r.end_ts + pd.Timedelta(days=1), color="#9ecae1" if r.major else "#deebf7", lw=0, zorder=0)
    sc = ax[0].scatter(t, ev.along_km, c=ev.dep, cmap="viridis_r", s=10, vmin=28, vmax=46, zorder=3, label="between-episode LFE")
    ax[0].scatter(tc, F.along_strike_km.reindex(evc.fam), marker="x", c="crimson", s=18, zorder=4, label="scrambled control")
    plt.colorbar(sc, ax=ax[0], label="Bostock depth (km)")
    for s, a in (("2006-07-01", "era→frozen"), ("2009-08-27", "S15"), ("2010-01-06", "env→PNSN tremor"), ("2014-01-01", "forward")):
        ax[0].axvline(pd.Timestamp(s), color="k", ls=":", lw=0.8); ax[0].text(pd.Timestamp(s), -52, " " + a, fontsize=7, rotation=90, va="bottom")
    ax[0].set_ylabel("along-strike (km, + NW)"); ax[0].legend(loc="lower left", fontsize=8)
    ax[0].set_title("SVI between-episode LFEs 2003–2026 (outside every listed episode ± 2 d, no tremor ≤ 30 km); blue = episodes")
    # group shares
    a = G[G.period == "2003-13"].set_index("group"); b = G[G.period == "2014-26"].set_index("group")
    x = np.arange(len(a)); w = 0.27
    ax[1].bar(x - w, a.BE_share, w, label="between 2003–13", color="#2c7fb8")
    ax[1].bar(x, b.BE_share, w, label="between 2014–26", color="#7fcdbb")
    ax[1].bar(x + w, b.ETS_share, w, label="ETS share 2014–26", color="#bdbdbd")
    ax[1].set_xticks(x); ax[1].set_xticklabels([f"G{g}\n{d:.0f} km" for g, d in zip(a.index, a.dep)], fontsize=7)
    ax[1].set_ylabel("share"); ax[1].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(f"{C.FIG}/r23_between_timeline.png", dpi=130); plt.close(fig)
    # burst sizes
    fig, ax = plt.subplots(1, 3, figsize=(13, 3.8))
    for u, col in (("2003-26:ERA", "#08519c"), ("2003-26:B2f+B3f", "#6baed6"), ("2003-26:S15f", "#31a354"), ("2003-26:FWD", "#fd8d3c")):
        bb = BT[BT.unit == u]
        if len(bb) == 0:
            continue
        for k, c in enumerate(("n", "size_k", "dur_h")):
            v = np.sort(bb[c].to_numpy())[::-1]
            ax[k].step(v, np.arange(1, len(v) + 1) / len(v), where="post", color=col, label=f"{u.split(':')[1]} ({len(v)})")
    for k, (c, lab) in enumerate((("n", "events per burst"), ("size_k", "Bostock-equivalent size (Σk)"), ("dur_h", "duration (h)"))):
        ax[k].set_xscale("symlog"); ax[k].set_yscale("log"); ax[k].set_xlabel(lab); ax[k].set_ylabel("share of bursts ≥ x")
    ax[0].legend(fontsize=8)
    fig.suptitle("Between-episode burst sizes by network segment (bursts = gaps < 24 h)")
    fig.tight_layout(); fig.savefig(f"{C.FIG}/r23_between_burst_sizes.png", dpi=130); plt.close(fig)


if __name__ == "__main__":
    main()
