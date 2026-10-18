"""Q8b: D2 RTR finder on PNSN tremor 2010-2026 (SVI, WA-NW, WA-SW boxes). Rules: q8_rtr/PLAN.md section Q8b."""
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_data as qd
import q8_rtr_lib as L
import q8_rtr_backzone as BZ
from q8_rtr_analyze import match

OUT = f"{L.OUT}/tables"
RNG = np.random.default_rng(88)
RES = {}
CATS = {"pnsn_svi": "SVI", "pnsn_nw": "NW", "pnsn_sw": "SW"}


def main():
    S = pd.read_parquet(f"{OUT}/streaks_all_pnsn_d2.parquet")
    M = pd.read_parquet(f"{OUT}/members_data_pnsn_d2.parquet")
    LD = {c: qd.pnsn(r) for c, r in CATS.items()}
    # absolute times
    for c, (ev, E) in LD.items():
        Ei = E.set_index("ep"); m = S.cat == c
        S.loc[m, "w0"] = S.loc[m, "ep"].map(Ei.w0).values; S.loc[m, "onset"] = S.loc[m, "ep"].map(Ei.onset).values
    S["w0"] = pd.to_datetime(S.w0); S["onset"] = pd.to_datetime(S.onset)
    S["start"] = S.w0 + pd.to_timedelta(S.t_first, unit="h"); S["end"] = S.w0 + pd.to_timedelta(S.t_last, unit="h")
    S["dur_h"] = S.t_last - S.t_first; S["days_after_onset"] = (S.start - S.onset).dt.total_seconds() / 86400
    S.to_parquet(f"{OUT}/streaks_all_abs_pnsn_d2.parquet", index=False)
    D = S[S["mode"] == "data"]
    rtr = lambda c: D[(D.cat == c) & (D.cls == "RTR")].copy()

    # ---------------- controls + a (counts, timing, kinematics)
    rows = []
    for c, (ev, E) in LD.items():
        d = D[D.cat == c]
        n1 = S[(S.cat == c) & (S["mode"] == "c1")]; n2 = S[(S.cat == c) & (S["mode"] == "c2")]
        R = rtr(c); nr = M[M.uid.isin(R.uid)].row.nunique()
        early = R[R.start < pd.Timestamp("2014-01-01")]; E_early = E[E.w0 < pd.Timestamp("2014-01-01")]
        rows.append(dict(box=CATS[c], n_eps=len(E), events=int(E.n.sum()), n_rtr=len(R), eps_with_rtr=R.ep.nunique(),
                         per_ep=len(R) / len(E), n_forward=int((d.cls == "forward").sum()), n_surge_only=int((d.cls == "surge_only").sum()),
                         c1_rtr=(n1.cls == "RTR").sum() / 5, c1_fwd=(n1.cls == "forward").sum() / 5,
                         c2_rtr=(n2.cls == "RTR").sum() / 5, c2_fwd=(n2.cls == "forward").sum() / 5,
                         frac_events_in_rtr=nr / E.n.sum(),
                         days_after_onset_med=R.days_after_onset.median(), speed_med=R.v_ev.abs().median(),
                         span_med=R.span_ev.median(), dur_med=R.dur_h.median(),
                         p_back_vs_fwd_binom=stats.binomtest(len(R), len(R) + int((d.cls == "forward").sum()), 0.5).pvalue
                         if len(R) + (d.cls == "forward").sum() > 0 else np.nan))
    A = pd.DataFrame(rows); A.to_csv(f"{OUT}/q8b_summary.csv", index=False); print(A.T.to_string())
    RES["summary"] = A.to_dict(orient="records")
    pe = []
    for c, (ev, E) in LD.items():
        R = rtr(c)
        for r in E.itertuples():
            q = R[R.ep == r.ep]
            pe.append(dict(box=CATS[c], ep=r.ep, s3_big=r.s3_big, n_windows=r.n, n_rtr=len(q),
                           days_after_onset=";".join(f"{v:.1f}" for v in sorted(q.days_after_onset)),
                           starts_PT=";".join(pd.Series(q.start.sort_values()).dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.strftime("%Y-%m-%d %H:%M").tolist())))
    PE = pd.DataFrame(pe); PE.to_csv(f"{OUT}/q8b_per_episode.csv", index=False)
    cols = ["cat", "ep", "start", "end", "days_after_onset", "dur_h", "v_ev", "span_ev", "x_start", "x_end", "front", "N_ev", "S"]
    pd.concat([rtr(c)[cols] for c in CATS]).to_csv(f"{OUT}/q8b_rtr_list.csv", index=False)

    # eras (one catalogue)
    era = {}
    for c in CATS:
        E = LD[c][1]; R = rtr(c)
        e1 = E.w0 < pd.Timestamp("2014-01-01")
        n1 = R.ep.isin(E[e1].ep).sum(); n2 = R.ep.isin(E[~e1].ep).sum(); k1, k2 = int(e1.sum()), int((~e1).sum())
        p = stats.binomtest(int(n1), int(n1 + n2), k1 / (k1 + k2)).pvalue if n1 + n2 > 0 else np.nan
        era[CATS[c]] = dict(eps_2010_13=k1, rtr_2010_13=int(n1), per_ep_2010_13=n1 / max(k1, 1), eps_2014_26=k2, rtr_2014_26=int(n2),
                            per_ep_2014_26=n2 / max(k2, 1), p_equal_rate=float(p))
    RES["eras"] = era; print("eras", era)
    # regions
    reg = {}
    for a_, b_ in [("pnsn_svi", "pnsn_nw"), ("pnsn_svi", "pnsn_sw"), ("pnsn_nw", "pnsn_sw")]:
        Ra, Rb = rtr(a_), rtr(b_)
        o = {}
        for col, f in [("speed", lambda d: d.v_ev.abs()), ("span", lambda d: d.span_ev), ("dur", lambda d: d.dur_h),
                       ("days_after_onset", lambda d: d.days_after_onset)]:
            if len(Ra) >= 3 and len(Rb) >= 3:
                o[col] = dict(med_a=float(f(Ra).median()), med_b=float(f(Rb).median()), p=float(stats.mannwhitneyu(f(Ra), f(Rb)).pvalue))
        ka, kb = len(LD[a_][1]), len(LD[b_][1])
        o["rate"] = dict(a=len(Ra) / ka, b=len(Rb) / kb, p=float(stats.binomtest(len(Ra), len(Ra) + len(Rb), ka / (ka + kb)).pvalue) if len(Ra) + len(Rb) else np.nan)
        reg[f"{CATS[a_]}_vs_{CATS[b_]}"] = o
    RES["regions"] = reg; print("regions", json.dumps(reg, indent=0))

    # ---------------- cross-check vs Bostock D2 RTRs and Royer (SVI, 2010-2013)
    B = pd.read_parquet(f"{OUT}/streaks_all_abs_d2.parquet")
    B = B[(B.cat == "bostock") & (B["mode"] == "data") & (B.cls == "RTR") & B.ep.isin(["2010-08", "2011-08", "2012-08", "2013-09"])]
    _, EB = qd.bostock(); EB = EB.set_index("ep")
    P = rtr("pnsn_svi")
    def overlap_rtr(Pl, Bl, sh=0.0):
        n = 0
        for ep in Bl.ep.unique():
            w0, w1 = EB.loc[ep, "w0"], EB.loc[ep, "w1"]
            a = Pl[(Pl.start >= w0) & (Pl.start < w1)].start.values + np.timedelta64(int(sh * 3600), "s")
            n += len(match(a, Bl[Bl.ep == ep].start.values, 6.0))
        return n
    nP = sum(((P.start >= EB.loc[e, "w0"]) & (P.start < EB.loc[e, "w1"])).sum() for e in B.ep.unique())
    mm = overlap_rtr(P, B)
    ch = [overlap_rtr(P, B, s) for s in (-72, -48, -24, 24, 48, 72)]
    xc = dict(n_bostock_rtr=len(B), n_pnsn_rtr_in_windows=int(nP), matched=mm, recall_of_bostock=mm / max(len(B), 1),
              precision_pnsn=mm / max(nP, 1), chance_matched_mean=float(np.mean(ch)), chance_recall=float(np.mean(ch)) / max(len(B), 1))
    Ro = pd.read_csv(f"{OUT}/royer2015_rtr_list.csv", parse_dates=["start_utc"])
    Ro = Ro[Ro.year.isin([2010, 2011, 2012])].rename(columns={"start_utc": "start"})
    Ro["ep"] = Ro.year.map({2010: "2010-08", 2011: "2011-08", 2012: "2012-08"})
    mr = overlap_rtr(P, Ro); chr_ = [overlap_rtr(P, Ro, s) for s in (-72, -48, -24, 24, 48, 72)]
    xc.update(royer_n=len(Ro), royer_matched=mr, royer_recall=mr / len(Ro), royer_chance=float(np.mean(chr_)) / len(Ro))
    # surge-level agreement: Bostock RTRs with any PNSN surge (any class) within +-6 h
    Ds = D[D.cat == "pnsn_svi"]
    xc["bostock_rtr_with_any_pnsn_surge_6h"] = float(np.mean([((Ds.start - b).abs() <= pd.Timedelta(hours=6)).any() for b in B.start]))
    xc["chance_any_surge"] = float(np.mean([np.mean([((Ds.start - b - pd.Timedelta(hours=s)).abs() <= pd.Timedelta(hours=6)).any() for b in B.start]) for s in (-48, -24, 24, 48)]))
    RES["crosscheck_2010_13"] = xc; print("xc", xc)

    # ---------------- b repeat sites
    bres = {}
    edges = np.arange(-80, 81, 10)
    for c, (ev, E) in LD.items():
        R = rtr(c).reset_index(drop=True)
        if R.ep.nunique() < 3:
            bres[CATS[c]] = dict(n=len(R), n_eps=int(R.ep.nunique())); continue
        pools = []
        for r in R.itertuples():
            g = ev[(ev.ep == r.ep) & ((ev.time - r.start).abs() <= pd.Timedelta(hours=12))]
            pools.append(g.x.values if len(g) else ev[ev.ep == r.ep].x.values)
        def nn(xx):
            ep = R.ep.values
            return np.mean([np.min(np.abs(xx[ep != ep[i]] - xx[i])) for i in range(len(xx))])
        def hc(xx):
            ep = R.ep.values
            H = np.array([np.histogram(xx[ep == e], edges)[0] for e in np.unique(ep)], float); H = H[H.sum(1) >= 2]
            if len(H) < 3:
                return np.nan
            C = np.corrcoef(H); iu = np.triu_indices(len(H), 1); return np.nanmean(C[iu])
        obs_nn, obs_hc = nn(R.x_start.values), hc(R.x_start.values)
        sims = [np.array([RNG.choice(p) for p in pools]) for _ in range(2000)]
        snn = np.array([nn(s) for s in sims]); shc = np.array([hc(s) for s in sims[:1000]])
        bres[CATS[c]] = dict(n=len(R), n_eps=int(R.ep.nunique()), nn_obs=float(obs_nn), nn_ctrl=float(snn.mean()),
                             p_nn=float((1 + (snn <= obs_nn).sum()) / 2001), hc_obs=float(obs_hc), hc_ctrl=float(np.nanmean(shc)),
                             p_hc=float((1 + (shc >= obs_hc).sum()) / 1001) if np.isfinite(obs_hc) else np.nan,
                             start_q10_50_90=R.x_start.quantile([.1, .5, .9]).round(1).tolist())
    RES["b_repeat"] = bres; print("b", bres)

    # ---------------- c start vs front
    cres = {}
    for c, (ev, E) in LD.items():
        R = rtr(c); pct, ds, de, dr = [], [], [], []
        for r in R.itertuples():
            if not np.isfinite(r.front):
                continue
            g = ev[ev.ep == r.ep]; t = ((g.time - r.w0).dt.total_seconds() / 3600).values; x = g.x.values
            on = BZ.onset_map(t, x)
            m = np.abs(t - r.t_first) <= 24
            onv = on.reindex(np.floor(x[m] / BZ.ONSET_BIN).astype(int)).values
            ref = np.abs(x[m][onv < t[m] - BZ.LAG_H] - r.front)
            if len(ref) < 20:
                continue
            d0 = abs(r.x_start - r.front); pct.append((ref < d0).mean()); ds.append(d0); de.append(abs(r.x_end - r.front)); dr.append(np.median(ref))
        if len(pct) >= 3:
            cres[CATS[c]] = dict(n=len(pct), start_km=float(np.median(ds)), end_km=float(np.median(de)), ref_km=float(np.median(dr)),
                                 mean_pct=float(np.mean(pct)), p=float(stats.wilcoxon(np.array(pct) - 0.5, alternative="less").pvalue),
                                 within5=float(np.mean(np.array(ds) <= 5)))
    RES["c_front"] = cres; print("c", cres)

    # ---------------- d tides (SVI)
    ev, E = LD["pnsn_svi"]; Ei = E.set_index("ep")
    tb = f"{L.ROOT}/svi_science/results/s5/tidal_stress_bostock_hourly.parquet"
    ts_ = f"{L.ROOT}/svi_science/results/s5_setup/tidal_stress_svi_hourly.parquet"
    def tide(Dset):
        res = []
        for r in Dset.itertuples():
            rows = M[M.uid == r.uid].row.values
            w = ev.loc[rows].tide_fid.value_counts()
            e = Ei.loc[r.ep]
            path = tb if e.w0 < pd.Timestamp("2014-01-15") else ts_
            T = pd.read_parquet(path, columns=["fid", "time", "tau_updip_Pa"],
                                filters=[("fid", "in", [int(f) for f in w.index]), ("time", ">=", e.w0), ("time", "<", e.w1 + pd.Timedelta(hours=2))])
            Pv = T.pivot(index="time", columns="fid", values="tau_updip_Pa"); Pv = Pv[[f for f in w.index if f in Pv.columns]]
            ww = w.reindex(Pv.columns).values.astype(float)
            ser = pd.Series(Pv.values @ ww / ww.sum(), index=Pv.index)
            grid = pd.date_range(e.w0, e.w1, freq="10min")
            tsv = np.interp(grid.astype("int64"), ser.index.astype("int64"), ser.values)
            n = max(int(round((r.end - r.start).total_seconds() / 600)), 1); i0 = int(np.searchsorted(grid, r.start))
            g = ev[ev.ep == r.ep]; days = g.time.dt.floor("D").value_counts(); act = days[days >= 0.01 * len(g)].index
            ok = np.isin(grid.floor("D"), act) & (np.arange(len(grid)) + n < len(grid))
            mins = pd.Series(tsv).rolling(n + 1).min().shift(-n).values; means = pd.Series(tsv).rolling(n + 1).mean().shift(-n).values
            res.append(dict(uid=r.uid, start=r.start, era="2010-13" if r.start < pd.Timestamp("2014-01-01") else "2014-26",
                            all_pos=bool((tsv[i0:i0 + n + 1] > 0).all()), mean_pos=bool(tsv[i0:i0 + n + 1].mean() > 0),
                            p_all=float(np.nanmean(mins[ok] > 0)), p_mean=float(np.nanmean(means[ok] > 0)), mean_udss=float(tsv[i0:i0 + n + 1].mean())))
        return pd.DataFrame(res)
    def summ(T):
        if len(T) == 0:
            return dict(n=0)
        s1 = (RNG.random((20000, len(T))) < T.p_all.values).sum(1); s2 = (RNG.random((20000, len(T))) < T.p_mean.values).sum(1)
        return dict(n=len(T), frac_all=float(T.all_pos.mean()), exp_all=float(T.p_all.mean()), p_all=float((1 + (s1 >= T.all_pos.sum()).sum()) / 20001),
                    frac_mean=float(T.mean_pos.mean()), exp_mean=float(T.p_mean.mean()), p_mean=float((1 + (s2 >= T.mean_pos.sum()).sum()) / 20001))
    dres = {}
    Dsvi = D[D.cat == "pnsn_svi"]
    for cls in ["RTR", "forward", "surge_only"]:
        T = tide(Dsvi[Dsvi.cls == cls]); T.to_csv(f"{OUT}/q8b_tides_{cls}.csv", index=False)
        dres[cls] = summ(T)
        if cls == "RTR":
            for e_ in ["2010-13", "2014-26"]:
                dres[f"RTR_{e_}"] = summ(T[T.era == e_])
    RES["d_tides_svi"] = dres; print("d", json.dumps(dres, indent=0))

    # ---------------- e up-dip (cells)
    eres = {}
    for c, (ev, E) in LD.items():
        R = rtr(c)
        if len(R) < 5:
            eres[CATS[c]] = dict(n=len(R)); continue
        isr = ev.index.isin(M[M.uid.isin(R.uid)].row.unique())
        n_all = ev.groupby("fam").size(); ok = n_all >= 50
        dd = ev.groupby("fam").dd.mean()
        n_r = ev[isr].groupby("fam").size().reindex(n_all.index, fill_value=0)
        rho = stats.spearmanr((n_r / n_all)[ok], dd[ok])[0]
        # random-time windows (same duration, slipped zone at that time, x range of the RTR)
        arr = {}
        for r in E.itertuples():
            g = ev[ev.ep == r.ep]; t = ((g.time - r.w0).dt.total_seconds() / 3600).values
            arr[r.ep] = dict(t=t, x=g.x.values, rows=g.index.values, T=(r.w1 - r.w0).total_seconds() / 3600, on=BZ.onset_map(t, g.x.values))
        rs = []
        for _ in range(200):
            rows_all = []
            for r in R.itertuples():
                a = arr[r.ep]; dur = r.t_last - r.t_first; xr = (min(r.x_start, r.x_end) - 3, max(r.x_start, r.x_end) + 3)
                for _t in range(20):
                    t0 = RNG.uniform(0, a["T"] - dur)
                    lo = np.searchsorted(a["t"], t0); hi = np.searchsorted(a["t"], t0 + dur); xx = a["x"][lo:hi]
                    onv = a["on"].reindex(np.floor(xx / BZ.ONSET_BIN).astype(int)).values
                    mk = (onv < t0 - BZ.LAG_H) & (xx >= xr[0]) & (xx <= xr[1])
                    if mk.sum() >= 5:
                        rows_all.append(a["rows"][lo:hi][mk]); break
            rr = np.concatenate(rows_all) if rows_all else np.array([], int)
            nr = ev.loc[rr].groupby("fam").size().reindex(n_all.index, fill_value=0)
            rs.append(stats.spearmanr((nr / n_all)[ok], dd[ok])[0])
        rs = np.array(rs)
        eres[CATS[c]] = dict(n_rtr=len(R), n_cells=int(ok.sum()), rho=float(rho), ctrl_mean=float(np.nanmean(rs)),
                             ctrl_95=[float(np.nanpercentile(rs, 2.5)), float(np.nanpercentile(rs, 97.5))],
                             p_more_updip=float((1 + (rs <= rho).sum()) / (len(rs) + 1)),
                             rtr_minus_all_downdip_km=float(ev[isr].dd.mean() - ev[ev.ep.isin(R.ep.unique())].dd.mean()))
    RES["e_updip"] = eres; print("e", eres)
    json.dump(RES, open(f"{L.OUT}/results_q8b.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
