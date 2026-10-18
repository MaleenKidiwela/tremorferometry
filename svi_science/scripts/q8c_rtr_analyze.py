"""Q8c analysis: D2 RTRs on all candidates. Rules: q8_rtr/PLAN.md section Q8c."""
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
RNG = np.random.default_rng(808)
RES = {}
LOAD = {"ac_svi": lambda v, s: qd.allcand_svi(v, s), "ac_back": lambda v, s: qd.back_2010(v, s),
        "ac_nw": lambda v, s: qd.allcand_wa("NW", v, s), "ac_sw": lambda v, s: qd.allcand_wa("SW", v, s)}
TB = f"{L.ROOT}/svi_science/results/s5/tidal_stress_bostock_hourly.parquet"
TS = f"{L.ROOT}/svi_science/results/s5_setup/tidal_stress_svi_hourly.parquet"


def main():
    S = pd.read_parquet(f"{OUT}/streaks_all_q8c_d2.parquet")
    M = pd.read_parquet(f"{OUT}/members_data_q8c_d2.parquet")
    LD = {}
    for c in LOAD:
        for v in ["raw", "once"]:
            for s in ["data", "c1"]:
                LD[(c, f"{v}|{s}")] = LOAD[c](v, s)
    _, EB = qd.bostock(); LD[("bostock", "all")] = (None, EB)
    # absolute times
    S["w0"] = pd.NaT; S["onset"] = pd.NaT
    for (c, var), (_, E) in LD.items():
        Ei = E.set_index("ep"); m = (S.cat == c) & (S["var"] == var)
        S.loc[m, "w0"] = S.loc[m, "ep"].map(Ei.w0).values; S.loc[m, "onset"] = S.loc[m, "ep"].map(Ei.onset).values
    S["w0"] = pd.to_datetime(S.w0); S["onset"] = pd.to_datetime(S.onset)
    S["start"] = S.w0 + pd.to_timedelta(S.t_first, unit="h"); S["end"] = S.w0 + pd.to_timedelta(S.t_last, unit="h")
    S["dur_h"] = S.t_last - S.t_first; S["days_after_onset"] = (S.start - S.onset).dt.total_seconds() / 86400
    S.to_parquet(f"{OUT}/streaks_all_abs_q8c_d2.parquet", index=False)

    def sel(c, var, mode="data", cls="RTR"):
        q = S[(S.cat == c) & (S["var"] == var) & (S["mode"] == mode)]
        return q[q.cls == cls] if cls else q

    # ---------------- 1. summary + pass rule
    rows = []
    for c in ["ac_svi", "ac_nw", "ac_sw", "ac_back"]:
        for v in ["raw", "once"]:
            var, varc = f"{v}|data", f"{v}|c1"
            E = LD[(c, var)][1]; Ec = LD[(c, varc)][1]
            R = sel(c, var); Rc = sel(c, varc)
            per = R.groupby("ep").size().reindex(E.ep, fill_value=0)
            perc = Rc.groupby("ep").size().reindex(E.ep, fill_value=0)
            diff = (per - perc).values
            pw = sel(c, var, "pw"); pwc = sel(c, varc, "pw")
            fwd = len(sel(c, var, cls="forward"))
            row = dict(set=c, variant=v, n_eps=len(E), cand_median=float(E.n.median()), scr_cand_median=float(Ec.n.median()),
                       rtr=len(R), per_ep=len(R) / len(E), forward=fwd, surge_only=len(sel(c, var, cls="surge_only")),
                       scr_rtr=len(Rc), scr_per_ep=len(Rc) / len(E), scr_forward=len(sel(c, varc, cls="forward")),
                       ratio_data_scr=len(R) / max(len(Rc), 0.5),
                       wilcoxon_p=float(stats.wilcoxon(diff, alternative="greater").pvalue) if np.any(diff != 0) else 1.0,
                       c1_rtr=len(sel(c, var, "c1")) / 5, c2_rtr=len(sel(c, var, "c2")) / 5, c2_fwd=len(sel(c, var, "c2", "forward")) / 5,
                       pw_rtr_mean=len(pw) / 5 if len(pw) else np.nan, pw_scr_rtr_mean=len(pwc) / 5 if len(pwc) else np.nan,
                       p_back_fwd=float(stats.binomtest(len(R), len(R) + fwd, 0.5).pvalue) if len(R) + fwd else np.nan,
                       days_after_onset=R.days_after_onset.median(), speed=R.v_ev.abs().median(), span=R.span_ev.median(), dur=R.dur_h.median())
            row["PASS"] = bool(row["ratio_data_scr"] >= 2 and row["wilcoxon_p"] < 0.05)
            rows.append(row)
    A = pd.DataFrame(rows); A.to_csv(f"{OUT}/q8c_summary.csv", index=False); print(A.T.to_string())
    RES["summary"] = A.to_dict(orient="records")
    pe = []
    for c in ["ac_svi", "ac_nw", "ac_sw", "ac_back"]:
        for v in ["raw", "once"]:
            E = LD[(c, f"{v}|data")][1]
            R = sel(c, f"{v}|data"); Rc = sel(c, f"{v}|c1")
            for r in E.itertuples():
                q = R[R.ep == r.ep]
                pe.append(dict(set=c, variant=v, ep=r.ep, n_cand=r.n, n_rtr=len(q), n_rtr_scrambled=int((Rc.ep == r.ep).sum()),
                               starts_PT=";".join(pd.Series(q.start.sort_values()).dt.tz_localize("UTC").dt.tz_convert("America/Los_Angeles").dt.strftime("%Y-%m-%d %H:%M").tolist())))
    pd.DataFrame(pe).to_csv(f"{OUT}/q8c_per_episode.csv", index=False)

    # ---------------- 2. thinned Bostock at the new densities
    full = S0 = pd.read_parquet(f"{OUT}/streaks_all_abs_d2.parquet")
    full = full[(full.cat == "bostock") & (full["mode"] == "data") & (full.cls == "RTR")]
    th = {}
    dens = pd.read_csv(f"{OUT}/q8c_thin_densities.csv", index_col=0).iloc[:, 0].to_dict()
    for lab, n in dens.items():
        q = sel("bostock", "all", f"thin@{lab}")
        prec = np.mean([len(full[(full.ep == r.ep) & (full.start <= r.end) & (full.end >= r.start)]) > 0 for r in q.itertuples()]) if len(q) else np.nan
        rec = []
        for rep in range(10):
            qr = q[q.rep == rep]
            rec.append(np.mean([len(qr[(qr.ep == f.ep) & (qr.start <= f.end) & (qr.end >= f.start)]) > 0 for f in full.itertuples()]))
        th[lab] = dict(n_per_ep=int(n), rtr_per_ep=len(q) / 100, precision=float(prec), recall=float(np.mean(rec)))
    th["full_density"] = dict(n_per_ep=19748, rtr_per_ep=len(full) / 10)
    RES["thinned_bostock"] = th; print("thin", th)

    # ---------------- 3. calibration 2010-13 (backward candidates)
    BR = full[full.ep.isin(["2010-08", "2011-08", "2012-08", "2013-09"])]
    Ro = pd.read_csv(f"{OUT}/royer2015_rtr_list.csv", parse_dates=["start_utc"])
    Ro = Ro[Ro.year.isin([2010, 2011, 2012])].rename(columns={"start_utc": "start"})
    Ro["ep"] = Ro.year.map({2010: "2010-08", 2011: "2011-08", 2012: "2012-08"})
    def mcount(our, ref, sh=0.0):
        n = 0
        for ep in ref.ep.unique():
            a = our[our.ep == ep].start.values + np.timedelta64(int(sh * 3600), "s")
            n += len(match(a, ref[ref.ep == ep].start.values, 6.0))
        return n
    cal = {}
    for v in ["raw", "once"]:
        for s in ["data", "c1"]:
            q = sel("ac_back", f"{v}|{s}")
            m1 = mcount(q, BR); ch1 = np.mean([mcount(q, BR, h) for h in (-72, -48, -24, 24, 48, 72)])
            m2 = mcount(q, Ro); ch2 = np.mean([mcount(q, Ro, h) for h in (-72, -48, -24, 24, 48, 72)])
            cal[f"{v}|{s}"] = dict(n_rtr=len(q), per_ep=len(q) / 4, bostock_rtr=len(BR), matched_bostock=m1, recall_bostock=m1 / len(BR),
                                   precision=m1 / max(len(q), 1), chance_recall_bostock=ch1 / len(BR),
                                   royer_n=len(Ro), matched_royer=m2, recall_royer=m2 / len(Ro), chance_recall_royer=ch2 / len(Ro))
    RES["calibration_2010_13"] = cal; print("cal", json.dumps(cal, indent=0))

    # ---------------- 4. questions per set
    def tide_tab(c, var, D):
        ev, E = LD[(c, var)]; Ei = E.set_index("ep"); res = []
        for r in D.itertuples():
            rows_ = M[M.uid == r.uid].row.values
            w = ev.loc[rows_].fam.value_counts()
            e = Ei.loc[r.ep]; path = TB if e.w0 < pd.Timestamp("2014-01-15") else TS
            T = pd.read_parquet(path, columns=["fid", "time", "tau_updip_Pa"],
                                filters=[("fid", "in", [int(f) for f in w.index]), ("time", ">=", e.w0), ("time", "<", e.w1 + pd.Timedelta(hours=2))])
            Pv = T.pivot(index="time", columns="fid", values="tau_updip_Pa"); Pv = Pv[[f for f in w.index if f in Pv.columns]]
            if Pv.shape[1] == 0:
                continue
            ww = w.reindex(Pv.columns).values.astype(float)
            ser = pd.Series(Pv.values @ ww / ww.sum(), index=Pv.index)
            grid = pd.date_range(e.w0, e.w1, freq="10min")
            tsv = np.interp(grid.astype("int64"), ser.index.astype("int64"), ser.values)
            n = max(int(round((r.end - r.start).total_seconds() / 600)), 1); i0 = int(np.searchsorted(grid, r.start))
            g = ev[ev.ep == r.ep]; days = g.time.dt.floor("D").value_counts(); act = days[days >= 0.01 * len(g)].index
            ok = np.isin(grid.floor("D"), act) & (np.arange(len(grid)) + n < len(grid))
            mins = pd.Series(tsv).rolling(n + 1).min().shift(-n).values; means = pd.Series(tsv).rolling(n + 1).mean().shift(-n).values
            res.append(dict(all_pos=bool((tsv[i0:i0 + n + 1] > 0).all()), mean_pos=bool(tsv[i0:i0 + n + 1].mean() > 0),
                            p_all=float(np.nanmean(mins[ok] > 0)), p_mean=float(np.nanmean(means[ok] > 0))))
        T = pd.DataFrame(res)
        if len(T) == 0:
            return dict(n=0)
        s1 = (RNG.random((20000, len(T))) < T.p_all.values).sum(1); s2 = (RNG.random((20000, len(T))) < T.p_mean.values).sum(1)
        return dict(n=len(T), frac_all=float(T.all_pos.mean()), exp_all=float(T.p_all.mean()), p_all=float((1 + (s1 >= T.all_pos.sum()).sum()) / 20001),
                    frac_mean=float(T.mean_pos.mean()), exp_mean=float(T.p_mean.mean()), p_mean=float((1 + (s2 >= T.mean_pos.sum()).sum()) / 20001))

    Q = {}
    for c in ["ac_svi", "ac_nw", "ac_sw", "ac_back"]:
        for v in ["raw", "once"]:
            var = f"{v}|data"; ev, E = LD[(c, var)]; R = sel(c, var).reset_index(drop=True)
            out = dict(n=len(R))
            if len(R) < 3:
                Q[f"{c}|{v}"] = out; continue
            # tides (SVI catalogues only)
            if c in ("ac_svi", "ac_back"):
                out["tides"] = {cls: tide_tab(c, var, sel(c, var, cls=cls)) for cls in ["RTR", "forward", "surge_only"]}
                out["tides"]["scrambled_RTR"] = tide_tab(c, f"{v}|c1", sel(c, f"{v}|c1")) if len(sel(c, f"{v}|c1")) else dict(n=0)
            # arrays per episode
            arr = {}
            for r in E.itertuples():
                g = ev[ev.ep == r.ep]; t = ((g.time - r.w0).dt.total_seconds() / 3600).values
                arr[r.ep] = dict(t=t, x=g.x.values, rows=g.index.values, T=(r.w1 - r.w0).total_seconds() / 3600, on=BZ.onset_map(t, g.x.values))
            # up-dip
            isr = ev.index.isin(M[M.uid.isin(R.uid)].row.unique())
            n_all = ev.groupby("fam").size(); ok = n_all >= 50; dd = ev.groupby("fam").dd.mean()
            n_r = ev[isr].groupby("fam").size().reindex(n_all.index, fill_value=0)
            if ok.sum() >= 8:
                rho = stats.spearmanr((n_r / n_all)[ok], dd[ok])[0]; rs = []
                for _ in range(200):
                    ra = []
                    for r in R.itertuples():
                        a = arr[r.ep]; dur = r.t_last - r.t_first; xr = (min(r.x_start, r.x_end) - 3, max(r.x_start, r.x_end) + 3)
                        for _t in range(20):
                            t0 = RNG.uniform(0, a["T"] - dur); lo = np.searchsorted(a["t"], t0); hi = np.searchsorted(a["t"], t0 + dur)
                            xx = a["x"][lo:hi]; onv = a["on"].reindex(np.floor(xx / BZ.ONSET_BIN).astype(int)).values
                            mk = (onv < t0 - BZ.LAG_H) & (xx >= xr[0]) & (xx <= xr[1])
                            if mk.sum() >= 5:
                                ra.append(a["rows"][lo:hi][mk]); break
                    rr = np.concatenate(ra) if ra else np.array([], int)
                    nr = ev.loc[rr].groupby("fam").size().reindex(n_all.index, fill_value=0)
                    rs.append(stats.spearmanr((nr / n_all)[ok], dd[ok])[0])
                rs = np.array(rs)
                out["updip"] = dict(n_fam=int(ok.sum()), rho=float(rho), ctrl=float(np.nanmean(rs)), ctrl95=[float(np.nanpercentile(rs, 2.5)), float(np.nanpercentile(rs, 97.5))],
                                    p_more_updip=float((1 + (rs <= rho).sum()) / (len(rs) + 1)),
                                    rtr_minus_all_km=float(ev[isr].dd.mean() - ev[ev.ep.isin(R.ep.unique())].dd.mean()))
            # repeat sites
            if R.ep.nunique() >= 3:
                pools = []
                for r in R.itertuples():
                    a = arr[r.ep]; m = np.abs(a["t"] - r.t_first) <= 12
                    pools.append(a["x"][m] if m.any() else a["x"])
                epv = R.ep.values
                nn = lambda xx: np.mean([np.min(np.abs(xx[epv != epv[i]] - xx[i])) for i in range(len(xx))])
                edges = np.arange(-80, 81, 10)
                def hc(xx):
                    H = np.array([np.histogram(xx[epv == e], edges)[0] for e in np.unique(epv)], float); H = H[H.sum(1) >= 2]
                    if len(H) < 3:
                        return np.nan
                    C = np.corrcoef(H); iu = np.triu_indices(len(H), 1); return np.nanmean(C[iu])
                sims = [np.array([RNG.choice(p) for p in pools]) for _ in range(1000)]
                o1, o2 = nn(R.x_start.values), hc(R.x_start.values)
                s1 = np.array([nn(z) for z in sims]); s2 = np.array([hc(z) for z in sims])
                out["repeat"] = dict(nn_obs=float(o1), nn_ctrl=float(s1.mean()), p_nn=float((1 + (s1 <= o1).sum()) / 1001),
                                     hc_obs=float(o2), hc_ctrl=float(np.nanmean(s2)), p_hc=float((1 + (s2 >= o2).sum()) / 1001) if np.isfinite(o2) else np.nan)
            # front
            pct, ds, de = [], [], []
            for r in R.itertuples():
                if not np.isfinite(r.front):
                    continue
                a = arr[r.ep]; m = np.abs(a["t"] - r.t_first) <= 24
                onv = a["on"].reindex(np.floor(a["x"][m] / BZ.ONSET_BIN).astype(int)).values
                ref = np.abs(a["x"][m][onv < a["t"][m] - BZ.LAG_H] - r.front)
                if len(ref) < 20:
                    continue
                d0 = abs(r.x_start - r.front); pct.append((ref < d0).mean()); ds.append(d0); de.append(abs(r.x_end - r.front))
            if len(pct) >= 3:
                out["front"] = dict(n=len(pct), start_km=float(np.median(ds)), end_km=float(np.median(de)), mean_pct=float(np.mean(pct)),
                                    p=float(stats.wilcoxon(np.array(pct) - 0.5, alternative="less").pvalue))
            Q[f"{c}|{v}"] = out
            print(c, v, json.dumps(out, default=str)[:1500], flush=True)
    RES["questions"] = Q
    # regions + eras
    reg = {}
    for v in ["raw", "once"]:
        for a_, b_ in [("ac_svi", "ac_nw"), ("ac_svi", "ac_sw")]:
            Ra, Rb = sel(a_, f"{v}|data"), sel(b_, f"{v}|data"); o = {}
            for col, f in [("speed", lambda d: d.v_ev.abs()), ("span", lambda d: d.span_ev), ("dur", lambda d: d.dur_h), ("timing", lambda d: d.days_after_onset)]:
                if len(Ra) >= 3 and len(Rb) >= 3:
                    o[col] = dict(a=float(f(Ra).median()), b=float(f(Rb).median()), p=float(stats.mannwhitneyu(f(Ra), f(Rb)).pvalue))
            reg[f"{v}:{a_}_vs_{b_}"] = o
    RES["regions"] = reg
    era = {}
    for v, lab in [("raw", "svi_raw"), ("once", "svi_once")]:
        n = len(sel("ac_svi", f"{v}|data")); k = len(LD[("ac_svi", f"{v}|data")][1]); exp = th[lab]["rtr_per_ep"] * k
        nb = len(sel("ac_back", f"{v}|data"))
        era[v] = dict(svi_2014_26=n, eps=k, per_ep=n / k, thinned_bostock_expect=exp, p_poisson_two_sided=float(min(1, 2 * min(stats.poisson.cdf(n, exp), stats.poisson.sf(n - 1, exp)))),
                      back_2010_13_per_ep=nb / 4, back_cand_median=float(LD[("ac_back", f"{v}|data")][1].n.median()),
                      svi_cand_median=float(LD[("ac_svi", f"{v}|data")][1].n.median()))
    RES["eras"] = era; print("eras", era, "regions", json.dumps(reg)[:2000])
    json.dump(RES, open(f"{L.OUT}/results_q8c.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
