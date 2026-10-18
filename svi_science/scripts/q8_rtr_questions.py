"""Q8 RTR questions a-e (PLAN.md §4). Needs q8_rtr_analyze.py outputs (results_step1.json, streaks_all_abs.parquet)."""
import json
import sys
import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_data as qd
import q8_rtr_lib as L
import q8_rtr_backzone as BZ

ROOT = L.ROOT
OUT = f"{L.OUT}/tables"
RES = {}
RNG = np.random.default_rng(81008)
SUF = sys.argv[1] if len(sys.argv) > 1 else "_d2"
NDRAW = 1000


def ep_arrays(ev, E):
    D = {}
    for r in E.itertuples():
        g = ev[ev.ep == r.ep]
        D[r.ep] = dict(t=((g.time - r.w0).dt.total_seconds() / 3600).values, x=g.x.values, rows=g.index.values,
                       T=(r.w1 - r.w0).total_seconds() / 3600, w0=r.w0, onset=r.onset)
    return D


def corridor(a, t0, D, v, x0):
    lo = np.searchsorted(a["t"], t0); hi = np.searchsorted(a["t"], t0 + D)
    xc = x0 + v * (a["t"][lo:hi] - t0)
    m = np.abs(a["x"][lo:hi] - xc) <= L.W
    return a["rows"][lo:hi][m]


def win_events(a, on, ta, tb, xr=None):
    """events in [ta, tb) inside the slipped zone at ta (bins with onset < ta - 12 h), optional x range."""
    lo = np.searchsorted(a["t"], ta); hi = np.searchsorted(a["t"], tb)
    xx = a["x"][lo:hi]
    onv = on.reindex(np.floor(xx / BZ.ONSET_BIN).astype(int)).values
    m = onv < ta - BZ.LAG_H
    if xr is not None:
        m &= (xx >= xr[0]) & (xx <= xr[1])
    return a["rows"][lo:hi][m]


def front_pos(onset, t0):
    b = [k for k, v in onset.items() if t0 - 24 <= v <= t0]
    return (np.median(b) + 0.5) * L.ONSET_BIN if b else np.nan


def line_class_geom(onset, t0, D, v, x0):
    xs = x0 + v * np.linspace(0, D, 41)
    po = L._path_onsets(onset, xs); have = ~np.isnan(po)
    if have.sum() == 0:
        return False
    behind = (po[have] < t0).mean()
    q1 = np.nanmedian(po[:11]) if np.isfinite(po[:11]).any() else np.nan
    q4 = np.nanmedian(po[-11:]) if np.isfinite(po[-11:]).any() else np.nan
    return bool(behind >= L.BEHIND_MIN and np.isfinite(q1 - q4) and (q1 - q4) >= L.BACK_H)


def main():
    st1 = json.load(open(f"{L.OUT}/results_step1{SUF}.json"))
    SST = st1["S_star_frozen"]; OWN = st1["S_star_own"]
    S = pd.read_parquet(f"{OUT}/streaks_all_abs{SUF}.parquet")
    M = pd.read_parquet(f"{OUT}/members_data{SUF}.parquet")
    sets = {"bostock": ("bostock", "all", qd.bostock()), "svi": ("svi", "head", qd.ours_svi("head")),
            "nw": ("nw", "head", qd.ours_nw("head")), "svi_p050": ("svi", "p050", qd.ours_svi("p050")),
            "nw_p050": ("nw", "p050", qd.ours_nw("p050")), "svi_raw075": ("svi", "raw075", qd.ours_svi("raw075")),
            "nw_raw075": ("nw", "raw075", qd.ours_nw("raw075"))}
    ARR = {k: ep_arrays(*v[2]) for k, v in sets.items()}
    ONB = {k: {ep: BZ.onset_map(a["t"], a["x"]) for ep, a in ARR[k].items()} for k in sets}

    def rtrs(k, thr="frozen"):
        c, v, _ = sets[k]
        s = SST if thr == "frozen" else OWN.get(f"{c}|{v}|data", np.inf)
        return S[(S.cat == c) & (S["var"] == v) & (S["mode"] == "data") & (S.cls == "RTR") & (S.S >= s)].copy()

    # ------------------------------------------------------------ (a) counts / timing / kinematics
    rows = []
    for k in sets:
        for thr in ("frozen", "own"):
            R = rtrs(k, thr); E = sets[k][2][1]
            nev_rtr = M[M.uid.isin(R.uid)].row.nunique()
            rows.append(dict(set=k, threshold=thr, S=SST if thr == "frozen" else OWN.get(f"{sets[k][0]}|{sets[k][1]}|data"),
                             n_rtr=len(R), n_eps=len(E), eps_with_rtr=R.ep.nunique(), per_ep_mean=len(R) / len(E),
                             frac_events_in_rtr=nev_rtr / E.n.sum(),
                             days_after_onset_med=R.days_after_onset.median(), speed_med=R.v_ev.abs().median(),
                             span_med=R.span_ev.median(), dur_med=R.dur_h.median(),
                             n_ev_med=R.N_ev.median(), n_forward=int(((S.cat == sets[k][0]) & (S["var"] == sets[k][1]) & (S["mode"] == "data")
                                                                     & (S.cls == "forward") & (S.S >= (SST if thr == "frozen" else OWN.get(f"{sets[k][0]}|{sets[k][1]}|data", np.inf)))).sum())))
    A = pd.DataFrame(rows); A.to_csv(f"{OUT}/a_summary{SUF}.csv", index=False); print(A.to_string())
    per_ep = []
    for k in ["bostock", "svi", "nw"]:
        R = rtrs(k); E = sets[k][2][1]
        for r in E.itertuples():
            q = R[R.ep == r.ep]
            per_ep.append(dict(set=k, ep=r.ep, n_events=r.n, n_rtr=len(q), days_after_onset=";".join(f"{d:.1f}" for d in sorted(q.days_after_onset))))
    pd.DataFrame(per_ep).to_csv(f"{OUT}/a_per_episode{SUF}.csv", index=False)
    cols = ["uid", "cat", "var", "ep", "start", "end", "days_after_onset", "dur_h", "v_ev", "span_ev", "x_start", "x_end",
            "front", "N_ev", "n_fam", "S", "rho", "bg_run", "n_run"]
    pd.concat([rtrs(k)[cols] for k in ["bostock", "svi", "nw"]]).to_csv(f"{OUT}/rtr_list_frozen{SUF}.csv", index=False)
    comp = {}
    for a_, b_ in [("bostock", "svi"), ("bostock", "nw"), ("svi", "nw")]:
        Ra, Rb = rtrs(a_), rtrs(b_)
        for col, f in [("speed", lambda d: d.v_ev.abs()), ("span", lambda d: d.span_ev), ("dur", lambda d: d.dur_h),
                       ("days_after_onset", lambda d: d.days_after_onset)]:
            if len(Ra) >= 3 and len(Rb) >= 3:
                comp[f"{a_}_vs_{b_}_{col}"] = dict(med_a=float(f(Ra).median()), med_b=float(f(Rb).median()),
                                                   p=float(stats.mannwhitneyu(f(Ra), f(Rb)).pvalue))
    RES["a_compare"] = comp
    # thinned Bostock: expected count at our density + recall of full-density RTRs
    Bf = rtrs("bostock")
    th = S[(S.cat == "bostock") & (S["mode"] == "thin") & (S.cls == "RTR")]
    own_thin = OWN.get("bostock|all|thin", np.inf)
    thin = {}
    for lab, s in [("frozen", SST), ("own_thin", own_thin)]:
        q = th[th.S >= s]
        cnt = q.groupby("rep").size().reindex(range(10), fill_value=0)
        rec = []
        for rep in range(10):
            qr = q[q.rep == rep]; hit = 0
            for f in Bf.itertuples():
                z = qr[(qr.ep == f.ep) & (qr.start <= f.end) & (qr.end >= f.start) & (np.sign(qr.v_ev) == np.sign(f.v_ev))]
                hit += len(z) > 0
            rec.append(hit / max(len(Bf), 1))
        thin[lab] = dict(S=s, rtr_per_10eps_mean=float(cnt.mean()), rtr_per_ep=float(cnt.mean() / 10), recall_of_full=float(np.mean(rec)),
                         precision_proxy=None)
    RES["a_thinned_bostock"] = thin
    print("thin", thin)

    # ------------------------------------------------------------ (b) same places
    bres = {}
    for k in ["bostock", "svi", "nw"]:
        R = rtrs(k).reset_index(drop=True)
        if R.ep.nunique() < 2:
            bres[k] = "fewer than 2 episodes with RTRs"; continue
        out = {}
        for which, xcol, tcol in [("start", "x_start", "t_first"), ("end", "x_end", "t_last")]:
            x = R[xcol].values; ep = R.ep.values
            def stat(xx):
                d = []
                for i in range(len(xx)):
                    o = xx[ep != ep[i]]
                    d.append(np.min(np.abs(o - xx[i])))
                return np.mean(d)
            obs = stat(x)
            pools = []
            for r in R.itertuples():
                a = ARR[k][r.ep]; tt = getattr(r, tcol)
                m = np.abs(a["t"] - tt) <= 12
                pools.append(a["x"][m] if m.any() else a["x"])
            sim = np.array([stat(np.array([RNG.choice(p) for p in pools])) for _ in range(NDRAW)])
            out[which] = dict(obs_mean_nn_km=float(obs), ctrl_mean=float(sim.mean()), ctrl_p05=float(np.percentile(sim, 5)),
                              p_closer=float((1 + (sim <= obs).sum()) / (NDRAW + 1)), n=len(R), n_eps=int(R.ep.nunique()))
        bres[k] = out
    RES["b_same_places"] = bres
    print("b", json.dumps(bres, indent=1))

    # ------------------------------------------------------------ (c) start at the main front (D2 front estimate)
    cres = {}
    for k in ["bostock", "svi", "nw"]:
        R = rtrs(k)
        pct, d_obs, d_ref, d_end = [], [], [], []
        for r in R.itertuples():
            if not np.isfinite(r.front):
                continue
            a = ARR[k][r.ep]; on = ONB[k][r.ep]
            lo = np.searchsorted(a["t"], r.t_first - 24); hi = np.searchsorted(a["t"], r.t_first + 24)
            xx = a["x"][lo:hi]; tt = a["t"][lo:hi]
            onv = on.reindex(np.floor(xx / BZ.ONSET_BIN).astype(int)).values
            m = onv < tt - BZ.LAG_H
            ref = np.abs(xx[m] - r.front)
            if len(ref) < 20:
                continue
            ds = abs(r.x_start - r.front)
            pct.append((ref < ds).mean()); d_obs.append(ds); d_ref.append(np.median(ref)); d_end.append(abs(r.x_end - r.front))
        if len(pct) < 3:
            cres[k] = dict(n=len(pct)); continue
        pct = np.array(pct)
        cres[k] = dict(n=len(pct), median_dist_start_to_front_km=float(np.median(d_obs)),
                       median_dist_end_to_front_km=float(np.median(d_end)),
                       ref_median_slipped_zone_km=float(np.median(d_ref)), mean_percentile=float(pct.mean()),
                       p_closer=float(stats.wilcoxon(pct - 0.5, alternative="less").pvalue),
                       share_start_within_5km=float(np.mean(np.array(d_obs) <= 5)),
                       share_start_within_10km=float(np.mean(np.array(d_obs) <= 10)))
    RES["c_front"] = cres
    print("c", json.dumps(cres, indent=1))

    # ------------------------------------------------------------ (d) tides (SVI only)
    dres = {}
    Ro = pd.read_csv(f"{OUT}/royer2015_rtr_list.csv", parse_dates=["start_utc", "end_utc"])
    evB, EB = sets["bostock"][2]
    evS, ES = sets["svi"][2]

    def tide_test(rt, tide_path, ev, E, members_fn):
        """rt: DataFrame with ep, start, end. members_fn(row) -> Series of fid counts."""
        res = []
        for r in rt.itertuples():
            w = members_fn(r)
            if len(w) == 0:
                continue
            e = E.set_index("ep").loc[r.ep]
            T = pd.read_parquet(tide_path, columns=["fid", "time", "tau_updip_Pa"],
                                filters=[("fid", "in", [int(f) for f in w.index]), ("time", ">=", e.w0), ("time", "<", e.w1 + pd.Timedelta(hours=2))])
            P = T.pivot(index="time", columns="fid", values="tau_updip_Pa")
            P = P[[c for c in w.index if c in P.columns]]
            ww = w.reindex(P.columns).values.astype(float)
            ser = pd.Series(P.values @ ww / ww.sum(), index=P.index)
            grid = pd.date_range(e.w0, e.w1, freq="10min")
            ts = np.interp(grid.astype("int64"), ser.index.astype("int64"), ser.values)
            dur = (r.end - r.start).total_seconds() / 3600
            n = max(int(round(dur * 6)), 1)
            i0 = int(np.searchsorted(grid, r.start))
            obs_all = bool((ts[i0:i0 + n + 1] > 0).all()); obs_mean = float(ts[i0:i0 + n + 1].mean())
            # random starts within active days
            g = ev[ev.ep == r.ep]; days = g.time.dt.floor("D").value_counts()
            act = days[days >= 0.01 * len(g)].index
            okstart = np.isin(grid.floor("D"), act) & (np.arange(len(grid)) + n < len(grid))
            mins = pd.Series(ts).rolling(n + 1).min().shift(-n).values
            means = pd.Series(ts).rolling(n + 1).mean().shift(-n).values
            p_all = float(np.nanmean(mins[okstart] > 0)); p_mean = float(np.nanmean(means[okstart] > 0))
            res.append(dict(ep=r.ep, start=r.start, dur_h=dur, all_pos=obs_all, mean_udss=obs_mean, mean_pos=obs_mean > 0,
                            p_all_pos=p_all, p_mean_pos=p_mean))
        return pd.DataFrame(res)

    def summarize(D, lab):
        if len(D) == 0:
            return dict(n=0)
        sim = (RNG.random((20000, len(D))) < D.p_all_pos.values).sum(1)
        sim2 = (RNG.random((20000, len(D))) < D.p_mean_pos.values).sum(1)
        o = int(D.all_pos.sum()); o2 = int(D.mean_pos.sum())
        D.to_csv(f"{OUT}/d_tides_{lab}{SUF}.csv", index=False)
        return dict(n=len(D), frac_all_positive=o / len(D), expected_all_positive=float(D.p_all_pos.mean()),
                    p_all=float((1 + (sim >= o).sum()) / 20001), frac_mean_positive=o2 / len(D),
                    expected_mean_positive=float(D.p_mean_pos.mean()), p_mean=float((1 + (sim2 >= o2).sum()) / 20001),
                    median_mean_udss_Pa=float(D.mean_udss.median()))
    tb = f"{ROOT}/svi_science/results/s5/tidal_stress_bostock_hourly.parquet"
    ts_ = f"{ROOT}/svi_science/results/s5_setup/tidal_stress_svi_hourly.parquet"
    ep_of_year = {2003: "2003-02", 2004: "2004-07", 2005: "2005-09", 2008: "2008-05", 2010: "2010-08", 2011: "2011-08", 2012: "2012-08"}
    Ro["ep"] = Ro.year.map(ep_of_year); Ro = Ro.rename(columns={"start_utc": "start", "end_utc": "end"})
    def royer_members(r):
        g = evB[(evB.time >= r.start) & (evB.time <= r.end)]
        return g.fam.value_counts()
    def mem_fn(ev):
        def f(r):
            rows = M[M.uid == r.uid].row.values
            return ev.loc[rows].fam.value_counts()
        return f
    dres["royer_list_2003_2012"] = summarize(tide_test(Ro, tb, evB, EB, royer_members), "royer")
    dres["ours_on_bostock_2003_2013"] = summarize(tide_test(rtrs("bostock"), tb, evB, EB, mem_fn(evB)), "bostock")
    dres["ours_svi_2014_2026"] = summarize(tide_test(rtrs("svi"), ts_, evS, ES, mem_fn(evS)), "svi")
    rs_own = rtrs("svi", "own")
    if len(rs_own):
        dres["ours_svi_2014_2026_ownS"] = summarize(tide_test(rs_own, ts_, evS, ES, mem_fn(evS)), "svi_own")
    # forward streaks on Bostock as a comparison class
    fw = S[(S.cat == "bostock") & (S["mode"] == "data") & (S.cls == "forward") & (S.S >= SST)]
    dres["bostock_forward_streaks"] = summarize(tide_test(fw, tb, evB, EB, mem_fn(evB)), "bostock_forward")
    RES["d_tides"] = dres
    print("d", json.dumps(dres, indent=1, default=str))

    # ------------------------------------------------------------ (e) sizes and depth
    eres = {}
    # attach S6 dM to our SVI events
    s6 = pd.read_parquet(f"{ROOT}/svi_science/results/s6/sizes.parquet", columns=["fam", "time_utc", "dM"])
    s6["time"] = s6.time_utc.astype("datetime64[ns]").dt.round("10ms")
    for k in ["svi", "svi_p050"]:
        ev = sets[k][2][0]
        ev["_tr"] = ev.time.dt.round("10ms")
        j = ev.reset_index().merge(s6[["fam", "time", "dM"]].rename(columns={"time": "_tr"}), on=["fam", "_tr"], how="left").drop_duplicates("index").set_index("index")
        ev["mag"] = j.dM.reindex(ev.index).values
    for k, lab in [("bostock", "mag"), ("svi", "dM"), ("svi_p050", "dM")]:
        ev, E = sets[k][2]
        R = rtrs(k)
        if len(R) < 3:
            eres[k] = dict(n=len(R)); continue
        mem_rows = set(M[M.uid.isin(R.uid)].row.values)
        isr = ev.index.isin(list(mem_rows))
        med = ev[~isr].groupby(["ep", "fam"]).mag.median()
        anom = ev.mag.values - med.reindex(list(zip(ev.ep, ev.fam))).values
        an = pd.Series(anom, index=ev.index)
        def per_rtr(rows):
            v = an.loc[rows].dropna()
            return v.mean() if len(v) >= 3 else np.nan
        obs = np.array([per_rtr(M[M.uid == u].row.values) for u in R.uid])
        # concurrent events outside the corridor
        conc = []
        for r in R.itertuples():
            a = ARR[k][r.ep]
            lo = np.searchsorted(a["t"], r.t_first); hi = np.searchsorted(a["t"], r.t_last, side="right")
            rows = np.setdiff1d(a["rows"][lo:hi], M[M.uid == r.uid].row.values)
            conc.append(per_rtr(rows) if len(rows) else np.nan)
        conc = np.array(conc)
        # shifted corridors
        sims = []
        for _ in range(NDRAW if k == "bostock" else 500):
            vals = []
            for r in R.itertuples():
                a = ARR[k][r.ep]
                dur = r.t_last - r.t_first; xr = (min(r.x_start, r.x_end) - 3, max(r.x_start, r.x_end) + 3)
                for _t in range(20):
                    t0 = RNG.uniform(0, a["T"] - dur)
                    rows = win_events(a, ONB[k][r.ep], t0, t0 + dur, xr)
                    if len(rows) >= 5:
                        vals.append(per_rtr(rows)); break
            sims.append(np.nanmean(vals))
        sims = np.array(sims)
        o = np.nanmean(obs)
        eres[k] = dict(measure=lab, n=int(np.isfinite(obs).sum()), rtr_minus_family_median=float(o),
                       concurrent_outside_corridor=float(np.nanmean(conc)),
                       rtr_minus_concurrent_paired=float(np.nanmean(obs - conc)),
                       paired_p=float(stats.wilcoxon((obs - conc)[np.isfinite(obs - conc)]).pvalue) if np.isfinite(obs - conc).sum() > 5 else np.nan,
                       shifted_ctrl_mean=float(np.nanmean(sims)), shifted_ctrl_p975=float(np.nanpercentile(sims, 97.5)),
                       p_larger_than_shifted=float((1 + (sims >= o).sum()) / (len(sims) + 1)))
        # depth: family RTR share vs down-dip
        if k == "bostock" or k == "svi":
            n_all = ev.groupby("fam").size(); n_r = ev[isr].groupby("fam").size().reindex(n_all.index, fill_value=0)
            dd = ev.groupby("fam").dd.first()
            ok = n_all >= (50 if k == "bostock" else 10)
            share = (n_r / n_all)[ok]
            rho = stats.spearmanr(share, dd[ok])[0]
            # control: shifted corridors -> shares
            rs = []
            for _ in range(200):
                rows_all = []
                for r in R.itertuples():
                    a = ARR[k][r.ep]
                    dur = r.t_last - r.t_first; xr = (min(r.x_start, r.x_end) - 3, max(r.x_start, r.x_end) + 3)
                    for _t in range(20):
                        t0 = RNG.uniform(0, a["T"] - dur)
                        rows = win_events(a, ONB[k][r.ep], t0, t0 + dur, xr)
                        if len(rows) >= 5:
                            rows_all.append(rows); break
                rows_all = np.concatenate(rows_all) if rows_all else np.array([], int)
                nr = ev.loc[rows_all].groupby("fam").size().reindex(n_all.index, fill_value=0)
                rs.append(stats.spearmanr((nr / n_all)[ok], dd[ok])[0])
            rs = np.array(rs)
            ddr = ev[isr].dd.mean() - ev[ev.ep.isin(R.ep.unique())].dd.mean()
            eres[k].update(depth_spearman_share_vs_downdip=float(rho), depth_ctrl_mean=float(np.nanmean(rs)),
                           depth_ctrl_p2_5=float(np.nanpercentile(rs, 2.5)), depth_ctrl_p97_5=float(np.nanpercentile(rs, 97.5)),
                           depth_p_more_updip=float((1 + (rs <= rho).sum()) / (len(rs) + 1)),
                           rtr_events_mean_downdip_minus_all_km=float(ddr), n_families=int(ok.sum()))
            fam_tab = pd.DataFrame(dict(n_all=n_all, n_rtr=n_r, share=n_r / n_all, downdip_km=dd))
            fam_tab.to_csv(f"{OUT}/e_family_rtr_share_{k}{SUF}.csv")
    RES["e_sizes_depth"] = eres
    print("e", json.dumps(eres, indent=1, default=str))
    json.dump(RES, open(f"{L.OUT}/results_questions{SUF}.json", "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
