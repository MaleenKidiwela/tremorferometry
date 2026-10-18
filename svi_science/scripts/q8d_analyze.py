"""Q8d analysis: WA P-weighted RTR test vs scrambled control; backward, up-dip, repeat, front, vs SVI."""
import json, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_data as qd, q8_rtr_lib as L, q8_rtr_backzone as BZ
OUT = f"{L.OUT}/tables"; RNG = np.random.default_rng(4)
S = pd.read_parquet(f"{OUT}/streaks_all_q8d_d2.parquet"); M = pd.read_parquet(f"{OUT}/members_data_q8d_d2.parquet")
RES = {}
svi = json.load(open(f"{L.OUT}/results_q8c_pw.json"))["SVI"]
ev_s, E_s = qd.allcand_svi("raw", "data"); svi_psum = float(ev_s.groupby("ep").P.sum().sum())
for cat, reg in [("ac_nw", "NW"), ("ac_sw", "SW")]:
    ev, E = qd.allcand_wa(reg, "raw", "data"); evc, Ec = qd.allcand_wa(reg, "raw", "c1")
    Ei = E.set_index("ep")
    d = S[(S.cat == cat) & (S["var"] == "raw|data")].copy()
    d["start_h"] = d.t_first; d["w0"] = d.ep.map(Ei.w0); d["onset"] = d.ep.map(Ei.onset)
    d["start"] = d.w0 + pd.to_timedelta(d.t_first, unit="h"); d["days_after_onset"] = (d.start - d.onset).dt.total_seconds() / 86400
    d["dur_h"] = d.t_last - d.t_first
    sc = S[(S.cat == cat) & (S["var"] == "raw|c1") & (S["mode"] == "pw")]
    pw = d[d["mode"] == "pw"]
    R = pw[pw.cls == "RTR"]
    per = R.groupby("ep").size().reindex(E.ep, fill_value=0) / 5
    perc = sc[sc.cls == "RTR"].groupby("ep").size().reindex(E.ep, fill_value=0) / 5
    diff = (per - perc).values
    rt, fw = len(R) / 5, (pw.cls == "forward").sum() / 5
    psum = float(ev.groupby("ep").P.sum().sum()); psum_c = float(evc.groupby("ep").P.sum().sum())
    r = dict(n_eps=len(E), P_sum_per_ep=psum / len(E), scr_P_sum_per_ep=psum_c / len(E),
             rtr_per_run=rt, per_ep=rt / len(E), scr_rtr_per_run=float(perc.sum()), ratio=rt / max(perc.sum(), 0.1),
             wilcoxon_p=float(stats.wilcoxon(diff, alternative="greater").pvalue) if np.any(diff != 0) else 1.0,
             forward_per_run=float(fw), scr_forward_per_run=(sc.cls == "forward").sum() / 5,
             p_back_fwd=float(stats.binomtest(int(len(R)), int(len(R) + (pw.cls == "forward").sum()), 0.5).pvalue) if rt + fw else None,
             time_shuffled_rtr=(d[(d["mode"] == "pw_c1") & (d.cls == "RTR")].shape[0]) / 5,
             pos_perm_rtr=(d[(d["mode"] == "pw_c2") & (d.cls == "RTR")].shape[0]) / 5, pos_perm_fwd=(d[(d["mode"] == "pw_c2") & (d.cls == "forward")].shape[0]) / 5,
             speed=float(R.v_ev.abs().median()) if len(R) else None, span=float(R.span_ev.median()) if len(R) else None,
             dur=float(R.dur_h.median()) if len(R) else None, days_after_onset=float(R.days_after_onset.median()) if len(R) else None)
    r["PASS"] = bool(r["ratio"] >= 2 and r["wilcoxon_p"] < 0.05)
    # vs SVI per episode and per 1000 P-sum
    n_svi, k_svi = svi["rtr_per_run"], 10
    r["vs_svi_per_ep_p"] = float(stats.binomtest(int(round(len(R))), int(round(len(R) + n_svi * 5)), len(E) / (len(E) + k_svi)).pvalue)
    r["per_1000_Psum"] = rt / psum * 1000; r["svi_per_1000_Psum"] = n_svi / svi_psum * 1000
    r["vs_svi_per_Psum_p"] = float(stats.binomtest(int(round(len(R))), int(round(len(R) + n_svi * 5)), psum / (psum + svi_psum)).pvalue)
    # up-dip, repeat, front (pooled draws for up-dip; draw 0 for repeat/front)
    if len(R) >= 5:
        arr = {}
        for e in E.itertuples():
            g = ev[ev.ep == e.ep]; t = ((g.time - e.w0).dt.total_seconds() / 3600).values
            arr[e.ep] = dict(t=t, x=g.x.values, rows=g.index.values, T=(e.w1 - e.w0).total_seconds() / 3600, on=BZ.onset_map(t, g.x.values))
        isr = ev.index.isin(M[M.uid.isin(R.uid)].row.unique())
        n_all = ev.groupby("fam").size(); ok = n_all >= 50; dd = ev.groupby("fam").dd.mean()
        n_r = ev[isr].groupby("fam").size().reindex(n_all.index, fill_value=0)
        rho = stats.spearmanr((n_r / n_all)[ok], dd[ok])[0]; rs = []
        R0 = R[R.rep == 0]
        for _ in range(200):
            ra = []
            for q in R.itertuples():
                a = arr[q.ep]; dur = q.t_last - q.t_first; xr = (min(q.x_start, q.x_end) - 3, max(q.x_start, q.x_end) + 3)
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
        r["updip"] = dict(n_fam=int(ok.sum()), rho=float(rho), ctrl=float(np.nanmean(rs)), p_more_updip=float((1 + (rs <= rho).sum()) / (len(rs) + 1)),
                          rtr_minus_all_km=float(ev[isr].dd.mean() - ev[ev.ep.isin(R.ep.unique())].dd.mean()))
        pct, ds, de = [], [], []
        for q in R.itertuples():
            if not np.isfinite(q.front):
                continue
            a = arr[q.ep]; m = np.abs(a["t"] - q.t_first) <= 24
            onv = a["on"].reindex(np.floor(a["x"][m] / BZ.ONSET_BIN).astype(int)).values
            ref = np.abs(a["x"][m][onv < a["t"][m] - BZ.LAG_H] - q.front)
            if len(ref) < 20:
                continue
            d0 = abs(q.x_start - q.front); pct.append((ref < d0).mean()); ds.append(d0); de.append(abs(q.x_end - q.front))
        if len(pct) >= 3:
            r["front"] = dict(n=len(pct), start_km=float(np.median(ds)), end_km=float(np.median(de)), mean_pct=float(np.mean(pct)),
                              p=float(stats.wilcoxon(np.array(pct) - 0.5, alternative="less").pvalue))
        if R0.ep.nunique() >= 3:
            epv = R0.ep.values; xs = R0.x_start.values
            pools = []
            for q in R0.itertuples():
                a = arr[q.ep]; m = np.abs(a["t"] - q.t_first) <= 12; pools.append(a["x"][m] if m.any() else a["x"])
            nn = lambda z: np.mean([np.min(np.abs(z[epv != epv[i]] - z[i])) for i in range(len(z))])
            o = nn(xs); sim = np.array([nn(np.array([RNG.choice(p) for p in pools])) for _ in range(1000)])
            r["repeat_draw0"] = dict(n=len(R0), nn_obs=float(o), nn_ctrl=float(sim.mean()), p=float((1 + (sim <= o).sum()) / 1001))
    per_ep_tab = pd.DataFrame(dict(ep=E.ep.values, P_sum=E.ep.map(ev.groupby("ep").P.sum()).values, rtr_per_run=per.values, scr_rtr_per_run=perc.values))
    per_ep_tab.to_csv(f"{OUT}/q8d_per_episode_{reg}.csv", index=False)
    RES[reg] = r
    print(reg, json.dumps(r, indent=1, default=str), flush=True)
RES["SVI_reference"] = dict(rtr_per_run=svi["rtr_per_run"], per_ep=svi["per_ep"], per_1000_Psum=svi["rtr_per_run"] / svi_psum * 1000,
                            scr_rtr_per_run=svi["scr_rtr_per_run"], forward_per_run=svi["forward_per_run"])
json.dump(RES, open(f"{L.OUT}/results_q8d.json", "w"), indent=1, default=str)
print(RES["SVI_reference"])
