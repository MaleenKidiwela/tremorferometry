"""Q8c add-on (pre-registered variant c): P-weighted pass rule per episode + tides/up-dip for SVI and NW."""
import json, sys
import numpy as np, pandas as pd
from scipy import stats
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_data as qd, q8_rtr_lib as L
OUT = f"{L.OUT}/tables"; RNG = np.random.default_rng(9)
S = pd.read_parquet(f"{OUT}/streaks_all_abs_q8c_d2.parquet"); M = pd.read_parquet(f"{OUT}/members_data_q8c_d2.parquet")
TB = f"{L.ROOT}/svi_science/results/s5/tidal_stress_bostock_hourly.parquet"; TS = f"{L.ROOT}/svi_science/results/s5_setup/tidal_stress_svi_hourly.parquet"
RES = {}
for c, reg in [("ac_svi", "SVI"), ("ac_nw", "NW"), ("ac_sw", "SW")]:
    E = (qd.allcand_svi("raw", "data") if c == "ac_svi" else qd.allcand_wa(reg, "raw", "data"))[1]
    d = S[(S.cat == c) & (S["var"] == "raw|data") & (S["mode"] == "pw")]
    sc = S[(S.cat == c) & (S["var"] == "raw|c1") & (S["mode"] == "pw")]
    pd_ = d[d.cls == "RTR"].groupby("ep").size().reindex(E.ep, fill_value=0) / 5
    ps_ = sc[sc.cls == "RTR"].groupby("ep").size().reindex(E.ep, fill_value=0) / 5
    fw = (d.cls == "forward").sum() / 5; rt = (d.cls == "RTR").sum() / 5
    r = dict(rtr_per_run=float(rt), per_ep=float(rt / len(E)), forward_per_run=float(fw), scr_rtr_per_run=float(ps_.sum()) if len(sc) else None,
             ratio=float(rt / max(ps_.sum(), 0.1)) if len(sc) else None,
             wilcoxon_p=float(stats.wilcoxon((pd_ - ps_).values, alternative="greater").pvalue) if len(sc) and np.any(pd_ != ps_) else None,
             speed=float(d[d.cls == "RTR"].v_ev.abs().median()) if rt else None, span=float(d[d.cls == "RTR"].span_ev.median()) if rt else None,
             dur=float(d[d.cls == "RTR"].dur_h.median()) if rt else None, days_after_onset=float(d[d.cls == "RTR"].days_after_onset.median()) if rt else None,
             p_back_fwd=float(stats.binomtest(int(round(rt * 5)), int(round((rt + fw) * 5)), 0.5).pvalue) if rt + fw else None)
    if r["ratio"] is not None:
        r["PASS"] = bool(r["ratio"] >= 2 and r["wilcoxon_p"] < 0.05)
    RES[reg] = r
print(json.dumps(RES, indent=1))
# tides for SVI pw RTRs, per rep
ev, E = qd.allcand_svi("raw", "data"); evc, Ec = qd.allcand_svi("raw", "c1")
def tides(D, ev, E):
    Ei = E.set_index("ep"); res = []
    for r in D.itertuples():
        w = ev.loc[M[M.uid == r.uid].row.values].fam.value_counts(); e = Ei.loc[r.ep]
        T = pd.read_parquet(TS, columns=["fid", "time", "tau_updip_Pa"], filters=[("fid", "in", [int(f) for f in w.index]), ("time", ">=", e.w0), ("time", "<", e.w1 + pd.Timedelta(hours=2))])
        Pv = T.pivot(index="time", columns="fid", values="tau_updip_Pa"); Pv = Pv[[f for f in w.index if f in Pv.columns]]
        ww = w.reindex(Pv.columns).values.astype(float); ser = pd.Series(Pv.values @ ww / ww.sum(), index=Pv.index)
        grid = pd.date_range(e.w0, e.w1, freq="10min"); tsv = np.interp(grid.astype("int64"), ser.index.astype("int64"), ser.values)
        n = max(int(round((r.end - r.start).total_seconds() / 600)), 1); i0 = int(np.searchsorted(grid, r.start))
        g = ev[ev.ep == r.ep]; days = g.time.dt.floor("D").value_counts(); act = days[days >= 0.01 * len(g)].index
        ok = np.isin(grid.floor("D"), act) & (np.arange(len(grid)) + n < len(grid))
        mins = pd.Series(tsv).rolling(n + 1).min().shift(-n).values
        res.append(dict(all_pos=bool((tsv[i0:i0 + n + 1] > 0).all()), p_all=float(np.nanmean(mins[ok] > 0))))
    T = pd.DataFrame(res)
    return T
out = {}
for lab, c_var, evx, Ex in [("data_RTR", "raw|data", ev, E), ("scrambled_RTR", "raw|c1", evc, Ec)]:
    for cls in (["RTR", "forward", "surge_only"] if lab == "data_RTR" else ["RTR"]):
        D = S[(S.cat == "ac_svi") & (S["var"] == c_var) & (S["mode"] == "pw") & (S.cls == cls)]
        T = tides(D, evx, Ex)
        if len(T) == 0:
            out[f"{lab}_{cls}"] = dict(n=0); continue
        # p from Poisson-binomial, with n divided by 5 (draws share events): conservative effective n
        k = len(T) / 5
        sim = (RNG.random((20000, len(T))) < T.p_all.values).mean(1)
        out[f"{lab}_{cls}"] = dict(n_pooled_5_draws=len(T), frac_all=float(T.all_pos.mean()), exp_all=float(T.p_all.mean()),
                                   p_pooled=float((1 + (sim >= T.all_pos.mean()).sum()) / 20001),
                                   p_effective_n_per_draw=float(stats.binomtest(int(round(T.all_pos.mean() * k)), max(int(round(k)), 1), T.p_all.mean(), alternative="greater").pvalue))
RES["tides_svi_pw"] = out
print(json.dumps(out, indent=1))
json.dump(RES, open(f"{L.OUT}/results_q8c_pw.json", "w"), indent=1)
