"""Q8 extra checks (post-hoc, D2): thinned precision, histogram-correlation repeat test, tides of surge classes."""
import json
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import q8_rtr_data as qd
import q8_rtr_lib as L

OUT = f"{L.OUT}/tables"
RNG = np.random.default_rng(5)
S = pd.read_parquet(f"{OUT}/streaks_all_abs_d2.parquet")
M = pd.read_parquet(f"{OUT}/members_data_d2.parquet")
RES = {}

# 1. thinned precision / recall
full = S[(S.cat == "bostock") & (S["mode"] == "data") & (S.cls == "RTR")]
allsurge = S[(S.cat == "bostock") & (S["mode"] == "data")]
th = S[(S.cat == "bostock") & (S["mode"] == "thin") & (S.cls == "RTR")]
hit_rtr = hit_any = 0
for r in th.itertuples():
    z = full[(full.ep == r.ep) & (full.start <= r.end) & (full.end >= r.start)]
    z2 = allsurge[(allsurge.ep == r.ep) & (allsurge.start <= r.end) & (allsurge.end >= r.start)]
    hit_rtr += len(z) > 0; hit_any += len(z2) > 0
RES["thin"] = dict(n_thin_rtr_total_10reps=len(th), per_ep=len(th) / 100, frac_overlapping_full_rtr=hit_rtr / max(len(th), 1),
                   frac_overlapping_any_full_surge=hit_any / max(len(th), 1),
                   thin_c2_per_ep=float(((S.cat == "bostock") & (S["mode"] == "thin_c2") & (S.cls == "RTR")).sum() / 100),
                   thin_c1_per_ep=float(((S.cat == "bostock") & (S["mode"] == "thin_c1") & (S.cls == "RTR")).sum() / 100))

# 2. repeat places: mean correlation of per-episode start histograms (10-km bins), control = random event of same episode within +-12 h
ev, E = qd.bostock()
R = full.reset_index(drop=True)
edges = np.arange(-60, 61, 10)
def hcorr(xs, eps):
    H = np.array([np.histogram(xs[eps == e], edges)[0] for e in np.unique(eps)], float)
    H = H[H.sum(1) >= 3]
    C = np.corrcoef(H); iu = np.triu_indices(len(H), 1)
    return np.nanmean(C[iu]), len(H)
pools = []
for r in R.itertuples():
    g = ev[(ev.ep == r.ep) & ((ev.time - r.start).abs() <= pd.Timedelta(hours=12))]
    pools.append(g.x.values)
out = {}
for col in ["x_start", "x_end"]:
    obs, nh = hcorr(R[col].values, R.ep.values)
    sim = np.array([hcorr(np.array([RNG.choice(p) for p in pools]), R.ep.values)[0] for _ in range(1000)])
    out[col] = dict(obs=float(obs), ctrl_mean=float(sim.mean()), ctrl_p95=float(np.percentile(sim, 95)),
                    p=float((1 + (sim >= obs).sum()) / 1001), n_eps=nh)
RES["b_hist_corr_bostock"] = out
# along-strike distribution of starts (for the write-up)
RES["bostock_start_quantiles_km"] = R.x_start.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).round(1).tolist()
RES["bostock_end_quantiles_km"] = R.x_end.quantile([0.1, 0.25, 0.5, 0.75, 0.9]).round(1).tolist()

# 3. tides for surge classes (Bostock): all-positive / mean-positive vs expectation
tb = f"{L.ROOT}/svi_science/results/s5/tidal_stress_bostock_hourly.parquet"
def tide_frac(D):
    res = []
    for r in D.itertuples():
        rows = M[M.uid == r.uid].row.values
        w = ev.loc[rows].fam.value_counts()
        e = E.set_index("ep").loc[r.ep]
        T = pd.read_parquet(tb, columns=["fid", "time", "tau_updip_Pa"],
                            filters=[("fid", "in", [int(f) for f in w.index]), ("time", ">=", e.w0), ("time", "<", e.w1 + pd.Timedelta(hours=2))])
        P = T.pivot(index="time", columns="fid", values="tau_updip_Pa")
        P = P[[c for c in w.index if c in P.columns]]
        ww = w.reindex(P.columns).values.astype(float)
        ser = pd.Series(P.values @ ww / ww.sum(), index=P.index)
        grid = pd.date_range(e.w0, e.w1, freq="10min")
        ts = np.interp(grid.astype("int64"), ser.index.astype("int64"), ser.values)
        n = max(int(round((r.end - r.start).total_seconds() / 600)), 1)
        i0 = int(np.searchsorted(grid, r.start))
        g = ev[ev.ep == r.ep]; days = g.time.dt.floor("D").value_counts(); act = days[days >= 0.01 * len(g)].index
        ok = np.isin(grid.floor("D"), act) & (np.arange(len(grid)) + n < len(grid))
        mins = pd.Series(ts).rolling(n + 1).min().shift(-n).values; means = pd.Series(ts).rolling(n + 1).mean().shift(-n).values
        res.append(dict(all_pos=(ts[i0:i0 + n + 1] > 0).all(), mean_pos=ts[i0:i0 + n + 1].mean() > 0,
                        p_all=np.nanmean(mins[ok] > 0), p_mean=np.nanmean(means[ok] > 0)))
    D2 = pd.DataFrame(res)
    return dict(n=len(D2), frac_all=float(D2.all_pos.mean()), exp_all=float(D2.p_all.mean()),
                frac_mean=float(D2.mean_pos.mean()), exp_mean=float(D2.p_mean.mean()))
RES["tides_surge_only_bostock"] = tide_frac(S[(S.cat == "bostock") & (S["mode"] == "data") & (S.cls == "surge_only")])
print(json.dumps(RES, indent=1))
json.dump(RES, open(f"{L.OUT}/results_extra_d2.json", "w"), indent=1)
