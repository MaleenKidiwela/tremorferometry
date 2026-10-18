"""C check (NOT pre-registered; added after seeing the IET figure): are the very short inter-event times (seconds to a
minute) double counts that survived 'counted once'? Then re-do Q2 with events merged within M seconds inside each patch.
Also: episode median patch duration vs front speed (post-hoc)."""
import json
import numpy as np
import pandas as pd
import sys
sys.path.insert(0, "/home/jovyan/tremorferometry/svi_science/scripts")
import c_decay as C

R = C.R; OUT = C.OUT
res = {}
Es = pd.read_parquet(f"{R}/s3/events_svi.parquet")
Es = Es[(Es.series == "data") & (Es.P_same >= 0.75) & Es.counted_once & Es.ep.isin(C.EPS["SVI"])].copy()
Es["patch"] = "G" + Es.group.astype(str); Es["src"] = Es.unit.astype(str)
En = pd.read_parquet(f"{R}/s3/events_nw.parquet")
En = En[(En.series == "data") & (En.P_same >= 0.75) & En.selected_v2 & En.counted_once & En.ep.isin(C.EPS["NW"])].copy()
En["dcl"] = C.dclass(En.dd_km.values, "NW")
En["patch"] = "B" + (np.floor(En.as_km / 10) * 10).astype(int).astype(str) + "_" + En.dcl
En["src"] = En.source
Ec = pd.read_parquet(f"{R}/s3/events_svi.parquet")
Ec = Ec[(Ec.series == "c1") & (Ec.P_same >= 0.5) & Ec.counted_once & Ec.ep.isin(C.EPS["SVI"])].copy()
Ec["patch"] = "G" + Ec.group.astype(str); Ec["src"] = Ec.unit.astype(str)

def pairs(E):
    out = []
    for (ep, pa), g in E.groupby(["ep", "patch"]):
        g = g.sort_values("time"); t = g.time.values.astype("datetime64[ns]").astype(np.int64) / 1e9
        d = np.diff(t); same = g.src.values[1:] == g.src.values[:-1]
        das = np.abs(np.diff(g.as_km.values))
        out.append(pd.DataFrame(dict(dt=d, same=same, das=das)))
    return pd.concat(out)

for nm, E in [("SVI", Es), ("NW", En), ("SVI_scrambled", Ec)]:
    P = pairs(E); r = {}
    for lim in [5, 10, 30, 60, 300]:
        m = P.dt < lim
        r[f"lt{lim}s"] = dict(n=int(m.sum()), frac_of_gaps=float(m.mean()), frac_same_source=float(P[m].same.mean()) if m.any() else None)
    r["n_gaps"] = int(len(P))
    res[f"short_iet_{nm}"] = r

res["merged"] = {}
for nm, E, reg in [("SVI", Es, "SVI"), ("NW", En, "NW")]:
    X = C.add_lag(C.prep(E[["time", "ep", "as_km", "dd_km", "patch", "P_same"]].assign(region=reg, qualC=False), reg), reg)
    for M in [10, 60]:
        keep = []
        for (ep, pa), g in X.groupby(["ep", "patch"]):
            t = g.time.values.astype(np.int64) / 1e9
            k = np.ones(len(g), bool); last = -1e18
            for i, x in enumerate(t):
                if x - last < M:
                    k[i] = False
                else:
                    last = x
            keep.append(g[k])
        Y = pd.concat(keep)
        C.NSUR = 200
        Cq, r2, sz, _ = C.q2_cells(Y, reg, f"merged{M}s", 30)
        r2["n_events"] = int(len(Y)); r2["n_cells"] = int(len(Cq))
        r2["median_frac_burst"] = float(Cq.frac_burst.median()); r2["median_frac_burst_Sday"] = float(Cq.fb_day_med.median())
        r2["median_frac_burst_hour"] = float(Cq.fb_hour_med.median())
        r2["median_mean_iet_min"] = float(Cq.mean_iet_min.median())
        r2["burst_sizes"] = C.size_fits(sz["obs"]); r2["burst_sizes_Sday"] = C.size_fits(sz["day"])
        res["merged"][f"{nm}_{M}s"] = r2

# post-hoc: episode median patch duration vs |front speed|
CT = pd.read_csv(f"{OUT}/cells_duration.csv")
res["duration_vs_speed"] = {}
for reg in ["SVI", "NW"]:
    m = CT[(CT.region == reg) & (CT.run == "headline")].dropna(subset=["D80_n30"]).groupby("ep").D80_n30.median()
    sp = C.K[C.K.region == reg].set_index("episode").lfe_speed.abs()
    x = pd.concat([m, sp], axis=1, join="inner").dropna()
    rho = C.stats.spearmanr(x.iloc[:, 0], x.iloc[:, 1])[0]
    nul = np.array([C.stats.spearmanr(C.rng.permutation(x.iloc[:, 0].values), x.iloc[:, 1])[0] for _ in range(10000)])
    res["duration_vs_speed"][reg] = dict(rho=float(rho), p=float((np.sum(np.abs(nul) >= abs(rho)) + 1) / 10001), n=int(len(x)))
json.dump(res, open(f"{OUT}/result_check.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
