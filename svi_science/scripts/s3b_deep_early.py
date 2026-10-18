"""S3b: do the deepest SVI families fire ahead of the (tremor-defined) slip front? Rules: results/s3b/PLAN_S3b.md."""
import numpy as np, pandas as pd, json
R = "/home/jovyan/tremorferometry/svi_science/results"
E = pd.read_parquet(f"{R}/s3/events_svi.parquet"); T = pd.read_parquet(f"{R}/s3/tremor_svi.parquet")
K = pd.read_csv(f"{R}/s3/episodes_kinematics.csv"); K = K[(K.region == "SVI")]
F = pd.read_csv(f"{R}/s0/families.csv")
q1, q2 = F.downdip_km.quantile([1/3, 2/3])
F["dclass"] = np.where(F.downdip_km >= q2, "deep", np.where(F.downdip_km <= q1, "shallow", "mid"))
E = E[(E.series == "data") & (E.P_same >= 0.75) & E.counted_once].copy()
E["fid"] = E.unit.astype(int); E = E.merge(F[["fid", "dclass", "downdip_km"]], on="fid")
eps = K[K.migrates == True].episode.tolist()
def arrival_map(tr):
    tr = tr.assign(b=np.floor(tr.as_km / 5) * 5)
    g = tr.groupby("b").time.apply(lambda s: s.sort_values().iloc[int(0.1 * (len(s) - 1))] if len(s) >= 20 else pd.NaT)
    return g.dropna()
maps = {e: arrival_map(T[T.ep == e]) for e in eps}
def lags(ev, amap):
    b = np.floor(ev.as_km / 5) * 5
    arr = b.map(amap)
    return (ev.time - arr).dt.total_seconds() / 86400
rows = []; allE = []
for e in eps:
    ev = E[E.ep == e].copy(); ev["lag"] = lags(ev, maps[e]); ev = ev.dropna(subset=["lag"]); ev["early"] = ev.lag < -1
    allE.append(ev)
    fr = ev.groupby("dclass").early.agg(["mean", "size"])
    rows.append(dict(episode=e, n=len(ev), early_deep=fr.loc["deep", "mean"] if "deep" in fr.index else np.nan,
                     n_deep=int(fr.loc["deep", "size"]) if "deep" in fr.index else 0,
                     early_shallow=fr.loc["shallow", "mean"] if "shallow" in fr.index else np.nan,
                     n_shallow=int(fr.loc["shallow", "size"]) if "shallow" in fr.index else 0,
                     early_mid=fr.loc["mid", "mean"] if "mid" in fr.index else np.nan,
                     median_lag_deep=ev[ev.dclass == "deep"].lag.median(), median_lag_shallow=ev[ev.dclass == "shallow"].lag.median()))
P = pd.DataFrame(rows); P["delta"] = P.early_deep - P.early_shallow
A = pd.concat(allE)
def delta(df):
    m = df.groupby("dclass").early.mean(); return m.get("deep", np.nan) - m.get("shallow", np.nan)
obs = delta(A)
rng = np.random.default_rng(1)
fam = A[["fid"]].drop_duplicates().merge(F[["fid", "dclass", "along_strike_km"]], on="fid")
fam["sbin"] = np.floor(fam.along_strike_km / 10)
null = []
for _ in range(5000):
    perm = fam.copy(); perm["dclass"] = perm.groupby("sbin").dclass.transform(lambda s: rng.permutation(s.values))
    null.append(delta(A.drop(columns="dclass").merge(perm[["fid", "dclass"]], on="fid")))
null = np.array(null); p = float((np.sum(null >= obs) + 1) / (len(null) + 1))
# mismatched fronts: arrival map of the next episode in the list
mm = []
for i, e in enumerate(eps):
    ev = E[E.ep == e].copy(); ev["lag"] = lags(ev, maps[eps[(i + 1) % len(eps)]] - (maps[eps[(i + 1) % len(eps)]].min() - maps[e].min()))
    ev = ev.dropna(subset=["lag"]); ev["early"] = ev.lag < -1; mm.append(ev)
obs_mm = delta(pd.concat(mm))
flag = ~P.episode.isin(["2017-03", "2025-10"])
A2 = A[A.ep.isin(P.episode[flag])]
out = dict(episodes=eps, n_events=int(len(A)), pooled_delta=obs, perm_p=p, null_q99=float(np.quantile(null, 0.99)),
           sign_pos_episodes=int((P.delta > 0).sum()), n_episodes=int(P.delta.notna().sum()),
           pooled_delta_without_flagged=delta(A2), mismatched_front_delta=obs_mm,
           early_frac_deep=float(A[A.dclass == "deep"].early.mean()), early_frac_shallow=float(A[A.dclass == "shallow"].early.mean()),
           median_lag_deep=float(A[A.dclass == "deep"].lag.median()), median_lag_shallow=float(A[A.dclass == "shallow"].lag.median()),
           downdip_thirds_km=[float(q1), float(q2)])
P.to_csv(f"{R}/s3b/per_episode.csv", index=False); json.dump(out, open(f"{R}/s3b/result.json", "w"), indent=1, default=float)
print(json.dumps(out, indent=1, default=float)); print(P.round(3).to_string())
