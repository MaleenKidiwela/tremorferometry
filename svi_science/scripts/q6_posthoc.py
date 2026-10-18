"""Q6 post-hoc (added 2026-10-08 after the main results): is the arrest effect only the component-linking rule?
Two components closer than 2 days / 15 km merge, so a bin that another component used within ~2 days can never be
'continued into'. Check: walk model restricted to rows with T_any >= 14 days (and >= 60 d), same episode-order shuffle."""
import importlib.util, json
import numpy as np, pandas as pd
ROOT = "/home/jovyan/tremorferometry"; OUT = f"{ROOT}/svi_science/results/q6_front_arrest"
spec = importlib.util.spec_from_file_location("t", f"{ROOT}/svi_science/scripts/q6_tests.py"); t = importlib.util.module_from_spec(spec); spec.loader.exec_module(t)
E = pd.read_csv(f"{OUT}/episodes_pnsn.csv"); be = pd.read_parquet(f"{OUT}/bin_episode.parquet")
W = pd.read_csv(f"{OUT}/walk_rows_pnsn.csv")
H = t.History(be, E); rng = np.random.default_rng(1)
real = np.zeros(E.comp.max() + 1); real[E.comp.values] = E.t_first.values
def stat(slot):
    Ta, Tb, T7 = H.T_all(slot, W.bin.values.astype(float), W.comp.values)
    out = {}
    for lim in (14, 60):
        m = Ta >= lim
        out[f"Tany_ge{lim}"] = t.beta_T(Ta, W.y.values, W.dist.values, W.lsize.values, W.prop.values, m)
        out[f"n_ge{lim}"] = int(m.sum())
    lt = np.log10(np.where(Ta >= 14, Ta, np.nan)); q = np.nanquantile(lt, [0.25, 0.5, 0.75]); k = np.digitize(lt, q); ok = np.isfinite(lt)
    out["stopfrac_q4_ge14"] = [float((W.y.values[ok & (k == j)] == 0).mean()) for j in range(4)]
    out["q_edges_days"] = (10 ** q).tolist()
    return out
r = stat(real); nl = [stat(t.perm_slots(E, rng)) for _ in range(1000)]
res = {"real": r}
for k in ("Tany_ge14", "Tany_ge60"):
    v = np.array([x[k] for x in nl]); res[k] = dict(beta=r[k], p_pos=t.pval(r[k], v, "+"), null_mean=float(np.nanmean(v)),
                                                 null_q05=float(np.nanquantile(v, .05)), null_q95=float(np.nanquantile(v, .95)))
sf = np.array([x["stopfrac_q4_ge14"] for x in nl]); res["stopfrac_null_q05"] = np.quantile(sf, .05, 0).tolist(); res["stopfrac_null_q95"] = np.quantile(sf, .95, 0).tolist()
json.dump(res, open(f"{OUT}/posthoc_ge14.json", "w"), indent=1); print(json.dumps(res, indent=1))
