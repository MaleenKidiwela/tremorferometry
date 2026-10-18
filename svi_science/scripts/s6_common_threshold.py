"""S6 post-hoc check (labelled as such in RESULTS_S6.md, not pre-registered): tests 1-3 with ONE completeness threshold
per family shared by the two compared classes. The pre-registered 'complete' flag uses each family-month's own threshold
(Mc_fm = Mc_z + noise), so a noisier class keeps only bigger events by construction. Here: T(f) = max over the two
classes of the median Mc_fm; events with dM >= T(f) from P>=0.75 counted-once rows; >= 3 events per class."""
import json
import numpy as np, pandas as pd

OUT = "/home/jovyan/tremorferometry/svi_science/results/s6"
S = pd.read_parquet(f"{OUT}/sizes_annotated.parquet")
S = S[(S.P_same >= 0.75) & S.keep_xfam8 & S.dM.notna() & S.Mc_fm.notna()].copy()
rng = np.random.default_rng(7)
NPERM = 300
EP = pd.read_csv("/home/jovyan/tremorferometry/lfe_extension_wa/results/episode_rates/episodes_svi_final.csv")
MJ = EP[EP.major]; big = dict(zip(MJ.episode, MJ.family_exposure_days > MJ.family_exposure_days.median()))


def cells(D, cls_col, a, b, unit):
    rows = []
    for key, g in D.groupby(unit, observed=True):
        ga, gb = g[g[cls_col] == a], g[g[cls_col] == b]
        if len(ga) < 3 or len(gb) < 3: continue
        T = max(ga.Mc_fm.median(), gb.Mc_fm.median())
        xa, xb = ga[ga.dM >= T].dM, gb[gb.dM >= T].dM
        if len(xa) >= 3 and len(xb) >= 3: rows.append((key, xa.median() - xb.median(), len(xa), len(xb)))
    return pd.DataFrame(rows, columns=["unit", "d", "na", "nb"])


def boot(C):
    fam = C.unit.map(lambda u: u[0] if isinstance(u, tuple) else u)
    g = {f: C.d.values[fam.values == f] for f in fam.unique()}; fs = list(g)
    bs = [np.mean(np.concatenate([g[fs[i]] for i in rng.choice(len(fs), len(fs))])) for _ in range(2000)]
    return list(np.percentile(bs, [2.5, 97.5]))


res = {}
# test 1: behind vs front, unit = family x episode, permutation of labels within cell
D1 = S[S.front_cls.notna()].copy()
C = cells(D1, "front_cls", "behind", "front", ["fam", "episode"])
null = []
for _ in range(NPERM):
    P = D1.copy(); P["front_cls"] = P.groupby(["fam", "episode"], observed=True).front_cls.transform(lambda s: rng.permutation(s.values))
    null.append(cells(P, "front_cls", "behind", "front", ["fam", "episode"]).d.mean())
res["t1_behind_minus_front"] = dict(n_cells=len(C), effect=C.d.mean(), ci=boot(C), p=float((np.sum(np.abs(null) >= abs(C.d.mean())) + 1) / (NPERM + 1)))
# test 2: big vs small, unit = family, permutation of episode labels
D2 = S[S.major == True].copy(); D2["cls"] = D2.episode.astype(str).map(big)
C = cells(D2, "cls", True, False, "fam")
eps = list(big); v = np.array([big[e] for e in eps]); null = []
for _ in range(NPERM):
    m = dict(zip(eps, rng.permutation(v))); P = D2.assign(cls=D2.episode.astype(str).map(m))
    null.append(cells(P, "cls", True, False, "fam").d.mean())
res["t2_big_minus_small"] = dict(n_fam=len(C), effect=C.d.mean(), ci=boot(C), p=float((np.sum(np.abs(null) >= abs(C.d.mean())) + 1) / (NPERM + 1)))
# test 3: outside all episodes vs major ETS (the 76 S7 events are too few for >=3 per family)
D3 = S[(S.major == True) | (S.in_episode == False)].copy(); D3["cls"] = np.where(S.loc[D3.index, "in_episode"] == False, "out", "ets")
C = cells(D3, "cls", "out", "ets", "fam")
null = []
for _ in range(NPERM):
    P = D3.copy(); P["cls"] = P.groupby("fam").cls.transform(lambda s: rng.permutation(s.values))
    null.append(cells(P, "cls", "out", "ets", "fam").d.mean())
res["t3_outside_minus_ETS"] = dict(n_fam=len(C), effect=C.d.mean(), ci=boot(C), p=float((np.sum(np.abs(null) >= abs(C.d.mean())) + 1) / (NPERM + 1)))
C3b = cells(S[(S.major == True) | S.between].assign(cls=lambda d: np.where(d.between, "be", "ets")), "cls", "be", "ets", "fam")
res["t3_S7_between_minus_ETS_n_fam"] = len(C3b)
json.dump(res, open(f"{OUT}/common_threshold_check.json", "w"), indent=1, default=float)
print(json.dumps(res, indent=1, default=float))
