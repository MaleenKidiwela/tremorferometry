#!/usr/bin/env python
"""s11_gate_a -- gate (a), kernel-free co-located test at the PGC site (PGC, B009, B010, B011).
Annual cycle per instrument (and per component when the S11 measurement exists) on series-demeaned episode values.
Part 1: existing S10-multi Z 2-4 s / 3-5 s episode stacks. Part 2 (if tables/episode_raw.parquet exists): S11 raw path
2-5 s, components 3C / Z / H. Output tables/gate_a.csv. eps = -dv/v; dv/v in %."""
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s11_common as C

U = C.episodes()
dep = pd.read_csv(f"{C.TB}/sensor_depths.csv"); dep["netsta"] = dep.Network + "." + dep.Station
dep = dep.set_index("netsta").sensor_depth_m
rows = []


def run(E, src, col, sigcol, extra_mask=None):
    E = E.copy(); E["dvv"] = -100 * E[col]; E["sig"] = 100 * E[sigcol]
    E = E[np.isfinite(E.dvv) & np.isfinite(E.sig) & (E.sig > 0)]
    if extra_mask is not None:
        E = E[extra_mask(E)]
    E["dvv"] = E.dvv - E.groupby(["fam", "netsta", "epoch"]).dvv.transform(lambda x: np.average(x, weights=1 / E.loc[x.index, "sig"] ** 2))
    E["mid"] = E.episode.map(U.mid)
    res = {}
    for ns in C.COLOC + ["PB.B926", "CN.VGZ", "CN.SNB"]:
        q = E[E.netsta == ns]
        # station-episode value: inverse-variance mean over families (cells), then annual fit over episodes
        g = q.groupby("episode").apply(lambda x: pd.Series(dict(v=np.average(x.dvv, weights=1 / x.sig ** 2),
                                                                 s=1 / np.sqrt(np.sum(1 / x.sig ** 2)), n=len(x))), include_groups=False)
        g = g[g.n >= 3]
        f = C.annual_fit(U.mid.reindex(g.index), g.v.to_numpy(), g.s.to_numpy())
        f.update(src=src, comp=col, netsta=ns, sensor_depth_m=float(dep.get(ns, np.nan)), n_cells=len(q), n_episodes=len(g))
        rows.append(f); res[ns] = f
    a = np.array([res[n]["amp"] for n in C.COLOC]); s = np.array([res[n]["amp_sig"] for n in C.COLOC])
    ok = np.isfinite(a) & np.isfinite(s)
    if ok.sum() >= 2:
        w = 1 / s[ok] ** 2; m = np.sum(w * a[ok]) / w.sum(); chi2 = float(np.sum(w * (a[ok] - m) ** 2))
        from scipy import stats
        rows.append(dict(src=src, comp=col, netsta="EQUALITY_4", amp=m, amp_sig=1 / np.sqrt(w.sum()), chi2red=chi2 / (ok.sum() - 1),
                         n=int(ok.sum()), p_equal=float(stats.chi2.sf(chi2, ok.sum() - 1))))


E = pd.read_parquet(f"{C.MTB}/episode_stacks.parquet")
MAJ = lambda q: q.episode.map(U.cls).eq("major")
for col in ("24", "35"):
    run(E, "s10multi_existing_all", f"eps{col}", f"sig{col}")
    run(E, "s10multi_existing_majors", f"eps{col}", f"sig{col}", extra_mask=MAJ)
p = f"{C.TB}/episode_raw.parquet"
if os.path.exists(p):
    E = pd.read_parquet(p)
    for col in ("eps3c_25", "epsZ_25", "epsH_25", "eps3c_2035", "eps3c_3560"):
        sc = col.replace("eps", "sig")
        if col in E and sc in E:
            run(E, "s11_raw_all", col, sc, extra_mask=lambda q: (q.n >= C.N_COL))
            run(E, "s11_raw_majors", col, sc, extra_mask=lambda q: (q.n >= C.N_COL) & MAJ(q))
T = pd.DataFrame(rows)
T.to_csv(f"{C.TB}/gate_a.csv", index=False)
pd.set_option("display.width", 200)
print(T[["src", "comp", "netsta", "sensor_depth_m", "amp", "amp_sig", "peak_doy", "chi2red", "n_episodes"] +
        (["p_equal"] if "p_equal" in T else [])].round(4).to_string())
