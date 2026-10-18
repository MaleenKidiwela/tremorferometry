#!/usr/bin/env python
"""s14_posthoc -- S14 POST HOC (after the C4 result was read, 2:35 pm PT): what is the shared (episode-common) part of the
CC residual? Per region: episode mean of the REG residual over all series (2-5 s), compared with the pooled dv/v anomaly,
the c2 noise-window CC residual, activity (cells, ln n) and class; cross-site correlation of site-episode means.
    python s14_posthoc.py
"""
import os, sys, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import s14_common as C


def site_of(region, ns):
    for g in C.COLOC[region]:
        if ns in g:
            return "+".join(x.split(".")[1] for x in g)
    return ns.split(".")[1]


def run(region):
    os.environ["S11_REGION"] = region
    import importlib, s14_tests as T
    importlib.reload(T)
    E = C.load_episode(region)
    des = T.Design(E, "z_25", T.COVS, "series")
    Q = des.D.assign(res=des.yr)
    Q["dv_a"] = Q.dvv25 - Q.groupby("series").dvv25.transform("median")
    dc = T.Design(E.assign(lnn_c2=np.log(E.n_ctrl2.clip(lower=1))), "cz_25", ["lnn_c2", "lnd"], "series")
    Q = Q.join(pd.Series(dc.yr, index=dc.D.index, name="c2res"))
    Q["site"] = [site_of(region, ns) for ns in Q.netsta]
    P = Q.groupby("episode").agg(date=("date", "median"), cls=("cls", "first"), res=("res", "mean"), dv=("dv_a", "median"),
                                 c2=("c2res", "mean"), ncell=("res", "size"), lnn=("lnn", "mean"), cc=("cc3c_25", "median")).reset_index()
    P = P[P.ncell >= 10].sort_values("date")
    P.to_csv(f"{C.TB}/posthoc_episode_common_{region}.csv", index=False)
    out = dict(n_episodes=len(P), sd_res=float(P.res.std()))
    for c in ("dv", "c2", "lnn", "ncell"):
        ok = P[c].notna() & P.res.notna()
        out[f"r_res_{c}"] = float(np.corrcoef(P.res[ok], P[c][ok])[0, 1])
    out["res_major_minus_burst"] = float(P[P.cls == "major"].res.mean() - P[P.cls == "burst"].res.mean())
    # share of the cell residual variance carried by the episode-common term
    Q = Q.join(P.set_index("episode").res.rename("epres"), on="episode")
    out["frac_var_episode_common"] = float(1 - (Q.res - Q.epres).var() / Q.res.var())
    # cross-site correlation of site-episode means (sites with >= 5 cells in an episode)
    S = Q.groupby(["site", "episode"]).filter(lambda q: len(q) >= 5).groupby(["site", "episode"]).res.mean().unstack(0)
    top = S.notna().sum().sort_values(ascending=False).index[:8]
    R = {}
    for i, a in enumerate(top):
        for b in top[i + 1:]:
            ok = S[a].notna() & S[b].notna()
            if ok.sum() >= 10:
                R[f"{a}~{b}"] = (round(float(np.corrcoef(S[a][ok], S[b][ok])[0, 1]), 3), int(ok.sum()))
    out["cross_site_r"] = R
    out["cross_site_r_median"] = float(np.median([v[0] for v in R.values()])) if R else np.nan
    # the few largest episode-common drops
    out["lowest_episodes"] = P.nsmallest(6, "res")[["episode", "cls", "res", "dv", "c2", "ncell", "lnn"]].round(3).to_dict("records")
    return P, out


def main():
    res = {}; Ps = {}
    for r in ("SVI", "WA_NW"):
        Ps[r], res[r] = run(r)
    json.dump(res, open(f"{C.TB}/posthoc_episode_common.json", "w"), indent=1, default=float)
    print(json.dumps(res, indent=1, default=float))
    fig, axs = plt.subplots(4, 1, figsize=(14, 11), sharex=True)
    for i, (r, P) in enumerate(Ps.items()):
        U = C.units(r)
        a, b = axs[2 * i], axs[2 * i + 1]
        for ax in (a, b):
            for _, u in U[U.cls == "major"].iterrows():
                ax.axvspan(u.start, u.end + pd.Timedelta(days=1), color="#f4b183", alpha=0.35, lw=0)
        a.plot(P.date, P.dv, "-", color="0.7", lw=0.6)
        sc = a.scatter(P.date, P.dv, c=P.cc, cmap="viridis", vmin=0.2, vmax=0.8, s=8 + 0.1 * P.ncell, lw=0.3, edgecolor="0.3")
        a.set_ylabel("dv/v (%)"); a.set_title(f"{r}: pooled per-episode dv/v (median of series anomalies)", fontsize=9)
        b.plot(P.date, P.res, "-", color="0.3", lw=0.7, label="CC residual after stack size + noise (mean over series)")
        b.scatter(P.date, P.res, c=P.cc, cmap="viridis", vmin=0.2, vmax=0.8, s=8 + 0.1 * P.ncell, lw=0.3, edgecolor="0.3")
        b.plot(P.date, P.c2, "-", color="#c00000", lw=0.6, alpha=0.7, label="c2 noise-window CC residual")
        b.axhline(0, color="k", lw=0.4); b.set_ylabel("Δ Fisher z"); b.legend(fontsize=7, loc="lower left")
        b.set_title(f"{r}: episode-common CC residual (POST HOC); r with dv/v {res[r]['r_res_dv']:+.2f}, with c2 {res[r]['r_res_c2']:+.2f}", fontsize=9)
    fig.colorbar(sc, ax=axs, shrink=0.4, label="median stretch CC")
    fig.savefig(f"{C.FIG}/figC1_posthoc_episode_common.png", dpi=120, bbox_inches="tight"); plt.close(fig)
    open(f"{C.OUT}/POSTHOC_DONE", "w").write(C.pt_now() + "\n")


if __name__ == "__main__":
    main()
