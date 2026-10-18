#!/usr/bin/env python
"""s11_ball2_locus -- S11 ballistic round 2, §4 source-locus nuisance (round2/PLAN.md + Amendment 3).
Pilot windows (families 10, 87, 17), 2003-02..2006-12, 11 POLARIS + 7 core stations. Per (family, station, epoch,
episode) S lag (and P3) against the leave-episode-out reference (round-1/2 machinery), demeaned per series. Per
(family, episode): lag_s = O + g_s . dx - T_S,s c (WLS). Wander test: chi2 of horizontal dx about zero + geometry-scramble
control (g permuted among stations, 1,000 draws).

    python s11_ball2_locus.py measure --workers 4
    python s11_ball2_locus.py fit
"""
from __future__ import annotations
import os, sys, json, argparse, time, io, urllib.request
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_multi_common as M
from s10_measure import pick_S
from s11_ball_measure import boot, pick_P
import s11_ball2_measure as B2

OUT = f"{S.ROOT}/svi_science/results/s11_depth/ballistic/round2"
TB, FIG, LOGD = f"{OUT}/tables", f"{OUT}/figures", f"{OUT}/logs"
POL = ["C8.KHVB", "C8.MGCB", "C8.SHDB", "C8.SHVB", "C8.TWBB", "CN.LZB", "PO.TWKB", "PO.TSJB", "PO.SILB", "PO.SSIB", "PO.TWGB"]
CORE = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011", "PB.B926", "CN.VGZ", "CN.SNB"]
STATIONS = POL + CORE
SITE = {"PB.B009": "PGC", "PB.B010": "PGC", "PB.B011": "PGC", "CN.PGC": "PGC"}
FAMS = (10, 87, 17)
END = "2007-01-01"
NMIN_LO, NMIN = 10, 20
SEED = 20261014
G_W, G_POS = {}, {}


def run_family(fam, RF, sp_map, A_map):
    rng = np.random.default_rng(SEED + fam)
    out, geo = [], []
    for ns in STATIONS:
        pos = G_POS.get(ns)
        if pos is None:
            continue
        r = RF[(RF.netsta == ns) & RF.row_id.isin(pos.index)].copy()
        r["ix"] = pos.reindex(r.row_id).to_numpy(); W = G_W[ns]
        ok = np.isfinite(W[r.ix.to_numpy(), :, 0]).all(1) if len(r) else np.zeros(0, bool)
        r = r[ok]; r["EP"] = [A_map.get((fam, d)) for d in r.date]
        for ep, g in r.groupby("epoch"):
            g = g.reset_index(drop=True)
            if len(g) < NMIN_LO + 6:
                continue
            Wr = W[g.ix.to_numpy()].astype(np.float64); w = g.P.to_numpy()
            Gs = np.einsum("n,ncs->cs", w, Wr); wt = w.sum(); G = Gs / wt
            tS = pick_S(G, 0.0); sp = float(sp_map.get((fam, ns), np.nan))
            tP = pick_P(G[0], tS, sp if np.isfinite(sp) else 6.0)
            gm = B2.geom2(G, tS, tP, sp)
            geo.append(dict(fam=fam, netsta=ns, epoch=ep, n_rows=len(g), tS=tS, tP=tP, sp_cat=sp, SmP_pick=tS - tP,
                            inc_deg=gm["inc_deg"]))
            for lab, ge in g[g.EP.notna()].groupby("EP"):
                if len(ge) < NMIN_LO or len(ge) > len(g) - 5:
                    continue
                m = (g.EP == lab).to_numpy()
                ref = (Gs - np.einsum("n,ncs->cs", w[m], Wr[m])) / (wt - w[m].sum())
                X = np.concatenate([np.einsum("n,ncs->cs", w[m], Wr[m])[None] / w[m].sum(), boot(Wr[m], w[m], rng)[0]])
                mm = B2.meas(X, ref, gm, ("S", "P3"))
                rr = dict(fam=fam, netsta=ns, epoch=ep, episode=lab, n=int(m.sum()), date=str(pd.to_datetime(ge.date).min().date()))
                for k in ("S", "P3"):
                    l, c, rl = mm[k]
                    rr[f"lag{k}"], rr[f"cc{k}"], rr[f"rail{k}"] = l[0] * 1e3, c[0], bool(rl[0])
                    rr[f"sig{k}"] = np.nanstd(l[1:], ddof=1) * 1e3
                out.append(rr)
    return out, geo


def _run(a):
    try:
        return a[0], run_family(*a), None
    except Exception:
        import traceback
        return a[0], ([], []), traceback.format_exc()[-1000:]


def measure(workers):
    os.makedirs(TB, exist_ok=True); os.makedirs(LOGD, exist_ok=True)
    log = S.Log(f"{LOGD}/locus.log")
    X = pd.read_parquet(f"{S.TB}/rows.parquet", columns=["row_id", "fam", "P", "series"])
    X = X[(X.series == "data") & (X.P >= 0.5)]
    RS = pd.read_parquet(f"{S.TB}/row_station.parquet", columns=["row_id", "netsta", "epoch", "date"])
    RS = RS[RS.netsta.isin(STATIONS) & (RS.date < END)].merge(X, on="row_id")
    Ad = pd.read_csv(f"{M.TB}/active_days.csv")
    A_map = {(int(f), d): e for f, d, e in zip(Ad.fam, Ad.date, Ad.episode)}
    O = pd.read_csv(f"{S.TB}/offsets_sp.csv"); sp_map = {(int(f), n): s for f, n, s in zip(O.fam, O.netsta, O.sp)}
    for ns in STATIONS:
        p = f"{S.WIN}/{ns}.npz"
        if not os.path.exists(p):
            log(f"{ns}: no windows"); continue
        z = np.load(p); rid = z["row_id"]; keep = np.isin(rid, RS.row_id.to_numpy())
        G_W[ns] = z["W"][keep]; G_POS[ns] = pd.Series(np.arange(keep.sum()), index=rid[keep])
        log(f"loaded {ns}: {keep.sum()} rows")
    import multiprocessing as mp
    res, geo = [], []
    with mp.get_context("fork").Pool(workers) as p:
        for f, (o, g), err in p.imap_unordered(_run, [(f, RS[RS.fam == f], sp_map, A_map) for f in FAMS]):
            if err:
                log(f"ERROR fam {f}: {err}")
            res += o; geo += g; log(f"fam {f}: {len(o)} station-episode values")
    pd.DataFrame(res).to_csv(f"{TB}/locus_lags.csv", index=False)
    pd.DataFrame(geo).to_csv(f"{TB}/locus_geom.csv", index=False)
    open(f"{OUT}/LOCUS_MEASURE_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


# ------------------------------------------------------------------------------------------------ geometry + fit
def station_coords():
    p = f"{TB}/locus_station_coords.csv"
    if os.path.exists(p):
        return pd.read_csv(p).set_index("netsta")
    rows = []
    for net in ("C8", "PO", "CN", "PB"):
        stas = ",".join(s.split(".")[1] for s in STATIONS if s.startswith(net + "."))
        u = f"https://service.earthscope.org/fdsnws/station/1/query?net={net}&sta={stas}&level=station&format=text"
        txt = urllib.request.urlopen(u, timeout=60).read().decode()
        D = pd.read_csv(io.StringIO(txt), sep="|")
        for r in D.itertuples():
            rows.append(dict(netsta=f"{r[1]}.{r[2]}", lat=r[3], lon=r[4], elev_m=r[5]))
    C = pd.DataFrame(rows).drop_duplicates("netsta")
    dep = pd.read_csv(f"{S.ROOT}/svi_science/results/s11_depth/tables/sensor_depths.csv")
    dep["netsta"] = dep.Network + "." + dep.Station
    C = C.merge(dep[["netsta", "sensor_depth_m"]], on="netsta", how="left").fillna({"sensor_depth_m": 0.0})
    C.to_csv(p, index=False)
    return C.set_index("netsta")


def vel_profile():
    V = pd.read_csv(f"{S.ROOT}/data/velocity_models/usgs_cascadia_v1.7/profiles/south_VI.csv")
    return V


def straight_TS(src, sta, V, nseg=400):
    """straight-ray S time through the 1-D profile (s)"""
    t = np.linspace(0, 1, nseg + 1); mid = 0.5 * (t[1:] + t[:-1])
    p = src[None] + (sta - src)[None] * mid[:, None]; L = np.linalg.norm(sta - src) / nseg
    z = np.clip(p[:, 2], 0, None)
    i = np.searchsorted(V.ztop_km.to_numpy(), z, side="right") - 1
    return float((L / V.vs.to_numpy()[i]).sum())


def fit():
    rng = np.random.default_rng(SEED)
    L = pd.read_csv(f"{TB}/locus_lags.csv"); Gm = pd.read_csv(f"{TB}/locus_geom.csv")
    C = station_coords(); V = vel_profile()
    fam = pd.read_csv(f"{S.ROOT}/data/catalogue/families_v1_3.csv", low_memory=False)
    fam = fam[fam.region == "S_Vancouver_Is"].copy(); fam["fam"] = pd.to_numeric(fam.family_id, errors="coerce")
    lat0, lon0 = 48.5, -123.6
    pj = lambda la, lo: np.array([(lo - lon0) * 111.19 * np.cos(np.radians(lat0)), (la - lat0) * 111.19])
    # clock: quartet instrument terms (round-2 EP table) where available
    kp = f"{TB}/clock_terms.csv"
    Kt = pd.read_csv(kp) if os.path.exists(kp) else pd.DataFrame(columns=["bin", "lab", "node", "K"])
    Kmap = {(r.lab, r.node): r.K for r in Kt[Kt.bin == "EP"].itertuples()}
    L = L[~L.railS].copy()
    L["kq"] = [Kmap.get((e, f"{ns}|{ep}"), 0.0) for e, ns, ep in zip(L.episode, L.netsta, L.epoch)]
    L["y"] = L.lagS - L.kq
    rows, res_all = [], []
    for nmin, tagn in ((NMIN, "n20"), (NMIN_LO, "n10")):
        d = L[L.n >= nmin].copy()
        d["K"] = d.groupby(["fam", "netsta", "epoch"]).y.transform("size"); d = d[d.K >= 3]
        d["r"] = d.y - d.groupby(["fam", "netsta", "epoch"]).y.transform("mean")
        d["sig"] = np.maximum(d.sigS, 0.5) * np.sqrt(d.K / (d.K - 1))
        geo = {}
        for f in FAMS:
            fr = fam[fam.fam == f].iloc[0]; src = np.r_[pj(fr.lat, fr.lon), fr.depth_km]
            i = np.searchsorted(V.ztop_km.to_numpy(), src[2], side="right") - 1; vs_src = V.vs.to_numpy()[i]
            for ns in STATIONS:
                if ns not in C.index:
                    continue
                c = C.loc[ns]; sta = np.r_[pj(c.lat, c.lon), c.sensor_depth_m / 1000.0]
                u = (sta - src) / np.linalg.norm(sta - src)
                geo[(f, ns)] = (-u / vs_src, straight_TS(src, sta, V), np.linalg.norm(sta - src))
        fits = []
        for (f, e), g in d.groupby(["fam", "episode"]):
            g = g[[(f, ns) in geo for ns in g.netsta]]
            nsite = len(set(SITE.get(ns, ns) for ns in g.netsta))
            if len(g) < 7 or nsite < 5:
                continue
            gv = np.array([geo[(f, ns)][0] for ns in g.netsta]); TS = np.array([geo[(f, ns)][1] for ns in g.netsta])
            y = g.r.to_numpy() / 1e3; s = g.sig.to_numpy() / 1e3
            fits.append(dict(fam=f, episode=e, date=g.date.min(), n_sta=len(g), n_site=nsite, gv=gv, TS=TS, y=y, s=s,
                             cc=float(g.ccS.median()), sta=list(g.netsta)))

        def solve(fi, gv=None, full=True):
            gv = fi["gv"] if gv is None else gv
            A = np.c_[np.ones(len(fi["y"])), gv, -fi["TS"]] if full else np.c_[np.ones(len(fi["y"])), -fi["TS"]]
            w = 1 / fi["s"]; Aw = A * w[:, None]; yw = fi["y"] * w
            N = Aw.T @ Aw; Ci = np.linalg.pinv(N); x = Ci @ Aw.T @ yw
            chi2 = float(((yw - Aw @ x) ** 2).sum()); dof = len(yw) - A.shape[1]
            return x, Ci, chi2, dof

        X2 = 0.0; out = []
        for fi in fits:
            x, Ci, chi2, dof = solve(fi); x0, C0, chi20, dof0 = solve(fi, full=False)
            h = x[1:3]; Ch = Ci[1:3, 1:3]; q = float(h @ np.linalg.pinv(Ch) @ h); X2 += q
            out.append(dict(nmin=tagn, fam=fi["fam"], episode=fi["episode"], date=fi["date"], n_sta=fi["n_sta"],
                            n_site=fi["n_site"], cc=fi["cc"], O_ms=x[0] * 1e3, dE_m=x[1] * 1e3, dN_m=x[2] * 1e3,
                            dZ_m=x[3] * 1e3, sE_m=np.sqrt(Ci[1, 1]) * 1e3, sN_m=np.sqrt(Ci[2, 2]) * 1e3,
                            sZ_m=np.sqrt(Ci[3, 3]) * 1e3, c_pct=x[4] * 100, s_c_pct=np.sqrt(Ci[4, 4]) * 100,
                            c_noloc_pct=x0[1] * 100, s_c_noloc_pct=np.sqrt(C0[1, 1]) * 100, chi2=chi2, dof=dof,
                            chi2_noloc=chi20, q_h=q, corr_Z_O=Ci[3, 0] / np.sqrt(Ci[3, 3] * Ci[0, 0]),
                            corr_c_O=Ci[4, 0] / np.sqrt(Ci[4, 4] * Ci[0, 0])))
        if not fits:
            rows.append(dict(nmin=tagn, n_fits=0)); continue
        from scipy.stats import chi2 as CH
        p_chi = float(CH.sf(X2, 2 * len(fits)))
        # geometry-scramble control: g rows permuted among the stations of each fit
        Xs = np.empty(1000)
        for b in range(1000):
            t = 0.0
            for fi in fits:
                perm = rng.permutation(len(fi["y"]))
                x, Ci, _, _ = solve(fi, gv=fi["gv"][perm]); h = x[1:3]
                t += float(h @ np.linalg.pinv(Ci[1:3, 1:3]) @ h)
            Xs[b] = t
        p_scr = float((1 + (Xs >= X2).sum()) / 1001)
        F = pd.DataFrame(out); res_all.append(F)
        rms_h = float(np.sqrt(np.mean(F.dE_m ** 2 + F.dN_m ** 2)))
        noise_h = float(np.sqrt(np.mean(F.sE_m ** 2 + F.sN_m ** 2)))
        motion = float(np.sqrt(max(rms_h ** 2 - noise_h ** 2, 0)))
        rows.append(dict(nmin=tagn, n_fits=len(fits), X2_h=X2, dof_h=2 * len(fits), p_chi2=p_chi, p_scramble=p_scr,
                         scr_X2_median=float(np.median(Xs)), detected=bool(p_chi < 0.01 and p_scr < 0.01),
                         rms_h_m=rms_h, formal_h_m=noise_h, motion_m=motion,
                         bound2s_m=float(np.sqrt(motion ** 2 + 2 * noise_h ** 2 / np.sqrt(len(fits)))),
                         chi2_dof_fit=float(F.chi2.sum() / F.dof.sum()),
                         c_mean_pct=float(np.average(F.c_pct, weights=1 / F.s_c_pct ** 2)),
                         c_noloc_mean_pct=float(np.average(F.c_noloc_pct, weights=1 / F.s_c_noloc_pct ** 2))))
    Rr = pd.DataFrame(rows); Rr.to_csv(f"{TB}/locus_summary.csv", index=False)
    Fa = pd.concat(res_all, ignore_index=True) if res_all else pd.DataFrame()
    Fa.to_csv(f"{TB}/locus_fits.csv", index=False)
    with pd.option_context("display.width", 250):
        print(Rr.round(3).to_string()); print(Fa.round(2).to_string())
    figure(Fa)


def figure(F):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import Normalize
    if not len(F):
        return
    F = F[F.nmin == "n20"] if (F.nmin == "n20").sum() >= 3 else F[F.nmin == "n10"]
    F = F.copy(); F["date"] = pd.to_datetime(F.date)
    norm = Normalize(0.5, 1.0); cmap = plt.get_cmap("viridis")
    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
    for f, g in F.groupby("fam"):
        g = g.sort_values("date")
        for ax, col in zip(axs[:2], ("dE_m", "dN_m")):
            ax.plot(g.date, g[col], color="0.6", lw=0.8)
            sc = ax.scatter(g.date, g[col], c=g.cc, cmap=cmap, norm=norm, s=30, lw=0)
            ax.text(g.date.iloc[-1], g[col].iloc[-1], f" fam {f}", fontsize=7)
        axs[2].plot(g.dE_m, g.dN_m, color="0.7", lw=0.6)
        axs[2].scatter(g.dE_m, g.dN_m, c=g.cc, cmap=cmap, norm=norm, s=30, lw=0)
    for ax, t in zip(axs, ("east shift (m)", "north shift (m)")):
        ax.set_ylabel(t); ax.axhline(0, color="0.5", lw=0.5); ax.grid(color="0.92")
    axs[2].set_xlabel("east shift (m)"); axs[2].set_ylabel("north shift (m)"); axs[2].axhline(0, color="0.5", lw=0.5)
    axs[2].axvline(0, color="0.5", lw=0.5); axs[2].set_aspect("equal", "datalim")
    fig.colorbar(sc, ax=axs, shrink=0.8, label="median S cc")
    fig.suptitle("Source-locus shift per family-episode (2003–06, POLARIS + core), relative to each family's mean", fontsize=10)
    fig.savefig(f"{FIG}/fig_locus.png", dpi=110, bbox_inches="tight"); plt.close(fig)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what"); ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    measure(a.workers) if a.what == "measure" else fit()
