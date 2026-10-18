#!/usr/bin/env python
"""s11_ball3_common -- S11 ballistic round 3, §3 region-wide (common-mode) change (round3/PLAN.md).
Per region: (1) S - P3 year series per station and pooled region-wide (dv/v units, random-effects error across stations);
(2) catalogue drift: absolute S and P3 lags (detector frame) vs S - P; (3) joint uniform a = dVs/Vs, b = dVp/Vp from
per-series decade contrasts of S - P, dS and dP3 (moveout) with within-series year permutation; (4) per ETS: region-wide
S - P per major (EP bins), chi2 about zero vs the scrambled control. Success rule 3.

    python s11_ball3_common.py [--tag _test] [--regions SVI,WA_NW]
"""
from __future__ import annotations
import os, sys, json, argparse
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_ball3_lib as L
import s11_ball2_longterm as LT

NPERM = 2000
SEED = 20261013
JREF = {"SVI": ("BH", ["PB.B009", "PB.B010", "PB.B011"]), "WA_NW": ("G1BH", ["PB.B005", "PB.B006", "PB.B007"])}
SPLIT_SETS = {"SVI": ["CN.PGC", "PB.B009", "PB.B010", "PB.B011"],
              "WA_NW": ["PB.B005", "PB.B006", "PB.B007", "UW.LRIV", "PB.B013", "PB.B943"]}
JREF2 = {"SVI": ("B011", ["PB.B011"]), "WA_NW": ("G2", ["PB.B013", "PB.B943"])}


class Contr:
    """per-series decade contrasts with within-series permutation (vectorised)"""

    def __init__(self, sid, y, s2, t):
        self.sid = np.asarray(sid); self.y = np.asarray(y, float); self.s2 = np.asarray(s2, float)
        self.late = np.asarray(t) >= L.SPLIT; self.ns = int(self.sid.max()) + 1 if len(self.sid) else 0
        self.nl = np.bincount(self.sid, self.late, self.ns); self.ne = np.bincount(self.sid, ~self.late, self.ns)
        self.ok = (self.nl >= 2) & (self.ne >= 2)
        self.order = np.argsort(self.sid, kind="stable"); self.sid_s = self.sid[self.order]

    def contrast(self, y=None, s2=None):
        y = self.y if y is None else y; s2 = self.s2 if s2 is None else s2
        C = np.bincount(self.sid, y * self.late, self.ns) / np.maximum(self.nl, 1) - \
            np.bincount(self.sid, y * ~self.late, self.ns) / np.maximum(self.ne, 1)
        v = np.bincount(self.sid, s2 * self.late, self.ns) / np.maximum(self.nl, 1) ** 2 + \
            np.bincount(self.sid, s2 * ~self.late, self.ns) / np.maximum(self.ne, 1) ** 2
        return np.where(self.ok, C, np.nan), np.where(self.ok, v, np.nan)

    def perm(self, rng):
        key = self.sid_s + rng.random(len(self.sid_s)); pr = self.order[np.argsort(key, kind="stable")]
        yp = np.empty_like(self.y); sp = np.empty_like(self.s2)
        yp[self.order] = self.y[pr]; sp[self.order] = self.s2[pr]
        return self.contrast(yp, sp)


def wmean(C, v, sc):
    x = C * sc; w = 1 / (v * sc ** 2); ok = np.isfinite(x) & np.isfinite(w)
    return (x[ok] * w[ok]).sum() / w[ok].sum() if ok.any() else np.nan


def station_series(st, kind, val_fn, sig, bin_="YR"):
    s = st[(st.kind == kind) & (st.bin == bin_)].copy()
    s["v"] = val_fn(s)
    d = LT.series_frame(s.assign(lab=s.lab), "v", sig, ["fam", "netsta", "epoch", "bin"], "ccP3")
    return d


def pooled_item(d, scale_fn, rng, nperm=NPERM, group_col="netsta"):
    """d series frame (sid, r, sig, t, fam, netsta) -> pooled decade contrast in scaled units, RE sigma, perm p"""
    if not len(d):
        return {}
    cc = Contr(d.sid, d.r, d.sig ** 2, d.t)
    C, v = cc.contrast()
    first = d.groupby("sid").first()
    sc = scale_fn(first).reindex(range(cc.ns)).to_numpy(dtype=float)
    grp = first[group_col].reindex(range(cc.ns)).to_numpy()
    ok = np.isfinite(C) & np.isfinite(sc)
    pr = L.pool_re((C * sc)[ok], (np.sqrt(v) * np.abs(sc))[ok], grp[ok])
    st_ = np.array([wmean(*cc.perm(rng), sc) for _ in range(nperm)])
    pr["perm_p"] = float((1 + (np.abs(st_) >= abs(pr["mean"])).sum()) / (nperm + 1)) if np.isfinite(pr["mean"]) else np.nan
    pr["perm_sd"] = float(np.nanstd(st_))
    pr["bound2s"] = abs(pr["mean"]) + 2 * max(pr["sig"], pr["perm_sd"])
    return pr


def year_curve(d, scale_fn):
    first = d.groupby("sid").first(); sc = scale_fn(first)
    x = d.r.to_numpy() * sc.reindex(d.sid).to_numpy(); w = 1 / (d.sig.to_numpy() * np.abs(sc.reindex(d.sid).to_numpy())) ** 2
    q = pd.DataFrame(dict(lab=d.lab.to_numpy(), x=x, w=w, cc=d.cc.to_numpy())).dropna()
    return q.groupby("lab").apply(lambda g: pd.Series(dict(value=(g.x * g.w).sum() / g.w.sum(), sig=1 / np.sqrt(g.w.sum()),
                                                           n=len(g), cc=g.cc.median()))).reset_index()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--tag", default=""); ap.add_argument("--regions", default="SVI,WA_NW")
    a = ap.parse_args(); tag = a.tag
    rng = np.random.default_rng(SEED)
    out, curves, joint, majors, splits, spreg = [], [], [], [], [], []
    for reg in a.regions.split(","):
        _, TB, _, _ = L.dirs(reg)
        st = L.load(reg, "lags_station", tag); SP = L.sp_table(reg, tag)
        spf = lambda f: pd.Series([SP.get((a_, b_), np.nan) for a_, b_ in zip(f.fam, f.netsta)], index=f.index, dtype=float)
        smp = lambda s: np.where(s.railS | s.railP3, np.nan, s.lagS - s.lagP3)
        scales = {"ms": lambda f: pd.Series(1.0, index=f.index), "constVpVs": lambda f: -100 / (spf(f) * 1e3),
                  "VsOnly": lambda f: -100 / (L.TS(spf(f)) * 1e3), "VpOnly": lambda f: 100 / (L.TP(spf(f)) * 1e3)}
        # ---- (1) S-P per station and region-wide; (2) absolute S, P3
        for kind in ("data", "ctrl"):
            items = {"S-P3": station_series(st, kind, smp, "sigSmP3"),
                     "absS": station_series(st, kind, lambda s: np.where(s.railS, np.nan, s.lagS), "sigS"),
                     "absP3": station_series(st, kind, lambda s: np.where(s.railP3, np.nan, s.lagP3), "sigP3")}
            if kind == "data":
                items["S-P3s"] = station_series(st, kind, lambda s: np.where(s.railS | s.railP3s, np.nan, s.lagS - s.lagP3s), "sigSmP3s")
                items["S3c-P3"] = station_series(st, kind, lambda s: np.where(s.railS3c | s.railP3, np.nan, s.lagS3c - s.lagP3), "sigS3cmP3")
            for obs, d in items.items():
                if not len(d):
                    continue
                units = ["ms", "constVpVs", "VsOnly", "VpOnly"] if obs.startswith(("S-", "S3c")) else ["ms"]
                for unit in units:
                    pr = pooled_item(d, scales[unit], rng, NPERM if (kind == "data" or unit == "ms") else 500)
                    out.append(dict(region=reg, obs=obs, kind=kind, X="all", unit=unit, **pr))
                    if kind == "data" and unit in ("ms", "constVpVs"):
                        cy = year_curve(d, scales[unit]); cy["region"], cy["obs"], cy["unit"] = reg, obs, unit; curves.append(cy)
                if kind == "data" and obs == "S-P3":
                    for ns, g in d.groupby("netsta"):
                        g = g.copy(); g["sid"] = g.groupby("sid").ngroup()
                        pr = pooled_item(g, scales["constVpVs"], rng, 1000)
                        out.append(dict(region=reg, obs=obs, kind=kind, X=ns, unit="constVpVs", **pr))
        # ---- (3) joint uniform (a, b): S-P, dS, dP3 per-series decade contrasts
        for jname, (rname, rs) in (("headline", JREF[reg]), ("column", JREF2[reg])):
            blocks = []
            d = station_series(st, "data", smp, "sigSmP3")
            if len(d):
                first = d.groupby("sid").first(); sp_ = spf(first)
                blocks.append(("S-P3", d, np.column_stack([-L.TS(sp_.to_numpy()), L.TP(sp_.to_numpy())]), first))
            for obs, fn in (("dS", "pair_corrected"), ("dP3", "pairP_corrected")):
                D_ = pd.read_parquet(f"{TB}/{fn}{tag}.parquet"); D_ = D_[D_.bin == "YR"]
                if obs == "dP3":
                    D_ = D_.assign(ccS_A=D_.ccP3_A, ccS_B=D_.ccP3_B)
                F = L.ref_frame(L.oriented(D_), rs, "rc")
                dd = LT.series_frame(F, "v", "sig", ["fam", "X", "epX", "bin"], "ccX")
                if not len(dd):
                    continue
                first = dd.groupby("sid").first()
                tsx = np.array([L.TS(SP.get((f, x), np.nan)) for f, x in zip(first.fam, first.X)])
                tsr = np.array([np.nanmean([L.TS(SP.get((f, m), np.nan)) for m in rs]) for f in first.fam])
                dT = tsx - tsr
                Dm = np.column_stack([-dT, np.zeros_like(dT)]) if obs == "dS" else np.column_stack([np.zeros_like(dT), -dT / L.VPVS])
                first = first.rename(columns={"X": "netsta"})
                blocks.append((obs, dd, Dm, first))
            ccs = [Contr(b[1].sid, b[1].r, b[1].sig ** 2, b[1].t) for b in blocks]

            def fit(Cv):
                Cs, Vs, Ds = [], [], []
                for (C, v), b in zip(Cv, blocks):
                    Dm = b[2]; ok = np.isfinite(C) & np.isfinite(Dm).all(1) & np.isfinite(v)
                    Cs.append(C[ok] / 1e3); Vs.append(v[ok] / 1e6); Ds.append(Dm[ok])
                C = np.concatenate(Cs); V = np.concatenate(Vs); Dm = np.concatenate(Ds); w = 1 / V
                N = (Dm * w[:, None]).T @ Dm; x = np.linalg.solve(N, (Dm * w[:, None]).T @ C)
                cov = np.linalg.inv(N); chi2 = float((w * (C - Dm @ x) ** 2).sum()); dof = max(len(C) - 2, 1)
                return x, cov, chi2 / dof, len(C)
            x, cov, cd, n = fit([c.contrast() for c in ccs])
            infl = max(1.0, np.sqrt(cd)); sa, sb = np.sqrt(np.diag(cov)) * infl
            sab = np.sqrt(cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]) * infl
            xp = np.array([fit([c.perm(rng) for c in ccs])[0] for _ in range(NPERM)])
            pa = float((1 + (np.abs(xp[:, 0]) >= abs(x[0])).sum()) / (NPERM + 1)); pb = float((1 + (np.abs(xp[:, 1]) >= abs(x[1])).sum()) / (NPERM + 1))
            pab = float((1 + (np.abs(xp[:, 1] - xp[:, 0]) >= abs(x[1] - x[0])).sum()) / (NPERM + 1))
            psd = np.std(xp, 0)
            joint.append(dict(region=reg, fit=jname, ref=rname, n_series=n, blocks=";".join(f"{b[0]}:{len(b[3])}" for b in blocks),
                              a_dVs_pct=x[0] * 100, a_sig_pct=sa * 100, a_perm_p=pa, a_perm_sd_pct=psd[0] * 100,
                              b_dVp_pct=x[1] * 100, b_sig_pct=sb * 100, b_perm_p=pb, b_perm_sd_pct=psd[1] * 100,
                              corr_ab=cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1]), chi2dof=cd,
                              bma_dVpVs_pct=(x[1] - x[0]) * 100, bma_sig_pct=sab * 100, bma_perm_p=pab,
                              a_bound2s_pct=(abs(x[0]) + 2 * max(sa, psd[0])) * 100, b_bound2s_pct=(abs(x[1]) + 2 * max(sb, psd[1])) * 100,
                              bma_bound2s_pct=(abs(x[1] - x[0]) + 2 * max(sab, np.std(xp[:, 1] - xp[:, 0]))) * 100))
        # ---- (4) per ETS: majors, region-wide S-P (constVpVs), demeaned over majors per series
        for kind in ("data", "ctrl"):
            s = st[(st.kind == kind) & (st.bin == "EP") & ~st.lab.str.startswith("B")].copy()
            s["v"] = smp(s)
            d = s[np.isfinite(s.v)].copy()
            d["K"] = d.groupby(["fam", "netsta", "epoch"]).v.transform("size"); d = d[d.K >= 3]
            d["r"] = d.v - d.groupby(["fam", "netsta", "epoch"]).v.transform("mean")
            d["sp"] = [SP.get((f, n), np.nan) for f, n in zip(d.fam, d.netsta)]
            d["x"] = -100 * d.r / (d.sp * 1e3); d["sx"] = 100 * np.maximum(d.sigSmP3, 0.5) * np.sqrt((d.K - 1) / d.K) / (d.sp * 1e3)
            dd_ = d[np.isfinite(d.x) & np.isfinite(d.sx)]
            for split, half in (("station", "A"), ("station", "B"), ("family", "A"), ("family", "B")):   # Amendments 3-4
                if split == "station":
                    hm = dd_.netsta.isin(SPLIT_SETS[reg]) if half == "A" else ~dd_.netsta.isin(SPLIT_SETS[reg])
                else:
                    hm = (dd_.fam % 2 == 0) if half == "A" else (dd_.fam % 2 == 1)
                for e, g in dd_[hm].groupby("lab"):
                    if len(g) < 3:
                        continue
                    pr = L.pool_re(g.x.to_numpy(), g.sx.to_numpy(), g.netsta.to_numpy())
                    splits.append(dict(region=reg, kind=kind, split=split, half=half, episode=e, value=pr["mean"], sig=pr["sig"],
                                       n=len(g), n_fam=g.fam.nunique()))
            if kind == "data":                          # Amendment 4 diagnostic: per (major, station) ms vs SP
                for e, g in dd_.groupby("lab"):
                    q = g.assign(rms=g.r, w=1 / np.maximum(g.sigSmP3, 0.5) ** 2).groupby("netsta").apply(
                        lambda x: pd.Series(dict(v=np.average(x.rms, weights=x.w), sp=np.average(x.sp, weights=x.w), W=x.w.sum())))
                    q = q[np.isfinite(q.v) & np.isfinite(q.sp)]
                    if len(q) >= 3 and q.sp.std() > 0.3:
                        X_ = np.column_stack([np.ones(len(q)), q.sp - q.sp.mean()]); W_ = q.W.to_numpy()
                        N_ = (X_ * W_[:, None]).T @ X_; bb = np.linalg.solve(N_, (X_ * W_[:, None]).T @ q.v.to_numpy())
                        cv_ = np.linalg.inv(N_)
                        spreg.append(dict(region=reg, episode=e, n_st=len(q), offset_ms=bb[0], offset_sig=np.sqrt(cv_[0, 0]),
                                          slope_ms_per_s=bb[1], slope_sig=np.sqrt(cv_[1, 1]), sp_mean=q.sp.mean()))
            for e, g in d.groupby("lab"):
                g = g[np.isfinite(g.x) & np.isfinite(g.sx)]
                if len(g) < 5:
                    continue
                pr = L.pool_re(g.x.to_numpy(), g.sx.to_numpy(), g.netsta.to_numpy())
                majors.append(dict(region=reg, kind=kind, episode=e, date=g.date.min(), n_values=len(g), n_fam=g.fam.nunique(),
                                   n_st=g.netsta.nunique(), dvv_constVpVs_pct=pr["mean"], sig_pct=pr["sig"], chi2dof=pr["chi2dof"],
                                   cc=float(g.ccP3.median())))
    R = pd.DataFrame(out); Cv = pd.concat(curves, ignore_index=True) if curves else pd.DataFrame()
    J = pd.DataFrame(joint); Mj = pd.DataFrame(majors)
    R.to_csv(f"{L.ROOT3}/tables/common_decade{tag}.csv", index=False); Cv.to_csv(f"{L.ROOT3}/tables/common_year_curves{tag}.csv", index=False)
    J.to_csv(f"{L.ROOT3}/tables/common_joint_ab{tag}.csv", index=False); Mj.to_csv(f"{L.ROOT3}/tables/common_majors{tag}.csv", index=False)
    # ---- Amendment 3: split-network coherence across majors
    Sp = pd.DataFrame(splits); coh = []
    from scipy.stats import spearmanr
    for (reg, kind, split), g in Sp.groupby(["region", "kind", "split"]):
        pv = g.pivot_table(index="episode", columns="half", values="value").dropna()
        ps = g.pivot_table(index="episode", columns="half", values="sig").reindex(pv.index)
        for sub in ("all", "sig<=0.05"):
            q = pv if sub == "all" else pv[(ps.A <= 0.05) & (ps.B <= 0.05)]
            if len(q) < 5 or not {"A", "B"} <= set(q.columns):
                continue
            r0 = float(np.corrcoef(q.A, q.B)[0, 1]); rr = np.array([np.corrcoef(rng.permutation(q.A.to_numpy()), q.B)[0, 1] for _ in range(10000)])
            coh.append(dict(region=reg, kind=kind, split=split, subset=sub, n_majors=len(q), r=r0,
                            perm_p_one_sided=float((1 + (rr >= r0).sum()) / 10001), spearman=float(spearmanr(q.A, q.B)[0]),
                            sd_A=float(q.A.std()), sd_B=float(q.B.std())))
    Sp.to_csv(f"{L.ROOT3}/tables/common_majors_split{tag}.csv", index=False)
    Coh = pd.DataFrame(coh); Coh.to_csv(f"{L.ROOT3}/tables/common_majors_coherence{tag}.csv", index=False)
    pd.DataFrame(spreg).to_csv(f"{L.ROOT3}/tables/common_majors_sp_regression{tag}.csv", index=False)
    print(Coh.round(4).to_string())
    # ---- calls
    v = {}
    for reg in a.regions.split(","):
        H = R[(R.region == reg) & (R.X == "all")]
        g = lambda obs, unit, kind="data": H[(H.obs == obs) & (H.unit == unit) & (H.kind == kind)].iloc[0].to_dict() \
            if len(H[(H.obs == obs) & (H.unit == unit) & (H.kind == kind)]) else {}
        smp_ms, aS, aP = g("S-P3", "ms"), g("absS", "ms"), g("absP3", "ms")
        drift = bool(aS and smp_ms and ((abs(aS["mean"]) >= 3 * aS["sig"]) or (abs(aP["mean"]) >= 3 * aP["sig"]))
                     and abs(smp_ms["mean"]) < 2 * smp_ms["sig"])
        dets = {}
        for unit in ("constVpVs", "VsOnly", "VpOnly"):
            x = g("S-P3", unit); c = g("S-P3", unit, "ctrl")
            dets[unit] = dict(mean=x.get("mean"), sig=x.get("sig"), perm_p=x.get("perm_p"), bound2s=x.get("bound2s"),
                              ctrl_p=c.get("perm_p"),
                              single_region_det=bool(x and abs(x["mean"]) >= 3 * x["sig"] and x["perm_p"] < 0.01 and not (c and c["perm_p"] < 0.01)))
        mj = Mj[(Mj.region == reg)]
        mchi = {k: float(((g_.dvv_constVpVs_pct / g_.sig_pct) ** 2).mean()) for k, g_ in mj.groupby("kind")}
        v[reg] = dict(drift_call=drift, absS=aS.get("mean"), absS_sig=aS.get("sig"), absP=aP.get("mean"), absP_sig=aP.get("sig"),
                      SmP_ms=smp_ms.get("mean"), SmP_ms_sig=smp_ms.get("sig"), S_P=dets, majors_chi2=mchi,
                      majors_median_sig_pct=float(mj[mj.kind == "data"].sig_pct.median()) if len(mj) else np.nan,
                      rule3_bound=bool((dets["VsOnly"]["bound2s"] or 9) < 0.03 or (dets["constVpVs"]["bound2s"] or 9) < 0.05))
    regs = a.regions.split(",")
    if len(regs) == 2:
        for unit in ("constVpVs", "VsOnly", "VpOnly"):
            x1, x2 = v[regs[0]]["S_P"][unit], v[regs[1]]["S_P"][unit]
            both = bool(x1["single_region_det"] and x2["single_region_det"] and np.sign(x1["mean"]) == np.sign(x2["mean"]))
            v[f"detection_{unit}"] = both
    v["rule3_PASS"] = bool(any(v[r]["rule3_bound"] for r in regs) or any(v.get(f"detection_{u}", False) for u in ("constVpVs", "VsOnly", "VpOnly")))
    json.dump(v, open(f"{L.ROOT3}/tables/rule3{tag}.json", "w"), indent=1, default=float)
    with pd.option_context("display.width", 250, "display.max_rows", 300):
        print(R[R.X == "all"].round(4).to_string()); print(J.round(4).T.to_string()); print(json.dumps(v, indent=1, default=float))


if __name__ == "__main__":
    main()
