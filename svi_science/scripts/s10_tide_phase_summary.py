#!/usr/bin/env python
"""s10_tide_phase_summary -- S10 follow-up part A, pooling / fits / controls / rules / bounds (PLAN.md A3-A6).
Reads tables/phase_series/*.npz (s10_tide_phase.py). Writes tables/phase_bins.csv, phase_fits.csv, phase_rules.csv,
phase_injection.csv, phase_population.csv, phase_bounds.csv, phase_result.json.
    python s10_tide_phase_summary.py [--indir DIR] [--tag _test]
"""
from __future__ import annotations
import os, sys, argparse, json, time
import numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_tide_common as T
import s10_tide_phase as TP

NB = TP.NB
CEN = np.radians(np.arange(NB) * 45.0)
XPH = np.c_[np.ones(NB), np.cos(CEN), np.sin(CEN)]
FAMSETS = ["10", "87", "17", "all"]


def fit(S1, S0, X):
    """batched WLS: S1, S0 (nd, nb); X (nb, p) -> coef (nd, p) (NaN if < p+1 bins with weight)."""
    S1 = np.asarray(S1, np.float64); S0 = np.asarray(S0, np.float64)
    N = np.einsum("db,bi,bj->dij", S0, X, X)
    r = np.einsum("db,bi->di", S1, X)
    coef = np.einsum("dij,dj->di", np.linalg.pinv(N), r)
    ok = (S0 > 0).sum(1) >= X.shape[1] + 1
    coef[~ok] = np.nan
    return coef


def load(indir):
    acc = {}       # (famset, level, key) -> sum array
    info = []
    for f in sorted(os.listdir(indir)):
        if not f.endswith(".npz"):
            continue
        z = np.load(f"{indir}/{f}")
        fam, ns = str(int(z["fam"])), str(z["netsta"])
        info.append(dict(fam=int(fam), netsta=ns, epoch=str(z["epoch"]), n=int(z["n"]), n_days=int(z["n_days"])))
        for k in z.files:
            if "~" not in k:
                continue
            a = z[k].astype(np.float64)
            for fs in (fam, "all"):
                for lev in (ns, "pooled"):
                    key = (fs, lev, k)
                    if key in acc:
                        acc[key] = acc[key] + a
                    else:
                        acc[key] = a.copy()
    return acc, pd.DataFrame(info)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--indir", default=TP.OUTD); ap.add_argument("--tag", default="")
    a = ap.parse_args()
    log = S.Log(f"{T.LOGD}/summary{a.tag}.log")
    meta = np.load(f"{T.TB}/phase_bins_meta.npz")
    fams_meta = [int(x) for x in meta["fams"]]
    acc, info = load(a.indir)
    info.to_csv(f"{T.TB}/phase_series_info{a.tag}.csv", index=False)
    log(f"loaded {len(info)} series, {len(acc)} accumulators")
    stations = sorted(info.netsta.unique())
    keys = sorted({k for (_, _, k) in acc})
    fit_rows, bin_rows = [], []

    def xval(var, fs):
        if fs == "all":
            return None
        return meta[f"xmean_{var}"][fams_meta.index(int(fs))]

    # ---------------- fits for every accumulator
    store = {}
    for (fs, lev, k), arr in acc.items():
        parts = k.split("~")
        if parts[0].startswith("inj") or parts[1] == "lbo":
            continue
        var = parts[0].replace("c2:", ""); ctrl = parts[0].startswith("c2:"); mode, des, win = parts[1], parts[2], parts[3]
        S1, S0, S3, NC = arr[:, 0], arr[:, 1], arr[:, 2], arr[:, 3]
        if var in TP.PHASE_VARS:
            co = fit(S1, S0, XPH)
            p_, q_ = co[:, 1], co[:, 2]
            A = np.hypot(p_, q_)
            store[(fs, lev, k)] = dict(kind="phase", p=p_, q=q_, A=A)
        else:
            if fs == "all":
                # joint fit: family intercepts + common slope, from the per-family accumulators at this level
                S1s, S0s, Xs = [], [], []
                for j, f in enumerate(("10", "87", "17")):
                    kk = (f, lev, k)
                    if kk not in acc:
                        continue
                    S1s.append(acc[kk][:, 0]); S0s.append(acc[kk][:, 1])
                    x = meta[f"xmean_{var}"][fams_meta.index(int(f))]
                    Xi = np.zeros((NB, 4)); Xi[:, j] = 1; Xi[:, 3] = x; Xs.append(Xi)
                X = np.vstack(Xs); S1c = np.hstack(S1s); S0c = np.hstack(S0s)
                used = X[:, :3].sum(0) > 0
                X = np.c_[X[:, :3][:, used], X[:, 3]]
                co = fit(S1c, S0c, X); slope = co[:, -1]
            else:
                co = fit(S1, S0, np.c_[np.ones(NB), xval(var, fs)]); slope = co[:, 1]
            store[(fs, lev, k)] = dict(kind="value", slope=slope)
        # bins (observed) with control SD
        with np.errstate(invalid="ignore", divide="ignore"):
            eb = S1 / np.where(S0 > 0, S0, np.nan)
            ccb = S3[0] / np.where(S0[0] > 0, S0[0], np.nan)
        for b in range(NB):
            bin_rows.append(dict(famset=fs, level=lev, var=("c2:" if ctrl else "") + var, mode=mode, design=des, win=win, bin=b,
                                 bin_center=(b * 45.0 if var in TP.PHASE_VARS else np.nan),
                                 x_mean=(xval(var, fs)[b] if (var in TP.VALUE_VARS and fs != "all") else np.nan),
                                 eps_pct=100 * eb[0, b], dvv_pct=-100 * eb[0, b], cc=ccb[b], w=S0[0, b], n_cells=NC[0, b],
                                 ctrl_sd_pct=100 * np.nanstd(eb[1:, b])))
    B = pd.DataFrame(bin_rows)
    B = B[B["mode"] == "rep"] if "mode" in B else B
    B.to_csv(f"{T.TB}/phase_bins{a.tag}.csv", index=False)

    # ---------------- rules per (famset, var, mode, design) ; ctrl prefix handled
    rules = []
    for (fs, lev, k), st in store.items():
        if lev != "pooled":
            continue
        parts = k.split("~")
        if parts[3] != "24":
            continue
        var = parts[0]; mode, des = parts[1], parts[2]
        k35 = "~".join(parts[:3] + ["35"])
        st35 = store.get((fs, lev, k35))
        r = dict(famset=fs, var=var, mode=mode, design=des)
        if st["kind"] == "phase":
            A = st["A"]; phi = np.arctan2(st["q"][0], st["p"][0])
            r.update(A_eps_pct=100 * A[0], phase_eps_max_deg=np.degrees(phi) % 360, phase_dvv_max_deg=(np.degrees(phi) + 180) % 360,
                     A_ctrl_med_pct=100 * np.nanmedian(A[1:]), A_ctrl99_pct=100 * np.nanpercentile(A[1:], 99),
                     p_A=(1 + np.sum(A[1:] >= A[0])) / (1 + np.isfinite(A[1:]).sum()),
                     sd_cos_ctrl_pct=100 * np.nanstd(st["p"][1:]))
            crit1 = A[0] > np.nanpercentile(A[1:], 99)
            u = np.array([np.cos(phi), np.sin(phi)])
            sites = set(); sta_list = []
            for ns in stations:
                s2 = store.get((fs, ns, k))
                if s2 is None or not np.isfinite(s2["p"][0]):
                    continue
                proj = s2["p"] * u[0] + s2["q"] * u[1]
                ok = proj[0] > np.nanpercentile(proj[1:], 95)
                sta_list.append(f"{ns}:{100 * proj[0]:+.3f}{'*' if ok else ''}")
                if ok:
                    sites.add(T.site(ns))
            r["stations_proj_pct"] = " ".join(sta_list)
            crit2 = len(sites) >= 2
            if st35 is not None and np.isfinite(st35["p"][0]):
                proj35 = st35["p"][0] * u[0] + st35["q"][0] * u[1]
                r.update(A35_eps_pct=100 * st35["A"][0], proj35_on_24_pct=100 * proj35,
                         p_A35=(1 + np.sum(st35["A"][1:] >= st35["A"][0])) / (1 + np.isfinite(st35["A"][1:]).sum()))
                crit3 = proj35 > 0
            else:
                crit3 = False
        else:
            s = st["slope"]
            sg = np.sign(s[0]) if np.isfinite(s[0]) else 1.0
            r.update(slope_eps_pct_per_unit=100 * s[0], slope_dvv_pct_per_unit=-100 * s[0],
                     slope_ctrl_sd_pct=100 * np.nanstd(s[1:]), slope_ctrl99_abs_pct=100 * np.nanpercentile(np.abs(s[1:]), 99),
                     p_slope=(1 + np.sum(np.abs(s[1:]) >= abs(s[0]))) / (1 + np.isfinite(s[1:]).sum()))
            crit1 = abs(s[0]) > np.nanpercentile(np.abs(s[1:]), 99)
            sites = set(); sta_list = []
            for ns in stations:
                s2 = store.get((fs, ns, k))
                if s2 is None or not np.isfinite(s2["slope"][0]):
                    continue
                pr = s2["slope"] * sg
                ok = pr[0] > np.nanpercentile(pr[1:], 95)
                sta_list.append(f"{ns}:{100 * s2['slope'][0]:+.4f}{'*' if ok else ''}")
                if ok:
                    sites.add(T.site(ns))
            r["stations_slope_pct"] = " ".join(sta_list)
            crit2 = len(sites) >= 2
            if st35 is not None and np.isfinite(st35["slope"][0]):
                r["slope35_eps_pct_per_unit"] = 100 * st35["slope"][0]
                crit3 = np.sign(st35["slope"][0]) == sg
            else:
                crit3 = False
        r.update(crit1_pooled_gt_ctrl99=bool(crit1), crit2_sites=len(sites), crit2_sites_list=";".join(sorted(sites)),
                 crit3_same_sign_35=bool(crit3), PASS=bool(crit1 and crit2 and crit3))
        rules.append(r)
    RU = pd.DataFrame(rules).sort_values(["famset", "design", "var", "mode"])
    RU.to_csv(f"{T.TB}/phase_rules{a.tag}.csv", index=False)

    # ---------------- injection recovery (family 10, pooled)
    inj = []
    for (fs, lev, k), arr in acc.items():
        if not (fs == "10" and lev in ("pooled",) + tuple(stations) and k.startswith("inj")):
            continue
        aa, des, win = k.split("~")
        amp = float(aa.replace("inj", ""))
        nat = acc.get((fs, lev, f"th_udss~rep~{des}~{win}"))
        if nat is None:
            continue
        ci = fit(arr[:1, 0], arr[:1, 1], XPH)[0]; cn = fit(nat[:1, 0], nat[:1, 1], XPH)[0]
        inj.append(dict(level=lev, a_pct=100 * amp, design=des, win=win, recovery=(ci[1] - cn[1]) / amp,
                        quad_leak=(ci[2] - cn[2]) / amp))
    INJ = pd.DataFrame(inj)
    INJ.to_csv(f"{T.TB}/phase_injection{a.tag}.csv", index=False)

    # ---------------- leave-bin-out D1 column (observed) fits
    lbo = []
    for (fs, lev, k), arr in acc.items():
        if lev != "pooled" or "~lbo~" not in k:
            continue
        var, _, des, win = k.split("~")
        if var in TP.PHASE_VARS:
            co = fit(arr[:1, 0], arr[:1, 1], XPH)[0]
            lbo.append(dict(famset=fs, var=var, win=win, A_eps_pct=100 * np.hypot(co[1], co[2]),
                            phase_dvv_max_deg=(np.degrees(np.arctan2(co[2], co[1])) + 180) % 360))
        elif fs != "all":
            co = fit(arr[:1, 0], arr[:1, 1], np.c_[np.ones(NB), xval(var, fs)])[0]
            lbo.append(dict(famset=fs, var=var, win=win, slope_dvv_pct_per_unit=-100 * co[1]))
    pd.DataFrame(lbo).to_csv(f"{T.TB}/phase_lbo{a.tag}.csv", index=False)

    # ---------------- bounds (PLAN A6)
    m2amp = pd.read_parquet(f"{T.TB}/rows_tide.parquet").groupby("fam").m2_udss_amp_Pa.mean() / 1e3
    bounds = []
    for _, r in RU[(RU["mode"] == "rep") & RU["var"].isin(["th_udss", "udss"])].iterrows():
        rec = INJ[(INJ.level == "pooled") & (INJ.design == r.design) & (INJ.win == "24") & (np.isclose(INJ.a_pct, 0.1))]
        rc = float(rec.recovery.iloc[0]) if len(rec) else np.nan
        if r["var"] == "th_udss":
            up = (r.A_eps_pct + 2.58 * r.sd_cos_ctrl_pct)
            amp_kpa = float(m2amp.get(int(r.famset), np.nan)) if r.famset != "all" else float(m2amp.mean())
            bounds.append(dict(famset=r.famset, design=r.design, quantity="tidal dv/v amplitude (M2 UDSS phase), %",
                               observed=r.A_eps_pct, bound99_raw=up, recovery=rc, bound99=up / rc if rc > 0 else np.nan,
                               per_kPa_bound=(up / rc / amp_kpa) if rc > 0 else np.nan, m2_udss_amp_kPa=amp_kpa))
        else:
            up = abs(r.slope_eps_pct_per_unit) + 2.58 * r.slope_ctrl_sd_pct
            bounds.append(dict(famset=r.famset, design=r.design, quantity="dv/v per kPa UDSS (value bins), %/kPa",
                               observed=r.slope_dvv_pct_per_unit, bound99_raw=up, recovery=rc,
                               bound99=up / rc if rc > 0 else np.nan))
    pd.DataFrame(bounds).to_csv(f"{T.TB}/phase_bounds{a.tag}.csv", index=False)

    # ---------------- population check (iv)
    try:
        pop = population(a)
        pop.to_csv(f"{T.TB}/phase_population{a.tag}.csv", index=False)
    except Exception as e:
        log(f"population check failed: {e}")
    log("summary done")
    if not a.tag:
        open(f"{T.OUT}/SUMMARY_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


def population(a):
    """per family, all series pooled: mean P, within-day log10 S amplitude, realign lag / cc per M2-UDSS phase bin and
    UDSS value bin; cos fit / slope vs the same 1,000 control draws."""
    R, codes, edges, xmean, fams = TP.build_codes(1000)
    pos = {r: i for i, r in enumerate(R.row_id)}
    sz = pd.read_parquet(f"{S.ROOT}/svi_science/results/s6/sizes.parquet", columns=["det_id", "dM"])
    R = R.merge(sz.drop_duplicates("det_id"), on="det_id", how="left")
    W = []
    for f in sorted(os.listdir(T.CACHE)):
        z = np.load(f"{T.CACHE}/{f}")
        ri = np.array([pos[r] for r in z["row_id"]])
        la = pd.Series(z["logamp"]).groupby(z["date"]).transform(lambda s: s - s.mean()).to_numpy()
        W.append(pd.DataFrame(dict(ri=ri, fam=int(z["fam"]), lag_ms=1e3 * z["lag"], rcc=z["rcc"], logamp_dm=la)))
    W = pd.concat(W, ignore_index=True)
    data = (R.series == "data").to_numpy()
    out = []
    for fam in T.FAMS:
        rows_f = np.flatnonzero(data & (R.fam == fam).to_numpy())
        wf = W[W.fam == fam]
        metrics = {"P": (rows_f, R.P.to_numpy()[rows_f]), "dM_S6": (rows_f, R.dM.to_numpy()[rows_f]),
                   "lag_ms": (wf.ri.to_numpy(), wf.lag_ms.to_numpy()), "realign_cc": (wf.ri.to_numpy(), wf.rcc.to_numpy()),
                   "logamp_withinday": (wf.ri.to_numpy(), wf.logamp_dm.to_numpy())}
        for var in ("th_udss", "udss"):
            CODES = codes[var]["rep"]
            for mn, (ri, x) in metrics.items():
                ok = np.isfinite(x); ri, x = ri[ok], x[ok]
                if len(x) < 20:
                    continue
                C = CODES[:, ri].astype(np.int64)
                off = (np.arange(C.shape[0])[:, None] * NB + C).ravel()
                s1 = np.bincount(off, weights=np.tile(x, C.shape[0]), minlength=C.shape[0] * NB).reshape(-1, NB)
                s0 = np.bincount(off, minlength=C.shape[0] * NB).reshape(-1, NB).astype(float)
                mb = s1 / np.where(s0 > 0, s0, np.nan)
                if var == "th_udss":
                    co = fit(s1, s0, XPH); stat = np.hypot(co[:, 1], co[:, 2])
                else:
                    co = fit(s1, s0, np.c_[np.ones(NB), xmean[var][fams.index(fam)]]); stat = np.abs(co[:, 1])
                out.append(dict(fam=fam, var=var, metric=mn, n=len(x), bin_means=" ".join(f"{v:.4g}" for v in mb[0]),
                                stat_obs=stat[0], stat_ctrl99=np.nanpercentile(stat[1:], 99),
                                p=(1 + np.sum(stat[1:] >= stat[0])) / (1 + np.isfinite(stat[1:]).sum()),
                                flag_population_effect=bool((1 + np.sum(stat[1:] >= stat[0])) / (1 + np.isfinite(stat[1:]).sum()) < 0.01)))
    return pd.DataFrame(out)


if __name__ == "__main__":
    main()
