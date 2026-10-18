#!/usr/bin/env python
"""s11_r3_zones -- S11 round 3, task 3 (PLAN §6.3): coarse along-strike deep zones per region.
Design-only forecast first (cell sigma + along-strike kernels, no dv/v): sigma of each zone's composite DeltaD over the
region's majors, per-major, per-year and zone contrasts; run on data only if every zone composite has forecast sigma <= 0.10 %
(3 zones; else 2 zones). When run: zone composites, contrasts, chi2 test of equal zones, scrambled control, leave-one-site-out,
waveform-level zone pulses (P store, raw path, ell40) and a common shallow pulse. dv/v in %.
    S11_REGION=... python s11_r3_zones.py [--workers 12] [--forecast-only]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json
import numpy as np
import pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_r3_common as R
import s11_r3_fit as LF
import s11_r3_inject as J
import s11_paths as P

SIG_RUN = 0.10


def wquant(x, w, q):
    o = np.argsort(x); x, w = np.asarray(x)[o], np.asarray(w)[o]
    c = np.cumsum(w) / w.sum()
    return float(np.interp(q, c, x))


def bins_table():
    B = [pd.read_parquet(f"{R.TB3}/kernel_asbins_r3.parquet")]
    if C.REGION == "SVI" and os.path.exists(f"{R.TB3}/kernel_asbins_r3_ext.parquet"):
        B.append(pd.read_parquet(f"{R.TB3}/kernel_asbins_r3_ext.parquet"))
    B = pd.concat(B, ignore_index=True)
    return B[B.win == C.HEAD].drop_duplicates("series").set_index("series")


def zone_coefs(cells, bounds, B, edges, setup=None):
    """a_d per zone for every cell: bins summed (fractional split of the bin straddling a boundary); open outer bins go to the
    end zones. setup given -> a_d(setup) x ell40 zone share."""
    nb = len(edges) + 1
    lo = np.r_[-np.inf, edges]; hi = np.r_[edges, np.inf]
    zl = np.r_[-np.inf, bounds]; zh = np.r_[bounds, np.inf]
    frac = np.zeros((nb, len(zl)))
    for j in range(nb):
        for z in range(len(zl)):
            if np.isinf(lo[j]) or np.isinf(hi[j]):
                mid = hi[j] - 1 if np.isinf(lo[j]) else lo[j] + 1
                frac[j, z] = float(zl[z] <= mid < zh[z])
            else:
                frac[j, z] = max(0.0, min(hi[j], zh[z]) - max(lo[j], zl[z])) / (hi[j] - lo[j])
    bcols = [f"b{j}" for j in range(nb)]
    key = cells.series if "series" in cells else cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch
    M = B.reindex(key)[bcols].to_numpy() @ frac                                  # (cells, zones)
    miss = ~np.isfinite(M).all(1)
    if miss.any():                                                              # same (family, station), any epoch
        alt = B.groupby(["fam", "netsta"])[bcols].mean()
        for i in np.flatnonzero(miss):
            k = (cells.fam.iloc[i], cells.netsta.iloc[i])
            if k in alt.index:
                M[i] = alt.loc[k].to_numpy() @ frac
    out = cells.copy()
    tot = M.sum(1)
    for z in range(M.shape[1]):
        out[f"z{z}"] = M[:, z] if setup is None else out.a_d.to_numpy() * M[:, z] / np.where(tot != 0, tot, np.nan)
    return out[np.isfinite(out[[f"z{z}" for z in range(M.shape[1])]]).all(1)].reset_index(drop=True)


def comps(F, majors, zs, wts=None):
    """IVW (weights fixed by the forecast) composite vectors per zone"""
    g, W = {}, {}
    for z in zs:
        num = np.zeros(len(F["th"])); den = 0.0; ws = {}
        for e in majors:
            ge = LF.dvec(F, e, z)
            if ge is None:
                continue
            w = (wts[z][e] if wts else 1.0 / float(ge @ F["C"] @ ge))
            ws[e] = w; num += w * ge; den += w
        g[z] = num / den if den > 0 else None; W[z] = ws
    return g, W


def forecast(cells, majors, nz, scale, lg):
    zs = [f"z{z}" for z in range(nz)]
    Lt = LF.temporal(cells, zcols=zs); F = Lt.delta(np.zeros(len(cells))); F["C"] = Lt.Ci * scale
    g, W = comps(F, majors, zs)
    rows = []
    for z in zs:
        sper = [float(np.sqrt(LF.dvec(F, e, z) @ F["C"] @ LF.dvec(F, e, z))) for e in majors if LF.dvec(F, e, z) is not None]
        yrs = {}
        for e in W[z]:
            yrs.setdefault(LF.U.mid[e].year, []).append(e)
        sy = []
        for y, ee in yrs.items():
            gy = sum(W[z][e] * LF.dvec(F, e, z) for e in ee) / sum(W[z][e] for e in ee)
            sy.append(float(np.sqrt(gy @ F["C"] @ gy)))
        rows.append(dict(n_zones=nz, zone=z, composite_sig=float(np.sqrt(g[z] @ F["C"] @ g[z])), per_major_sig_median=float(np.median(sper)),
                         per_year_sig_median=float(np.median(sy)), n_years=len(sy),
                         a_dz_median=float(cells[z].median()), a_dz_sum_weight=float((cells[z] ** 2 / cells.sig ** 2).sum())))
    for i in range(nz):
        for j in range(i + 1, nz):
            gd = g[zs[i]] - g[zs[j]]
            rows.append(dict(n_zones=nz, zone=f"{zs[i]}-{zs[j]}", composite_sig=float(np.sqrt(gd @ F["C"] @ gd))))
    T = pd.DataFrame(rows)
    for r in rows:
        lg("forecast " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    return T, W


def zstats(F, g, zs):
    c = np.array([float(g[z] @ F["th"]) for z in zs])
    Gm = np.array([g[z] for z in zs]); S = Gm @ F["C"] @ Gm.T
    Si = np.linalg.inv(S); one = np.ones(len(zs))
    mu = float(one @ Si @ c / (one @ Si @ one)); r = c - mu
    chi = float(r @ Si @ r); p = float(stats.chi2.sf(chi, len(zs) - 1))
    return c, np.sqrt(np.diag(S)), S, mu, chi, p


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=12); ap.add_argument("--forecast-only", action="store_true")
    a = ap.parse_args()
    lg = R.log("zones"); t0 = time.time()
    majors = R.majors()
    ax = json.load(open(f"{R.TB3}/along_strike_axis.json")); edges = np.array(ax["edges"]); fx = {int(k): v for k, v in ax["fam_x"].items()}
    B = bins_table()
    cells = R.head_cells()
    cells["series"] = cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch
    fw = (1 / cells.sig ** 2).groupby(cells.fam).sum()
    x = np.array([fx[f] for f in fw.index]); w = fw.to_numpy()
    scale = 1.0                                     # design only: formal sigma from the cell errors (PLAN §6.3), no dv/v used
    res = {}; FT = []
    for nz in (3, 2):
        bounds = [wquant(x, w, q) for q in np.arange(1, nz) / nz]
        cz = zone_coefs(cells, bounds, B, edges)
        lg(f"{nz} zones: bounds {np.round(bounds, 1).tolist()} km along strike (az {ax['az_deg']:.0f}); families per zone "
           f"{np.bincount(np.searchsorted(bounds, x), minlength=nz).tolist()}; cells {len(cz)}; zone a_d medians "
           f"{[round(float(cz[f'z{z}'].median()), 3) for z in range(nz)]}")
        T, W = forecast(cz, majors, nz, scale, lg)
        FT.append(T)
        res[nz] = dict(bounds=bounds, W=W, cells=cz, run=bool(T[T.zone.str.fullmatch(r"z\d")].composite_sig.max() <= SIG_RUN))
        if res[nz]["run"]:
            break
    FTT = pd.concat(FT, ignore_index=True); FTT.to_csv(f"{R.TB3}/zones_forecast.csv", index=False)
    nz = 3 if res[3]["run"] else (2 if 2 in res and res[2]["run"] else None)
    R.dump(dict(axis_az=ax["az_deg"], bounds={k: v["bounds"] for k, v in res.items()}, run_zones=nz, sig_rule=SIG_RUN,
                fam_x_weighted=True, written=R.pt_now()), f"{R.TB3}/zones_def.json")
    lg(f"forecast decision: run with {nz} zones" if nz else "forecast decision: NO zone run (forecast sigma > 0.10 % for 3 and 2 zones)")
    if nz is None or a.forecast_only:
        open(f"{R.R3}/ZONES_DONE", "w").write(R.pt_now() + f"\nrun_zones {nz}\n"); return
    zs = [f"z{z}" for z in range(nz)]; cz = res[nz]["cells"]; Wf = res[nz]["W"]; bounds = res[nz]["bounds"]
    # ---------------- data
    Ft = LF.temporal(cz, zcols=zs).solve(); g, _ = comps(Ft, majors, zs, Wf)
    c, s, S, mu, chi, p = zstats(Ft, g, zs)
    rows = [dict(item="data", zone=z, composite=float(c[i]), sig=float(s[i])) for i, z in enumerate(zs)]
    for i in range(nz):
        for j in range(i + 1, nz):
            gd = g[zs[i]] - g[zs[j]]
            rows.append(dict(item="data", zone=f"{zs[i]}-{zs[j]}", composite=float(gd @ Ft["th"]), sig=float(np.sqrt(gd @ Ft["C"] @ gd))))
    rows.append(dict(item="data_equal_zones_test", chi2=chi, p=p, mu=mu, dof=nz - 1, chi2_fit=Ft["chi2"]))
    # per-year column
    py = []
    for z in zs:
        yrs = {}
        for e in Wf[z]:
            yrs.setdefault(LF.U.mid[e].year, []).append(e)
        for y, ee in sorted(yrs.items()):
            gy = sum(Wf[z][e] * LF.dvec(Ft, e, z) for e in ee) / sum(Wf[z][e] for e in ee)
            py.append(dict(zone=z, year=y, majors=",".join(ee), dD=float(gy @ Ft["th"]), sig=float(np.sqrt(gy @ Ft["C"] @ gy))))
    pd.DataFrame(py).to_csv(f"{R.TB3}/zones_peryear.csv", index=False)
    G = []
    # scrambled control
    cc = R.head_cells(control=True); cc["series"] = cc.fam.astype(str) + "|" + cc.netsta + "|" + cc.epoch
    ccz = zone_coefs(cc, bounds, B, edges)
    Fc = LF.temporal(ccz, zcols=zs).solve(); gc_, _ = comps(Fc, majors, zs, Wf)
    c2, s2, _, _, chi2c, pc = zstats(Fc, gc_, zs)
    for i, z in enumerate(zs):
        rows.append(dict(item="control", zone=z, composite=float(c2[i]), sig=float(s2[i])))
    G.append(dict(gate="control", within2sig=bool(np.all(np.abs(c2) < 2 * s2)), chi2=chi2c, p=pc,
                  passed=bool(np.all(np.abs(c2) < 2 * s2) and pc > 0.05), n_cells=len(ccz)))
    # leave-one-site-out
    lo_ok, sign_ok, lrows = True, True, []
    for st in sorted(cz.site.unique()):
        sub = cz[cz.site != st].reset_index(drop=True)
        if sub.episode.nunique() < 5:
            continue
        Fs = LF.temporal(sub, zcols=zs).solve(); gs_, _ = comps(Fs, majors, zs, Wf)
        if any(gs_[z] is None for z in zs):
            continue
        cs_, ss_, _, _, chs, ps = zstats(Fs, gs_, zs)
        ok = bool(np.all(np.abs(cs_ - c) < 2 * ss_)); lo_ok &= ok
        dsg = [np.sign(cs_[i] - cs_[j]) == np.sign(c[i] - c[j]) for i in range(nz) for j in range(i + 1, nz)]
        sign_ok &= bool(all(dsg))
        lrows.append(dict(site=st, **{f"{z}": float(cs_[i]) for i, z in enumerate(zs)}, **{f"{z}_sig": float(ss_[i]) for i, z in enumerate(zs)},
                          chi2=chs, p=ps, within2sig=ok, contrast_signs_kept=bool(all(dsg))))
    pd.DataFrame(lrows).to_csv(f"{R.TB3}/zones_loso.csv", index=False)
    G.append(dict(gate="loso", passed=bool(lo_ok), contrast_signs_kept_all=bool(sign_ok), n_sites=len(lrows)))
    # waveform-level zone pulses (P store, raw path)
    sk = cz.drop_duplicates("series")
    coef = {k: (a_s, dict(zip(zs, [r[z] for z in zs]))) for k, a_s, r in zip(sk.series, sk.a_s, sk.to_dict("records"))}

    def runs_fn(key, st, unit):
        if key not in coef:
            return {"nat": None}
        a_s, adz = coef[key]
        ep = np.array([unit.get(d, None) for d in st["days"]], dtype=object); mall = np.isin(ep, majors)
        out = {"nat": None}
        if mall.any():
            for z in zs:
                out[f"Z|{z}"] = -(0.003 * adz[z]) * mall
            out["Psall"] = -(0.001 * a_s) * mall
        return out
    unit = P.unit_map()
    vals = J.run(set(sk.series), runs_fn, J.est_raw, J.store_paths_P(), a.workers, lg, unit=unit)
    dn = J.cell_vector(vals, cz, "nat"); ok = np.isfinite(dn)
    lg(f"natural-from-store vs measured: max |diff| {np.nanmax(np.abs(dn - cz.d.to_numpy())):.1e}")
    cs2 = cz[ok].reset_index(drop=True); Lt = LF.temporal(cs2, zcols=zs); F0 = Lt.delta(np.zeros(len(cs2)))
    g0, _ = comps(F0, majors, zs, Wf)
    inj = []
    rec_ok, leak_ok = True, True
    for zi in zs:
        di = J.cell_vector(vals, cz, f"Z|{zi}")[ok]; xx = np.where(np.isfinite(di), di - dn[ok], 0.0); th = Lt.Op @ xx
        for zj in zs:
            v = float(g0[zj] @ th) / 0.3
            inj.append(dict(pulse=f"deep_{zi}", zone=zj, ratio=v))
            if zi == zj:
                rec_ok &= bool(0.7 <= v <= 1.3)
            else:
                leak_ok &= bool(abs(v) <= 0.30)
    di = J.cell_vector(vals, cz, "Psall")[ok]; xx = np.where(np.isfinite(di), di - dn[ok], 0.0); th = Lt.Op @ xx
    sh_ok = True
    for zj in zs:
        inj.append(dict(pulse="shallow_all", zone=zj, ratio=float(g0[zj] @ th) / 0.1))
    for i in range(nz):
        for j in range(i + 1, nz):
            v = float((g0[zs[i]] - g0[zs[j]]) @ th) / 0.1; inj.append(dict(pulse="shallow_all", zone=f"{zs[i]}-{zs[j]}", ratio=v))
            sh_ok &= bool(abs(v) <= 0.35)
    pd.DataFrame(inj).to_csv(f"{R.TB3}/zones_injection.csv", index=False)
    G.append(dict(gate="zone_pulse_recovery", passed=bool(rec_ok))); G.append(dict(gate="zone_pulse_crossleak", passed=bool(leak_ok)))
    G.append(dict(gate="shallow_pulse_contrast_leak", passed=bool(sh_ok)))
    allg = all(x["passed"] for x in G if x["passed"] is not None)
    G.append(dict(gate="WORDING_differs", passed=bool(p < 0.01 and allg and sign_ok), p_data=p))
    # columns: task-1 Q3 winner kernel, adopted task-2 estimator
    vj = f"{R.TB3}/kerr_verdict.json"
    if os.path.exists(vj):
        v = json.load(open(vj)); w3 = v["winner"].get("q3")
        if w3:
            import s11_r3_kerr as KR
            st_ = KR.setups(v["lapse"])[w3]
            cw = zone_coefs(R.with_setup(cells, st_), bounds, B, edges, setup=st_)
            Fw = LF.temporal(cw, zcols=zs).solve(); gw, _ = comps(Fw, majors, zs, Wf); cwv, swv, _, _, chw, pw = zstats(Fw, gw, zs)
            for i, z in enumerate(zs):
                rows.append(dict(item=f"kernel_{w3}", zone=z, composite=float(cwv[i]), sig=float(swv[i])))
            rows.append(dict(item=f"kernel_{w3}_equal_zones_test", chi2=chw, p=pw))
    qv = f"{R.TB3}/r3q3_verdict.json"
    if os.path.exists(qv) and json.load(open(qv)).get("adopted"):
        import s11_r3_q3 as Q
        cfg = json.load(open(qv))["adopted"]; ce, _ = Q.cells_for(cfg, Q.r3_tables())
        cez = zone_coefs(ce, bounds, B, edges)
        Fe = LF.temporal(cez, zcols=zs).solve(); ge, _ = comps(Fe, majors, zs, Wf); cev, sev, _, _, che, pe = zstats(Fe, ge, zs)
        for i, z in enumerate(zs):
            rows.append(dict(item=f"estimator_{cfg}", zone=z, composite=float(cev[i]), sig=float(sev[i])))
        rows.append(dict(item=f"estimator_{cfg}_equal_zones_test", chi2=che, p=pe))
    pd.DataFrame(rows).to_csv(f"{R.TB3}/zones_result.csv", index=False)
    pd.DataFrame(G).to_csv(f"{R.TB3}/zones_gates.csv", index=False)
    for r in rows + G:
        lg("zones " + ", ".join(f"{k} {v:.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in r.items()))
    lg(f"zones done {time.time() - t0:.0f}s")
    open(f"{R.R3}/ZONES_DONE", "w").write(R.pt_now() + f"\nrun_zones {nz}\n")


if __name__ == "__main__":
    main()
