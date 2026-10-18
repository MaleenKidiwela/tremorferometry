#!/usr/bin/env python
"""s11_r7_zones -- S11 round 7, task 2 (PLAN.md §11.2): along-strike deep zones with the E4 model, one major at a time,
pooled, and slipped-vs-not-slipped zones from S3 patch coverage. Zone shares from the round-3 ell40 along-strike bins;
time-of-day scrambled control cells; model-consistent zone pulses and the round-4 waveform pulses.
    S11_REGION=... python s11_r7_zones.py [--test N]
"""
from __future__ import annotations
import os
for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys, argparse, time, json, itertools
import numpy as np
import pandas as pd
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_r7_common as Q
import s11_r5_common as R5
import s11_r4_common as R4
import s11_r3_common as R3
import s11_r3_zones as RZ

RES = f"{C.ROOT}/svi_science/results"
NPERM = 5000


# ------------------------------------------------------------------------------ zone shares
def zone_share(cells, bounds, B, edges):
    """(cells, nz) share of the ell40 deep coefficient per zone; NaN where the series has no bins"""
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
    M = B.reindex(cells.series.to_numpy())[bcols].to_numpy() @ frac
    tot = M.sum(1, keepdims=True)
    with np.errstate(invalid="ignore", divide="ignore"):
        Sh = M / np.where(tot != 0, tot, np.nan)
    return Sh


def shares_for(cs, bounds, B, edges, lg=None):
    Sh = zone_share(cs, bounds, B, edges)
    ok = np.isfinite(Sh).all(1)
    fam_mean = pd.DataFrame(Sh[ok]).groupby(cs.fam.to_numpy()[ok]).mean()
    nfill = 0
    for i in np.flatnonzero(~ok):
        f = cs.fam.iloc[i]
        if f in fam_mean.index:
            Sh[i] = fam_mean.loc[f].to_numpy(); nfill += 1
    ok2 = np.isfinite(Sh).all(1)
    if lg:
        lg(f"zone shares: {int(ok.sum())} cells from bins, {nfill} filled with the family mean (far), {int((~ok2).sum())} dropped")
    return Sh, ok2


# ------------------------------------------------------------------------------ slipped labels
def axis_proj():
    ax = json.load(open(f"{R3.TB3}/along_strike_axis.json"))
    st = np.array(ax["strike"]); c0 = np.array(ax["centroid_km"])
    if C.REGION == "SVI":
        lat0, lon0 = 48.5, -123.6
    else:
        fl = pd.read_csv(f"{C.MTB}/families.csv"); lat0, lon0 = float(fl.lat.mean()), float(fl.lon.mean())
    def f(lat, lon):
        xy = np.column_stack([(np.asarray(lon) - lon0) * 111.19 * np.cos(np.radians(lat0)), (np.asarray(lat) - lat0) * 111.19])
        return (xy - c0) @ st
    return f, ax


def s3_labels(majors, bounds, lg):
    """zone z slipped in S3 episode e iff >= 50 % of its patches have >= 3 S3 headline events -> {major: [bool per zone]}"""
    import s11_r7_front as FR
    f, _ = axis_proj()
    EP = pd.read_csv(f"{RES}/s3/episodes_all.csv", parse_dates=["w0", "w1"])
    if C.REGION == "SVI":
        ev = pd.read_parquet(f"{RES}/s3/events_svi.parquet")
        ev = ev[(ev.series.astype(str) == "data") & (ev.P_same >= 0.75) & ev.keep_xfam8 & ev.counted_once].copy()
        ev["patch"] = ev.group.astype(int); EP = EP[EP.region == "SVI"]
    else:
        ev = pd.read_parquet(f"{RES}/s3/events_nw.parquet")
        ev = ev[(ev.series.astype(str) == "data") & ev.selected_v2 & ev.counted_once].copy()
        ev["patch"] = np.floor(ev.as_km / 10.0).astype(int)
        sh = ev.patch.value_counts() / len(ev); ev = ev[ev.patch.isin(sh.index[sh >= 0.01])]
        EP = EP[EP.region == "NW"]
    ev = ev[np.isfinite(ev.lat) & np.isfinite(ev.lon)]
    ev["x"] = f(ev.lat.to_numpy(), ev.lon.to_numpy())
    pz = ev.groupby("patch").x.median()
    zone_of_patch = pd.Series(np.searchsorted(bounds, pz.to_numpy()), index=pz.index)
    nz = len(bounds) + 1
    EPi = EP.set_index("episode")
    out, rows = {}, []
    for e3, g in ev.groupby("ep"):
        if e3 not in EPi.index:
            continue
        e = FR.s11_major_for(e3, EPi.loc[e3, "w0"], EPi.loc[e3, "w1"], majors)
        if e is None:
            continue
        cnt = g.patch.value_counts()
        lab = []
        for z in range(nz):
            pp = zone_of_patch.index[zone_of_patch == z]
            if not len(pp):
                lab.append(None); continue
            cov = float(np.mean([cnt.get(p, 0) >= 3 for p in pp]))
            lab.append(bool(cov >= 0.5))
            rows.append(dict(source="S3", s3_episode=e3, major=e, zone=f"z{z}", n_patches=len(pp), coverage=cov, slipped=cov >= 0.5,
                             n_events=int(g.patch.isin(pp).sum())))
        out[e] = lab
    lg(f"S3 slipped labels: {len(out)} majors; patches per zone {[int((zone_of_patch == z).sum()) for z in range(nz)]}")
    return out, rows


def fam_labels(majors, bounds, fams, lg):
    ax = json.load(open(f"{R3.TB3}/along_strike_axis.json")); fx = {int(k): v for k, v in ax["fam_x"].items()}
    A = pd.read_csv(f"{C.MTB}/active_days.csv")
    if C.REGION == "SVI" and os.path.exists(f"{C.TB}/ext/active_days.csv"):
        A = pd.concat([A, pd.read_csv(f"{C.TB}/ext/active_days.csv")], ignore_index=True).drop_duplicates(["fam", "date"])
    fz = {f: int(np.searchsorted(bounds, fx[f])) for f in fams if f in fx}
    pe = A[A.episode.isin(majors)].groupby(["episode", "fam"]).sumP.sum()
    nz = len(bounds) + 1
    out, rows = {}, []
    for e in majors:
        lab = []
        for z in range(nz):
            fz_ = [f for f, zz in fz.items() if zz == z]
            if not fz_:
                lab.append(None); continue
            act = [pe.get((e, f), 0.0) >= 1 for f in fz_]
            cov = float(np.mean(act)); lab.append(bool(cov >= 0.5))
            rows.append(dict(source="F", major=e, zone=f"z{z}", n_patches=len(fz_), coverage=cov, slipped=cov >= 0.5))
        out[e] = lab
    return out, rows


# ------------------------------------------------------------------------------ fits and statistics
def fit(cs, Sh, nz, d=None):
    cols = [("D", f"z{k}", cs.a_d.to_numpy(float) * Sh[:, k]) for k in range(nz)]
    L = Q.e4_model(cs, cols)
    F = L.solve(cs.d.to_numpy() if d is None else d)
    return L, F


def zst(F, gs):
    c = np.array([float(g @ F["th"]) for g in gs]); Gm = np.array(gs); S = Gm @ F["C"] @ Gm.T
    Si = np.linalg.pinv(S); one = np.ones(len(gs))
    mu = float(one @ Si @ c / (one @ Si @ one)); r = c - mu
    chi = float(r @ Si @ r); p = float(stats.chi2.sf(chi, len(gs) - 1))
    return c, np.sqrt(np.diag(S)), chi, p


def per_major(L, F, majors, nz):
    rows, gz = [], {}
    for e in majors:
        gs = [Q.dvec_g(F, L.eps, e, "D", f"z{k}") for k in range(nz)]
        if any(g is None for g in gs):
            continue
        gz[e] = gs
        c, s, chi, p = zst(F, gs)
        r = dict(major=e, chi2=chi, p=p)
        for k in range(nz):
            r[f"z{k}"] = c[k]; r[f"z{k}_sig"] = s[k]
        for i, j in itertools.combinations(range(nz), 2):
            gd = gs[i] - gs[j]; v, sg = Q.val_sig(F, gd); r[f"z{i}-z{j}"] = v; r[f"z{i}-z{j}_sig"] = sg
        rows.append(r)
    return pd.DataFrame(rows), gz


def pooled(F, gz, nz):
    """IVW composite per zone over majors, contrasts, equal-zones test"""
    comp = []
    for k in range(nz):
        num = np.zeros(len(F["th"])); den = 0.0
        for e, gs in gz.items():
            g = gs[k]; w = 1.0 / float(g @ F["C"] @ g); num += w * g; den += w
        comp.append(num / den)
    c, s, chi, p = zst(F, comp)
    r = dict(chi2=chi, p=p, n_majors=len(gz))
    for k in range(nz):
        r[f"z{k}"] = c[k]; r[f"z{k}_sig"] = s[k]
    for i, j in itertools.combinations(range(nz), 2):
        v, sg = Q.val_sig(F, comp[i] - comp[j]); r[f"z{i}-z{j}"] = v; r[f"z{i}-z{j}_sig"] = sg
    return r, comp


def slip_contrast(F, gz, labels, rng, nperm=NPERM):
    """per major with both classes: mean ΔD(slipped) - mean ΔD(not); pooled IVW; label permutation within majors"""
    rows, gl = [], []
    for e, gs in gz.items():
        lab = labels.get(e)
        if lab is None:
            continue
        idx = [k for k, l in enumerate(lab) if l is not None]
        sl = [k for k in idx if lab[k]]; ns = [k for k in idx if not lab[k]]
        if not sl or not ns:
            continue
        g = np.mean([gs[k] for k in sl], 0) - np.mean([gs[k] for k in ns], 0)
        v, s = Q.val_sig(F, g)
        rows.append(dict(major=e, slipped=";".join(f"z{k}" for k in sl), not_slipped=";".join(f"z{k}" for k in ns), contrast=v, sig=s))
        gl.append((e, g, [gs[k] for k in idx], len(sl)))
    if not gl:
        return pd.DataFrame(rows), None
    w = np.array([1.0 / float(g @ F["C"] @ g) for _, g, _, _ in gl])
    gp = sum(wi * g for wi, (_, g, _, _) in zip(w, gl)) / w.sum()
    v, s = Q.val_sig(F, gp)
    # permutation: labels permuted among the labelled zones of each major (count of slipped kept)
    zv = [np.array([float(g @ F["th"]) for g in gl_[2]]) for gl_ in gl]
    ws = w / w.sum()
    perm = np.zeros(nperm)
    for t in range(nperm):
        acc = 0.0
        for i, (_, _, gg, nsl) in enumerate(gl):
            o = rng.permutation(len(gg)); sl = o[:nsl]; ns = o[nsl:]
            acc += ws[i] * (zv[i][sl].mean() - zv[i][ns].mean())
        perm[t] = acc
    p = float((np.sum(np.abs(perm) >= abs(v)) + 1) / (nperm + 1))
    return pd.DataFrame(rows), dict(pooled=v, pooled_sig=s, z=v / s, perm_p=p, n_majors=len(gl))


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--test", type=int, default=0)
    a = ap.parse_args(); test = a.test > 0
    TB = Q.tbdir(test); lg = Q.log("zones_test" if test else "zones"); t0 = time.time()
    H, inj = R5.load_H(lg, limit=a.test)
    majors = Q.majors()
    zd = json.load(open(f"{R3.TB3}/zones_def.json")); nz = int(zd["run_zones"]); bounds = list(zd["bounds"][str(nz)])
    ax = json.load(open(f"{R3.TB3}/along_strike_axis.json")); edges = np.array(ax["edges"])
    B = RZ.bins_table()
    lg(f"zones: {nz}, bounds {np.round(bounds, 1).tolist()} km")
    LS3, rs3 = s3_labels(majors, bounds, lg)
    LF_, rf = fam_labels(majors, bounds, set(int(x) for x in H.fam.unique()), lg)
    pd.DataFrame(rs3 + rf).to_csv(f"{TB}/r7_zones_labels.csv", index=False)
    PM, PO, SL, SLR, INJ, CT, LO = [], [], [], [], [], [], []
    verdict = {}
    for setup in (Q.LSEL, Q.NOM):
        cs0, good = R5.setup_cells(H, setup)
        good = good & np.isfinite(inj["nat"])
        cs = cs0[good].reset_index(drop=True)
        Sh, ok = shares_for(cs, bounds, B, edges, lg)
        sel = np.flatnonzero(good)[ok]                      # rows of H used
        cs = cs[ok].reset_index(drop=True); Sh = Sh[ok]
        cs, _ = R4.identifiable(cs) if False else (cs, 0)
        L, F = fit(cs, Sh, nz)
        T, gz = per_major(L, F, majors, nz); T["setup"] = setup; PM.append(T)
        P, comp = pooled(F, gz, nz); P["setup"] = setup; PO.append(P)
        lg(f"[{setup}] cells {len(cs)}, majors {len(gz)}; pooled " + ", ".join(f"{k} {v:+.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in P.items()))
        for lname, labs in (("S3", LS3), ("F", LF_)):
            rows, pc = slip_contrast(F, gz, labs, np.random.default_rng(Q.SEED))
            if len(rows):
                rows["setup"] = setup; rows["labels"] = lname; SLR.append(rows)
            if pc:
                pc.update(setup=setup, labels=lname); SL.append(pc)
                lg(f"[{setup}] slipped − not ({lname}): {pc['pooled']:+.4f} ± {pc['pooled_sig']:.4f} % (z {pc['z']:+.2f}, permutation p {pc['perm_p']:.3f}, {pc['n_majors']} majors)")
        # control
        cc = R4.control_cells(cs)
        key = cs.series + "|" + cs.episode
        ShC = Sh[pd.Series(np.arange(len(cs)), index=key).reindex(cc.series + "|" + cc.episode).to_numpy().astype(int)]
        Lc, Fc = fit(cc, ShC, nz)
        Tc, gzc = per_major(Lc, Fc, majors, nz)
        con = [c for c in Tc.columns if "-" in c and not c.endswith("_sig")]
        fr = float(np.mean([abs(r[c]) < 2 * r[c + "_sig"] for _, r in Tc.iterrows() for c in con])) if len(Tc) else np.nan
        Pc, _ = pooled(Fc, gzc, nz) if gzc else ({}, None)
        pool_ok = bool(all(abs(Pc[c]) < 2 * Pc[c + "_sig"] for c in con)) if Pc else None
        ctrl = dict(setup=setup, cells=len(cc), majors=len(gzc), frac_within2=fr, pooled_within2=pool_ok,
                    passed=bool(fr >= 0.9 and pool_ok), **{f"pooled_{k}": v for k, v in Pc.items()})
        for lname, labs in (("S3", LS3), ("F", LF_)):
            _, pcc = slip_contrast(Fc, gzc, labs, np.random.default_rng(Q.SEED), nperm=200)
            if pcc:
                ctrl[f"slip_{lname}"] = pcc["pooled"]; ctrl[f"slip_{lname}_sig"] = pcc["pooled_sig"]
                ctrl[f"slip_{lname}_ok"] = bool(abs(pcc["pooled"]) < 2 * pcc["pooled_sig"])
        CT.append(ctrl); lg(f"[{setup}] control: " + ", ".join(f"{k} {v:+.4f}" if isinstance(v, float) else f"{k} {v}" for k, v in ctrl.items()))
        # injections
        rec, leak = [], []
        for e in gz:
            m = (cs.episode == e).to_numpy()
            for k in range(nz):
                x = np.where(m, cs.a_d.to_numpy() * Sh[:, k] * 0.3, 0.0)
                th = L.Op @ x; Fx = dict(F, th=th)
                vals = [float(g @ th) for g in gz[e]]
                rec.append(vals[k] / 0.3); leak.append(max(abs(vals[j]) for j in range(nz) if j != k) / 0.3)
        dn = inj["nat"][sel]
        pdl, psl = [], []
        for e in gz:
            for run, lst, amp in ((f"Pd|{setup}|{e}", pdl, 0.3), (f"Ps|{setup}|{e}", psl, 0.1)):
                if run not in inj:
                    continue
                v = inj[run][sel]; x = np.where(np.isfinite(v), v - dn, 0.0); th = L.Op @ x
                for i, j in itertools.combinations(range(nz), 2):
                    lst.append(abs(float((gz[e][i] - gz[e][j]) @ th)) / amp)
        I = dict(setup=setup, zone_pulse_rec_med=float(np.median(rec)), zone_pulse_rec_lo=float(np.min(rec)), zone_pulse_rec_hi=float(np.max(rec)),
                 cross_leak_med=float(np.median(leak)), cross_leak_max=float(np.max(leak)),
                 Pd_contrast_leak_med=float(np.median(pdl)) if pdl else np.nan, Ps_contrast_leak_med=float(np.median(psl)) if psl else np.nan,
                 Ps_contrast_leak_max=float(np.max(psl)) if psl else np.nan)
        I["passed"] = bool(0.8 <= I["zone_pulse_rec_med"] <= 1.1 and I["cross_leak_med"] <= 0.2 and (not psl or I["Ps_contrast_leak_med"] <= 0.35))
        INJ.append(I); lg(f"[{setup}] injections: " + ", ".join(f"{k} {v:.3f}" if isinstance(v, float) else f"{k} {v}" for k, v in I.items()))
        # leave-one-site-out (pooled contrasts)
        for s in sorted(cs.site.unique()):
            m = (cs.site != s).to_numpy()
            sub = cs[m].reset_index(drop=True)
            sub2, _ = R4.identifiable(sub)
            keep = (sub.series + "|" + sub.episode).isin(sub2.series + "|" + sub2.episode).to_numpy()
            if sub2.episode.nunique() < 3:
                continue
            try:
                Ls, Fs = fit(sub2, Sh[m][keep], nz)
                _, gzs = per_major(Ls, Fs, majors, nz)
                if not gzs:
                    continue
                Ps_, _ = pooled(Fs, gzs, nz)
            except Exception as ex:
                lg(f"loso {s}: {ex}"); continue
            r = dict(setup=setup, left_out=s, **Ps_)
            r["within2"] = bool(all(abs(Ps_[c] - P[c]) < 2 * Ps_[c + "_sig"] for c in P if "-" in c and not c.endswith("_sig")))
            LO.append(r)
        if setup == Q.LSEL:
            nmj = len(T)
            pmin = float(T.p.min()) if len(T) else np.nan
            con = [c for c in P if "-" in c and not c.endswith("_sig")]
            zmax = max(abs(P[c]) / P[c + "_sig"] for c in con)
            slp = [x for x in SL if x["setup"] == setup and x["labels"] == "S3"]
            verdict = dict(region=C.REGION, n_zones=nz, n_majors=nmj, min_major_p=pmin, bonferroni=0.01 / max(nmj, 1),
                           episode_differs=bool(len(T) and pmin < 0.01 / nmj and ctrl["passed"]),
                           episodes_below_bonf=";".join(T.major[T.p < 0.01 / nmj]) if len(T) else "",
                           pooled_p=P["p"], pooled_max_z=zmax,
                           pooled_differs=bool(P["p"] < 0.01 and zmax >= 3 and ctrl["passed"] and I["passed"]),
                           sum_chi2=float(T.chi2.sum()), sum_chi2_dof=int(nmj * (nz - 1)),
                           sum_chi2_p=float(stats.chi2.sf(T.chi2.sum(), nmj * (nz - 1))) if nmj else np.nan,
                           control_passed=ctrl["passed"], injection_passed=I["passed"],
                           pooled_contrasts={c: [P[c], P[c + "_sig"]] for c in con},
                           bound2_contrasts={c: abs(P[c]) + 2 * P[c + "_sig"] for c in con})
            if slp:
                s_ = slp[0]
                verdict.update(slip_S3=s_["pooled"], slip_S3_sig=s_["pooled_sig"], slip_S3_perm_p=s_["perm_p"], slip_S3_majors=s_["n_majors"],
                               slip_locked=bool(abs(s_["z"]) >= 3 and s_["perm_p"] < 0.01 and ctrl.get("slip_S3_ok", True)),
                               slip_bound2=abs(s_["pooled"]) + 2 * s_["pooled_sig"])
            slf = [x for x in SL if x["setup"] == setup and x["labels"] == "F"]
            if slf:
                verdict.update(slip_F=slf[0]["pooled"], slip_F_sig=slf[0]["pooled_sig"], slip_F_perm_p=slf[0]["perm_p"], slip_F_majors=slf[0]["n_majors"])
    pd.concat(PM, ignore_index=True).to_csv(f"{TB}/r7_zones_permajor.csv", index=False)
    pd.DataFrame(PO).to_csv(f"{TB}/r7_zones_pooled.csv", index=False)
    pd.DataFrame(SL).to_csv(f"{TB}/r7_zones_slip.csv", index=False)
    (pd.concat(SLR, ignore_index=True) if SLR else pd.DataFrame()).to_csv(f"{TB}/r7_zones_slip_permajor.csv", index=False)
    pd.DataFrame(INJ).to_csv(f"{TB}/r7_zones_inject.csv", index=False)
    pd.DataFrame(CT).to_csv(f"{TB}/r7_zones_ctrl.csv", index=False)
    pd.DataFrame(LO).to_csv(f"{TB}/r7_zones_loso.csv", index=False)
    verdict["written"] = Q.pt_now()
    Q.dump(verdict, f"{TB}/r7_zones_verdict.json")
    lg("VERDICT " + json.dumps(verdict, default=str))
    lg(f"done {time.time() - t0:.0f}s")
    if not test:
        open(f"{Q.R7}/ZONES_{C.REGION}_DONE", "w").write(Q.pt_now() + "\n")


if __name__ == "__main__":
    main()
