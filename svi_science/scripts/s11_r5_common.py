"""s11_r5_common -- shared helpers for S11 round 5 (PLAN.md §9). Region via S11_REGION (SVI default, WA_NW).
Read-only on every earlier S11 table; outputs under s11_depth/round5/<svi|wa_nw>/. eps = -dv/v; cell d = dv/v in %.
Contents: paths; units (one shallow term per site-station, basis size by cell count); the general linear model Model5
(episode or smooth deep term; Ua / U1 unit terms; optional Dp projection; M, B; Sc per episode or smooth spline; ridge per
group or exact basis with null-space constraints); the round-4 injections re-aligned to H; the injection suite (T1-T5) and
the gates (G-ctrl, G-loso, G-fam, G-site, G-prior r5) for any model factory."""
from __future__ import annotations
import os, sys, json
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C
import s11_v2 as V
import s11_r3_common as R
import s11_r3_fit as LF
import s11_r4_common as R4

BASE = R4.BASE
SUB = R4.SUB
R5 = f"{BASE}/round5"
OUT5 = f"{R5}/{SUB}"
TB5, LOG5, FIG5 = f"{OUT5}/tables", f"{OUT5}/logs", f"{R5}/figures"
for _d in (TB5, LOG5, FIG5):
    os.makedirs(_d, exist_ok=True)
V6 = R4.V6
GEN = V6 + ["AVG6"]
U = C.episodes()
STERMS = R4.STERMS                       # c1 s1 c2 s2 lin h2007 h2011 h2015 h2019 h2023
SMALL = ["c1", "s1", "lin"]
RT = R4.RT                               # slow + c1 s1
PROJ = ["lin", "h2007", "h2011", "h2015", "h2019", "h2023", "c1", "s1"]
UMIN, UFULL, USMALL = 20, 40, 10         # PLAN §9.0 unit rules
PLAUS = {"SVI": ["two3_20", "two5_40", "ell40"], "WA_NW": ["two3_20", "two5_40"]}   # PLAN §9.5 (fixed)
HNOSC_R4 = {("SVI", "ell40"): 0.111, ("SVI", "two3_20"): 0.082, ("WA_NW", "ell40"): 0.128, ("WA_NW", "two3_20"): 0.124}


def log(name):
    return C.Log(f"{LOG5}/{name}.log")


def dump(obj, path):
    json.dump(obj, open(path, "w"), indent=1, default=lambda x: x.item() if hasattr(x, "item") else str(x))


def pt_now():
    return C.pt_now()


def fmt(r):
    return R4.fmt(r)


# ------------------------------------------------------------------------------ units
def units_of(cells, umin=UMIN):
    """-> (unit label per cell, {unit: basis terms}) — PLAN §9.0"""
    cnt = cells.netsta.value_counts()
    lab = np.array([ns if cnt[ns] >= umin else f"rest:{s}" for ns, s in zip(cells.netsta, cells.site)], dtype=object)
    n = pd.Series(lab).value_counts()
    terms = {}
    for u, k in n.items():
        terms[u] = STERMS if k >= UFULL else (SMALL if k >= USMALL else [])
    return lab, terms


# ------------------------------------------------------------------------------ smooth regional spline (task 2)
def spline_basis(mids_dec, doy, k_years):
    """cubic B-splines with interior knots every k years over the episode span, orthogonalised (unweighted, at the episode
    mid-dates) against the unit basis STERMS -> (n_ep, r) columns"""
    from scipy.interpolate import BSpline
    t0, t1 = float(mids_dec.min()) - 0.01, float(mids_dec.max()) + 0.01
    inner = np.arange(t0 + k_years, t1 - 0.5 * k_years + 1e-9, k_years)
    kn = np.r_[[t0] * 4, inner, [t1] * 4]
    nb = len(kn) - 4
    Bm = np.column_stack([BSpline(kn, np.eye(nb)[j], 3, extrapolate=False)(mids_dec) for j in range(nb)])
    Bm = np.nan_to_num(Bm)
    bb = V.basis(mids_dec, doy, STERMS)
    Um = np.column_stack([bb[t] for t in STERMS])
    Rr = Bm - Um @ np.linalg.lstsq(Um, Bm, rcond=None)[0]
    u, s, vt = np.linalg.svd(Rr, full_matrices=False)
    keep = s > 1e-6 * s.max()
    return u[:, keep], dict(n_spline=nb, n_kept=int(keep.sum()), knots=inner.round(2).tolist())


# ------------------------------------------------------------------------------ the model
def _demean(cells, G, w):
    key = cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch
    codes, uniq = pd.factorize(key)
    W = np.bincount(codes, w)
    Gm = np.zeros((len(uniq), G.shape[1])); np.add.at(Gm, codes, w[:, None] * G)
    return G - (Gm / W[:, None])[codes], codes, W, len(uniq)


class Model5:
    """d = o_series + a_d D + [unit terms] + [a_m M + a_b B] + [Sc] (+ round-4 S_site / R_sta for reproductions).
    kind 'episode': D(e) per episode (D0 sum-to-zero; Dp + orthogonal to PROJ at episode mids); 'smooth': D = a_d x RT basis.
    unit: None | 'as' | 'one'. Sc: None | 'episode' | ('spline', k). ridge: None (exact) | float | {group: sigma} on every
    non-episode group (groups: U, M, B, Scs, S, R, Dsm); ridge_D (smooth kind only): prior on the smooth D (round-4 form).
    r4_site / r4_rec: round-4 S per site (a_s x STERMS) / R per station (c1 s1 lin, coefficient 1, stations >= 20 cells)."""
    def __init__(self, cells, kind="episode", unit="as", Dproj=False, MB=True, Sc=None, ridge=None, ridge_D=None,
                 r4_site=False, r4_rec=False, apparent=False, reduce=True, pinv_tol=1e-10, umin=UMIN, fold_b=False, mid=None):
        c = cells.reset_index(drop=True); self.cells = c
        mid = c.episode.map(U.mid) if mid is None else pd.Series(np.asarray(mid)); t, doy = V.tdec(mid)
        n = len(c); r = np.arange(n)
        one = np.ones(n)
        ad, as_, am, ab = (one, one, one, one) if apparent else tuple(c[k].to_numpy(float) for k in ("a_d", "a_s", "a_m", "a_b"))
        if fold_b and not apparent:              # Amendment 1: the > 45 km layer joins the deep term (no B term)
            ad = ad + ab; MBb = False
        else:
            MBb = True
        eps_ = sorted(c.episode.unique(), key=lambda e: U.mid[e])
        blocks = []                              # (group, names, columns (n, k), reducible)
        info = {}
        if kind == "episode":
            Gd = np.zeros((n, len(eps_))); ie = {e: i for i, e in enumerate(eps_)}
            Gd[r, [ie[e] for e in c.episode]] = ad
            blocks.append(("D", [("D", "", e) for e in eps_], Gd, False))
        else:
            b = V.basis(t, doy, RT)
            blocks.append(("Dsm", [("D", "", tm) for tm in RT], np.column_stack([ad * b[tm] for tm in RT]), False))
        if Sc == "episode":
            Gs = np.zeros((n, len(eps_))); ie = {e: i for i, e in enumerate(eps_)}
            Gs[r, [ie[e] for e in c.episode]] = as_
            blocks.append(("Sc", [("Sc", "", e) for e in eps_], Gs, False))
        elif isinstance(Sc, tuple) and Sc[0] == "spline":
            em = U.mid.reindex(eps_); te, de = V.tdec(em.to_numpy())
            Bs, si = spline_basis(te, de, float(Sc[1])); info["spline"] = si
            ie = {e: i for i, e in enumerate(eps_)}; rows = np.array([ie[e] for e in c.episode])
            blocks.append(("Scs", [("Scs", "", j) for j in range(Bs.shape[1])], as_[:, None] * Bs[rows], True))
        if unit is not None:
            lab, terms = units_of(c, umin); self.unit_lab = lab; self.unit_terms = terms
            bS = V.basis(t, doy, STERMS)
            coef = as_ if unit == "as" else one
            for u_, tl in sorted(terms.items()):
                if not tl:
                    continue
                m = lab == u_
                G = np.zeros((n, len(tl)))
                for j, tm in enumerate(tl):
                    G[m, j] = coef[m] * bS[tm][m]
                blocks.append(("U", [("U", u_, tm) for tm in tl], G, True))
        if r4_site:
            bS = V.basis(t, doy, STERMS); site = c.site.to_numpy()
            for s in sorted(c.site.unique()):
                m = site == s; G = np.zeros((n, len(STERMS)))
                for j, tm in enumerate(STERMS):
                    G[m, j] = as_[m] * bS[tm][m]
                blocks.append(("S", [("S", s, tm) for tm in STERMS], G, False))
        if r4_rec:
            bR = V.basis(t, doy, ["c1", "s1", "lin"]); cnt = c.netsta.value_counts(); ns = c.netsta.to_numpy()
            for s in sorted([s for s in cnt.index if cnt[s] >= R4.RMIN]):
                m = ns == s; G = np.zeros((n, 3))
                for j, tm in enumerate(["c1", "s1", "lin"]):
                    G[m, j] = bR[tm][m]
                blocks.append(("R", [("R", s, tm) for tm in ["c1", "s1", "lin"]], G, False))
        if MB:
            bR = V.basis(t, doy, RT)
            blocks.append(("M", [("M", "", tm) for tm in RT], np.column_stack([am * bR[tm] for tm in RT]), True))
            if MBb:
                blocks.append(("B", [("B", "", tm) for tm in RT], np.column_stack([ab * bR[tm] for tm in RT]), True))
        self.w = 1.0 / c.sig.to_numpy() ** 2
        Graw = np.column_stack([b[2] for b in blocks])
        Gt, self.codes, self.Wg, self.nser = _demean(c, Graw, self.w)
        # exact-null reduction inside reducible blocks (does not change any fit)
        cols, names, groups = [], [], []
        o = 0; sw = np.sqrt(self.w)[:, None]; ndrop = 0
        for g, nm, G, red in blocks:
            k = G.shape[1]; B = Gt[:, o:o + k]; o += k
            if red and reduce and k > 0:
                u_, s_, vt = np.linalg.svd(B * sw, full_matrices=False)
                keep = s_ > 1e-6 * max(s_.max(), 1e-300)
                ndrop += int((~keep).sum())
                if not keep.any():
                    continue
                B = B @ vt[keep].T
                nm = [(nm[0][0], nm[0][1], f"v{j}") for j in range(int(keep.sum()))]
            cols.append(B); names += nm; groups += [g] * B.shape[1]
        self.Gt = np.column_stack(cols); self.names = names; self.groups = np.array(groups)
        self.idx = {nm: i for i, nm in enumerate(names)}
        info["n_null_dropped"] = ndrop
        self.N = self.Gt.T @ (self.Gt * self.w[:, None])
        P = len(names)
        prec = np.zeros(P)
        if ridge is not None:
            for g in set(groups):
                if g in ("D", "Sc"):
                    continue
                if g == "Dsm":
                    if ridge_D is not None:
                        prec[self.groups == g] = 1.0 / ridge_D ** 2
                    continue
                sp = ridge.get(g) if isinstance(ridge, dict) else ridge
                if sp is not None and np.isfinite(sp) and sp > 0:
                    prec[self.groups == g] = 1.0 / sp ** 2
        elif ridge_D is not None:
            prec[self.groups == "Dsm"] = 1.0 / ridge_D ** 2
        Mx = self.N + np.diag(prec)
        A = []
        if kind == "episode":
            iD = [self.idx[("D", "", e)] for e in eps_]
            a = np.zeros(P); a[iD] = 1.0; A.append(a)
            if Dproj:
                em = U.mid.reindex(eps_); te, de = V.tdec(em.to_numpy()); bb = V.basis(te, de, PROJ)
                for tm in PROJ:
                    if np.any(bb[tm] != 0):
                        a = np.zeros(P); a[iD] = bb[tm]; A.append(a)
            if Sc == "episode":
                a = np.zeros(P); a[[self.idx[("Sc", "", e)] for e in eps_]] = 1.0; A.append(a)
        self.singular = False
        if A:
            from scipy.linalg import null_space
            Z = null_space(np.array(A))
            Mz = Z.T @ Mx @ Z
        else:
            Z = None; Mz = Mx
        self.cond = float(np.linalg.cond(Mz))
        if not np.isfinite(self.cond) or self.cond > 1e12:
            self.singular = True; Mzi = np.linalg.pinv(Mz, rcond=pinv_tol, hermitian=True)
        else:
            Mzi = np.linalg.inv(Mz)
        self.Ci = Z @ Mzi @ Z.T if Z is not None else Mzi
        self.Op = self.Ci @ (self.Gt.T * self.w[None, :])
        self.peff = float(np.trace(self.Ci @ self.N))
        self.eps = eps_ if kind == "episode" else None
        self.extra = dict(eps=eps_ if kind == "episode" else None, U=U, ncell=c.groupby("episode").size(), zs=[""])
        self.info = info; self.kind = kind

    def demean_d(self, d):
        dm = np.bincount(self.codes, self.w * d) / self.Wg
        return d - dm[self.codes]

    def solve(self, d=None):
        d = self.cells.d.to_numpy() if d is None else d
        th = self.Op @ d
        r = self.demean_d(d) - self.Gt @ th
        chi2 = float(np.sum(self.w * r ** 2) / max(len(d) - self.nser - self.peff, 1))
        F = dict(th=th, C=self.Ci * max(chi2, 1.0), idx=self.idx, names=self.names, chi2=chi2, n=len(d), nser=self.nser, resid=r,
                 groups=self.groups)
        F.update(self.extra)
        return F

    def delta(self, x):
        F = dict(th=self.Op @ x, C=self.Ci, idx=self.idx, names=self.names, groups=self.groups); F.update(self.extra); return F


# ------------------------------------------------------------------------------ forms (PLAN §9.1/9.2)
FORMS = {
    "r4_HnoSc": dict(unit=None, r4_site=True, r4_rec=True, ridge=1.0, reduce=False),
    "r4_H": dict(unit=None, r4_site=True, r4_rec=True, ridge=1.0, reduce=False, Sc="episode"),
    "E1": dict(unit="as", Dproj=False), "E2": dict(unit="as", Dproj=True),
    "E3": dict(unit="one", Dproj=False), "E4": dict(unit="one", Dproj=True),
    "E2b": dict(unit="as", Dproj=True, fold_b=True), "E4b": dict(unit="one", Dproj=True, fold_b=True),   # Amendment 1
    "E2m": dict(unit="as", Dproj=True, fold_b=True, MB=False),                                         # column only
}
ELIGIBLE = ["E1", "E2", "E3", "E4", "E2b", "E4b"]
TIE = ["E2", "E2b", "E4", "E4b", "E1", "E3"]


def factory(form_kw, **over):
    kw = dict(form_kw); kw.update(over)
    return lambda cells, **o: Model5(cells, **{**kw, **o})


def q3(F, majors):
    return LF.q3_table(F, majors)


def summary(F, T, comp, label=""):
    return R4.summary(F, T, comp, label)


# ------------------------------------------------------------------------------ cells + round-4 injections
def load_H(lg, limit=0):
    import s11_r4_combined as CB
    H = CB.load_H(lg)
    z = np.load(f"{R4.TB4}/r4_injected.npz", allow_pickle=True)
    ser, ep = z["series"].astype(str), z["episode"].astype(str)
    if len(ser) != len(H) or not (np.all(ser == H.series.astype(str).to_numpy()) and np.all(ep == H.episode.astype(str).to_numpy())):
        raise RuntimeError("round-4 injections are not aligned with H")
    inj = {str(k): v for k, v in zip(z["runs"], z["d"])}
    if limit:                                    # code test only: first `limit` series
        keep = list(dict.fromkeys(H.series))[:limit]; m = H.series.isin(keep).to_numpy()
        H = H[m].reset_index(drop=True); inj = {k: v[m] for k, v in inj.items()}
    dn = inj["nat"]; ok = np.isfinite(dn)
    lg(f"H {len(H)} cells; round-4 injections {len(inj)} runs aligned; natural-from-store vs cells max |diff| "
       f"{np.nanmax(np.abs(dn[ok] - H.d.to_numpy()[ok])):.1e} %, missing {int((~ok).sum())}")
    return H, inj


def setup_cells(H, setup):
    cs = R4.with_setup(H, setup) if setup != "ell40" else H.copy()
    good = (np.isfinite(cs.a_d) & np.isfinite(cs.a_s) & np.isfinite(cs.a_m) & np.isfinite(cs.a_b)).to_numpy()
    return cs, np.array(good, bool)


# ------------------------------------------------------------------------------ suite (T1-T5 from the round-4 injections)
def resp(L, F0, x, majors):
    th = L.Op @ x
    gm = {e: LF.dvec(F0, e) for e in majors}
    gm = {e: g for e, g in gm.items() if g is not None}
    gc = LF.comp_vec(F0, majors)
    return {e: float(g @ th) for e, g in gm.items()}, (float(gc @ th) if gc is not None else np.nan)


def evaluate(H, inj, setup, majors, fac, sites_w, label, lg):
    cs, good = setup_cells(H, setup)
    good = good & np.isfinite(inj["nat"])
    cs = cs[good].reset_index(drop=True)
    L = fac(cs); F0 = L.delta(np.zeros(len(cs)))
    dn = inj["nat"][good]

    def X(run, mask=None):
        v = inj.get(run, inj["nat"])[good]
        return np.where(np.isfinite(v), v - dn, 0.0) * (1.0 if mask is None else mask)
    site = cs.site.to_numpy(); rows = []
    for g in GEN:
        rec, lk = [], []
        for e in majors:
            if f"Pd|{g}|{e}" not in inj:
                continue
            r_, _ = resp(L, F0, X(f"Pd|{g}|{e}"), [e])
            if e in r_:
                rec.append(r_[e] / 0.3)
            r_, _ = resp(L, F0, X(f"Ps|{g}|{e}"), [e])
            if e in r_:
                lk.append(abs(r_[e]) / 0.1)
        _, cps = resp(L, F0, X(f"Psall|{g}"), majors)
        rb, cb = resp(L, F0, X(f"B|{g}"), majors)
        ra, ca = resp(L, F0, X(f"A|{g}"), majors)
        rows.append(dict(label=label, setup=setup, gen=g, T1_Pd_rec=float(np.median(rec)) if rec else np.nan,
                         T1_lo=float(np.min(rec)) if rec else np.nan, T1_hi=float(np.max(rec)) if rec else np.nan,
                         T2_Ps_leak=float(np.median(lk)) if lk else np.nan, T2_hi=float(np.max(lk)) if lk else np.nan,
                         Psall_comp_leak=cps / 0.1, T5_B_comp=abs(cb) / 0.1,
                         B_major_med=float(np.median(np.abs(list(rb.values())))) / 0.1 if rb else np.nan, A_comp=ca / 0.3,
                         n_major=len(rec)))
    T = pd.DataFrame(rows)
    gin = setup if setup in GEN else "ell40"
    ons = []
    sites = sites_w.index[sites_w >= 0.05] if sites_w is not None else []
    for s in sites:
        m = (site == s).astype(float); lk = []
        for e in majors:
            if f"Ps|{gin}|{e}" not in inj:
                continue
            r_, _ = resp(L, F0, X(f"Ps|{gin}|{e}", m), [e])
            if e in r_:
                lk.append(abs(r_[e]) / 0.1)
        rr, cr = resp(L, F0, X("Rann", m), majors)
        ons.append(dict(label=label, setup=setup, site=s, weight_share=float(sites_w[s]), T3_onesite_leak=float(np.median(lk)) if lk else np.nan,
                        T3_hi=float(np.max(lk)) if lk else np.nan, T4_rann_comp=abs(cr) / 0.1,
                        T4_rann_major_med=float(np.median(np.abs(list(rr.values())))) / 0.1 if rr else np.nan))
    rr, cr = resp(L, F0, X("Rann"), majors)
    ons.append(dict(label=label, setup=setup, site="ALL", weight_share=1.0, T3_onesite_leak=np.nan, T4_rann_comp=abs(cr) / 0.1,
                    T4_rann_major_med=float(np.median(np.abs(list(rr.values())))) / 0.1 if rr else np.nan))
    O = pd.DataFrame(ons)
    ins = T[T.gen == gin].iloc[0]
    lg(f"  [{label}|{setup}] in-sample T1 {ins.T1_Pd_rec:.3f} T2 {ins.T2_Ps_leak:.3f} T5 {ins.T5_B_comp:.3f}; T1 over gen "
       f"{T.T1_Pd_rec.min():.2f}-{T.T1_Pd_rec.max():.2f}; T4 all {O[O.site == 'ALL'].T4_rann_comp.iloc[0]:.3f}/"
       f"{O[O.site == 'ALL'].T4_rann_major_med.iloc[0]:.3f}")
    return T, O


def rules(T, O, setup):
    gin = setup if setup in GEN else "ell40"
    ins = T[T.gen == gin].iloc[0]
    v = dict(T1=bool(0.8 <= ins.T1_Pd_rec <= 1.1), T1_val=float(ins.T1_Pd_rec), T2=bool(ins.T2_Ps_leak <= 0.35),
             T2_val=float(ins.T2_Ps_leak), T2_fixed=bool(ins.T2_Ps_leak <= 0.10), T5=bool(ins.T5_B_comp <= 0.35), T5_val=float(ins.T5_B_comp))
    o = O[O.site != "ALL"]
    v["T3"] = bool((o.T3_onesite_leak <= 0.35).all()) if len(o) else None
    v["T3_max"] = float(o.T3_onesite_leak.max()) if len(o) else np.nan
    v["T4"] = bool(((O.T4_rann_comp <= 0.35) & (O.T4_rann_major_med <= 0.35)).all())
    a = O[O.site == "ALL"].iloc[0]; v["T4_all_comp"] = float(a.T4_rann_comp); v["T4_all_major"] = float(a.T4_rann_major_med)
    v["T1_range_over_gen"] = [float(T.T1_Pd_rec.min()), float(T.T1_Pd_rec.max())]
    return v


# ------------------------------------------------------------------------------ gates
def _shift(Tm, comp, Tf, cf, mj):
    Tf = Tf.set_index("episode")
    e_ = [e for e in mj.episode if e in Tf.index]
    mv = float(np.max(np.abs((Tf.dD[e_] - Tm.dD[e_]) / Tm.sig[e_]))) if e_ else np.nan
    cv = abs(cf["composite"] - comp["composite"]) / comp["composite_sig"] if (cf and comp) else np.nan
    return mv, cv


def gates(H, setup, majors, fac, label, lg, exact=True, ridge_kw=None, ctrl_b1_fn=None, out=None):
    """G-site, G-prior (r5), G-ctrl (B2 + B1 column), G-loso, G-fam for one model factory and setup"""
    cs, good = setup_cells(H, setup); cs = cs[good].reset_index(drop=True)
    L = fac(cs); F = L.solve(); T, comp = q3(F, majors); Tm = T.set_index("episode"); mj = T[T.major]
    G = summary(F, T, comp, label); G.update(setup=setup, cond=L.cond, singular=L.singular, n_null_dropped=L.info.get("n_null_dropped", 0),
                                             n_params=len(L.names))
    if "spline" in L.info:
        G["spline_kept"] = L.info["spline"]["n_kept"]
    G["acceptance"] = bool(G["n_major_sig_le_0p1"] >= 5)
    w = 1 / cs.sig ** 2; sw = (w.groupby(cs.site).sum() / w.sum()).sort_values(ascending=False)
    G["top_site"], G["top_site_share"] = sw.index[0], float(sw.iloc[0]); G["G_site_pass"] = bool(sw.iloc[0] <= 0.5)
    pv = []
    if exact:
        for sp in (1.0, 2.0):
            Tf, cf = q3(fac(cs, ridge=sp).solve(), majors); mv, cv = _shift(Tm, comp, Tf, cf, mj)
            G[f"prior_ridge{sp:g}"] = mv; G[f"prior_ridge{sp:g}_comp"] = cv; pv += [mv, cv]
        for g in sorted(set(L.groups) - {"D", "Sc", "Dsm"}):          # Amendment 1: one group at a time (diagnostic)
            for sp in (1.0, 2.0):
                Tf, cf = q3(fac(cs, ridge={g: sp}).solve(), majors); mv, cv = _shift(Tm, comp, Tf, cf, mj)
                G[f"diag_ridge_{g}{sp:g}"] = max(mv, cv) if np.isfinite(cv) else mv
        if L.singular:
            for tol in (1e-8, 1e-12):
                Tf, cf = q3(fac(cs, pinv_tol=tol).solve(), majors); mv, cv = _shift(Tm, comp, Tf, cf, mj)
                G[f"prior_tol{tol:g}"] = mv; G[f"prior_tol{tol:g}_comp"] = cv; pv += [mv, cv]
    else:
        base = ridge_kw["ridge"]
        for f_ in (2.0, 0.5):
            rr = {k: v * f_ for k, v in base.items()} if isinstance(base, dict) else base * f_
            Tf, cf = q3(fac(cs, ridge=rr).solve(), majors); mv, cv = _shift(Tm, comp, Tf, cf, mj)
            G[f"prior_x{f_:g}"] = mv; G[f"prior_x{f_:g}_comp"] = cv; pv += [mv, cv]
    G["prior_max"] = float(np.nanmax(pv)) if pv else np.nan
    G["G_prior_pass"] = bool(G["prior_max"] < 0.1)
    for kind in ("g", "b1"):
        src = cs if kind == "g" else (ctrl_b1_fn(cs) if ctrl_b1_fn else None)
        if src is None:
            continue
        cc = R4.control_cells(src); cc, _ = R4.identifiable(cc)
        if len(cc) < 20 or cc.episode.nunique() < 3:
            G[f"ctrl_{kind}_pass"] = None; G[f"ctrl_{kind}_cells"] = len(cc); continue
        Tc, ccomp = q3(fac(cc).solve(), majors); Tc = Tc.set_index("episode")
        meas = [e for e in mj.episode if e in Tc.index and Tc.n_cells[e] >= 5]
        frac = float(np.mean([abs(Tc.dD[e]) < 2 * Tc.sig[e] for e in meas])) if meas else np.nan
        G.update({f"ctrl_{kind}_cells": len(cc), f"ctrl_{kind}_n": len(meas), f"ctrl_{kind}_frac": frac,
                  f"ctrl_{kind}_comp": ccomp["composite"] if ccomp else np.nan, f"ctrl_{kind}_comp_sig": ccomp["composite_sig"] if ccomp else np.nan,
                  f"ctrl_{kind}_pass": bool(frac >= 0.9 and (ccomp is None or abs(ccomp["composite"]) < 2 * ccomp["composite_sig"])) if meas else None})
        if out:
            Tc.reset_index().to_csv(f"{out}_ctrl_{kind}.csv", index=False)
    worst, ws_ = 1.0, ""; lo = []
    for s in sorted(cs.site.unique()):
        sub, _ = R4.identifiable(cs[cs.site != s])
        if sub.episode.nunique() < 3:
            continue
        Ts, _ = q3(fac(sub).solve(), majors); Ts = Ts.set_index("episode")
        ee = [e for e in mj.episode if e in Ts.index]
        if not ee:
            continue
        fr = float(np.mean([abs(Ts.dD[e] - Tm.dD[e]) < 2 * Ts.sig[e] for e in ee]))
        lo.append(Ts.reset_index().assign(left_out=s))
        if fr < worst:
            worst, ws_ = fr, s
    if lo and out:
        pd.concat(lo, ignore_index=True).to_csv(f"{out}_loso.csv", index=False)
    G.update(loso_worst=worst, loso_worst_site=ws_, G_loso_pass=bool(worst >= 0.9))
    oo, ee_ = (R4.identifiable(cs[cs.fam % 2 == k])[0] for k in (1, 0))
    co = ce = None
    if oo.episode.nunique() >= 3 and ee_.episode.nunique() >= 3:
        To, co = q3(fac(oo).solve(), majors); Te, ce = q3(fac(ee_).solve(), majors)
    if co and ce:
        J_ = To.set_index("episode").join(Te.set_index("episode"), lsuffix="_o", rsuffix="_e", how="inner"); J_ = J_[J_.major_o]
        fr = float(np.mean(np.abs(J_.dD_o - J_.dD_e) < 2 * np.hypot(J_.sig_o, J_.sig_e))) if len(J_) else np.nan
        G.update(famsplit_frac=fr, famsplit_odd=co["composite"], famsplit_odd_sig=co["composite_sig"], famsplit_even=ce["composite"],
                 famsplit_even_sig=ce["composite_sig"],
                 G_fam_pass=bool(fr >= 0.9 and abs(co["composite"] - ce["composite"]) < 2 * np.hypot(co["composite_sig"], ce["composite_sig"])))
    if out:
        T.to_csv(f"{out}_dD.csv", index=False)
    lg(f"gates [{label}|{setup}]: " + fmt({k: v for k, v in G.items() if k not in ("label",)}))
    return G, T, sw, F, L, cs


def coef_rms(F, groups=("U", "M", "B", "Scs")):
    out = {}
    for g in groups:
        m = F["groups"] == g
        if m.any():
            out[g] = float(np.sqrt(np.mean(F["th"][m] ** 2)))
    return out
