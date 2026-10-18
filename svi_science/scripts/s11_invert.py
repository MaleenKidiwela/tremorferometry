#!/usr/bin/env python
"""s11_invert -- S11 step 5, inversion D4 (rules: svi_science/results/s11_depth/PLAN.md).
d (dv/v, %) per cell (family, station, epoch, episode) = o_series + a_s S_site(e) + a_d D(e)
- offsets per (family, station, epoch) series free (eliminated exactly by weighted de-meaning = Schur complement);
- S per site (PGC = PGC/B009/B010/B011, B926, VGZ, SNB) per episode; D per episode (25-45 km layer, SVI);
- weights 1/sigma^2, sigma inflated per series by sqrt(max(1, chi2_red of its scrambled-control cells));
- Gaussian prior sigma 1 % on S and D (conditioning only); variant S = 2 annual harmonics + per-year step per site.
Reports formal and chi2-scaled sigma, resolution diagonal, low-eigenvalue directions (data-only sigma > 0.1 %), the prior
x2 / x0.5 test, and the projection of an episode-varying uniform 5-25 km (and > 45 km) change onto D.

    python s11_invert.py [--paths raw,svd] [--tag ""]    -> tables/inv_*.csv, tables/inv_summary.csv, INVERT_DONE
"""
from __future__ import annotations
import os, sys, argparse, json, time
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s11_common as C

PRIOR = 1.0          # % on S and D
LOWSIG = 0.1         # % data-only sigma along an eigen-direction = "below the noise level"
# Instrument vetoes (added 2026-10-14, S14 follow-up): rows of svi_science/instrument_vetoes.csv with decision == "veto"
# drop cells whose date_mid falls in [start, end] at that station. Default ON; S11_APPLY_VETO=0 reproduces the pre-veto
# numbers on disk (all S11 tables written before 2026-10-14 are pre-veto).
VETO_FILE = f"{C.ROOT}/svi_science/instrument_vetoes.csv"
APPLY_VETO = os.environ.get("S11_APPLY_VETO", "1") != "0"


def veto_mask(netsta, date_mid):
    """True where a cell (station, mid date) falls inside an instrument veto window"""
    m = np.zeros(len(netsta), bool)
    if not os.path.exists(VETO_FILE):
        return m
    Vt = pd.read_csv(VETO_FILE)
    Vt = Vt[Vt.decision.str.strip() == "veto"]
    d = pd.to_datetime(pd.Series(date_mid)).astype("datetime64[ns]").to_numpy()
    ns = np.asarray(netsta)
    for r in Vt.itertuples():
        m |= (ns == r.netsta) & (d >= np.datetime64(pd.Timestamp(r.start))) & (d <= np.datetime64(pd.Timestamp(r.end)))
    return m


# ------------------------------------------------------------------------------ data assembly
def kernels(variant="ell40", tag=""):
    K = pd.read_parquet(f"{C.TB}/kernel_coef{tag}.parquet")
    pe = f"{C.TB}/kernel_coef_ext.parquet"           # Revision 2: kernels of the added families (distinct series keys)
    if os.path.exists(pe) and not tag:
        K = pd.concat([K, pd.read_parquet(pe)], ignore_index=True)
    K = K[(K.variant == variant) & (K.comp == "3c")]
    return K


def attach_kernels(cells, K, win):
    k = K[K.win == win][["series", "fam", "netsta", "a_s", "a_m", "a_d", "a_b"]]
    cells = cells.copy(); cells["series"] = cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch
    out = cells.merge(k.drop(columns=["fam", "netsta"]), on="series", how="left")
    # series absent from the kernel table (no S10 grand stack): same (family, station), any epoch
    miss = out.a_d.isna()
    if miss.any():
        alt = k.groupby(["fam", "netsta"])[["a_s", "a_m", "a_d", "a_b"]].mean()
        fill = out.loc[miss, ["fam", "netsta"]].merge(alt, left_on=["fam", "netsta"], right_index=True, how="left")
        out.loc[miss, ["a_s", "a_m", "a_d", "a_b"]] = fill[["a_s", "a_m", "a_d", "a_b"]].to_numpy()
    return out


def load_cells(path, win, tag="", control=False, nmin=C.N_HEAD, ccmin=C.CC_HEAD, realign_rule=False, causal=1.5):
    """headline cells of one path/window -> DataFrame(fam, netsta, epoch, episode, site, d [%], sig [%], cc, n)"""
    E = pd.read_parquet(f"{C.TB}/episode_raw{tag}.parquet")
    ser = pd.read_csv(f"{C.TB}/series{tag}.csv")
    ser = ser[ser.n_data > 0]
    keep = ser.causal25 > causal
    if realign_rule:
        keep &= ser.realign_cc_med >= 0.5
    ok_series = set(zip(ser.fam[keep], ser.netsta[keep], ser.epoch[keep]))
    base = E[["fam", "netsta", "epoch", "episode", "n", "pre_rel"]].copy()
    base["instr_veto"] = veto_mask(E.netsta, E.date_mid) if "date_mid" in E else False
    if path == "raw":
        pre = "c2_" if control else ""
        V = E[["fam", "netsta", "epoch", "episode", f"{pre}eps3c_{win}", f"cc3c_{win}"]].copy()
        sigc = f"{pre}sig3c_{win}" if f"{pre}sig3c_{win}" in E else None
        V["sig"] = E[sigc] if sigc else np.nan
        if control and sigc is None:                  # control sigma only at the headline window
            V["sig"] = np.nan
        V = V.rename(columns={f"{pre}eps3c_{win}": "eps", f"cc3c_{win}": "cc"})
        if control:
            V["cc_ctrl"] = E[f"c2_cc3c_{win}"]
    else:
        f = f"{C.TB}/episode_svd{'_ctrl' if control else ''}{tag}.parquet"
        S_ = pd.read_parquet(f)
        V = S_[["fam", "netsta", "epoch", "episode", f"eps3c_{win}", f"cc3c_{win}", f"sig3c_{win}"]].rename(
            columns={f"eps3c_{win}": "eps", f"cc3c_{win}": "cc", f"sig3c_{win}": "sig"})
        if control:
            # data-cell cc decides which cells exist (same headline cells as the data path)
            Sd = pd.read_parquet(f"{C.TB}/episode_svd{tag}.parquet")[["fam", "netsta", "epoch", "episode", f"cc3c_{win}"]]
            V = V.drop(columns="cc").merge(Sd.rename(columns={f"cc3c_{win}": "cc"}), on=["fam", "netsta", "epoch", "episode"],
                                           how="inner")
    V = V.merge(base, on=["fam", "netsta", "epoch", "episode"], how="inner")
    V = V[[(f, n, e) in ok_series for f, n, e in zip(V.fam, V.netsta, V.epoch)]]
    V = V[(V.n >= nmin)]
    if APPLY_VETO:
        V = V[~V.instr_veto.astype(bool)]
    if path == "raw" and not control:
        V = V[V.cc >= ccmin]
    if path == "raw" and control:
        # same headline cells as the data path: data cc >= ccmin
        V = V[V.cc >= ccmin]
    V = V[np.isfinite(V.eps) & np.isfinite(V.sig) & (V.sig > 0)]
    V["d"] = -100 * V.eps; V["sig"] = 100 * V.sig
    V["site"] = V.netsta.map(C.SITE)
    return V.reset_index(drop=True)


def control_inflation(path, win, tag=""):
    """per series sqrt(max(1, chi2_red)) of scrambled-control cells about zero (station median if < 3 cells)"""
    Cc = load_cells(path, C.HEAD, tag, control=True, nmin=0, ccmin=-9)
    Cc = Cc[Cc.sig >= 0.01]                     # 2026-10-12 fix: a degenerate bootstrap sigma (1.8e-16 %) overflowed B010
    Cc["z2"] = (Cc.d / Cc.sig) ** 2
    g = Cc.groupby(["fam", "netsta", "epoch"]).z2.agg(["mean", "size"])
    gs = g[g["size"] >= 3].reset_index()
    st = gs.groupby("netsta")["mean"].median()
    st0 = Cc.groupby("netsta").z2.median()
    out = {}
    for (f, n, e), r in g.iterrows():
        chi = r["mean"] if r["size"] >= 3 else st.get(n, st0.get(n, 1.0))
        out[(f, n, e)] = float(np.sqrt(max(1.0, chi)))
    return out, st.to_dict()


# ------------------------------------------------------------------------------ core inversion
def design(cells, mode="episode"):
    """parameter index + dense design (without offsets)"""
    U = C.episodes()
    eps_ = sorted(cells.episode.unique(), key=lambda e: U.mid.get(e, pd.Timestamp("2100-01-01")))
    D_idx = {e: i for i, e in enumerate(eps_)}
    p = len(eps_)
    names = [("D", "", e) for e in eps_]
    if mode == "episode":
        se = sorted(set(zip(cells.site, cells.episode)), key=lambda t: (t[0], D_idx[t[1]]))
        S_idx = {t: p + i for i, t in enumerate(se)}
        names += [("S", s, e) for s, e in se]
        p += len(se)
    else:                                         # 2 annual harmonics + per-year step per site
        t = cells.episode.map(U.mid)
        cells = cells.assign(year=t.dt.year, ph=2 * np.pi * (t.dt.dayofyear - 1) / 365.25)
        S_idx = {}
        for s in sorted(cells.site.unique()):
            for h in ("c1", "s1", "c2", "s2"):
                S_idx[(s, h)] = p; names.append(("Sh", s, h)); p += 1
            for y in sorted(cells.year[cells.site == s].unique()):
                S_idx[(s, int(y))] = p; names.append(("Sy", s, int(y))); p += 1
    G = np.zeros((len(cells), p))
    r = np.arange(len(cells))
    G[r, [D_idx[e] for e in cells.episode]] = cells.a_d.to_numpy()
    if mode == "episode":
        G[r, [S_idx[(s, e)] for s, e in zip(cells.site, cells.episode)]] = cells.a_s.to_numpy()
    else:
        a_s = cells.a_s.to_numpy(); ph = cells.ph.to_numpy()
        for j, (s, y) in enumerate(zip(cells.site, cells.year)):
            G[j, S_idx[(s, "c1")]] = a_s[j] * np.cos(ph[j]); G[j, S_idx[(s, "s1")]] = a_s[j] * np.sin(ph[j])
            G[j, S_idx[(s, "c2")]] = a_s[j] * np.cos(2 * ph[j]); G[j, S_idx[(s, "s2")]] = a_s[j] * np.sin(2 * ph[j])
            G[j, S_idx[(s, int(y))]] = a_s[j]
    return G, names, D_idx


def demean(cells, G, d, w):
    key = cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch
    codes, uniq = pd.factorize(key)
    W = np.bincount(codes, w)
    Gm = np.zeros((len(uniq), G.shape[1])); np.add.at(Gm, codes, w[:, None] * G)
    dm = np.bincount(codes, w * d)
    return G - (Gm / W[:, None])[codes], d - (dm / W)[codes], len(uniq)


def invert(cells, mode="episode", prior=PRIOR, extras=True, d_override=None):
    """cells with d, sig (inflated), a_s, a_d, site, episode -> dict(params DataFrame, summary)"""
    G, names, D_idx = design(cells, mode)
    d = cells.d.to_numpy() if d_override is None else d_override
    w = 1.0 / cells.sig.to_numpy() ** 2
    Gt, dt, nser = demean(cells, G, d, w)
    Nd = Gt.T @ (Gt * w[:, None]); b = Gt.T @ (w * dt)
    P = np.eye(len(names)) / prior ** 2
    Ci = np.linalg.inv(Nd + P)
    th = Ci @ b
    res = dt - Gt @ th
    Rm = Ci @ Nd
    peff = float(np.trace(Rm))
    dof = max(len(d) - nser - peff, 1.0)
    chi2 = float(np.sum(w * res ** 2) / dof)
    sf = np.sqrt(np.diag(Ci)); ss = sf * np.sqrt(max(chi2, 1.0))
    par = pd.DataFrame(names, columns=["kind", "site", "episode"])
    par["value"] = th; par["sig_formal"] = sf; par["sig"] = ss; par["resolution"] = np.diag(Rm)
    cnt = cells.groupby("episode").size()
    ccw = cells.groupby("episode").apply(lambda q: np.average(q.cc, weights=1 / q.sig ** 2), include_groups=False)
    par["n_cells"] = [cnt.get(e, 0) if k == "D" else int(((cells.site == s) & (cells.episode == e)).sum()) if k == "S" else np.nan
                      for k, s, e in zip(par.kind, par.site, par.episode)]
    par["cc_w"] = [ccw.get(e, np.nan) if k in ("D", "S") else np.nan for k, e in zip(par.kind, par.episode)]
    U = C.episodes()
    par["mid"] = [U.mid.get(e, pd.NaT) if k in ("D", "S") else pd.NaT for k, e in zip(par.kind, par.episode)]
    summ = dict(mode=mode, n_cells=int(len(d)), n_series=int(nser), n_par=len(names), p_eff=peff, chi2red=chi2,
                n_D=int((par.kind == "D").sum()))
    if extras:
        lam, V = np.linalg.eigh(Nd)
        low = lam < 1.0 / LOWSIG ** 2
        isD = (par.kind == "D").to_numpy()
        summ["n_low_eig"] = int(low.sum())
        summ["frac_D_in_low"] = float(np.mean([np.sum(V[isD, i] ** 2) for i in np.flatnonzero(low)])) if low.any() else 0.0
        if mode == "episode" and low.any():
            # overlap of the low directions with D = -mean(S): per episode corr of v_D with mean_site v_S
            Si = par.index[par.kind == "S"]; ep_of = par.episode.to_numpy()
            cors = []
            for i in np.flatnonzero(low)[:20]:
                vD = pd.Series(V[isD, i], index=ep_of[isD]); vS = pd.Series(V[Si, i], index=ep_of[Si]).groupby(level=0).mean()
                j = pd.concat([vD, vS], axis=1).dropna()
                if len(j) > 3 and j.std().min() > 0:
                    cors.append(float(np.corrcoef(j.iloc[:, 0], j.iloc[:, 1])[0, 1]))
            summ["low_eig_corr_D_meanS_median"] = float(np.median(cors)) if cors else np.nan
        # prior test
        Dv = th[isD]
        for f in (2.0, 0.5):
            Pf = np.eye(len(names)) / (prior * f) ** 2
            thf = np.linalg.solve(Nd + Pf, b)
            summ[f"prior_x{f:g}_change"] = float(np.sqrt(np.mean((thf[isD] - Dv) ** 2)) / max(np.sqrt(np.mean(Dv ** 2)), 1e-12))
        summ["prior_pass"] = bool(max(summ["prior_x2_change"], summ["prior_x0.5_change"]) < 0.10)
        # projection of an episode-varying uniform 5-25 km (and > 45 km) change onto D
        rng = np.random.default_rng(C.SEED)
        eps_list = list(D_idx)
        for grp in ("a_m", "a_b"):
            sgn = rng.choice([-1.0, 1.0], len(eps_list)); dm = pd.Series(sgn, index=eps_list)
            ds = cells[grp].to_numpy() * cells.episode.map(dm).to_numpy()
            Gt2, dt2, _ = demean(cells, G, ds, w)
            th2 = Ci @ (Gt2.T @ (w * dt2))
            summ[f"proj_{grp[2:]}_onto_D"] = float(np.polyfit(dm.to_numpy(), th2[isD], 1)[0])
        summ["D_rms"] = float(np.sqrt(np.mean(Dv ** 2))); summ["D_sig_median"] = float(np.median(ss[isD]))
    return par, summ


def prepare(path, win, tag="", variant="ell40", nmin=C.N_HEAD, realign_rule=False, inflate=True, cells=None, ccmin=C.CC_HEAD):
    K = kernels(variant)
    if cells is None:
        cells = load_cells(path, win, tag, nmin=nmin, realign_rule=realign_rule, ccmin=ccmin)
    cells = attach_kernels(cells, K, win)
    cells = cells[np.isfinite(cells.a_d) & np.isfinite(cells.a_s)]
    if inflate:
        inf, _ = control_inflation(path, win, tag)
        cells["infl"] = [inf.get((f, n, e), 1.0) for f, n, e in zip(cells.fam, cells.netsta, cells.epoch)]
        cells["sig_raw"] = cells.sig; cells["sig"] = cells.sig * cells.infl
    # series with >= 3 cells
    key = cells.fam.astype(str) + "|" + cells.netsta + "|" + cells.epoch
    cells = cells[key.map(key.value_counts()) >= 3].reset_index(drop=True)
    return cells


RUNS = [  # label, win, variant, nmin, realign_rule, mode
    ("head", C.HEAD, "ell40", C.N_HEAD, False, "episode"),
    ("w2035", "2035", "ell40", C.N_HEAD, False, "episode"),
    ("w3560", "3560", "ell40", C.N_HEAD, False, "episode"),
    ("w36", "36", "ell40", C.N_HEAD, False, "episode"),
    ("n20", C.HEAD, "ell40", C.N_COL, False, "episode"),
    ("realign", C.HEAD, "ell40", C.N_HEAD, True, "episode"),
    ("ell20", C.HEAD, "ell20", C.N_HEAD, False, "episode"),
    ("ell80", C.HEAD, "ell80", C.N_HEAD, False, "episode"),
    ("two5_40", C.HEAD, "two5_40", C.N_HEAD, False, "episode"),
    ("seed2", C.HEAD, "ell40_s2", C.N_HEAD, False, "episode"),
    ("harm", C.HEAD, "ell40", C.N_HEAD, False, "harm"),
    # post hoc (added 2026-10-12 after the headline numbers): the raw cc3c >= 0.8 cell cut keeps ~20 % of n >= 50 cells
    ("cc0", C.HEAD, "ell40", C.N_HEAD, False, "episode", -9.0),
    ("cc0_w2035", "2035", "ell40", C.N_HEAD, False, "episode", -9.0),
    ("cc0_w3560", "3560", "ell40", C.N_HEAD, False, "episode", -9.0),
    ("cc0_harm", C.HEAD, "ell40", C.N_HEAD, False, "harm", -9.0),
]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--paths", default="raw,svd"); ap.add_argument("--tag", default="")
    ap.add_argument("--only", default=""); ap.add_argument("--out", default="")
    a = ap.parse_args()
    log = C.Log(f"{C.LOGD}/invert{a.tag}.log"); t0 = time.time()
    rows = []
    for path in a.paths.split(","):
        for run in RUNS:
            label, win, var, nmin, rr, mode = run[:6]
            ccmin = run[6] if len(run) > 6 else C.CC_HEAD
            if a.only and label not in a.only.split(","):
                continue
            try:
                cells = prepare(path, win, a.tag, var, nmin, rr, ccmin=ccmin)
                par, summ = invert(cells, mode)
            except Exception as e:
                import traceback
                log(f"ERROR {path} {label}: {e} {traceback.format_exc()[-500:]}"); continue
            par.to_csv(f"{C.TB}/inv_{path}_{label}{a.tag}.csv", index=False)
            cells.to_csv(f"{C.TB}/cells_{path}_{label}{a.tag}.csv", index=False)
            summ.update(path=path, label=label, win=win, variant=var, nmin=nmin, realign_rule=rr)
            rows.append(summ)
            log(f"{path} {label}: cells {summ['n_cells']} series {summ['n_series']} D {summ['n_D']} chi2 {summ['chi2red']:.2f} "
                f"D rms {summ.get('D_rms', np.nan):.4f} % sigD med {summ.get('D_sig_median', np.nan):.4f} % low-eig {summ.get('n_low_eig')} "
                f"prior x2 {summ.get('prior_x2_change', np.nan):.3f} x0.5 {summ.get('prior_x0.5_change', np.nan):.3f} "
                f"proj mid {summ.get('proj_m_onto_D', np.nan):.3f} below {summ.get('proj_b_onto_D', np.nan):.3f}")
    pd.DataFrame(rows).to_csv(f"{C.TB}/inv_summary{a.tag}{a.out}.csv", index=False)
    log(f"invert done {time.time() - t0:.0f}s")
    if not a.only:
        open(f"{C.OUT}/INVERT{a.tag}_DONE", "w").write(time.strftime("%Y-%m-%d %H:%M:%S\n"))


if __name__ == "__main__":
    main()
