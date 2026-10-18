#!/usr/bin/env python
"""SH_KMRFIX task 3 (2026-10-01, CPU): confidence model v2 with NC.KMR removed from EVERY feature.

Why: from 2025-09-17 the forward catalogue is the 16-station run mf_NK (KMR records no ground motion after its swap),
so the model must see 16-station inputs in training too.  Raw continuous data are deleted, so the no-KMR features of
the training rows are rebuilt from saved per-candidate values (results/sh_conf/candidates_conf.parquet + the per
family-day station noise tables results/sh_conf/feat/<run>/<date>_noise.parquet):
  network  ccsum' = ccsum - cc_KMR * n_KMR (channel-mean cc of KMR at the candidate x its live channels: 3, 2 on
           2025-08-03..25 when HHE was masked dead); nlive' = nlive - n_KMR; meancc' = ccsum'/nlive'; nchan', nsta';
           net_mad' = (ccsum' - med') / mad' with med = ccsum - net_mad * MAD (MAD of the day's network sum, stored),
           med' = med - n_KMR * median(KMR), sigma'^2 = (MAD/0.6745)^2 - (n_KMR * sigma_KMR)^2 (sigma_KMR = spread of
           KMR's channel-mean cc at 1000 random origins of that family-day).  Checked against the real 16-station run
           on 2025-09-18 .. 2026-09-28 (B/C rows vs mf_NK at the same peaks).
  site     every site feature recomputed from the per-station cc_* / z_* columns without KMR (KBU + B049 = one site:
           live-channel-weighted cc; z from the two stations' noise, independent-noise approximation).
  PNET     per-site max of the per-station PNET columns without KMR; pnet_rel_max shifted by the change of pnet_max.
  HF       hfr_rel set missing where the family's HF station was NC.KMR.
  context  S1 flags (pass_9x, pass_mcc, pass_floor, 4-s winner, in_S1_9x) recomputed from the v2 network values;
           densities unchanged except for the mf_NK rows (recomputed from the mf_NK >= 8x candidates).
Forward rows 2025-09-17..: the mf_NK >= 8x candidates (candidates8x_v2.parquet) with network values from mf_NK and
site/PNET/HF values from the cut-window candidate at the SAME peak (same family, |dt| <= 0.1 s); others use a fallback
model without site/PNET/HF features (FALLBACK_v2, fitted the same way).
Model, folds, fit days, labels, evaluation sets: exactly as sh_conf_fit.py (HGB + isotonic, stored 5 day-folds).
Outputs results/sh_kmrfix/{conf_v2_rows.parquet, t3_*.csv, t3_*.json}, results/cut/sh_ext/conf_P_v2_by_det_id.parquet."""
import os, sys, json, time, glob, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings("ignore")
import sh_conf_common as C
import sh_conf_fit as FT
from sh_v3_common import apply_rules, decluster
RES = C.RES; OUT = f"{RES}/sh_kmrfix"; CUT = f"{RES}/cut/sh_ext"; EXT = C.EXT
KMR = "NC.KMR"; STAS = [s for s in C.STAS if s != KMR]
NCOMP = {s: (1 if s in ("NC.KPR", "NC.KBU") else 3) for s in C.STAS}
SITE = C.sites()
SW1 = pd.Timestamp("2025-09-17")
MAP_DT = 0.1


def log(s): print(f"[{time.strftime('%H:%M:%S')}] {s}", flush=True)


def noise_table():
    p = f"{OUT}/t3_noise_pivot.parquet"
    if os.path.exists(p): return pd.read_parquet(p)
    fs = sorted(glob.glob(f"{RES}/sh_conf/feat/*/*_noise.parquet"))
    N = pd.concat([pd.read_parquet(f, columns=["run", "date", "fam", "netsta", "noise_med", "noise_mad"]) for f in fs], ignore_index=True)
    N["fam"] = N.fam.astype(str); N["pres"] = 1.0
    P = N.pivot_table(index=["run", "date", "fam"], columns="netsta", values=["pres", "noise_med", "noise_mad"], aggfunc="first")
    P.columns = [f"{a}_{b}" for a, b in P.columns]; P = P.reset_index()
    P.to_parquet(p, index=False); return P


# ------------------------------------------------------------------------------------------------ site / PNET features
def site_feats(R, stations):
    groups = {}
    for s in stations: groups.setdefault(SITE[s], []).append(s)
    n = len(R); K = len(groups)
    SC = np.full((n, K), np.nan); Z = np.full((n, K), np.nan); MD = np.full((n, K), np.nan); PR = np.zeros((n, K), bool); SS = np.zeros((n, K))
    g = lambda c: R[c].to_numpy(float) if c in R else np.full(n, np.nan)
    for k, (sid, sts) in enumerate(sorted(groups.items())):
        if len(sts) == 1:
            s = sts[0]; c = g(f"cc_{s}"); pr = np.nan_to_num(g(f"pres_{s}")) > 0
            SC[:, k] = c; Z[:, k] = g(f"z_{s}"); MD[:, k] = g(f"noise_mad_{s}"); PR[:, k] = pr; SS[:, k] = np.nan_to_num(c) * NCOMP[s]
        else:
            ss = np.zeros(n); sl = np.zeros(n); wm = np.zeros(n); wv = np.zeros(n); wp = np.zeros(n); pr = np.zeros(n, bool)
            for s in sts:
                c = g(f"cc_{s}"); w = NCOMP[s]; p = np.nan_to_num(g(f"pres_{s}")) > 0
                ss += w * np.nan_to_num(c); sl += w * np.isfinite(c)
                wm += np.where(p, w * np.nan_to_num(g(f"noise_med_{s}")), 0); wv += np.where(p, (w * np.nan_to_num(g(f"noise_mad_{s}"))) ** 2, 0)
                wp += np.where(p, w, 0); pr |= p
            with np.errstate(invalid="ignore", divide="ignore"):
                SC[:, k] = np.where(sl > 0, ss / sl, np.nan)
                medc = np.where(wp > 0, wm / wp, np.nan); madc = np.where(wp > 0, np.sqrt(wv) / wp, np.nan)
                Z[:, k] = (SC[:, k] - medc) / np.where(madc > 0, madc, np.nan)
            MD[:, k] = madc; PR[:, k] = pr; SS[:, k] = ss
    SC[~PR] = np.nan; Z[~PR] = np.nan; MD[~PR] = np.nan; SS[~PR] = 0.0
    Cs = np.sort(np.nan_to_num(SC, nan=-99), 1)[:, ::-1]; Zs = np.sort(np.nan_to_num(Z, nan=-99), 1)[:, ::-1]
    top = lambda A, k: np.where(A[:, k] <= -98, np.nan, A[:, k])
    nsite = PR.sum(1)
    with np.errstate(invalid="ignore", divide="ignore", all="ignore"):
        tot = SS.sum(1); share = np.where(PR, SS / np.where(tot > 0, tot, np.nan)[:, None], np.nan)
        out = dict(n_sites=nsite.astype(float), n_sites_live=np.isfinite(SC).sum(1).astype(float),
                   site_cc_max=np.where(nsite >= 1, top(Cs, 0), np.nan), site_cc_2nd=np.where(nsite >= 2, top(Cs, 1), np.nan),
                   site_z_max=np.where(nsite >= 1, top(Zs, 0), np.nan), site_z_2nd=np.where(nsite >= 2, top(Zs, 1), np.nan),
                   site_z_min=np.nanmin(np.where(np.isfinite(Z), Z, np.nan), 1),
                   n_z2=(np.nan_to_num(Z, nan=-9) >= 2).sum(1).astype(float), n_z3=(np.nan_to_num(Z, nan=-9) >= 3).sum(1).astype(float),
                   n_z5=(np.nan_to_num(Z, nan=-9) >= 5).sum(1).astype(float), dom_top=np.nanmax(share, 1),
                   site_noise_mad_med=np.nanmedian(MD, 1))
    out["frac_z3"] = out["n_z3"] / np.clip(out["n_sites"], 1, None)
    return out


def pnet_feats(R, stations, sfx):
    groups = {}
    for s in stations: groups.setdefault(SITE[s], []).append(s)
    n = len(R); V = np.full((n, len(groups)), -1.0)
    for k, (sid, sts) in enumerate(sorted(groups.items())):
        for s in sts:
            c = f"pn{sfx}_{s}"
            if c in R: V[:, k] = np.fmax(V[:, k], np.nan_to_num(R[c].to_numpy(float), nan=-1.0))
    V[V < 0] = np.nan
    Vs = np.sort(np.nan_to_num(V, nan=-1.0), 1)[:, ::-1]
    nfin = np.isfinite(V).sum(1)
    return dict(pnet_max=np.where(Vs[:, 0] >= 0, Vs[:, 0], np.nan), pnet_2nd=np.where((Vs[:, 1] >= 0) & (nfin >= 2), Vs[:, 1], np.nan),
                n_sites_pnet03=(np.nan_to_num(V, nan=0) > 0.3).sum(1).astype(float), n_sites_pnet05=(np.nan_to_num(V, nan=0) > 0.5).sum(1).astype(float))


def recon_check(A):
    """reproduce the stored (with-KMR) site / PNET features from the per-station columns."""
    rows = []
    sf = site_feats(A, C.STAS)
    for c in ["n_sites", "n_sites_live", "site_cc_max", "site_cc_2nd", "site_z_max", "site_z_2nd", "site_z_min", "n_z3", "dom_top", "site_noise_mad_med"]:
        a = A[c].to_numpy(float); b = sf[c]; m = np.isfinite(a) & np.isfinite(b)
        rows.append(dict(feature=c, n=int(m.sum()), frac_both_nan_or_finite=float((np.isfinite(a) == np.isfinite(b)).mean()),
                         frac_within_1e3=float((np.abs(a[m] - b[m]) <= 1e-3 * np.maximum(1, np.abs(a[m]))).mean()), r=float(np.corrcoef(a[m], b[m])[0, 1])))
    for sfx in ("w", "f"):
        pf = pnet_feats(A, C.STAS, sfx)
        for c in ["pnet_max", "pnet_2nd", "n_sites_pnet03", "n_sites_pnet05"]:
            a = A[f"{c}_{sfx}"].to_numpy(float); b = pf[c]; m = np.isfinite(a) & np.isfinite(b)
            rows.append(dict(feature=f"{c}_{sfx}", n=int(m.sum()), frac_both_nan_or_finite=float((np.isfinite(a) == np.isfinite(b)).mean()),
                             frac_within_1e3=float((np.abs(a[m] - b[m]) <= 1e-3 * np.maximum(1, np.abs(a[m]))).mean()), r=float(np.corrcoef(a[m], b[m])[0, 1])))
    # network: sum over stations of cc * ncomp vs stored ccsum
    ss = np.zeros(len(A))
    for s in C.STAS: ss += np.nan_to_num(A[f"cc_{s}"].to_numpy(float)) * NCOMP[s]
    m = A.series.eq("data").to_numpy()
    rows.append(dict(feature="ccsum_from_stations", n=int(m.sum()), frac_within_1e3=float((np.abs(ss[m] - A.ccsum.to_numpy()[m]) <= 1e-3 * np.maximum(1, np.abs(A.ccsum.to_numpy()[m]))).mean()),
                     r=float(np.corrcoef(ss[m], A.ccsum.to_numpy()[m])[0, 1])))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ network approximation
def kmr_nlive_day(run, date):
    d = np.asarray(date)
    return np.where((np.asarray(run) == "fwd") & (d >= "2025-08-03") & (d <= "2025-08-25"), 2, 3)


def net_approx(R):
    """no-KMR network features for rows whose network sum included KMR."""
    n = len(R)
    pres = np.nan_to_num(R[f"pres_{KMR}"].to_numpy(float)) > 0 if f"pres_{KMR}" in R else np.zeros(n, bool)
    if "seg" in R: pres &= ~R.seg.isin(["NK"]).to_numpy()
    ccK = R[f"cc_{KMR}"].to_numpy(float)
    nd = kmr_nlive_day(R.run.to_numpy(), R.date.to_numpy()).astype(float)
    nl = np.where(pres & np.isfinite(ccK), nd, 0.0)
    mad = np.where(R.series.eq("data").to_numpy(), R.fam_mad.to_numpy(float), R.fam_mad_c.to_numpy(float))
    med = R.ccsum.to_numpy(float) - R.net_mad.to_numpy(float) * mad
    nmed = np.nan_to_num(R[f"noise_med_{KMR}"].to_numpy(float)) if f"noise_med_{KMR}" in R else np.zeros(n)
    nmad = np.nan_to_num(R[f"noise_mad_{KMR}"].to_numpy(float)) if f"noise_mad_{KMR}" in R else np.zeros(n)
    ndp = np.where(pres, nd, 0.0)
    cs = R.ccsum.to_numpy(float) - np.where(nl > 0, ccK, 0.0) * nl
    nlive = R.nlive.to_numpy(float) - nl
    sig2 = np.sqrt(np.clip((mad / 0.6745) ** 2 - (ndp * nmad) ** 2, ((mad / 0.6745) * 0.5) ** 2, None))
    mad2 = 0.6745 * sig2
    med2 = med - ndp * nmed
    return dict(ccsum=cs, nlive=nlive, meancc=cs / np.clip(nlive, 1, None), net_mad=(cs - med2) / mad2,
                nchan=R.nchan.to_numpy(float) - ndp, nsta=R.nsta.to_numpy(float) - pres.astype(float), mad_v2=mad2, kmr_in_net=pres)


def v2_transform(R, exact_net=False):
    out = {}
    if exact_net:
        for c in ["ccsum", "nlive", "meancc", "net_mad", "nchan", "nsta"]: out[c] = R[c].to_numpy(float)
    else:
        na = net_approx(R)
        for c in ["ccsum", "nlive", "meancc", "net_mad", "nchan", "nsta"]: out[c] = na[c]
        out["mad_v2"] = na["mad_v2"]; out["kmr_in_net"] = na["kmr_in_net"]
    out["live_frac"] = out["nlive"] / np.clip(out["nchan"], 1, None)
    out.update(site_feats(R, STAS))
    for sfx in ("w", "f"):
        pf = pnet_feats(R, STAS, sfx)
        old = R[f"pnet_max_{sfx}"].to_numpy(float); rel = R[f"pnet_rel_max_{sfx}"].to_numpy(float)
        shift = np.where(np.isfinite(old) & np.isfinite(pf["pnet_max"]), old - pf["pnet_max"], 0.0)
        for c, v in pf.items(): out[f"{c}_{sfx}"] = v
        out[f"pnet_rel_max_{sfx}"] = np.where(np.isfinite(pf["pnet_max"]), rel - shift, np.nan)
    hs = R.hf_sta.fillna("").to_numpy() if "hf_sta" in R else np.array([""] * len(R))
    out["hfr_rel"] = np.where(hs == KMR, np.nan, R.hfr_rel.to_numpy(float))
    return pd.DataFrame(out, index=R.index)


def context_v2(F):
    """S1 flags from v2 network values per (run, series); densities recomputed for rows flagged dens_recompute."""
    F["pass_9x"] = F.net_mad >= 9.0; F["pass_mcc"] = F.meancc >= 0.12; F["pass_floor"] = F.ccsum >= 2.75
    F["win4s_raw"] = False; F["in_S1_9x"] = False
    for (run, ser), g in F.groupby(["run", "series"]):
        D = g[["fam", "ccsum", "net_mad", "meancc"]].copy(); D["t"] = pd.to_datetime(g.t_epoch, unit="s").astype("datetime64[ns]")
        D = D.reset_index()
        F.loc[D["index"].to_numpy(), "win4s_raw"] = decluster(D.drop(columns="index"))
        E = apply_rules(D, 9); F.loc[E["index"].to_numpy(), "in_S1_9x"] = True
    G = C.fam_groups()
    sel = F.index[F.dens_recompute]
    if len(sel):
        S = F.loc[sel]
        for d, gi in S.groupby("date").groups.items():
            g = F.loc[gi]; fam = g.fam.to_numpy(); te = g.t_epoch.to_numpy()
            TT = {f: np.sort(te[fam == f]) for f in np.unique(fam)}; allt = np.sort(te)
            for f in TT:
                pos = gi[fam == f]; tq = F.loc[pos, "t_epoch"].to_numpy()
                cw = FT.count_within
                F.loc[pos, "dens_own_1h"] = cw(tq, TT[f], 3600) - 1; F.loc[pos, "dens_own_10m"] = cw(tq, TT[f], 600) - 1
                F.loc[pos, "dens_grp_1h"] = sum(cw(tq, TT[x], 3600) for x in G[f] if x in TT) - 1
                F.loc[pos, "dens_grp_10m"] = sum(cw(tq, TT[x], 600) for x in G[f] if x in TT) - 1
                F.loc[pos, "n_oth_4s"] = cw(tq, allt, 4.0 - 1e-9) - cw(tq, TT[f], 4.0 - 1e-9)
    for c in FT.S1: F[c] = F[c].astype(bool)
    return F


# ------------------------------------------------------------------------------------------------ main
def main():
    t00 = time.time()
    A = pd.read_parquet(f"{RES}/sh_conf/candidates_conf.parquet"); A["fam"] = A.fam.astype(str)
    NP = noise_table()
    A = A.merge(NP, on=["run", "date", "fam"], how="left")
    log(f"loaded {len(A)} rows; noise pivot {len(NP)}")
    # (0) reconstruction check
    rc = recon_check(A); rc.to_csv(f"{OUT}/t3_reconstruction_check.csv", index=False); log("\n" + rc.round(4).to_string())
    # det_id for forward data rows
    Dc = pd.read_parquet(f"{CUT}/detections.parquet", columns=["det_id", "fam", "t_epoch"]); Dc["fam"] = Dc.fam.astype(str)
    Dc["key"] = np.round(Dc.t_epoch * 20).astype(np.int64)
    A["key"] = np.round(A.t_epoch * 20).astype(np.int64)
    A = A.merge(Dc[["fam", "key", "det_id"]], on=["fam", "key"], how="left")
    A.loc[~((A.run == "fwd") & (A.series == "data")), "det_id"] = np.nan
    log(f"det_id joined {int(A.det_id.notna().sum())} of {int(((A.run == 'fwd') & (A.series == 'data')).sum())} forward data rows")
    # (1) v2 features for every stored row
    V = v2_transform(A, exact_net=False)
    NK_sw = A.seg.eq("NK").to_numpy()
    if NK_sw.any():                 # switch-day rows are already the no-KMR run: keep their network values
        for c in ["ccsum", "nlive", "meancc", "net_mad", "nchan", "nsta"]:
            V.loc[NK_sw, c] = A.loc[NK_sw, c].to_numpy(float)
        V.loc[NK_sw, "live_frac"] = V.loc[NK_sw, "nlive"] / V.loc[NK_sw, "nchan"].clip(lower=1)
    keep = ["series", "i", "t_epoch", "run", "date", "fam", "seg", "role", "build_any", "year", "fam_code", "build_own", "fold",
            "m_same", "m_any", "det_id", "dens_own_1h", "dens_own_10m", "dens_grp_1h", "dens_grp_10m", "n_oth_4s",
            "P_same_FULL", "P_any_FULL", "in_S1_9x"]
    B = pd.concat([A[keep].rename(columns={"in_S1_9x": "in_S1_9x_v1feat"}), V], axis=1)
    B["dens_recompute"] = False; B["feat_src"] = np.where(NK_sw, "NK_switchday_exact", "approx_noKMR")
    # (2) approximation check on 2025-09-18..: B/C rows (with-KMR network) vs the real 16-station run at the same peaks
    NKr = pd.read_parquet(f"{EXT}/mf_NK.parquet"); NKr["fam"] = NKr.fam.astype(str); NKr["key"] = np.round(NKr.t_epoch * 20).astype(np.int64)
    bc = A[A.seg.isin(["B", "C"]) & A.series.eq("data")].index
    J = B.loc[bc, ["fam", "date", "t_epoch", "net_mad", "meancc", "ccsum", "mad_v2"]].assign(key=A.loc[bc, "key"].to_numpy(),
            net_mad_withKMR=A.loc[bc, "net_mad"].to_numpy(), meancc_withKMR=A.loc[bc, "meancc"].to_numpy(), mad_withKMR=A.loc[bc, "fam_mad"].to_numpy())
    J = J.merge(NKr[["fam", "key", "net_mad", "meancc", "ccsum"]].rename(columns={"net_mad": "net_mad_NK", "meancc": "meancc_NK", "ccsum": "ccsum_NK"}), on=["fam", "key"], how="inner")
    GN = pd.read_csv(f"{EXT}/mf_NK_diag.csv"); GN["fam"] = GN.fam.astype(str); GN["date"] = GN.day.astype(str)
    J = J.merge(GN[["fam", "date", "mad"]].rename(columns={"mad": "mad_NK"}), on=["fam", "date"], how="left")
    ac = {}
    for a, b, name in [("net_mad", "net_mad_NK", "net_mad"), ("net_mad_withKMR", "net_mad_NK", "net_mad_uncorrected"), ("meancc", "meancc_NK", "meancc"),
                       ("meancc_withKMR", "meancc_NK", "meancc_uncorrected"), ("ccsum", "ccsum_NK", "ccsum"), ("mad_v2", "mad_NK", "mad"), ("mad_withKMR", "mad_NK", "mad_uncorrected")]:
        x = J[a].to_numpy(float); y = J[b].to_numpy(float); m = np.isfinite(x) & np.isfinite(y)
        ac[name] = dict(n=int(m.sum()), median_ratio=float(np.median(x[m] / y[m])), p10_ratio=float(np.percentile(x[m] / y[m], 10)),
                        p90_ratio=float(np.percentile(x[m] / y[m], 90)), r=float(np.corrcoef(x[m], y[m])[0, 1]))
    xs = J.net_mad.to_numpy(); ys = J.net_mad_NK.to_numpy()
    ac["net_mad_cross_9x"] = dict(both=int(((xs >= 9) & (ys >= 9)).sum()), approx_only=int(((xs >= 9) & (ys < 9)).sum()), NK_only=int(((xs < 9) & (ys >= 9)).sum()))
    json.dump(ac, open(f"{OUT}/t3_network_approx_check.json", "w"), indent=1); log("approx check " + json.dumps(ac))
    # (3) the forward 16-station rows (mf_NK >= 8x, 2025-09-17 ..)
    R = pd.read_parquet(f"{OUT}/candidates8x_v2.parquet"); R = R[R.seg == "NK"].copy()
    R["date"] = R.t.dt.strftime("%Y-%m-%d"); R["fam"] = R.fam.astype(str)
    useF = R.map_type.isin(["exact", "same_fam"]) & (R.map_dt_s.abs() <= MAP_DT)
    Arow = A[(A.run == "fwd") & (A.series == "data") & A.det_id.notna()].set_index("det_id")
    src = Arow.reindex(R.det_id.where(useF, -1).astype(float).to_numpy())
    src.index = R.index
    for c in ["ccsum", "nlive", "meancc", "net_mad", "nchan", "nsta"]: src[c] = R[c].to_numpy(float)
    src["run"] = "fwd"; src["series"] = "data"; src["date"] = R.date.to_numpy(); src["fam"] = R.fam.to_numpy(); src["seg"] = "NK"
    VN = v2_transform(src, exact_net=True)
    netc = ["ccsum", "nlive", "meancc", "net_mad", "nchan", "nsta", "live_frac"]
    VN.loc[~useF.to_numpy(), [c for c in VN.columns if c not in netc]] = np.nan
    N = pd.DataFrame(dict(series="data", i=-1, t_epoch=R.t_epoch.to_numpy(), run="fwd", date=R.date.to_numpy(), fam=R.fam.to_numpy(), seg="NK",
                          role="fwd_ext", build_any=False, year=R.t.dt.year.to_numpy(), fam_code=R.fam.astype(int).to_numpy() - 1, build_own=False,
                          fold=-1, m_same=np.nan, m_any=np.nan, det_id=R.det_id.where(R.map_type != "none", np.nan).to_numpy(float),
                          det_id_v2=R.det_id_v2.to_numpy(), P_same_FULL=np.nan, P_any_FULL=np.nan, in_S1_9x_v1feat=False), index=R.index)
    for c in ["dens_own_1h", "dens_own_10m", "dens_grp_1h", "dens_grp_10m", "n_oth_4s"]: N[c] = 0
    N = pd.concat([N, VN], axis=1); N["dens_recompute"] = True
    N["feat_src"] = np.where(useF, "NK_net_exact+site_from_cut_peak", "NK_net_exact_only(fallback)")
    # (4) assemble: drop the with-KMR B/C DATA rows from the v2 forward series (kept apart for the consistency check)
    isBCd = B.seg.isin(["B", "C"]) & B.series.eq("data")
    drop = B.seg.isin(["B", "C", "NK"]) & B.series.eq("data")      # replaced by the mf_NK rows (switch days included)
    BC = B[isBCd].copy(); BC["run"] = "fwdBC"
    F = pd.concat([B[~drop], N], ignore_index=True)
    F = context_v2(F)
    BC["dens_recompute"] = False; BC = context_v2(BC)
    log(f"assembled v2 rows {len(F)} ({F.groupby(['run', 'series']).size().to_dict()}); BC check rows {len(BC)}; {time.time() - t00:.0f}s")
    # (5) fit: same design as sh_conf_fit (stored folds), FULL_v2 and FALLBACK_v2
    PNW = FT.PNW; PNF = FT.PNF
    FS = {"FULL_v2": FT.MF + FT.CTX + FT.S1 + PNW + PNF + FT.FAM,
          "FALLBACK_v2": ["net_mad", "meancc", "ccsum", "nlive", "live_frac", "nchan", "nsta"] + [c for c in FT.CTX if c != "hfr_rel"] + FT.S1 + FT.FAM}
    data = F.series.eq("data").to_numpy(); role = F.role.to_numpy(); dates = F.date.to_numpy(); fold = F.fold.to_numpy()
    fitm = data & (role == "fit")
    MODELS = {}
    for fs, feats in FS.items():
        X = F[feats].to_numpy(float); Xb = BC[feats].to_numpy(float)
        for lab, lc in FT.LABELS.items():
            y = F[lc].to_numpy(); pred = np.full(len(F), np.nan)
            for k in range(FT.NFOLD):
                tr = fitm & (fold != k)
                m = FT.fit_cal(feats, X[tr], y[tr].astype(int), dates[tr])
                sc = (role == "fit") & (fold == k); pred[sc] = m.predict_proba(X[sc])[:, 1]
            m = FT.fit_cal(feats, X[fitm], y[fitm].astype(int), dates[fitm])
            sc = role != "fit"; pred[sc] = m.predict_proba(X[sc])[:, 1]
            F[f"P_{lab}_{fs}"] = pred; BC[f"P_{lab}_{fs}"] = m.predict_proba(Xb)[:, 1]; MODELS[f"{fs}_{lab}"] = m
            log(f"fit {fs} {lab}: n_fit {int(fitm.sum())} pos {int(np.nansum(y[fitm]))}; {time.time() - t00:.0f}s")
    import joblib; joblib.dump(dict(models=MODELS, features=FS), f"{OUT}/models_v2.joblib")
    # final P_v2: FULL_v2 where site features were available, else FALLBACK_v2
    fb = F.feat_src.eq("NK_net_exact_only(fallback)").to_numpy()
    for lab in FT.LABELS:
        F[f"P_{lab}_v2"] = np.where(fb, F[f"P_{lab}_FALLBACK_v2"], F[f"P_{lab}_FULL_v2"])
    F["P_model_v2"] = np.where(fb, "FALLBACK_v2", "FULL_v2")
    # (6) evaluation vs the original FULL model on the same rows
    rows = []
    for es in FT.EVAL:
        em = (role == ("fit" if es == "fit_oof" else es))
        Dd = F[em & data]; Cc = F[em & ~data]
        for lab, lc in FT.LABELS.items():
            y = Dd[lc].to_numpy(float)
            for mname in ["FULL", "FULL_v2", "FALLBACK_v2"]:
                p = Dd[f"P_{lab}_{mname}"].to_numpy(float); ok = np.isfinite(p) & np.isfinite(y)
                if not ok.any(): continue
                auc, ap = FT.auc_ap(y[ok], p[ok]); sl, ic, _ = FT.calib(y[ok], p[ok]); pc = Cc[f"P_{lab}_{mname}"].to_numpy(float)
                rows.append(dict(eval_set=es, label=lab, model=mname, n=int(ok.sum()), n_days=Dd.date.nunique(), base_rate=float(y[ok].mean()), auc=auc, ap=ap,
                                 slope=sl, icpt=ic, matched=int(y[ok].sum()), sumP=float(p[ok].sum()), obs_over_exp=float(y[ok].sum() / p[ok].sum()),
                                 brier=float(np.mean((p[ok] - y[ok]) ** 2)), data_fracP_gt05=float((p[ok] > 0.5).mean()),
                                 ctl_fracP_gt05=float(np.mean(pc[np.isfinite(pc)] > 0.5)) if np.isfinite(pc).any() else np.nan,
                                 ctl_over_data_sumP=float(np.nansum(pc) / p[ok].sum())))
    EV = pd.DataFrame(rows); EV.to_csv(f"{OUT}/t3_eval.csv", index=False); log("\n" + EV.round(3).to_string())
    # (7) consistency: v2 model on the B/C (approx) rows vs on the real 16-station rows, same days, per family
    Fn = F[(F.seg == "NK") & (F.date > "2025-09-17") & (F.date != "2026-03-12")]
    cons = pd.DataFrame(dict(sumP_same_v2_NK=Fn.groupby("fam").P_same_v2.sum(), sumP_same_v2_BCapprox=BC.groupby("fam").P_same_FULL_v2.sum(),
                             sumP_same_v1_BC=A[A.seg.isin(["B", "C"]) & A.series.eq("data")].groupby("fam").P_same_FULL.sum(),
                             n_NK=Fn.groupby("fam").size(), n_BC=BC.groupby("fam").size()))
    cons["NK_over_BCapprox"] = cons.sumP_same_v2_NK / cons.sumP_same_v2_BCapprox; cons["NK_v2_over_BC_v1"] = cons.sumP_same_v2_NK / cons.sumP_same_v1_BC
    cons.loc["all"] = cons.sum(numeric_only=True); cons.loc["all", "NK_over_BCapprox"] = cons.loc["all", "sumP_same_v2_NK"] / cons.loc["all", "sumP_same_v2_BCapprox"]
    cons.loc["all", "NK_v2_over_BC_v1"] = cons.loc["all", "sumP_same_v2_NK"] / cons.loc["all", "sumP_same_v1_BC"]
    cons.to_csv(f"{OUT}/t3_consistency_NK_vs_BC.csv"); log("\n" + cons.round(3).tail(5).to_string())
    # (8) save rows + joinable P table
    sv = ["series", "run", "date", "fam", "seg", "role", "build_any", "t_epoch", "det_id", "det_id_v2", "feat_src", "m_same", "m_any",
          "net_mad", "meancc", "ccsum", "nlive", "nchan", "nsta", "in_S1_9x", "P_same_FULL", "P_any_FULL", "P_same_FULL_v2", "P_any_FULL_v2",
          "P_same_FALLBACK_v2", "P_any_FALLBACK_v2", "P_same_v2", "P_any_v2", "P_model_v2"]
    F[sv].rename(columns={"in_S1_9x": "in_S1_9x_v2feat"}).to_parquet(f"{OUT}/conf_v2_rows.parquet", index=False)
    BC[["date", "fam", "seg", "t_epoch", "P_same_FULL", "P_same_FULL_v2", "P_any_FULL_v2"]].to_parquet(f"{OUT}/conf_v2_BCcheck_rows.parquet", index=False)
    log(f"done {time.time() - t00:.0f}s")


if __name__ == "__main__":
    main()
