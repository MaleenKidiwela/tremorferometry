#!/usr/bin/env python
"""Track D, D1(b) -- DAILY SENSITIVITY: a per-family-day completeness magnitude Mc(f, d) at both tiers
(Merlin Step 6, component 2), calibrated on the referee year and VALIDATED out of sample on 2010.

MODEL (probit, one fit per station-set/epoch run and per tier T)

    P(detect | family f, day d, magnitude M) = Phi( b_f + alpha*(M - 1.8) + c1*u(f,d) + c2*v(f,d) )

  b_f     per-family intercept (ridge-shrunk toward the pooled value; the ridge weight is chosen by
          5-fold CV on the referee events).  It carries BOTH the family's intrinsic detectability and
          the station set of the epoch -- that is the same quantity D1(a) measures directly as recall,
          so component 1 and component 2 are one model, not two.
  alpha   magnitude slope on the probit scale.  E4 measured an AMPLITUDE slope of 10^(0.38 dM) within
          family; the recall slope is correspondingly shallow, so Mc is a wide, soft threshold here and
          is reported as such.
  u(f,d)  log10 of the family-day network MAD scale divided by that family's median over the same
          station set  -- a noisier cc series needs a bigger event.  Stored by Track C as `mad_scale`.
  v(f,d)  log10 of the SNR-weighted noise AMPLITUDE factor over the stations that carried a template
          that day: v = 0.5*log10( sum_s w(f,s)^2 * e8_s(d)/e8_ref(s,epoch) / sum_s w(f,s)^2 ), with
          e8 the daily median 8-s window energy (3-component mean) of `daily_noise_*.parquet` and
          e8_ref the median over that station's INSTRUMENT epoch (so a gain change cancels: the 20 Hz
          chain removes no response, and cc is amplitude-normalised, so only the RELATIVE noise matters).
          Physical anchor: in the low-SNR regime cc is proportional to 1/noise, so c2 should come out near -1
          in the units of the fitted statistic.

  Mc(f,d,T) = 1.8 - (b_f + c1*u + c2*v)/alpha     (the magnitude at which P = 0.5)

VALIDATION (the pre-registered bar): fit on the 2012 ablation whose station set equals the 2010 network
(9 stations, tag a2_abl2012_7p), predict 2010 recall per magnitude bin, compare with measured.  Bar: |predicted
- measured| <= 0.05 in every bin."""
from __future__ import annotations
import argparse, json, os
import numpy as np, pandas as pd
from scipy.stats import norm
from scipy.optimize import minimize
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d_tiers  # noqa: E402

IN = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
EXT = "/home/jovyan/tremorferometry/lfe_extension"
OUT = f"{EXT}/results"
TD = f"{OUT}/trackD"
MAGBINS = [(-9, 1.6, "M<1.6"), (1.6, 1.8, "M1.6-1.8"), (1.8, 2.0, "M1.8-2.0"), (2.0, 9, "M>=2.0")]
NETSTA = dict(PGC="CN.PGC", VGZ="CN.VGZ", NLLB="CN.NLLB", SNB="CN.SNB", LRIV="UW.LRIV", B001="PB.B001",
              B003="PB.B003", B004="PB.B004", B005="PB.B005", B006="PB.B006", B007="PB.B007",
              B009="PB.B009", B010="PB.B010", B011="PB.B011", B926="PB.B926")
RUNS = {  # station-set id -> (tag, referee year, stations)
    "S15": ("a1_mergeB_L12", 2012, "PGC,VGZ,NLLB,SNB,LRIV,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926"),
    "S14": ("d1_abl2012_E2", 2012, "PGC,VGZ,NLLB,SNB,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926"),
    "S13": ("d1_abl2012_E3", 2012, "PGC,VGZ,NLLB,B001,B003,B004,B005,B006,B007,B009,B010,B011,B926"),
    "P9": ("a2_abl2012_7p", 2012, "PGC,VGZ,SNB,NLLB,B926,LRIV,B009,B010,B011"),
    "S15_2010": ("a2_xfer2010", 2010, "PGC,VGZ,SNB,NLLB,B926,LRIV,B009,B010,B011"),
}
SET_OF_RUN = {"S15": "S15", "S14": "S14", "S13": "S13", "P9": "S15", "S15_2010": "S15"}
HF_CUT = 3.6384


# ----------------------------------------------------------------- data
def one2one(a, b_sorted, tol=2.0):
    a = np.asarray(a, float); out = np.full(len(a), np.inf)
    if len(a) == 0 or len(b_sorted) == 0:
        return out
    j = 0; nb = len(b_sorted)
    for i in np.argsort(a, kind="stable"):
        ti = a[i]
        while j < nb and b_sorted[j] < ti - tol:
            j += 1
        if j < nb and abs(b_sorted[j] - ti) <= tol:
            out[i] = abs(b_sorted[j] - ti); j += 1
    return out


def load_det(tag, L=12.0):
    det = pd.read_parquet(f"{OUT}/e5r_detections_{tag}.parquet"); det = det[det.L == L]
    hp = f"{OUT}/a1_hf_ratio_{tag}.parquet"
    if os.path.exists(hp):
        H = pd.read_parquet(hp); H = H[H.L == L][["doy", "fam", "t", "hfr_PGC"]]
        det = det.merge(H, on=["doy", "fam", "t"], how="left")
        v = det.hfr_PGC.to_numpy(); det = det[~(np.isfinite(v) & (v > HF_CUT))]
    return det.reset_index(drop=True)


def gt(year):
    if year == 2012:
        g = pd.read_parquet(f"{IN}/e5r_gt_bostock2012.parquet").reset_index(drop=True)
        g = g[(g.half == "test") & (~g.guard)].copy()
        g["time"] = g.time.astype("datetime64[ns]")
        g["doy"] = pd.DatetimeIndex(g.time).dayofyear
    else:
        g = pd.read_parquet(f"{OUT}/a2_gt_bostock2010.parquet")
        off = pd.read_csv(f"{IN}/e3_bostock_family_offsets.csv"); off = off[off.accept.astype(bool)]
        om = dict(zip(off.bkey.astype(int), off.offset_s.astype(float)))
        g = g[g.bkey.astype(int).isin(om)].copy()
        g["time"] = g.time.astype("datetime64[ns]")
        g["t_origin"] = g.time.astype("int64") / 1e9 + g.bkey.astype(int).map(om).to_numpy()
        g["doy"] = pd.DatetimeIndex(g.time).dayofyear
    g["bkey"] = g.bkey.astype(int)
    return g.reset_index(drop=True)


def noise_ref(D):
    """per (netsta, instrument epoch) median of the 3-component mean e8 -> the gain-free reference."""
    D = D[D.status.isin(["cached", "built"])].copy()
    D["e8"] = D[["e8_Z", "e8_N", "e8_E"]].mean(axis=1)
    return D


def day_terms(runid, det, ev, DN, W, TPL):
    """u(f,d) and v(f,d) for the reference events of one run.  The family-day MAD scale comes from the
    BUDGET file (`madn`), which has a row for every scored family-day, not only the ones with a
    detection."""
    tag, year, stas = RUNS[runid]
    stal = stas.split(",")
    b = pd.read_csv(f"{OUT}/e5_budget_{tag}.csv"); b = b[b.L == 12.0]
    b = b[b["T"] == sorted(b["T"].unique())[0]]
    ms = b[["fam", "doy", "madn"]].drop_duplicates(["fam", "doy"]).rename(columns={"madn": "mad_scale"})
    ms["bkey"] = ms.fam.astype(int)
    medf = ms.groupby("bkey").mad_scale.median().rename("mad_ref")
    ms = ms.merge(medf, on="bkey")
    ms["u"] = np.log10(ms.mad_scale / ms.mad_ref)
    # v: SNR-weighted relative noise amplitude over the template-carrying stations present that day
    N = DN[(DN.year == year) & (DN.netsta.isin([NETSTA[s] for s in stal]))].copy()
    ref = N.groupby("netsta").e8.median()
    N["r"] = N.e8 / N.netsta.map(ref)
    rmap = {(r.netsta, int(r.doy)): float(r.r) for r in N.itertuples()}
    v = {}
    for bk in ev.bkey.unique():
        ws = {s: W.get((float(bk), s), np.nan) for s in stal if s in TPL.get(int(bk), set())}
        ws = {s: w for s, w in ws.items() if np.isfinite(w)}
        if not ws:
            continue
        for d in sorted(ev[ev.bkey == bk].doy.unique()):
            num = den = 0.0
            for s, w in ws.items():
                r = rmap.get((NETSTA[s], int(d)))
                if r is None or not np.isfinite(r):
                    continue
                num += w * w * r; den += w * w
            v[(int(bk), int(d))] = 0.5 * np.log10(num / den) if den > 0 else np.nan
    return ms.set_index(["bkey", "doy"]).u.to_dict(), v


def build_xy(runid, T, DN, W, TPL):
    tag, year, stas = RUNS[runid]
    det = load_det(tag); ev = gt(year)
    D = det[det.net_mad >= T]
    hit = np.zeros(len(ev), bool)
    for bk, g in ev.groupby("bkey"):
        td = np.sort(D[D.fam == float(bk)].t.to_numpy())
        hit[ev.index.get_indexer(g.index)] = one2one(g.t_origin.to_numpy(), td) <= 2.0
    umap, vmap = day_terms(runid, det, ev, DN, W, TPL)
    ev = ev.copy(); ev["y"] = hit.astype(float)
    ev["u"] = [umap.get((int(b), int(d)), np.nan) for b, d in zip(ev.bkey, ev.doy)]
    ev["v"] = [vmap.get((int(b), int(d)), np.nan) for b, d in zip(ev.bkey, ev.doy)]
    ev["u"] = ev.u.fillna(0.0); ev["v"] = ev.v.fillna(0.0)
    return ev


# ----------------------------------------------------------------- probit with a ridge on the family effects
def fit_probit(ev, lam, fams=None):
    fams = fams if fams is not None else np.array(sorted(ev.bkey.unique()))
    fidx = {f: i for i, f in enumerate(fams)}
    k = np.array([fidx[b] for b in ev.bkey])
    M = (ev.mag.to_numpy() - 1.8); U = ev.u.to_numpy(); V = ev.v.to_numpy(); y = ev.y.to_numpy()
    nf = len(fams)

    def nll(p):
        b0, alpha, c1, c2 = p[0], p[1], p[2], p[3]
        b = p[4:]
        eta = b0 + b[k] + alpha * M + c1 * U + c2 * V
        lp = norm.logcdf(eta); lq = norm.logcdf(-eta)
        f = -(y * lp + (1 - y) * lq).sum() + lam * (b * b).sum()
        g = np.zeros_like(p)
        w = y * np.exp(norm.logpdf(eta) - lp) - (1 - y) * np.exp(norm.logpdf(eta) - lq)
        g[0] = -w.sum(); g[1] = -(w * M).sum(); g[2] = -(w * U).sum(); g[3] = -(w * V).sum()
        gb = -np.bincount(k, weights=w, minlength=nf) + 2 * lam * b
        g[4:] = gb
        return f, g

    p0 = np.zeros(4 + nf); p0[0] = norm.ppf(np.clip(y.mean(), 1e-3, 1 - 1e-3))
    r = minimize(nll, p0, jac=True, method="L-BFGS-B", options=dict(maxiter=2000, maxfun=4000))
    return dict(b0=r.x[0], alpha=r.x[1], c1=r.x[2], c2=r.x[3],
                b={int(f): float(r.x[4 + i]) for i, f in enumerate(fams)}, fams=fams, nll=float(r.fun),
                converged=bool(r.success))


def predict(fit, bkey, mag, u, v):
    b = np.array([fit["b"].get(int(f), 0.0) for f in bkey])
    eta = fit["b0"] + b + fit["alpha"] * (np.asarray(mag) - 1.8) + fit["c1"] * np.asarray(u) + fit["c2"] * np.asarray(v)
    return norm.cdf(eta)


def cv_lambda(ev, lams, k=5, seed=7):
    rng = np.random.default_rng(seed)
    fold = rng.integers(0, k, len(ev))
    fams = np.array(sorted(ev.bkey.unique()))
    out = []
    for lam in lams:
        ll = 0.0
        for f in range(k):
            tr, te = ev[fold != f], ev[fold == f]
            fit = fit_probit(tr, lam, fams)
            p = np.clip(predict(fit, te.bkey, te.mag, te.u, te.v), 1e-6, 1 - 1e-6)
            ll += float((te.y * np.log(p) + (1 - te.y) * np.log(1 - p)).sum())
        out.append((lam, ll))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lams", default="0.3,1,3,10,30,100")
    a = ap.parse_args()
    os.makedirs(TD, exist_ok=True)
    EPOCHS, TIERS, TSRC = d_tiers.load()
    TH_SET = {ss: d_tiers.thresholds_for_set(ss, EPOCHS, TIERS) for ss in ("S15", "S14", "S13")}
    print(f"[tiers] {TSRC}")
    DN = noise_ref(pd.read_parquet(f"{TD}/daily_noise_ref.parquet"))
    w = pd.read_csv(f"{IN}/e5r_weights_snr_2012.csv")
    W = {(float(r.fam), r.sta): float(r.w) for r in w.itertuples()}
    qc = pd.read_csv(f"{IN}/e5r_template_qc_2012.csv")
    if "usable" in qc.columns:
        qc = qc[qc.usable.astype(bool)]
    TPL = qc.groupby("bkey").sta.apply(lambda s: set(s)).to_dict()

    # export axis = Track C's EPOCHS, each fitted at ITS OWN re-derived thresholds through the referee
    # run of its station set; the validation axis stays the station set (P9 -> held-out 2010).
    JOBS = [(e["epoch_id"], e["station_set"], e["epoch_id"]) for e in EPOCHS] + \
           [("P9", "P9", None), ("S15_2010", "S15_2010", None)]
    fits, rows, valrows = {}, [], []
    for nm in ("T_ext", "T_core"):
        for key, runid, eid in JOBS:
            tag = RUNS[runid][0]
            if not os.path.exists(f"{OUT}/e5r_detections_{tag}.parquet"):
                print(f"[skip] {key}: {tag} not on disk"); continue
            T = TIERS[eid][nm] if eid else TH_SET[SET_OF_RUN[runid]][nm]
            ev = build_xy(runid, T, DN, W, TPL)
            lams = [float(x) for x in a.lams.split(",")]
            cv = cv_lambda(ev, lams) if key in (EPOCHS[0]["epoch_id"], "P9") else None
            lam = max(cv, key=lambda z: z[1])[0] if cv else fits[(nm, EPOCHS[0]["epoch_id"])]["lam"]
            fit = fit_probit(ev, lam); fit["lam"] = lam; fit["cv"] = cv; fit["n"] = len(ev)
            fit["recall_meas"] = float(ev.y.mean()); fit["recall_pred"] = float(predict(fit, ev.bkey, ev.mag, ev.u, ev.v).mean())
            fits[(nm, key)] = fit
            fit["u_p01"], fit["u_p99"] = [float(x) for x in np.percentile(ev.u, [1, 99])]
            fit["v_p01"], fit["v_p99"] = [float(x) for x in np.percentile(ev.v, [1, 99])]
            print(f"[{nm} {key}] run {runid} T={T} n={len(ev)} lam={lam} alpha={fit['alpha']:.3f} c1={fit['c1']:.3f} "
                  f"c2={fit['c2']:.3f} recall meas {fit['recall_meas']:.4f} pred {fit['recall_pred']:.4f}", flush=True)
            for f_, b_ in fit["b"].items():
                rows.append(dict(tier=nm, T=T, epoch_id=key, station_set=runid, family_id=f_, b_family=b_,
                                 b0=fit["b0"], alpha=fit["alpha"], c1=fit["c1"], c2=fit["c2"], lam=lam))
        # ---- VALIDATION: fit on the 2012 9-station ablation, predict the held-out 2010 year
        if (nm, "P9") in fits and (nm, "S15_2010") in fits:
            T10 = TH_SET["S15"][nm]
            ev10 = build_xy("S15_2010", T10, DN, W, TPL)
            fitP9 = fits[(nm, "P9")]
            p = predict(fitP9, ev10.bkey, ev10.mag, ev10.u, ev10.v)
            for lo, hi, lab in MAGBINS + [(-9, 9, "all")]:
                k = (ev10.mag >= lo) & (ev10.mag < hi)
                if k.sum() < 20:
                    continue
                valrows.append(dict(tier=nm, bin=lab, n=int(k.sum()), measured=float(ev10.y[k].mean()),
                                    predicted=float(p[k.to_numpy()].mean())))
            # in-sample control on 2012 with the same station set
            evP9 = build_xy("P9", TH_SET["S15"][nm], DN, W, TPL)
            pin = predict(fitP9, evP9.bkey, evP9.mag, evP9.u, evP9.v)
            for lo, hi, lab in MAGBINS + [(-9, 9, "all")]:
                k = (evP9.mag >= lo) & (evP9.mag < hi)
                if k.sum() < 20:
                    continue
                valrows.append(dict(tier=nm, bin=lab + " (2012 in-sample)", n=int(k.sum()),
                                    measured=float(evP9.y[k].mean()), predicted=float(pin[k.to_numpy()].mean())))
    P = pd.DataFrame(rows)
    rng = {(k[0], k[1]): (v["u_p01"], v["u_p99"], v["v_p01"], v["v_p99"]) for k, v in fits.items()
           if "u_p01" in v}
    P["u_calib_lo"] = [rng.get((t, e), (np.nan,) * 4)[0] for t, e in zip(P.tier, P.epoch_id)]
    P["u_calib_hi"] = [rng.get((t, e), (np.nan,) * 4)[1] for t, e in zip(P.tier, P.epoch_id)]
    P["v_calib_lo"] = [rng.get((t, e), (np.nan,) * 4)[2] for t, e in zip(P.tier, P.epoch_id)]
    P["v_calib_hi"] = [rng.get((t, e), (np.nan,) * 4)[3] for t, e in zip(P.tier, P.epoch_id)]
    P.to_csv(f"{TD}/mc_model_params.csv", index=False)
    V = pd.DataFrame(valrows)
    V["delta"] = V.predicted - V.measured
    V["within_0.05"] = V.delta.abs() <= 0.05
    V.to_csv(f"{TD}/mc_model_validation.csv", index=False)
    print("\n[VALIDATION -- fit on 2012 (9-station ablation), predicted on the held-out 2010 year]")
    print(V.round(4).to_string(index=False))
    json.dump({f"{k[0]}|{k[1]}": {kk: vv for kk, vv in v.items() if kk not in ("b", "fams")}
               for k, v in fits.items()}, open(f"{TD}/mc_model_fits.json", "w"), indent=1, default=str)
    print(f"[write] {TD}/mc_model_params.csv, mc_model_validation.csv, mc_model_fits.json")


if __name__ == "__main__":
    main()
