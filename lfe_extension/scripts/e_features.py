#!/usr/bin/env python
"""Track E step 2 — assemble the per-detection FEATURE table (labels are NEVER features).

  (1) independent-site breadth   n_sites_above (z_site >= Z_ABOVE), frac_sites_above, z_site_med,
                                 z_site_min, z_ex_top (mean site z EXCLUDING the heaviest site)
  (2) independent confirmation   the N_LOW sites with the SMALLEST template-weight share for that
                                 family — the ones that contributed least to triggering: conf_z
                                 (their own cc in units of their own daily MAD) and conf_ze (the
                                 S-window energy excess at the predicted moveout, in units of that
                                 site's own daily distribution)
  (3) burst clustering           n_same10 / n_same60 (same family, +-10 / +-60 min), n_other5
                                 (OTHER families, +-5 min), and the excess of n_same60 over the
                                 Poisson expectation of that family-day's own rate
  (4) physics                    net_mad, cc, rel_amp, hfr, n_sites, wsh_pgc
Context carried but NOT given to the model: tremor_flag, quiet_flag, season, epoch, tier, date.

Burst features are computed on the WHOLE record at once (never per month), so no detection loses
neighbours at a month boundary, and separately inside each "world" (data / site-rolled control).

  python scripts/e_features.py --prod --out results/trackE/feat_prod.parquet
  python scripts/e_features.py --ref 2012 --det results/e5r_detections_a1_mergeB_L12.parquet \
      --persite results/trackE/persite_ref2012.parquet --out results/trackE/feat_ref2012.parquet
"""
from __future__ import annotations
import os, sys, glob, argparse
import numpy as np
import pandas as pd

Z_ABOVE = 3.0       # "exceeds its own site MAD" — fixed by declaration, never tuned on the labels
N_LOW = 2           # how many lowest-weight sites define the independent confirmation


def _same_fam_counts(t, fam, w):
    """count of same-family detections within +-w s, self excluded."""
    n = len(t)
    out = np.zeros(n)
    o = np.lexsort((t, fam))
    ts = t[o]; fs = fam[o]
    b = np.flatnonzero(np.r_[True, fs[1:] != fs[:-1], True])
    for i in range(len(b) - 1):
        lo, hi = b[i], b[i + 1]
        tt = ts[lo:hi]
        out[o[lo:hi]] = (np.searchsorted(tt, tt + w, "right") -
                         np.searchsorted(tt, tt - w, "left") - 1)
    return out


def burst_features(D):
    t = D.t.to_numpy(float); fam = D.fam.to_numpy()
    n10 = _same_fam_counts(t, fam, 600.0)
    n60 = _same_fam_counts(t, fam, 3600.0)
    n5s = _same_fam_counts(t, fam, 300.0)
    o = np.argsort(t); ts = t[o]
    allc = (np.searchsorted(ts, ts + 300.0, "right") -
            np.searchsorted(ts, ts - 300.0, "left") - 1)
    n_all5 = np.empty(len(t)); n_all5[o] = allc
    nfd = D.groupby(["fam", "date"]).t.transform("size").to_numpy(float)
    mu = (nfd - 1) * (7200.0 / 86400.0)
    return pd.DataFrame(dict(n_same10=n10, n_same60=n60, n_other5=n_all5 - n5s,
                             burst_mu=mu, burst_excess=n60 - mu,
                             burst_z=(n60 - mu) / np.sqrt(mu + 1.0)), index=D.index)


def persite_features(PS):
    key = ["fam", "date", "t", "ctrl"]
    PS = PS.sort_values(key + ["wsh"], kind="stable").reset_index(drop=True)
    g = PS.groupby(key, sort=False)
    # a normalised correlation cannot exceed 1: on 18 station-days one station (mostly UW.LRIV
    # before it was dropped) has a near-zero-energy interval that slips past the energy floor and
    # its cc blows up.  Track C's own stored net_cc carries the same blow-up, so this is a data-QC
    # flag on those detections, not a discrepancy with production.  They are never scored.
    PS["_bad"] = (~np.isfinite(PS.cc_site.to_numpy())) | (np.abs(PS.cc_site.to_numpy()) > 1.0)
    F = g.agg(n_sites_obs=("site", "size"), persite_bad=("_bad", "max"),
              z_site_med=("z_site", "median"),
              z_site_min=("z_site", "min"), z_site_max=("z_site", "max"),
              ze_site_med=("ze_site", "median"), ze_site_min=("ze_site", "min")).reset_index()
    ab = (PS.z_site.to_numpy() >= Z_ABOVE).astype(float)
    F["n_sites_above"] = pd.Series(ab).groupby(g.ngroup().to_numpy(), sort=False).sum().to_numpy()
    F["frac_sites_above"] = F.n_sites_above / F.n_sites_obs.clip(lower=1)
    r = g.cumcount().to_numpy()
    nq = g.size().to_numpy()
    rep = np.repeat(nq, nq)
    low = PS[r < N_LOW]
    Fl = low.groupby(key, sort=False).agg(conf_z=("z_site", "mean"), conf_ze=("ze_site", "mean"),
                                          conf_wsh=("wsh", "mean"), n_conf=("site", "size")).reset_index()
    top = PS[r == rep - 1]
    Ft = top.groupby(key, sort=False).agg(top_z=("z_site", "mean"), top_wsh=("wsh", "mean")).reset_index()
    F = F.merge(Fl, on=key, how="left").merge(Ft, on=key, how="left")
    zsum = PS.groupby(key, sort=False).z_site.sum().reset_index().rename(columns={"z_site": "z_sum"})
    F = F.merge(zsum, on=key, how="left")
    F["z_ex_top"] = (F.z_sum - F.top_z) / (F.n_sites_obs - 1).clip(lower=1)
    return F.drop(columns=["z_sum"])


CARRY = ["fam", "date", "doy", "t", "net_mad", "net_cc", "rel_amp", "n_sta", "mad_scale", "hfr",
         "hf_pass", "tremor_flag", "quiet_flag", "n_tremor_day", "eq_flag", "tier_ext", "tier_core",
         "cat_ext", "cat_core", "epoch", "match_own", "match_any"]
# n_sites and wsh_pgc are deliberately NOT taken from the detection file: the referee runs spell the
# site columns differently (wsh_PGC) and the control has none, so both are re-derived from the
# per-site table, identically in every world.


def load_det(paths, year=0, record_thr=8.0):
    fr = []
    for p in paths:
        d = pd.read_parquet(p)
        if "L" in d.columns:
            d = d[d.L == 12.0]
        d = d.rename(columns={c: c.replace("hfr_PGC", "hfr") for c in d.columns if c == "hfr_PGC"})
        if "amp" in d.columns and "rel_amp" not in d.columns:
            d = d.rename(columns={"amp": "rel_amp"})
        if "date" not in d.columns:
            d["date"] = (pd.Timestamp(f"{year}-01-01") +
                         pd.to_timedelta(d.doy.astype(int) - 1, unit="D")).dt.strftime("%Y-%m-%d")
        fr.append(d[[c for c in CARRY if c in d.columns]])
    D = pd.concat(fr, ignore_index=True)
    D = D[D.net_mad >= record_thr].copy()
    D["fam"] = D.fam.astype(int)
    D["ctrl"] = 0
    return D


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prod", action="store_true")
    ap.add_argument("--ref", type=int, default=0)
    ap.add_argument("--det", default="")
    ap.add_argument("--persite", default="")
    ap.add_argument("--persite-det", default="")
    ap.add_argument("--out", required=True)
    ap.add_argument("--record-thr", type=float, default=8.0)
    a = ap.parse_args()
    RES = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results")

    if a.prod:
        dets = sorted(glob.glob(f"{RES}/trackC/detections_*.parquet"))
        pss = sorted(p for p in glob.glob(f"{RES}/trackE/persite_[0-9]*.parquet")
                     if "_det" not in p)
        D = load_det(dets, record_thr=a.record_thr)
        PS = pd.concat([pd.read_parquet(p) for p in pss], ignore_index=True)
        frames = [D]
    else:
        D = load_det([a.det], year=a.ref, record_thr=a.record_thr)
        PS = pd.read_parquet(a.persite)
        frames = [D]
        pdt = a.persite_det or a.persite.replace(".parquet", "_det.parquet")
        if os.path.exists(pdt):
            PD = pd.read_parquet(pdt)
            Cc = PD[PD.ctrl == 1].copy()
            if len(Cc):
                Cc["fam"] = Cc.fam.astype(int)
                Cc["rel_amp"] = np.nan; Cc["hfr"] = np.nan; Cc["wsh_pgc"] = np.nan
                for c in ("tremor_flag", "quiet_flag", "eq_flag", "match_own", "match_any"):
                    Cc[c] = False
                frames.append(Cc[[c for c in D.columns if c in Cc.columns]])
    PS["fam"] = PS.fam.astype(int)
    F = persite_features(PS)
    B = []
    for fr in frames:
        fr = fr.reset_index(drop=True)
        B.append(pd.concat([fr, burst_features(fr)], axis=1))
    X = pd.concat(B, ignore_index=True)
    X["ctrl"] = X.ctrl.astype(np.int8)
    F["ctrl"] = F.ctrl.astype(np.int8)
    X = X.merge(F, on=["fam", "date", "t", "ctrl"], how="left")
    W = PS[PS.site == "PGC"][["fam", "date", "t", "ctrl", "wsh"]].rename(columns={"wsh": "wsh_pgc"})
    W["ctrl"] = W.ctrl.astype(np.int8)
    X = X.merge(W.drop_duplicates(["fam", "date", "t", "ctrl"]), on=["fam", "date", "t", "ctrl"],
                how="left")
    # HF ratio, made epoch-stable: the PGC anti-alias roll-off moves the whole distribution at the
    # BH->HH change (Track C README caveat), so the model is given the ratio RELATIVE to its own
    # month's median rather than the absolute one.
    lh = np.log10(X.hfr.to_numpy(float)) if "hfr" in X.columns else np.full(len(X), np.nan)
    X["log_hfr"] = lh
    ym = pd.Index(X.date.astype(str).str[:7])
    med = pd.Series(lh).groupby(ym).transform("median").to_numpy()
    X["hfr_rel"] = lh - med
    X.to_parquet(a.out, index=False)
    miss = int(X.n_sites_obs.isna().sum())
    print(f"[feat] {a.out}  {len(X)} rows ({int((X.ctrl==0).sum())} data / {int((X.ctrl==1).sum())} "
          f"control); per-site features missing on {miss} ({100*miss/max(len(X),1):.2f} %)", flush=True)


if __name__ == "__main__":
    main()
