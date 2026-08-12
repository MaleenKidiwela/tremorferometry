#!/usr/bin/env python
"""B011 (borehole) vs PGC (broadband, ~290 m) co-located validation of the broadband dv/v (Merlin spec).
The raw dv/v correlation is confounded by the shared noise field (65% artifact at 290 m), so the verdict lives
in MIRROR-CORRECTED, deseasoned residuals, referenced against a PGC-mirror-vs-B011-coda NULL (same triggers/
days/noise, zero LFE energy). Replicates cm()/des()/beta from plot_mirror_corrected_dvv.py exactly.

Per era e in {BHZ 2010-2017, HHZ 2018-2026}:
  F_e = (PGC-era causality-cert) ∩ (B011 fwd_vs_rev ratio>1.5).  Build cm() medians for PGC/B011 x coda/mirror.
  R_X = des(coda_X) - beta_X*des(mirror_X).   S1 = spearman(R_PGC, R_B011).  N1 = spearman(des(PGC mirror), R_B011).
  Block-bootstrap (block=60d) CIs on S1, slope b, and PAIRED (S1-N1).
  PASS(e): S1 CI excludes 0 AND (S1-N1) CI excludes 0 AND b CI excludes 0 with point in [0.3,3]; S1<0.2 -> INCONCLUSIVE.
Program validated iff BOTH eras PASS.  Plus per-family matched-vs-mismatched Delta, and the B011 seam-step test.
Usage: python scripts/b011_pgc_validation.py
"""
import os, sys, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr

B011_LA, B011_LO = 48.65, -123.45
ERAS = [("BHZ", "PGCbhz", "2010-01-01", "2017-09-01"),
        ("HHZ", "PGChhz", "2017-09-01", "2026-12-31")]
SEAM = pd.Timestamp("2017-08-15")


# ---- protocol replicated verbatim from plot_mirror_corrected_dvv.py ----
def cm(csv, fams):
    d = pd.read_csv(csv); d = d[d.patch.isin(fams) & (d.cc_max >= 0.6)]; d["date"] = pd.to_datetime(d.date)
    m = d.groupby("date").dvv.median(); n = d.groupby("date").patch.nunique()
    return (m[n >= 3].rolling(15, center=True, min_periods=5).median()) * 100

def des(s):
    s = s.dropna()
    if len(s) < 90: return s
    t = (s.index - s.index[0]).days.values.astype(float); yr = 365.25
    X = np.column_stack([np.ones_like(t), t, np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr),
                         np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
    b, *_ = np.linalg.lstsq(X, s.values, rcond=None); return pd.Series(s.values - X@b, index=s.index)

def corrected(coda_csv, mir_csv, fams, lo=None, hi=None):
    """des(coda) - beta*des(mirror) on common days, restricted to [lo,hi]. Returns (R, des_coda, des_mir, beta)."""
    fwd = des(cm(coda_csv, fams)); rev = des(cm(mir_csv, fams))
    j = pd.concat([fwd.rename("f"), rev.rename("r")], axis=1).dropna()
    if lo: j = j[j.index >= pd.Timestamp(lo)]
    if hi: j = j[j.index <= pd.Timestamp(hi)]
    if len(j) < 90: return None
    beta = np.polyfit(j.r, j.f, 1)[0]; R = j.f - beta*j.r
    return R, j.f, j.r, beta


# ---- block bootstrap ----
def _mbb(n, block, rng):
    nb = int(np.ceil(n/block)); starts = rng.integers(0, max(1, n-block+1), size=nb)
    return np.concatenate([np.arange(s, s+block) for s in starts])[:n]

def boot(fn, n, block=60, nboot=2000, seed=0):
    rng = np.random.default_rng(seed); v = [fn(_mbb(n, block, rng)) for _ in range(nboot)]
    v = np.array([x for x in v if x is not None and np.isfinite(x)])
    return np.nanpercentile(v, [2.5, 97.5])

def tls_slope(x, y):
    xy = np.column_stack([x - x.mean(), y - y.mean()])
    w, V = np.linalg.eigh(np.cov(xy, rowvar=False)); pc = V[:, np.argmax(w)]
    return pc[1] / pc[0] if pc[0] != 0 else np.nan


def run_era(name, tag, lo, hi):
    b011_fams = set(pd.read_csv("data/b011_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)
    pgc_cert = pd.read_csv(f"data/{tag.lower()}_causality_cert.csv")
    pgc_fams = set(pgc_cert[pgc_cert.reliable].fam)
    F = b011_fams & pgc_fams
    pc, pm = f"data/daily_dvv_{tag}_Z_2to4.csv", f"data/daily_dvv_{tag}_MIRROR_2to4.csv"
    bc, bm = "data/daily_dvv_B011p90f40_Z_2to4.csv", "data/daily_dvv_B011p90f40_MIRROR_2to4.csv"
    for f in (pc, pm, bc, bm):
        if not os.path.exists(f):
            print(f"  [{name}] MISSING {f} -> mirror build not finished; skipping era"); return None

    P = corrected(pc, pm, F, lo, hi); B = corrected(bc, bm, F, lo, hi)
    if P is None or B is None:
        print(f"  [{name}] insufficient overlap"); return None
    Rp, fp, mp, betaP = P; Rb, fb, mb, betaB = B
    # documentation-only raw / deseasoned coda-vs-coda correlations
    jraw = pd.concat([cm(pc, F).rename("p"), cm(bc, F).rename("b")], axis=1).dropna()
    jraw = jraw[(jraw.index >= pd.Timestamp(lo)) & (jraw.index <= pd.Timestamp(hi))]
    raw_r = pearsonr(jraw.p, jraw.b)[0] if len(jraw) > 30 else np.nan
    jdes = pd.concat([fp.rename("p"), fb.rename("b")], axis=1).dropna()
    des_r = pearsonr(jdes.p, jdes.b)[0] if len(jdes) > 30 else np.nan

    # align verdict series on common days: R_PGC, R_B011, des(PGC mirror)
    J = pd.concat([Rp.rename("Rp"), Rb.rename("Rb"), mp.rename("mp")], axis=1).dropna()
    n = len(J)
    if n < 200:
        print(f"  [{name}] only {n} common days (<200) -> underpowered");
    x, y, c = J.Rp.values, J.Rb.values, J.mp.values
    S1 = spearmanr(x, y)[0]; N1 = spearmanr(c, y)[0]
    b_ols = np.polyfit(y, x, 1)[0]; b_tls = tls_slope(y, x)
    s1_ci = boot(lambda i: spearmanr(x[i], y[i])[0], n)
    diff_ci = boot(lambda i: spearmanr(x[i], y[i])[0] - spearmanr(c[i], y[i])[0], n)
    b_ci = boot(lambda i: np.polyfit(y[i], x[i], 1)[0], n)

    passes = (s1_ci[0] > 0) and (diff_ci[0] > 0) and (b_ci[0] > 0) and (0.3 <= b_ols <= 3)
    verdict = "PASS" if passes else "FAIL"
    if passes and S1 < 0.2: verdict = "INCONCLUSIVE (S1<0.2)"

    print(f"\n  ===== {name} era ({tag}) =====")
    print(f"  families F_e = PGCcert ∩ B011cert = {len(F)}   common days = {n}")
    print(f"  beta_PGC={betaP:.2f}  beta_B011={betaB:.2f}")
    print(f"  raw coda-vs-coda r      = {raw_r:.3f}   (confounded by shared noise field — documentation only)")
    print(f"  deseasoned coda r       = {des_r:.3f}   (still shares non-seasonal noise)")
    print(f"  S1 corrected-resid r    = {S1:.3f}   95% CI [{s1_ci[0]:.3f}, {s1_ci[1]:.3f}]")
    print(f"  N1 mirror-null r        = {N1:.3f}   (PGC mirror vs B011 corrected)")
    print(f"  S1 - N1 (paired)        = {S1-N1:.3f}  95% CI [{diff_ci[0]:.3f}, {diff_ci[1]:.3f}]")
    print(f"  slope b (OLS/TLS)       = {b_ols:.2f} / {b_tls:.2f}   95% CI [{b_ci[0]:.2f}, {b_ci[1]:.2f}]")
    print(f"  --> {name}: {verdict}")

    return dict(era=name, tag=tag, nfam=len(F), ndays=int(n), betaP=float(betaP), betaB=float(betaB),
                raw_r=float(raw_r), des_r=float(des_r), S1=float(S1), S1_ci=[float(s1_ci[0]), float(s1_ci[1])],
                N1=float(N1), diff=float(S1-N1), diff_ci=[float(diff_ci[0]), float(diff_ci[1])],
                b_ols=float(b_ols), b_tls=float(b_tls), b_ci=[float(b_ci[0]), float(b_ci[1])], verdict=verdict)


def per_family(name, tag, lo, hi, min_days=400, n_mismatch=20, seed=1):
    """Matched (fam i PGC vs fam i B011) vs mismatched (fam i PGC vs random fam j!=i B011) Spearman r."""
    b011_fams = set(pd.read_csv("data/b011_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)
    pgc_cert = pd.read_csv(f"data/{tag.lower()}_causality_cert.csv")
    F = sorted(b011_fams & set(pgc_cert[pgc_cert.reliable].fam))
    dp = pd.read_csv(f"data/daily_dvv_{tag}_Z_2to4.csv"); dp = dp[(dp.cc_max >= 0.6)]; dp["date"] = pd.to_datetime(dp.date)
    db = pd.read_csv("data/daily_dvv_B011p90f40_Z_2to4.csv"); db = db[(db.cc_max >= 0.6)]; db["date"] = pd.to_datetime(db.date)
    dp = dp[(dp.date >= pd.Timestamp(lo)) & (dp.date <= pd.Timestamp(hi))]
    db = db[(db.date >= pd.Timestamp(lo)) & (db.date <= pd.Timestamp(hi))]
    Pg = {f: des(g.set_index("date").dvv.sort_index()) for f, g in dp[dp.patch.isin(F)].groupby("patch")}
    Bb = {f: des(g.set_index("date").dvv.sort_index()) for f, g in db[db.patch.isin(F)].groupby("patch")}
    rng = np.random.default_rng(seed); matched, mismatched = [], []
    Bkeys = [f for f in F if f in Bb]
    for f in F:
        if f not in Pg or f not in Bb: continue
        j = pd.concat([Pg[f].rename("p"), Bb[f].rename("b")], axis=1).dropna()
        if len(j) < min_days: continue
        matched.append(spearmanr(j.p, j.b)[0])
        others = [g for g in Bkeys if g != f]
        for g in rng.choice(others, size=min(n_mismatch, len(others)), replace=False):
            jm = pd.concat([Pg[f].rename("p"), Bb[g].rename("b")], axis=1).dropna()
            if len(jm) >= min_days: mismatched.append(spearmanr(jm.p, jm.b)[0])
    matched, mismatched = np.array(matched), np.array(mismatched)
    if len(matched) < 5:
        print(f"  [{name}] per-family: only {len(matched)} families with >= {min_days} days -> skip"); return None
    delta = np.median(matched) - np.median(mismatched)
    # permutation p: pool and reshuffle labels
    pool = np.concatenate([matched, mismatched]); nm = len(matched)
    perm = np.array([np.median(p[:nm]) - np.median(p[nm:]) for p in [rng.permutation(pool) for _ in range(2000)]])
    pval = (np.sum(perm >= delta) + 1) / (len(perm) + 1)
    print(f"  [{name}] per-family: matched med r={np.median(matched):.3f} (n={len(matched)}), "
          f"mismatched med r={np.median(mismatched):.3f} (n={len(mismatched)}), Delta={delta:.3f}, p={pval:.4f}")
    return dict(era=name, n_matched=len(matched), matched_med=float(np.median(matched)),
                mismatched_med=float(np.median(mismatched)), delta=float(delta), pval=float(pval))


def seam_test():
    """B011 corrected-residual step at 2017-08-15 vs a random-date null (is there a REAL site velocity step?)."""
    fams = set(pd.read_csv("data/b011_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)
    R = corrected("data/daily_dvv_B011p90f40_Z_2to4.csv", "data/daily_dvv_B011p90f40_MIRROR_2to4.csv", fams)
    if R is None: print("  seam: B011 corrected residual unavailable"); return None
    R = R[0]
    def step(t0):
        pre = R[(R.index >= t0 - pd.Timedelta("90D")) & (R.index < t0)]
        post = R[(R.index >= t0) & (R.index < t0 + pd.Timedelta("90D"))]
        return (post.mean() - pre.mean(), len(pre), len(post)) if (len(pre) > 20 and len(post) > 20) else (np.nan, len(pre), len(post))
    obs, npre, npost = step(SEAM)
    print(f"\n  ===== B011 seam step @ {SEAM.date()} =====")
    print(f"  corrected-residual R: {len(R)} days total; window at seam has pre={npre}, post={npost}")
    if not np.isfinite(obs):
        # NaN obs would collapse the null p-value to its floor -> report UNMEASURABLE, never a spurious step
        print(f"  --> UNMEASURABLE: too few corrected-residual days in the +-90d seam window (need >20 each).")
        print(f"      Retry on the deseasoned coda (not mirror-corrected) to get a step estimate with a valid null...")
        Rc = des(cm("data/daily_dvv_B011p90f40_Z_2to4.csv", fams))
        def step2(t0):
            pre = Rc[(Rc.index >= t0 - pd.Timedelta("90D")) & (Rc.index < t0)]
            post = Rc[(Rc.index >= t0) & (Rc.index < t0 + pd.Timedelta("90D"))]
            return (post.mean() - pre.mean()) if (len(pre) > 20 and len(post) > 20) else np.nan
        obs = step2(SEAM)
        lo, hi = Rc.index.min() + pd.Timedelta("180D"), Rc.index.max() - pd.Timedelta("180D")
        rng = np.random.default_rng(2); cand = pd.date_range(lo, hi)
        null = np.array([step2(d) for d in rng.choice(cand, size=min(200, len(cand)), replace=False)])
        null = null[np.isfinite(null)]
        src = "deseasoned coda (mirror-corrected residual too sparse at seam)"
    else:
        lo, hi = R.index.min() + pd.Timedelta("180D"), R.index.max() - pd.Timedelta("180D")
        rng = np.random.default_rng(2); cand = pd.date_range(lo, hi)
        null = np.array([step(d)[0] for d in rng.choice(cand, size=min(200, len(cand)), replace=False)])
        null = null[np.isfinite(null)]
        src = "mirror-corrected residual"
    if not np.isfinite(obs) or len(null) < 20:
        print(f"  --> UNMEASURABLE even on deseasoned coda (insufficient data at seam); cannot test the Aug-2017 step.");
        return dict(seam_step=None, pval=None, n_null=int(len(null)), note="unmeasurable")
    pval = (np.sum(np.abs(null) >= abs(obs)) + 1) / (len(null) + 1)
    is_corrected = src.startswith("mirror")
    print(f"  step = {obs:+.4f}%  (post90d - pre90d, {src})   null p = {pval:.3f} (N={len(null)} random dates)")
    if not is_corrected:
        # fallback is mirror-UNCORRECTED -> can be the shared-noise artifact; a low p here does NOT establish a real step
        note = "INCONCLUSIVE: estimate is mirror-UNCORRECTED (artifact-prone); no validated step — treat Aug-2017 site step as untested"
    elif pval < 0.1:
        note = "REAL site step near Aug-2017 (flag it)"
    else:
        note = "no detectable real velocity step at the site"
    print(f"  --> {note}")
    return dict(seam_step=float(obs), pval=float(pval), n_null=int(len(null)), source=src, verdict=note)


if __name__ == "__main__":
    print("=" * 70); print("B011 (borehole) vs PGC (broadband) co-located validation"); print("=" * 70)
    era_res = [run_era(*e) for e in ERAS]
    fam_res = [per_family(*e) for e in ERAS]
    seam = seam_test()
    era_ok = [r for r in era_res if r]
    program = "VALIDATED" if (len(era_ok) == 2 and all(r["verdict"] == "PASS" for r in era_ok)) else \
              "PARTIAL/PENDING" if era_ok else "NOT YET (mirror build incomplete)"
    print("\n" + "=" * 70)
    print(f"BROADBAND PROGRAM: {program}  (requires BOTH eras PASS)")
    print("Cross-seam dv/v STEP at PGC is UNMEASURABLE by construction (per-era references); B011 is the sole site authority.")
    print("=" * 70)
    out = dict(eras=era_res, per_family=fam_res, seam=seam, program=program)
    with open("data/b011_pgc_validation.json", "w") as f: json.dump(out, f, indent=2, default=str)
    print("-> data/b011_pgc_validation.json")
