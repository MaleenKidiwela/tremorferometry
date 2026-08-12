#!/usr/bin/env python
"""Can we identify reliable families WITHOUT full reverse densify? Evidence for Merlin.
Uses existing full fwd + rev daily-stack npz. Ground truth = full-reverse fwd/rev ratio>1.5.
Tests: (A) spread of reversed coda floor across families (constant=global scalar, cheap);
       (B) sampled reverse: K random days -> ratio -> agreement with full-reverse reliable set;
       (C) causality-only (fwd coda/mirror, NO reverse) -> agreement + correlation with full ratio.
"""
import sys, numpy as np, pandas as pd

TAG = sys.argv[1] if len(sys.argv) > 1 else "B926p90f40"
rng = np.random.RandomState(0)


def load(rev=False):
    d = np.load(f"data/long_window_daily_{TAG}{'rev' if rev else ''}_Z.npz", allow_pickle=True)
    return d["t"], d["patches"], d["dates"], d["n_det"].astype(float), d["stacks"]


t, p, dt, nd, S = load(False)
_, pr, dtr, ndr, Sr = load(True)
coda = (t >= 2) & (t <= 4); mir = (t >= -2) & (t < 0)
rms = lambda x: np.sqrt(np.mean(x**2, axis=-1))
fams = pd.unique(p)


def grand(P, D, ND, ST, fam, days=None):
    m = P == fam
    if days is not None:
        keep = np.isin(D[m], days); rows = ST[m][keep]; w = ND[m][keep]
    else:
        rows = ST[m]; w = ND[m]
    return (rows * w[:, None]).sum(0) / w.sum() if w.sum() > 0 else None


# --- forward coda + mirror (causality), full reverse coda (ground truth) ---
rec = []
for fam in fams:
    gf = grand(p, dt, nd, S, fam); gr = grand(pr, dtr, ndr, Sr, fam)
    if gf is None or gr is None: continue
    fc = float(rms(gf[coda])); fm = float(rms(gf[mir])); rc = float(rms(gr[coda]))
    rec.append(dict(fam=fam, fwd_coda=fc, fwd_mir=fm, rev_coda=rc,
                    ratio=fc/rc if rc else np.nan, caus=fc/fm if fm else np.nan))
R = pd.DataFrame(rec).dropna()
truth = R.ratio > 1.5
print(f"=== {TAG}: {len(R)} families, reliable (full fwd/rev>1.5) = {truth.sum()} ===")

# (A) reversed floor spread
rc = R.rev_coda
print(f"\n(A) reversed coda floor across families: median {rc.median():.4f}, "
      f"IQR {rc.quantile(.25):.4f}-{rc.quantile(.75):.4f}, CoV {rc.std()/rc.mean():.2f}  "
      f"(low CoV => ~global scalar, no per-family reverse needed)")

# (B) sampled reverse: K random days
alldays = np.unique(dtr)
print(f"\n(B) sampled reverse ({len(alldays)} reversed days total):")
for K in [5, 10, 30, 100, 365]:
    if K > len(alldays): continue
    days = rng.choice(alldays, K, replace=False)
    rk = []
    for fam in R.fam:
        gr = grand(pr, dtr, ndr, Sr, fam, days=days)
        rk.append(rms(gr[coda]) if gr is not None else np.nan)
    ratio_k = R.fwd_coda.values / np.array(rk)
    pred = ratio_k > 1.5
    tp = (pred & truth.values).sum(); fp = (pred & ~truth.values).sum(); fn = (~pred & truth.values).sum()
    prec = tp/(tp+fp) if tp+fp else 0; rec_ = tp/(tp+fn) if tp+fn else 0
    print(f"   K={K:4d} days: reliable set precision {prec:.2f} recall {rec_:.2f}  "
          f"(agree {(pred==truth.values).mean()*100:.0f}%)")

# (C) causality-only (NO reverse)
print(f"\n(C) causality-only (fwd coda/mirror, ZERO reverse):")
print(f"   corr(causality, full ratio) = {R.caus.corr(R.ratio):.2f}")
for thr in [1.3, 1.5, 1.7, 2.0]:
    pred = R.caus > thr
    tp = (pred & truth).sum(); fp = (pred & ~truth).sum(); fn = (~pred & truth).sum()
    prec = tp/(tp+fp) if tp+fp else 0; rec_ = tp/(tp+fn) if tp+fn else 0
    print(f"   caus>{thr}: precision {prec:.2f} recall {rec_:.2f} agree {(pred==truth).mean()*100:.0f}% (n_pred {pred.sum()})")
