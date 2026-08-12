#!/usr/bin/env python
"""Mirror-window reference vs reversed-densify reference (Merlin step 2 acceptance).
For B926/B011: does the FREE mirror dv/v track the reversed (noise) dv/v, and does correcting with
mirror give the same residual as correcting with reverse? Accept: r(mirror-cm,rev-cm)>=0.8 AND
r(resid_mirror,resid_rev)>=0.9 with matching std -> mirror reference is free, reverse stays retired."""
import numpy as np, pandas as pd

def cm(csv, fams):
    d = pd.read_csv(csv); d = d[d.patch.isin(fams) & (d.cc_max >= 0.6)]; d["date"] = pd.to_datetime(d.date)
    m = d.groupby("date").dvv.median(); n = d.groupby("date").patch.nunique()
    return (m[n >= 3].rolling(15, center=True, min_periods=5).median()) * 100

def des(s):
    s = s.dropna(); t = (s.index - s.index[0]).days.values.astype(float); yr = 365.25
    X = np.column_stack([np.ones_like(t), t, np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr), np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
    b, *_ = np.linalg.lstsq(X, s.values, rcond=None); return pd.Series(s.values - X@b, index=s.index)

def rr(a, b):
    j = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna(); return j.a.corr(j.b), j

for S in ["B926", "B011"]:
    TAG = f"{S}p90f40"; s = S.lower()
    fams = set(pd.read_csv(f"data/{s}_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)
    fwd = des(cm(f"data/daily_dvv_{TAG}_Z_2to4.csv", fams))
    rev = des(cm(f"data/daily_dvv_{TAG}rev_Z_2to4.csv", fams))
    mir = des(cm(f"data/daily_dvv_{TAG}_MIRROR_2to4.csv", fams))
    r_mr, _ = rr(mir, rev)
    _, jm = rr(fwd, mir); bm = np.polyfit(jm.b, jm.a, 1)[0]; res_m = jm.a - bm*jm.b
    _, jv = rr(fwd, rev); bv = np.polyfit(jv.b, jv.a, 1)[0]; res_v = jv.a - bv*jv.b
    r_res, _ = rr(res_m, res_v)
    print(f"\n===== {S} =====")
    print(f"  mirror std {mir.std():.3f}%  rev std {rev.std():.3f}%")
    print(f"  r(mirror-cm, rev-cm) = {r_mr:.2f}   (accept >=0.8)")
    print(f"  residual std: via-mirror {res_m.std():.3f}%  via-rev {res_v.std():.3f}%")
    print(f"  r(resid_mirror, resid_rev) = {r_res:.2f}   (accept >=0.9)")
    ok = (r_mr >= 0.8) and (r_res >= 0.9)
    print(f"  -> {'PASS: mirror reference is FREE, keep reverse retired' if ok else 'FAIL: mirror != reverse, reverse needed for correction'}")
print("\nMIRROR_TEST_DONE")
