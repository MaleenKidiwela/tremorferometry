#!/usr/bin/env python
"""DECISIVE control (Merlin): is the dv/v real coda-borne signal, or a stacked-noise-field artifact?
Run the SAME stretch measurement on the REVERSED (noise-triggered, no real coda) daily stacks and
correlate the reversed common-mode with the forward common-mode. High corr => artifact.
Common-mode = per-date median over the certified families (cc_max>=0.6), 15-day centered rolling median.
Reports raw AND deseasoned (remove trend+annual+semiannual) correlation.
Acceptance (Z): deseasoned r(fwdZ,revZ) <= 0.2 -> Z payload real. >= 0.4 -> STOP, Z contaminated.
Assumes the rev dv/v CSVs already exist (daily_dvv_<TAG>rev_<C>_2to4.csv)."""
import sys, numpy as np, pandas as pd

def common_mode(csv, fams, cc_min=0.6):
    df = pd.read_csv(csv); df = df[df.patch.isin(fams) & (df.cc_max >= cc_min)]
    df["date"] = pd.to_datetime(df.date)
    med = df.groupby("date").dvv.median(); n = df.groupby("date").patch.nunique()
    med = med[n >= 3]
    return med.rolling(15, center=True, min_periods=5).median()

def deseason(s):
    s = s.dropna(); t = (s.index - s.index[0]).days.values.astype(float)
    yr = 365.25
    X = np.column_stack([np.ones_like(t), t,
                         np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr),
                         np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
    beta, *_ = np.linalg.lstsq(X, s.values, rcond=None)
    return pd.Series(s.values - X @ beta, index=s.index), 1 - np.var(s.values - X@beta)/np.var(s.values)

def corr(a, b):
    j = pd.concat([a.rename("a"), b.rename("b")], axis=1).dropna()
    return (float(j.a.corr(j.b)), len(j)) if len(j) > 30 else (np.nan, len(j))

for STA in ["B926", "B011"]:
    TAG = f"{STA}p90f40"; s = STA.lower()
    zfams = set(pd.read_csv(f"data/{s}_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)
    try:
        h2fams = set(pd.read_csv(f"data/{s}_h2_pass_families.csv").fam)
    except Exception:
        h2fams = set()
    print(f"\n===== {STA} =====")
    for C, fams in [("Z", zfams), ("H2", h2fams), ("H1", zfams)]:
        if not fams: continue
        try:
            fwd = common_mode(f"data/daily_dvv_{TAG}_{C}_2to4.csv", fams)
            rev = common_mode(f"data/daily_dvv_{TAG}rev_{C}_2to4.csv", fams)
        except FileNotFoundError as e:
            print(f"  {C}: missing {e.filename}"); continue
        r_raw, n = corr(fwd, rev)
        fd, _ = deseason(fwd); rd, _ = deseason(rev)
        r_des, _ = corr(fd, rd)
        flag = ""
        if C == "Z":
            flag = "  <-- Z REAL (<=0.2)" if (r_des <= 0.2) else ("  <-- ** STOP: Z CONTAMINATED (>=0.4) **" if r_des >= 0.4 else "  <-- Z ambiguous (0.2-0.4)")
        print(f"  {C}: r(fwd,rev) raw {r_raw:+.2f}  deseasoned {r_des:+.2f}  (n={n}){flag}")
