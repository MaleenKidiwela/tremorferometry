#!/usr/bin/env python
"""Re-plot gated 3-comp dv/v using the FROZEN certification (no re-densify/re-certify).
Reads the existing data/<sta>_fwd_vs_rev_coda.csv (Z reliable, ratio>1.5) + data/<sta>_h2_pass_families.csv
and gates the (possibly refreshed) daily_dvv_<TAG>_{Z,H2}_2to4.csv. For gap-fills where the reliable set
is already decided. Usage: replot_gated_dvv.py <STA> <TAG>"""
import sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

STA, TAG = sys.argv[1], sys.argv[2]; sta = STA.lower()
zcert = set(pd.read_csv(f"data/{sta}_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)
h2c = set(pd.read_csv(f"data/{sta}_h2_pass_families.csv").fam) if __import__("os").path.exists(f"data/{sta}_h2_pass_families.csv") else set()


def cert_med(csv, fams, cc_min=0.6):
    df = pd.read_csv(csv); df = df[df.patch.isin(fams) & (df.cc_max >= cc_min)]
    df["date"] = pd.to_datetime(df.date)
    med = df.groupby("date").dvv.median(); n = df.groupby("date").patch.nunique()
    med = med[n >= 3]
    return med.rolling(15, center=True, min_periods=5).median()


zmed = cert_med(f"data/daily_dvv_{TAG}_Z_2to4.csv", zcert)
fig, ax = plt.subplots(figsize=(13, 4.5))
ax.plot(zmed.index, zmed.values*100, lw=1.4, color="#1a1a2e", label=f"Z  ({len(zcert)} certified families)")
summ = dict(station=STA, tag=TAG, z_certified=len(zcert), z_std_pct=round(float(zmed.std()*100), 3), h2_pass=len(h2c))
if len(h2c) >= 10:
    h2med = cert_med(f"data/daily_dvv_{TAG}_H2_2to4.csv", h2c)
    ax.plot(h2med.index, h2med.values*100, lw=1.1, color="#c1440e", alpha=0.85, label=f"H2 ({len(h2c)} families, earned)")
    j = pd.concat([zmed.rename("Z"), h2med.rename("H2")], axis=1).dropna()
    summ["h2_std_pct"] = round(float(h2med.std()*100), 3)
    summ["z_h2_corr"] = round(float(j.Z.corr(j.H2)), 3) if len(j) > 30 else None
ax.axhline(0, color="0.6", lw=0.6); ax.set_ylim(-1, 1); ax.grid(alpha=0.2)
ax.set_ylabel("dv/v (%)"); ax.legend(loc="upper right", frameon=False)
ax.set_title(f"{STA} certified coda dv/v (2-4 s, origin-anchored) — Z{' + H2' if len(h2c)>=10 else ''}  [2017-filled]")
fig.tight_layout(); fig.savefig(f"figures/{sta}_certified_3comp_gated_dvv.png", dpi=130)
json.dump(summ, open(f"data/{sta}_3comp_summary.json", "w"), indent=2)
print(f"[{STA}] replot (frozen cert): Z {len(zcert)} std {summ['z_std_pct']}% | H2 {len(h2c)}"
      + (f" std {summ.get('h2_std_pct')}%" if len(h2c)>=10 else "") + f" -> figures/{sta}_certified_3comp_gated_dvv.png")
