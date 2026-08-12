#!/usr/bin/env python
"""General per-station 3-component finalization (CPU): certify Z (fwd-vs-rev), run the free H-gate,
produce gated Z + H2 dv/v and a figure. Encodes the Merlin-validated borehole method (see RESUME).
Usage: python finalize_3comp_dvv.py <STA> <TAG>   e.g. B011 B011p90f40
Needs: long_window_daily_<TAG>_{Z,H1,H2}.npz, long_window_daily_<TAG>rev_{Z,H1,H2}.npz,
       daily_dvv_<TAG>_{Z,H2}_2to4.csv, catalogs/pnsn_tremor_cascadia_full.csv.
       (Episode-day region is derived from family locations, NOT station coords — families can be offset.)
Writes: data/<sta>_fwd_vs_rev_coda.csv, data/<sta>_h2_pass_families.csv,
        data/<sta>_3comp_summary.json, figures/<sta>_certified_3comp_gated_dvv.png"""
import sys, json, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

STA, TAG = sys.argv[1], sys.argv[2]
sta = STA.lower()
EPI_HALF, EPI_MIN = 0.5, 3


def _famloc(f):
    p = str(f).split("_")
    try:
        return float(p[0]), float(p[1])
    except Exception:
        return np.nan, np.nan


# Episode-day region from the FAMILY locations (station-coord-independent; families can sit far from the
# sensor — e.g. B011 families are ~0.8deg west of the station, so a station-centered box misses them).
# Robust center = median family lat/lon +- 0.5deg (ETS episodes are broad/segment-wide, so this captures
# regional tremor days = high-purity days for the whole family set).
_pz0 = np.load(f"data/long_window_daily_{TAG}_Z.npz", allow_pickle=True)["patches"]
_ll = np.array([_famloc(f) for f in pd.unique(_pz0)])
_ll = _ll[~np.isnan(_ll).any(1)]
CLAT, CLON = float(np.median(_ll[:, 0])), float(np.median(_ll[:, 1]))
cat = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
cat = cat[cat.lat.between(CLAT-EPI_HALF, CLAT+EPI_HALF) & cat.lon.between(CLON-EPI_HALF, CLON+EPI_HALF)]
cat["d"] = pd.to_datetime(cat.time).dt.strftime("%Y-%m-%d")
g = cat.groupby("d").size(); epi = set(g[g >= EPI_MIN].index)
print(f"[{STA}] episode box: family centroid {CLAT:.2f},{CLON:.2f} +-{EPI_HALF} -> {len(epi)} episode days")


def load(ch, rev=False):
    p = f"data/long_window_daily_{TAG}{'rev' if rev else ''}_{ch}.npz"
    d = np.load(p, allow_pickle=True)
    return d["t"], d["patches"], d["dates"], d["n_det"].astype(float), d["stacks"]


def rms(x): return float(np.sqrt(np.mean(x ** 2)))


def grand(p, dt, nd, S, fam, epi_only):
    m = p == fam
    if epi_only:
        dm = np.array([x in epi for x in dt[m]]); rows = S[m][dm]; w = nd[m][dm]
    else:
        rows = S[m]; w = nd[m]
    return ((rows * w[:, None]).sum(0) / w.sum(), w.sum()) if w.sum() > 0 else (None, 0)


# ---- Z certification: fwd/rev grand coda ratio ----
tz, pz, dz, nz, Sz = load("Z"); _, pzr, dzr, nzr, Szr = load("Z", rev=True)
coda = (tz >= 2) & (tz <= 4); mir = (tz >= -2) & (tz < 0)
zc = []
for fam in pd.unique(pz):
    gf, _ = grand(pz, dz, nz, Sz, fam, False)
    gr, _ = grand(pzr, dzr, nzr, Szr, fam, False)
    if gf is None or gr is None or rms(gr[coda]) <= 0: continue
    zc.append(dict(fam=fam, fwd=rms(gf[coda]), rev=rms(gr[coda]), ratio=rms(gf[coda]) / rms(gr[coda])))
ZC = pd.DataFrame(zc); ZC.to_csv(f"data/{sta}_fwd_vs_rev_coda.csv", index=False)
zcert = set(ZC.query("ratio>1.5").fam)
print(f"[{STA}] Z certified (fwd/rev>1.5): {len(zcert)}/{len(ZC)}")

# ---- H2 gate (causality + fwd/rev on episode days), base-rate QC ----
def h_gate(ch):
    t, p, dt, nd, S = load(ch); _, pr, dtr, ndr, Sr = load(ch, rev=True)
    rows = []
    for fam in pd.unique(p):
        gf, n = grand(p, dt, nd, S, fam, True)
        if gf is None or n < 200: continue
        gr, _ = grand(pr, dtr, ndr, Sr, fam, True)
        caus = rms(gf[coda]) / (rms(gf[mir]) + 1e-30)
        fr = rms(gf[coda]) / (rms(gr[coda]) + 1e-30) if gr is not None else np.nan
        rows.append(dict(fam=fam, caus=caus, fwd_rev=fr, cert=fam in zcert))
    return pd.DataFrame(rows)

H = h_gate("H2")
Hpass = H[(H.caus > 1.5) & (H.fwd_rev > 1.5) & H.cert]
Hpass.fam.to_csv(f"data/{sta}_h2_pass_families.csv", index=False)
br_c = ((H.cert) & (H.caus > 1.5) & (H.fwd_rev > 1.5)).sum() / max(1, H.cert.sum())
br_n = ((~H.cert) & (H.caus > 1.5) & (H.fwd_rev > 1.5)).sum() / max(1, (~H.cert).sum())
print(f"[{STA}] H2 gate: {len(Hpass)} certified families pass | base-rate cert {br_c*100:.0f}% vs non-cert {br_n*100:.0f}%")
h2_real = len(Hpass) >= 10 and br_c > 2 * max(br_n, 0.02)

# ---- gated dv/v medians + figure ----
def cert_med(csv, fams, cc_min=0.6):
    df = pd.read_csv(csv); df = df[df.patch.isin(fams) & (df.cc_max >= cc_min)]
    df["date"] = pd.to_datetime(df.date)
    med = df.groupby("date").dvv.median(); n = df.groupby("date").patch.nunique()
    med = med[n >= 3]
    return med.rolling(15, center=True, min_periods=5).median()

zmed = cert_med(f"data/daily_dvv_{TAG}_Z_2to4.csv", zcert)
fig, ax = plt.subplots(figsize=(13, 4.5))
ax.plot(zmed.index, zmed.values * 100, lw=1.4, color="#1a1a2e", label=f"Z  ({len(zcert)} certified families)")
summ = dict(station=STA, tag=TAG, z_certified=len(zcert), z_std_pct=round(float(zmed.std()*100), 3),
            h2_pass=int(len(Hpass)), h2_real=bool(h2_real), br_cert=round(br_c, 3), br_noncert=round(br_n, 3))
if h2_real:
    h2med = cert_med(f"data/daily_dvv_{TAG}_H2_2to4.csv", set(Hpass.fam))
    ax.plot(h2med.index, h2med.values * 100, lw=1.1, color="#c1440e", alpha=0.85, label=f"H2 ({len(Hpass)} families, earned)")
    j = pd.concat([zmed.rename("Z"), h2med.rename("H2")], axis=1).dropna()
    summ["h2_std_pct"] = round(float(h2med.std()*100), 3)
    summ["z_h2_corr"] = round(float(j.Z.corr(j.H2)), 3) if len(j) > 30 else None
    print(f"[{STA}] H2 REAL -> H2 dv/v std {summ['h2_std_pct']}%, r(Z,H2)={summ['z_h2_corr']}")
else:
    print(f"[{STA}] H2 does NOT clear gate -> Z-only (H1 anti-causal at boreholes; honest negative)")
ax.axhline(0, color="0.6", lw=0.6); ax.set_ylim(-1, 1); ax.grid(alpha=0.2)
ax.set_ylabel("dv/v (%)"); ax.legend(loc="upper right", frameon=False)
ax.set_title(f"{STA} certified coda dv/v (2-4 s, origin-anchored) — Z{' + H2' if h2_real else ' only'}")
fig.tight_layout(); fig.savefig(f"figures/{sta}_certified_3comp_gated_dvv.png", dpi=130)
json.dump(summ, open(f"data/{sta}_3comp_summary.json", "w"), indent=2)
print(f"[{STA}] -> figures/{sta}_certified_3comp_gated_dvv.png ; data/{sta}_3comp_summary.json")
print("SUMMARY", json.dumps(summ))
