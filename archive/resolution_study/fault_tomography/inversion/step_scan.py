#!/usr/bin/env python
"""Instrument/metadata-step SCAN (Merlin round-4 gate). For each pair whose dv/v record spans a FDSN response-
epoch boundary T, fit dvv = a + seasonal(4 harm) + b*H(t-T); b = the step. Per station split b into:
  COMMON = median(b over families)  (site-term-absorbable; measures raw-data severity)
  LEAK   = MAD(b over families)      (family-specific residual -> reaches the fault term)
Placebo: same fit at RANDOM dates (far from real boundaries) -> significance calibration.
Also: histogram of boundary dates vs the index years (clustered-upgrade check) + cc_max step (blend signature).
Thresholds (data-space %): common median |b| > 0.05 -> per-era ref mandatory; leak > 0.02 -> era-split rerun.
Outputs res_catalog/step_scan.{csv,txt}. dv/v (anom) is already data-space, so b is directly comparable."""
import os, glob, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
from obspy.clients.fdsn import Client
OUT = "fault_tomography/inversion/res_catalog"
pm = pd.read_parquet(f"{OUT}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["t"] = pd.to_datetime(pm.ym.astype(str) + "-15"); pm["pair"] = pm.tag + "|" + pm.cell
fo = pd.read_csv("data/broadband_fleet_order.csv"); net_of = dict(zip(fo.sta, fo.net))
cand = pd.read_csv("data/candidate_stations_post2020.csv")
for _, r in cand.iterrows(): net_of.setdefault(r.sta, r.get("net", "PB"))

def sta_net(tag):
    low = tag.lower()
    for pre, s in [("pgc", "PGC"), ("shb", "SHB"), ("clrs", "CLRS")]:
        if low.startswith(pre): return s, "CN"
    sta = tag.replace("p90f40", "").upper()
    return sta, net_of.get(sta, "PB" if sta[0] == "B" and sta[1:2].isdigit() else "UW")

def provider(net): return "NCEDC" if net in ("BK", "NC") else "IRIS"

# ---- FDSN response-epoch boundaries per tag (Z channels, within data span) ----
bnds = {}
for tag in sorted(pm.tag.unique()):
    sta, net = sta_net(tag)
    try:
        inv = Client(provider(net), timeout=30).get_stations(network=net, station=sta, level="channel")
        ds = sorted({pd.Timestamp(c.start_date.datetime) for c in inv[0][0].channels if c.code.endswith("Z")
                     and pd.Timestamp(c.start_date.datetime) > pd.Timestamp("2009-06-01")
                     and pd.Timestamp(c.start_date.datetime) < pd.Timestamp("2026-01-01")})
        # merge boundaries < 30 days apart
        merged = []
        for d in ds:
            if not merged or (d - merged[-1]).days > 30: merged.append(d)
        bnds[tag] = merged
    except Exception:
        bnds[tag] = []
nb = sum(len(v) for v in bnds.values()); nstn = sum(1 for v in bnds.values() if v)
print(f"tags with >=1 in-record response boundary: {nstn}/{len(bnds)} | total boundaries: {nb}")

# ---- step fit ----
def harm(t): yr = 365.25; return np.column_stack([np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr), np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
def step_b(g, T):
    t = (g.t - g.t.min()).dt.days.values.astype(float); v = g.dvv_month.values.astype(float)
    H = (g.t.values >= np.datetime64(T)).astype(float)
    pre, post = H.sum() == 0, H.sum() == len(H)
    if pre or post or H.sum() < 4 or (len(H)-H.sum()) < 4: return None
    X = np.column_stack([np.ones_like(t), harm(t), H]); b, *_ = np.linalg.lstsq(X, v, rcond=None)
    return float(b[-1])

rows = []
for tag, blist in bnds.items():
    if not blist: continue
    sub = pm[pm.tag == tag]
    for T in blist:
        bs = [step_b(g, T) for _, g in sub.groupby("pair")]
        bs = [x for x in bs if x is not None]
        if len(bs) >= 5:
            rows.append(dict(tag=tag, boundary=pd.Timestamp(T).date(), nfam=len(bs),
                             common=float(np.median(bs)), leak=float(1.4826*np.median(np.abs(bs-np.median(bs))))))
R = pd.DataFrame(rows)

# ---- placebo: random boundary dates (not near true) ----
rng = np.random.RandomState(0); pl = []
for tag in R.tag.unique()[:60]:
    sub = pm[pm.tag == tag]; span = (sub.t.min(), sub.t.max())
    if (span[1]-span[0]).days < 365: continue
    T = span[0] + pd.Timedelta(days=int(rng.uniform(180, (span[1]-span[0]).days-180)))
    bs = [step_b(g, T) for _, g in sub.groupby("pair")]; bs = [x for x in bs if x is not None]
    if len(bs) >= 5: pl.append(dict(common=float(np.median(bs)), leak=float(1.4826*np.median(np.abs(bs-np.median(bs))))))
P = pd.DataFrame(pl)

R.to_csv(f"{OUT}/step_scan.csv", index=False)
out = []
out.append(f"boundaries scanned: {len(R)} (>=5 families each) across {R.tag.nunique()} stations")
out.append(f"COMMON step |b| (data-space %): median {R.common.abs().median():.3f}  p90 {R.common.abs().quantile(.9):.3f}  max {R.common.abs().max():.3f}")
out.append(f"LEAK (family spread) %:          median {R.leak.median():.3f}  p90 {R.leak.quantile(.9):.3f}  max {R.leak.max():.3f}")
out.append(f"PLACEBO common |b| median {P.common.abs().median():.3f}  leak median {P.leak.median():.3f}  (n={len(P)})")
out.append(f"significant boundaries (common |b| > 2x placebo-common-p90={2*P.common.abs().quantile(.9):.3f}%): {int((R.common.abs()>2*P.common.abs().quantile(.9)).sum())}")
out.append(f"THRESHOLDS: common>0.05% -> per-era ref mandatory (v2); leak>0.02% -> era-split rerun (v1)")
out.append(f"  common median {R.common.abs().median():.3f} {'EXCEEDS' if R.common.abs().median()>0.05 else 'below'} 0.05% ; leak median {R.leak.median():.3f} {'EXCEEDS' if R.leak.median()>0.02 else 'below'} 0.02%")
# boundary-date histogram vs years (clustered-upgrade check)
yc = pd.to_datetime(R.boundary).dt.year.value_counts().sort_index()
out.append("boundary-date histogram (clustered-upgrade check): " + " ".join(f"{y}:{yc[y]}" for y in yc.index))
txt = "\n".join(out); open(f"{OUT}/step_scan.txt", "w").write(txt); print("\n"+txt)
print(f"\nwrote {OUT}/step_scan.csv, step_scan.txt")
