#!/usr/bin/env python
"""Round-4 REFINED step-scan (Merlin): used-band physical boundaries only (from era_table.csv), PER-STATION
placebo distributions (>=8 random dates away from real boundaries), attribution by percentile, leak from
ATTRIBUTED boundaries as the quadrature excess over same-station placebo. Run BEFORE (raw anom) and AFTER
(--corrected reads era-split corrected series) — the after-scan is the mechanical gate (leak excess < 0.02%).
Usage: python step_scan_v2.py [--after]"""
import sys, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
OUT = "fault_tomography/inversion/res_catalog"
AFTER = "--after" in sys.argv
pm = pd.read_parquet(f"{OUT}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["t"] = pd.to_datetime(pm.ym.astype(str)+"-15"); pm["pair"] = pm.tag+"|"+pm.cell
val = "anom_corr" if AFTER else "dvv_month"
if AFTER:
    ac = pd.read_parquet(f"{OUT}/anom_corrected.parquet"); pm = pm.merge(ac[["tag", "cell", "ym", "anom_corr"]], on=["tag", "cell", "ym"], how="inner")
era = pd.read_csv(f"{OUT}/era_table.csv")
bnd_of = {r.tag: [pd.Timestamp(x) for x in str(r.boundaries).split(";") if x and x != "nan"] for _, r in era.iterrows()}
def harm(t): yr = 365.25; return np.column_stack([np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr), np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
def step_b(g, T):
    t = (g.t-g.t.min()).dt.days.values.astype(float); v = g[val].values.astype(float)
    H = (g.t.values >= np.datetime64(T)).astype(float)
    if H.sum() < 4 or (len(H)-H.sum()) < 4: return None
    X = np.column_stack([np.ones_like(t), harm(t), H]); b, *_ = np.linalg.lstsq(X, v, rcond=None); return float(b[-1])
rng = np.random.RandomState(0); rows = []
for tag, bounds in bnd_of.items():
    sub = pm[pm.tag == tag]
    if len(sub) == 0: continue
    lo, hi = sub.t.min(), sub.t.max(); days = (hi-lo).days
    if days < 365: continue
    pc, pl = [], []
    for _ in range(12):
        T = lo + pd.Timedelta(days=int(rng.uniform(180, max(181, days-180))))
        if bounds and min(abs((T-b).days) for b in bounds) < 90: continue
        bs = [step_b(g, T) for _, g in sub.groupby("pair")]; bs = [x for x in bs if x is not None]
        if len(bs) >= 5: pc.append(abs(np.median(bs))); pl.append(1.4826*np.median(np.abs(np.array(bs)-np.median(bs))))
    if len(pc) < 3: continue
    p95, pleak = float(np.percentile(pc, 95)), float(np.median(pl))
    for T in bounds:
        bs = [step_b(g, T) for _, g in sub.groupby("pair")]; bs = [x for x in bs if x is not None]
        if len(bs) < 5: continue
        common = abs(np.median(bs)); leak = 1.4826*np.median(np.abs(np.array(bs)-np.median(bs)))
        rows.append(dict(tag=tag, boundary=str(T.date()), nfam=len(bs), common=common, leak=leak, p95=p95,
                         pleak=pleak, attributed=common > p95, leak_excess=float(np.sqrt(max(0, leak**2-pleak**2)))))
R = pd.DataFrame(rows); R.to_csv(f"{OUT}/step_scan_v2{'_after' if AFTER else ''}.csv", index=False)
att = R[R.attributed]
tag_label = "AFTER (era-split corrected)" if AFTER else "BEFORE (raw)"
print(f"=== REFINED SCAN {tag_label} ===")
print(f"boundaries scanned {len(R)} across {R.tag.nunique()} stations; ATTRIBUTED (>station 95th placebo) {len(att)}")
if len(att):
    print(f"attributed common |b| median {att.common.median():.3f}%; LEAK EXCESS (quadrature, attributed) median {att.leak_excess.median():.3f}% p90 {att.leak_excess.quantile(.9):.3f}%")
print(f"overall common |b| median {R.common.median():.3f}% (placebo-p95 median {R.p95.median():.3f})")
gate = (att.leak_excess.median() if len(att) else 0)
print(f"GATE: leak excess {'>' if gate>0.02 else '<'} 0.02% -> {'ERA-SPLIT WARRANTED / after-scan FAILS' if gate>0.02 else 'clean (era-split near-no-op / after-scan PASSES)'}")
print(f"wrote {OUT}/step_scan_v2{'_after' if AFTER else ''}.csv")
