#!/usr/bin/env python
"""Lightweight Lin-coincidence validation of B926's picker-LFE families.
Non-circular ground-truth check: do the families' detection times coincide with
known Lin (2023) LFE origin times near B926? Compares to a random-offset null.
Single process, no waveforms — safe to run alongside the pipeline."""
import numpy as np, pandas as pd

STA_LAT, STA_LON = 48.8202, -124.1312
WIN = 15.0          # +-s coincidence window
OFFSET_LO, OFFSET_HI = 2.0, 20.0   # B926 envelope-peak ~ Lin_OT + S-traveltime (2-20 s)

print("loading Lin catalog...", flush=True)
lin = pd.read_csv("data/raw_lfe/lin2023_lfe.csv", usecols=["OT", "lat", "lon"])
dkm = np.sqrt(((lin.lat - STA_LAT)*111.0)**2 +
              ((lin.lon - STA_LON)*111.0*np.cos(np.radians(STA_LAT)))**2)
lin = lin[dkm <= 80].copy()
lin_ot = np.sort(pd.to_datetime(lin.OT).astype("int64").values / 1e9)  # epoch s
print(f"  Lin LFEs within 80 km of B926: {len(lin_ot):,}")

mem = pd.read_parquet("data/b926_disc_p70_2010_2026_m3.members.parquet")
pk = pd.read_csv("data/family_picker_b926_p70_m3.csv")
lfe_fams = set(pk[pk.pred == "LFE"].fam)
mem = mem[mem.family_id.isin(lfe_fams)].copy()
mem["ep"] = pd.to_datetime(mem.time, utc=True).astype("int64").values / 1e9
print(f"  530-set member detections: {len(mem):,}")


def coinc_frac(det_ep, shift=0.0):
    """fraction of detections whose (t - offset) lands within WIN of a Lin OT,
    scanning the expected S-traveltime offset band."""
    t = det_ep + shift
    hit = np.zeros(len(t), bool)
    for off in (OFFSET_LO, (OFFSET_LO+OFFSET_HI)/2, OFFSET_HI):
        idx = np.searchsorted(lin_ot, t - off)
        for cand in (idx-1, idx):
            cand = np.clip(cand, 0, len(lin_ot)-1)
            hit |= np.abs(lin_ot[cand] - (t - off)) <= WIN
    return hit.mean()

rng_shift = 86400 * 37.0   # 37-day shift -> destroys real coincidence, keeps tremor rate
rows = []
for fam, g in mem.groupby("family_id"):
    ep = g.ep.values
    rows.append((fam, len(ep), coinc_frac(ep), coinc_frac(ep, rng_shift)))
d = pd.DataFrame(rows, columns=["fam", "n", "coinc", "null"])
d["excess"] = d.coinc - d.null
print()
print(f"=== B926 530-family Lin coincidence ===")
print(f"  families evaluated: {len(d)}")
print(f"  median real coincidence: {d.coinc.median():.2f} | median null: {d.null.median():.2f}")
print(f"  families with real >> null (excess>0.15): {int((d.excess>0.15).sum())} ({100*(d.excess>0.15).mean():.0f}%)")
print(f"  families with coinc>0.5: {int((d.coinc>0.5).sum())}")
print(f"  pooled: {d.coinc.mean():.3f} real vs {d.null.mean():.3f} null (chance)")
d.to_csv("data/lin_coinc_b926_byfamily.csv", index=False)
print("  wrote data/lin_coinc_b926_byfamily.csv")
