#!/usr/bin/env python
"""Merlin Step 1: purity-stratified horizontal grand-stack test (B926) — GO/NO-GO for 3-comp build.
Zero GPU. Uses existing fwd+rev daily-stack npz. Restrict to the 110 certified real-coda families
(b926_fwd_vs_rev_coda ratio>1.5) AND PNSN tremor-episode days near the station. Per family+channel:
  causality = RMS(coda 2-4 s) / RMS(mirror -2..0 s)   [real scattered coda causal; ringing symmetric]
  fwd/rev   = RMS(fwd coda 2-4 s) / RMS(rev coda 2-4 s) [real > noise-match floor]
GATE (Merlin): >=10 families with causality>1.5 AND fwd/rev>1.5 on a horizontal -> per-event H signal
exists -> build 3-comp. Else horizontals dead -> honest Z-only. Z is the positive control.
"""
import numpy as np, pandas as pd

STA_LAT, STA_LON = 48.82, -124.131
EPI_BOX = 0.4          # deg half-width for tremor-near-station
EPI_MIN = 3            # >= this many PNSN tremor dets in box that day => episode day
CERT = "data/b926_fwd_vs_rev_coda.csv"

# --- episode days from PNSN master catalog (tremor near B926) ---
cat = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time","lat","lon"])
cat = cat[(cat.lat.between(STA_LAT-EPI_BOX, STA_LAT+EPI_BOX)) &
          (cat.lon.between(STA_LON-EPI_BOX, STA_LON+EPI_BOX))]
cat["d"] = pd.to_datetime(cat.time).dt.strftime("%Y-%m-%d")
epi = cat.groupby("d").size()
epi_days = set(epi[epi >= EPI_MIN].index)
print(f"episode days near B926 (>= {EPI_MIN} tremor/day in +-{EPI_BOX}deg): {len(epi_days)}")

cert = set(pd.read_csv(CERT).query("ratio > 1.5").fam)
print(f"certified families (ratio>1.5): {len(cert)}")


def grand(stacks, patches, dates, ndet, fam, daymask_set):
    m = (patches == fam)
    if daymask_set is not None:
        dm = np.array([d in daymask_set for d in dates[m]])
        rows = stacks[m][dm]; w = ndet[m][dm].astype(float)
    else:
        rows = stacks[m]; w = ndet[m].astype(float)
    if w.sum() <= 0 or len(rows) == 0:
        return None, 0
    g = (rows * w[:, None]).sum(0) / w.sum()
    return g, int(w.sum())


def rms(x):
    return float(np.sqrt(np.mean(x**2)))


results = {}
for ch in ["Z", "H1", "H2"]:
    f = np.load(f"data/long_window_daily_B926p90f40_{ch}.npz", allow_pickle=True)
    r = np.load(f"data/long_window_daily_B926p90f40rev_{ch}.npz", allow_pickle=True)
    t = f["t"]
    coda = (t >= 2) & (t <= 4)
    mirror = (t >= -2) & (t < 0)
    fp, fd, fn, fs_ = f["patches"], f["dates"], f["n_det"], f["stacks"]
    rp, rd, rn, rs_ = r["patches"], r["dates"], r["n_det"], r["stacks"]
    rows = []
    for fam in cert:
        gf, nf = grand(fs_, fp, fd, fn, fam, epi_days)
        gr, nr = grand(rs_, rp, rd, rn, fam, epi_days)
        if gf is None or nf < 200:      # need enough episode-day detections to stack
            continue
        caus = rms(gf[coda]) / (rms(gf[mirror]) + 1e-30)
        fr = rms(gf[coda]) / (rms(gr[coda]) + 1e-30) if gr is not None else np.nan
        rows.append(dict(fam=fam, n_epi_det=nf, caus=caus, fwd_rev=fr))
    R = pd.DataFrame(rows)
    passed = R[(R.caus > 1.5) & (R.fwd_rev > 1.5)]
    results[ch] = R
    print(f"\n=== {ch} ===  families with enough episode-day dets: {len(R)}")
    if len(R):
        print(f"  causality  median {R.caus.median():.2f}  (>1.5: {(R.caus>1.5).sum()})")
        print(f"  fwd/rev    median {R.fwd_rev.median():.2f}  (>1.5: {(R.fwd_rev>1.5).sum()})")
        print(f"  PASS BOTH (caus>1.5 & fwd/rev>1.5): {len(passed)} families")
    del f, r, fs_, rs_

print("\n" + "="*60)
h_pass = max(len(results["H1"][(results["H1"].caus>1.5)&(results["H1"].fwd_rev>1.5)]),
             len(results["H2"][(results["H2"].caus>1.5)&(results["H2"].fwd_rev>1.5)]))
z_pass = len(results["Z"][(results["Z"].caus>1.5)&(results["Z"].fwd_rev>1.5)])
print(f"Z control passes: {z_pass} families (sanity: should be many)")
print(f"BEST horizontal passes: {h_pass} families")
print("VERDICT:", "PROCEED to 3-comp build (>=10 H families)" if h_pass >= 10
      else "HORIZONTALS DEAD at B926 -> honest Z-only (report gate stats)")
for ch in results:
    results[ch].to_csv(f"data/b926_purity_h_test_{ch}.csv", index=False)
print("wrote data/b926_purity_h_test_{Z,H1,H2}.csv")
