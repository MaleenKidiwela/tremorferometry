#!/usr/bin/env python
"""Merlin resonance guards on the H2 families that passed the purity gate (B926).
Decides real scattered coda vs event-excited resonance / Rayleigh-coupled noise BEFORE any GPU.

Guard A (station-artifact): pairwise cross-FAMILY |cc| of coda (2-4 s). Scattered coda is
  family-specific (different sources light different scatterers) -> LOW cross-family cc. Resonance/
  ringing is family-invariant -> HIGH. Compare H2 passers vs Z passers (Z = known-real reference).
Guard B (Rayleigh check, Merlin exp 2): corr(H1 grand stack, Hilbert(Z grand stack)) on FAIL
  families. High -> H1's anti-causal shape IS Z-selected microseism (90-deg Rayleigh Z-H coupling).
Guard C (narrowband): coda spectral peak/median ratio, H2 passers vs Z. Resonance = narrowband line.
"""
import numpy as np, pandas as pd
from scipy.signal import hilbert

STA_LAT, STA_LON, EPI_BOX, EPI_MIN = 48.82, -124.131, 0.4, 3
cat = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time","lat","lon"])
cat = cat[cat.lat.between(STA_LAT-EPI_BOX,STA_LAT+EPI_BOX) & cat.lon.between(STA_LON-EPI_BOX,STA_LON+EPI_BOX)]
cat["d"] = pd.to_datetime(cat.time).dt.strftime("%Y-%m-%d")
g = cat.groupby("d").size(); epi = set(g[g>=EPI_MIN].index)
cert = set(pd.read_csv("data/b926_fwd_vs_rev_coda.csv").query("ratio>1.5").fam)


def load_grands(ch):
    d = np.load(f"data/long_window_daily_B926p90f40_{ch}.npz", allow_pickle=True)
    t, p, dt, nd, S = d["t"], d["patches"], d["dates"], d["n_det"], d["stacks"]
    out = {}
    for fam in cert:
        m = p == fam
        dm = np.array([x in epi for x in dt[m]])
        rows = S[m][dm]; w = nd[m][dm].astype(float)
        if w.sum() < 200: continue
        out[fam] = (rows*w[:,None]).sum(0)/w.sum()
    return t, out

t, Zg = load_grands("Z")
_, H1g = load_grands("H1")
_, H2g = load_grands("H2")
coda = (t>=2)&(t<=4); fs=40.0

# passers from the gate test
def passers(ch):
    r = pd.read_csv(f"data/b926_purity_h_test_{ch}.csv")
    return list(r[(r.caus>1.5)&(r.fwd_rev>1.5)].fam)
H2p, Zp = passers("H2"), passers("Z")

def cross_family_cc(grands, fams):
    C = np.array([grands[f][coda] for f in fams if f in grands])
    C = C - C.mean(1, keepdims=True); C /= (np.linalg.norm(C,axis=1,keepdims=True)+1e-30)
    M = C @ C.T; iu = np.triu_indices(len(C),1)
    return np.abs(M[iu])

ccz = cross_family_cc(Zg, Zp); cch = cross_family_cc(H2g, H2p)
print("=== Guard A: cross-family coda |cc| (LOW=real family-specific, HIGH=resonance/artifact) ===")
print(f"  Z passers (n={len(Zp)}):  median |cc| = {np.median(ccz):.2f}  90th = {np.quantile(ccz,.9):.2f}")
print(f"  H2 passers (n={len(H2p)}): median |cc| = {np.median(cch):.2f}  90th = {np.quantile(cch,.9):.2f}")

# Guard B: H1 vs Hilbert(Z) on FAIL families
fail = [f for f in cert if f not in set(passers("H2")) and f in H1g and f in Zg]
rho = []
for f in fail:
    z = Zg[f][coda]; h = H1g[f][coda]
    zh = np.imag(hilbert(Zg[f]))[coda]   # 90-deg shifted Z
    z=z-z.mean(); h=h-h.mean(); zh=zh-zh.mean()
    r1 = abs(np.corrcoef(h,z)[0,1]); r2 = abs(np.corrcoef(h,zh)[0,1])
    rho.append((r1,r2))
rho=np.array(rho)
print("\n=== Guard B: H1 coda vs Z (Rayleigh Z-H coupling test, FAIL fams) ===")
print(f"  |corr(H1, Z)|      median {np.median(rho[:,0]):.2f}")
print(f"  |corr(H1, Hilb Z)| median {np.median(rho[:,1]):.2f}  (>Z-corr => 90deg-coupled microseism)")

# Guard C: narrowband peak/median of coda spectrum
def pkmed(grands, fams):
    v=[]
    for f in fams:
        if f not in grands: continue
        x = grands[f][coda]; X = np.abs(np.fft.rfft(x*np.hanning(len(x))))
        fr = np.fft.rfftfreq(len(x),1/fs); b=(fr>=2)&(fr<=8)
        if b.sum()>2: v.append(X[b].max()/(np.median(X[b])+1e-30))
    return np.array(v)
pz, ph = pkmed(Zg,Zp), pkmed(H2g,H2p)
print("\n=== Guard C: coda spectral peak/median (HIGH=narrowband resonance) ===")
print(f"  Z passers  median {np.median(pz):.2f}  95th {np.quantile(pz,.95):.2f}")
print(f"  H2 passers median {np.median(ph):.2f}  95th {np.quantile(ph,.95):.2f}")

print("\n" + "="*60)
verdict_A = "REAL (family-specific)" if np.median(cch) < 0.5 and np.median(cch) < np.median(ccz)+0.15 else "ARTIFACT (family-invariant coda)"
print(f"Guard A verdict for H2: {verdict_A}")
