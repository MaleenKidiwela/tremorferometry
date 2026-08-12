#!/usr/bin/env python
"""Merlin decisive tests (all CPU) for the H2 go/no-go at B926.
T1 base-rate: joint gate (caus>1.5 & fwd/rev>1.5) pass rate for CERTIFIED (110) vs NON-cert (190) H2.
             If equal -> gate uninformative about H2 -> Z-only.
T2 split-half within-vs-cross: odd/even episode-day half grand stacks A_f,B_f. within=cc(A_f,B_f);
             cross=cc(A_f,B_g), f!=g, families >0.1deg apart. Raw 2-4s WAVEFORM. Real: within>>cross.
             Calibrate on Z (real ref) and rev stacks (zero). Compare H2 gap to half the Z gap.
"""
import numpy as np, pandas as pd

STA_LAT, STA_LON, EPI_BOX, EPI_MIN = 48.82, -124.131, 0.4, 3
cat = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time","lat","lon"])
cat = cat[cat.lat.between(STA_LAT-EPI_BOX,STA_LAT+EPI_BOX) & cat.lon.between(STA_LON-EPI_BOX,STA_LON+EPI_BOX)]
cat["d"] = pd.to_datetime(cat.time).dt.strftime("%Y-%m-%d")
g = cat.groupby("d").size(); epi = set(g[g>=EPI_MIN].index)
cdf = pd.read_csv("data/b926_fwd_vs_rev_coda.csv")
cert = set(cdf.query("ratio>1.5").fam); allfam = set(cdf.fam)


def load(tag):
    d = np.load(f"data/long_window_daily_{tag}.npz", allow_pickle=True)
    return d["t"], d["patches"], d["dates"], d["n_det"].astype(float), d["stacks"]


def wmean(rows, w):
    return (rows*w[:,None]).sum(0)/w.sum() if w.sum() > 0 else None


def rms(x): return float(np.sqrt(np.mean(x**2)))


def famloc(f):
    p = f.split("_"); return float(p[0]), float(p[1])


# ---------- T1 base-rate + build half-stacks for T2 ----------
def analyze(ch):
    t, p, dt, nd, S = load(f"B926p90f40_{ch}")
    tr, pr, dtr, ndr, Sr = load(f"B926p90f40rev_{ch}")
    coda = (t>=2)&(t<=4); mir = (t>=-2)&(t<0)
    rows, halves = [], {}
    for fam in allfam:
        m = p==fam; dm = np.array([x in epi for x in dt[m]])
        rw = S[m][dm]; w = nd[m][dm]
        if w.sum() < 200: continue
        gf = wmean(rw, w)
        mr = pr==fam; dmr = np.array([x in epi for x in dtr[mr]])
        gr = wmean(Sr[mr][dmr], ndr[mr][dmr]) if (mr.any() and np.array([x in epi for x in dtr[mr]]).any()) else None
        caus = rms(gf[coda])/(rms(gf[mir])+1e-30)
        fr = rms(gf[coda])/(rms(gr[coda])+1e-30) if gr is not None else np.nan
        rows.append(dict(fam=fam, caus=caus, fwd_rev=fr, cert=fam in cert))
        # odd/even day half stacks for T2
        days = dt[m][dm]
        parity = np.array([int(x.replace("-",""))%2 for x in days])
        for par, key in [(1,"A"),(0,"B")]:
            sel = parity==par
            if nd[m][dm][sel].sum() > 0:
                halves.setdefault(fam, {})[key] = wmean(rw[sel], w[sel])
    R = pd.DataFrame(rows)
    return R, halves, coda


def base_rate(R, tag):
    j = (R.caus>1.5)&(R.fwd_rev>1.5)
    c = R[R.cert]; n = R[~R.cert]
    print(f"[{tag}] joint-gate pass:  CERT {j[R.cert].mean()*100:.0f}% ({j[R.cert].sum()}/{len(c)})"
          f"   NON-CERT {j[~R.cert].mean()*100:.0f}% ({j[~R.cert].sum()}/{len(n)})")
    return j[R.cert].mean(), j[~R.cert].mean()


def split_half(halves, coda, fams_use):
    def cc(a,b):
        a=a[coda]-a[coda].mean(); b=b[coda]-b[coda].mean()
        return float(a@b/((np.linalg.norm(a)*np.linalg.norm(b))+1e-30))
    fu = [f for f in fams_use if f in halves and "A" in halves[f] and "B" in halves[f]]
    within = np.array([cc(halves[f]["A"], halves[f]["B"]) for f in fu])
    cross = []
    for i,f in enumerate(fu):
        lf = famloc(f)
        for gg in fu[i+1:i+40]:   # sample cross pairs
            lg = famloc(gg)
            if abs(lf[0]-lg[0])>0.1 or abs(lf[1]-lg[1])>0.1:
                cross.append(cc(halves[f]["A"], halves[gg]["B"]))
    cross = np.array(cross)
    return np.median(within), np.median(cross), len(fu)


print("="*64)
res = {}
for ch in ["Z","H2","H1"]:
    R, halves, coda = analyze(ch)
    base_rate(R, ch)
    w, x, n = split_half(halves, coda, list(cert))       # split-half on certified set
    gap = w - x
    print(f"[{ch}] split-half (cert, n={n}): within {w:.2f}  cross {x:.2f}  GAP {gap:.2f}")
    res[ch] = dict(within=w, cross=x, gap=gap)
    del R, halves

print("="*64)
zg = res["Z"]["gap"]; h2g = res["H2"]["gap"]
print(f"Z gap {zg:.2f} (real ref) | H2 gap {h2g:.2f} | half-Z {zg/2:.2f}")
print("T2 verdict H2:", "REAL family-specific coda (gap>=half-Z)" if h2g >= zg/2
      else "ARTIFACT/shared-waveform (within~cross, gap<half-Z)")
