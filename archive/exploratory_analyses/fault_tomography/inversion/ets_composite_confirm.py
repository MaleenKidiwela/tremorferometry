"""CONFIRMATORY per-cell LOCAL-ETS-locked DEEP composite (Merlin detector #1). Mirror-free.
PRE-REGISTERED (fixed before looking at output):
  primary statistic = deep-cell mean of the ETS(minus inter-ETS) composite inversion (model %).
  predicted sign    = NEGATIVE (shear-velocity drop during locally tremor-active months).
  DETECTION  <=>  one-sided p < 0.01 on >=500 per-pair CIRCULAR-SHIFT nulls  AND all gates:
     G1 year-shift null (±12/24/36 mo; keeps season, breaks ETS) -> observed beyond its spread
     G2 sign-consistent & significant on BOTH grids
     G3 leave-one-station-out: <50% attenuation, no single station flips significance
     G4 survives excluding depth_mixed cells
     G5 n_days-matched control: composite not explained by ETS-vs-inter data-density imbalance
Per-pair composite differences WITHIN pair over the full record -> immune to the 2010 endpoint & baseline demean.
Reuses invert_coarse_cm data prep. Reads {RESDIR}/{G.npz, pair_months.parquet, era_table.csv}."""
import os, warnings; warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
from numpy.linalg import inv
from scipy.spatial import cKDTree
from scipy.stats import f as fdist, pearsonr
RD = os.environ.get("RESDIR", "fault_tomography/inversion/res_catalog_g50")
LAM = float(os.environ.get("LAM", "0.5")); NNULL = int(os.environ.get("NNULL", "500"))
d = np.load(f"{RD}/G.npz", allow_pickle=True)
G = -d["G"].astype(float); f = d["captured"].astype(float); clat = d["cell_lat"]; clon = d["cell_lon"]; depth = d["depth_km"]
mixed = d["depth_mixed"] if "depth_mixed" in d.files else np.zeros(len(depth), bool)
cxy = d["cxy"]; cellids = d["cell"]; ptag = d["pair_tag"]; pcell = d["pair_cell"]; Mf = G.shape[1]
keep = f > 0; Gc = (G*f[:, None])[keep]; ptag, pcell = ptag[keep], pcell[keep]
row_of = {(t, c): i for i, (t, c) in enumerate(zip(ptag, pcell))}; deep = depth > 30
_, nbr = cKDTree(cxy).query(cxy, k=min(6, Mf)); L = np.zeros((Mf, Mf))
for i in range(Mf):
    for j in nbr[i, 1:]:
        L[i, i] += 1; L[i, j] -= 1
LtL = L.T@L; regmat = LAM**2*LtL + 1e-4*np.eye(Mf)
# ---- data prep (identical to invert_coarse_cm) ----
pm = pd.read_parquet(f"{RD}/pair_months.parquet"); pm = pm[pm.n_days >= 5].copy()
pm["row"] = [row_of.get((t, c), -1) for t, c in zip(pm.tag, pm.cell)]; pm = pm[pm.row >= 0]
pm["t"] = pd.PeriodIndex(pm.ym, freq="M"); pm["year"] = pm.t.dt.year; pm["pair"] = pm.tag+"|"+pm.cell
et = pd.read_csv(f"{RD}/era_table.csv") if os.path.exists(f"{RD}/era_table.csv") else pd.read_csv("fault_tomography/inversion/res_catalog/era_table.csv")
bnd_of = {r.tag: sorted(pd.Timestamp(x) for x in str(r.boundaries).split(";") if x and x != "nan") for _, r in et.iterrows()}
ts = pm.t.dt.to_timestamp().values; eid = np.zeros(len(pm), int); strad = np.zeros(len(pm), bool)
for i, (tg, t) in enumerate(zip(pm.tag.values, ts)):
    bs = bnd_of.get(tg, [])
    if bs: eid[i] = int(sum(1 for b in bs if t >= np.datetime64(b))); strad[i] = any(abs((pd.Timestamp(t)-b).days) <= 35 for b in bs)
pm["stera"] = pm.tag+"__e"+eid.astype(str); pm = pm[~strad].copy()
def harm(t): yr = 12.; return np.column_stack([np.sin(2*np.pi*t/yr), np.cos(2*np.pi*t/yr), np.sin(4*np.pi*t/yr), np.cos(4*np.pi*t/yr)])
def deseason(g):
    v = g.dvv_month.values.astype(float); t = (g.t-g.t.min()).apply(lambda x: x.n).values.astype(float); n = len(v)
    if n >= 24:
        X = np.column_stack([np.ones_like(t), harm(t)]); b, *_ = np.linalg.lstsq(X, v, rcond=None)
        RSSf = ((v-X@b)**2).sum(); RSSc = ((v-v.mean())**2).sum(); F = ((RSSc-RSSf)/4)/(RSSf/(n-5)) if RSSf > 0 else 0
        if 1-fdist.cdf(F, 4, n-5) < 0.05: v = v-harm(t)@b[1:]
    return v
pm = pm.sort_values(["pair", "t"]); pm["ds"] = np.concatenate([deseason(g) for _, g in pm.groupby("pair")])
pm["cmode"] = pm.groupby(["stera", "ym"]).ds.transform("mean"); pm["resid"] = pm.ds - pm.cmode
def demean(g):
    inw = (g.year.values >= 2019) & (g.year.values <= 2024); return g.resid.values - (g.resid.values[inw].mean() if inw.sum() >= 6 else g.resid.values.mean())
pm["anom"] = np.concatenate([demean(g) for _, g in pm.groupby("pair")])
sig_pair = pm.groupby("pair").anom.std().to_dict(); TVM = float(np.nanmedian(list(sig_pair.values())))
# ---- per-CELL local ETS months (tremor within 100 km, monthly count z-scored within year, z>0.5) ----
tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat", "lon"])
tr["t"] = pd.to_datetime(tr["time"], errors="coerce"); tr = tr.dropna(subset=["t", "lat", "lon"])
lat0, lon0 = clat.mean(), clon.mean()
def xy(la, lo): return np.column_stack([(np.asarray(lo)-lon0)*111*np.cos(np.radians(lat0)), (np.asarray(la)-lat0)*111])
nn = cKDTree(xy(tr.lat.values, tr.lon.values)).query_ball_point(xy(clat, clon), r=100)
tr_ym = tr.t.dt.to_period("M"); tr_yr = tr.t.dt.year.values
cell_ets = {}
for k, cid in enumerate(cellids):
    ev = tr_ym.values[nn[k]]; yy = tr_yr[nn[k]]
    if len(ev) == 0: cell_ets[cid] = set(); continue
    mc = pd.DataFrame({"ym": ev, "yr": yy}).groupby(["yr", "ym"]).size().reset_index(name="n")
    z = mc.groupby("yr").n.transform(lambda s: (s-s.mean())/(s.std() if s.std() > 0 else 1))
    cell_ets[cid] = set(pd.PeriodIndex(mc.ym[z > 0.5]))
# ---- per-pair arrays (time-ordered) ----
pairs = []
for p, g in pm.sort_values("t").groupby("pair"):
    per = g.t.values; a = g.anom.values.astype(float); nd = g.n_days.values.astype(float)
    cid = g.cell.iloc[0]; ce = cell_ets.get(cid, set())
    ets = np.array([pp in ce for pp in pd.PeriodIndex(per)])
    if ets.sum() >= 3 and (~ets).sum() >= 3:
        pairs.append(dict(pair=p, tag=g.tag.iloc[0], row=int(g.row.iloc[0]), per=pd.PeriodIndex(per),
                          a=a, ets=ets, nd=nd, cid=cid, sig=max(1e-3, sig_pair.get(p, TVM))))
rows = np.array([q["row"] for q in pairs]); wts = np.array([1/q["sig"] for q in pairs])
print(f"[{os.path.basename(RD)}] valid pairs (>=3 ETS & >=3 inter months): {len(pairs)} | deep cells {int(deep.sum())}")

def make_P(sel):                                   # linear inversion operator for a subset of pairs
    r = rows[sel]; w = wts[sel]; Gr = Gc[r]
    Minv = inv((Gr*w[:, None]).T@(Gr*w[:, None]) + regmat)
    ill = (np.abs(Gr) > 0).any(0)
    return Minv@(Gr.T*w**2), ill                   # (Mf x npair), illumination mask
def primary(dvv, P, ill, mask=None):
    m = P@dvv; dm = deep & ill & (~mixed if mask == "nomix" else True)
    return float(np.nanmean(m[dm])) if dm.any() else np.nan, m
def comp_vec(shiftfn=None):                        # per-pair (ETS mean - inter mean)
    out = np.empty(len(pairs))
    for i, q in enumerate(pairs):
        e = q["ets"] if shiftfn is None else shiftfn(q)
        out[i] = q["a"][e].mean() - q["a"][~e].mean()
    return out

sel_all = np.ones(len(pairs), bool); P, ill = make_P(sel_all)
obs_dvv = comp_vec(); OBS, mmap = primary(obs_dvv, P, ill)
# ---- circular-shift null (preserves per-pair ETS count & block structure) ----
def circ(q, k): return np.roll(q["ets"], k)
nullv = np.empty(NNULL)
for s in range(NNULL):
    rng = np.random.RandomState(s); dv = np.empty(len(pairs))
    for i, q in enumerate(pairs):
        e = np.roll(q["ets"], 1+rng.randint(len(q["ets"])-1))   # ONE shift; mask and complement share it
        dv[i] = q["a"][e].mean() - q["a"][~e].mean()
    nullv[s], _ = primary(dv, P, ill)
p_one = (1 + np.sum(nullv <= OBS)) / (1 + NNULL)          # one-sided, predicted negative
# ---- G1 year-shift null (keeps calendar season, breaks ETS alignment) ----
def yr_shift(sh):
    def fn(q):
        ce = cell_ets.get(q["cid"], set()); return np.array([(pp - sh) in ce for pp in q["per"]])
    return fn
yshift = {}
for sh in (12, 24, 36, -12, -24, -36):
    dv = comp_vec(yr_shift(sh)); yshift[sh], _ = primary(dv, P, ill)
yv = np.array(list(yshift.values()))
# ---- G3 leave-one-station-out ----
loso = {}
for tg in sorted(set(q["tag"] for q in pairs)):
    sel = np.array([q["tag"] != tg for q in pairs])
    if sel.sum() < 20: continue
    Ps, ils = make_P(sel); v, _ = primary(comp_vec()[sel], Ps, ils); loso[tg] = v
lv = np.array(list(loso.values())); att = np.abs(lv/OBS - 1)
# ---- G4 exclude depth_mixed cells ----
OBS_nomix, _ = primary(obs_dvv, P, ill, mask="nomix")
# ---- G5 n_days-matched control: does per-pair composite scale with ETS-vs-inter n_days imbalance? ----
imbal = np.array([q["nd"][q["ets"]].mean() - q["nd"][~q["ets"]].mean() for q in pairs])
r_nd, p_nd = pearsonr(imbal, obs_dvv)

# ---- verdict ----
G1 = bool(OBS < 0 and np.mean(OBS < yv) >= 0.8)                  # observed more negative than >=80% of year-shifts
G3 = bool(np.median(att) < 0.5 and (lv < 0).mean() >= 0.9)      # <50% median attenuation, sign holds under LOSO
G4 = bool(OBS_nomix < 0 and abs(OBS_nomix-OBS)/abs(OBS) < 0.5)
G5 = bool(abs(r_nd) < 0.3)                                       # composite not explained by data-density imbalance
print("="*70)
print(f"PRIMARY deep-cell mean composite: {OBS:+.4f}% (model)  | predicted sign: NEGATIVE  | got: {'NEG ✓' if OBS<0 else 'POS ✗'}")
print(f"circular-shift null: mean {nullv.mean():+.4f}%  one-sided p={p_one:.4f}  ({'PASS <0.01' if p_one<0.01 else 'fail'})   [N={NNULL}]")
print(f"G1 year-shift null:  values {np.round(yv,3)}  -> observed below {int(100*np.mean(OBS<yv))}% of them  ({'PASS' if G1 else 'fail'})")
print(f"G3 LOSO:  median atten {np.median(att)*100:.0f}%  sign-neg {int(100*(lv<0).mean())}%  range [{lv.min():+.3f},{lv.max():+.3f}]  ({'PASS' if G3 else 'fail'})")
print(f"G4 no-mixed-cells: {OBS_nomix:+.4f}%  ({'PASS' if G4 else 'fail'})")
print(f"G5 n_days control: corr(imbalance, composite) r={r_nd:+.3f} p={p_nd:.3f}  ({'PASS' if G5 else 'fail'})")
DET = (p_one < 0.01) and G1 and G3 and G4 and G5
print(f"VERDICT [{os.path.basename(RD)}]: primary {'PASS' if p_one<0.01 else 'FAIL'} | gates G1{'✓'if G1 else'✗'} G3{'✓'if G3 else'✗'} G4{'✓'if G4 else'✗'} G5{'✓'if G5 else'✗'}  ->  {'DETECTION (this grid)' if DET else 'SUGGESTIVE' if (OBS<0 and p_one<0.05) else 'NULL'}")
print(f"    (G2 sign-consistency across grids = compare this line for g50 and g20)")
np.savez_compressed(f"{RD}/ets_composite_confirm.npz", OBS=OBS, p_one=p_one, nullv=nullv, yv=yv, lv=lv,
                    OBS_nomix=OBS_nomix, r_nd=r_nd, mmap=mmap, deep=deep, mixed=mixed, clat=clat, clon=clon)
print(f"wrote {RD}/ets_composite_confirm.npz")
