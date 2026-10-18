#!/usr/bin/env python
"""Space-time connected slip patches in the SSE model (user 2026-09-09): thrust slip accumulated within the trailing
30-d window (t-30, t] (our dv/v frame; thrust-positive) on the model's own 742-triangle mesh, thresholded at THR mm per
30 d; patch-days above threshold are linked when (a) two mesh triangles share a vertex on the same day or (b) the same
triangle is above threshold on consecutive days. Each connected component = one temporally inter-connected slip patch.
Per component: start/end/duration, peak rate (mm/30 d) + date, footprint (n triangles, area), centroid, along-strike
extent, depth range, migration (fit of daily active-area latitude vs time), cumulative slip on the footprint over the
active span (raw_cum[end] - raw_cum[start-30]) -> potency, M0 (mu 30 GPa), Mw, mean slip.
Bridging (env GAP days, BRIDGE_KM km; suffix _gap{G}_r{R}): components are linked across temporal gaps <= GAP d and spatial
gaps <= BRIDGE_KM km (morphological dilation before linking; footprints unchanged), so one migrating event that briefly dips
below THR is not split.  SWEEP=1 prints a parameter table only.
Outputs (res_catalog_sse/): slip_patches_{THR}mm.csv, slip_patches_{THR}mm.npz (labels T x 742, inc30 mm, mesh),
figures/sse_slip_patches_{THR}mm.png.  Env: THR (mm, default 5); SIGN=-1 -> mirror catalog of back-slip < -THR (suffix _mirror), the noise-floor reference.
"""
import os, sys, numpy as np, pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
REPO = "/home/jovyan/cascadia-sse-strain"; sys.path.insert(0, f"{REPO}/src")
from sse_strain import load_fault_ascii, load_solution
from sse_strain.geodesy import llh2localxy
SD = "fault_tomography/inversion/res_catalog_sse"; THR = float(os.environ.get("THR", "5")); WIN = 30; MU = 30e9
SIGN = float(os.environ.get("SIGN", "1")); SUF = "" if SIGN > 0 else "_mirror"      # SIGN=-1: catalog BACK-slip < -THR (mirror null for the noise floor)
fault = load_fault_ascii(f"{REPO}/data/cascadia_mesh.txt"); sol = load_solution(f"{REPO}/data/2026-08-18", fault); sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); A_raw = np.nan_to_num(zc["A_raw"].astype(float)); T = len(days)
raw_cum = -(Ld @ A_raw.T)                                    # (742, T) m, thrust-positive (sign as sse_trailing_products.py)
inc = np.full_like(raw_cum, np.nan); inc[:, WIN:] = (raw_cum[:, WIN:] - raw_cum[:, :-WIN]) * 1e3      # mm within (t-30, t]
NP = raw_cum.shape[0]; cll = fault.centroid_llh; plon, plat, pdep = cll[:, 0], cll[:, 1], -cll[:, 2]
V = fault.vertices_m()                                        # (NP, 3, 3) m
area = 0.5 * np.linalg.norm(np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0]), axis=1) / 1e6       # km2
# triangle adjacency: shared vertex (rounded to 1 m)
key = {}; adj = [set() for _ in range(NP)]
for p in range(NP):
    for v in np.round(V[p]).astype(np.int64): key.setdefault(tuple(v), []).append(p)
for ps in key.values():
    for a in ps:
        for b in ps:
            if a != b: adj[a].add(b)
B = SIGN * np.nan_to_num(inc) > THR                           # (NP, T); SIGN=-1 -> mirror (back-slip) components
GAP = int(os.environ.get("GAP", "0")); BRIDGE_KM = float(os.environ.get("BRIDGE_KM", "0"))     # bridge temporal gaps <= GAP d and spatial gaps <= BRIDGE_KM km
SUF += (f"_gap{GAP}_r{BRIDGE_KM:g}" if (GAP or BRIDGE_KM) else "")
x_km, y_km = llh2localxy(plat, plon, fault.origin); XY = np.column_stack([x_km, y_km])
EA = np.array(sorted({(a, b) for a in range(NP) for b in adj[a] if a < b}))                     # vertex-sharing triangle pairs
from scipy.ndimage import binary_dilation
from scipy.sparse import csr_matrix
def label_components(B, gap, bridge_km):
    """Connected components of B (NP x T) under same-day vertex adjacency + same-triangle consecutive days, computed on a
    morphologically dilated mask (time half-width ceil(gap/2), space radius bridge_km/2) so that gaps <= gap d / <= bridge_km km
    are bridged; labels are then read back at the ORIGINAL nodes (footprints unchanged)."""
    D = B.copy()
    if bridge_km > 0:
        from scipy.spatial import cKDTree
        pr = cKDTree(XY).query_pairs(bridge_km / 2, output_type="ndarray"); W = csr_matrix((np.ones(len(pr) * 2 + NP), (np.r_[pr[:, 0], pr[:, 1], np.arange(NP)], np.r_[pr[:, 1], pr[:, 0], np.arange(NP)])), shape=(NP, NP))
        D = (W @ D.astype(np.int32)) > 0
    if gap > 0: D = binary_dilation(D, structure=np.ones((1, 2 * int(np.ceil(gap / 2)) + 1), bool))
    nodes = np.flatnonzero(D.ravel()); nid = -np.ones(NP * T, np.int64); nid[nodes] = np.arange(len(nodes))
    pt, tt_ = np.where(D[:, :-1] & D[:, 1:]); r1, c1 = nid[pt * T + tt_], nid[pt * T + tt_ + 1]
    pa, ta = np.where(D[EA[:, 0]] & D[EA[:, 1]]); r2, c2 = nid[EA[pa, 0] * T + ta], nid[EA[pa, 1] * T + ta]
    Gr = coo_matrix((np.ones(len(r1) + len(r2)), (np.r_[r1, r2], np.r_[c1, c2])), shape=(len(nodes), len(nodes)))
    _, labD = connected_components(Gr, directed=False)
    L = np.zeros((NP, T), np.int32); pp, tt = np.where(B); L[pp, tt] = labD[nid[pp * T + tt]] + 1
    u, inv_ = np.unique(L[L > 0], return_inverse=True); L[L > 0] = inv_ + 1                    # renumber 1..ncomp over ORIGINAL nodes
    return L, len(u)
if os.environ.get("SWEEP"):
    gi_ = pd.Series(np.arange(T), index=days); probes = [("2011-06-13", 42.2), ("2011-07-01", 45.3), ("2011-08-15", 48.0)]
    print("gap_d bridge_km | ncomp  non-minor(>=3 tri,>=30 pd)  maxdur_d  max_strike_km | 2011 probes (labels: S-blip, C, N)")
    for gap in (0, 3, 7, 14):
        for br in (0, 50, 100, 200):
            L, nc_ = label_components(B, gap, br); ids, cnt = np.unique(L[L > 0], return_counts=True)
            npatch = np.array([len(np.unique(np.where(L == c)[0])) for c in ids]); nonminor = int(((npatch >= 3) & (cnt >= 30)).sum())
            dur = np.array([np.ptp(np.where(L == c)[1]) + 1 for c in ids]); ext = np.array([np.ptp(plat[np.unique(np.where(L == c)[0])]) * 111 for c in ids])
            lab = []
            for d_, la in probes:
                t = int(gi_[pd.Timestamp(d_)]); act = np.where(B[:, t])[0]; lab.append(int(L[act[np.argmin(np.abs(plat[act] - la))], t]) if len(act) else -1)
            print(f"{gap:5d} {br:9.0f} | {nc_:5d}  {nonminor:4d}                        {dur.max():5d}     {ext.max():6.0f}      | {lab}")
    sys.exit(0)
L, ncomp = label_components(B, GAP, BRIDGE_KM); nodes = np.flatnonzero(B.ravel())
print(f"threshold {THR:g} mm/30 d: {len(nodes)} patch-days above, {ncomp} connected components", flush=True)
rows = []
for c in range(1, ncomp + 1):
    sel = L == c; pc, tc = np.where(sel); ts, te = tc.min(), tc.max(); fp = np.unique(pc)
    k = np.argmax(np.where(sel, SIGN * np.nan_to_num(inc), -np.inf)); pk, tk = divmod(k, T)
    # daily active-area centroid latitude -> migration
    tl = np.arange(ts, te + 1); lat_t = np.array([np.average(plat[sel[:, t]], weights=area[sel[:, t]]) if sel[:, t].any() else np.nan for t in tl])
    ok = np.isfinite(lat_t); mig = np.polyfit(tl[ok] - ts, lat_t[ok], 1)[0] * 111 if ok.sum() > 3 and te > ts else np.nan     # km/day, + = northward
    t0 = max(ts - WIN, 0); dslip = SIGN * (raw_cum[fp, te] - raw_cum[fp, t0]); pot = float((area[fp] * 1e6 * dslip).sum())     # m^3
    M0 = MU * pot; Mw = (2 / 3) * (np.log10(max(M0, 1.0)) - 9.1) if M0 > 0 else np.nan
    rows.append(dict(id=c, start=days[ts].date(), end=days[te].date(), duration_d=int(te - ts + 1), peak_date=days[tk].date(), peak_mm_per30d=float(inc[pk, tk]),
                     n_patches=len(fp), area_km2=float(area[fp].sum()), patch_days=int(sel.sum()), lat_centroid=float(np.average(plat[fp], weights=area[fp])), lon_centroid=float(np.average(plon[fp], weights=area[fp])),
                     lat_min=float(plat[fp].min()), lat_max=float(plat[fp].max()), strike_extent_km=float((plat[fp].max() - plat[fp].min()) * 111 + 28), depth_min_km=float(pdep[fp].min()), depth_max_km=float(pdep[fp].max()),
                     migration_km_per_d=float(mig), mean_slip_mm=float(1e3 * pot / (area[fp].sum() * 1e6)), max_cum_slip_mm=float(1e3 * dslip.max()), potency_km3=pot / 1e9, Mw=float(Mw),
                     minor=bool((len(fp) < 3) or (sel.sum() < 30)),                          # tiny: <3 triangles or <30 patch-days
                     broad_pulse=bool((te - ts + 1 < 14) and ((plat[fp].max() - plat[fp].min()) * 111 > 300))))   # <14 d yet >300 km along strike: margin-wide blip, suspect model/GNSS common mode
cat = pd.DataFrame(rows).sort_values("start"); cat.to_csv(f"{SD}/slip_patches_{THR:g}mm{SUF}.csv", index=False)
np.savez_compressed(f"{SD}/slip_patches_{THR:g}mm{SUF}.npz", label=L.T, inc30_mm=inc.T.astype(np.float32), days=days.astype(str).values, patch_lat=plat, patch_lon=plon, patch_depth_km=pdep, patch_area_km2=area, thr_mm=THR)
maj = cat[~cat.minor]
print(f"components: {len(cat)} total, {len(maj)} non-minor (>=3 triangles & >=30 patch-days, {maj.patch_days.sum()/cat.patch_days.sum()*100:.0f}% of patch-days); broad pulses (<14 d, >300 km): {cat.broad_pulse.sum()}; Mw>=6.0: {(cat.Mw>=6).sum()}, Mw>=6.5: {(cat.Mw>=6.5).sum()}")
big = maj[(maj.Mw >= 6.2) & ~maj.broad_pulse]; print(f"listed below: Mw>=6.2 and not broad pulse ({len(big)})")
print(big.assign(start=maj.start.astype(str), end=maj.end.astype(str)).round(2)[["id","start","end","duration_d","peak_mm_per30d","n_patches","area_km2","lat_min","lat_max","depth_min_km","depth_max_km","migration_km_per_d","mean_slip_mm","Mw"]].to_string(index=False))
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
fig, ax = plt.subplots(2, 1, figsize=(16, 9), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
rng = np.random.RandomState(0); cols = rng.rand(ncomp + 1, 3) * 0.8
pcs, tcs = np.where(L > 0); ax[0].scatter(days[tcs], plat[pcs], c=cols[L[pcs, tcs]], s=4, marker="s", lw=0)
for _, r in big.iterrows(): ax[0].text(pd.Timestamp(r.start), r.lat_max + 0.1, str(r.id), fontsize=6, ha="center")
ax[0].set_ylabel("latitude of mesh triangle"); ax[0].set_title(f"SSE model: {'thrust slip >' if SIGN > 0 else 'MIRROR: back-slip < -'}{THR:g} mm within the trailing 30 d, space-time connected components (colour = component; labels = Mw>=6.2, not broad pulses)", fontsize=10); ax[0].grid(alpha=.3)
ax[1].bar(pd.to_datetime(cat.start), cat.Mw, width=cat.duration_d, color=cols[cat.id], align="edge"); ax[1].set_ylabel("Mw of component"); ax[1].set_ylim(5, 7.2); ax[1].grid(alpha=.3)
ax[1].set_xlim(pd.Timestamp("2007-01-01"), pd.Timestamp("2026-09-01")); fig.savefig(f"figures/sse_slip_patches_{THR:g}mm{SUF}.png", dpi=120, bbox_inches="tight")
print(f"wrote {SD}/slip_patches_{THR:g}mm{SUF}.csv/.npz, figures/sse_slip_patches_{THR:g}mm{SUF}.png")
