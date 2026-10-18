#!/usr/bin/env python
"""Null A for the block-gated 3-5 s footprint mean (lean): gated per-(cell,day) series built once (rail rows removed,
flagged station-episode blocks removed by masking those stations' rows for that episode's cells is not possible per
episode in a single grid, so the gate here removes rail rows fleet-wide and drops the flagged station-episode
contributions at evaluation time via a per-episode station exclusion list). 300 rigid shifts per segment."""
import numpy as np, pandas as pd, glob, os
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; rng = np.random.RandomState(29)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; M = len(cells); gi = pd.Series(np.arange(G), index=days); cpos = {c: i for i, c in enumerate(cells)}
ep = pd.read_csv(f"{SD}/episodes.csv"); ep["key"] = ep.seg + "_" + ep.peak; okep = ep[ep.ok & ep.seg.isin(["N", "C"])]; zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True)
EP = [dict(key=r.key, seg=r.seg, ion=int(gi[pd.Timestamp(r.fp_onset)]), fp=zf[f"fp_{r.key}"]) for _, r in okep.iterrows()]
E0, E1 = int(gi[pd.Timestamp("2011-04-01")]), int(gi[pd.Timestamp("2026-04-01")]); span = E1 - E0; wrap = lambda i: E0 + (i - E0) % span
B = pd.read_csv(f"{SD}/lapse_gated_blocks.csv"); flag = {k: set(B[(B.key == k) & B.flagged].tag) for k in B.key.unique()}
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "fam", "cell"]).drop_duplicates(); fd = fd[fd.cell.isin(cpos)]; famset = {t: set(g.fam) for t, g in fd.groupby("tag")}
# per-station gated cell-day means (rail rows removed), stored per station so flagged blocks can be excluded per episode
month = days.month.values; ST = {}
for f in glob.glob("data/daily_dvv_*_Z_3to5.csv"):
    tag = os.path.basename(f).split("daily_dvv_")[1].split("_Z_3to5")[0]
    if tag not in famset: continue
    x = pd.read_csv(f, usecols=["patch", "date", "dvv", "cc_max"]); x = x[(x.cc_max > 0.7) & (x.dvv.abs() < 0.0297) & x.patch.isin(famset[tag])]
    x = x.merge(fd[fd.tag == tag][["fam", "cell"]].rename(columns={"fam": "patch"}), on="patch"); x["di"] = gi.reindex(pd.to_datetime(x.date)).values; x = x.dropna(subset=["di"]); x["v"] = -x.dvv * 100
    g = x.groupby([x.cell.map(cpos), x.di.astype(int)]).v.mean(); ST[tag] = g
print(f"stations {len(ST)}", flush=True)
# assemble a per-station sparse (cell, day) table -> for speed build one dense grid per station only for footprint cells
allfp = sorted(set(np.concatenate([e["fp"] for e in EP]))); fpi = {c: k for k, c in enumerate(allfp)}
grids = {}
for tag, g in ST.items():
    A = np.full((G, len(allfp)), np.nan, np.float32); idx = g.index.get_level_values(0).values; di = g.index.get_level_values(1).values; keep = np.isin(idx, allfp)
    A[di[keep], [fpi[c] for c in idx[keep]]] = g.values[keep]
    # monthly deseason (episode-masked)
    mask = np.zeros(G, bool)
    for e in EP: mask[max(e["ion"] - 30, 0):e["ion"] + 121] = True
    for m in range(1, 13): sel = (month == m) & ~mask; A[month == m] -= np.nanmean(A[sel], 0)[None]
    grids[tag] = A
def gated_mean(shift):
    out = []
    for e in EP:
        i = wrap(e["ion"] + shift.get(e["seg"], 0)) if shift else e["ion"]; on = slice(max(i + 15, 0), i + 61); pre = slice(max(i - 90, 0), i - 29); cols = [fpi[c] for c in e["fp"]]
        num = np.zeros(len(cols)); cnt = np.zeros(len(cols))
        for tag, A in grids.items():
            if tag in flag.get(e["key"], set()): continue
            a = np.nanmean(A[on][:, cols], 0); b = np.nanmean(A[pre][:, cols], 0); va = np.isfinite(A[on][:, cols]).mean(0) > .6; vb = np.isfinite(A[pre][:, cols]).mean(0) > .6; ok = va & vb & np.isfinite(a) & np.isfinite(b)
            num[ok] += (a - b)[ok]; cnt[ok] += 1
        v = num[cnt > 0] / cnt[cnt > 0]; out.append(np.mean(v) if len(v) >= 3 else np.nan)
    out = np.array(out); return np.nanmean(out), int(np.nansum(out < 0)), int(np.isfinite(out).sum())
m0, K0, N0 = gated_mean({}); print(f"gated pooled 3-5 s footprint mean = {m0:+.3f}% (K={K0}/{N0})", flush=True)
nul = np.array([gated_mean({g: int(rng.randint(90, span - 90)) for g in ("N", "C")})[0] for _ in range(300)])
print(f"Null A (300 rigid shifts): p = {np.mean(np.abs(nul) >= abs(m0)):.3f} (sd {np.nanstd(nul):.3f})"); print("GATED_NULL_DONE")
