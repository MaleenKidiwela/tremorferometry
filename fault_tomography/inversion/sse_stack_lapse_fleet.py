#!/usr/bin/env python
"""Fleet-wide DIRECT stack test of the ETS-locked change (2026-09-04): for every station, per certified family, the mean
daily stack during ON windows [+15,+60] of admitted slip episodes (footprint within 150 km) vs the mean stack during
PRE windows [-90,-30]; stretch ON vs PRE about the S (t_s=+1.0 s) in the 2-4 s and 3-5 s stack windows (S-pinned:
t_min/t_max = window - 1.0). Lapse consistency per station: median ratio e35/e24 and sign agreement; a near-station
velocity change gives same sign, ratio ~0.6-1.2; a waveform/noise-composition change gives sign flips or ratio ~0.
Bootstrap (resample ON/PRE days, 40 draws) -> per-station SE; time-half replication (episodes <=2018 vs >=2019).
Also per family: coda/direct energy change, pre-arrival RMS ratio ON/PRE (noise-field proxy).
Output: res_catalog_sse/stack_lapse_fleet.csv (station rows) + stack_lapse_fleet_families.csv; figure later.
"""
import numpy as np, pandas as pd, os, sys, glob, time
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"
sys.path.insert(0, "src"); from tremorferometry.dvv import stretch_dvv
rng = np.random.RandomState(5)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); gi = pd.Series(np.arange(G), index=days)
clat, clon = zc["cell_lat"], zc["cell_lon"]; km = lambda la, lo: np.column_stack([(lo + 124) * 111 * np.cos(np.radians(47)), (la - 47) * 111]); XY = km(clat, clon)
ep = pd.read_csv(f"{SD}/episodes.csv"); ep["key"] = ep.seg + "_" + ep.peak; okep = ep[ep.ok & ep.seg.isin(["N", "C", "S"])]; zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True)
sll = pd.read_csv(f"{RD}/pairs.csv").groupby("tag")[["sta_lat", "sta_lon"]].first(); fams = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "fam"]).drop_duplicates()
def st(wb, wa, fs, w, emax, n):
    i0 = int(np.searchsorted(T_AX, 1.0)); return -stretch_dvv(wb[i0:], wa[i0:], fs=fs, t_min=w[0] - 1.0, t_max=w[1] - 1.0, eps_max=emax, n_eps=n).dvv * 100
SHARD, NSHARD = int(os.environ.get("SHARD", "0")), int(os.environ.get("NSHARD", "1"))
OUT_S = f"{SD}/stack_lapse_fleet_s{SHARD}.csv" if NSHARD > 1 else f"{SD}/stack_lapse_fleet.csv"
rows, frows = [], []
for k_, tag in enumerate(sorted(fams.tag.unique())):
    if k_ % NSHARD != SHARD: continue
    fn = f"data/long_window_daily_{tag}_Z.npz"
    if not os.path.exists(fn) or tag not in sll.index: continue
    t0 = time.time()
    try:
        z = np.load(fn, allow_pickle=True); T_AX = np.asarray(z["t"], float); fs = float(z["fs"]); P = np.asarray(z["patches"]).astype(str); D = pd.to_datetime(np.asarray(z["dates"]).astype(str))
        sel = np.isin(P, list(set(fams[fams.tag == tag].fam))); idx = np.where(sel)[0]
        if len(idx) < 100: continue
        S = np.asarray(z["stacks"][idx], np.float32); P = P[idx]; D = D[idx]; del z
    except Exception as e:
        print(f"{tag}: load failed {e}", flush=True); continue
    sxy = km(np.array([sll.sta_lat[tag]]), np.array([sll.sta_lon[tag]]))[0]
    eps = [(int(gi[pd.Timestamp(r.fp_onset)]), int(r.peak[:4])) for _, r in okep.iterrows() if np.sqrt(((XY[zf[f"fp_{r.key}"]] - sxy) ** 2).sum(1)).min() <= 150]
    if len(eps) < 4: continue
    def masks(eplist):
        on = np.zeros(len(D), bool); pre = np.zeros(len(D), bool)
        for i, _ in eplist: on |= (D >= days[i] + pd.Timedelta(days=15)) & (D <= days[min(i + 60, G - 1)]); pre |= (D >= days[max(i - 90, 0)]) & (D <= days[max(i - 30, 0)])
        return on, pre
    on, pre = masks(eps); on1, pre1 = masks([e for e in eps if e[1] <= 2018]); on2, pre2 = masks([e for e in eps if e[1] >= 2019])
    w24, w35, dr, pn = (T_AX >= 2) & (T_AX <= 4), (T_AX >= 3) & (T_AX <= 5), (T_AX >= 0.5) & (T_AX <= 1.5), (T_AX >= -2) & (T_AX < 0)
    e24, e35, e24_1, e24_2, e35_1, e35_2, cdr, nz, boots = [], [], [], [], [], [], [], [], []
    for fam in np.unique(P):
        m = P == fam; a = S[m & on]; b = S[m & pre]
        if len(a) < 15 or len(b) < 15: continue
        wa = a.mean(0).astype(float); wb = b.mean(0).astype(float); wa /= np.abs(wa).max(); wb /= np.abs(wb).max()
        v24 = st(wb, wa, fs, (2, 4), 0.02, 201); v35 = st(wb, wa, fs, (3, 5), 0.03, 301); e24.append(v24); e35.append(v35)
        cdr.append(np.log((wa[w24] ** 2).mean() / (wa[dr] ** 2).mean()) - np.log((wb[w24] ** 2).mean() / (wb[dr] ** 2).mean())); nz.append(np.sqrt((wa[pn] ** 2).mean() / (wb[pn] ** 2).mean()))
        for k in range(4):                                  # bootstrap over days (pooled later)
            aa = a[rng.randint(0, len(a), len(a))].mean(0).astype(float); bb = b[rng.randint(0, len(b), len(b))].mean(0).astype(float); aa /= np.abs(aa).max(); bb /= np.abs(bb).max()
            boots.append(st(bb, aa, fs, (2, 4), 0.02, 201))
        a1, b1, a2, b2 = S[m & on1], S[m & pre1], S[m & on2], S[m & pre2]
        if len(a1) >= 15 and len(b1) >= 15:
            x = a1.mean(0).astype(float) / np.abs(a1.mean(0)).max(); y = b1.mean(0).astype(float) / np.abs(b1.mean(0)).max(); e24_1.append(st(y, x, fs, (2, 4), 0.02, 201)); e35_1.append(st(y, x, fs, (3, 5), 0.03, 301))
        if len(a2) >= 15 and len(b2) >= 15:
            x = a2.mean(0).astype(float) / np.abs(a2.mean(0)).max(); y = b2.mean(0).astype(float) / np.abs(b2.mean(0)).max(); e24_2.append(st(y, x, fs, (2, 4), 0.02, 201)); e35_2.append(st(y, x, fs, (3, 5), 0.03, 301))
        frows.append(dict(tag=tag, fam=fam, e24=v24, e35=v35, cdr=cdr[-1], noise_ratio=nz[-1], n_on=len(a), n_pre=len(b)))
    if len(e24) < 5: continue
    e24, e35 = np.array(e24), np.array(e35); ratio = np.median(e35 / e24) if np.median(np.abs(e24)) > 0.02 else np.nan
    rows.append(dict(tag=tag, lat=float(sll.sta_lat[tag]), lon=float(sll.sta_lon[tag]), borehole="p90f40" in tag, n_fam=len(e24), n_ep=len(eps),
                     e24=float(np.median(e24)), e24_fpos=float(np.mean(e24 > 0)), e24_se=float(np.std(boots) / np.sqrt(max(len(boots) / 4, 1))), e35=float(np.median(e35)), e35_fpos=float(np.mean(e35 > 0)),
                     ratio35_24=float(ratio), sign_agree=float(np.mean(np.sign(e24) == np.sign(e35))), cdr=float(np.median(cdr)), noise_ratio=float(np.median(nz)),
                     e24_early=float(np.median(e24_1)) if e24_1 else np.nan, e24_late=float(np.median(e24_2)) if e24_2 else np.nan,
                     e35_early=float(np.median(e35_1)) if e35_1 else np.nan, e35_late=float(np.median(e35_2)) if e35_2 else np.nan))
    r = rows[-1]; print(f"{tag:11s} fams {r['n_fam']:3d} ep {r['n_ep']:2d} | e24 {r['e24']:+.3f} (fpos {r['e24_fpos']:.2f}, se {r['e24_se']:.3f}) e35 {r['e35']:+.3f} ratio {r['ratio35_24']:+.2f} agree {r['sign_agree']:.2f} | cdr {r['cdr']:+.3f} noise {r['noise_ratio']:.3f} | early/late e24 {r['e24_early']:+.3f}/{r['e24_late']:+.3f} | {time.time() - t0:.0f}s", flush=True)
    pd.DataFrame(rows).to_csv(OUT_S, index=False); pd.DataFrame(frows).to_csv(OUT_S.replace(".csv", "_families.csv"), index=False)
print("FLEET_STACK_DONE")
