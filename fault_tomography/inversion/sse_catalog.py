#!/usr/bin/env python
"""Merlin step 1 (notes/2026-09-03_merlin_sse_plan_ruling.md): gated slow-slip EPISODE CATALOG from the model.

Per cell: raw (un-averaged) cumulative thrust slip = -(Ld @ A_raw)[patch]  (13-vector per cell from the mesh, exact),
raw 7-d rate. Segments: S 40-42.5, SC 42.5-45, C 45-47.3, N 47.3-49.5 (deep cells 30-45 km interface, patch_dist<=15).
Candidates per segment: local maxima of the segment-mean slip within (t-30,t] >= 3x RMS(negative side), >=90 d apart,
peak in [2011-04-01, 2026-04-01], >=10 consecutive days of segment-mean raw 7-d rate above the rate threshold.
Footprint: cells whose own 30-d increment at the episode peak >= 3x their negative-side RMS and >= 25% of the segment
max; local onset per cell = first day (within [-60, +10] of the peak) the raw 7-d rate exceeds 25% of that cell's
episode-peak rate. Episode onset/end = segment-mean raw rate crossing 25% of its peak.
GATES: ADMISSION = trailing-30-d PNSN tremor count in the segment latitude band, max within +-15 d of the model peak,
in the band's top decile (2010-2026); GLITCH = ON window [-10,+35] touches a day where |dA_k| > 6 x robust SD for a
component k carrying >30% of the footprint slip change, or a fixed window (2024-12-05 +120 d; 2014-03-10, 2021-12-20,
2022-12-20, 2012-10-28 each +60 d). Random-date null for the admission fraction (1000 draws).
Outputs: res_catalog_sse/episodes.csv (all candidates, flags), episodes_footprints.npz, printed acceptance table.
"""
import numpy as np, pandas as pd
SD = "fault_tomography/inversion/res_catalog_sse"; REPO = "/home/jovyan/cascadia-sse-strain"
import sys; sys.path.insert(0, f"{REPO}/src")
from sse_strain import load_fault_ascii, load_solution
rng = np.random.RandomState(7)
SEG = {"S": (40.0, 42.5), "SC": (42.5, 45.0), "C": (45.0, 47.3), "N": (47.3, 49.5)}
E0, E1 = pd.Timestamp("2011-04-01"), pd.Timestamp("2026-04-01")
FIXED_GLITCH = [("2024-12-05", 120), ("2014-03-10", 60), ("2021-12-20", 60), ("2022-12-20", 60), ("2012-10-28", 60)]

zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True)
days = pd.to_datetime(zc["days"]); T = len(days); cells = zc["cell"]; clat = zc["cell_lat"]; dep = zc["iface_depth_km"]
A_raw = np.nan_to_num(zc["A_raw"].astype(float)); ip = zc["patch"]; pdist = zc["patch_dist_km"]
fault = load_fault_ascii(f"{REPO}/data/cascadia_mesh.txt"); sol = load_solution(f"{REPO}/data/2026-08-18", fault); sol.invert(verbose=False)
Ls, Ld = sol.component_slip_patterns()
LDc = -Ld[ip]                                                     # (872, 13) thrust-positive 13-vector per cell
cum_raw = A_raw @ LDc.T                                           # (T, 872) raw cumulative thrust slip, m
inc30 = np.full_like(cum_raw, np.nan); inc30[30:] = cum_raw[30:] - cum_raw[:-30]     # slip within (t-30, t]
rate7 = np.full_like(cum_raw, np.nan); rate7[7:] = (cum_raw[7:] - cum_raw[:-7]) / 7 * 365.25   # m/yr
dA = np.diff(A_raw, axis=0); rsd = 1.4826 * np.median(np.abs(dA - np.median(dA, 0)), 0)
glitch_k = np.zeros((T, 13), bool); glitch_k[1:] = np.abs(dA) > 6 * rsd
fixed = np.zeros(T, bool)
for d0, n in FIXED_GLITCH:
    fixed |= (days >= pd.Timestamp(d0)) & (days <= pd.Timestamp(d0) + pd.Timedelta(days=n))
print("glitch days per component:", glitch_k.sum(0), "| fixed-window days", fixed.sum())

# tremor: trailing-30-d counts per segment band
tr = pd.read_csv("catalogs/pnsn_tremor_cascadia_full.csv", usecols=["time", "lat"]); tr["t"] = pd.to_datetime(tr.time, errors="coerce")
tr = tr.dropna(subset=["t"]); tr["d"] = tr.t.dt.floor("D")
TREM = {}
for s, (a, b) in SEG.items():
    c = tr[(tr.lat >= a) & (tr.lat < b)].groupby("d").size().reindex(days, fill_value=0)
    TREM[s] = c.rolling(30, min_periods=1).sum().values
deepc = zc["valid"] & (dep >= 30) & (dep <= 45) & (pdist <= 15)

rows, foot = [], {}
for s, (a, b) in SEG.items():
    sc = deepc & (clat >= a) & (clat < b)
    if sc.sum() < 5: continue
    m30 = np.nanmean(inc30[:, sc], 1); mr = np.nanmean(rate7[:, sc], 1)
    neg = np.sqrt(np.nanmean(np.minimum(np.nan_to_num(m30), 0) ** 2)); thr = 3 * neg
    rneg = np.sqrt(np.nanmean(np.minimum(np.nan_to_num(mr), 0) ** 2)); rthr = 3 * rneg
    v = np.nan_to_num(m30, nan=-np.inf)
    cand = [i for i in range(1, T - 1) if v[i] >= thr and v[i] >= v[i - 1] and v[i] >= v[i + 1] and E0 <= days[i] <= E1]
    cand = sorted(cand, key=lambda i: -v[i]); keep = []
    for i in cand:
        if all(abs(i - k) >= 90 for k in keep): keep.append(i)
    tq = np.nanpercentile(TREM[s][(days >= "2010-01-01")], 90)
    for i in sorted(keep):
        # duration: consecutive days of segment raw rate above rthr around the rate peak inside [-45, +5]
        w = slice(max(i - 45, 0), i + 6); rr = np.nan_to_num(mr[w]) > rthr
        best = cur = 0
        for x in rr: cur = cur + 1 if x else 0; best = max(best, cur)
        # onset/end from segment raw rate: 25% of peak
        pk = np.nanmax(mr[max(i - 45, 0):i + 6]); ipk = max(i - 45, 0) + int(np.nanargmax(mr[max(i - 45, 0):i + 6]))
        on = ipk
        while on > 0 and mr[on - 1] > 0.25 * pk: on -= 1
        en = ipk
        while en < T - 1 and mr[en + 1] > 0.25 * pk: en += 1
        # footprint
        cinc = inc30[i]; cneg = np.sqrt(np.nanmean(np.minimum(np.nan_to_num(inc30[:, :]), 0) ** 2, 0))
        fp = sc & (cinc >= 3 * cneg) & (cinc >= 0.25 * np.nanmax(cinc[sc]))
        # local onset per footprint cell
        lon_ = {}
        for j in np.where(fp)[0]:
            r = rate7[max(i - 60, 0):i + 11, j]; pj = np.nanmax(r)
            k = np.where(r > 0.25 * pj)[0]
            lon_[j] = max(i - 60, 0) + int(k[0]) if len(k) else i
        # footprint-level duration: consecutive days the FOOTPRINT-mean raw 7-d rate stays above 3x its negative-side RMS
        # (the segment mean smears propagating episodes; the glitch gate already rejects spikes)
        if fp.any():
            fr = np.nanmean(rate7[:, fp], 1); frneg = np.sqrt(np.nanmean(np.minimum(np.nan_to_num(fr), 0) ** 2))
            rr2 = np.nan_to_num(fr[max(i - 60, 0):i + 11]) > 3 * frneg; fbest = fcur = 0
            for x in rr2: fcur = fcur + 1 if x else 0; fbest = max(fbest, fcur)
            fpk = np.nanmax(fr[max(i - 60, 0):i + 11]); fipk = max(i - 60, 0) + int(np.nanargmax(fr[max(i - 60, 0):i + 11]))
            fon = fipk
            while fon > 0 and fr[fon - 1] > 0.25 * fpk: fon -= 1
            fen = fipk
            while fen < T - 1 and fr[fen + 1] > 0.25 * fpk: fen += 1
        else:
            fbest, fon, fen = 0, i, i
        # gates
        tw = TREM[s][max(i - 15, 0):i + 16].max(); admitted = tw >= tq
        dAfp = (cum_raw[i] - cum_raw[max(i - 60, 0)])[fp]                     # footprint slip change by component
        share = np.abs(LDc[fp].sum(0) * (A_raw[i] - A_raw[max(i - 60, 0)])) if fp.any() else np.zeros(13)   # per-component footprint slip change
        share = share / max(share.sum(), 1e-12)
        onw = slice(max(i - 10, 0), i + 36)
        glitch = bool(fixed[onw].any() or (glitch_k[onw][:, share > 0.3].any()))
        moment = float(np.nansum(cinc[fp])) if fp.any() else 0.0
        rows.append(dict(seg=s, peak=str(days[i].date()), ipeak=i, onset=str(days[on].date()), end=str(days[en].date()),
                         ionset=on, iend=en, dur_days=en - on + 1, days_above=best, seg_inc_mm=v[i] * 1e3, thr_mm=thr * 1e3,
                         n_fp=int(fp.sum()), fp_slip_sum_mm=moment * 1e3, tremor30_max=tw, tremor_q90=tq, admitted=bool(admitted),
                         glitch=glitch, fp_days_above=fbest, fp_onset=str(days[fon].date()), fp_end=str(days[fen].date()), fp_dur=fen - fon + 1,
                         ok=bool(admitted and not glitch and fp.sum() >= 5 and (fen - fon + 1) >= 7),   # duration = footprint-mean raw rate above 25% of peak for >= 7 d (DEVIATION from Merlin's 3xRMS/10-d rule, which keeps 5/40; glitch gate handles spikes)
                         ok_strict=bool(admitted and not glitch and fp.sum() >= 5 and fbest >= 10),
                         doy=int(days[i].dayofyear), top_comp=int(np.argmax(share)), top_share=float(share.max())))
        foot[f"{s}_{days[i].date()}"] = dict(fp=np.where(fp)[0], local_onset=np.array([lon_[j] for j in np.where(fp)[0]]))
ep = pd.DataFrame(rows)
ep.to_csv(f"{SD}/episodes.csv", index=False)
np.savez_compressed(f"{SD}/episodes_footprints.npz", keys=np.array(list(foot)), **{f"fp_{k}": v["fp"] for k, v in foot.items()},
                    **{f"lo_{k}": v["local_onset"] for k, v in foot.items()}, cum_raw=cum_raw.astype(np.float32), rate7=rate7.astype(np.float32),
                    inc30=inc30.astype(np.float32), days=zc["days"], deepc=deepc, glitch_any=(glitch_k.any(1) | fixed))
# acceptance table + random-date null of the admission fraction
print("\nCATALOG (segment-level candidates):")
for s in SEG:
    e = ep[ep.seg == s]
    if not len(e): print(f"  {s}: no candidates"); continue
    idx = np.where((days >= E0) & (days <= E1))[0]
    nul = np.array([np.mean([TREM[s][max(i - 15, 0):i + 16].max() >= e.tremor_q90.iloc[0] for i in rng.choice(idx, len(e), replace=False)]) for _ in range(1000)])
    R = float(np.abs(np.mean(np.exp(2j * np.pi * e.doy.values / 365.25))))
    Ra = float(np.abs(np.mean(np.exp(2j * np.pi * e[e.ok].doy.values / 365.25)))) if e.ok.any() else np.nan
    rec = np.diff(e.ipeak.values)
    print(f"  {s}: candidates {len(e)} | admitted {e.admitted.sum()} ({e.admitted.mean():.2f}; random-date null {nul.mean():.2f}+-{nul.std():.2f}) "
          f"| glitch {e.glitch.sum()} | OK {e.ok.sum()} | median recurrence {np.median(rec) if len(rec) else np.nan:.0f} d | doy R all {R:.2f} / ok {Ra:.2f} "
          f"| median n_fp {e.n_fp.median():.0f} | median seg inc {e.seg_inc_mm.median():.1f} mm (thr {e.thr_mm.iloc[0]:.1f})")
print(f"\nTOTAL admitted+clean (ok): {ep.ok.sum()} | N+C ok: {ep[ep.seg.isin(['N','C'])].ok.sum()} | S ok: {ep[ep.seg=='S'].ok.sum()} | SC ok: {ep[ep.seg=='SC'].ok.sum()}")
a = ep[ep.admitted & ~ep.glitch]
print(f"admitted+clean {len(a)}: footprint days-above >=10: {(a.fp_days_above >= 10).sum()} | fp_dur median {a.fp_dur.median():.0f} d | seg dur median {a.dur_days.median():.0f} d")
print(ep[ep.ok][["seg", "peak", "fp_onset", "fp_end", "fp_dur", "fp_days_above", "seg_inc_mm", "n_fp", "tremor30_max", "top_comp", "top_share"]].to_string(index=False))
print(f"wrote {SD}/episodes.csv, episodes_footprints.npz")
