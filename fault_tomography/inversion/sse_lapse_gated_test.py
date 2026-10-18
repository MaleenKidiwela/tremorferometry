#!/usr/bin/env python
"""Merlin's decisive test for the 3-5 s footprint decrease (audit 2026-09-04, prescription step 2), plus step 1 (July-2018
per-station paired table) and step 3 (per-episode 3-5 s table saved under its own name).
Family-level paired lapse test on FOOTPRINT rows of the 21 N+C tremor-admitted episodes: per family-episode, the
ON[+15,+60]-PRE[-90,-30] contrast in 3-5 s vs 2-4 s from the same stacks (cc>0.7 each window; monthly deseason).
Variant B (gated): exclude rail rows (|eps| >= 0.99 eps_max of the window) and any station-episode block whose 3-5 s ON
rail fraction > 3x PRE or whose 3-5 s cc median drops > 0.1 ON vs PRE. Statistics: Theil-Sen slope c35 on c24 with a
station-resampling bootstrap, sign agreement, pooled gated 3-5 s mean with Null A (rigid shift).
RULE: slope in [0.6, 1.2] (CI excl. 0) AND sign agreement >= 70% AND gated mean <= -0.04% at p < 0.01 -> lapse-consistent
candidate; else waveform-change class -> BOUND, close.  Outputs: res_catalog_sse/lapse_gated_{blocks,pairs,episodes35}.csv
"""
import numpy as np, pandas as pd, glob, os
from scipy.stats import theilslopes
SD = "fault_tomography/inversion/res_catalog_sse"; RD = "fault_tomography/inversion/res_catalog"; rng = np.random.RandomState(29)
zc = np.load(f"{SD}/components_at_cells_trailing30.npz", allow_pickle=True); days = pd.to_datetime(zc["days"]); G = len(days); cells = zc["cell"]; gi = pd.Series(np.arange(G), index=days)
ep = pd.read_csv(f"{SD}/episodes.csv"); ep["key"] = ep.seg + "_" + ep.peak; okep = ep[ep.ok & ep.seg.isin(["N", "C"])]; zf = np.load(f"{SD}/episodes_footprints.npz", allow_pickle=True)
EP = [dict(key=r.key, seg=r.seg, ion=int(gi[pd.Timestamp(r.fp_onset)]), fp=set(cells[zf[f"fp_{r.key}"]])) for _, r in okep.iterrows()]
E0, E1 = int(gi[pd.Timestamp("2011-04-01")]), int(gi[pd.Timestamp("2026-04-01")]); span = E1 - E0; wrap = lambda i: E0 + (i - E0) % span
fd = pd.read_parquet(f"{RD}/family_days.parquet", columns=["tag", "fam", "cell"]).drop_duplicates(); fam2cell = fd.set_index(["tag", "fam"]).cell
allfp = set().union(*[e["fp"] for e in EP]); fd = fd[fd.cell.isin(allfp)]; famset = {t: set(g.fam) for t, g in fd.groupby("tag")}
month = days.month.values
def load(tag, w):
    fn = f"data/daily_dvv_{tag}_Z_{w}.csv"
    if not os.path.exists(fn): return None
    x = pd.read_csv(fn, usecols=["patch", "date", "dvv", "cc_max"]); x = x[x.patch.isin(famset[tag])].rename(columns={"patch": "fam"}); x["date"] = pd.to_datetime(x.date)
    emax = 0.02 if w == "2to4" else 0.03; x["rail"] = x.dvv.abs() >= 0.99 * emax; x["v"] = -x.dvv * 100; return x
pairs, blocks = [], []
for tag in sorted(famset):
    a24, a35 = load(tag, "2to4"), load(tag, "3to5")
    if a24 is None or a35 is None: continue
    for w, x in (("24", a24), ("35", a35)):
        x["di"] = gi.reindex(x.date).values
    for x in (a24, a35):
        x.dropna(subset=["di"], inplace=True); x["di"] = x.di.astype(int)
    def piv(x, gate):
        y = x[(x.cc_max > 0.7) & (~x.rail if gate else True)]; p = y.pivot_table(index="di", columns="fam", values="v").reindex(range(G)); p = p - p.groupby(month).transform("mean"); return p
    P24, P35, P35g = piv(a24, False), piv(a35, False), piv(a35, True)
    for e in EP:
        fams_e = [f for f in P24.columns if fam2cell.get((tag, f), None) in e["fp"]]
        fams_e = [f for f in fams_e if f in P35.columns]
        if not fams_e: continue
        i = e["ion"]; on = slice(max(i + 15, 0), i + 61); pre = slice(max(i - 90, 0), i - 29)
        b35 = a35[a35.fam.isin(fams_e)]; bon = b35[(b35.di >= i + 15) & (b35.di <= i + 60)]; bpre = b35[(b35.di >= i - 90) & (b35.di <= i - 30)]
        if len(bon) < 30 or len(bpre) < 30: continue
        rail_on, rail_pre = bon.rail.mean(), bpre.rail.mean(); cc_on, cc_pre = bon.cc_max.median(), bpre.cc_max.median()
        b24 = a24[a24.fam.isin(fams_e)]; r24_on = b24[(b24.di >= i + 15) & (b24.di <= i + 60)].rail.mean()
        flagged = bool(rail_on > 3 * max(rail_pre, 0.005) or (cc_pre - cc_on) > 0.1)
        blocks.append(dict(tag=tag, key=e["key"], n_fam=len(fams_e), rail35_on=rail_on, rail35_pre=rail_pre, rail24_on=r24_on, cc35_on=cc_on, cc35_pre=cc_pre, flagged=flagged))
        for f in fams_e:
            c24 = np.nanmean(P24[f].values[on]) - np.nanmean(P24[f].values[pre]); c35 = np.nanmean(P35[f].values[on]) - np.nanmean(P35[f].values[pre]); c35g = np.nanmean(P35g[f].values[on]) - np.nanmean(P35g[f].values[pre]) if f in P35g.columns else np.nan
            v24 = np.isfinite(P24[f].values[on]).mean() > .6 and np.isfinite(P24[f].values[pre]).mean() > .6; v35 = np.isfinite(P35[f].values[on]).mean() > .6 and np.isfinite(P35[f].values[pre]).mean() > .6
            if v24 and v35: pairs.append(dict(tag=tag, key=e["key"], fam=f, c24=c24, c35=c35, c35_gated=c35g, flagged=flagged))
B = pd.DataFrame(blocks); Pp = pd.DataFrame(pairs); B.to_csv(f"{SD}/lapse_gated_blocks.csv", index=False); Pp.to_csv(f"{SD}/lapse_gated_pairs.csv", index=False)
print(f"station-episode blocks {len(B)} | flagged (rail ON > 3x PRE or cc drop > 0.1): {B.flagged.sum()} | family-episode pairs {len(Pp)}")
# step 1: July-2018 per-station table
j = B[B.key == "N_2018-07-15"].merge(Pp[Pp.key == "N_2018-07-15"].groupby("tag").agg(c24=("c24", "mean"), c35=("c35", "mean"), fneg35=("c35", lambda v: (v < 0).mean())).reset_index(), on="tag")
print("\nJULY-2018 per-station paired table (footprint families):"); print(j.sort_values("c35")[["tag", "n_fam", "c24", "c35", "fneg35", "rail35_on", "rail35_pre", "rail24_on", "cc35_on", "cc35_pre", "flagged"]].round(3).to_string(index=False))
# step 2: paired lapse test, all rows and gated
def summarize(name, d, col):
    d = d.dropna(subset=["c24", col]); sl = theilslopes(d[col].values, d.c24.values)[0]
    tags = d.tag.unique(); bs = []
    for _ in range(300):
        pick = rng.choice(tags, len(tags), replace=True); dd = pd.concat([d[d.tag == t] for t in pick]); bs.append(theilslopes(dd[col].values, dd.c24.values)[0])
    agree = np.mean(np.sign(d.c24) == np.sign(d[col]))
    print(f"{name}: n={len(d):,} family-episodes, {len(tags)} stations | Theil-Sen slope c35 on c24 = {sl:.2f} (station bootstrap 95% [{np.percentile(bs, 2.5):.2f}, {np.percentile(bs, 97.5):.2f}]) | sign agreement {100*agree:.0f}% | mean c24 {d.c24.mean():+.3f} c35 {d[col].mean():+.3f}")
    return sl, agree
summarize("ALL rows", Pp, "c35"); sl, agree = summarize("GATED (rail rows + flagged blocks removed)", Pp[~Pp.flagged], "c35_gated")
# gated pooled footprint mean with Null A: rebuild per-episode footprint means from the gated family contrasts (cell-mean then footprint-mean)
Pg = Pp[~Pp.flagged].dropna(subset=["c35_gated"]).copy(); Pg["cell"] = [fam2cell.get((t, f)) for t, f in zip(Pg.tag, Pg.fam)]
s = Pg.groupby(["key", "cell"]).c35_gated.mean().groupby("key").mean().reindex([e["key"] for e in EP])
print(f"\ngated pooled 3-5 s footprint mean (from the paired table): {s.mean():+.3f}% (K={int((s<0).sum())}/{int(s.notna().sum())} decrease)")
# Null A for the gated statistic: recompute the per-episode footprint contrasts from gated pivots at shifted times (station-level cache)
print("Null A for the gated mean: computing from gated pivots ...", flush=True)
cache = {}
for tag in sorted(famset):
    a35 = load(tag, "3to5")
    if a35 is None: continue
    a35["di"] = gi.reindex(a35.date).values; a35 = a35.dropna(subset=["di"]); a35["di"] = a35.di.astype(int)
    y = a35[(a35.cc_max > 0.7) & ~a35.rail]; p = y.pivot_table(index="di", columns="fam", values="v").reindex(range(G)); p = p - p.groupby(month).transform("mean")
    cache[tag] = (p.values.astype(np.float32), list(p.columns))
flagset = set(zip(B[B.flagged].tag, B[B.flagged].key))
def gated_mean(shift):
    out = []
    for e in EP:
        i = wrap(e["ion"] + shift.get(e["seg"], 0)) if shift else e["ion"]; on = slice(max(i + 15, 0), i + 61); pre = slice(max(i - 90, 0), i - 29); cellvals = {}
        for tag, (arr, cols) in cache.items():
            if (tag, e["key"]) in flagset: continue
            for k, f in enumerate(cols):
                c = fam2cell.get((tag, f))
                if c not in e["fp"]: continue
                a, b = arr[on, k], arr[pre, k]
                if np.isfinite(a).mean() > .6 and np.isfinite(b).mean() > .6: cellvals.setdefault(c, []).append(np.nanmean(a) - np.nanmean(b))
        if len(cellvals) >= 3: out.append(np.mean([np.mean(v) for v in cellvals.values()]))
        else: out.append(np.nan)
    return np.nanmean(out), int(np.nansum(np.array(out) < 0)), int(np.isfinite(out).sum())
m0, K0, N0 = gated_mean({}); nul = [gated_mean({g: int(rng.randint(90, span - 90)) for g in ("N", "C")})[0] for _ in range(300)]
p0 = float(np.mean(np.abs(nul) >= abs(m0)))
print(f"GATED pooled 3-5 s footprint mean = {m0:+.3f}% (K={K0}/{N0}), Null A p = {p0:.3f} (sd {np.std(nul):.3f})")
verdict = "LAPSE-CONSISTENT CANDIDATE" if (0.6 <= sl <= 1.2 and agree >= 0.7 and m0 <= -0.04 and p0 < 0.01) else "WAVEFORM-CHANGE CLASS -> BOUND, close"
print(f"\nRULE: slope in [0.6,1.2] & agreement >= 70% & gated mean <= -0.04% at p<0.01 -> candidate. RESULT: {verdict}")
pd.DataFrame(dict(key=[e["key"] for e in EP], c35_gated_fp=s.values)).to_csv(f"{SD}/lapse_gated_episodes35.csv", index=False)
