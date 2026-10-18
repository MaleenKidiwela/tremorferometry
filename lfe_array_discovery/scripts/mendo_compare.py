#!/usr/bin/env python
"""MENDO (2026-09-29): real families of the blind Royer-chain run vs Shelly et al. (2025)'s 27 families.

Real family  = auto_score.csv n_sites_3km >= 3 (>= 3 independent sites with an impulsive 2-8 Hz arrival,
               stations < 3 km apart = one site; ab_nomerge_score.py), as NW_lag15 / SW_lag15 / NVI_lag15.
Three "known" tests (a family is NEW only if it fails all three):
 (a) LOCATION  headline Hyp2000 epicentre (hypo/loc_headline.csv) within LOC_KM (10 km, the nw_union.py
               candidate radius) of any Shelly Table S1 family (lfe_south/results/families_south.csv =
               the Table S1 / data-release locations).  Unlocated families cannot pass (a); the depth
               difference is a column.
 (b) WAVEFORM  duplicate of any Shelly v3b it0 template (lfe_south/results/sh_v3b/tpl_it0.npz, 6 s,
               20 sps -> resampled to 40 sps) by the nw_union.py rule: >= 2 common stations with
               cc >= 0.5 on 2-8 Hz (same channel letter Z/N/E; PB 1/2 horizontals are unrotated in v3b,
               so only Z there), and spread <= 0.15 s of the implied origin offset
               (lag of the template inside our 44 s 2-8 Hz stack minus the template's own start offset
               tt_s + best_tau - pre_s relative to Shelly's catalogue time).
 (c) TIMES     shared detection times with the Shelly catalogue (LFEcat_MTJ_2018-2024_shelly2025.csv, all
               27 families): our detection times t vs his, differences within +-SEARCH_S; the family
               offset is unknown (different time references), so the peak 4 s bin (+-2 s) of the
               difference histogram is taken; shared if (peak - median bin) >= 0.10 of our detections AND
               peak >= 2 x the median bin (median bin = background coincidence rate, which is high
               during bursts: a synthetic copy of his family 5 reads peak 1.0 / background 0.16).
-> <run>/compare_shelly/families_vs_shelly.csv (+ printed counts)
  python scripts/mendo_compare.py results/ab_royer/MENDO
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.signal import butter, correlate, resample_poly, sosfiltfilt

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ab_common as A  # noqa: E402

SOUTH = f"{A.ROOT}/lfe_south"
CAT = f"{A.ROOT}/lfe_catalog_repair/results/external/other/LFEcat_MTJ_2018-2024_shelly2025.csv"
LOC_KM = 10.0
SEARCH_S = 20.0
PRE = 12.0
SOS = butter(4, [2.0 / 20.0, 8.0 / 20.0], "band", output="sos")


def km(a, b, c, d):
    a, b, c, d = map(np.radians, (a, b, c, d))
    q = np.sin((c - a) / 2) ** 2 + np.cos(a) * np.cos(c) * np.sin((d - b) / 2) ** 2
    return 6371 * 2 * np.arcsin(np.sqrt(q))


def load_ours(npz):
    z = np.load(npz)
    S = {}
    for k in z.files:
        g, s = k.split("_", 1)
        S.setdefault(int(g), {})[s] = z[k].astype(float)
    return S


def load_shelly_tpl():
    T = pd.read_csv(f"{SOUTH}/results/sh_v3b/tpl_it0.csv")
    T = T[(T.source == "Shelly2025") & (T.status == "ok")]
    z = np.load(f"{SOUTH}/results/sh_v3b/tpl_it0.npz", allow_pickle=True)
    D = {}
    for r in T.itertuples():
        if r.key not in z.files:
            continue
        x = z[r.key].astype(float)
        x = resample_poly(x, 2, 1, axis=-1)
        x = sosfiltfilt(SOS, x, axis=-1)
        chans = str(r.chans).split(",")
        comp = {}
        for i, c in enumerate(chans[:x.shape[0]]):
            L = c[-1]
            if L in "ZNE":
                comp["ZNE".index(L)] = x[i]
        off = float(r.tt_s) + float(r.best_tau) - float(r.pre_s)
        D.setdefault(int(r.fam), {})[r.sta] = (comp, off)
    return D


def wave_dup(ours, his):
    """nw_union.py rule adapted to a 6 s template inside our 44 s stack."""
    est, ccs = [], []
    for s in set(ours) & set(his):
        comp, off = his[s]
        best = (0.0, np.nan)
        for c, y in comp.items():
            x = ours[s][c]
            if not (x ** 2).sum() > 0 or not (y ** 2).sum() > 0:
                continue
            n = len(y)
            num = correlate(x, y, "valid")
            e = np.convolve(x * x, np.ones(n), "valid")
            cc = num / np.sqrt(np.maximum(e, 1e-30) * (y @ y))
            k = int(np.argmax(cc))
            if cc[k] > best[0]:
                best = (float(cc[k]), k / A.FS - off)
        if best[0] >= 0.5:
            ccs.append(best[0]); est.append(best[1])
    spread = float(np.ptp(est)) if len(est) > 1 else np.nan
    return (len(ccs) >= 2 and spread <= 0.15), len(ccs), spread


def main():
    run = sys.argv[1]
    out = f"{run}/compare_shelly"
    os.makedirs(out, exist_ok=True)
    sc = pd.read_csv(f"{run}/hypo_28/review/auto_score.csv")
    real = sorted(sc[sc.n_sites_3km >= 3].grp)
    L = pd.read_csv(f"{run}/hypo/loc_headline.csv").set_index("grp")
    det = pd.read_parquet(sorted([f"{run}/grow/{f}" for f in os.listdir(f"{run}/grow")
                                  if f.startswith("det_it")],
                                 key=lambda p: int(p.split("det_it")[1].split(".")[0]))[-1])
    F = pd.read_csv(f"{SOUTH}/results/families_south.csv")
    F = F[F.source == "Shelly2025"]
    C = pd.read_csv(CAT)
    C.columns = [c.strip("%") for c in C.columns]
    ts = (pd.to_datetime(C[["year", "month", "day"]]).astype("datetime64[ns]").astype("int64") / 1e9
          + C.s_of_day.astype(float)).to_numpy(copy=True)
    ts = np.sort(ts)
    ours = load_ours(f"{run}/hypo_28/stacks_loc.npz")
    his = load_shelly_tpl()
    rows = []
    for g in real:
        r = dict(grp=g, n_sites_3km=int(sc.set_index("grp").n_sites_3km[g]))
        # (a)
        if g in L.index and np.isfinite(L.lat.get(g, np.nan)):
            d = km(L.lat[g], L.lon[g], F.lat.to_numpy(), F.lon.to_numpy())
            j = int(np.argmin(d))
            r.update(lat=L.lat[g], lon=L.lon[g], dep=L.dep[g], erh=L.erh[g], erz=L.erz[g],
                     boot_km=L.get("boot_km", pd.Series(dtype=float)).get(g, np.nan),
                     near_fam=int(F.fam.iloc[j]), near_km=float(d[j]),
                     near_ddep=float(L.dep[g] - F.depth.iloc[j]), a_loc=bool(d[j] <= LOC_KM))
        else:
            r.update(a_loc=False)
        # (b)
        hit = ("", 0, np.nan)
        for f_, H in his.items():
            ok, n, sp = wave_dup(ours.get(g, {}), H)
            if ok:
                hit = (f_, n, sp)
                break
        r.update(b_wave=bool(hit[0] != ""), b_fam=hit[0], b_nsta=hit[1], b_spread=hit[2])
        # (c)
        t = np.sort(det[det.grp == g].t.to_numpy(float))
        diffs = []
        for x in t:
            i0, i1 = np.searchsorted(ts, [x - SEARCH_S, x + SEARCH_S])
            diffs += list(ts[i0:i1] - x)
        h, _ = np.histogram(diffs, bins=np.arange(-SEARCH_S, SEARCH_S + 1e-9, 1.0))
        h4 = np.convolve(h, np.ones(4), "valid")
        pk = float(h4.max()) / max(len(t), 1) if len(h4) else 0.0
        bg = float(np.median(h4)) / max(len(t), 1) if len(h4) else 0.0
        r.update(n_det=len(t), c_peak_frac=pk, c_bg_frac=bg,
                 c_times=bool(pk - bg >= 0.10 and pk >= 2 * max(bg, 1e-9)))
        r["known"] = bool(r["a_loc"] or r["b_wave"] or r["c_times"])
        r["new"] = not r["known"]
        rows.append(r)
    R = pd.DataFrame(rows)
    R.to_csv(f"{out}/families_vs_shelly.csv", index=False)
    if len(R):
        print(f"real {len(R)}; located {int(R.get('lat', pd.Series(dtype=float)).notna().sum())}; "
              f"(a) loc<= {LOC_KM:g} km {int(R.a_loc.sum())}; (b) waveform dup {int(R.b_wave.sum())}; "
              f"(c) shared times {int(R.c_times.sum())}; NEW {int(R.new.sum())} "
              f"(located new {int((R.new & R.get('lat', pd.Series(np.nan, index=R.index)).notna()).sum())})")
    else:
        print("no real families")


if __name__ == "__main__":
    main()
