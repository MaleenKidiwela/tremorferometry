#!/usr/bin/env python
"""P4-0b seeds -- both arms, plus the seed-z band that TRUE hp events occupy.

  arm GIVEN  3 hp events (p_real >= 0.8) per eligible kept family, spread over episodes and days,
             ranked by the frozen catalogue's own net_mad (Bostock's construction starts from a LOUD
             window).  Moveout = the family's own cell.
  arm FREE   Merlin's CORRECTED seed rule (P4-2 step 2): 12-s windows at 6-s hop on PNSN tremor
             hours (>= 5 epicentres within 50 km of the region centre), episode days only; per-site
             3-C energy z at the REGION moveout, MAD-normalised per station-day on the LOWER HALF of
             the day's window distribution; >= 3 admitted sites individually z >= 3 and
             4 <= z_med <= 25; no catalogued earthquake within +-60 s; HF ratio E(10-19)/E(2-8)
             < 3.5 at the two nearest stations (NATIVE rate); <= 3 seeds/hour, >= 120 s apart;
             budget <= 1500, ordered by z_med within the band.

  python scripts/v_seeds.py --arm given|free|hpband
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import v_common as V                                                        # noqa: E402
import s_common as S                                                        # noqa: E402

OUT = V.OUT


# ------------------------------------------------------------------ energy z --------------------
def energy_z_day(date, stas, t_abs, cst, idx, win_s=V.SEED_WIN_S):
    """(len(t_abs), len(stas)) z of the 3-C energy in the win_s window at each station's moveout.

    Background = median and LOWHALF MAD of the day's own 1-s-stride window-energy distribution."""
    n = int(round(win_s * V.FS))
    Z = np.full((len(t_abs), len(stas)), np.nan, np.float32)
    y, doy = int(date.year), int(date.dayofyear)
    for si, s in enumerate(stas):
        v = S.load_day_zne(s, y, doy, idx)
        if v is None:
            continue
        t0, X = v
        e = (np.asarray(X, np.float64) ** 2).sum(0)
        cs = np.concatenate([[0.0], np.cumsum(e)])
        st = np.arange(0, len(e) - n, int(V.FS))
        bgw = cs[st + n] - cs[st]
        m = float(np.median(bgw))
        q = float(np.quantile(bgw, 0.25))
        mad = max(1.4826 * (m - q), 1e-9)
        i0 = np.rint((np.asarray(t_abs, float) + cst[s] - t0) * V.FS).astype(np.int64)
        ok = (i0 >= 0) & (i0 + n < len(e))
        w = np.full(len(t_abs), np.nan)
        w[ok] = cs[i0[ok] + n] - cs[i0[ok]]
        Z[:, si] = (w - m) / mad
    return Z


def site_stats(Z, stas):
    """(z_med over sites, n sites with z >= SEED_SITE_Z, n sites scored)."""
    sites = sorted({V.SITE_OF[s] for s in stas})
    P = np.full((Z.shape[0], len(sites)), np.nan, np.float32)
    for qi, q in enumerate(sites):
        cols = [i for i, s in enumerate(stas) if V.SITE_OF[s] == q]
        with np.errstate(all="ignore"):
            P[:, qi] = np.nanmean(Z[:, cols], axis=1)
    with np.errstate(all="ignore"):
        zmed = np.nanmedian(P, axis=1)
        n3 = np.nansum(P >= V.SEED_SITE_Z, axis=1)
        nsc = np.isfinite(P).sum(axis=1)
    return zmed, n3.astype(int), nsc.astype(int), P, sites


# ------------------------------------------------------------------ HF ratio ---------------------
def hf_ratio(netsta, date, t_abs, win_s=V.SEED_WIN_S, band=(10.0, 19.0)):
    """E(10-19)/E(2-8) over the window, from the NATIVE-rate day trace (o_run_month.hf_ratio_day)."""
    import o_common as C
    import o_chanorder as chanorder
    from scipy.signal import butter, sosfiltfilt
    out = np.full(len(t_abs), np.nan)
    p = C.trace_path(netsta, int(date.year), int(date.dayofyear))
    if p is None or not len(t_abs):
        return out, np.nan
    try:
        from obspy import read
        st = read(p)
        for tr in st:
            fr = round(float(tr.stats.sampling_rate))
            if fr > 0 and float(tr.stats.sampling_rate) != fr and \
               abs(float(tr.stats.sampling_rate) - fr) / fr < 1e-4:
                tr.stats.sampling_rate = float(fr)
        st.merge(method=1, fill_value=0)
        chans = sorted({tr.stats.channel for tr in st})
        order = chanorder.resolve_order(chans, station=netsta, strict=True)
        keep = [c for c in order if c is not None]
        npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in keep)
        X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in keep])
        fs = float(st[0].stats.sampling_rate)
        t0 = float(st.select(channel=keep[0])[0].stats.starttime.timestamp)

        def bp(x, f1, f2):
            return sosfiltfilt(butter(4, [f1 / (fs / 2), min(f2, 0.95 * fs / 2) / (fs / 2)],
                                      btype="band", output="sos"), x, axis=-1)
        lo = (bp(X, 2.0, 8.0) ** 2).sum(0)
        hi = (bp(X, band[0], band[1]) ** 2).sum(0)
        clo = np.concatenate([[0.0], np.cumsum(lo)])
        chi = np.concatenate([[0.0], np.cumsum(hi)])
        n = int(round(win_s * fs))
        i0 = np.rint((np.asarray(t_abs, float) - t0) * fs).astype(np.int64)
        g = (i0 >= 0) & (i0 + n <= len(lo))
        den = clo[i0[g] + n] - clo[i0[g]]
        num = chi[i0[g] + n] - chi[i0[g]]
        out[g] = np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)
        return out, fs
    except Exception as e:
        V.log(f"  [hf] {netsta} {date.date()} FAILED: {type(e).__name__}: {e}")
        return out, np.nan


# ------------------------------------------------------------------ arms ------------------------
def arm_given():
    keep, D, _ = V.eligible_families()
    A = V.availability()
    usable = {t for t, u in zip(A.date, A.usable) if u}
    ED = V.episode_days()
    D = D[[pd.Timestamp(d) in usable and pd.Timestamp(d) in ED for d in D.dt]]
    rows = []
    for b in keep:
        g = D[D.bkey == b].sort_values("net_mad", ascending=False)
        picked, days, eps = [], set(), set()
        for pref_new_ep in (True, False):
            for r in g.itertuples():
                if len(picked) >= V.GIVEN_SEEDS_PER_FAM:
                    break
                if r.dt in days:
                    continue
                if pref_new_ep and r.episode in eps:
                    continue
                picked.append(r); days.add(r.dt); eps.add(r.episode)
        for j, r in enumerate(picked):
            rows.append(dict(seed=f"G{b}_{j}", arm="given", bkey=int(b), t=float(r.t),
                             date=r.dt.strftime("%Y-%m-%d"), episode=r.episode,
                             net_mad=float(r.net_mad), p_real=float(r.p_real)))
    Sd = pd.DataFrame(rows)
    Sd.to_parquet(f"{OUT}/seeds_given.parquet", index=False)
    meta = dict(n_families=len(keep), families=keep, n_seeds=int(len(Sd)),
                per_family=Sd.groupby("bkey").size().to_dict(),
                rule=f"3 hp events per eligible kept family (>= {V.GIVEN_MIN_HP} hp in S18+S21+S22), "
                     "distinct days, distinct episodes preferred, ranked by frozen net_mad",
                n_hp_by_family=D.groupby("bkey").size().to_dict())
    json.dump(meta, open(f"{OUT}/seeds_given.json", "w"), indent=1, default=str)
    V.log(f"[given] {len(keep)} eligible families -> {len(Sd)} seeds")
    return Sd


def _candidate_windows():
    lat, lon, dep, _ = V.region_centre()
    tr = pd.read_csv(V.PNSN, usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    tr = tr.dropna(subset=["time"])
    tr = tr[S.epi_km(tr.lat.values, tr.lon.values, lat, lon) <= V.SEED_EPI_KM]
    ED = V.episode_days()
    A = V.availability()
    usable = {t for t, u in zip(A.date, A.usable) if u}
    ok_days = {t for t in ED if t in usable}
    tr = tr[[pd.Timestamp(t.date()) in ok_days for t in tr.time]]
    h = tr.time.dt.floor("h").value_counts()
    hrs = sorted(h[h >= V.SEED_MIN_EPI_HOUR].index)
    nper = int(round(3600.0 / V.SEED_HOP_S))
    C = []
    for hh in hrs:
        t0 = hh.value / 1e9
        C.append(pd.DataFrame(dict(t=t0 + np.arange(nper) * V.SEED_HOP_S,
                                   hour=hh, date=pd.Timestamp(hh.date()),
                                   episode=ED[pd.Timestamp(hh.date())])))
    C = pd.concat(C, ignore_index=True)
    return C.sort_values("t").reset_index(drop=True), len(hrs)


def arm_free():
    t00 = time.time()
    lat, lon, dep, dep_s2 = V.region_centre()
    cst = V.cell_tt(lat, lon, dep)
    C, nhours = _candidate_windows()
    V.log(f"[free] {nhours} tremor hours -> {len(C)} candidate windows")
    idx = V.cache_idx()
    stas = V.DET_STAS
    Z = np.full((len(C), len(stas)), np.nan, np.float32)
    for k, (d, g) in enumerate(C.groupby("date")):
        av = [s for s in stas if (s, int(d.year), int(d.dayofyear)) in idx]
        zz = energy_z_day(d, av, g.t.to_numpy(), cst, idx)
        for si, s in enumerate(av):
            Z[g.index.to_numpy(), stas.index(s)] = zz[:, si]
        if k % 10 == 0:
            V.log(f"  [free] day {k} {d.date()} ({time.time()-t00:.0f}s)")
    zmed, n3, nsc, P, sites = site_stats(Z, stas)
    C = C.assign(z_med=zmed, n_site_z3=n3, n_sites_scored=nsc)
    np.save(f"{OUT}/free_candidate_persite_z.npy", P.astype(np.float32))
    C.to_parquet(f"{OUT}/free_candidates.parquet", index=False)

    band = C[(C.n_site_z3 >= V.SEED_MIN_SITES) & (C.z_med >= V.SEED_ZMED_LO)
             & (C.z_med <= V.SEED_ZMED_HI)].copy()
    V.log(f"[free] in band: {len(band)} of {int(np.isfinite(C.z_med).sum())} scored")
    eq = pd.read_csv(V.EQCAT)
    eq["ot"] = pd.to_datetime(eq.ot, errors="coerce")
    teq = np.sort(eq.ot.dropna().astype("datetime64[ns]").astype("int64").to_numpy() / 1e9)
    j = np.searchsorted(teq, band.t.to_numpy())
    dmin = np.full(len(band), np.inf)
    for off in (-1, 0):
        k = np.clip(j + off, 0, len(teq) - 1)
        dmin = np.minimum(dmin, np.abs(band.t.to_numpy() - teq[k]))
    band = band[dmin > V.SEED_EQ_VETO_S].copy()
    V.log(f"[free] after +-{V.SEED_EQ_VETO_S:.0f}s earthquake veto: {len(band)}")

    band = band.sort_values("z_med", ascending=False).reset_index(drop=True)
    per_hour, last = {}, []
    keep = []
    for r in band.itertuples():
        if per_hour.get(r.hour, 0) >= V.SEED_MAX_PER_HOUR:
            continue
        if any(abs(r.t - u) < V.SEED_MIN_SEP_S for u in last):
            continue
        keep.append(r.Index)
        per_hour[r.hour] = per_hour.get(r.hour, 0) + 1
        last.append(r.t)
        if len(last) > 4000:
            last = last[-4000:]
    band = band.loc[keep].reset_index(drop=True)
    V.log(f"[free] after <=3/hour and >=120 s separation: {len(band)}")

    ST = V.station_table()
    dist = sorted(((float(S.epi_km(lat, lon, ST.lat[s], ST.lon[s])), s) for s in stas))
    near2 = [s for _, s in dist[:2]]
    V.log(f"[free] HF veto stations (nearest to centre): {near2}")
    cand = band.head(min(len(band), 3 * V.SEED_BUDGET)).copy()
    hf = np.full((len(cand), 2), np.nan)
    nat = {}
    for si, s in enumerate(near2):
        ns = str(ST.netsta[s])
        for d, g in cand.groupby("date"):
            r, fs = hf_ratio(ns, d, g.t.to_numpy() + cst[s])
            hf[g.index.to_numpy(), si] = r
            nat[s] = fs
    cand["hf1"], cand["hf2"] = hf[:, 0], hf[:, 1]
    with np.errstate(all="ignore"):
        hfmax = np.nanmax(hf, axis=1)
    cand["hf_max"] = hfmax
    cand["hf_native"] = all(np.isfinite(list(nat.values()))) if nat else False
    pre = len(cand)
    cand = cand[~(cand.hf_max >= V.SEED_HF_MAX)].copy()      # NaN (unmeasurable) is NEVER vetoed
    V.log(f"[free] after HF < {V.SEED_HF_MAX} at the 2 nearest stations: {len(cand)} of {pre}")
    Sd = cand.sort_values("z_med", ascending=False).head(V.SEED_BUDGET).copy()
    Sd["seed"] = [f"F{i}" for i in range(len(Sd))]
    Sd["arm"] = "free"
    Sd = Sd.sort_values("t").reset_index(drop=True)
    Sd.to_parquet(f"{OUT}/seeds_free.parquet", index=False)

    hh = Sd.hour.value_counts()
    dd = Sd.date.value_counts()
    meta = dict(n_tremor_hours=int(nhours), n_candidates=int(len(C)),
                n_scored=int(np.isfinite(C.z_med).sum()), n_in_band=int(len(band)),
                n_after_hf=int(len(cand)), n_seeds=int(len(Sd)),
                z_med_range=[float(Sd.z_med.min()), float(Sd.z_med.max())] if len(Sd) else None,
                max_seeds_per_hour=int(hh.max()) if len(hh) else 0,
                top10_hours_fraction=float(hh.head(10).sum() / max(len(Sd), 1)),
                top5_days_fraction=float(dd.head(5).sum() / max(len(Sd), 1)),
                n_distinct_hours=int(len(hh)), n_distinct_days=int(len(dd)),
                hf_stations=near2, hf_native_fs={k: float(v) for k, v in nat.items()},
                centre=dict(lat=lat, lon=lon, dep_used=dep, dep_slab2=dep_s2),
                cache_band="1-8 Hz (the frozen 20 Hz cache); the 2-8 Hz seed energy is therefore "
                           "measured in the cache's 1-8 Hz band -- FLAGGED",
                runtime_s=time.time() - t00)
    json.dump(meta, open(f"{OUT}/seeds_free.json", "w"), indent=1, default=str)
    V.log(f"[free] {len(Sd)} seeds; max/hour {meta['max_seeds_per_hour']}, "
          f"top-10-hours {meta['top10_hours_fraction']:.1%}")
    return Sd


def hp_band():
    """The seed-z band TRUE hp events occupy -- the number that fixes the Puget seed rule."""
    lat, lon, dep, _ = V.region_centre()
    cst = V.cell_tt(lat, lon, dep)
    keep, D, _ = V.eligible_families()
    hid, _, _ = S.holdouts()
    All = S.hp_events(min_p=0.8, episodes=V.EPISODES)
    All = All[~All.bkey.isin(hid)]
    idx = V.cache_idx()
    stas = V.DET_STAS
    Z = np.full((len(All), len(stas)), np.nan, np.float32)
    All = All.reset_index(drop=True)
    # the hp event time is on the family ORIGIN axis; the region-moveout window is placed at
    # t + cst[s], exactly as a template-free seed would be.
    for d, g in All.groupby("dt"):
        av = [s for s in stas if (s, int(d.year), int(d.dayofyear)) in idx]
        if not av:
            continue
        zz = energy_z_day(d, av, g.t.to_numpy(), cst, idx)
        for si, s in enumerate(av):
            Z[g.index.to_numpy(), stas.index(s)] = zz[:, si]
    zmed, n3, nsc, P, sites = site_stats(Z, stas)
    All = All.assign(z_med=zmed, n_site_z3=n3, n_sites_scored=nsc)
    All[["det_uid", "bkey", "t", "date", "episode", "p_real", "net_mad",
         "z_med", "n_site_z3", "n_sites_scored"]].to_parquet(f"{OUT}/hp_seed_z.parquet", index=False)
    ok = All[All.n_sites_scored >= V.SEED_MIN_SITES]
    q = np.nanpercentile(ok.z_med.dropna(), [5, 10, 25, 50, 75, 90, 95]) if len(ok) else []
    meta = dict(n_hp=int(len(All)), n_scored_ge3sites=int(len(ok)),
                z_med_percentiles=dict(zip([5, 10, 25, 50, 75, 90, 95], [float(x) for x in q])),
                frac_in_band=float(((ok.z_med >= V.SEED_ZMED_LO) & (ok.z_med <= V.SEED_ZMED_HI)
                                    & (ok.n_site_z3 >= V.SEED_MIN_SITES)).mean()) if len(ok) else None,
                frac_zmed_below_4=float((ok.z_med < 4).mean()) if len(ok) else None,
                frac_zmed_above_25=float((ok.z_med > 25).mean()) if len(ok) else None,
                frac_ge3_sites_z3=float((ok.n_site_z3 >= 3).mean()) if len(ok) else None,
                note="median-over-sites z of the 12-s region-moveout energy window at TRUE hp events")
    json.dump(meta, open(f"{OUT}/hp_seed_z.json", "w"), indent=1, default=str)
    V.log(f"[hpband] {len(ok)} hp events scored; z_med median "
          f"{meta['z_med_percentiles'].get(50)}; in-band {meta['frac_in_band']}")
    return All


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", required=True, choices=["given", "free", "hpband"])
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    {"given": arm_given, "free": arm_free, "hpband": hp_band}[a.arm]()


if __name__ == "__main__":
    main()
