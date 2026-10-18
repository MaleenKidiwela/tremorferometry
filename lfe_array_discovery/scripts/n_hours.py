#!/usr/bin/env python
"""STAGE 1 -- the 2008 (and 2007) TREMOR-HOUR AUDIT.  Prefix `n_`.

There is no tremor catalogue before 2010-01-06, and the hour flag that the 2008 site probe used
(network-median hourly energy over a fixed p95) is **86 % daytime** (`results/WA_V2_INPUTS.md`
sec.2.2) -- i.e. it is a cultural-noise detector.  This module replaces it with a statistic that
cultural noise cannot fake, in three parts:

  (1) DIURNAL CYCLE REMOVED PER STATION.  For every admitted station s, every day d and every
      UTC hour h, E(s,d,h) = median over the hour of the 6-s-window mean square of the 2-8 Hz
      3-C waveform (median over components).  The reference is that STATION's OWN median of the
      same quantity at the SAME hour-of-day over the literature-quiet days, M(s,h).  The ratio
      r = E/M therefore has the station's site response, its instrument gain AND its own diurnal
      cycle divided out.  A station that is simply loud at 16 UTC every day has r ~ 1 at 16 UTC.

  (2) SIMULTANEITY ACROSS THE ARRAY.  An hour is a CANDIDATE when at least `--min-frac` (0.60) of
      the reporting near-field stations -- admitted stations within 60 km of an active cell of the
      sub-region, health-veto applied -- are simultaneously above THEIR OWN p95 of r, where that
      p95 is taken over the same station's literature-quiet days.  Per-station thresholds mean a
      noisy station cannot carry the hour and a quiet one cannot veto it.

  (3) SPATIAL COHERENCE.  Tremor is a source at depth: the excess must DECAY with distance from
      one place.  Cultural noise is either site-local or basin-wide.  For each candidate hour the
      log-excess u_s = log10 r(s,h) of the reporting stations is fitted, over every active cell c,
      with u = a - b log10(max(d_s(c), 10 km)), b >= 0, and compared with the uniform model u = a.
      The statistic is the variance reduction VR of the best cell.  Its null is the SAME fit with
      the STATION LABELS SHUFFLED (`--n-shuf`, default 200), which preserves both the geometry and
      the distribution of u and destroys only their pairing.  An hour passes when its VR beats the
      shuffled distribution at p <= `--p-coh` (0.05).

A TREMOR HOUR = (2) AND (3).  Both the (2)-only and the (2)+(3) sets are reported with their
UTC-hour and weekday/weekend histograms, so the diurnal question is answered with numbers.

For the 2007 confirmation bed there is no literature-quiet day inside the window, but tranche T5
carries 15 out-of-episode days (2006-11-xx, 2007-04-xx) which are used as M(s,h)'s day set in
exactly the same way.

  python scripts/n_hours.py --stage energy        # per station-day hourly 2-8 Hz energy
  python scripts/n_hours.py --stage audit         # the audit itself + histograms
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import n_common as N                                                        # noqa: E402

FS = 20.0
WIN_S = 6.0
FMIN, FMAX = 2.0, 8.0
MIN_REPORT = 20          # an hour is only judged when this many near stations report
MIN_FRAC = 0.60          # ... and this fraction of them must be above their own p95
P_COH = 0.05
N_SHUF = 200
D_FLOOR_KM = 10.0        # the near-field flattening of any amplitude-decay law


# ----------------------------------------------------------------- (1) the energy table --------
def one_file(job):
    """(netsta, date, 24 hourly 2-8 Hz energies) for one cached station-day."""
    path, netsta, date = job
    from scipy.signal import butter, sosfiltfilt
    try:
        z = np.load(path)
        x = np.asarray(z["x"], np.float64)
    except Exception as e:                                                  # noqa: BLE001
        return dict(netsta=netsta, date=date, status=f"err:{type(e).__name__}")
    if x.ndim != 2 or x.shape[0] != 3:
        return dict(netsta=netsta, date=date, status="not3c")
    sos = butter(4, [FMIN / (FS / 2), min(FMAX, FS / 2 - 0.1) / (FS / 2)],
                 btype="band", output="sos")
    bad = ~np.isfinite(x)
    if bad.any():
        x = np.where(bad, 0.0, x)
    y = sosfiltfilt(sos, x, axis=-1)
    W = int(round(WIN_S * FS))
    per_h = int(round(3600 * FS)) // W          # 600 windows per hour
    nw = min(y.shape[1] // W, 24 * per_h)
    p = (y[:, :nw * W] ** 2).reshape(3, nw, W).mean(-1)
    out = dict(netsta=netsta, date=date, status="ok")
    hh = np.arange(nw) // per_h
    for h in range(24):
        m = hh == h
        if m.sum() < per_h // 2:
            out[f"h{h:02d}"] = np.nan
            continue
        v = np.median(p[:, m], axis=1)          # per component, median over the hour
        v = v[v > 0]
        out[f"h{h:02d}"] = float(np.median(v)) if v.size else np.nan
    return out


def stage_energy(a):
    os.makedirs(N.OUTROOT, exist_ok=True)
    idx = N.cache_index()
    ns = N.netsta()
    rev = {v: k for k, v in ns.items()}
    want = set(N.bed_days())
    jobs = []
    for (sta, yr, doy), p in idx.items():
        d = (pd.Timestamp(year=yr, month=1, day=1)
             + pd.Timedelta(days=doy - 1)).strftime("%Y-%m-%d")
        if d in want:
            jobs.append((p, ns[sta], d))
    jobs.sort(key=lambda j: (j[2], j[1]))
    if a.days:
        keep = set(a.days.split(","))
        jobs = [j for j in jobs if j[2] in keep]
    N.log(f"[energy {N.BED}] {len(jobs)} station-days ({len(rev)} stations), "
          f"{a.workers} workers -> {a.out}")
    rows = []
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for k, r in enumerate(ex.map(one_file, jobs, chunksize=4), 1):
            rows.append(r)
            if k % 500 == 0:
                N.log(f"   {k}/{len(jobs)}")
    E = pd.DataFrame(rows)
    E.to_parquet(a.out, index=False)
    N.log(f"[energy {N.BED}] {len(E)} rows ({int((E.status == 'ok').sum())} ok) -> {a.out}")


# ----------------------------------------------------------------- (2)+(3) the audit -----------
def ratio_table(E, quiet):
    """(long DataFrame date, hour, sta, r) plus {sta: p95 of r on the quiet days}."""
    hc = [f"h{h:02d}" for h in range(24)]
    E = E[E.status == "ok"].copy()
    L = E.melt(id_vars=["netsta", "date"], value_vars=hc, var_name="hh", value_name="e")
    L["hour"] = L.hh.str[1:].astype(int)
    L = L.drop(columns=["hh"]).dropna(subset=["e"])
    L = L[L.e > 0]
    Q = L[L.date.isin(set(quiet))]
    M = Q.groupby(["netsta", "hour"]).e.median().rename("m").reset_index()
    L = L.merge(M, on=["netsta", "hour"], how="inner")
    L["r"] = L.e / L.m
    P = (L[L.date.isin(set(quiet))].groupby("netsta").r
         .quantile(0.95).rename("p95").reset_index())
    L = L.merge(P, on="netsta", how="inner")
    return L, M, P


def coherence(u, dkm, n_shuf, rng):
    """Best single-cell amplitude-decay fit against the uniform model.

    u    (n_sta,)          log10 excess
    dkm  (n_cell, n_sta)   epicentral distance station-cell
    -> (VR, p, best cell index, slope b)"""
    X = np.log10(np.maximum(dkm, D_FLOOR_KM))                  # (n_cell, n_sta)
    n = u.size
    ss_tot = float(((u - u.mean()) ** 2).sum())
    if ss_tot <= 0:
        return 0.0, 1.0, -1, 0.0

    def best(v):
        xm = X.mean(1, keepdims=True)
        vm = v.mean()
        sxy = ((X - xm) * (v - vm)).sum(1)
        sxx = ((X - xm) ** 2).sum(1)
        b = np.where(sxx > 0, sxy / np.maximum(sxx, 1e-12), 0.0)
        b = np.minimum(b, 0.0)                                  # u must DECAY: slope <= 0 in log d
        res = (v[None, :] - vm) - b[:, None] * (X - xm)
        ss = (res ** 2).sum(1)
        i = int(np.argmin(ss))
        return 1.0 - float(ss[i]) / ss_tot, i, float(-b[i])

    vr, i, bb = best(u)
    nul = np.array([best(u[rng.permutation(n)])[0] for _ in range(n_shuf)])
    return vr, float((nul >= vr).mean()), i, bb


def stage_audit(a):
    import s_common as S
    os.makedirs(N.OUTROOT, exist_ok=True)
    E = pd.read_parquet(a.energy)
    quiet = N.quiet_days()
    near_sta = N.near_stations(a.near_km)
    ns = N.netsta()
    near_ns = {ns[s] for s in near_sta}
    ST = N.station_table()
    AC = N.active_cells()
    veto = N.health_veto()
    N.log(f"[audit {N.BED}] {len(near_sta)} near-field stations (<= {a.near_km:.0f} km of an "
          f"active cell) of {len(N.stations())} admitted; {len(quiet)} literature-quiet days; "
          f"{len(AC)} active cells")

    L, M, P = ratio_table(E, quiet)
    L = L[L.netsta.isin(near_ns)].copy()
    L["sta"] = [n.split(".")[1] for n in L.netsta]
    # health veto (per station-hour) applied to the audit exactly as to the pair stage
    L["ok"] = [(r.sta not in veto.get((r.date, int(r.hour)), set())) for r in L.itertuples()]
    L = L[L.ok].drop(columns=["ok"])
    L["above"] = L.r > L.p95
    P.to_csv(f"{N.OUTROOT}/station_p95_{N.BED}.csv", index=False)

    cla, clo = AC.lat.to_numpy(), AC.lon.to_numpy()
    DK = {}
    for s in near_sta:
        DK[s] = S.epi_km(float(ST.lat[s]), float(ST.lon[s]), cla, clo)

    rng = np.random.default_rng(7)
    rows = []
    for (date, hour), g in L.groupby(["date", "hour"]):
        nrep = len(g)
        nab = int(g.above.sum())
        frac = nab / max(nrep, 1)
        rec = dict(date=date, hour=int(hour), n_report=nrep, n_above=nab, frac_above=frac,
                   med_r=float(g.r.median()), quiet_day=bool(date in set(quiet)))
        cand = (nrep >= a.min_report) and (frac >= a.min_frac)
        rec["clause2"] = bool(cand)
        if cand:
            u = np.log10(g.r.to_numpy())
            dk = np.vstack([DK[s] for s in g.sta]).T          # (n_cell, n_sta)
            vr, p, ic, bb = coherence(u, dk, a.n_shuf, rng)
            rec.update(vr=vr, p_coh=p, slope_b=bb,
                       fit_lat=float(AC.lat.iloc[ic]) if ic >= 0 else np.nan,
                       fit_lon=float(AC.lon.iloc[ic]) if ic >= 0 else np.nan)
            rec["clause3"] = bool(p <= a.p_coh)
        else:
            rec.update(vr=np.nan, p_coh=np.nan, slope_b=np.nan,
                       fit_lat=np.nan, fit_lon=np.nan, clause3=False)
        rec["tremor_hour"] = bool(rec["clause2"] and rec["clause3"])
        rows.append(rec)
    A = pd.DataFrame(rows).sort_values(["date", "hour"]).reset_index(drop=True)
    A.to_csv(f"{N.OUTROOT}/hour_audit_{N.BED}.csv", index=False)

    W = A[~A.quiet_day]                                   # in-window hours only for selection
    sel = W[W.tremor_hour]
    fallback = ""
    if len(sel) < a.min_hours:
        alt = W[W.clause2].sort_values("p_coh")
        take = alt[alt.p_coh <= a.p_coh_relaxed]
        fallback = (f"only {len(sel)} hours pass clause 2 AND clause 3 at p<={a.p_coh}; "
                    f"the relaxed coherence bar p<={a.p_coh_relaxed} is used and RECORDED")
        N.log(f"[audit {N.BED}] FALLBACK: {fallback}")
        sel = take
    sel = sel.sort_values(["date", "hour"])
    sel[["date", "hour"]].to_csv(f"{N.OUTROOT}/tremor_hours_{N.BED}.csv", index=False)
    open(f"{N.OUTROOT}/tremor_dates_{N.BED}.txt", "w").write(",".join(sorted(set(sel.date))))
    open(f"{N.OUTROOT}/quiet_dates_{N.BED}.txt", "w").write(",".join(quiet))

    def hist(D):
        wk = pd.to_datetime(D.date).dt.weekday if len(D) else pd.Series(dtype=int)
        return dict(n=int(len(D)), n_days=int(D.date.nunique()) if len(D) else 0,
                    by_utc_hour={int(h): int(v) for h, v in
                                 D.hour.value_counts().sort_index().items()},
                    frac_daytime_15_22utc=float(D.hour.between(15, 22).mean()) if len(D) else 0.0,
                    frac_weekday=float((wk < 5).mean()) if len(D) else 0.0,
                    dates=sorted(set(D.date)))

    # SENSITIVITY DIAGNOSTIC.  The per-station p95 is pooled over all quiet-day hours, so the
    # day-to-day spread of r -- which is larger in the cultural daytime -- makes a daytime hour
    # harder to trigger than a night hour.  The quiet-day frac_above BY UTC HOUR measures exactly
    # that asymmetry: if it is flat, the selected hours' clock distribution is a property of the
    # signal; if it is night-heavy, part of it is the threshold.  Reported either way.
    QD = A[A.quiet_day]
    q_by_hour = {int(h): float(v) for h, v in QD.groupby("hour").frac_above.mean().items()}
    # ... and the array's own DIURNAL NOISE FLOOR, the network median of M(s,h) by UTC hour,
    # normalised to its minimum.  A fixed tremor amplitude produces a SMALLER ratio r when the
    # floor is high, so the clock distribution of the selected hours is bounded from below by
    # this curve: where the floor is 5x, the detectable tremor is 5x larger.  The selected hours
    # are therefore a DETECTION-LIMITED sample of the episode, not a claim that tremor is
    # nocturnal.
    Mn = M[M.netsta.isin(near_ns)].groupby("hour").m.median()
    floor = {int(h): float(v / Mn.min()) for h, v in Mn.items()}
    rep = dict(bed=N.BED, window=[N.DAY0, N.DAY1],
               quiet_frac_above_by_utc_hour=q_by_hour,
               diurnal_noise_floor_by_utc_hour=floor,
               quiet_frac_above_night_00_12=float(QD[QD.hour <= 12].frac_above.mean())
               if len(QD) else np.nan,
               quiet_frac_above_day_13_23=float(QD[QD.hour >= 13].frac_above.mean())
               if len(QD) else np.nan,
               quiet_hours_clause2=int(QD.clause2.sum()) if len(QD) else 0,
               n_near_stations=len(near_sta), near_stations=near_sta,
               n_quiet_days=len(quiet), quiet_days=quiet,
               min_report=a.min_report, min_frac=a.min_frac, p_coh=a.p_coh, n_shuf=a.n_shuf,
               n_hours_judged=int(len(W)),
               clause2_only=hist(W[W.clause2]),
               clause2_and_3=hist(W[W.tremor_hour]),
               selected=hist(sel), fallback=fallback,
               quiet_day_hours_passing_clause2=int(A[A.quiet_day & A.clause2].shape[0]),
               quiet_day_hours_passing_both=int(A[A.quiet_day & A.tremor_hour].shape[0]),
               n_quiet_day_hours_judged=int(A[A.quiet_day].shape[0]))
    json.dump(rep, open(f"{N.OUTROOT}/hour_audit_{N.BED}.json", "w"), indent=1, default=float)
    N.log(json.dumps({k: v for k, v in rep.items()
                      if k not in ("near_stations", "quiet_days")}, indent=1, default=float))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["energy", "audit"], required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--energy", default="")
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--days", default="")
    ap.add_argument("--near-km", type=float, default=60.0)
    ap.add_argument("--min-report", type=int, default=MIN_REPORT)
    ap.add_argument("--min-frac", type=float, default=MIN_FRAC)
    ap.add_argument("--p-coh", type=float, default=P_COH)
    ap.add_argument("--p-coh-relaxed", type=float, default=0.20)
    ap.add_argument("--min-hours", type=int, default=24)
    ap.add_argument("--n-shuf", type=int, default=N_SHUF)
    a = ap.parse_args()
    a.out = a.out or f"{N.OUTROOT}/hour_energy28_{N.BED}.parquet"
    a.energy = a.energy or f"{N.OUTROOT}/hour_energy28_{N.BED}.parquet"
    if a.stage == "energy":
        stage_energy(a)
    else:
        stage_audit(a)


if __name__ == "__main__":
    main()
