#!/usr/bin/env python
"""P4-0 step 4 - per station-day site-noise and registration-proxy statistics.

For every cached probe day:
  * VELOCITY PSD.  The 20 Hz cache is counts that have already been through the frozen chain's
    zero-phase Butterworth-4 1-8 Hz at the NATIVE rate, so the full instrument response |R(f)|
    (StationXML, evalresp, output=VEL, counts/(m/s); scalar instrument sensitivity as the
    fallback) is multiplied by |H(f)|^2 - squared because sosfiltfilt runs the filter forward AND
    backward - and the median-over-segments Welch PSD is divided by that amplitude squared.  The
    reported number is the MEDIAN over the 2-8 Hz bins, per component, in (m/s)^2/Hz.
    Segments containing a gap-fill zero are dropped, so a half-dead day is not reported as quiet.
  * INCOHERENT 6-s S-WINDOW ENERGY.  Median over the day of the mean square in non-overlapping
    6 s windows, per component (u_common.med_win_power).  No alignment and no stacking: this is
    the registration proxy, and its tremor-day / quiet-day ratio is scale-invariant, so a wrong
    sensitivity cannot manufacture a contrast.
  * COVERAGE = fraction of the day whose 3-C sample sum is non-zero.

DIP is honoured here: candidates.csv carries zsign = -1 for the stations whose vertical is
declared at dip +90 (positive DOWN); that sign is applied to Z on read, exactly as s_common.zsign
does for the P3D roster.  (It cannot change an energy, but the contract is kept and recorded.)

Output: results/p4/station_day_stats.csv
"""
from __future__ import annotations
import os
import sys
import time
import argparse
import warnings
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import u_common as U                                                   # noqa: E402

NPERSEG = 256                                    # 12.8 s at 20 Hz -> 0.078 Hz bins
_RESP = {}                                       # (netsta, chan, date) -> |R(f)| counts/(m/s)
_META = {}


def _init(resp, meta):
    global _RESP, _META
    _RESP, _META = resp, meta


def build_response_table(need, log):
    """need: DataFrame(net, sta, loc, band, date).  -> {(netsta, chan, date): |R(f)| array}."""
    from obspy import read_inventory, UTCDateTime
    f = np.fft.rfftfreq(NPERSEG, 1.0 / U.FS_CACHE)
    fq = f.copy()
    fq[0] = 1e-6
    out, mode = {}, {}
    xmldir = os.path.join(U.P4, "stationxml")
    invs = {}
    for p in sorted(os.listdir(xmldir)):
        if p.endswith(".xml"):
            invs[p[:-4]] = read_inventory(os.path.join(xmldir, p))
    chans = {}
    for name, inv in invs.items():
        for n in inv:
            for s in n:
                for ch in s:
                    chans.setdefault((f"{n.code}.{s.code}", ch.location_code or "", ch.code), []).append(ch)
    memo = {}
    for r in need.itertuples():
        ns = f"{r.net}.{r.sta}"
        for suf in ("Z", "N", "1", "E", "2"):
            key = (ns, (r.loc or ""), f"{r.band}{suf}")
            cand = chans.get(key)
            if not cand:
                continue
            t = UTCDateTime(r.date) + 43200
            ep = [c for c in cand
                  if (c.start_date is None or c.start_date <= t)
                  and (c.end_date is None or c.end_date >= t)]
            if not ep:
                continue
            ch = ep[0]
            mk = (key, str(ch.start_date))
            if mk not in memo:
                amp, how = None, "none"
                try:
                    R = ch.response.get_evalresp_response_for_frequencies(fq, output="VEL")
                    amp = np.abs(R)
                    if np.all(np.isfinite(amp[1:])) and np.nanmedian(amp[1:]) > 0:
                        how = "evalresp"
                    else:
                        amp = None
                except Exception:
                    amp = None
                if amp is None:
                    try:
                        s0 = float(ch.response.instrument_sensitivity.value)
                        u0 = str(ch.response.instrument_sensitivity.input_units or "").upper()
                        if u0.replace(" ", "") in ("M/S", "M/SEC", "M/S**1"):
                            amp = np.full(f.size, s0)
                            how = "sensitivity"
                    except Exception:
                        amp = None
                memo[mk] = (amp, how)
            amp, how = memo[mk]
            out[(ns, f"{r.band}{suf}", r.date)] = amp
            mode[(ns, f"{r.band}{suf}", r.date)] = how
    log(f"  [resp] {len(memo)} channel epochs evaluated; modes "
        f"{pd.Series([v[1] for v in memo.values()]).value_counts().to_dict()}")
    return out, mode, f


def one(job):
    netsta, date_s, y, doy, band, loc, fs_native, chain, zsign, day_type, role, site, sensor = job
    rec = dict(netsta=netsta, net=netsta.split(".")[0], sta=netsta.split(".")[1], date=date_s,
               day_type=day_type, role=role, site=site, sensor_type=sensor, band=band, loc=loc,
               fs_native=fs_native, chain=chain, zsign=zsign, err="")
    try:
        p = U.cache_path(netsta, y, doy)
        z = np.load(p)
        X = z["x"].astype(np.float64)
        order = [str(c) for c in z["order"]]
        X[0] *= float(zsign)                                    # dip honoured (positive-up Z)
        rec["coverage"] = float((np.abs(X).sum(0) > 0).mean())
        rec["order"] = ",".join(order)
        amp = []
        for c in order:
            a = _RESP.get((netsta, c, date_s))
            amp.append(a)
        if any(a is None for a in amp):
            rec["err"] = "no_response"
            amp = None
        else:
            Hsq = U.butter_band_amp(_META["freqs"], fs_native) if np.isfinite(fs_native) else 1.0
            amp = np.vstack([np.asarray(a, float) * Hsq for a in amp])
        bm, _, _, nseg = U.band_median_psd(X, U.FS_CACHE, *U.PSD_BAND, nperseg=NPERSEG,
                                           resp_amp=amp)
        rec["n_seg"] = nseg
        for i, lab in enumerate("ZNE"):
            rec[f"psd_{lab}"] = float(bm[i]) if i < len(bm) else np.nan
        rec["psd_3c"] = float(np.nanmedian(bm))
        # scalar-sensitivity velocity scaling for the (scale-invariant) contrast statistic
        if amp is not None:
            sc = np.nanmedian(amp[:, (_META["freqs"] >= 2) & (_META["freqs"] <= 8)], axis=1)
            Xv = X / np.maximum(sc[:, None], 1e-30)
        else:
            Xv = X
        e6, nwin = U.med_win_power(Xv, int(round(U.E_WIN_S * U.FS_CACHE)))
        rec["n_win"] = nwin
        for i, lab in enumerate("ZNE"):
            rec[f"e6_{lab}"] = float(e6[i]) if i < len(e6) else np.nan
        rec["e6_3c"] = float(np.nanmedian(e6))
        rec["status"] = "ok"
    except Exception as e:
        rec["status"] = "error"
        rec["err"] = f"{type(e).__name__}: {e}"
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=24)
    ap.add_argument("--out", default=os.path.join(U.P4, "station_day_stats.csv"))
    a = ap.parse_args()
    warnings.filterwarnings("ignore")

    def log(m):
        print(m, flush=True)

    Q = pd.read_csv(os.path.join(U.P4, "probe_queue.csv"), keep_default_na=False)
    PD = pd.read_csv(os.path.join(U.P4, "probe_days.csv"))
    dtype = dict(zip(PD.date, PD.day_type))
    C = pd.read_csv(os.path.join(U.P4, "candidates.csv"), keep_default_na=False)
    ref = os.path.join(U.P4, "reference_stations.csv")
    if os.path.exists(ref):
        RS = pd.read_csv(ref, keep_default_na=False)
        for c in ("site", "sensor_type"):
            if c not in RS:
                RS[c] = RS.sta if c == "site" else "reference"
        C = pd.concat([C, RS], ignore_index=True).drop_duplicates(["net", "sta"], keep="first")
    meta = {(r.net, r.sta): r for r in C.itertuples()}
    M = pd.read_csv(os.path.join(U.P4, "download_manifest.csv"), keep_default_na=False)
    M = M[M.status == "data"].drop_duplicates(["net", "sta", "date"], keep="last")
    fsn = {(r.net, r.sta, r.date): (float(r.fs) if str(r.fs) not in ("", "nan") else np.nan)
           for r in M.itertuples()}
    RESP, MODE, freqs = build_response_table(Q, log)
    jobs = []
    for r in Q.itertuples():
        ns = f"{r.net}.{r.sta}"
        d = pd.Timestamp(r.date)
        if not os.path.exists(U.cache_path(ns, d.year, int(d.dayofyear))):
            continue
        m = meta.get((r.net, r.sta))
        fs = fsn.get((r.net, r.sta, r.date), np.nan)
        if not np.isfinite(fs) and m is not None:
            fs = float(m.fs)
        chain = "frozen" if (np.isfinite(fs) and abs(fs / U.FS_CACHE - round(fs / U.FS_CACHE)) < 1e-9) else "rational"
        jobs.append((ns, r.date, d.year, int(d.dayofyear), r.band, (r.loc or ""), fs, chain,
                     float(getattr(m, "zsign", 1.0)) if m is not None else 1.0,
                     dtype.get(r.date, "?"), r.role,
                     str(getattr(m, "site", r.sta)) if m is not None else r.sta,
                     str(getattr(m, "sensor_type", "")) if m is not None else ""))
    log(f"[stats] {len(jobs)} cached station-days, {a.workers} workers")
    t0 = time.time()
    import multiprocessing as mp
    rows = []
    with mp.Pool(a.workers, initializer=_init,
                 initargs=(RESP, {"freqs": freqs}), maxtasksperchild=40) as pool:
        for i, r in enumerate(pool.imap_unordered(one, jobs, chunksize=1), 1):
            rows.append(r)
            if i % 100 == 0 or i == len(jobs):
                log(f"  {i}/{len(jobs)} {time.time()-t0:.0f}s")
    D = pd.DataFrame(rows)
    D["resp_mode"] = [MODE.get((r.netsta, f"{r.band}Z", r.date), "none") for r in D.itertuples()]
    D = D.sort_values(["netsta", "date"])
    D.to_csv(a.out, index=False)
    log(D.groupby(["status"]).size().to_string())
    for r in D[D.status == "error"].head(15).itertuples():
        log(f"  ERROR {r.netsta} {r.date}: {r.err}")
    log(f"[stats] -> {a.out} in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
