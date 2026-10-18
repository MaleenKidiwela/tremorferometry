#!/usr/bin/env python
"""DUC_FIX data stage -- plan, fetch and cache one MONTH of HIS stations with the epoch-specific
channel group (results/duc_fix/station_epochs.csv).

  plan   every (station, day) inside an epoch of a station his families_permanent.txt names
  fetch  reuse traces/<NET.STA>/<yyyy>/<doy>.mseed (shared, all channels) when present, or its
         .nd marker; else download ONLY the epoch's loc + channels into traces_duc/ (<= 8 threads)
  cache  Ducellier's chain verbatim (t6_duc_cache.py = his get_data.py): detrend, 5 s hann taper,
         remove_response VEL pre_filt (0.2,0.5,10,15) water_level 80, bandpass 1.5-9 zerophase,
         lanczos a=10 interpolate + decimate by Fraction(fs, 20), per native rate; VALID mask of the
         pre-fill trace segments.  Two speed-ups that do not change a sample: the inventory is read
         once per worker, and evalresp is memoised per (response, dt, nfft) (obspy's remove_response
         inverts the returned array in place, so a copy is handed out).
  python scripts/duc_fix_data.py --month 2012-03 [--workers 10] [--fetch-workers 8]
"""
import os, sys, time, argparse, warnings, threading
from fractions import Fraction
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from duc_fix_common import *

FILT = (1.5, 9.0)


def month_plan(month):
    E = epochs(); fam = families()
    used = set(",".join(fam.stations).split(","))
    E = E[E.sta.isin(used)]
    m0 = pd.Timestamp(month + "-01"); m1 = m0 + pd.offsets.MonthBegin(1)
    rows = []
    for r in E.itertuples():
        a = max(r.start, m0); b = min(r.end, m1)
        for d in pd.date_range(a, b, freq="D", inclusive="left"):
            rows.append((r.net, r.sta, str(d.date()), r.loc, r.chans, r.epoch, r.status))
    return pd.DataFrame(rows, columns=["net", "sta", "day", "loc", "chans", "epoch", "status"])


# ------------------------------------------------------------------------------------------ fetch
_cl = {}; _lk = threading.Lock()
def _client(prov):
    from obspy.clients.fdsn import Client
    k = (prov, threading.get_ident())
    if k not in _cl: _cl[k] = Client(prov, timeout=180)
    return _cl[k]


def src_path(net, sta, day):
    p = day_path(TR, net, sta, day)
    if os.path.exists(p): return p
    p2 = day_path(TR2, net, sta, day)
    if os.path.exists(p2): return p2
    return None


def fetch_one(job):
    net, sta, day, loc, chans = job
    if src_path(net, sta, day): return "have", 0
    p = day_path(TR2, net, sta, day)
    if os.path.exists(p + ".nd") or os.path.exists(day_path(TR, net, sta, day) + ".nd"): return "nodata", 0
    from obspy import UTCDateTime
    t0 = UTCDateTime(day); msg = ""
    for att in range(3):
        try:
            st = _client(PROV.get(net, "IRIS")).get_waveforms(net, sta, loc if loc else "--", chans, t0 - 60, t0 + 86400 + 60)
            if len(st) == 0: raise ValueError("empty")
            os.makedirs(os.path.dirname(p), exist_ok=True)
            tmp = p + f".tmp{os.getpid()}_{threading.get_ident()}"
            st.write(tmp, format="MSEED"); os.replace(tmp, p)
            return "data", os.path.getsize(p)
        except Exception as e:
            msg = str(e)[:120]
            if "No data" in msg or "empty" in msg or "204" in msg:
                os.makedirs(os.path.dirname(p), exist_ok=True); open(p + ".nd", "w").close()
                return "nodata", 0
            time.sleep(3 * (att + 1))
    return "error:" + msg, 0


# ------------------------------------------------------------------------------------------ cache
_INV = {}; _ER = {}
PREF4 = 15.0          # his pre_filt (0.2, 0.5, 10, 15): the taper is exactly 0 at f >= 15 Hz
def _patch_evalresp():
    """memoise evalresp per (response, dt, nfft) and evaluate it only where his pre_filt taper is
    non-zero (f < 15 Hz).  Above 15 Hz the data spectrum is multiplied by an exact 0, so the value
    there is irrelevant; it is filled with the 15 Hz value so the water-level maximum is unaffected
    in practice (water level 80 dB only clips bins < 1e-4 x max, none of which are in 0.2-15 Hz).
    evalresp evaluates every frequency independently, so the in-band values are the same numbers."""
    from obspy.core.inventory.response import Response
    if getattr(Response, "_duc_memo", False): return
    def memo(self, t_samp, nfft, output="VEL", start_stage=None, end_stage=None, **kw):
        k = (id(self), float(t_samp), int(nfft), output, start_stage, end_stage)
        if k not in _ER:
            if len(_ER) > 8: _ER.pop(next(iter(_ER)))
            fy = 1 / (t_samp * 2.0)
            freqs = np.linspace(0, fy, int(nfft // 2) + 1, dtype=np.float64)
            n = int(np.searchsorted(freqs, PREF4 + 0.5)) + 1
            n = min(n, freqs.size)
            r = self.get_evalresp_response_for_frequencies(freqs[:n], output=output, start_stage=start_stage,
                                                           end_stage=end_stage, **kw)
            R = np.empty(freqs.size, dtype=np.complex128); R[:n] = r; R[n:] = r[-1]
            _ER[k] = (R, freqs)
        r, f = _ER[k]
        return r.copy(), f.copy()
    Response.get_evalresp_response = memo; Response._duc_memo = True


_GPUFFT = None
def _patch_fft():
    """remove_response's day-long rfft/irfft (~17 M points, ~1 s each on one CPU core) on the GPU (cuFFT,
    double precision; identical to round-off, then stored as float32).  Falls back to numpy."""
    global _GPUFFT
    if _GPUFFT is not None or os.environ.get("DUC_GPUFFT", "1") != "1": return
    try:
        import cupy as cp
        cp.get_default_memory_pool().set_limit(size=int(1.5 * 1024 ** 3))
        # DUC_DIAG fix (2026-09-30): cuFFT plans are cached per FFT length and each holds a work area
        # allocated from the (1.5 GB-limited) pool; free_all_blocks() cannot release it.  Day lengths vary
        # (gaps, 40/50/100 Hz), so the plan cache filled the pool and 53,329 of 77,820 station-days of the
        # 2008-2012 scan died with "Out of memory" -> those stations were silently absent from the MF.
        # Now: no plan cache, and any GPU failure falls back to numpy (identical to round-off, 1e-10).
        try: cp.fft.config.get_plan_cache().set_size(0)
        except Exception: pass
        _nr, _ni = np.fft.rfft, np.fft.irfft
        def _gpu(fn, a, n, dt):
            try:
                r = cp.asnumpy(fn(cp.asarray(a, dtype=dt), n=n))
            except Exception:
                r = None
            cp.get_default_memory_pool().free_all_blocks()
            return r
        def rfft(a, n=None, axis=-1, norm=None):
            if axis != -1 or norm is not None or np.asarray(a).ndim != 1 or np.asarray(a).size < 100000:
                return _nr(a, n=n, axis=axis, norm=norm)
            r = _gpu(cp.fft.rfft, a, n, cp.float64)
            return r if r is not None else _nr(a, n=n)
        def irfft(a, n=None, axis=-1, norm=None):
            if axis != -1 or norm is not None or np.asarray(a).ndim != 1 or np.asarray(a).size < 100000:
                return _ni(a, n=n, axis=axis, norm=norm)
            r = _gpu(cp.fft.irfft, a, n, cp.complex128)
            return r if r is not None else _ni(a, n=n)
        np.fft.rfft, np.fft.irfft = rfft, irfft
        _GPUFFT = True
    except Exception:
        _GPUFFT = False


def inventory(net, sta):
    k = f"{net}.{sta}"
    if k not in _INV:
        from obspy import read_inventory
        if len(_INV) > 4: _INV.pop(next(iter(_INV)))
        _INV[k] = read_inventory(f"{RESP}/{k}.xml", format="STATIONXML")
    return _INV[k]


def his_resample(tr, fr, ratio):
    """== tr.interpolate(ratio.denominator*fr, 'lanczos', a=10) followed by decimate(ratio.numerator,
    no_filter=True), evaluated only at the samples decimate keeps.  Integer ratio: the lanczos kernel
    at integer offsets is the identity, so the result is data[::numerator] (to ~1e-16).  Rational
    ratio: obspy's own lanczos_interpolation (the function interpolate calls) evaluated directly on the
    20 Hz grid -- the same kernel sums at the same times, 1/numerator of the work."""
    x = np.require(tr.data, dtype=np.float64)
    if ratio.denominator == 1:
        return x[::ratio.numerator].copy()
    from obspy.signal.interpolation import lanczos_interpolation
    import math
    st = tr.stats.starttime.timestamp; en = tr.stats.endtime.timestamp
    dt_up = 1.0 / (ratio.denominator * fr)
    npts_up = int(math.floor((en - st) / dt_up)) + 1
    n_out = (npts_up - 1) // ratio.numerator + 1
    return np.atleast_1d(lanczos_interpolation(x, st, tr.stats.delta, st, dt_up * ratio.numerator, n_out,
                                               a=10, window="lanczos"))


def cache_one(job):
    net, sta, day, loc, chans = job
    d = pd.Timestamp(day); dst = cache_path(net, sta, d)
    want = chans.split(",")
    if os.path.exists(dst):
        try:
            z = np.load(dst)
            if set(want) <= set(str(c) for c in z["chans"]): return "cached"
        except Exception:
            pass
    src = src_path(net, sta, day)
    if src is None: return "nosrc"
    try:
        _patch_evalresp(); _patch_fft()
        from obspy import read, UTCDateTime
        from obspy.core.stream import Stream
        inv = inventory(net, sta)
        D = read(src)
        D = Stream([tr for tr in D if tr.stats.channel in want and tr.stats.location == loc])
        if len(D) == 0: return "nochan"
        t0 = UTCDateTime(day) - 60; t1 = UTCDateTime(day) + 86400 + 60
        D.trim(t0, t1)
        D = Stream([tr for tr in D if tr.stats.npts > 1])
        if len(D) == 0: return "nochan"
        segs = {}
        for tr in D: segs.setdefault(tr.stats.channel, []).append((tr.stats.starttime.timestamp, tr.stats.endtime.timestamp))
        D.detrend(type="linear")
        D.taper(type="hann", max_percentage=None, max_length=5.0)
        D.remove_response(inventory=inv, output="VEL", pre_filt=(0.2, 0.5, 10.0, 15.0), water_level=80.0)
        D.filter("bandpass", freqmin=FILT[0], freqmax=FILT[1], zerophase=True)
        out = Stream(); chains = set()
        for fr in sorted({float(tr.stats.sampling_rate) for tr in D}):
            G = Stream([tr for tr in D if float(tr.stats.sampling_rate) == fr])
            ratio = Fraction(int(fr), int(FS))
            chains.add("frozen" if ratio.denominator == 1 else "rational")
            for tr in G:
                tr.data = his_resample(tr, fr, ratio)
                tr.stats.sampling_rate = FS
            out += G
        D = out
        D.merge(method=1, fill_value=0)
        tA = UTCDateTime(day); tB = tA + 86400
        D.trim(tA, tB, pad=True, fill_value=0)
        cs = [c for c in want if len(D.select(channel=c))]
        N = int(round(86400 * FS))
        X = np.zeros((len(cs), N), np.float32); V = np.zeros((len(cs), N), bool)
        for i, c in enumerate(cs):
            tr = D.select(channel=c)[0]
            n = min(tr.stats.npts, N); X[i, :n] = tr.data[:n]
            for (a0, a1) in segs.get(c, []):
                i0 = max(0, int(round((a0 - tA.timestamp) * FS))); i1 = min(N, int(round((a1 - tA.timestamp) * FS)))
                if i1 > i0: V[i, i0:i1] = True
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        tmp = dst + f".tmp{os.getpid()}.npz"
        np.savez(tmp, x=X, validp=np.packbits(V, axis=1), chans=np.array(cs), t0=np.float64(tA.timestamp),
                 fs=np.float64(FS), chain=np.array("+".join(sorted(chains))), loc=np.array(loc))
        os.replace(tmp, dst)
        return "built"
    except Exception as e:
        return "error:" + str(e)[:100]


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--month", required=True)
    ap.add_argument("--stage", default="all", choices=["fetch", "cache", "all"])
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--fetch-workers", type=int, default=8)
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    os.makedirs(f"{FIX}/data", exist_ok=True)
    pp = f"{FIX}/data/plan_{a.month}.csv"
    if a.stage in ("fetch", "all"):
        P = month_plan(a.month)
        jobs = [(r.net, r.sta, r.day, r.loc, r.chans) for r in P.itertuples()]
        t0 = time.time()
        from concurrent.futures import ThreadPoolExecutor
        with ThreadPoolExecutor(a.fetch_workers) as ex:
            F = list(ex.map(fetch_one, jobs))
        P["fetch"] = [f[0] for f in F]; P["bytes"] = [f[1] for f in F]
        log(f"[data {a.month}] fetch {P.fetch.str.split(':').str[0].value_counts().to_dict()} "
            f"{P.bytes.sum() / 1e9:.2f} GB new, {time.time() - t0:.0f}s")
        P.to_csv(pp, index=False)
    if a.stage in ("cache", "all"):
        P = pd.read_csv(pp, dtype={"loc": str}).fillna({"loc": ""})
        t1 = time.time()
        jobs2 = sorted([(r.net, r.sta, r.day, r.loc, r.chans) for r in P.itertuples()], key=lambda j: (j[0], j[1], j[2]))
        import multiprocessing as mp
        with mp.Pool(a.workers, maxtasksperchild=200) as pool:
            C = pool.map(cache_one, jobs2, chunksize=max(1, len(jobs2) // (a.workers * 4)))
        cm = dict(zip([(j[1], j[2]) for j in jobs2], C))
        P["cache"] = [cm[(r.sta, r.day)] for r in P.itertuples()]
        log(f"[data {a.month}] cache {P.cache.str.split(':').str[0].value_counts().to_dict()} {time.time() - t1:.0f}s")
        P.to_csv(pp, index=False)


if __name__ == "__main__":
    main()
