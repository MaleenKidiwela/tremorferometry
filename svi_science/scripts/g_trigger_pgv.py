#!/usr/bin/env python
"""g_trigger_pgv (2026-10-11) -- Study G step 2: measured peak ground velocity / dynamic stress at CN.PGC from raw
files on disk (PLAN.md §2).  8 workers max.

    python svi_science/scripts/g_trigger_pgv.py   -> results/g_trigger/tables/pgv_pgc.csv
"""
import os, glob, warnings
from concurrent.futures import ProcessPoolExecutor
import numpy as np
import pandas as pd
from obspy import read, read_inventory, UTCDateTime, Stream

warnings.filterwarnings("ignore")
ROOT = "/home/jovyan/tremorferometry"
G = f"{ROOT}/svi_science/results/g_trigger"
DIRS = [f"{ROOT}/lfe_extension/traces/CN.PGC", f"{ROOT}/lfe_extension/traces_backward/CN.PGC"]
INV = None


def files_for(t):
    out = []
    for d in DIRS:
        p = f"{d}/{t.year}/{t.julday:03d}.mseed"
        for q in (p, p + ".nd"):
            if os.path.exists(q):
                out.append(q); break
        if out:
            break
    return out


def load(t1, t2):
    st = Stream()
    day = UTCDateTime(t1.year, t1.month, t1.day)
    while day < t2:
        for f in files_for(day):
            try:
                st += read(f, starttime=t1, endtime=t2)
            except Exception:
                pass
        day += 86400
    # R7: backward day files overlap the next day by ~1 min and the next file starts with a digitiser transient;
    # keep the earlier piece where they overlap (obspy merge method 1 would keep the later one).
    out = Stream()
    for cid in sorted(set(tr.id for tr in st)):
        prev_end = None
        for tr in sorted(st.select(id=cid), key=lambda x: x.stats.starttime):
            if prev_end is not None and tr.stats.starttime <= prev_end:
                if tr.stats.endtime <= prev_end:
                    continue
                tr = tr.slice(prev_end + tr.stats.delta, tr.stats.endtime)
            out += tr
            prev_end = tr.stats.endtime
    return out


def one(row):
    global INV
    if INV is None:
        INV = read_inventory(f"{ROOT}/codex-LFE/results/full_history_metadata_v1/CN.xml").select(station="PGC")
    t0 = UTCDateTime(pd.Timestamp(row["t0_utc"]).to_pydatetime())
    tP = t0 + row["SVI_tP_s"]; tE = t0 + row["SVI_tR_s"] + 1800
    pad = 1800
    res = dict(id=row["id"], pgc_ok=False, pgc_reason="")
    try:
        st = load(tP - pad, tE + pad)
        if len(st) == 0:
            res["pgc_reason"] = "no_file"; return res
        for tr in st:
            tr.stats.sampling_rate = round(tr.stats.sampling_rate, 3)   # 100.0000076 Hz jitter blocks merge
            tr.data = tr.data.astype(np.float64) - np.median(tr.data)   # R7: per-piece offset (day-file splice step)
        st.merge(method=1, fill_value=None)
        comps = {}
        for tr in st:
            if tr.stats.channel[:2] not in ("BH", "HH"):
                continue
            comps.setdefault(tr.stats.channel[-1], tr)
        if len(comps) < 3:
            res["pgc_reason"] = f"comps_{len(comps)}"; return res
        pv, pv_hf, rms_hf, cov = [], [], [], []
        for c, tr in comps.items():
            seg = tr.slice(tP, tE)
            n_need = (tE - tP) * tr.stats.sampling_rate
            dd = seg.data
            good = (~np.ma.getmaskarray(dd)).sum() if np.ma.isMaskedArray(dd) else len(dd)
            cov.append(min(1.0, good / n_need))
            if np.ma.isMaskedArray(tr.data):
                tr.data = tr.data.filled(0)
            # R7: 1-s BHE glitch at 00:00:00-00:00:01 UTC at the end of day files -> blank midnight +-3 s
            m0 = UTCDateTime(tr.stats.starttime.year, tr.stats.starttime.month, tr.stats.starttime.day) + 86400
            while m0 < tr.stats.endtime:
                i0 = int((m0 - 3 - tr.stats.starttime) * tr.stats.sampling_rate); i1 = int((m0 + 3 - tr.stats.starttime) * tr.stats.sampling_rate)
                tr.data[max(i0, 0):max(i1, 0)] = 0
                m0 += 86400
            seg = tr.slice(tP, tE)
            a_ = np.abs(np.ma.filled(seg.data, 0))       # R6: flat-top test (counts exceed 2^23 unclipped)
            if (a_ == a_.max()).sum() >= 3:
                res["pgc_reason"] = "clip"; return res
            trv = tr.copy()
            trv.detrend("demean"); trv.taper(0.05)
            if trv.stats.sampling_rate > 20:
                trv.filter("lowpass", freq=8.0, corners=4, zerophase=True)
                trv.decimate(int(trv.stats.sampling_rate // 20), no_filter=True)
            # CN.xml carries the overall sensitivity only (no stages): flat-velocity broadband in 0.02-2 Hz
            sens = INV.get_response(tr.id, tP).instrument_sensitivity.value
            trv.data = trv.data / sens
            a = trv.copy().filter("bandpass", freqmin=0.02, freqmax=2.0, corners=4, zerophase=True).slice(tP, tE)
            b = trv.copy().filter("bandpass", freqmin=2.0, freqmax=8.0, corners=4, zerophase=True).slice(tP, tE)
            pv.append(np.abs(a.data).max()); pv_hf.append(np.abs(b.data).max()); rms_hf.append(np.sqrt(np.mean(b.data ** 2)))
        res["pgc_cov"] = float(min(cov))
        if min(cov) < 0.95:
            res["pgc_reason"] = "gap"; res["pgv_pgc"] = float(max(pv)); return res
        res.update(pgc_ok=True, pgv_pgc=float(max(pv)), pgv_hf_pgc=float(max(pv_hf)), rms_hf_pgc=float(np.mean(rms_hf)))
        return res
    except Exception as e:
        res["pgc_reason"] = f"err:{type(e).__name__}:{str(e)[:60]}"
        return res


def main():
    E = pd.read_csv(f"{G}/tables/eq_list_raw.csv", parse_dates=["t0_utc"])
    E = E[E.t0_utc >= "2003-01-01"]
    rows = E.to_dict("records")
    with ProcessPoolExecutor(8) as ex:
        R = list(ex.map(one, rows, chunksize=4))
    P = pd.DataFrame(R)
    P["eps_pgc"] = P.pgv_pgc / 3500.0
    P["sig_pgc_kpa"] = 30e9 * P.eps_pgc / 1e3
    P.loc[~P.pgc_ok.astype(bool), "sig_pgc_kpa"] = np.nan
    P.to_csv(f"{G}/tables/pgv_pgc.csv", index=False)
    print(P.pgc_ok.value_counts(), P.pgc_reason.value_counts().head(10))


if __name__ == "__main__":
    main()
