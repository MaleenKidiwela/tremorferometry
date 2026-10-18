#!/usr/bin/env python
"""DUC_FIN high-frequency ratio per candidate (SH_CONF 'hfr_rel', WA step-1 design) from the RAW miniSEED.

hfr = E(10-19 Hz) / E(2-8 Hz) summed over the epoch's channels, 12-s window starting 1 s before the expected S at the
family's HF station (S = detection time + duc_fin_common.arrival_offsets 'template' offset), native rate, raw counts
(median removed; the ratio is unit-free).  hfr_rel = log10(hfr) - median log10 hfr of 300 random windows of the same
station-day.  HF station of a family = its template stations ranked by HIS template's energy-envelope peak / median
(fixed per family), first one with a raw file that day and >= 40 sps.  Scrambled-control candidates are read at their
site's rolled (original-time) position.

  python scripts/duc_fin_hf.py --workers 10 [--days ...]   -> OUT/hf/<date>.parquet, OUT/status/hf_0.csv
"""
import os, sys, json, time, argparse, warnings, traceback, zlib
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duc_fin_common as F
import duc_fix_common as DC
HF_BAND = (10.0, 19.0); LF_BAND = (2.0, 8.0); WIN_S = 12.0; LEAD_S = 1.0; NR = 300
NDAY = int(round(86400 * DC.FS))
_G = {}


def hf_rank():
    T = DC.his_templates(); no = F.netof(); rank = {}
    for f, d in T.items():
        q = {}
        for (sta, cha), tp in d.items():
            e = np.convolve(tp * tp, np.ones(20) / 20, "same"); q.setdefault(sta, []).append(e.max() / max(np.median(e), 1e-30))
        rank[f] = [f"{no[s]}.{s}" for s, v in sorted(q.items(), key=lambda kv: -np.median(kv[1]))]
    return rank


def _init():
    O = pd.read_csv(f"{F.OUT}/arrival_offsets.csv", dtype={"fam": str})
    O = O[O.offset_source == "template"]
    _G.update(rank=hf_rank(), off={(r.fam, r.netsta): r.offset_s for r in O.itertuples()}, SITE=F.sites(),
              look=DC.epoch_lookup(), plan=F.plan_status())


def bp(x, fs, f1, f2):
    from scipy.signal import butter, sosfiltfilt
    return sosfiltfilt(butter(4, [f1 / (fs / 2), min(f2, 0.95 * fs / 2) / (fs / 2)], btype="band", output="sos"), x, axis=-1)


def station_energy(ns, date):
    from obspy import read, UTCDateTime
    import duc_fix_data as DD
    net, sta = ns.split(".")
    p = DD.src_path(net, sta, date)
    if p is None:
        st_ = _G["plan"].get((ns, date))
        return None, ("nodata" if st_ in (None, "nosrc") else "missing")
    ep = _G["look"](sta, pd.Timestamp(date))
    if ep is None: return None, "nodata"
    try:
        st = read(p)
    except Exception:
        return None, "error"
    try:
        loc = ep.loc if isinstance(ep.loc, str) else ""
        st = st.select(location=loc)
        st = type(st)([tr for tr in st if tr.stats.channel in ep.chans.split(",")])
        if not len(st): return None, "nochan"
        fs0 = float(st[0].stats.sampling_rate)
        if fs0 < 40: return None, "lowrate"
        D0 = UTCDateTime(date); E28 = E19 = None
        for ch in sorted({tr.stats.channel for tr in st}):
            s2 = st.select(channel=ch).copy(); s2.merge(method=1, fill_value=0); tr = s2[0]
            fs_ = float(tr.stats.sampling_rate)
            if abs(fs_ - fs0) > 1e-3: continue
            x = tr.data.astype(np.float64); x -= np.median(x)
            i0 = int(round((tr.stats.starttime - D0) * fs0)); N = int(86400 * fs0)
            y = np.zeros(N); a, b = max(i0, 0), min(i0 + len(x), N)
            if b <= a: continue
            y[a:b] = x[a - i0:b - i0]
            e1 = bp(y, fs0, *LF_BAND) ** 2; e2 = bp(y, fs0, *HF_BAND) ** 2
            E28 = e1 if E28 is None else E28 + e1; E19 = e2 if E19 is None else E19 + e2
        if E28 is None: return None, "short"
        return (np.concatenate([[0], np.cumsum(E28)]), np.concatenate([[0], np.cumsum(E19)]), fs0), "ok"
    except Exception:
        return None, "error"


def day_job(date):
    t00 = time.time(); status = []
    try:
        meta = json.load(open(f"{F.OUT}/feat/{date}.json"))
        pq = f"{F.OUT}/feat/{date}.parquet"
        Fr = pd.read_parquet(pq, columns=["series", "k", "fam", "g", "t_epoch"]) if os.path.getsize(pq) > 0 else pd.DataFrame()
        out = Fr[["series", "k"]].copy() if len(Fr) else pd.DataFrame(columns=["series", "k"])
        out["hfr"] = np.nan; out["hfr_rel"] = np.nan; out["hf_sta"] = ""
        if not len(Fr): return date, out, status, time.time() - t00
        E = {}; tref = pd.Timestamp(date).timestamp(); avail = set(meta["stations"])
        for f, gi in Fr.groupby("fam").groups.items():
            g = Fr.loc[gi]
            for ns in _G["rank"].get(str(f), []):
                if ns not in avail or (str(f), ns) not in _G["off"]: continue
                if ns not in E:
                    E[ns] = station_energy(ns, date); status.append((date, ns, E[ns][1], ""))
                if E[ns][0] is None: continue
                c1, c2, fs0 = E[ns][0]
                off = _G["off"][(str(f), ns)]
                tdet = g.t_epoch.to_numpy().copy()
                mc = (g.series == "c1").to_numpy()
                if mc.any():
                    sh = int(meta["shifts"][str(_G["SITE"][ns])])
                    tdet[mc] = tref + ((g.g.to_numpy()[mc] % NDAY + sh) % NDAY) / DC.FS
                ts = tdet + off - LEAD_S - tref
                nwn = int(round(WIN_S * fs0)); j = np.round(ts * fs0).astype(np.int64)
                ok = (j >= 0) & (j + nwn < len(c1))
                jj = np.clip(j, 0, len(c1) - nwn - 1)
                s28 = c1[jj + nwn] - c1[jj]; s19 = c2[jj + nwn] - c2[jj]
                with np.errstate(divide="ignore", invalid="ignore"):
                    hr = np.where(ok & (s28 > 0), s19 / s28, np.nan)
                rng = np.random.default_rng(int(date.replace("-", "")) + zlib.crc32(ns.encode()) % 1000)
                jr = rng.integers(0, len(c1) - nwn - 1, 4 * NR)
                r28 = c1[jr + nwn] - c1[jr]; r19 = c2[jr + nwn] - c2[jr]; kk = (r28 > 0) & (r19 > 0)
                ref = float(np.median(np.log10(r19[kk][:NR] / r28[kk][:NR]))) if kk.sum() >= 50 else np.nan
                out.loc[gi, "hfr"] = hr
                with np.errstate(divide="ignore", invalid="ignore"):
                    out.loc[gi, "hfr_rel"] = np.log10(hr) - ref
                out.loc[gi, "hf_sta"] = ns
                break
        return date, out, status, time.time() - t00
    except Exception as e:
        try: names = json.load(open(f"{F.OUT}/feat/{date}.json"))["stations"]
        except Exception: names = ["ALL"]
        return date, None, status + [(date, ns, "error", f"day job: {str(e)[:100]} {traceback.format_exc()[-200:]}") for ns in names], time.time() - t00


def main():
    from concurrent.futures import ProcessPoolExecutor, as_completed
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=10); ap.add_argument("--days", default="")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    days = a.days.split(",") if a.days else list(F.day_table().date)
    todo = [d for d in days if os.path.exists(f"{F.OUT}/feat/{d}.json") and not os.path.exists(f"{F.OUT}/hf/{d}.parquet")]
    F.log(f"[hf] {len(todo)} days, {a.workers} workers")
    sp = f"{F.OUT}/status/hf_0.csv"
    status = [tuple(x) for x in pd.read_csv(sp).fillna("").itertuples(index=False)] if os.path.exists(sp) else []
    t0 = time.time(); od = f"{F.OUT}/hf"; os.makedirs(od, exist_ok=True)
    with ProcessPoolExecutor(a.workers, initializer=_init) as ex:
        fut = [ex.submit(day_job, d) for d in todo]
        for n, fu in enumerate(as_completed(fut)):
            date, out, st, secs = fu.result(); status += st
            if out is not None:
                out.to_parquet(f"{od}/{date}.parquet.tmp", index=False); os.replace(f"{od}/{date}.parquet.tmp", f"{od}/{date}.parquet")
            else:
                F.log(f"  ERROR {date}: {st[-1][-1]}")
            if n % 50 == 0 or n < 3:
                F.log(f"  {date}: {0 if out is None or not len(out) else out.hfr_rel.notna().mean():.2f} defined, {secs:.0f}s; {n + 1}/{len(todo)} {time.time() - t0:.0f}s")
                F.write_status("hf", "0", status)
    F.write_status("hf", "0", status)
    F.log(f"[hf] done {len(todo)} in {time.time() - t0:.0f}s")
    F.guard("hf")


if __name__ == "__main__":
    main()
