#!/usr/bin/env python
"""SH_CONF stage 2 (CPU): high-frequency ratio per candidate (WA step-1 'hfr_rel', lfe_extension_wa/scripts/step1_features.py).

hfr = E(10-19 Hz) / E(2-8 Hz) summed over the band's channels, 12-s window starting 1 s before the expected S
(= the MF template window start), native rate, from the RAW day file (traces/<NET.STA>/<yyyy>/<doy>.mseed, band
preference HH > BH > EH > SH as the cache, >= 40 sps).  hfr_rel = log10(hfr) - median log10 of 300 random windows of
the same station-day.  HF station of a family = the highest-pick_q station of its it0 templates that has a raw file
that day (ranking fixed per family; fallback down the list).  Scrambled-control candidates are read at their site's
rolled position (sh_conf_common.s_index).

  python scripts/sh_conf_hf.py --workers 6  -> OUT/hf/<run>/<date>.parquet, OUT/status/hf_0.csv
"""
import os, sys, json, time, argparse, warnings, traceback, zlib
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sh_conf_common as C
BANDPREF = ["HH", "BH", "EH", "SH"]
HF_BAND = (10.0, 19.0); LF_BAND = (2.0, 8.0); WIN_S = 12.0; LEAD_S = 1.0; NR = 300
_G = {}


def _init():
    M, _ = C.load_tpl(C.TPL["v3b"])
    rank = {f: g.sort_values("pick_q", ascending=False).netsta.tolist() for f, g in M.groupby("fam")}
    _G.update(rank=rank, SITE=C.sites(), TPLS={})


def bp(x, fs, f1, f2):
    from scipy.signal import butter, sosfiltfilt
    return sosfiltfilt(butter(4, [f1 / (fs / 2), f2 / (fs / 2)], btype="band", output="sos"), x, axis=-1)


def station_energy(ns, date):
    """-> (cumsum E28, cumsum E1019, fs) over the day from midnight, or status"""
    from obspy import read, UTCDateTime
    net, sta = ns.split(".")
    p = C.day_path(C.TR, net, sta, date)
    if not os.path.exists(p):
        return None, "nofile"
    try:
        st = read(p)
    except Exception as e:
        return None, "error"
    for tr in st:
        fs_ = float(tr.stats.sampling_rate); fr = round(fs_)
        if fr > 0 and fs_ != fr and abs(fs_ - fr) / fr < 1e-4: tr.stats.sampling_rate = float(fr)
    bands = {}
    for tr in st: bands.setdefault(tr.stats.channel[:2], []).append(tr)
    band = next((b for b in BANDPREF if b in bands and float(bands[b][0].stats.sampling_rate) >= 40), None)
    if band is None:
        return None, "lowrate"
    try:
        s1 = st.select(channel=band + "?")
        locs = sorted({tr.stats.location for tr in s1}); s1 = s1.select(location="" if "" in locs else locs[0])
        D0 = UTCDateTime(date)
        E28 = E19 = None; fs0 = None
        for ch in sorted({tr.stats.channel for tr in s1}):
            s2 = s1.select(channel=ch).copy(); s2.merge(method=1, fill_value=0); tr = s2[0]
            fs0 = float(tr.stats.sampling_rate)
            x = tr.data.astype(np.float64); x -= np.median(x)
            i0 = int(round((tr.stats.starttime - D0) * fs0)); N = int(86400 * fs0)
            y = np.zeros(N); a, b = max(i0, 0), min(i0 + len(x), N)
            if b <= a: continue
            y[a:b] = x[a - i0:b - i0]
            e1 = bp(y, fs0, *LF_BAND) ** 2; e2 = bp(y, fs0, *HF_BAND) ** 2
            E28 = e1 if E28 is None else E28 + e1; E19 = e2 if E19 is None else E19 + e2
        if E28 is None:
            return None, "short"
        return (np.concatenate([[0], np.cumsum(E28)]), np.concatenate([[0], np.cumsum(E19)]), fs0), "ok"
    except Exception:
        return None, "error"


def day_job(args):
    run, date, tpl = args
    t00 = time.time()
    status = []
    try:
        meta = json.load(open(f"{C.OUT}/feat/{run}/{date}.json"))
        F = pd.read_parquet(f"{C.OUT}/feat/{run}/{date}.parquet", columns=["series", "k", "fam", "i"])
        if tpl not in _G["TPLS"]: _G["TPLS"][tpl] = C.load_tpl(tpl)[0]
        M = _G["TPLS"][tpl]; SITE = _G["SITE"]
        out = F[["series", "k", "fam"]].copy(); out["hfr"] = np.nan; out["hfr_rel"] = np.nan; out["hf_sta"] = ""
        if not len(F):
            return run, date, out, status, time.time() - t00
        E = {}; tref = pd.Timestamp(date).timestamp(); L = meta["L"]; pre_n = meta["pre_n"]; PRE = pre_n / C.FS
        for f, gi in F.groupby("fam").groups.items():
            g = F.loc[gi]
            for ns in _G["rank"].get(f, []):
                if ns not in meta["stations"]: continue
                r = M[(M.fam == f) & (M.netsta == ns)]
                if not len(r): continue
                if ns not in E:
                    E[ns] = station_energy(ns, date)
                    status.append((run, date, ns, E[ns][1], ""))
                if E[ns][0] is None: continue
                c1, c2, fs0 = E[ns][0]
                info = meta["stations"][ns]; nw = info["N"] - L + 1
                sh = int(round((float(r.tt_s.iloc[0]) + float(r.best_tau.iloc[0]) - PRE - (info["t0"] - tref)) * C.FS))
                shift = int(meta["shifts"][str(SITE[ns])])
                jd = np.full(len(g), -1, np.int64)
                for ser, sft in (("data", 0), ("c1", shift)):
                    mm = (g.series == ser).to_numpy()
                    if mm.any(): jd[mm] = C.s_index(g.i.to_numpy()[mm], sh, sft, nw, pre_n)
                ts = info["t0"] + jd / C.FS - LEAD_S - tref          # window start, s after midnight
                nwn = int(round(WIN_S * fs0)); j = np.round(ts * fs0).astype(np.int64)
                ok = (jd >= 0) & (j >= 0) & (j + nwn < len(c1))
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
        return run, date, out, status, time.time() - t00
    except Exception as e:
        msg = f"day job: {str(e)[:100]} {traceback.format_exc()[-300:]}"
        return run, date, None, status + [(run, date, ns, "error", msg) for ns in C.STAS], time.time() - t00


def main():
    from concurrent.futures import ProcessPoolExecutor, as_completed
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=6); ap.add_argument("--days", default="")
    ap.add_argument("--runs", default="v3b,fwd")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    D = pd.read_csv(f"{C.OUT}/day_table.csv", dtype={"date": str})
    D = D[D.run.isin(a.runs.split(","))]
    if a.days: D = D[D.date.isin(a.days.split(","))]
    todo = [(r.run, r.date, r.tpl) for r in D.itertuples() if os.path.exists(f"{C.OUT}/feat/{r.run}/{r.date}.json")
            and not os.path.exists(f"{C.OUT}/hf/{r.run}/{r.date}.parquet")]
    C.log(f"[hf] {len(todo)} run-days, {a.workers} workers")
    sp = f"{C.OUT}/status/hf_0.csv"
    status = [tuple(x) for x in pd.read_csv(sp).fillna("").itertuples(index=False)] if os.path.exists(sp) else []
    t0 = time.time()
    with ProcessPoolExecutor(a.workers, initializer=_init) as ex:
        fut = [ex.submit(day_job, j) for j in todo]
        for n, fu in enumerate(as_completed(fut)):
            run, date, out, st, secs = fu.result()
            status += st
            if out is not None:
                od = f"{C.OUT}/hf/{run}"; os.makedirs(od, exist_ok=True)
                out.to_parquet(f"{od}/{date}.parquet.tmp", index=False); os.replace(f"{od}/{date}.parquet.tmp", f"{od}/{date}.parquet")
            else:
                C.log(f"  ERROR {run} {date}: {st[-1][-1]}")
            if n % 25 == 0 or n < 3:
                C.log(f"  {run} {date}: {0 if out is None else out.hfr_rel.notna().mean():.2f} defined, {secs:.0f}s; {n + 1}/{len(todo)} {time.time() - t0:.0f}s")
                C.write_status("hf", "0", status)
    C.write_status("hf", "0", status)
    C.log(f"[hf] done {len(todo)} in {time.time() - t0:.0f}s")
    C.guard("hf")


if __name__ == "__main__":
    main()
