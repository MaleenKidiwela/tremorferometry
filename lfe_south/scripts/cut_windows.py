#!/usr/bin/env python
"""CUT REPEAT WINDOWS (2026-09-30) -- keep raw-count windows around EVERY family detection so the continuous
raw miniSEED can later be deleted (user decision 09-30: deletion only after the user confirms an exact list).

For every candidate detection (>= 8 x MAD, all of them, with their statistics and catalogue flags) and EVERY
station with a raw day file at that time (the run's stations + any other station in the raw folder(s) that
day), cut -15 s .. +30 s around the family's expected arrival at that station:
    arrival = detection time + that station's template offset        (offset_source = template:<file>)
    arrival = detection time (offset 0)  for a station without one    (offset_source = detection_time)
The window keeps the S arrival at +15 s, a 15 s pre-arrival (mirror / causality) slice and 30 s after it
(coda 2-4 s after S and beyond).  Raw counts, native sample rate, no filtering, every channel of the
station's preferred band (HH > BH > EH > SH, >= 20 sps -- the t3_cache band rule), all location codes.

Storage (deduplicated): windows of one station overlap heavily (27 families detect the same event, 4 s
declustering), so each channel stores the UNION of its windows as contiguous segments:
  results/cut/<run>/windows_<YYYY-MM>.h5   one resizable 1-D dataset per channel  /data/<NET.STA.LOC.CHA>
                                           (int32 raw counts, or float32 when the miniSEED is float), gzip+shuffle
  results/cut/<run>/segments_<YYYY-MM>.parquet  seg_id, seed_id, dataset, h5_start, n, t0 (first sample, epoch s),
                                           fs, gaps (JSON [[i0, i1], ...] zero-filled sample runs), src path, resp_epoch
  results/cut/<run>/index_<YYYY-MM>.parquet    one row per detection x station x channel: det_id, fam, t_det,
                                           netsta, loc, cha, seed_id, fs, t0 (window first sample), offset_s,
                                           offset_source, n_samples, n_expected, n_gap, partial, h5_file, dataset,
                                           h5_start, seg_id, resp_epoch, catalogue flags, detection statistics
  results/cut/<run>/detections.parquet     the full candidate table with det_id and flags
  results/cut/<run>/stationdays.parquet    per raw station-day file: path, bytes, status, windows cut
  results/cut/<run>/channel_epochs.csv     FDSN channel epochs (resp_epoch id = SEED@epoch-start)
  results/cut/<run>/cut_summary.csv        per month: detections, windows, station-days, raw GB, cut GB
read_window(index_row) returns (data, t0, fs) for one window.

  python scripts/cut_windows.py cut        --run sh_ext [--workers 6] [--months 2024-01,2024-02]
  python scripts/cut_windows.py verify     --run sh_ext [--n 20]
  python scripts/cut_windows.py stackcheck --run sh_ext
  python scripts/cut_windows.py deletelist --run sh_ext
Runs: sh_ext (Shelly forward, 2023-12-01..), sh_back (Shelly backward BK2, 2008-2017), duc_fix (Ducellier 2008-2026).
Nothing is deleted by this script.
"""
import os, sys, json, glob, time, argparse, warnings, collections
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from t_common import WS, RES, TR, CACHE, log, day_path

TR2 = f"{WS}/traces_duc"
PRE_S, POST_S = 15.0, 30.0
BANDPREF = ["HH", "BH", "EH", "SH"]
PROV = {'NC': 'NCEDC', 'BK': 'NCEDC', 'PB': 'EARTHSCOPE', 'UW': 'EARTHSCOPE', 'CC': 'EARTHSCOPE', 'US': 'EARTHSCOPE', 'UO': 'EARTHSCOPE'}


# ================================================================================================ run profiles
def _flag_match(D, C, tol=0.05, same_fam=True):
    """True where (fam, t) of D has an entry in C within tol s (same family if same_fam)."""
    out = np.zeros(len(D), bool)
    if C is None or not len(C): return out
    te = D.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    ce = C.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    if same_fam:
        for f, idx in D.groupby("fam").indices.items():
            c = np.sort(ce[C.fam.values == f])
            if not c.size: continue
            x = te[idx]; j = np.clip(np.searchsorted(c, x), 1, c.size - 1) if c.size > 1 else np.zeros(x.size, int)
            d = np.minimum(np.abs(c[j] - x), np.abs(c[np.maximum(j - 1, 0)] - x))
            out[idx] = d <= tol
    else:
        c = np.sort(ce); j = np.clip(np.searchsorted(c, te), 1, max(c.size - 1, 1))
        d = np.minimum(np.abs(c[np.minimum(j, c.size - 1)] - te), np.abs(c[np.maximum(j - 1, 0)] - te))
        out = d <= tol
    return out


def prof_sh_ext():
    E = os.environ.get("CUT_SH_EXT_DIR", f"{RES}/sh_ext2026")       # override only for smoke tests
    D = pd.read_parquet(f"{E}/candidates_raw8x.parquet")
    D["fam"] = D.fam.astype(str); D["t"] = D.t.astype("datetime64[ns]")
    D = D.rename(columns={"net_mad": "xmad"})
    CAT = pd.read_parquet(f"{E}/cat_S1_9x_internal.parquet"); CAT["fam"] = CAT.fam.astype(str); CAT["t"] = CAT.t.astype("datetime64[ns]")
    D["in_cat_S1_9x"] = _flag_match(D, CAT, 0.01)
    from sh_v3_common import his_catalogue
    S = his_catalogue(); S = S[S.t >= pd.Timestamp("2023-11-30")]
    D["his_same_fam_2s"] = _flag_match(D, S, 2.0); D["his_any_fam_2s"] = _flag_match(D, S, 2.0, same_fam=False)
    D.loc[D.t >= pd.Timestamp("2024-12-01"), ["his_same_fam_2s", "his_any_fam_2s"]] = False   # his catalogue ends 2024-11-30
    M = pd.read_csv(f"{RES}/sh_v3b/tpl_it0.csv"); M["fam"] = M.fam.astype(str)
    OFF = pd.DataFrame(dict(fam=M.fam, netsta=M.netsta, offset_s=M.tt_s + M.best_tau, offset_source="template:sh_v3b/tpl_it0 (tt_s+best_tau = S arrival)"))
    import sh_ext2026 as X
    return dict(det=D, off=OFF, stations=list(X.STAS), t0="2023-12-01", t1=str((X.T_END).date()),
                roots=[TR, TR2], flags=["in_cat_S1_9x", "his_same_fam_2s", "his_any_fam_2s"],
                stats=["xmad", "ccsum", "meancc", "nlive", "nchan", "nsta", "seg"], tpl=f"{RES}/sh_v3b/tpl_it0",
                done="SH_EXT_DONE", failed="SH_EXT_FAILED", label="shelly_forward")


def prof_sh_back():
    B = f"{RES}/sh_back"
    D = pd.read_parquet(f"{B}/candidates_bk2_raw8x.parquet"); D["fam"] = D.fam.astype(str); D["t"] = D.t.astype("datetime64[ns]")
    D = D.rename(columns={"net_mad": "xmad"})
    CAT = pd.read_parquet(f"{B}/catalogue_bk2_S1_9x.parquet"); CAT["fam"] = CAT.fam.astype(str); CAT["t"] = CAT.t.astype("datetime64[ns]")
    D["in_cat_S1_9x"] = _flag_match(D, CAT, 0.01)
    M = pd.read_csv(f"{B}/tpl_S6.csv"); M["fam"] = M.fam.astype(str)
    OFF = pd.DataFrame(dict(fam=M.fam, netsta=M.netsta, offset_s=M.tt_s + M.best_tau, offset_source="template:sh_back/tpl_S6 (tt_s+best_tau)"))
    return dict(det=D, off=OFF, stations=sorted(M.netsta.unique()), t0="2008-01-01", t1="2018-01-01", roots=[TR, TR2],
                flags=["in_cat_S1_9x"], stats=["xmad", "ccsum", "meancc", "nlive", "nchan", "nsta"], tpl=f"{B}/tpl_S6",
                done="SH_BACK_DONE", failed="SH_BACK_FAILED", label="shelly_backward")


def prof_duc_fix():
    F = f"{RES}/duc_fix"
    ps = {os.path.basename(p): p for p in glob.glob(f"{F}/mf1/*.parquet")}
    ps.update({os.path.basename(p): p for p in glob.glob(f"{F}/mf2/*.parquet")})      # pass-2 re-scans replace pass 1
    D = pd.concat([pd.read_parquet(p).assign(pass_=("mf2" if "/mf2/" in p else "mf1")) for p in sorted(ps.values())], ignore_index=True)
    D["fam"] = D.fam.astype(str); D["t"] = D.t.astype("datetime64[ns]"); D["xmad"] = D.cc / D.mad
    D["pass_his"] = D.pass_his.astype(bool)
    import duc_fix_common as DC
    T = DC.his_templates(); rows = []
    for f, d in T.items():
        for (sta, cha), tp in d.items():
            env = np.convolve(np.abs(tp), np.ones(5) / 5, "same"); rows.append((f, sta, cha, float(np.argmax(env)) / DC.FS))
    O = pd.DataFrame(rows, columns=["fam", "sta", "cha", "pk"])
    H = DC.his_rows(); netof = dict(zip(H.sta, H.net))
    O = O.groupby(["fam", "sta"]).pk.median().reset_index()
    OFF = pd.DataFrame(dict(fam=O.fam, netsta=[f"{netof.get(s, 'NC')}.{s}" for s in O.sta], offset_s=O.pk,
                            offset_source="template:ducellier templates_v1 (envelope peak after template start)"))
    stations = sorted({f"{r.net}.{r.sta}" for r in H.itertuples()})
    return dict(det=D, off=OFF, stations=stations, t0="2008-01-01", t1="2026-10-01", roots=[TR, TR2],
                flags=["pass_his"], stats=["xmad", "cc", "nchannel", "nlive", "mad", "ccsum", "thr_perm", "pass_"], tpl=None,
                done="DUC_FIX_DONE", failed="DUC_FIX_FAILED", label="ducellier")


PROFILES = dict(sh_ext=prof_sh_ext, sh_back=prof_sh_back, duc_fix=prof_duc_fix)


# ================================================================================================ channel epochs
def channel_epochs(netstas, t0, t1, path):
    have = pd.read_csv(path, dtype={"loc": str}) if os.path.exists(path) else pd.DataFrame(columns=["netsta"])
    todo = [ns for ns in netstas if ns not in set(have.netsta)]
    rows = []
    if todo:
        from obspy.clients.fdsn import Client
        from obspy import UTCDateTime
        for ns in todo:
            n, s = ns.split("."); got = False
            for prov in [PROV.get(n, "EARTHSCOPE")] + [p for p in ("NCEDC", "EARTHSCOPE") if p != PROV.get(n, "EARTHSCOPE")]:
                try:
                    inv = Client(prov, timeout=120).get_stations(network=n, station=s, channel="HH?,BH?,EH?,SH?", level="channel",
                                                                 starttime=UTCDateTime(t0) - 86400, endtime=UTCDateTime(t1) + 86400)
                except Exception as e:
                    continue
                for net in inv:
                    for st in net:
                        for ch in st:
                            sens = ch.response.instrument_sensitivity if ch.response is not None else None
                            rows.append(dict(netsta=ns, loc=ch.location_code, cha=ch.code, start=str(ch.start_date)[:19],
                                             end=str(ch.end_date)[:19] if ch.end_date else "", sr=ch.sample_rate,
                                             sensor=(ch.sensor.description or ch.sensor.model or "") if ch.sensor else "",
                                             gain=sens.value if sens else np.nan, az=ch.azimuth, dip=ch.dip, prov=prov))
                got = True; break
            if not got: rows.append(dict(netsta=ns, loc="", cha="", start="", end="", sr=np.nan, sensor="FDSN query failed", prov=""))
        have = pd.concat([have, pd.DataFrame(rows)], ignore_index=True)
        have.to_csv(path, index=False)
    have["loc"] = have["loc"].fillna("").astype(str)
    return have


def epoch_lookup(E):
    """{seed_id: [(t0_epoch, t1_epoch, epoch_id)]}"""
    L = {}
    for r in E.itertuples():
        if not isinstance(r.cha, str) or not r.cha: continue
        seed = f"{r.netsta.split('.')[0]}.{r.netsta.split('.')[1]}.{r.loc}.{r.cha}"
        a = pd.Timestamp(r.start).timestamp()
        en = str(r.end) if isinstance(r.end, str) else ""
        b = pd.Timestamp(en).timestamp() if en and not en.startswith("3000") and en != "nan" else 4102444800.0
        L.setdefault(seed, []).append((a, b, f"{seed}@{str(r.start)[:19]}"))
    return L


# ================================================================================================ cutting
_G = {}


def _init(G):
    warnings.filterwarnings("ignore"); _G.update(G)


def union(iv, join=1.0):
    iv = sorted(iv); out = []
    for a, b in iv:
        if out and a <= out[-1][1] + join: out[-1][1] = max(out[-1][1], b)
        else: out.append([a, b])
    return out


def build_segment(trs, a, b):
    """samples of one seed id covering [a, b] from the (unmerged) traces; returns (data, t0, fs, gaps) or None.
    The grid is the first overlapping trace's own sample grid; other traces are placed by rounding (gaps -> 0)."""
    ov = [t for t in trs if t.stats.endtime.timestamp >= a and t.stats.starttime.timestamp <= b]
    if not ov: return None
    T = ov[0]; fs = float(T.stats.sampling_rate); ts = T.stats.starttime.timestamp
    k0 = int(np.ceil((a - ts) * fs - 1e-6)); k0 = max(k0, 0)
    t0 = ts + k0 / fs
    n = int(np.floor((b - t0) * fs + 1e-6)) + 1
    if n <= 0: return None
    dt = np.float32 if T.data.dtype.kind == "f" else np.int32
    out = np.zeros(n, dt); filled = np.zeros(n, bool)
    for t in ov:
        f2 = float(t.stats.sampling_rate)
        if abs(f2 - fs) / fs > 1e-4: continue
        s2 = t.stats.starttime.timestamp
        j0 = max(int(np.ceil((t0 - s2) * f2 - 0.01)), 0)          # t0 lies on this grid for the anchor trace
        p0 = int(round((s2 + j0 / f2 - t0) * fs))
        m = min(t.stats.npts - j0, n - p0)
        if m <= 0 or p0 >= n: continue
        if p0 < 0: j0 -= p0; m += p0; p0 = 0
        if m <= 0: continue
        x = np.asarray(t.data[j0:j0 + m])
        if np.ma.isMaskedArray(x): x = x.filled(0)
        out[p0:p0 + m] = x.astype(dt); filled[p0:p0 + m] = True
    if not filled.any(): return None
    d = np.diff(np.concatenate([[1], filled.astype(np.int8), [1]]))
    gaps = list(zip(np.flatnonzero(d == -1).tolist(), np.flatnonzero(d == 1).tolist()))
    return out, t0, fs, gaps


def cut_day(job):
    """one day: every station file present in the raw roots x every detection of the day."""
    day, dets = job                      # dets: array (det_id, fam_code, t_epoch)
    from obspy import read
    G = _G; d = pd.Timestamp(day); segs = []; idx = []; sdays = []
    for ns in G["allsta"]:
        net, sta = ns.split(".")
        src = None
        for r in G["roots"]:
            p = day_path(r, net, sta, d)
            if os.path.exists(p): src = p; break
        if src is None: continue
        sd = dict(day=day, netsta=ns, path=src, bytes=os.path.getsize(src), status="ok", n_win=0, n_partial=0, n_chan=0, band="")
        try:
            st = read(src)
        except Exception as e:
            sd.update(status="read_error:" + str(e)[:80]); sdays.append(sd); continue
        bands = {}
        for tr in st: bands.setdefault(tr.stats.channel[:2], []).append(tr)
        band = next((b_ for b_ in BANDPREF if b_ in bands and max(float(t.stats.sampling_rate) for t in bands[b_]) >= 20), None)
        if band is None:
            sd.update(status="no_band"); sdays.append(sd); continue
        sd["band"] = band
        byseed = {}
        for tr in bands[band]: byseed.setdefault(tr.id, []).append(tr)
        sd["n_chan"] = len(byseed)
        # windows of this station
        offs = G["off"].get(ns, {})
        W = []
        for det_id, fc, te in dets:
            f = G["fams"][fc]
            if f in offs: o, osrc = offs[f], 1
            else: o, osrc = 0.0, 0
            arr = te + o
            W.append((det_id, fc, te, o, osrc, arr - PRE_S, arr + POST_S))
        if not W:
            sdays.append(sd); continue
        U = union([(w[5], w[6]) for w in W])
        for seed, trs in sorted(byseed.items()):
            trs = sorted(trs, key=lambda t: t.stats.starttime)
            ss = []
            for a, b in U:
                r = build_segment(trs, a, b)
                if r is None: continue
                ss.append(r)
            if not ss: continue
            base = len(segs)
            for (x, t0, fs, gaps) in ss:
                segs.append(dict(seed_id=seed, data=x, t0=t0, fs=fs, gaps=gaps, src=src))
            st0 = np.array([s_[1] for s_ in ss]); fs_ = ss[0][2]
            ends = np.array([s_[1] + (s_[0].size - 1) / s_[2] for s_ in ss])
            loc, cha = seed.split(".")[2], seed.split(".")[3]
            nexp = int(round((PRE_S + POST_S) * fs_))                 # samples in [arrival-15 s, arrival+30 s)
            for det_id, fc, te, o, osrc, wa, wb in W:
                k = int(np.searchsorted(st0, wa + 0.5 / fs_, "right") - 1)
                if k < 0 or ends[k] < wa - 0.5 / fs_:
                    k2 = int(np.searchsorted(st0, wa))          # window starts in a data gap: take the next segment if it overlaps
                    if k2 >= len(ss) or st0[k2] > wb: continue
                    k = k2
                x, t0, fs, gaps = ss[k]
                i0 = int(np.ceil((wa - t0) * fs - 1e-6)); i1 = i0 + nexp
                j0, j1 = max(i0, 0), min(i1, x.size)
                if j1 <= j0: continue
                ng = sum(max(0, min(g1, j1) - max(g0, j0)) for g0, g1 in gaps)
                part = (j0 != i0) or (j1 != i1)
                idx.append((det_id, ns, loc, cha, seed, fs, t0 + j0 / fs, o, osrc, j1 - j0, nexp, ng, part, base + k, j0))
                sd["n_win"] += 1; sd["n_partial"] += int(part)
        sdays.append(sd)
    if G.get("out"):                     # day-file mode: this worker compresses + writes windows_<day>.h5 itself
        import h5py
        fn = f"{G['out']}/windows_{day}.h5"
        with h5py.File(fn + ".tmp", "w") as h:
            h.attrs.update(day=day, pre_s=PRE_S, post_s=POST_S, created=str(pd.Timestamp.now("UTC")))
            for s_ in segs:
                x = s_.pop("data"); dn = "data/" + s_["seed_id"] + ("" if x.dtype == np.int32 else ".f32")
                if dn not in h:
                    h.create_dataset(dn, shape=(0,), maxshape=(None,), dtype=x.dtype, chunks=(16384,), compression="gzip",
                                     compression_opts=4, shuffle=True)
                ds = h[dn]; st_ = ds.shape[0]; ds.resize((st_ + x.size,)); ds[st_:] = x
                s_.update(dataset=dn, h5_start=st_, n=int(x.size), dtype=str(x.dtype))
        if segs: os.replace(fn + ".tmp", fn)
        else: os.remove(fn + ".tmp")
    return day, segs, idx, sdays


def month_of(day): return day[:7]


def cut(run, workers, months=None, outdir=None, max_days=0):
    import h5py, multiprocessing as mp
    P = PROFILES[run](); OUT = outdir or f"{RES}/cut/{run}"; os.makedirs(OUT, exist_ok=True)
    D = P["det"]
    D = D[(D.t >= pd.Timestamp(P["t0"])) & (D.t < pd.Timestamp(P["t1"]))].sort_values(["t", "fam"]).reset_index(drop=True)
    D.insert(0, "det_id", np.arange(len(D), dtype=np.int64))
    D["t_epoch"] = D.t.values.astype("datetime64[ns]").astype(np.int64) / 1e9
    D.to_parquet(f"{OUT}/detections.parquet")
    fams = sorted(D.fam.unique()); fcode = {f: i for i, f in enumerate(fams)}
    off = {}
    for r in P["off"].itertuples():
        if np.isfinite(r.offset_s): off.setdefault(r.netsta, {})[str(r.fam)] = float(r.offset_s)
    P["off"].to_csv(f"{OUT}/offsets_used.csv", index=False)
    osrc_txt = P["off"].offset_source.iloc[0] if len(P["off"]) else ""
    allsta = set(P["stations"])
    for r in P["roots"]:
        if os.path.isdir(r): allsta |= {s for s in os.listdir(r) if "." in s and os.path.isdir(f"{r}/{s}")}
    allsta = sorted(allsta)
    days = pd.date_range(P["t0"], P["t1"], freq="D", inclusive="left")
    days = [str(x.date()) for x in days]
    if months: days = [x for x in days if month_of(x) in set(months)]
    done = {os.path.basename(p)[5:12] for p in glob.glob(f"{OUT}/done_*.json")}
    days = [x for x in days if month_of(x) not in done]
    if max_days: days = days[:max_days]                                   # smoke tests only
    # stations with at least one file in range (for the channel-epoch query)
    ystas = sorted({ns for ns in allsta for r in P["roots"] for y in {x[:4] for x in days} if os.path.isdir(f"{r}/{ns}/{y}")})
    log(f"[cut] run={run} {len(D)} detections {P['t0']}..{P['t1']}; {len(days)} days to cut ({len(done)} months done); "
        f"{len(allsta)} stations known, {len(ystas)} with files in range; workers={workers}")
    E = channel_epochs(ystas, P["t0"], P["t1"], f"{OUT}/channel_epochs.csv"); EL = epoch_lookup(E)
    def resp_id(seed, t):
        for a, b, i in EL.get(seed, []):
            if a <= t < b: return i
        return f"{seed}@unknown"
    dday = D.t.dt.floor("D").dt.strftime("%Y-%m-%d").values
    grp = {k: v for k, v in pd.Series(np.arange(len(D))).groupby(dday).groups.items()}
    jobs = []
    for x in days:
        ii = np.asarray(grp.get(x, []), dtype=int)
        jobs.append((x, [(int(D.det_id.values[i]), fcode[D.fam.values[i]], float(D.t_epoch.values[i])) for i in ii]))
    detcols = ["det_id", "fam", "t"] + [c for c in P["stats"] if c in D.columns] + [c for c in P["flags"] if c in D.columns]
    DI = D[detcols].set_index("det_id")
    DAYFILES = os.environ.get("CUT_DAYFILES", "1") == "1"   # 2026-10-07: workers write windows_<day>.h5 (single gzip writer was the bottleneck)
    G = dict(allsta=allsta, roots=P["roots"], off=off, fams=fams, out=OUT if DAYFILES else "")

    # The index is built per day and streamed to index_<m>.parquet.tmp (one row group per day): holding a busy month's
    # ~100 M index rows as Python tuples until close_month was the next pod-OOM risk (2026-10-06).
    import pyarrow as pa, pyarrow.parquet as pq
    ICOLS = ["det_id", "netsta", "loc", "cha", "seed_id", "fs", "t0", "offset_s", "offset_is_template", "n_samples", "n_expected", "n_gap",
             "partial", "seg_id", "i_in_seg"]
    cur = dict(m=None); ISCHEMA = []
    def open_month(m):
        cur.update(m=m, t=time.time(), h5=None if DAYFILES else h5py.File(f"{OUT}/windows_{m}.h5.tmp", "w"), segrows=[], sdrows=[],
                   nseg=0, iw=None, nwin=0, nfull=0)
        if not DAYFILES: cur["h5"].attrs.update(run=run, month=m, pre_s=PRE_S, post_s=POST_S, created=str(pd.Timestamp.now("UTC")))
    def write_day_index(idx, segmeta, h5name):
        m = cur["m"]
        I = pd.DataFrame(idx, columns=ICOLS)
        I = I.merge(pd.DataFrame(segmeta, columns=["seg_id", "dataset", "h5_start", "resp_epoch"]), on="seg_id", how="left")
        I["h5_start"] = I.h5_start + I.i_in_seg; I = I.drop(columns=["i_in_seg"])
        I["offset_source"] = np.where(I.offset_is_template, osrc_txt, "detection_time"); I = I.drop(columns=["offset_is_template"])
        I["h5_file"] = h5name
        I = I.join(DI, on="det_id"); I["t_det"] = I.pop("t")
        if not ISCHEMA:
            sc = pa.Table.from_pandas(I, preserve_index=False).schema
            ISCHEMA.append(pa.schema([f.with_type(pa.large_string()) if pa.types.is_null(f.type) else f for f in sc], metadata=sc.metadata))
        if cur["iw"] is None: cur["iw"] = pq.ParquetWriter(f"{OUT}/index_{m}.parquet.tmp", ISCHEMA[0])
        cur["iw"].write_table(pa.Table.from_pandas(I, schema=ISCHEMA[0], preserve_index=False))
        cur["nwin"] += len(I); cur["nfull"] += int((~I.partial).sum())
    def close_month():
        m = cur["m"]
        if not DAYFILES: cur["h5"].close()
        S = pd.DataFrame(cur["segrows"])
        if cur["iw"] is not None:
            cur["iw"].close(); os.replace(f"{OUT}/index_{m}.parquet.tmp", f"{OUT}/index_{m}.parquet")
        else:
            pd.DataFrame(columns=[c for c in ICOLS if c not in ("offset_is_template", "i_in_seg")]).to_parquet(f"{OUT}/index_{m}.parquet")
        SD = pd.DataFrame(cur["sdrows"])
        if not DAYFILES: os.replace(f"{OUT}/windows_{m}.h5.tmp", f"{OUT}/windows_{m}.h5")
        S.drop(columns=[], errors="ignore").to_parquet(f"{OUT}/segments_{m}.parquet")
        SD.to_parquet(f"{OUT}/stationdays_{m}.parquet")
        gb = (sum(os.path.getsize(f) for f in glob.glob(f"{OUT}/windows_{m}-??.h5")) if DAYFILES
              else os.path.getsize(f"{OUT}/windows_{m}.h5")) / 1e9
        summ = dict(month=m, n_det=int(D.t.dt.strftime("%Y-%m").eq(m).sum()), n_windows=cur["nwin"], n_full=cur["nfull"],
                    n_station_days=len(SD), n_read_ok=int((SD.status == "ok").sum()) if len(SD) else 0,
                    raw_GB=round(SD.bytes.sum() / 1e9, 3) if len(SD) else 0, cut_GB=round(gb, 3),
                    n_segments=len(S), seg_hours=round(float((S.n / S.fs).sum() / 3600), 2) if len(S) else 0,
                    minutes=round((time.time() - cur["t"]) / 60, 1))
        json.dump(summ, open(f"{OUT}/done_{m}.json", "w"), indent=1)
        log(f"[cut] month {m}: {summ}")
        cur["m"] = None

    def bounded_imap(pool, fn, it, depth):
        # ordered like pool.imap, but at most `depth` days in flight: imap reads ahead without limit, so when the
        # single h5 writer here is slower than the workers the raw-day backlog filled RAM (pod OOM 2026-10-06, 2x).
        # Busy days hold ~6 GB of waveforms each: 2 x workers in flight reached 153 GB on 2010-09 -> default workers + 1.
        q = collections.deque()
        for j in it:
            q.append(pool.apply_async(fn, (j,)))
            if len(q) >= depth: yield q.popleft().get()
        while q: yield q.popleft().get()

    t00 = time.time()
    with mp.Pool(workers, initializer=_init, initargs=(G,)) as pool:
        for n, (day, segs, idx, sdays) in enumerate(bounded_imap(pool, cut_day, jobs, int(os.environ.get("CUT_INFLIGHT", workers + 1))), 1):
            m = month_of(day)
            if cur["m"] != m:
                if cur["m"] is not None: close_month()
                open_month(m)
            h5 = cur["h5"]; h5name = f"windows_{day}.h5" if DAYFILES else f"windows_{m}.h5"
            loc_ids = []; segmeta = []
            for s in segs:
                if DAYFILES:                     # the worker already wrote the samples
                    dn, st_, nx, dts = s["dataset"], s["h5_start"], s["n"], s["dtype"]
                else:
                    x = s["data"]; nm = s["seed_id"] + ("" if x.dtype == np.int32 else ".f32")
                    dn = f"data/{nm}"
                    if dn not in h5:
                        h5.create_dataset(dn, shape=(0,), maxshape=(None,), dtype=x.dtype, chunks=(16384,), compression="gzip",
                                          compression_opts=4, shuffle=True)
                    ds = h5[dn]; st_ = ds.shape[0]; ds.resize((st_ + x.size,)); ds[st_:] = x; nx, dts = x.size, str(x.dtype)
                sid = cur["nseg"]; cur["nseg"] += 1; loc_ids.append(sid)
                cur["segrows"].append(dict(seg_id=sid, seed_id=s["seed_id"], dataset=dn, h5_start=st_, n=nx, t0=s["t0"], fs=s["fs"],
                                           dtype=dts, gaps=json.dumps(s["gaps"]), src=s["src"],
                                           resp_epoch=resp_id(s["seed_id"], s["t0"]), h5_file=h5name))
                segmeta.append((sid, dn, st_, cur["segrows"][-1]["resp_epoch"]))
            if idx:
                write_day_index([r[:13] + (loc_ids[r[13]],) + tuple(r[14:]) for r in idx], segmeta, h5name)
            cur["sdrows"] += sdays
            if n % 20 == 0: log(f"  {day} {n}/{len(jobs)} days {time.time() - t00:.0f}s")
    if cur["m"] is not None: close_month()
    # combined summary + index
    Sm = pd.DataFrame([json.load(open(p)) for p in sorted(glob.glob(f"{OUT}/done_*.json"))])
    Sm.to_csv(f"{OUT}/cut_summary.csv", index=False)
    SD = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{OUT}/stationdays_*.parquet"))], ignore_index=True)
    SD.to_parquet(f"{OUT}/stationdays.parquet")
    log(f"[cut] done: {int(Sm.n_windows.sum())} windows, cut {Sm.cut_GB.sum():.2f} GB from raw {Sm.raw_GB.sum():.1f} GB "
        f"({len(Sm)} months) in {time.time() - t00:.0f}s")


# ================================================================================================ reading
def read_window(row, outdir):
    import h5py
    with h5py.File(f"{outdir}/{row['h5_file']}", "r") as h:
        x = h[row["dataset"]][int(row["h5_start"]):int(row["h5_start"]) + int(row["n_samples"])]
    return x, float(row["t0"]), float(row["fs"])


def load_index(outdir, columns=None, filters=None):
    return pd.concat([pd.read_parquet(p, columns=columns, filters=filters) for p in sorted(glob.glob(f"{outdir}/index_*.parquet"))],
                     ignore_index=True)


def load_index_sample(outdir, n, seed=0):
    """n random index rows over all months, reading only the row groups that hold them (a full load_index of a
    production store is billions of rows)."""
    import pyarrow.parquet as pq
    fs_ = sorted(glob.glob(f"{outdir}/index_*.parquet")); cnt = np.array([pq.ParquetFile(p).metadata.num_rows for p in fs_], dtype=np.int64)
    if not cnt.sum(): return pd.DataFrame()
    g = np.sort(np.random.default_rng(seed).choice(int(cnt.sum()), size=min(n, int(cnt.sum())), replace=False))
    off = np.r_[0, np.cumsum(cnt)]; out = []
    for k in np.unique(np.searchsorted(off, g, side="right") - 1):
        rows = g[(g >= off[k]) & (g < off[k + 1])] - off[k]
        pf = pq.ParquetFile(fs_[k]); rg = np.r_[0, np.cumsum([pf.metadata.row_group(i).num_rows for i in range(pf.num_row_groups)])]
        for i in np.unique(np.searchsorted(rg, rows, side="right") - 1):
            out.append(pf.read_row_group(i).to_pandas().iloc[rows[(rows >= rg[i]) & (rows < rg[i + 1])] - rg[i]])
    return pd.concat(out, ignore_index=True)


# ================================================================================================ verify
def verify(run, n=20, outdir=None, seed=0):
    from obspy import read, UTCDateTime
    OUT = outdir or f"{RES}/cut/{run}"
    I = load_index_sample(OUT, n, seed)        # sampled rows only (was: whole index, then I.sample)
    I["month"] = I.h5_file.str[8:15]                 # month files windows_<m>.h5 or day files windows_<day>.h5
    S = pd.concat([pd.read_parquet(f"{OUT}/segments_{m}.parquet", columns=["seg_id", "src"]).assign(month=m)
                   for m in sorted(I.month.unique())], ignore_index=True)
    pick = I.merge(S, on=["month", "seg_id"], how="left")
    rows = []
    for r in pick.to_dict("records"):
        x, t0, fs = read_window(r, OUT)
        st = read(r["src"]).select(id=r["seed_id"])
        # independent path: obspy slice of the raw traces at the window's first-sample time
        ref = np.zeros(x.size, x.dtype); have = np.zeros(x.size, bool)
        for tr in st:
            s = tr.slice(UTCDateTime(t0) - 0.25 / fs, UTCDateTime(t0 + (x.size - 1) / fs) + 0.25 / fs)
            if s.stats.npts == 0: continue
            p0 = int(round((s.stats.starttime.timestamp - t0) * fs)); m = min(s.stats.npts, x.size - p0)
            if m <= 0: continue
            ref[p0:p0 + m] = s.data[:m]; have[p0:p0 + m] = True
        eq = bool(np.array_equal(ref[have], x[have])); zeros_where_missing = bool(np.all(x[~have] == 0))
        rows.append(dict(det_id=r["det_id"], fam=r["fam"], seed_id=r["seed_id"], t0=str(UTCDateTime(t0)), fs=fs, n=x.size,
                         n_compared=int(have.sum()), n_gap_index=int(r["n_gap"]), n_missing_raw=int((~have).sum()),
                         identical=eq, gap_zero=zeros_where_missing, max_abs_diff=float(np.max(np.abs(ref[have].astype(float) - x[have]))) if have.any() else np.nan,
                         src=r["src"]))
    V = pd.DataFrame(rows); V.to_csv(f"{OUT}/verify_samples.csv", index=False)
    ok = bool(V.identical.all() and V.gap_zero.all() and (V.n_missing_raw == V.n_gap_index).all())
    msg = (f"{'PASS' if ok else 'FAIL'}: {len(V)} random windows re-read from HDF5 vs the raw miniSEED (obspy slice): "
           f"{int(V.identical.sum())} sample-identical, {int(V.n_compared.sum())} samples compared, "
           f"gap counts consistent {int((V.n_missing_raw == V.n_gap_index).sum())}/{len(V)}")
    open(f"{OUT}/VERIFY.txt", "w").write(msg + "\n"); log("[verify] " + msg + "\n" + V.drop(columns=["src"]).to_string())
    return ok


# ================================================================================================ stack check
def tcc(a, b, lag=5):
    out = []
    for c in range(min(a.shape[0], b.shape[0])):
        x, y = a[c].astype(float), b[c].astype(float)
        if not np.any(x) or not np.any(y): continue
        best = -1
        for s in range(-lag, lag + 1):
            u = x[max(0, s):len(x) + min(0, s)]; v = y[max(0, -s):len(y) + min(0, -s)]
            u = u - u.mean(); v = v - v.mean()
            best = max(best, float((u * v).sum() / (np.sqrt((u * u).sum() * (v * v).sum()) + 1e-30)))
        out.append(best)
    return float(np.mean(out)) if out else np.nan


def to_template_domain(x, fs, pre):
    """t3_cache chain (bp 1-8 Hz at native rate, resample_poly to 20 sps) + sh_v2 hp 2 Hz; returns the
    template slice [arrival - pre, arrival - pre + 6 s] (arrival at +15 s in the window)."""
    from scipy.signal import butter, sosfiltfilt, resample_poly
    from fractions import Fraction
    y = sosfiltfilt(butter(4, [1.0 / (fs / 2), 8.0 / (fs / 2)], btype="band", output="sos"), x.astype(float) - x.mean())
    fr = Fraction(20.0 / fs).limit_denominator(1000)
    if fr != 1: y = resample_poly(y, fr.numerator, fr.denominator)
    y = sosfiltfilt(butter(4, 2.0 / 10.0, btype="highpass", output="sos"), y)
    i0 = int(round((PRE_S - pre) * 20)); return y[i0:i0 + 120], y[int(round(1.0 * 20)):int(round(1.0 * 20)) + 120]


def stackcheck(run, outdir=None, fam=None, netsta=None, month=None):
    OUT = outdir or f"{RES}/cut/{run}"; P = PROFILES[run]()
    if not P["tpl"]:
        log("[stackcheck] no Shelly-style template set for this run; skipped"); return True
    M = pd.read_csv(P["tpl"] + ".csv"); M["fam"] = M.fam.astype(str); Z = np.load(P["tpl"] + ".npz"); pre = float(M.pre_s.iloc[0])
    I = load_index(OUT); I["fam"] = I.fam.astype(str); I["month"] = I.h5_file.str[8:15]
    flag = "in_cat_S1_9x" if "in_cat_S1_9x" in I.columns else P["flags"][0]
    J = I[(I[flag]) & (~I.partial) & (I.n_gap == 0) & (I.offset_source != "detection_time")]
    if fam is None or netsta is None or month is None:
        g = J.groupby(["fam", "netsta", "month"]).det_id.nunique().reset_index(name="n")
        g = g.merge(M[["fam", "netsta", "cc50"]], on=["fam", "netsta"]); g = g[g.n >= 30]
        g["score"] = g.cc50 * np.sqrt(g.n); b = g.sort_values("score").iloc[-1]
        fam, netsta, month = b.fam, b.netsta, b.month
    K = J[(J.fam == fam) & (J.netsta == netsta) & (J.month == month)]
    key = M[(M.fam == fam) & (M.netsta == netsta)].key.iloc[0]; T = Z[key]
    tchans = M[M.key == key].chans.iloc[0].split(",")
    stk = {c: np.zeros(120) for c in tchans}; mir = {c: np.zeros(120) for c in tchans}; nst = {c: 0 for c in tchans}; single = []
    for det, g in K.groupby("det_id"):
        per = {}
        for r in g.to_dict("records"):
            if r["cha"] not in tchans: continue
            x, t0, fs = read_window(r, OUT)
            if x.size < int(round(45 * fs)): continue
            w, m_ = to_template_domain(x, fs, pre); per[r["cha"]] = (w, m_)
        if not per: continue
        nrm = np.sqrt(np.mean([np.mean(w * w) for w, _ in per.values()])) + 1e-30
        for c, (w, m_) in per.items():
            stk[c] += w / nrm; mir[c] += m_ / nrm; nst[c] += 1
        A = np.vstack([per[c][0] if c in per else np.zeros(120) for c in tchans]); single.append(tcc(A, T))
    S_ = np.vstack([stk[c] for c in tchans]); Mi = np.vstack([mir[c] for c in tchans])
    cc = tcc(S_, T); ccm = tcc(Mi, T)
    res = dict(run=run, fam=fam, netsta=netsta, month=month, template=key, n_windows=int(max(nst.values())),
               cc_stack_vs_template=round(cc, 3), cc_single_median=round(float(np.nanmedian(single)), 3),
               cc_prearrival_stack_vs_template=round(ccm, 3), template_cc50=float(M[M.key == key].cc50.iloc[0]))
    res["pass"] = bool(cc >= 0.5 and cc > ccm + 0.2)
    json.dump(res, open(f"{OUT}/STACKCHECK.json", "w"), indent=1)
    log(f"[stackcheck] {res}")
    np.savez(f"{OUT}/stackcheck_arrays.npz", stack=S_, prearrival=Mi, template=T, chans=np.array(tchans))
    return res["pass"]


# ================================================================================================ delete-candidate list
def _duc_needs(days):
    """(net, sta, day) the Ducellier queue still reads from raw (months not yet cached)."""
    need = set(); FX = f"{RES}/duc_fix"
    if os.path.exists(f"{RES}/DUC_FIX_DONE"): return need, "DUC_FIX_DONE present: no raw reads left"
    try:
        import duc_fix_data as DD
    except Exception as e:
        return None, f"cannot import duc_fix_data ({e}); treat every file as needed"
    for m in sorted({x[:7] for x in days}):
        if os.path.exists(f"{FX}/data/cached_{m}"): continue
        P = DD.month_plan(m)
        need |= set(zip(P.net, P.sta, P.day))
    return need, "months without cached_<M> marker"


def _back_needs(days):
    need = set(); B = f"{RES}/sh_back"
    if os.path.exists(f"{RES}/SH_BACK_DONE"): return need, "SH_BACK_DONE present"
    p = f"{B}/download_plan_bk2.csv"
    if os.path.exists(p):
        P = pd.read_csv(p); need |= set(zip(P.net, P.sta, P.day.astype(str)))
    else:   # BK2 not planned yet: its window is 2008-01-01..2017-12-31 for the S6 stations (sh_back_queue.sh header)
        S6 = ["NC.KMPB", "NC.KMR", "PB.B049", "PB.B932", "PB.B933", "PB.B934"]
        for x in days:
            if "2008-01-01" <= x < "2018-01-01":
                for ns in S6: need.add((*ns.split("."), x))
    return need, "BK2 download plan (or its 2008-2017 S6 window when not yet planned)"


def _fwd_needs(days):
    need = set()
    if os.path.exists(f"{RES}/SH_EXT_DONE"): return need, "SH_EXT_DONE present"
    P = pd.read_csv(f"{RES}/sh_ext2026/download_plan.csv"); need |= set(zip(P.net, P.sta, P.day.astype(str)))
    return need, "forward queue still running"


def deletelist(run, outdir=None):
    OUT = outdir or f"{RES}/cut/{run}"; P = PROFILES[run]()
    SD = pd.read_parquet(f"{OUT}/stationdays.parquet")
    Sm = pd.read_csv(f"{OUT}/cut_summary.csv")
    days = [str(x.date()) for x in pd.date_range(P["t0"], P["t1"], freq="D", inclusive="left")]
    dset = set(days)
    # every raw file of the date range in the roots (also ones the cut did not see)
    files = []
    for r in P["roots"]:
        if not os.path.isdir(r): continue
        for ns in os.listdir(r):
            for y in sorted({x[:4] for x in days}):
                dd = f"{r}/{ns}/{y}"
                if not os.path.isdir(dd): continue
                for fn in os.listdir(dd):
                    if not fn.endswith(".mseed"): continue
                    day = str((pd.Timestamp(f"{y}-01-01") + pd.Timedelta(days=int(fn[:3]) - 1)).date())
                    if day in dset: files.append((f"{dd}/{fn}", ns, day))
    F = pd.DataFrame(files, columns=["path", "netsta", "day"])
    F["bytes"] = [os.path.getsize(p) for p in F.path]
    cutok = SD[SD.status == "ok"].set_index("path")
    needs = {}
    for name, fn in (("duc_fix", _duc_needs), ("sh_back", _back_needs), ("sh_ext", _fwd_needs)):
        if name == run: continue
        needs[name] = fn(days)
    rows = []; keep = []
    months_done = set(Sm.month)
    for r in F.itertuples():
        net, sta = r.netsta.split(".")
        why = []
        if r.day[:7] not in months_done: why.append("month not cut")
        elif r.path not in cutok.index:
            s = SD[SD.path == r.path]
            why.append("not fully cut (" + (s.status.iloc[0] if len(s) else "file not seen by the cut: another root holds this station-day") + ")")
        for name, (need, note) in needs.items():
            if need is None: why.append(f"{name}: {note}")
            elif (net, sta, r.day) in need: why.append(f"needed by {name} queue ({note})")
        if why:
            keep.append(dict(path=r.path, bytes=r.bytes, reason="; ".join(why)))
        else:
            c = cutok.loc[r.path]
            rows.append(dict(path=r.path, size_bytes=r.bytes, size_GB=round(r.bytes / 1e9, 4),
                             reason=f"fully cut ({int(c.n_win)} windows, {int(c.n_chan)} channels {c.band}); not read by duc_fix / sh_back / sh_ext queues"))
    C = pd.DataFrame(rows); K = pd.DataFrame(keep)
    lab = P["label"]
    C.to_csv(f"{OUT}/../DELETE_CANDIDATES_{lab}.csv", index=False); K.to_csv(f"{OUT}/../KEEP_{lab}.csv", index=False)
    tb = C.size_bytes.sum() / 1e12 if len(C) else 0.0; tbk = K.bytes.sum() / 1e12 if len(K) else 0.0
    # regenerable caches of this run (derived from the raw files; regenerable only while raw exists or by re-download)
    R = []
    if run == "sh_ext":
        Cx = pd.read_parquet(f"{RES}/sh_ext2026/cache_ext.parquet")
        for r in Cx[Cx.status.isin(["built", "cached"])].itertuples():
            n, s = r.netsta.split("."); p = f"{CACHE}/{r.netsta}/{pd.Timestamp(r.day).year}/{pd.Timestamp(r.day).dayofyear:03d}.npz"
            if os.path.exists(p):
                R.append(dict(path=p, size_bytes=os.path.getsize(p), kind="cache20hz (1-8 Hz, 20 sps)",
                              regenerate="python scripts/t3_cache.py results/sh_ext2026/download_plan.csv (needs the raw day file)"))
        for p in sorted(glob.glob(f"{RES}/sh_ext2026/prep/*.npz")):
            R.append(dict(path=p, size_bytes=os.path.getsize(p), kind="sh_v2 prep masks/refs",
                          regenerate="python scripts/sh_ext2026.py prep --days-file results/sh_ext2026/days/all.json (needs cache20hz)"))
    elif run == "sh_back":
        for p in sorted(glob.glob(f"{RES}/sh_back/prep_S6/*.npz")):
            R.append(dict(path=p, size_bytes=os.path.getsize(p), kind="sh_v2 prep", regenerate="sh_back_queue.sh prep_bk2 stage"))
    elif run == "duc_fix":
        for r in F.drop_duplicates(["netsta", "day"]).itertuples():
            p = f"{WS}/cache_duc/{r.netsta}/{r.day[:4]}/{pd.Timestamp(r.day).dayofyear:03d}.npz"
            if os.path.exists(p): R.append(dict(path=p, size_bytes=os.path.getsize(p), kind="cache_duc", regenerate="duc_fix_data.py --stage cache"))
    Rg = pd.DataFrame(R); Rg.to_csv(f"{OUT}/../REGENERABLE_CACHES_{lab}.csv", index=False)
    msg = (f"[deletelist] {run}: {len(C)} delete-candidate raw files, {tb:.3f} TB -> results/cut/DELETE_CANDIDATES_{lab}.csv; "
           f"kept {len(K)} files {tbk:.3f} TB (reasons in KEEP_{lab}.csv: "
           f"{K.reason.str.split(';').str[0].value_counts().head(6).to_dict() if len(K) else {}}); "
           f"regenerable caches {len(Rg)} files {Rg.size_bytes.sum() / 1e12 if len(Rg) else 0:.3f} TB -> REGENERABLE_CACHES_{lab}.csv; "
           f"cut store {Sm.cut_GB.sum():.1f} GB. NOTHING DELETED.")
    open(f"{OUT}/../DELETE_CANDIDATES_{lab}.txt", "w").write(msg + "\n" + "\n".join(f"  {k}: {v[1]}" for k, v in needs.items()) + "\n")
    log(msg)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("what", choices=["cut", "verify", "stackcheck", "deletelist"])
    ap.add_argument("--run", required=True, choices=list(PROFILES)); ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--months", default=""); ap.add_argument("--out", default=""); ap.add_argument("--n", type=int, default=20)
    ap.add_argument("--max-days", type=int, default=0, help="smoke test: cut only the first N days")
    ap.add_argument("--fam", default=None); ap.add_argument("--netsta", default=None); ap.add_argument("--month", default=None)
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    o = a.out or None
    if a.what == "cut": cut(a.run, a.workers, [m for m in a.months.split(",") if m] or None, o, a.max_days)
    elif a.what == "verify": sys.exit(0 if verify(a.run, a.n, o) else 5)
    elif a.what == "stackcheck": sys.exit(0 if stackcheck(a.run, o, a.fam, a.netsta, a.month) else 6)
    elif a.what == "deletelist": deletelist(a.run, o)
