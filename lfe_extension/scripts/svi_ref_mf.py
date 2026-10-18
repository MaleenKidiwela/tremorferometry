#!/usr/bin/env python
"""SVI_REF feature matched filter (GPU) -- the FROZEN Track C detector re-run with per-site / per-station gathers and
ONE scrambled control per day.

The detection path is c_run_month.py's measurement core copied line for line (same cache, templates, placement,
mergeB site merge, robust energy floor + dilation, per family-day network MAD, declustering, float32 order), so on
the same data and station set the 'data' rows reproduce results/trackC/detections_<ym>.parquet (t, fam, net_mad,
net_cc) exactly; smoke-tested in SVI_REF before launch.  What is added never touches a detection:

  data   every >= 8 x MAD network peak with: net_mad, net_cc, mad_scale, med_scale, rel_amp, con_pgc, wsh_pgc,
         n_sites, n_sta; per SITE (mergeB sites) the w-weighted site cc and its z against that site's own series at
         NRAND random valid origins of the same family-day (site_cc_max/2nd, site_z_max/2nd/min, n_z2/3/5, frac_z3,
         dom_top = largest site share of the network sum, site_noise_mad_med); per station cc_<STA>, z_<STA> and the
         S-window energy excess ez_<STA> (log 3-C energy over [S, S+4 s] at the predicted S = template start + 0.5 s,
         in MADs of the same station's values at the random origins); the Track E low-weight confirmation of the 2
         lowest-weight sites (conf_z, conf_ze; KEPT AS A VALIDATION INSTRUMENT, never a model feature); PGC HF ratio
         E(10-19)/E(2-8 Hz) from the native-rate trace (hfr) and hfr_rel = hfr / that station-day's median at 2000
         random times; cleaning flags (n_sta_clbad = contributing stations whose window touches a flat/muted sample
         of svi_ref_prep, n_sta_lines = contributing stations with notched lines that day); the station-set columns.
  c1     scrambled control: every site's cc series (and validity) rolled by that site's shift of the day
         (svi_ref_common.site_shifts), network sum rebuilt with the same weights, median / MAD, >= 8 x MAD, the
         same declustering; the same gathers at each site's ROLLED position (HF at the PGC site's rolled time).
  --clean 1   the southern cleaning is APPLIED (prep masks zero the data and mask every window that touches an
              invalid sample; lines notched); --clean 0 keeps the frozen data path and stores the flags only.
Station-day status (ok / nodata / short / missing / error) -> OUT/status/mf_<run>_<shard>.csv (guarded in the queue).

  python scripts/svi_ref_mf.py --run prod --kind prod --sset eff --days-file OUT/days/prod.json --shard 0 --nshard 3
  python scripts/svi_ref_mf.py --run t2011_S14 --kind truth --year 2011 --sset S14 --days-file ...
"""
from __future__ import annotations
import os, sys, json, time, argparse, warnings, traceback
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C
import c_run_month as CR
import svi_ref_common as R

FS = C.FS
PRE = C.PRE
E_WIN_S = 4.0
N_LOW = 2


def hf_prep(netsta, year, doy, note="", band=(10.0, 19.0)):
    """native-rate PGC day -> (t0, fs, prefix E(2-8), prefix E(HF)) -- CR.hf_ratio_day's chain, kept as prefix sums
    so the ratio can be read at data, control and random times.  None when unavailable."""
    p = C.trace_path(netsta, year, doy, note)
    if p is None:
        return None
    try:
        from obspy import read
        st = read(p)
        for tr in st:
            fr = round(float(tr.stats.sampling_rate))
            if fr > 0 and float(tr.stats.sampling_rate) != fr and abs(float(tr.stats.sampling_rate) - fr) / fr < 1e-4:
                tr.stats.sampling_rate = float(fr)
        st.merge(method=1, fill_value=0)
        chans = sorted({tr.stats.channel for tr in st})
        order = CR.chanorder.resolve_order(chans, station=f"{netsta}/{year}.{doy}", strict=True)
        keep = [c for c in order if c is not None]
        npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in keep)
        X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in keep])
        fs = float(st[0].stats.sampling_rate)
        t0 = float(st.select(channel=keep[0])[0].stats.starttime.timestamp)
        lo = (CR.bp_native(X, fs, 2.0, 8.0) ** 2).sum(0)
        hi = (CR.bp_native(X, fs, band[0], min(band[1], 0.95 * fs / 2)) ** 2).sum(0)
        return t0, fs, np.concatenate([[0.0], np.cumsum(lo)]), np.concatenate([[0.0], np.cumsum(hi)])
    except Exception as e:
        print(f"  [hfr] {netsta} {year}.{doy:03d} FAILED: {type(e).__name__}: {e}", flush=True)
        return None


def hf_at(H, tstarts, L):
    out = np.full(len(tstarts), np.nan)
    if H is None or not len(tstarts):
        return out
    t0, fs, clo, chi = H
    n = int(round(L * fs))
    i0 = np.rint((np.asarray(tstarts, float) - t0) * fs).astype(np.int64)
    g = (i0 >= 0) & (i0 + n < len(clo))
    den = clo[i0[g] + n] - clo[i0[g]]; num = chi[i0[g] + n] - chi[i0[g]]
    out[g] = np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)
    return out


def robust(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 20:
        return np.nan, np.nan
    m = float(np.median(v)); s = 1.4826 * float(np.median(np.abs(v - m)))
    return m, (s if s > 0 else np.nan)


def expect_data(mrow):
    """manifest says the raw day exists with >= 50 % coverage (c_month_ready's 'contributes' rule)."""
    if mrow is None or mrow[0] not in ("data", "present_elsewhere"):
        return False
    try:
        c = float(mrow[1])
    except (TypeError, ValueError):
        return True
    return (not np.isfinite(c)) or c >= 0.5


def load_excl(path):
    """R2 effective-set exclusions: {(sta, 'YYYY-MM'): reason}."""
    if not path or not os.path.exists(path):
        return {}
    E = pd.read_csv(path, dtype=str).fillna("")
    return {(r.sta, r.month): r.reason for r in E.itertuples()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True); ap.add_argument("--kind", default="prod", choices=["prod", "truth"])
    ap.add_argument("--year", type=int, default=0); ap.add_argument("--sset", default="eff")
    ap.add_argument("--clean", type=int, default=0); ap.add_argument("--days-file", default=""); ap.add_argument("--days", default="")
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--excl-file", default=""); ap.add_argument("--prep-dir", default="")
    ap.add_argument("--cache-dirs", default=""); ap.add_argument("--fam-chunk", type=int, default=24)
    ap.add_argument("--energy-floor", type=float, default=1e-3); ap.add_argument("--cc-device", default="auto")
    ap.add_argument("--gpu-free-gb", type=float, default=20.0); ap.add_argument("--gpu-need-gb", type=float, default=8.0)
    ap.add_argument("--gpu-wait-s", type=float, default=7200.0); ap.add_argument("--min-sites", type=int, default=2)
    ap.add_argument("--no-hfr", action="store_true"); ap.add_argument("--hf-sta", default="PGC")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    import svi_ref_prep as P
    days = json.load(open(a.days_file)) if a.days_file else [d for d in a.days.split(",") if d]
    days = [d for j, d in enumerate(sorted(days)) if j % a.nshard == a.shard]
    OD = f"{R.OUT}/feat/{a.run}"; R.ensure(OD, f"{R.OUT}/status")
    todo = [d for d in days if not os.path.exists(f"{OD}/{d}.json")]

    fz = C.frozen(); L = float(fz["template_length_s"])
    tpl_path, tpl_tag = C.resolve_templates(fz, "")
    tpl, sp, F, wgt, mv, ampref = C.load_frozen_tables(fz, tpl_path, L)
    stas_tbl = C.station_table()
    off = C.placement(tpl, sp, F, mv, stas_tbl, L)
    STAS = C.order_stations({k[1] for k in tpl})
    ep = C.epochs(); man = C.manifest()
    POLICY = C.epoch_policy()
    TSTA, LBL, DROPPED = {}, {}, {}
    def build_set(s, es, post):                     # c_run_month.build_set verbatim (no per-epoch npz in use)
        pol = POLICY.get((s, es))
        if pol and pol["drop"]:
            TSTA[(s, es)] = {}; LBL[(s, es)] = "dropped"; DROPPED[(s, es)] = True; return
        f_ovr = C.epoch_template_file(s, es, POLICY)
        if f_ovr is not None:
            raise SystemExit(f"[mf] per-epoch template file {f_ovr} in use: not supported by SVI_REF (v1.1 has none)")
        TSTA[(s, es)] = {k[0]: (tpl[k][L], off[k]) for k in tpl if k[1] == s and k in off}
        LBL[(s, es)] = (pol["label"] if (pol and pol["label"]) else ("2012_provisional" if post else "2012_final"))
    for s in STAS:
        seen = set()
        for r in ep[ep.sta == s].itertuples():
            for d0 in (r.start, pd.Timestamp(C.POST_2012_BOUNDARY[s]) if s in C.POST_2012_BOUNDARY else None):
                if d0 is None or not (r.start <= d0 <= r.end): continue
                es = C.eff_epoch_start(s, d0, ep)
                if es in seen: continue
                seen.add(es); build_set(s, es, C.is_post_boundary(s, d0))
        if not seen: build_set(s, "base", False)
    NETSTA = {s: (ep[ep.sta == s].netsta.iloc[0] if (ep.sta == s).any() else str(stas_tbl.loc[s, "netsta"])) for s in STAS}
    SITE = C.site_map(STAS)
    share_ok = C.share_ok_matrix()
    if a.cache_dirs: CACHES = [d.rstrip("/") for d in a.cache_dirs.split(",") if d]
    elif a.kind == "truth": CACHES = list(R.TRUTH[a.year]["caches"])
    else: CACHES = [R.CACHE_R1, C.CACHE]
    SSET = None if a.sset == "eff" else set(R.SSETS[a.sset])
    EXCL = load_excl(a.excl_file) if a.sset == "eff" else {}
    n = int(round(L * FS)); Lmax = n; nE = int(round(E_WIN_S * FS)); PREn = int(round(PRE * FS))
    mann = {(r.netsta, r.date): (r.status, r.coverage_frac, r.note if isinstance(r.note, str) else "") for r in man.itertuples()}
    R.log(f"[mf {a.run} {a.shard}/{a.nshard}] {len(todo)} of {len(days)} days; kind={a.kind} sset={a.sset} clean={a.clean} "
          f"caches={CACHES} excl={len(EXCL)} stations={len(STAS)} templates={tpl_tag}")

    import torch
    import torch.nn.functional as Fnn
    from concurrent.futures import ThreadPoolExecutor
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    torch.set_num_threads(4)
    CCDEV = [dev]
    if dev == "cuda":
        free_b, tot_b = torch.cuda.mem_get_info()
        if a.cc_device == "cpu" or (a.cc_device == "auto" and free_b / 2 ** 30 < a.gpu_free_gb): CCDEV[0] = "cpu"
    FCH = [a.fam_chunk]
    EX = ThreadPoolExecutor(1)

    def wait_gpu():
        if dev != "cuda": return
        t0w = time.time()
        while True:
            torch.cuda.empty_cache(); free_b = torch.cuda.mem_get_info()[0]
            if free_b / 2 ** 30 >= a.gpu_need_gb or time.time() - t0w > a.gpu_wait_s:
                if a.cc_device == "auto":
                    CCDEV[0] = "cuda" if free_b / 2 ** 30 >= a.gpu_free_gb else "cpu"
                return
            time.sleep(20)

    sp_ = f"{R.OUT}/status/mf_{a.run}_s{a.shard}.csv"
    status = [tuple(x) for x in pd.read_csv(sp_, dtype=str).fillna("").itertuples(index=False)] if os.path.exists(sp_) else []
    t_all = time.time()
    for nday, ds in enumerate(todo):
        day = pd.Timestamp(ds); doy = int(day.dayofyear); year = int(day.year); ym = ds[:7]
        t_day = time.time()
        try:
            wait_gpu()
            al = C.alive(ep, day)
            alive_sta = {r.sta: r.epoch_id for r in al.itertuples()}
            EK = {s: C.eff_epoch_start(s, day, ep) for s in STAS if s in alive_sta}
            TS = {s: TSTA.get((s, EK[s]), TSTA.get((s, "base"), {})) for s in EK}
            hfut = None
            if not a.no_hfr and a.hf_sta in NETSTA and a.hf_sta in alive_sta and (SSET is None or a.hf_sta in SSET):
                hfut = EX.submit(hf_prep, NETSTA[a.hf_sta], year, doy, mann.get((NETSTA[a.hf_sta], ds), ("", 0, ""))[2])
            DP = P.DayPrep(ds, a.prep_dir or None)
            if a.clean and not DP.ok():
                raise RuntimeError(f"--clean 1 but no prep file for {ds}")
            D, T0, NPT, EPID, CLINV, NLINES = {}, {}, {}, {}, {}, {}
            for s in STAS:
                if s not in alive_sta: continue
                ns = NETSTA[s]
                if SSET is not None and s not in SSET: continue
                if DROPPED.get((s, EK[s])): continue                     # epoch policy (LRIV >= 2019-05-30)
                if (s, ym) in EXCL:
                    status.append((ds, ns, "excluded", EXCL[(s, ym)])); continue
                p = next((q for q in (f"{c}/{ns}/{year}/{doy:03d}.npz" for c in CACHES) if os.path.exists(q)), None)
                mrow = mann.get((ns, ds))
                expect = expect_data(mrow)
                if p is None:
                    status.append((ds, ns, "missing" if (expect and a.kind == "prod") else "nodata",
                                   "no cache file" + (f" (manifest {mrow[0]})" if mrow else ""))); continue
                try:
                    z = np.load(p); x = z["x"]; t0s = float(z["t0"])
                    if x.ndim != 2 or x.shape[0] != 3 or x.shape[1] < 0.5 * 86400 * FS:
                        status.append((ds, ns, "short", f"cache shape {x.shape}")); continue
                except Exception as e:
                    status.append((ds, ns, "error", f"cache unreadable: {str(e)[:60]}")); continue
                inv = None
                if DP.ok():
                    if a.clean:
                        x, inv = DP.apply(ns, x, t0s)
                    else:
                        inv = DP.mask(ns, t0s, x.shape[1])
                    NLINES[s] = len(DP.lines(ns))
                D[s] = torch.from_numpy(np.ascontiguousarray(x)).to(dev)
                T0[s] = t0s; NPT[s] = D[s].shape[1]; EPID[s] = alive_sta[s]; CLINV[s] = inv
                status.append((ds, ns, "ok", ""))
            meta = dict(date=ds, run=a.run, kind=a.kind, sset=a.sset, clean=a.clean, stations=sorted(D),
                        prep=DP.ok(), excluded=sorted(s for s in STAS if (s, ym) in EXCL))
            if not D:
                pd.DataFrame().to_parquet(f"{OD}/{ds}.parquet"); json.dump(meta, open(f"{OD}/{ds}.json", "w")); continue
            net_epoch = "|".join(f"{s}:{EPID[s].split(':')[1]}" for s in sorted(D))
            Nmin = min(NPT.values())
            NF = int(2 ** np.ceil(np.log2(Nmin + Lmax)))
            NV = Nmin - n
            DH, EN = {}, {}
            for s in D:
                x = D[s][:, :Nmin]
                DH[s] = torch.fft.rfft(x, n=NF)
                xd = x.double()
                EN[s] = torch.cumsum(torch.cat([torch.zeros((x.shape[0], 1), device=dev, dtype=torch.float64), xd * xd], 1), 1)
            fam_list = sorted({f for s in D for f in TS.get(s, {})})
            S3, VALID, E3W, E4W, CLW = {}, {}, {}, {}, {}
            for s in D:
                fs_ = [f for f in fam_list if f in TS.get(s, {})]
                if not fs_: continue
                try:
                    out3 = torch.empty((len(fs_), NV), device=CCDEV[0])
                except torch.OutOfMemoryError:
                    torch.cuda.empty_cache(); CCDEV[0] = "cpu"; out3 = torch.empty((len(fs_), NV), device="cpu")
                ccdev = CCDEV[0]
                ew = (EN[s][:, n:Nmin + 1] - EN[s][:, :NV + 1])[:, :NV]
                med_all = ew.median(); med_ch = ew.median(dim=1).values
                e3 = ew.sum(0); med3 = e3.median()
                ref_ch = torch.maximum(med_ch, 1e-2 * med3)
                bad_w = ((ew < a.energy_floor * ref_ch[:, None]).any(0) | (e3 < a.energy_floor * med3) | (ew.min(0).values <= 0))
                if float(med3) <= 0: bad_w = torch.ones_like(bad_w)
                if bool(bad_w.any()):
                    bad_w = (Fnn.max_pool1d(bad_w.float()[None, None], kernel_size=2 * n + 1, stride=1, padding=n)[0, 0] > 0)
                # cleaning: window touches an invalid sample (exact, no extra dilation) -- applied (clean 1) or flag only
                clw = None
                if CLINV.get(s) is not None:
                    inv = CLINV[s][:Nmin]
                    ci = np.concatenate([[0], np.cumsum(inv.astype(np.int64))])
                    clw = torch.from_numpy((ci[n:n + NV] - ci[:NV]) > 0).to(dev)
                    if a.clean:
                        bad_w = bad_w | clw
                CLW[s] = clw
                efloor = 1e-10 * med_all
                e = ew.clamp_min(efloor).sqrt().float()
                ch = FCH[0]; i0 = 0
                while i0 < len(fs_):
                    ff = fs_[i0:i0 + ch]
                    try:
                        W = torch.from_numpy(np.stack([TS[s][f][0] for f in ff])).to(dev)
                        WH = torch.fft.rfft(W, n=NF)
                        num = torch.fft.irfft(WH.conj() * DH[s][None], n=NF)[:, :, :NV]
                        out3[i0:i0 + len(ff)] = (num / e[None]).mean(1).to(ccdev)
                        del W, WH, num
                    except torch.OutOfMemoryError:
                        torch.cuda.empty_cache()
                        if ch == 1: raise
                        ch = max(1, ch // 2); FCH[0] = ch; continue
                    i0 += len(ff)
                if bool(bad_w.any()): out3[:, bad_w] = 0.0
                S3[s] = (fs_, out3); VALID[s] = ~bad_w; E3W[s] = e3.float()
                e4 = (EN[s][:, nE:Nmin + 1] - EN[s][:, :Nmin + 1 - nE]).sum(0)
                E4W[s] = torch.log(e4.clamp_min(1e-30)).float()
                del ew, e3, e, bad_w, e4
                del DH[s], EN[s]
                if dev == "cuda": torch.cuda.empty_cache()
            DH.clear(); EN.clear()
            H = hfut.result() if hfut is not None else None
            hf_fs = H[1] if H is not None else np.nan
            hf_med = np.nan
            if H is not None:
                rr = np.random.default_rng(7 + doy + 1000 * year)
                tr_ = H[0] + rr.uniform(0, (len(H[2]) - 1) / H[1] - L, 2000)
                hf_med = float(np.nanmedian(hf_at(H, tr_, L)))
            SH = R.site_shifts(ds, NV)
            meta["shifts"] = SH
            out_rows, fd_rows = [], []
            for f in sorted({f2 for s in S3 for f2 in TS.get(s, {})}):
                ss = [s for s in D if s in S3 and f in TS.get(s, {})]
                if not ss: continue
                groups = {}
                for s in ss: groups.setdefault(SITE[s], []).append(s)
                base_fd = dict(fam=int(f), date=ds, n_sites=len(groups), n_sta=len(ss), sites=",".join(sorted(groups)))
                if len(groups) < a.min_sites:
                    fd_rows.append(dict(series="data", observable=False, **base_fd)); continue
                o = {s: TS[s][f][1] for s in ss}
                tref0 = min(T0[s] for s in ss)
                k = {s: int(round((tref0 - T0[s] + o[s]) * FS)) for s in ss}
                b = max(0, -min(k.values())); k = {s: k[s] + b for s in ss}
                tref = tref0 + b / FS
                NVo = NV - max(k.values())
                if NVo <= 1000:
                    fd_rows.append(dict(series="data", observable=False, **base_fd)); continue
                wm = {s: 3.0 * float(wgt.get((f, s), 1.0)) for s in ss}
                Wsite = {q: float(np.sqrt(sum(wm[s] ** 2 for s in mem))) for q, mem in groups.items()}
                acc_n = torch.zeros(NVo, device=dev)
                vmask = torch.ones(NVo, dtype=torch.bool, device=dev)
                weff, wtot, SEG, CS, VM = {}, 0.0, {}, {}, {}
                for q, mem in groups.items():
                    wsum_m = sum(wm[s] for s in mem)
                    cs = torch.zeros(NVo, device=dev)
                    vm = torch.ones(NVo, dtype=torch.bool, device=dev)
                    for s in mem:
                        fs_, cc3 = S3[s]; i = fs_.index(f)
                        seg = cc3[i, k[s]:k[s] + NVo]
                        if seg.device.type != dev: seg = seg.to(dev, non_blocking=True)
                        SEG[s] = seg
                        cs += (wm[s] / wsum_m) * seg
                        vm &= VALID[s][k[s]:k[s] + NVo]
                        weff[s] = Wsite[q] * wm[s] / wsum_m
                    acc_n += Wsite[q] * cs
                    vmask &= vm
                    wtot += Wsite[q]
                    CS[q] = cs; VM[q] = vm
                acc_n /= wtot
                am = acc_n[vmask]
                if am.numel() < 1000:
                    fd_rows.append(dict(series="data", observable=False, **base_fd)); del acc_n, vmask; continue
                medn = am.median(); madn = (am - medn).abs().median().clamp_min(1e-6)
                idx = CR.decluster_idx(acc_n, medn + R.REC_THR * madn, n, vmask)
                wsum = sum(weff[s] for s in ss)
                sites = sorted(groups); Wv = np.array([Wsite[q] for q in sites])
                low = [sites[j] for j in np.argsort(Wv / Wv.sum(), kind="stable")[:N_LOW]] if len(sites) > N_LOW else []
                # ---- random valid origins -> per-site / per-station background
                vidx = torch.nonzero(vmask).squeeze(1).cpu().numpy()
                RI = R.rand_idx(ds, f, vidx)
                RIt = torch.from_numpy(RI).to(dev)
                bg_site = {q: robust(CS[q][RIt].cpu().numpy()) for q in sites}
                bg_sta = {s: robust(SEG[s][RIt].cpu().numpy()) for s in ss}
                bg_e = {}
                for s in ss:
                    j = RIt + k[s] + PREn
                    ok = j < E4W[s].numel()
                    bg_e[s] = robust(E4W[s][j[ok]].cpu().numpy())
                nl_sta = sum(1 for s in ss if NLINES.get(s, 0) > 0)

                def gather(II, series, nmad, vals, med_s, mad_s):
                    """II: {site: index array (len m)} -> feature DataFrame."""
                    m = len(next(iter(II.values())))
                    cols = {}
                    SC = np.full((len(sites), m), np.nan); SZ = np.full((len(sites), m), np.nan)
                    SHR = np.full((len(sites), m), np.nan); SE = np.full((len(sites), m), np.nan)
                    amp = np.zeros(m); cpg = np.zeros(m); clb = np.zeros(m, np.int16)
                    for j_, q in enumerate(sites):
                        it = torch.from_numpy(II[q]).to(dev)
                        cv = CS[q][it].cpu().numpy().astype(float)
                        SC[j_] = cv; SHR[j_] = Wsite[q] * cv
                        mq, sq = bg_site[q]; SZ[j_] = (cv - mq) / sq
                        eq_ = []
                        for s in groups[q]:
                            ccv = SEG[s][it]
                            js = it + k[s]
                            amp += (weff[s] * ccv * torch.sqrt(E3W[s][js] / ampref[(s, L)])).cpu().numpy()
                            if SITE[s] == "PGC": cpg += (weff[s] * ccv).cpu().numpy()
                            c_ = ccv.cpu().numpy().astype(float); ms_, ss_ = bg_sta[s]
                            cols[f"cc_{s}"] = c_.astype(np.float32); cols[f"z_{s}"] = ((c_ - ms_) / ss_).astype(np.float32)
                            je = (js + PREn).clamp_max(E4W[s].numel() - 1)
                            ev = E4W[s][je].cpu().numpy().astype(float); me_, se_ = bg_e[s]
                            ez = (ev - me_) / se_; cols[f"ez_{s}"] = ez.astype(np.float32); eq_.append(ez)
                            if CLW.get(s) is not None:
                                clb += CLW[s][js.clamp_max(CLW[s].numel() - 1)].cpu().numpy().astype(np.int16)
                        SE[j_] = np.nanmean(np.vstack(eq_), 0) if eq_ else np.nan
                    with np.errstate(invalid="ignore", divide="ignore"):
                        tot = np.nansum(SHR, 0)
                        dom = np.nanmax(SHR, 0) / np.where(tot > 0, tot, np.nan)
                        Cs = np.sort(np.nan_to_num(SC, nan=-99), 0)[::-1]; Zs = np.sort(np.nan_to_num(SZ, nan=-99), 0)[::-1]
                        top = lambda A, r: np.where(A[r] > -98, A[r], np.nan) if A.shape[0] > r else np.full(m, np.nan)
                        zmin = np.nanmin(np.where(np.isfinite(SZ), SZ, np.inf), 0); zmin[~np.isfinite(zmin)] = np.nan
                        lo_i = [sites.index(q) for q in low]
                        conf_z = np.nanmean(SZ[lo_i], 0) if lo_i else np.full(m, np.nan)
                        conf_ze = np.nanmean(SE[lo_i], 0) if lo_i else np.full(m, np.nan)
                    tvals = tref + II["__net"] / FS
                    Dd = pd.DataFrame(dict(series=series, fam=int(f), date=ds, t=tvals, i=II["__net"].astype(np.int64),
                                           net_mad=nmad, net_cc=vals, mad_scale=float(mad_s), med_scale=float(med_s),
                                           rel_amp=amp / wsum, con_pgc=cpg / wsum,
                                           wsh_pgc=float(sum(weff[s] for s in ss if SITE[s] == "PGC") / wsum),
                                           n_sites=len(groups), n_sta=len(ss), n_sta_day=len(D), epoch=net_epoch,
                                           has_LRIV="LRIV" in ss, has_SNB="SNB" in ss,
                                           site_cc_max=top(Cs, 0), site_cc_2nd=top(Cs, 1), site_z_max=top(Zs, 0),
                                           site_z_2nd=top(Zs, 1), site_z_min=zmin,
                                           n_z2=(np.nan_to_num(SZ, nan=-9) >= 2).sum(0), n_z3=(np.nan_to_num(SZ, nan=-9) >= 3).sum(0),
                                           n_z5=(np.nan_to_num(SZ, nan=-9) >= 5).sum(0), dom_top=dom,
                                           site_noise_mad_med=float(np.nanmedian([bg_site[q][1] for q in sites])),
                                           conf_z=conf_z, conf_ze=conf_ze, conf_sites=",".join(low),
                                           n_sta_clbad=clb, n_sta_lines=nl_sta))
                    Dd["frac_z3"] = Dd.n_z3 / max(len(sites), 1)
                    for c_, v_ in cols.items(): Dd[c_] = v_
                    # HF at PGC (station PGC itself, its own window = site PGC index + k[PGC])
                    if "PGC" in ss and H is not None:
                        ip = II[SITE["PGC"]]
                        tst = tref + ip / FS + o["PGC"]
                        hv = hf_at(H, tst, L)
                    else:
                        hv = np.full(m, np.nan)
                    Dd["hfr"] = hv; Dd["hfr_rel"] = hv / hf_med if np.isfinite(hf_med) and hf_med > 0 else np.nan
                    Dd["hf_fs"] = hf_fs; Dd["prep_ok"] = DP.ok()
                    return Dd

                if len(idx):
                    ii = torch.from_numpy(idx).to(dev); vals = acc_n[ii]
                    nmad = ((vals - medn) / madn).cpu().numpy()
                    II = {q: idx for q in sites}; II["__net"] = idx
                    out_rows.append(gather(II, "data", nmad, vals.cpu().numpy(), float(medn), float(madn)))
                fd_rows.append(dict(series="data", observable=True, mad=float(madn), med=float(medn), tref=tref, NVo=int(NVo),
                                    n_valid_s=float(am.numel() / FS), n_rec=int(len(idx)), **base_fd))
                # ---- scrambled control: each site rolled by its own shift of the day
                acc_c = torch.zeros(NVo, device=dev); vmc = torch.ones(NVo, dtype=torch.bool, device=dev)
                shq = {q: int(SH[q]) % NVo for q in sites}
                for q in groups:
                    acc_c += Wsite[q] * torch.roll(CS[q], -shq[q])
                    vmc &= torch.roll(VM[q], -shq[q])
                acc_c /= wtot
                amc = acc_c[vmc]
                if amc.numel() >= 1000:
                    medc = amc.median(); madc = (amc - medc).abs().median().clamp_min(1e-6)
                    idc = CR.decluster_idx(acc_c, medc + R.REC_THR * madc, n, vmc)
                    if len(idc):
                        vc = acc_c[torch.from_numpy(idc).to(dev)]
                        II = {q: (idc + shq[q]) % NVo for q in sites}; II["__net"] = idc
                        out_rows.append(gather(II, "c1", ((vc - medc) / madc).cpu().numpy(), vc.cpu().numpy(),
                                               float(medc), float(madc)))
                    fd_rows.append(dict(series="c1", observable=True, mad=float(madc), med=float(medc), tref=tref, NVo=int(NVo),
                                        n_valid_s=float(amc.numel() / FS), n_rec=int(len(idc)), **base_fd))
                else:
                    fd_rows.append(dict(series="c1", observable=False, **base_fd))
                del acc_n, vmask, am, SEG, CS, VM, acc_c, vmc, amc
            del S3, VALID, E3W, E4W, D, CLW
            if dev == "cuda": torch.cuda.empty_cache()
            Fo = pd.concat(out_rows, ignore_index=True) if out_rows else pd.DataFrame()
            if len(Fo):
                Fo["keep_xfam8"] = False
                for ser, g in Fo.groupby("series"):
                    Fo.loc[g.index, "keep_xfam8"] = C.local_max_keep_cc(g.t.to_numpy(), g.net_cc.to_numpy().astype(float),
                                                                        g.fam.to_numpy().astype(int), share_ok, 2.0)
                Fo["sset"] = a.sset; Fo["run"] = a.run; Fo["k"] = np.arange(len(Fo))
            FD = pd.DataFrame(fd_rows)
            Fo.to_parquet(f"{OD}/{ds}.parquet.tmp", index=False); os.replace(f"{OD}/{ds}.parquet.tmp", f"{OD}/{ds}.parquet")
            FD.to_parquet(f"{OD}/{ds}_fd.parquet.tmp", index=False); os.replace(f"{OD}/{ds}_fd.parquet.tmp", f"{OD}/{ds}_fd.parquet")
            meta.update(n_data=int((Fo.series == "data").sum()) if len(Fo) else 0, n_c1=int((Fo.series == "c1").sum()) if len(Fo) else 0,
                        hf_fs=hf_fs, hf_med=hf_med, NV=NV, secs=round(time.time() - t_day, 1))
            json.dump(meta, open(f"{OD}/{ds}.json.tmp", "w"), default=str); os.replace(f"{OD}/{ds}.json.tmp", f"{OD}/{ds}.json")
            if nday < 3 or nday % 20 == 0:
                R.log(f"  [{a.run} {a.shard}] {ds}: data {meta['n_data']} c1 {meta['n_c1']} stations {len(meta['stations'])} "
                      f"{meta['secs']}s; {nday + 1}/{len(todo)} {time.time() - t_all:.0f}s")
        except Exception as e:
            R.log(f"  ERROR {ds}: {e}\n{traceback.format_exc()[-1500:]}")
            status += [(ds, NETSTA[s], "error", f"day job: {str(e)[:100]}") for s in STAS if s in C.alive(ep, day).sta.values]
            if dev == "cuda": torch.cuda.empty_cache()
        if nday % 10 == 0: R.write_status(f"mf_{a.run}", f"s{a.shard}", status)
    R.write_status(f"mf_{a.run}", f"s{a.shard}", status)
    R.log(f"[mf {a.run} {a.shard}] done {len(todo)} days {time.time() - t_all:.0f}s")


if __name__ == "__main__":
    main()
