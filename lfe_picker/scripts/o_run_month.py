#!/usr/bin/env python
"""Track C step 3 — the PRODUCTION network matched filter for one month, on the FROZEN chain.

This is e4_template_test.py's measurement core with everything that only exists for the referee
experiment (Bostock ground truth, halves, random controls, rolled controls, Z-only series, moveout
pass 1, event scoring) removed, and nothing else changed.  Every numerical step below is the same
code path as the frozen a1_mergeB_L12 run (scripts/a1_run_mergeB_L12.sh):

  day file       cache20hz/<NET.STA>/<YYYY>/<DOY>.npz  (20 Hz, 1-8 Hz Butterworth-4 zero-phase,
                 component order Z / N|1 / E|2, not rotated)
  Nmin           min npts over the stations present that day; NF = 2^ceil(log2(Nmin + Lmax))
  cc             full-day normalised correlation per channel by rFFT, denominator from a FLOAT64
                 prefix sum of squares (the E4 bugfix), station series = mean over its 3 channels
  valid          3-C window energy > 1e-6 x the day's median AND every channel energy > 0
  placement      template start = origin + ttP(Crosson 1-D) + (S-P) - 0.5 s + empirical moveout
  network        mergeB: within-site w-weighted MEAN cc for {PGC,B009,B010,B011} and
                 {LRIV,B005,B006,B007}; site weight sqrt(sum w^2); w = 3 x SNR weight; sum over
                 sites / total weight
  background     median and MAD of THAT family-day's own network series over the valid samples
  detections     declustered peaks (max of a template-length window) above med + T_rec x MAD
  amplitude      SNR-weighted mean of cc_s x sqrt(E3_s / E3_ref_s)  (units of the station's 2012
                 quiet-day window rms)

Production rules added on top (Merlin Step 7):
  - >= 2 SITES with a usable template AND data, else the family-day is `observable = False` and is
    NOT scored (recorded in family_days_<YYYY-MM>.parquet);
  - both tiers recorded, never re-run: every detection carries net_mad and the family-day MAD scale,
    so T_core (8.75) and T_ext (8.25) are filters on one parquet;
  - frozen cross-family assignment: among families sharing >= 50 % of their (2012) stations a
    detection is dropped if a station-sharing family has a higher RAW weighted cc within +-2 s;
    the max-MAD assignment is stored alongside as a column, never as the default;
  - PNSN tremor-in-box flag, catalogued-earthquake flag (60 s), HF-ratio column at the PGC site.

Outputs (results/trackC/):
  detections_<YYYY-MM>.parquet, family_days_<YYYY-MM>.parquet, stacks_<YYYY-MM>.npz, run_<YYYY-MM>.json

  python scripts/c_run_month.py 2015-12
"""
from __future__ import annotations
import os, sys, json, time, argparse, warnings
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import o_common as C
import o_chanorder as chanorder

FS = C.FS
PRE = C.PRE


def decluster_idx(x, thr, n, valid):
    """declustered peak indices of a 1-D series (max of a template-length window, masked samples out)."""
    import torch
    import torch.nn.functional as Fnn
    mx = Fnn.max_pool1d(x[None, None], kernel_size=2 * (n // 2) + 1, stride=1, padding=n // 2)[0, 0]
    return torch.nonzero((x >= mx - 1e-12) & valid & (x >= thr)).squeeze(1).cpu().numpy()


def bp_native(x, fs, f1, f2):
    from scipy.signal import butter, sosfiltfilt
    return sosfiltfilt(butter(4, [f1 / (fs / 2), f2 / (fs / 2)], btype="band", output="sos"), x, axis=-1)


def hf_ratio_day(netsta, year, doy, tstarts, L, band=(10.0, 19.0), note=""):
    """E(HF)/E(2-8 Hz) over the scored window at the HF station, from the NATIVE-rate day trace
    (identical to a1_hf_veto.py; the 20 Hz cache cannot carry the band)."""
    out = np.full(len(tstarts), np.nan)
    fs = np.nan
    p = C.trace_path(netsta, year, doy, note)
    if p is None or len(tstarts) == 0:
        return out, fs
    try:
        from obspy import read
        st = read(p)
        for tr in st:                       # same CN fragment-rate snap as c_cache_month (see README)
            fr = round(float(tr.stats.sampling_rate))
            if fr > 0 and float(tr.stats.sampling_rate) != fr and \
               abs(float(tr.stats.sampling_rate) - fr) / fr < 1e-4:
                tr.stats.sampling_rate = float(fr)
        st.merge(method=1, fill_value=0)
        chans = sorted({tr.stats.channel for tr in st})
        order = chanorder.resolve_order(chans, station=f"{netsta}/{year}.{doy}", strict=True)
        keep = [c for c in order if c is not None]
        npts = min(min(tr.stats.npts for tr in st.select(channel=c)) for c in keep)
        X = np.vstack([st.select(channel=c)[0].data[:npts].astype(np.float64) for c in keep])
        fs = float(st[0].stats.sampling_rate)
        t0 = float(st.select(channel=keep[0])[0].stats.starttime.timestamp)
        lo = (bp_native(X, fs, 2.0, 8.0) ** 2).sum(0)
        hi = (bp_native(X, fs, band[0], min(band[1], 0.95 * fs / 2)) ** 2).sum(0)
        clo = np.concatenate([[0.0], np.cumsum(lo)]); chi = np.concatenate([[0.0], np.cumsum(hi)])
        n = int(round(L * fs))
        i0 = np.rint((np.asarray(tstarts, float) - t0) * fs).astype(np.int64)
        g = (i0 >= 0) & (i0 + n <= len(lo))
        den = clo[i0[g] + n] - clo[i0[g]]; num = chi[i0[g] + n] - chi[i0[g]]
        out[g] = np.where(den > 0, num / np.maximum(den, 1e-30), np.nan)
    except Exception as e:
        print(f"  [hfr] {netsta} {year}.{doy:03d} FAILED: {type(e).__name__}: {e}", flush=True)
    return out, fs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ym")
    ap.add_argument("--templates", default="", help="override the template npz (default: frozen / TRACK_A_RESULTS.md)")
    ap.add_argument("--record-thr", type=float, default=8.0,
                    help="lowest network MAD written to the parquet; both tiers are filters on it")
    ap.add_argument("--fam-chunk", type=int, default=24)
    ap.add_argument("--cc-device", default="auto", choices=["auto", "gpu", "cpu"],
                    help="where the per-station cc series live between the correlation and the network sum. "
                         "cpu = modest GPU footprint (~2 GB) on a shared card; the float32 values are identical.")
    ap.add_argument("--gpu-free-gb", type=float, default=20.0,
                    help="auto: keep the cc series on the GPU only if at least this much is free")
    ap.add_argument("--gpu-need-gb", type=float, default=6.0,
                    help="wait at the start of each day until at least this much GPU memory is free "
                         "(the card is shared with other agents)")
    ap.add_argument("--gpu-wait-s", type=float, default=3600.0, help="give up waiting after this long")
    ap.add_argument("--min-sites", type=int, default=2)
    ap.add_argument("--stack-pre-s", type=float, default=2.0)
    ap.add_argument("--stack-len-s", type=float, default=16.0)
    ap.add_argument("--no-stacks", action="store_true")
    ap.add_argument("--no-hfr", action="store_true")
    ap.add_argument("--hf-sta", default="PGC")
    ap.add_argument("--doy-from", type=int, default=0)
    ap.add_argument("--doy-to", type=int, default=999)
    ap.add_argument("--tag", default="")
    ap.add_argument("--out-dir", default="")
    ap.add_argument("--cache-dirs", default="", help="extra 20 Hz cache roots searched BEFORE cache20hz/")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    OUT = a.out_dir or C.TRC
    os.makedirs(OUT, exist_ok=True)
    tag = a.tag or a.ym
    t_start = time.time()

    fz = C.frozen()
    L = float(fz["template_length_s"])
    T_ext = float(fz["T_ext"]); T_core = float(fz["T_core"])
    # PER-NETWORK-EPOCH tiers (coordinator 2026-09-17): the false-alarm rate at a fixed MAD is a
    # property of the network, so the two tiers are re-derived per epoch by the same FA-only rule
    # (scripts/c_tiers_by_epoch.py).  Falls back to the frozen pair when the epoch is not derived yet.
    TIERS = C.tiers_by_epoch()
    TEPOCH = C.epoch_of_month(a.ym)
    if TEPOCH in TIERS:
        T_core, T_ext = TIERS[TEPOCH]
        print(f"[tiers] {a.ym} is network epoch {TEPOCH}: T_core {T_core} / T_ext {T_ext} "
              f"(frozen anchor {fz['T_core']} / {fz['T_ext']})", flush=True)
    hf = fz.get("hf_veto", {})
    hf_adopted = bool(hf.get("adopted", False)); hf_cut = float(hf.get("cut", np.nan))
    HFC = C.hf_cuts()                      # results/a4_hf_veto_cuts.csv, per-epoch cuts (may be empty)
    if HFC:
        print(f"[hf] {len(HFC)} per-epoch veto cuts from results/a4_hf_veto_cuts.csv", flush=True)
    tpl_path, tpl_tag = C.resolve_templates(fz, a.templates)
    tpl, sp, F, wgt, mv, ampref = C.load_frozen_tables(fz, tpl_path, L)
    stas_tbl = C.station_table()
    off = C.placement(tpl, sp, F, mv, stas_tbl, L)
    STAS = C.order_stations({k[1] for k in tpl})
    ep = C.epochs(); man = C.manifest()
    # ---- per-(station, epoch) template sets -------------------------------------------------------
    # Base = the frozen set.  A station-epoch AT OR AFTER its post-2012 instrument boundary is not
    # covered by 2012 instruments: label 2012_provisional, unless Track A's A-2 agent has dropped
    # results/a4_templates_<STA>_<epoch>.npz in (then that file is used and labelled by its basename).
    # TSTA[(sta, epoch_start)] = {fam: (template (3,n), placement offset s)}; LBL[(sta, epoch)] = label.
    TSTA, LBL = {}, {}
    POLICY = C.epoch_policy()
    DROPPED = {}                      # (sta, epoch) -> True when the policy excludes the station
    if POLICY:
        print(f"[tpl] epoch policy: {len(POLICY)} rows from results/trackC/templates_by_epoch.csv", flush=True)

    def build_set(s, es, post):
        pol = POLICY.get((s, es))
        if pol and pol["drop"]:
            # Track A-2: LRIV must not carry 2012 templates past 2019-05-30 (stack cc goes NEGATIVE).
            # The station is EXCLUDED from the network for this epoch; the loss is recorded per
            # family-day in `excluded_sta` so the completeness model can absorb it.
            TSTA[(s, es)] = {}
            LBL[(s, es)] = "dropped"
            DROPPED[(s, es)] = True
            return
        f_ovr = C.epoch_template_file(s, es, POLICY)
        if f_ovr is not None:
            z = np.load(f_ovr, allow_pickle=True)
            t2 = {k: v for k, v in z["tpl"].item().items() if k[1] == s and L in v}
            sp2 = z["sp"].item()
            off2 = C.placement(t2, sp2, F, mv, stas_tbl, L)
            TSTA[(s, es)] = {k[0]: (t2[k][L], off2[k]) for k in t2 if k in off2}
            LBL[(s, es)] = os.path.basename(f_ovr).replace(".npz", "")
        else:
            TSTA[(s, es)] = {k[0]: (tpl[k][L], off[k]) for k in tpl if k[1] == s and k in off}
            LBL[(s, es)] = (pol["label"] if (pol and pol["label"])
                            else ("2012_provisional" if post else "2012_final"))

    for s in STAS:
        seen = set()
        for r in ep[ep.sta == s].itertuples():
            for d0 in (r.start, pd.Timestamp(C.POST_2012_BOUNDARY[s]) if s in C.POST_2012_BOUNDARY else None):
                if d0 is None or not (r.start <= d0 <= r.end):
                    continue
                es = C.eff_epoch_start(s, d0, ep)
                if es in seen:
                    continue
                seen.add(es)
                build_set(s, es, C.is_post_boundary(s, d0))
        if not seen:
            build_set(s, "base", False)
    print("[tpl] station-epoch template sets: "
          + "; ".join(f"{s}@{es}={LBL[(s, es)]}({len(TSTA[(s, es)])})" for (s, es) in sorted(TSTA)), flush=True)
    NETSTA = {s: (ep[ep.sta == s].netsta.iloc[0] if (ep.sta == s).any()
                  else str(stas_tbl.loc[s, "netsta"])) for s in STAS}
    SITE = C.site_map(STAS)
    share_ok = C.share_ok_matrix()
    CACHES = [d.rstrip("/") for d in a.cache_dirs.split(",") if d] + [C.CACHE]
    n = int(round(L * FS)); Lmax = n
    print(f"[cfg] {a.ym}  L={L}s  T_ext={T_ext} T_core={T_core}  record>={a.record_thr}  "
          f"templates={tpl_tag} ({len(tpl)} (fam,sta), {len({k[0] for k in tpl})} families)  "
          f"stations={len(STAS)}  hf_veto_adopted={hf_adopted} cut={hf_cut:.4f}", flush=True)

    # ---- PNSN tremor in the family box, and the regional EQ catalog ----
    lat0, lat1, lon0, lon1 = C.family_box(F)
    tr = pd.read_csv(C.PNSN, usecols=["time", "lat", "lon"])
    tr["time"] = pd.to_datetime(tr.time, errors="coerce")
    tr = tr.dropna(subset=["time"])
    tr = tr[tr.lat.between(lat0, lat1) & tr.lon.between(lon0, lon1)]
    t_trem = np.sort(tr.time.astype("datetime64[ns]").astype("int64").to_numpy() / 1e9)
    trem_day = tr.time.dt.strftime("%Y-%m-%d").value_counts().to_dict()
    eqp = f"{C.TRC}/eq_catalog.csv"
    if os.path.exists(eqp):
        eq = pd.read_csv(eqp)
        eq["ot"] = pd.to_datetime(eq.ot, errors="coerce")
        t_eq = np.sort(eq.ot.dropna().astype("datetime64[ns]").astype("int64").to_numpy() / 1e9)
    else:
        print(f"[eq] WARNING {eqp} missing -> eq_flag is NaN; run scripts/c_eq_catalog.py", flush=True)
        t_eq = np.zeros(0)

    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    torch.set_num_threads(4)
    CCDEV = [dev]
    if dev == "cuda":
        free_b, tot_b = torch.cuda.mem_get_info()
        if a.cc_device == "cpu" or (a.cc_device == "auto" and free_b / 2 ** 30 < a.gpu_free_gb):
            CCDEV[0] = "cpu"
        print(f"[env] {dev} {torch.cuda.get_device_name(0)}  free {free_b/2**30:.1f}/{tot_b/2**30:.1f} GB "
              f"-> cc series on {CCDEV[0].upper()}", flush=True)
    else:
        print("[env] cpu", flush=True)

    days = [d for d in C.month_days(a.ym) if a.doy_from <= d.dayofyear <= a.doy_to]
    note = {(r.netsta, r.date): (r.note if isinstance(r.note, str) else "")
            for r in man[man.date.str[:7] == a.ym].itertuples()}
    det_rows, fd_rows = [], []
    stack_keys, stack_arr = [], []
    nsp = int(round(a.stack_pre_s * FS)); nsl = int(round(a.stack_len_s * FS))
    n_unobs = 0
    per_day_log = []
    FCH = [a.fam_chunk]      # adaptive family chunk: only ever shrinks, after a CUDA OOM on a shared card

    def wait_gpu():
        """The H100 is shared: never start a day unless the working set fits, and never crash the month
        because another agent grabbed the card between two days."""
        if dev != "cuda":
            return
        t0w = time.time()
        said = False
        while True:
            torch.cuda.empty_cache()
            free_b = torch.cuda.mem_get_info()[0]
            if free_b / 2 ** 30 >= a.gpu_need_gb or time.time() - t0w > a.gpu_wait_s:
                if a.cc_device == "auto":
                    want = "cuda" if free_b / 2 ** 30 >= a.gpu_free_gb else "cpu"
                    if want != CCDEV[0]:
                        print(f"  [gpu] {free_b/2**30:.1f} GB free -> cc series on {want.upper()}", flush=True)
                        CCDEV[0] = want
                if said:
                    print(f"  [gpu] {free_b/2**30:.1f} GB free after {time.time()-t0w:.0f}s, continuing", flush=True)
                return
            if not said:
                print(f"  [gpu] only {free_b/2**30:.1f} GB free, waiting for {a.gpu_need_gb} GB ...", flush=True)
                said = True
            time.sleep(20)

    for day in days:
        doy = int(day.dayofyear); year = int(day.year); ds = day.strftime("%Y-%m-%d")
        t_day = time.time()
        wait_gpu()
        al = C.alive(ep, day)
        alive_sta = {r.sta: r.epoch_id for r in al.itertuples()}
        EK = {s: C.eff_epoch_start(s, day, ep) for s in STAS if s in alive_sta}
        TS = {}
        for s in EK:
            key = (s, EK[s])
            if key not in TSTA and (s, "base") not in TSTA:
                raise SystemExit(f"[tpl] no template set for {s} epoch {EK[s]} — station_epochs.csv and "
                                 f"the template-set build disagree; refusing to run a station silently blind")
            TS[s] = TSTA.get(key, TSTA.get((s, "base"), {}))
        # stations the epoch policy EXCLUDES today, and which had data (so the exclusion is a real loss)
        DROP_TODAY = {s for s in EK if DROPPED.get((s, EK[s]))
                      and any(os.path.exists(f"{c}/{NETSTA[s]}/{year}/{doy:03d}.npz") for c in CACHES)}
        D, T0, NPT, EPID = {}, {}, {}, {}
        for s in STAS:
            if s not in alive_sta or s in DROP_TODAY:
                continue
            p = next((q for q in (f"{c}/{NETSTA[s]}/{year}/{doy:03d}.npz" for c in CACHES)
                      if os.path.exists(q)), None)
            if p is None:
                continue
            try:
                z = np.load(p)
            except Exception as e:
                print(f"  [skip] {NETSTA[s]} {year}.{doy:03d}: {e}", flush=True)
                continue
            D[s] = torch.from_numpy(z["x"]).to(dev)
            T0[s] = float(z["t0"]); NPT[s] = D[s].shape[1]; EPID[s] = alive_sta[s]
        if not D:
            per_day_log.append(dict(date=ds, n_sta=0, n_fam=0, n_det=0, secs=0.0))
            continue
        net_epoch = "|".join(f"{s}:{EPID[s].split(':')[1]}" for s in sorted(D))
        Nmin = min(NPT.values())
        NF = int(2 ** np.ceil(np.log2(Nmin + Lmax)))
        DH, EN = {}, {}
        for s in D:
            x = D[s][:, :Nmin]
            DH[s] = torch.fft.rfft(x, n=NF)
            xd = x.double()
            EN[s] = torch.cumsum(torch.cat([torch.zeros((x.shape[0], 1), device=dev, dtype=torch.float64),
                                            xd * xd], 1), 1)
        NV = Nmin - n
        fam_list = sorted({f for s in D for f in TS.get(s, {})})
        S3, VALID, E3W = {}, {}, {}
        for s in D:
            fs_ = [f for f in fam_list if f in TS.get(s, {})]
            if not fs_:
                continue
            try:
                out3 = torch.empty((len(fs_), NV), device=CCDEV[0])
            except torch.OutOfMemoryError:
                torch.cuda.empty_cache(); CCDEV[0] = "cpu"
                print(f"  [oom] cc series moved to CPU for the rest of the run", flush=True)
                out3 = torch.empty((len(fs_), NV), device="cpu")
            ccdev = CCDEV[0]
            ew = (EN[s][:, n:Nmin + 1] - EN[s][:, :NV + 1])[:, :NV]
            efloor = 1e-10 * ew.median()
            e = ew.clamp_min(efloor).sqrt().float()
            ch = FCH[0]
            i0 = 0
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
                    if ch == 1:
                        raise
                    ch = max(1, ch // 2); FCH[0] = ch
                    print(f"  [oom] {s}: fam-chunk -> {ch} (kept for the rest of the run)", flush=True)
                    continue
                i0 += len(ff)
            S3[s] = (fs_, out3)
            e3 = ew.sum(0)
            VALID[s] = (e3 > 1e-6 * e3.median()) & (ew.min(0).values > 0)
            E3W[s] = e3.float()
            del ew, e3, e
            del DH[s], EN[s]
            if dev == "cuda":
                torch.cuda.empty_cache()
        DH.clear(); EN.clear()

        ntr = int(trem_day.get(ds, 0))
        day_det = []
        for f in sorted({f2 for s in S3 for f2 in TS.get(s, {})}):
            ss = [s for s in D if s in S3 and f in TS.get(s, {})]
            if not ss:
                continue
            groups = {}
            for s in ss:
                groups.setdefault(SITE[s], []).append(s)
            n_sites_f = len(groups)
            sites_f = ",".join(sorted(groups))
            excl_f = ",".join(sorted(s2 for s2 in DROP_TODAY if (float(f), s2) in off))
            if n_sites_f < a.min_sites:
                n_unobs += 1
                fd_rows.append(dict(fam=int(f), date=ds, doy=doy, observable=False, n_sites=n_sites_f,
                                    n_sta=len(ss), mad_scale=np.nan, med_scale=np.nan,
                                    n_tremor_day=ntr, n_valid_s=np.nan, n_det_rec=0, n_det_ext=0,
                                    n_det_core=0, epoch=net_epoch, sites=sites_f, excluded_sta=excl_f))
                continue
            o = {s: TS[s][f][1] for s in ss}
            tref0 = min(T0[s] for s in ss)
            k = {s: int(round((tref0 - T0[s] + o[s]) * FS)) for s in ss}
            b = max(0, -min(k.values()))
            k = {s: k[s] + b for s in ss}
            tref = tref0 + b / FS
            NVo = NV - max(k.values())
            if NVo <= 1000:
                fd_rows.append(dict(fam=int(f), date=ds, doy=doy, observable=False, n_sites=n_sites_f,
                                    n_sta=len(ss), mad_scale=np.nan, med_scale=np.nan,
                                    n_tremor_day=ntr, n_valid_s=np.nan, n_det_rec=0, n_det_ext=0,
                                    n_det_core=0, epoch=net_epoch, sites=sites_f, excluded_sta=excl_f))
                continue
            wm = {s: 3.0 * float(wgt.get((f, s), 1.0)) for s in ss}
            Wsite = {q: float(np.sqrt(sum(wm[s] ** 2 for s in mem))) for q, mem in groups.items()}  # mergeB
            acc_n = torch.zeros(NVo, device=dev)
            vmask = torch.ones(NVo, dtype=torch.bool, device=dev)
            weff, wtot, SEG = {}, 0.0, {}
            for q, mem in groups.items():
                wsum_m = sum(wm[s] for s in mem)
                cs = torch.zeros(NVo, device=dev)
                vm = torch.ones(NVo, dtype=torch.bool, device=dev)
                for s in mem:
                    fs_, cc3 = S3[s]
                    i = fs_.index(f)
                    seg = cc3[i, k[s]:k[s] + NVo]
                    if seg.device.type != dev:
                        seg = seg.to(dev, non_blocking=True)
                    SEG[s] = seg
                    cs += (wm[s] / wsum_m) * seg
                    vm &= VALID[s][k[s]:k[s] + NVo]
                    weff[s] = Wsite[q] * wm[s] / wsum_m
                acc_n += Wsite[q] * cs
                vmask &= vm
                wtot += Wsite[q]
            acc_n /= wtot
            am = acc_n[vmask]
            if am.numel() < 1000:
                fd_rows.append(dict(fam=int(f), date=ds, doy=doy, observable=False, n_sites=n_sites_f,
                                    n_sta=len(ss), mad_scale=np.nan, med_scale=np.nan,
                                    n_tremor_day=ntr, n_valid_s=float(am.numel() / FS), n_det_rec=0,
                                    n_det_ext=0, n_det_core=0, epoch=net_epoch, sites=sites_f, excluded_sta=excl_f))
                del acc_n, vmask
                continue
            medn = am.median(); madn = (am - medn).abs().median().clamp_min(1e-6)
            idx = decluster_idx(acc_n, medn + a.record_thr * madn, n, vmask)
            n_ext = n_core = 0
            if len(idx):
                ii = torch.from_numpy(idx).to(dev)
                vals = acc_n[ii]
                nmad = ((vals - medn) / madn).cpu().numpy()
                amp_acc = torch.zeros(len(ii), device=dev); wsum = 0.0
                con_pgc = torch.zeros(len(ii), device=dev)
                for s in ss:
                    w = weff[s]
                    js = ii + k[s]
                    ccv = SEG[s][ii]                       # = cc3[i, ii + k[s]] (seg starts at k[s])
                    amp_acc += w * ccv * torch.sqrt(E3W[s][js] / ampref[(s, L)])
                    if SITE[s] == "PGC":
                        con_pgc += w * ccv
                    wsum += w
                td = tref + idx / FS
                n_ext = int((nmad >= T_ext).sum()); n_core = int((nmad >= T_core).sum())
                labs = {LBL.get((s, EK[s]), "2012_final") for s in ss}
                if labs == {"2012_final"}:
                    tlab = "2012_final"
                elif labs <= {"2012_final", "2012_provisional"}:
                    tlab = "2012_provisional"
                else:
                    tlab = "per_epoch:" + ",".join(sorted(l for l in labs if l.startswith("a4")))
                day_det.append(pd.DataFrame(dict(
                    fam=int(f), date=ds, doy=doy, t=td, net_mad=nmad,
                    net_cc=vals.cpu().numpy(), rel_amp=(amp_acc / wsum).cpu().numpy(),
                    n_sites=n_sites_f, n_sta=len(ss),
                    con_pgc=(con_pgc / wsum).cpu().numpy(),
                    wsh_pgc=float(sum(weff[s] for s in ss if SITE[s] == "PGC") / wsum),
                    mad_scale=float(madn), med_scale=float(medn), epoch=net_epoch,
                    epoch_template=tlab,
                    n_sta_provisional=sum(LBL.get((s, EK[s]), "") == "2012_provisional" for s in ss),
                    n_sta_post_boundary=sum(C.is_post_boundary(s, day) for s in ss),
                    n_sta_excluded=len([x for x in excl_f.split(",") if x]))))
                # ---- per-station daily stacks of the T_core detections ----
                if not a.no_stacks and n_core:
                    sel = np.where(nmad >= T_core)[0]
                    for s in ss:
                        st_i0 = idx[sel] + k[s] - nsp
                        ok = (st_i0 >= 0) & (st_i0 + nsl <= Nmin)
                        if not ok.any():
                            continue
                        cols = torch.from_numpy(st_i0[ok]).to(dev)[:, None] + \
                            torch.arange(nsl, device=dev)[None]
                        segs = D[s][:, :Nmin][:, cols].permute(1, 0, 2).double()   # (K,3,nsl)
                        r = segs.pow(2).mean(dim=(1, 2)).sqrt()
                        good = torch.isfinite(r) & (r > 0)
                        if not bool(good.any()):
                            continue
                        segs = segs[good] / r[good][:, None, None]
                        stack_keys.append((int(f), s, doy, int(good.sum())))
                        stack_arr.append(segs.mean(0).float().cpu().numpy())
                        del segs
            fd_rows.append(dict(fam=int(f), date=ds, doy=doy, observable=True, n_sites=n_sites_f,
                                n_sta=len(ss), mad_scale=float(madn), med_scale=float(medn),
                                n_tremor_day=ntr, n_valid_s=float(am.numel() / FS),
                                n_det_rec=int(len(idx)), n_det_ext=n_ext, n_det_core=n_core,
                                epoch=net_epoch, sites=sites_f, excluded_sta=excl_f))
            del acc_n, vmask, am, SEG
        del S3, VALID, E3W, D
        if dev == "cuda":
            torch.cuda.empty_cache()
        nd = 0
        if day_det:
            Dd = pd.concat(day_det, ignore_index=True)
            nd = len(Dd)
            # ---- HF ratio at the PGC site, from the native-rate trace of THIS day ----
            if not a.no_hfr and a.hf_sta in NETSTA:
                hs = a.hf_sta
                tstart = Dd.t.to_numpy() + np.array([TS.get(hs, {}).get(float(fm), (None, np.nan))[1]
                                                     for fm in Dd.fam.to_numpy()])
                ok = np.isfinite(tstart)
                r = np.full(len(Dd), np.nan); hfs = np.nan
                if ok.any():
                    r[ok], hfs = hf_ratio_day(NETSTA[hs], year, doy, tstart[ok], L,
                                              note=note.get((NETSTA[hs], ds), ""))
                Dd["hfr"] = r
                Dd["hf_fs"] = hfs          # the HF station's NATIVE rate: the frozen cut 3.6384 was
                                           # calibrated on PGC BH 40 Hz, and the ratio is band-limited
                                           # by the anti-alias filter, so it is epoch-dependent
            else:
                Dd["hfr"] = np.nan; Dd["hf_fs"] = np.nan
            det_rows.append(Dd)
        per_day_log.append(dict(date=ds, n_sta=len(alive_sta), n_fam=len(fam_list), n_det=nd,
                                secs=time.time() - t_day))
        print(f"  [{ds}] stations {len(T0)} families {len(fam_list)} detections>={a.record_thr} {nd} "
              f"({time.time()-t_day:.1f}s)", flush=True)

    # ---------------- month outputs ----------------
    FD = pd.DataFrame(fd_rows)
    if len(FD):
        # quiet day = no PNSN tremor in the box that day AND >= 3 days from any tremor day (Merlin Step 6.3)
        td = pd.Series({d: c for d, c in trem_day.items()})
        tdays = pd.to_datetime(pd.Index(td[td > 0].index))
        dts = pd.to_datetime(FD.date.unique())
        near = {}
        for d in dts:
            near[d.strftime("%Y-%m-%d")] = bool(len(tdays) and (np.abs((tdays - d).days) <= 3).any())
        FD["quiet_flag"] = (FD.n_tremor_day == 0) & (~FD.date.map(near).astype(bool))
        FD.to_parquet(f"{OUT}/family_days_{tag}.parquet", index=False)
    DET = pd.concat(det_rows, ignore_index=True) if det_rows else pd.DataFrame(
        columns=["fam", "date", "doy", "t", "net_mad", "net_cc", "rel_amp", "n_sites", "n_sta",
                 "con_pgc", "wsh_pgc", "mad_scale", "med_scale", "epoch", "epoch_template",
                 "n_sta_provisional", "n_sta_post_boundary", "n_sta_excluded", "hfr", "hf_fs"])
    if len(DET):
        DET = DET.sort_values("t").reset_index(drop=True)
        DET["time"] = pd.to_datetime(DET.t, unit="s")
        DET["tier_ext"] = DET.net_mad >= T_ext
        DET["tier_core"] = DET.net_mad >= T_core
        DET["T_core_used"] = T_core; DET["T_ext_used"] = T_ext; DET["tier_epoch"] = TEPOCH
        if HFC:
            _hc = [C.hf_cut_for(HFC, a.hf_sta, d_, f_, hf_cut)
                   for d_, f_ in zip(DET.date.to_numpy(), DET.hf_fs.to_numpy())]
            DET["hf_cut"] = [x[0] for x in _hc]
            DET["hf_cut_src"] = [x[1] for x in _hc]
            for u in sorted(set(DET.hf_cut_src)):
                print(f"[hf] cut {DET.hf_cut[DET.hf_cut_src == u].iloc[0]:.4f} from {u}", flush=True)
        else:
            DET["hf_cut"] = hf_cut
            DET["hf_cut_src"] = "frozen_operating_point"
        DET["hf_pass"] = ~(np.isfinite(DET.hfr.to_numpy()) & (DET.hfr.to_numpy() > DET.hf_cut.to_numpy()))
        DET["hf_veto_adopted"] = hf_adopted
        # tremor / EQ flags
        tv = DET.t.to_numpy()
        DET["tremor_flag"] = (np.abs(tv - t_trem[np.clip(np.searchsorted(t_trem, tv), 0, len(t_trem) - 1)]) <= 7200) | \
                             (np.abs(tv - t_trem[np.clip(np.searchsorted(t_trem, tv) - 1, 0, len(t_trem) - 1)]) <= 7200) \
                             if len(t_trem) else False
        DET["n_tremor_day"] = DET.date.map(lambda d: int(trem_day.get(d, 0)))
        if len(t_eq):
            j = np.searchsorted(t_eq, tv)
            dlo = np.abs(tv - t_eq[np.clip(j - 1, 0, len(t_eq) - 1)])
            dhi = np.abs(tv - t_eq[np.clip(j, 0, len(t_eq) - 1)])
            DET["eq_dt"] = np.minimum(dlo, dhi)
            DET["eq_flag"] = DET.eq_dt <= 60.0
        else:
            DET["eq_dt"] = np.nan; DET["eq_flag"] = False
        # ---- FROZEN cross-family assignment, evaluated at each tier ----
        fam_i = DET.fam.to_numpy().astype(int)
        for lbl, T in (("ext", T_ext), ("core", T_core)):
            sel = (DET.net_mad >= T).to_numpy()
            kc = np.zeros(len(DET), bool); km = np.zeros(len(DET), bool)
            if sel.any():
                kc[sel] = C.local_max_keep_cc(DET.t.to_numpy()[sel], DET.net_cc.to_numpy()[sel].astype(float),
                                              fam_i[sel], share_ok, 2.0)
                km[sel] = C.local_max_keep(DET.t.to_numpy()[sel], DET.net_mad.to_numpy()[sel], 2.0)
            DET[f"keep_xfam_{lbl}"] = kc
            DET[f"keep_mad_{lbl}"] = km
        # quiet-day flag per detection (Merlin Step 6.3 / Track A step 3b: the floor lambda is computed
        # per family per season from these, and ~22 % of quiet-day T_core detections are real inter-ETS
        # repeats, so quiet != false)
        if len(FD):
            qmap = {(int(r.fam), r.date): bool(r.quiet_flag) for r in FD.itertuples()}
            DET["quiet_flag"] = [qmap.get((int(f_), d_), False)
                                 for f_, d_ in zip(DET.fam.to_numpy(), DET.date.to_numpy())]
        else:
            DET["quiet_flag"] = False
        # THE CATALOG COLUMNS: tier AND the frozen HF veto AND the frozen cross-family assignment.
        # Vetoed / suppressed rows are KEPT (every column above is preserved) so the false-alarm floor
        # can be re-derived; they are simply not members of the catalog.
        DET["cat_ext"] = DET.tier_ext & DET.hf_pass & DET.keep_xfam_ext
        DET["cat_core"] = DET.tier_core & DET.hf_pass & DET.keep_xfam_core
        DET.to_parquet(f"{OUT}/detections_{tag}.parquet", index=False)
    if not a.no_stacks and stack_keys:
        K = np.array([[k[0], 0, k[2], k[3]] for k in stack_keys], np.int64)
        np.savez_compressed(f"{OUT}/stacks_{tag}.npz",
                            fam=K[:, 0], sta=np.array([k[1] for k in stack_keys]),
                            doy=K[:, 2], n=K[:, 3], stack=np.stack(stack_arr).astype(np.float32),
                            pre_s=np.float64(a.stack_pre_s), len_s=np.float64(a.stack_len_s),
                            fs=np.float64(FS), tier=np.str_("T_core"))
    meta = dict(ym=a.ym, tag=tag, templates=tpl_path, template_set=tpl_tag, L=L,
                T_ext=T_ext, T_core=T_core, tier_epoch=TEPOCH,
                T_frozen=(float(fz["T_core"]), float(fz["T_ext"])),
                record_thr=a.record_thr, min_sites=a.min_sites,
                hf_veto_adopted=hf_adopted, hf_cut=hf_cut, hf_station=a.hf_sta,
                hf_cuts_per_epoch=len(HFC),
                epoch_policy_rows=len(POLICY),
                n_days=len(days), n_detections=int(len(DET)),
                n_det_ext=int(DET.tier_ext.sum()) if len(DET) else 0,
                n_det_core=int(DET.tier_core.sum()) if len(DET) else 0,
                n_det_core_hfpass=int((DET.tier_core & DET.hf_pass).sum()) if len(DET) else 0,
                n_cat_ext=int(DET.cat_ext.sum()) if len(DET) else 0,
                n_cat_core=int(DET.cat_core.sum()) if len(DET) else 0,
                epoch_template_counts=(DET.epoch_template.value_counts().to_dict() if len(DET) else {}),
                n_family_days=int(len(FD)),
                n_family_days_observable=int(FD.observable.sum()) if len(FD) else 0,
                n_family_days_unobservable=int((~FD.observable).sum()) if len(FD) else 0,
                n_stacks=len(stack_keys), runtime_s=time.time() - t_start,
                per_day=per_day_log)
    with open(f"{OUT}/run_{tag}.json", "w") as fh:
        json.dump(meta, fh, indent=1, default=str)
    print(f"[done] {a.ym}: {len(DET)} detections >= {a.record_thr} MAD "
          f"({meta['n_det_ext']} T_ext, {meta['n_det_core']} T_core), "
          f"{meta['n_family_days_observable']}/{len(FD)} family-days observable, "
          f"{len(stack_keys)} stacks, {meta['runtime_s']:.0f}s -> {OUT}/detections_{tag}.parquet")


if __name__ == "__main__":
    main()
