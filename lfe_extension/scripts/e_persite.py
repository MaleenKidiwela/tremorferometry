#!/usr/bin/env python
"""Track E step 1 — PER-SITE recomputation of the frozen matched filter at stored detection times.

The production run (scripts/c_run_month.py) stores only the NETWORK score of each detection plus the
PGC site's contribution.  Track E needs, for every detection, what each INDEPENDENT SITE saw on its
own: the site's weighted-mean cc at the detection time in units of THAT site's own daily MAD, and the
S-window energy excess at the predicted moveout against the same station's own distribution that day.

This script is a copy of c_run_month.py's measurement core (same cache, same templates, same
placement, same mergeB site merge, same float32 path) with the network sum replaced by the per-site
series and with the family loop restricted to the (family, day) pairs that actually carry a stored
detection.  `c_run_month.py` is NEVER modified or imported for anything but its shared config module.

Self-check: the weighted sum of the per-site series at the stored detection index must reproduce the
stored `net_cc` and `net_mad`; the max deviation over the month is written to the meta json.

Modes
  --mode prod --ym YYYY-MM     detections from results/trackC/detections_<ym>.parquet, production
                               station/epoch policy from c_common.
  --mode ref  --year 2012 ...  detections from a Track A frozen-run parquet (--det), fixed station
                               list (--stations), no epoch policy.  Also generates the SITE-ROLLED
                               chance control (one distinct shift per site, 300-3000 s) and writes
                               its detections with the same per-site features -> the negative class.

Outputs (results/trackE/):
  persite_<tag>.parquet   long format, one row per (detection, site)
  persite_<tag>_det.parquet  one row per detection (control detections included, ctrl=1)
  persite_<tag>_meta.json
"""
from __future__ import annotations
import os, sys, json, time, argparse, warnings
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import c_common as C

FS = C.FS
OUT = f"{C.RES}/trackE"


def decluster_idx(x, thr, n, valid):
    import torch, torch.nn.functional as Fnn
    mx = Fnn.max_pool1d(x[None, None], kernel_size=2 * (n // 2) + 1, stride=1, padding=n // 2)[0, 0]
    return torch.nonzero((x >= mx - 1e-12) & valid & (x >= thr)).squeeze(1).cpu().numpy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="prod", choices=["prod", "ref"])
    ap.add_argument("--ym", default="")
    ap.add_argument("--year", type=int, default=0)
    ap.add_argument("--doy-from", type=int, default=0)
    ap.add_argument("--doy-to", type=int, default=999)
    ap.add_argument("--det", default="", help="ref mode: the Track A detections parquet")
    ap.add_argument("--stations", default="", help="ref mode: station list")
    ap.add_argument("--netsta", default="", help="ref mode: STA=NET.STA map")
    ap.add_argument("--cache-dirs", default="", help="extra 20 Hz cache roots searched BEFORE cache20hz/")
    ap.add_argument("--tag", default="")
    ap.add_argument("--record-thr", type=float, default=8.0)
    ap.add_argument("--control", action="store_true", help="also generate site-rolled control detections")
    ap.add_argument("--fam-chunk", type=int, default=24)
    ap.add_argument("--gpu-need-gb", type=float, default=6.0)
    ap.add_argument("--e-win-s", type=float, default=4.0, help="S / pre-arrival energy window")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    os.makedirs(OUT, exist_ok=True)
    tag = a.tag or (a.ym if a.mode == "prod" else f"ref{a.year}")
    t_start = time.time()

    fz = C.frozen()
    L = float(fz["template_length_s"])
    tpl_path, tpl_tag = C.resolve_templates(fz, "")
    tpl, sp, F, wgt, mv, ampref = C.load_frozen_tables(fz, tpl_path, L)
    stas_tbl = C.station_table()
    off = C.placement(tpl, sp, F, mv, stas_tbl, L)
    n = int(round(L * FS))
    nE = int(round(a.e_win_s * FS))

    # ---------------- detections to score ----------------
    if a.mode == "prod":
        dp = f"{C.TRC}/detections_{a.ym}.parquet"
        DET = pd.read_parquet(dp, columns=["fam", "date", "doy", "t", "net_mad", "net_cc",
                                           "mad_scale", "med_scale", "n_sites", "n_sta"])
        days = [d for d in C.month_days(a.ym)]
        year_of = {d.strftime("%Y-%m-%d"): int(d.year) for d in days}
        ep = C.epochs()
        STAS = C.order_stations({k[1] for k in tpl})
        NETSTA = {s: (ep[ep.sta == s].netsta.iloc[0] if (ep.sta == s).any()
                      else str(stas_tbl.loc[s, "netsta"])) for s in STAS}
        POLICY = C.epoch_policy()
    else:
        DET = pd.read_parquet(a.det)
        DET = DET[(DET.L == L) if "L" in DET.columns else slice(None)]
        DET = DET[(DET.doy >= a.doy_from) & (DET.doy <= a.doy_to)].copy()
        DET["date"] = [pd.Timestamp(f"{a.year}-01-01") + pd.Timedelta(days=int(d) - 1) for d in DET.doy]
        DET["date"] = DET.date.dt.strftime("%Y-%m-%d")
        DET = DET.rename(columns={"medn": "med_scale", "madn": "mad_scale"})
        days = sorted(pd.to_datetime(DET.date.unique()))
        year_of = {d.strftime("%Y-%m-%d"): a.year for d in days}
        STAS = [s for s in a.stations.split(",") if s]
        NETSTA = dict(kv.split("=") for kv in a.netsta.split(",") if kv)
        ep, POLICY = None, {}
    DET = DET[DET.net_mad >= a.record_thr].reset_index(drop=True)
    SITE = C.site_map(STAS)
    CACHES = [d.rstrip("/") for d in a.cache_dirs.split(",") if d] + [C.CACHE]
    print(f"[cfg] {tag} mode={a.mode} detections={len(DET)} days={len(days)} stations={len(STAS)} "
          f"templates={tpl_tag}", flush=True)
    if not len(DET):
        pd.DataFrame().to_parquet(f"{OUT}/persite_{tag}.parquet")
        json.dump(dict(tag=tag, n_det=0), open(f"{OUT}/persite_{tag}_meta.json", "w"))
        return

    # per-(station, epoch) template sets — production only; ref mode uses the frozen set everywhere
    TSTA, LBL, DROPPED = {}, {}, {}
    if a.mode == "prod":
        def build_set(s, es):
            pol = POLICY.get((s, es))
            if pol and pol["drop"]:
                TSTA[(s, es)] = {}; DROPPED[(s, es)] = True; return
            f_ovr = C.epoch_template_file(s, es, POLICY)
            if f_ovr is not None:
                z = np.load(f_ovr, allow_pickle=True)
                t2 = {k: v for k, v in z["tpl"].item().items() if k[1] == s and L in v}
                off2 = C.placement(t2, z["sp"].item(), F, mv, stas_tbl, L)
                TSTA[(s, es)] = {k[0]: (t2[k][L], off2[k]) for k in t2 if k in off2}
            else:
                TSTA[(s, es)] = {k[0]: (tpl[k][L], off[k]) for k in tpl if k[1] == s and k in off}
        for s in STAS:
            seen = set()
            for r in ep[ep.sta == s].itertuples():
                for d0 in (r.start, pd.Timestamp(C.POST_2012_BOUNDARY[s])
                           if s in C.POST_2012_BOUNDARY else None):
                    if d0 is None or not (r.start <= d0 <= r.end):
                        continue
                    es = C.eff_epoch_start(s, d0, ep)
                    if es in seen:
                        continue
                    seen.add(es); build_set(s, es)
            if not seen:
                build_set(s, "base")
    else:
        for s in STAS:
            TSTA[(s, "base")] = {k[0]: (tpl[k][L], off[k]) for k in tpl if k[1] == s and k in off}

    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    torch.set_num_threads(4)
    if dev == "cuda":
        while torch.cuda.mem_get_info()[0] / 2 ** 30 < a.gpu_need_gb:
            torch.cuda.empty_cache(); time.sleep(20)

    rows_site, rows_det = [], []
    chk = []
    rng_day = np.random.default_rng(97)
    FCH = [a.fam_chunk]

    for day in days:
        ds = pd.Timestamp(day).strftime("%Y-%m-%d")
        Dd = DET[DET.date == ds]
        if not len(Dd) and not a.control:
            continue
        doy = int(pd.Timestamp(day).dayofyear); year = year_of[ds]
        if a.mode == "prod":
            al = C.alive(ep, pd.Timestamp(day))
            alive_sta = {r.sta: r.epoch_id for r in al.itertuples()}
            EK = {s: C.eff_epoch_start(s, pd.Timestamp(day), ep) for s in STAS if s in alive_sta}
            TS = {s: TSTA.get((s, EK[s]), TSTA.get((s, "base"), {})) for s in EK}
            DROP_TODAY = {s for s in EK if DROPPED.get((s, EK[s]))}
        else:
            alive_sta = {s: "ref" for s in STAS}
            TS = {s: TSTA[(s, "base")] for s in STAS}
            DROP_TODAY = set()
        D, T0, NPT = {}, {}, {}
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
                print(f"  [skip] {NETSTA[s]} {year}.{doy:03d}: {e}", flush=True); continue
            D[s] = torch.from_numpy(z["x"]).to(dev); T0[s] = float(z["t0"]); NPT[s] = D[s].shape[1]
        if not D:
            continue
        Nmin = min(NPT.values()); NF = int(2 ** np.ceil(np.log2(Nmin + n)))
        DH, EN = {}, {}
        for s in D:
            x = D[s][:, :Nmin]
            DH[s] = torch.fft.rfft(x, n=NF)
            xd = x.double()
            EN[s] = torch.cumsum(torch.cat([torch.zeros((x.shape[0], 1), device=dev,
                                                        dtype=torch.float64), xd * xd], 1), 1)
        NV = Nmin - n
        fams_day = sorted(set(Dd.fam.astype(int).tolist())) if len(Dd) else []
        if a.control:
            fams_day = sorted({int(f) for s in D for f in TS.get(s, {})})
        if not fams_day:
            del DH, EN, D
            if dev == "cuda":
                torch.cuda.empty_cache()
            continue
        # ---- cc series + energy-ratio series, only for the families needed today ----
        S3, VALID, LRz = {}, {}, {}
        for s in D:
            fs_ = [float(f) for f in fams_day if float(f) in TS.get(s, {})]
            if not fs_:
                continue
            out3 = torch.empty((len(fs_), NV), device=dev)
            ew = (EN[s][:, n:Nmin + 1] - EN[s][:, :NV + 1])[:, :NV]
            efloor = 1e-10 * ew.median()
            e = ew.clamp_min(efloor).sqrt().float()
            ch = FCH[0]; i0 = 0
            while i0 < len(fs_):
                ff = fs_[i0:i0 + ch]
                try:
                    W = torch.from_numpy(np.stack([TS[s][f][0] for f in ff])).to(dev)
                    WH = torch.fft.rfft(W, n=NF)
                    num = torch.fft.irfft(WH.conj() * DH[s][None], n=NF)[:, :, :NV]
                    out3[i0:i0 + len(ff)] = (num / e[None]).mean(1)
                    del W, WH, num
                except torch.OutOfMemoryError:
                    torch.cuda.empty_cache()
                    if ch == 1:
                        raise
                    ch = max(1, ch // 2); FCH[0] = ch; continue
                i0 += len(ff)
            S3[s] = (fs_, out3)
            e3 = ew.sum(0)
            VALID[s] = (e3 > 1e-6 * e3.median()) & (ew.min(0).values > 0)
            # ---- S-window / pre-arrival energy ratio at every template-start index j ----
            cs = EN[s].sum(0)                                  # (Nmin+1,) cumulative 3-C energy
            jmax = min(NV, Nmin - nE)
            eS = (cs[nE:jmax + nE] - cs[:jmax])
            ePre = torch.full((jmax,), float("nan"), device=dev, dtype=torch.float64)
            ePre[nE:] = cs[nE:jmax] - cs[:jmax - nE]
            r = torch.log10(eS.clamp_min(1e-30) / ePre.clamp_min(1e-30))
            g = torch.isfinite(r) & VALID[s][:jmax]
            if int(g.sum()) > 100:
                md = r[g].median(); sd = (r[g] - md).abs().median().clamp_min(1e-6) * 1.4826
                zz = torch.full((NV,), float("nan"), device=dev, dtype=torch.float64)
                zz[:jmax] = (r - md) / sd
                zz[~torch.isfinite(zz)] = float("nan")
            else:
                zz = torch.full((NV,), float("nan"), device=dev, dtype=torch.float64)
            LRz[s] = zz.float()
            del ew, e, e3, cs, eS, ePre, r, DH[s], EN[s]
            if dev == "cuda":
                torch.cuda.empty_cache()
        DH.clear(); EN.clear()

        # one distinct site roll per day, 300-3000 s
        shifts = {}
        if a.control:
            qs = sorted({SITE[s] for s in D})
            sh = rng_day.choice(np.arange(int(300 * FS), int(3000 * FS)), size=len(qs), replace=False)
            shifts = {q: int(v) for q, v in zip(qs, sh)}

        for f in fams_day:
            ff = float(f)
            ss = [s for s in D if s in S3 and ff in TS.get(s, {})]
            if not ss:
                continue
            groups = {}
            for s in ss:
                groups.setdefault(SITE[s], []).append(s)
            o = {s: TS[s][ff][1] for s in ss}
            tref0 = min(T0[s] for s in ss)
            k = {s: int(round((tref0 - T0[s] + o[s]) * FS)) for s in ss}
            b = max(0, -min(k.values()))
            k = {s: k[s] + b for s in ss}
            tref = tref0 + b / FS
            NVo = NV - max(k.values())
            if NVo <= 1000:
                continue
            wm = {s: 3.0 * float(wgt.get((ff, s), 1.0)) for s in ss}
            Wsite = {q: float(np.sqrt(sum(wm[s] ** 2 for s in mem))) for q, mem in groups.items()}
            QS = sorted(groups)
            nq = len(QS)
            SC = torch.zeros((nq, NVo), device=dev)          # per-site weighted-mean cc
            SZ = torch.zeros((nq, NVo), device=dev)          # per-site energy-excess z
            vmask = torch.ones(NVo, dtype=torch.bool, device=dev)
            VMq = torch.ones((nq, NVo), dtype=torch.bool, device=dev)
            for qi, q in enumerate(QS):
                mem = groups[q]
                wsum_m = sum(wm[s] for s in mem)
                for s in mem:
                    fs_, cc3 = S3[s]
                    i = fs_.index(ff)
                    SC[qi] += (wm[s] / wsum_m) * cc3[i, k[s]:k[s] + NVo]
                    SZ[qi] += (wm[s] / wsum_m) * LRz[s][k[s]:k[s] + NVo]
                    VMq[qi] &= VALID[s][k[s]:k[s] + NVo]
                vmask &= VMq[qi]
            Wq = torch.tensor([Wsite[q] for q in QS], device=dev, dtype=torch.float32)
            wtot = float(Wq.sum())
            acc = (Wq[:, None] * SC).sum(0) / wtot
            nval = int(vmask.sum())
            if nval < 1000:
                del SC, SZ, VMq, acc, vmask
                continue
            V = SC[:, vmask]                                   # (nq, nval)
            med_q = V.median(dim=1).values
            mad_q = (V - med_q[:, None]).abs().median(dim=1).values.clamp_min(1e-9)
            del V
            jobs = []
            dsub = Dd[Dd.fam.astype(int) == int(f)]
            if len(dsub):
                idx = np.rint((dsub.t.to_numpy() - tref) * FS).astype(np.int64)
                ok = (idx >= 0) & (idx < NVo)
                if ok.any():
                    jobs.append((0, idx[ok], dsub[ok].reset_index(drop=True), None))
            if a.control:
                shq = torch.tensor([shifts[q] % NVo for q in QS], device=dev)
                SCr = torch.stack([torch.roll(SC[qi], int(shq[qi])) for qi in range(nq)])
                VMr = torch.stack([torch.roll(VMq[qi], int(shq[qi])) for qi in range(nq)])
                accs = (Wq[:, None] * SCr).sum(0) / wtot
                vms = VMr.all(0)
                del SCr, VMr
                if int(vms.sum()) >= 1000:
                    ams = accs[vms]
                    med_s = ams.median(); mad_s = (ams - med_s).abs().median().clamp_min(1e-6)
                    ic = decluster_idx(accs, med_s + a.record_thr * mad_s, n, vms)
                    if len(ic):
                        vals = accs[torch.from_numpy(ic).to(dev)]
                        nm = ((vals - med_s) / mad_s).cpu().numpy()
                        cdf = pd.DataFrame(dict(fam=int(f), date=ds, doy=doy,
                                                t=tref + ic / FS, net_mad=nm,
                                                net_cc=vals.cpu().numpy(),
                                                mad_scale=float(mad_s), med_scale=float(med_s),
                                                n_sites=nq, n_sta=len(ss)))
                        jobs.append((1, ic, cdf, [int(shq[qi]) for qi in range(nq)]))
                del accs, vms
            for ctrl, ii, dsub2, shf in jobs:
                it = torch.from_numpy(np.asarray(ii)).to(dev)
                if ctrl:
                    J = (it[None, :] - torch.tensor(shf, device=dev)[:, None]) % NVo
                    ccq = torch.gather(SC, 1, J)
                    zeq = torch.gather(SZ, 1, J)
                else:
                    ccq = SC[:, it]
                    zeq = SZ[:, it]
                zq = (ccq - med_q[:, None]) / mad_q[:, None]
                netv = (Wq[:, None] * ccq).sum(0) / wtot
                blk = torch.stack([ccq, zq, zeq]).cpu().numpy()      # ONE sync per job
                nc = netv.cpu().numpy()
                nd_ = blk.shape[2]
                rows_site.append(dict(fam=int(f), date=ds, t=np.repeat(dsub2.t.to_numpy(), 1),
                                      ctrl=ctrl, sites=QS, nq=nq,
                                      nsta=[len(groups[q]) for q in QS],
                                      wsh=[Wsite[q] / wtot for q in QS], blk=blk))
                d2 = dsub2.copy()
                d2["ctrl"] = ctrl
                d2["net_cc_chk"] = nc
                d2["n_sites_obs"] = nq
                rows_det.append(d2)
                if not ctrl and "net_cc" in dsub2.columns:
                    chk.append(float(np.abs(nc - dsub2.net_cc.to_numpy()).max()))
            del SC, SZ, VMq, acc, vmask
        del S3, VALID, LRz, D
        if dev == "cuda":
            torch.cuda.empty_cache()

    if rows_site:
        cols = dict(fam=[], date=[], t=[], ctrl=[], site=[], nsta_site=[], wsh=[],
                    cc_site=[], z_site=[], ze_site=[])
        for r in rows_site:
            nd_ = r["blk"].shape[2]
            for qi, q in enumerate(r["sites"]):
                cols["fam"].append(np.full(nd_, r["fam"], np.int32))
                cols["date"].append(np.full(nd_, r["date"], object))
                cols["t"].append(r["t"])
                cols["ctrl"].append(np.full(nd_, r["ctrl"], np.int8))
                cols["site"].append(np.full(nd_, q, object))
                cols["nsta_site"].append(np.full(nd_, r["nsta"][qi], np.int8))
                cols["wsh"].append(np.full(nd_, r["wsh"][qi], np.float32))
                cols["cc_site"].append(r["blk"][0, qi])
                cols["z_site"].append(r["blk"][1, qi])
                cols["ze_site"].append(r["blk"][2, qi])
        PS = pd.DataFrame({k: np.concatenate(v) for k, v in cols.items()})
    else:
        PS = pd.DataFrame()
    PD = pd.concat(rows_det, ignore_index=True) if rows_det else pd.DataFrame()
    PS.to_parquet(f"{OUT}/persite_{tag}.parquet", index=False)
    PD.to_parquet(f"{OUT}/persite_{tag}_det.parquet", index=False)
    meta = dict(tag=tag, mode=a.mode, n_det_in=int(len(DET)), n_det_out=int((PD.ctrl == 0).sum()) if len(PD) else 0,
                n_ctrl=int((PD.ctrl == 1).sum()) if len(PD) else 0, n_site_rows=int(len(PS)),
                max_abs_net_cc_dev=float(max(chk)) if chk else None,
                record_thr=a.record_thr, e_win_s=a.e_win_s, runtime_s=time.time() - t_start)
    json.dump(meta, open(f"{OUT}/persite_{tag}_meta.json", "w"), indent=1)
    print(f"[done] {tag} {meta}", flush=True)


if __name__ == "__main__":
    main()
