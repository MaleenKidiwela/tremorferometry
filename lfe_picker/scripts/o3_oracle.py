#!/usr/bin/env python
"""P0-O item 3 -- the PER-STATION HEADROOM ORACLE on the 2012 development year.

Question it answers, before any GPU-day is spent on a learned model: how much detection headroom is
there in a SOURCE-AGNOSTIC single-station statistic, relative to the template cc the frozen detector
already uses?  Two statistics are computed on exactly the same station-day series, the same cache,
the same placement and the same site merge:

  (i)  ENERGY  -- the E4-style incoherent 3-C S-window energy excess:
       log10 E(S-1 .. S+5 s) / E(the same 6-s window 20 s earlier), 1-8 Hz, 20 Hz cache,
       MAD-normalised over the valid samples of that station-day;
  (ii) CC      -- the 3-C normalised correlation with that family's 2012 template (the frozen
       measurement core of c_run_month / e_persite), MAD-normalised over the valid samples of that
       (family, station, day).

Both are read at the SAME predicted template start (origin + ttP + (S-P) - 0.5 s + moveout), merged
within the PGC cluster {PGC,B009,B010,B011} and the LRIV cluster {LRIV,B005,B006,B007} by the mean of
their MAD units, and SUMMED over the 9 sites -- "sum of MAD units at the predicted S per site".

Positives   Bostock's 2012 TEST-half, non-guard events (Track A's split, e5r_gt_bostock2012.parquet).
Negatives   (a) quiet-day windows, 2012 doy 228-242 (0 Bostock events, no PNSN tremor in the box);
            (b) ETS-day windows >= 120 s from ANY Bostock event (any family).
Rates       network-wide quiet-day FALSE EVENTS/day = declustered peaks per template-day x the median
            number of observable families per day; recall is read at 12/day and 116/day for BOTH
            statistics, i.e. at matched false-alarm rate, never at a matched threshold.

Decision (PRE-REGISTERED, results/merlin_review_plan_2026-09-17.md):
  energy-network recall < 1.25 x template  -> P2 is DISCOVERY-ONLY
  energy-network recall >= 1.50 x template -> a sensitivity claim is allowed
  in between                               -> report both readings.
"""
from __future__ import annotations
import os, sys, json, time, argparse, warnings
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import o_common as C

FS = C.FS
HERE = os.path.dirname(os.path.abspath(__file__))
PICK = os.path.dirname(HERE)
RES = f"{PICK}/results"
IN = C.IN
QUIET = set(range(228, 243))
S15 = ["PGC", "VGZ", "NLLB", "SNB", "LRIV", "B001", "B003", "B004", "B005", "B006", "B007",
       "B009", "B010", "B011", "B926"]
NETSTA = {"PGC": "CN.PGC", "VGZ": "CN.VGZ", "NLLB": "CN.NLLB", "SNB": "CN.SNB", "LRIV": "UW.LRIV",
          "B001": "PB.B001", "B003": "PB.B003", "B004": "PB.B004", "B005": "PB.B005",
          "B006": "PB.B006", "B007": "PB.B007", "B009": "PB.B009", "B010": "PB.B010",
          "B011": "PB.B011", "B926": "PB.B926"}


def mad_norm(x, valid):
    """(x - median) / (1.4826 * MAD) over the valid samples.  torch, float32."""
    import torch
    v = x[valid]
    if v.numel() < 1000:
        return None
    md = v.median()
    sd = (v - md).abs().median() * 1.4826
    if not torch.isfinite(sd) or sd <= 0:
        return None
    return (x - md) / sd


def decluster_peaks(x, valid, n):
    """values of the declustered local maxima (max over a window of n samples) among valid samples."""
    import torch, torch.nn.functional as Fnn
    mx = Fnn.max_pool1d(x[None, None], kernel_size=2 * (n // 2) + 1, stride=1, padding=n // 2)[0, 0]
    m = (x >= mx - 1e-12) & valid
    return x[m]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--year", type=int, default=2012)
    ap.add_argument("--doy-from", type=int, default=0)
    ap.add_argument("--doy-to", type=int, default=999)
    ap.add_argument("--e-win-s", type=float, default=6.0, help="S-1 .. S+5 s")
    ap.add_argument("--e-pre-s", type=float, default=20.0, help="reference window this much earlier")
    ap.add_argument("--n-neg-per-famday", type=int, default=40)
    ap.add_argument("--fam-chunk", type=int, default=24)
    ap.add_argument("--tag", default="oracle2012")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    t_start = time.time()
    import torch
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    torch.set_num_threads(8)

    fz = C.frozen()
    L = float(fz["template_length_s"])
    n = int(round(L * FS))
    tpl_path, tpl_tag = C.resolve_templates(fz, "")
    tpl, sp, F, wgt, mv, ampref = C.load_frozen_tables(fz, tpl_path, L)
    stas_tbl = C.station_table()
    off = C.placement(tpl, sp, F, mv, stas_tbl, L)
    TS = {s: {f: (tpl[(f, s)][L], off[(f, s)]) for (f, s2) in tpl if s2 == s and (f, s) in off}
          for s in S15}
    SITE = C.site_map(S15)
    SITES = sorted(set(SITE.values()))
    print(f"[cfg] templates={tpl_tag} sites={SITES}", flush=True)

    gt = pd.read_parquet(f"{IN}/e5r_gt_bostock2012.parquet").reset_index(drop=True)
    gt["doy"] = pd.DatetimeIndex(gt.time).dayofyear
    per_day = gt.groupby("doy").size()
    ETS = set(per_day[per_day >= 100].index)
    POS = gt[(gt.half == "test") & (~gt.guard)].copy().reset_index(drop=True)
    t_all = np.sort(gt.t_origin.to_numpy())
    print(f"[gt] {len(gt)} events, {len(POS)} test-half non-guard positives, "
          f"{len(ETS)} ETS days, {len(QUIET)} quiet days", flush=True)

    nE = int(round(a.e_win_s * FS))
    nPre = int(round(a.e_pre_s * FS))
    HALF = int(round(0.5 * FS))          # template start = S - 0.5 s

    CACHES = [C.REP_CACHE, C.EXT_CACHE, C.CACHE]
    days = sorted({int(d) for d in gt.doy.unique()} | QUIET)
    days = [d for d in days if a.doy_from <= d <= a.doy_to]

    rng = np.random.default_rng(11)
    pos_rows, neg_rows = [], []          # per (event/negative, station) single-station values
    net_pos, net_neg = [], []            # per event / negative, the 2 network statistics
    fa_peaks = {"E": [], "C": []}        # declustered NETWORK peak values on quiet family-days
    n_quiet_famday = 0
    sta_peaks = {"E": {s: [] for s in S15}, "C": {s: [] for s in S15}}
    n_quiet_staday = {s: 0 for s in S15}
    n_quiet_tpl_staday = {s: 0 for s in S15}
    nfam_obs = []

    for doy in days:
        t_day = time.time()
        D, T0, NPT = {}, {}, {}
        for s in S15:
            p = next((q for q in (f"{c}/{NETSTA[s]}/{a.year}/{doy:03d}.npz" for c in CACHES)
                      if os.path.exists(q)), None)
            if p is None:
                continue
            z = np.load(p)
            D[s] = torch.from_numpy(z["x"]).to(dev); T0[s] = float(z["t0"]); NPT[s] = D[s].shape[1]
        if not D:
            continue
        Nmin = min(NPT.values()); NF = int(2 ** np.ceil(np.log2(Nmin + n)))
        ZE, VAL, S3 = {}, {}, {}
        for s in list(D):
            x = D[s][:, :Nmin]
            xd = x.double()
            EN = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64),
                                         xd * xd], 1), 1)
            NV = Nmin - n
            ew = (EN[:, n:Nmin + 1] - EN[:, :NV + 1])[:, :NV]        # 3-C template-window energy
            e3 = ew.sum(0)
            VAL[s] = (e3 > 1e-6 * e3.median()) & (ew.min(0).values > 0)
            # ---- ENERGY statistic: 6-s S window over the same window 20 s earlier ----
            cs = EN.sum(0)                                            # (Nmin+1,)
            m = Nmin - nE + 1
            eS = (cs[nE:nE + m] - cs[:m])                             # energy of [i, i+nE)
            rat = torch.full((m,), float("nan"), device=dev, dtype=torch.float64)
            rat[nPre:] = torch.log10(eS[nPre:].clamp_min(1e-30) / eS[:m - nPre].clamp_min(1e-30))
            # index by TEMPLATE-START j: the S-1 window starts at i = j - HALF
            zr = torch.full((NV,), float("nan"), device=dev, dtype=torch.float64)
            hi = min(NV, m + HALF)
            zr[HALF:hi] = rat[:hi - HALF]
            g = torch.isfinite(zr) & VAL[s]
            zz = mad_norm(zr.float(), g)
            ZE[s] = zz if zz is not None else torch.full((NV,), float("nan"), device=dev)
            ZE[s] = torch.nan_to_num(ZE[s], nan=-99.0)
            # ---- CC statistic ----
            XH = torch.fft.rfft(x, n=NF)
            fs_ = sorted(TS[s])
            out3 = torch.empty((len(fs_), NV), device=dev)
            e = ew.clamp_min(1e-10 * ew.median()).sqrt().float()
            i0 = 0; ch = a.fam_chunk
            while i0 < len(fs_):
                ff = fs_[i0:i0 + ch]
                W = torch.from_numpy(np.stack([TS[s][f][0] for f in ff])).to(dev)
                WH = torch.fft.rfft(W, n=NF)
                num = torch.fft.irfft(WH.conj() * XH[None], n=NF)[:, :, :NV]
                out3[i0:i0 + len(ff)] = (num / e[None]).mean(1)
                del W, WH, num
                i0 += len(ff)
            md = out3[:, VAL[s]].median(1).values
            sd = (out3 - md[:, None]).abs()[:, VAL[s]].median(1).values * 1.4826
            S3[s] = (fs_, (out3 - md[:, None]) / sd.clamp_min(1e-9)[:, None])
            del XH, out3, EN, ew, e3, e, cs, eS, rat, zr, xd, x, D[s]
            torch.cuda.empty_cache() if dev == "cuda" else None
            # ---- single-station declustered peaks on QUIET days (peaks per station-day) ----
            if doy in QUIET:
                sta_peaks["E"][s].append(decluster_peaks(ZE[s], VAL[s], n).cpu().numpy())
                n_quiet_staday[s] += 1
                fs2, cc2 = S3[s]
                pk = []
                for i in range(len(fs2)):
                    pk.append(decluster_peaks(cc2[i], VAL[s], n).cpu().numpy())
                sta_peaks["C"][s].append(np.concatenate(pk) if pk else np.zeros(0))
                n_quiet_tpl_staday[s] += len(fs2)
        if not S3:
            continue
        NV = min(min(v.shape[0] for v in ZE.values()), min(v[1].shape[1] for v in S3.values()))
        # ---------------- per-family network combination ----------------
        POSd = POS[POS.doy == doy]
        nfam_here = 0
        for f in sorted({fl for s in S3 for fl in TS[s]}):
            ss = [s for s in S3 if f in TS[s]]
            groups = {}
            for s in ss:
                groups.setdefault(SITE[s], []).append(s)
            if len(groups) < 2:
                continue
            nfam_here += 1
            tref0 = min(T0[s] for s in ss)
            k = {s: int(round((tref0 - T0[s] + TS[s][f][1]) * FS)) for s in ss}
            b = max(0, -min(k.values()))
            k = {s: k[s] + b for s in ss}
            tref = tref0 + b / FS               # origin time of network index 0
            NVo = NV - max(k.values())
            if NVo <= 1000:
                continue
            accE = torch.zeros(NVo, device=dev)
            accC = torch.zeros(NVo, device=dev)
            vmask = torch.ones(NVo, dtype=torch.bool, device=dev)
            for q, mem in groups.items():
                sE = torch.zeros(NVo, device=dev); sC = torch.zeros(NVo, device=dev)
                for s in mem:
                    fs2, cc2 = S3[s]
                    sE += ZE[s][k[s]:k[s] + NVo]
                    sC += cc2[fs2.index(f)][k[s]:k[s] + NVo]
                    vmask &= VAL[s][k[s]:k[s] + NVo]
                accE += sE / len(mem); accC += sC / len(mem)
            # ---- positives: max of the network series within +-2 s of the event origin ----
            ev = POSd[POSd.bkey == int(f)]
            if len(ev):
                j = np.rint((ev.t_origin.to_numpy() - tref) * FS).astype(np.int64)
                w = int(round(2.0 * FS))
                ok = (j - w >= 0) & (j + w < NVo)
                for idx, jj, good in zip(ev.index.to_numpy(), j, ok):
                    if not good:
                        continue
                    sl = slice(int(jj - w), int(jj + w + 1))
                    if not bool(vmask[sl].any()):
                        continue
                    net_pos.append(dict(idx=int(idx), fam=int(f), doy=doy,
                                        mag=float(gt.mag.iloc[idx]),
                                        E=float(accE[sl].max()), C=float(accC[sl].max()),
                                        n_sites=len(groups)))
                    # per-station values AT the predicted S (no +-2 s search)
                    for s in ss:
                        jjs = int(jj) + k[s]
                        if 0 <= jjs < NV and bool(VAL[s][jjs]):
                            fs2, cc2 = S3[s]
                            pos_rows.append(dict(idx=int(idx), sta=s, fam=int(f), doy=doy,
                                                 zE=float(ZE[s][jjs]),
                                                 zC=float(cc2[fs2.index(f)][jjs])))
            # ---- FA: declustered network peaks on QUIET family-days ----
            if doy in QUIET:
                pE = decluster_peaks(accE, vmask, n).cpu().numpy()
                pC = decluster_peaks(accC, vmask, n).cpu().numpy()
                fa_peaks["E"].append(pE[pE > 0]); fa_peaks["C"].append(pC[pC > 0])
                n_quiet_famday += 1
            # ---- negatives: sampled origin times (quiet days, and ETS days >= 120 s from any event)
            if doy in QUIET or doy in ETS:
                cand = rng.integers(0, NVo, size=a.n_neg_per_famday * 6)
                tt = tref + cand / FS
                if doy in ETS:
                    ii = np.searchsorted(t_all, tt)
                    lo = t_all[np.clip(ii - 1, 0, len(t_all) - 1)]
                    hi = t_all[np.clip(ii, 0, len(t_all) - 1)]
                    far = np.minimum(np.abs(tt - lo), np.abs(tt - hi)) >= 120.0
                    cand = cand[far]
                cand = cand[:a.n_neg_per_famday]
                cls = "quiet" if doy in QUIET else "ets"
                for jj in cand:
                    jj = int(jj)
                    if not bool(vmask[jj]):
                        continue
                    net_neg.append(dict(cls=cls, fam=int(f), doy=doy,
                                        E=float(accE[jj]), C=float(accC[jj])))
                    for s in ss:
                        jjs = jj + k[s]
                        if 0 <= jjs < NV and bool(VAL[s][jjs]):
                            fs2, cc2 = S3[s]
                            neg_rows.append(dict(cls=cls, sta=s, fam=int(f), doy=doy,
                                                 zE=float(ZE[s][jjs]),
                                                 zC=float(cc2[fs2.index(f)][jjs])))
        nfam_obs.append(nfam_here)
        del S3, ZE, VAL, D
        if dev == "cuda":
            torch.cuda.empty_cache()
        print(f"  [{a.year}.{doy:03d}] fams {nfam_here}  pos {len(net_pos)}  neg {len(net_neg)} "
              f"({time.time()-t_day:.1f}s)", flush=True)

    os.makedirs(f"{RES}/oracle", exist_ok=True)
    P = pd.DataFrame(net_pos); N = pd.DataFrame(net_neg)
    PS = pd.DataFrame(pos_rows); NS = pd.DataFrame(neg_rows)
    P.to_parquet(f"{RES}/oracle/o3_net_pos_{a.tag}.parquet", index=False)
    N.to_parquet(f"{RES}/oracle/o3_net_neg_{a.tag}.parquet", index=False)
    PS.to_parquet(f"{RES}/oracle/o3_sta_pos_{a.tag}.parquet", index=False)
    NS.to_parquet(f"{RES}/oracle/o3_sta_neg_{a.tag}.parquet", index=False)
    np.savez_compressed(f"{RES}/oracle/o3_fa_peaks_{a.tag}.npz",
                        E=np.concatenate(fa_peaks["E"]) if fa_peaks["E"] else np.zeros(0),
                        C=np.concatenate(fa_peaks["C"]) if fa_peaks["C"] else np.zeros(0),
                        n_quiet_famday=n_quiet_famday,
                        nfam_median=float(np.median(nfam_obs)) if nfam_obs else np.nan)
    np.savez_compressed(
        f"{RES}/oracle/o3_sta_peaks_{a.tag}.npz",
        **{f"E_{s}": (np.concatenate(sta_peaks["E"][s]) if sta_peaks["E"][s] else np.zeros(0)) for s in S15},
        **{f"C_{s}": (np.concatenate(sta_peaks["C"][s]) if sta_peaks["C"][s] else np.zeros(0)) for s in S15},
        n_quiet_staday=json.dumps(n_quiet_staday), n_quiet_tpl_staday=json.dumps(n_quiet_tpl_staday))
    json.dump(dict(tag=a.tag, year=a.year, n_days=len(days), n_pos=len(P), n_neg=len(N),
                   n_quiet_famday=n_quiet_famday,
                   nfam_median=float(np.median(nfam_obs)) if nfam_obs else np.nan,
                   n_quiet_staday=n_quiet_staday, n_quiet_tpl_staday=n_quiet_tpl_staday,
                   e_win_s=a.e_win_s, e_pre_s=a.e_pre_s, runtime_s=time.time() - t_start,
                   sites=SITES),
              open(f"{RES}/oracle/o3_meta_{a.tag}.json", "w"), indent=1)
    print(f"[done] {len(P)} positives, {len(N)} negatives, {n_quiet_famday} quiet family-days, "
          f"{time.time()-t_start:.0f}s", flush=True)


if __name__ == "__main__":
    main()
