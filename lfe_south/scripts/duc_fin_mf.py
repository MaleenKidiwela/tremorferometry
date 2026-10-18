#!/usr/bin/env python
"""DUC_FIN feature matched filter (GPU) -- Ducellier's detector re-run with per-station gathers + ONE scrambled control.

Per day it runs the SAME detector as scripts/duc_fix_mf.py (its DayData, decluster and constants are imported; the
per-family loop below is a line-for-line copy of duc_fix_mf.run with extra gathers that never touch the detection):
his templates for 'his' epochs + the DUC_FIN re-templated new epochs (OUT/retpl/templates_new, keep == 1).
  data   every >= 8 x MAD hourly detection (= duc_fix mf1 where no new template is active; the 2011-10..2012-07 days
         with new templates are the pass-2 re-scan), with: cc (channel mean), nchannel, nlive, mad, ccsum, thr_perm,
         pass_his, the number of templates (n_tpl) and of new-epoch templates (n_new_tpl); per station: cc_<NS> (mean cc
         of its participating channels at the detection lag), nc_<NS> (participating channels), z_<NS> (station cc vs
         that station's cc at NRAND random origins of the same family-day); site features as SH_CONF (each of his 47
         stations is its own <= 3-km site): site_cc_max/2nd, site_z_max/2nd/min, n_z2/3/5, dom_top, site_noise_mad_med.
  c1     scrambled control: every site's per-channel cc series (and its hour participation) rolled by its own random
         shift of the same day (duc_fin_common.site_shifts); network sum rebuilt; per-sample channel count; hourly MAD,
         8 x MAD, his declustering; same gathers at the rolled positions.
Station-day status (ok / nodata / short / missing / error) -> OUT/status/mf_<shard>.csv; guard in the queue.

  python scripts/duc_fin_mf.py --shard 0 --nshard 3 [--days 2009-06-01,...]   -> OUT/feat/<date>.parquet + .json
"""
import os, sys, json, time, argparse, warnings, traceback
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import duc_fin_common as F
import duc_fix_common as DC
import duc_fix_mf as MF
from duc_fix_mf import THR, FREQ0, NH, NDAY, TAIL, NREF, decluster
import duc_fix_clean as CL
XP = MF.XP; GPU = MF.GPU
if GPU:
    try: MF.cp.fft.config.get_plan_cache().set_size(4)       # the day uses ONE fft length; keep the cache tiny
    except Exception: pass
FS = DC.FS


def new_templates():
    p = f"{F.OUT}/retpl/templates_new.csv"
    if not os.path.exists(p): return {}
    M = pd.read_csv(p, dtype={"fam": str}); M = M[M.keep == 1]
    Z = np.load(f"{F.OUT}/retpl/templates_new.npz")
    return {(r.fam, r.epoch, r.cha): Z[r.key].astype(np.float64) for r in M.itertuples() if r.key in Z.files}


def site_feats(stS, stL, sts, SITE, noise_S, noise_L):
    """sh_conf_mf.site_feats (copied): stS/stL (n_st, n) station channel-sum cc and channel count at the candidates;
    noise_* (n_st, NRAND) at random origins (data)."""
    sites = sorted({SITE[s] for s in sts})
    si = np.array([sites.index(SITE[s]) for s in sts])
    n = stS.shape[1]
    SS = np.zeros((len(sites), n)); SL = np.zeros((len(sites), n))
    RS = np.zeros((len(sites), noise_S.shape[1])); RL = np.zeros_like(RS)
    for k in range(len(sites)):
        m = si == k
        SS[k] = stS[m].sum(0); SL[k] = stL[m].sum(0)
        RS[k] = noise_S[m].sum(0); RL[k] = noise_L[m].sum(0)
    with np.errstate(invalid="ignore", divide="ignore"):
        SC = np.where(SL > 0, SS / SL, np.nan)
        RC = np.where(RL > 0, RS / RL, np.nan)
        med = np.nanmedian(RC, 1) if RC.size else np.zeros(0)
        mad = 1.4826 * np.nanmedian(np.abs(RC - med[:, None]), 1)
        mad = np.where(mad > 0, mad, np.nan)
        Z = (SC - med[:, None]) / mad[:, None]
        tot = SS.sum(0)
        share = SS / np.where(tot > 0, tot, np.nan)
        stc = np.where(stL > 0, stS / np.maximum(stL, 1), np.nan)
        RSt = np.where(noise_L > 0, noise_S / np.maximum(noise_L, 1), np.nan)
        smed = np.nanmedian(RSt, 1); smad = 1.4826 * np.nanmedian(np.abs(RSt - smed[:, None]), 1)
        smad = np.where(smad > 0, smad, np.nan)
        stz = (stc - smed[:, None]) / smad[:, None]
    Cs = np.sort(np.nan_to_num(SC, nan=-99), 0)[::-1]
    Zs = np.sort(np.nan_to_num(Z, nan=-99), 0)[::-1]
    def top(A, k):
        if A.shape[0] <= k: return np.full(n, np.nan)
        v = A[k].astype(float); v[v <= -98] = np.nan; return v
    with np.errstate(invalid="ignore"):
        zmin = np.nanmin(np.where(np.isfinite(Z), Z, np.inf), 0) if len(sites) else np.full(n, np.nan)
        dom = np.nanmax(np.where(np.isfinite(share), share, -np.inf), 0) if len(sites) else np.full(n, np.nan)
    out = dict(n_sites=len(sites), n_sites_live=(SL > 0).sum(0),
               site_cc_max=top(Cs, 0), site_cc_2nd=top(Cs, 1), site_z_max=top(Zs, 0), site_z_2nd=top(Zs, 1),
               site_z_min=np.where(np.isfinite(zmin), zmin, np.nan),
               n_z2=(np.nan_to_num(Z, nan=-9) >= 2).sum(0), n_z3=(np.nan_to_num(Z, nan=-9) >= 3).sum(0),
               n_z5=(np.nan_to_num(Z, nan=-9) >= 5).sum(0),
               dom_top=np.where(np.isfinite(dom), dom, np.nan),
               site_noise_mad_med=float(np.nanmedian(mad)) if np.isfinite(mad).any() else np.nan)
    return out, stc, stz


def station_status(d, want, look, DD, plan):
    rows = []; ds = str(d.date())
    loaded = {k[0] for k in DD.ch}
    for ns in want:
        net, sta = ns.split(".")
        if look(sta, d) is None: continue
        st = plan.get((ns, ds))
        if sta in loaded: rows.append((ds, ns, "ok", "")); continue
        if st not in ("built", "cached"): rows.append((ds, ns, "nodata", str(st))); continue
        p = DC.cache_path(net, sta, d)
        if not os.path.exists(p): rows.append((ds, ns, "missing", "plan says on disk, no cache file")); continue
        try:
            z = np.load(p); x = z["x"]
            if x.ndim != 2 or x.shape[1] < 0.5 * 86400 * FS or not np.any(x):
                rows.append((ds, ns, "short", "cache < 12 h or all zero")); continue
        except Exception as e:
            rows.append((ds, ns, "error", f"cache unreadable: {str(e)[:60]}")); continue
        rows.append((ds, ns, "missing", "cache readable but not loaded (prep / t0)"))
    return rows


def day_job(d, G, cache):
    t00 = time.time(); d = pd.Timestamp(d); date = str(d.date())
    fams, want, look, orient, HT, NT, famst, dur, thr, SITE = (G[k] for k in
        ["fams", "want", "look", "orient", "HT", "NT", "famst", "dur", "thr", "SITE"])
    for k in [k for k in cache if k[-1] < d]: cache.pop(k)
    DD = MF.DayData(d, want, look, orient, cache)
    status = station_status(d, want, look, DD, G["plan"])
    meta = dict(date=date, stations=sorted({f"{G['NETOF'][k[0]]}.{k[0]}" for k in DD.ch}), n_channels=len(DD.ch))
    if not DD.ch:
        return pd.DataFrame(), meta, status
    nfft = int(2 ** np.ceil(np.log2(NDAY + TAIL + 1000)))
    SH = F.site_shifts(date, list(SITE.values()), NDAY)
    meta["shifts"] = {str(k): v for k, v in SH.items()}
    R = F.rand_origins(date, NDAY)
    Rg = XP.asarray(R); Rh = R // NH
    Gc = {}
    def gch(k):
        if k in Gc: return Gc[k]
        x, m, r, ref, ep = DD.ch[k]
        xg = XP.asarray(x)
        cr = np.concatenate([[0], np.cumsum(~r)])
        Gc[k] = dict(F=XP.fft.rfft(xg, nfft), x=xg, x2=xg * xg, m=XP.asarray(m), cr=cr, ref=ref, ep=ep, S={})
        return Gc[k]
    PH = XP.asarray(np.minimum(np.arange(NDAY) // NH, 23))
    hidx_np = (np.arange(24)[:, None] * NH + np.arange(NH + 1)[None, :])
    hidx = XP.asarray(hidx_np); hidx_c = XP.asarray(hidx_np % NDAY)
    out = []; lastL = None; gpu_peak = 0; NTPL = {}
    for f in sorted(fams, key=lambda f_: dur[f_]):
        L = int(round(dur[f] * FS)) + 1
        if L != lastL:
            for g_ in Gc.values(): g_["S"].clear()
            lastL = L
        TT = []
        for sta in famst[f]:
            for (s2, c), (x, m, r, ref, ep) in DD.ch.items():
                if s2 != sta: continue
                if ep.status == "his": tp = HT.get(f, {}).get((sta, c))
                else: tp = NT.get((f, ep.epoch, c))
                if tp is None or tp.size != L or not np.all(np.isfinite(tp)) or not np.any(tp): continue
                TT.append(((sta, c), tp))
        if not TT: continue
        NTPL[f] = [len(TT), int(sum(1 for k_, _ in TT if DD.ch[k_][4].status == "new"))]
        nw = NDAY + 1
        SUM = XP.zeros((24, NH + 1)); LIVE = XP.zeros((24, NH + 1)); NCH = np.zeros(24)
        SUMc = XP.zeros(NDAY); NCc = XP.zeros(NDAY); LVc = XP.zeros(NDAY)
        ngt1 = 0; CHI = []
        for k, tp in TT:
            g = gch(k)
            if g["ref"] <= 0: continue
            part = np.array([g["cr"][h * NH + NH + L] - g["cr"][h * NH] == 0
                             if h * NH + NH + L <= NDAY + TAIL else False for h in range(24)])
            if not part.any(): continue
            if L not in g["S"]:
                starts = XP.arange(nw)
                S1 = CL.box_sum(g["x"][None, :], L, starts, XP)[0]
                S2 = CL.box_sum(g["x2"][None, :], L, starts, XP)[0]
                inv = (~g["m"]).astype(XP.int32)
                cf = XP.concatenate([XP.zeros(1, dtype=XP.int64), XP.cumsum(inv)])
                ninv = cf[L:L + nw] - cf[:nw]
                var = L * S2 - S1 * S1
                ok = (ninv == 0) & (S2 > CL.GUARD_REL * g["ref"] * L / NREF) & (var > 0)
                g["S"][L] = (S1, XP.sqrt(XP.where(ok, var, 1.0)), ok)
            S1, sd, ok = g["S"][L]
            tg = XP.asarray(tp)
            st = float(np.sqrt(max(L * np.sum(tp * tp) - np.sum(tp) ** 2, 0.0)))
            if st <= 0: continue
            dot = XP.fft.irfft(g["F"] * XP.conj(XP.fft.rfft(tg, nfft)), nfft)[:nw]
            cc = (L * dot - float(np.sum(tp)) * S1) / (st * sd)
            okc = ok & XP.isfinite(cc)
            bad = okc & (XP.abs(cc) > 1 + 1e-6)
            ngt1 += int(bad.sum()); okc = okc & ~bad
            v = XP.where(okc, cc, 0.0)[hidx]; lv = okc.astype(XP.float64)[hidx]
            pg = XP.asarray(part.astype(np.float64))[:, None]
            SUM += v * pg; LIVE += lv * pg; NCH += part
            # ---- extra (never touches the detection above): keep the channel, add it to the rolled control
            vfull = XP.where(okc, cc, 0.0)
            pf = XP.asarray(part)[PH].astype(XP.float64)                  # per-sample hour participation (data frame)
            sta = k[0]; ns = f"{G['NETOF'][sta]}.{sta}"; s_ = SH[SITE[ns]]
            SUMc += XP.roll(vfull[:NDAY] * pf, -s_); NCc += XP.roll(pf, -s_); LVc += XP.roll(okc[:NDAY].astype(XP.float64) * pf, -s_)
            CHI.append((ns, vfull.astype(XP.float32), okc, part, pf.astype(XP.bool_), s_,
                        int(DD.ch[k][4].status == "new")))
            del dot, cc, bad
        if not NCH.any() or not CHI: continue
        SUMn = MF.to_np(SUM); LIVEn = MF.to_np(LIVE)
        thrp = float(thr.get(f, 0.0))
        det = []
        for h in range(24):
            if NCH[h] == 0: continue
            mcc = SUMn[h] / NCH[h]
            mad = float(np.median(np.abs(mcc - mcc.mean())))
            idx = np.where(mcc >= THR * mad)[0]
            if idx.size == 0 or mad <= 0: continue
            keep = decluster(idx, mcc, int(np.floor(1.0 / ((1 / FS) * FREQ0))))
            for i in keep:
                c_ = float(mcc[i]); nc = int(NCH[h])
                det.append((h, int(i), c_, nc, int(LIVEn[h, i]), mad, c_ * nc, thrp, int(thrp > 0 and c_ * nc >= thrp)))
        # ---- control detection (rolled frame, per-sample channel count)
        Sc = MF.to_np(SUMc[hidx_c]); Nc = MF.to_np(NCc[hidx_c]); Lc = MF.to_np(LVc[hidx_c])
        detc = []
        for h in range(24):
            if not (Nc[h] > 0).any(): continue
            mcc = np.where(Nc[h] > 0, Sc[h] / np.maximum(Nc[h], 1), 0.0)
            mad = float(np.median(np.abs(mcc - mcc.mean())))
            idx = np.where(mcc >= THR * mad)[0]
            if idx.size == 0 or mad <= 0: continue
            keep = decluster(idx, mcc, int(np.floor(1.0 / ((1 / FS) * FREQ0))))
            for i in keep:
                c_ = float(mcc[i]); nc = int(Nc[h, i])
                detc.append((h, int(i), c_, nc, int(Lc[h, i]), mad, c_ * nc, thrp, int(thrp > 0 and c_ * nc >= thrp)))
        # ---- gathers: station channel-sum and channel count at the candidates (data / control) and random origins
        sts = sorted({c[0] for c in CHI}); si = {s: j for j, s in enumerate(sts)}
        def gather(hh, ii, rolled):
            n = len(hh)
            S_ = np.zeros((len(sts), n)); P_ = np.zeros((len(sts), n)); L_ = np.zeros((len(sts), n))
            if not n: return S_, P_, L_
            g0 = np.asarray(hh, np.int64) * NH + np.asarray(ii, np.int64)
            for ns, v32, okc, part, pf, s_, _ in CHI:
                j = si[ns]
                if rolled:
                    go = XP.asarray((g0 % NDAY + s_) % NDAY)
                    pp = MF.to_np(pf[go]).astype(float)
                else:
                    go = XP.asarray(np.minimum(g0, NDAY))
                    pp = part[np.asarray(hh)].astype(float)
                S_[j] += MF.to_np(v32[go]).astype(float) * pp; P_[j] += pp; L_[j] += MF.to_np(okc[go]).astype(float) * pp
            return S_, P_, L_
        nS, nP, _ = gather(Rh, R % NH, False)
        for series, D_ in (("data", det), ("c1", detc)):
            if not D_: continue
            A = np.array(D_, dtype=float)
            hh = A[:, 0].astype(int); ii = A[:, 1].astype(int)
            S_, P_, L_ = gather(hh, ii, series == "c1")
            sf, stc, stz = site_feats(S_, P_, sts, SITE, nS, nP)
            g0 = hh * NH + ii
            Rw = pd.DataFrame(dict(series=series, fam=f, h=hh, i=ii, g=g0, t_epoch=d.timestamp() + g0 / FS,
                                   cc=A[:, 2], nchannel=A[:, 3].astype(int), nlive=A[:, 4].astype(int), mad=A[:, 5],
                                   ccsum=A[:, 6], thr_perm=A[:, 7], pass_his=A[:, 8].astype(int),
                                   n_tpl=len(TT), n_new_tpl=int(sum(c[6] for c in CHI)), nsta=len(sts), **sf))
            for ns in sts:
                j = si[ns]
                Rw[f"cc_{ns}"] = stc[j].astype(np.float32); Rw[f"nc_{ns}"] = P_[j].astype(np.int8); Rw[f"z_{ns}"] = stz[j].astype(np.float32)
            out.append(Rw)
        if GPU:
            gpu_peak = max(gpu_peak, MF.cp.get_default_memory_pool().total_bytes())
        del CHI, SUMc, NCc, LVc
    del Gc
    if GPU:
        gpu_peak = max(gpu_peak, MF.cp.get_default_memory_pool().total_bytes())
        MF.cp.get_default_memory_pool().free_all_blocks()
    Fo = pd.concat(out, ignore_index=True) if out else pd.DataFrame()
    if len(Fo): Fo["k"] = np.arange(len(Fo))
    meta["ntpl"] = NTPL
    meta.update(n_data=int((Fo.series == "data").sum()) if len(Fo) else 0, n_c1=int((Fo.series == "c1").sum()) if len(Fo) else 0,
                n_pass_data=int(Fo[Fo.series == "data"].pass_his.sum()) if len(Fo) else 0,
                n_pass_c1=int(Fo[Fo.series == "c1"].pass_his.sum()) if len(Fo) else 0,
                gpu_pool_peak_GB=round(gpu_peak / 1e9, 2), secs=round(time.time() - t00, 1))
    return Fo, meta, status


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--days", default="")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    if a.days: days = [x for j, x in enumerate(sorted(a.days.split(","))) if j % a.nshard == a.shard]
    else:
        ms = [m for j, m in enumerate(F.MONTHS) if j % a.nshard == a.shard]       # contiguous months per shard
        days = [x for x in F.day_table().date if x[:7] in set(ms)]
    od = f"{F.OUT}/feat"; os.makedirs(od, exist_ok=True)
    todo = [x for x in days if not os.path.exists(f"{od}/{x}.json")]
    fp = DC.families(); H = DC.his_rows(); NETOF = dict(zip(H.sta, H.net))
    famst = {r.family: r.stations.split(",") for r in fp.itertuples()}
    G = dict(fams=list(fp.family), famst=famst, dur={r.family: float(r.duration) for r in fp.itertuples()},
             thr=DC.thresholds(), HT=DC.his_templates(), NT=new_templates(), orient=MF.template_orient(), look=DC.epoch_lookup(),
             NETOF=NETOF, SITE=F.sites(), plan=F.plan_status())
    G["want"] = sorted({f"{NETOF[s]}.{s}" for f in G["fams"] for s in famst[f] if s in NETOF})
    F.log(f"[fin-mf {a.shard}/{a.nshard}] {len(todo)} of {len(days)} days to do; GPU={GPU}; new templates {len(G['NT'])}; "
          f"DUC_MUTE={os.environ.get('DUC_MUTE')} (mute on={CL.MUTE_ON})")
    sp = f"{F.OUT}/status/mf_{a.shard}.csv"
    status = [tuple(x) for x in pd.read_csv(sp).fillna("").itertuples(index=False)] if os.path.exists(sp) else []
    cache = {}; t0 = time.time()
    for n, ds in enumerate(todo):
        try:
            Fo, meta, st = day_job(ds, G, cache)
            status += st
            Fo.to_parquet(f"{od}/{ds}.parquet.tmp", index=False); os.replace(f"{od}/{ds}.parquet.tmp", f"{od}/{ds}.parquet")
            json.dump(meta, open(f"{od}/{ds}.json.tmp", "w")); os.replace(f"{od}/{ds}.json.tmp", f"{od}/{ds}.json")
            if n < 3 or n % 15 == 0:
                F.log(f"  [{a.shard}] {ds}: data {meta.get('n_data')} (pass {meta.get('n_pass_data')}), c1 {meta.get('n_c1')} "
                      f"(pass {meta.get('n_pass_c1')}), stations {len(meta['stations'])}, gpu {meta.get('gpu_pool_peak_GB')} GB, "
                      f"{meta.get('secs')}s; {n + 1}/{len(todo)} {time.time() - t0:.0f}s")
        except Exception as e:
            F.log(f"  ERROR {ds}: {e}\n{traceback.format_exc()[-1500:]}")
            status += [(ds, ns, "error", f"day job: {str(e)[:100]}") for ns in G["want"]]
            cache.clear()
            if GPU: MF.cp.get_default_memory_pool().free_all_blocks()
        if n % 10 == 0: F.write_status("mf", str(a.shard), status)
    F.write_status("mf", str(a.shard), status)
    F.log(f"[fin-mf {a.shard}] done {len(todo)} days {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
