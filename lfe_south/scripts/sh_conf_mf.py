#!/usr/bin/env python
"""SH_CONF stage 1 (GPU): per-candidate matched-filter features + the site-rolled scrambled control.

For every day of a run (sh_conf_common.day_table) this re-runs the validated detector (t11_mf_v2: its Chan class,
cc formula, validity guard, moveout shift, daily MAD and fixed 4-s bins -- imported, not copied) on the SAME
cleaned data (sh_v2_common prep of that run) with the SAME template set of that day's segment, and in one pass:
  data     re-accumulates the network sum (a check: net_mad / ccsum at every >= 8 x MAD candidate of the stored run
           must be reproduced; max |error| written to the day JSON) and gathers, for every stored candidate and for
           NRAND random origins, each station's channel-summed cc and live-channel count at its aligned position;
  control  ONE scrambled control 'c1': every <= 3-km site's per-channel cc series is rolled (circularly) by its own
           random shift of the same day (sh_conf_common.site_shifts), the network sum rebuilt, daily MAD + 4-s bins
           -> control candidates >= 8 x MAD, and the same per-station gathers at the rolled positions.
Per candidate (data and control) it writes the MF features, per-station cc / z, and the site features
(site_cc_max/2nd, site_z_max/2nd/min, n_z2/3/5, dom_top, site_noise_mad_med); z = (site cc - median) / (1.4826 MAD)
of that site's cc at the NRAND random origins of the same family-day (data, unrolled).

  python scripts/sh_conf_mf.py --shard 0 --nshard 4 [--days 2022-06-30,...] [--runs v3b,fwd]
  -> OUT/feat/<run>/<date>.parquet, OUT/feat/<run>/<date>.json, OUT/feat/<run>/<date>_noise.parquet, OUT/status/mf_<shard>.csv
"""
import os, sys, json, time, argparse, warnings, traceback
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sh_conf_common as C
import sh_v2_common as C2
import t11_mf_v2 as MF
XP = MF.XP; GPU = MF.GPU


def load_station(P, d, ns):
    """DayPrep.load with explicit error capture -> ((t0, X, mc, ref) | None, status, msg)"""
    if ns not in P.stas:
        return None, "noprep", ""
    net, sta = ns.split(".")
    p = C.npz_path(net, sta, d)
    if not os.path.exists(p):
        return None, "nofile", ""
    try:
        z = np.load(p); x = z["x"].astype(np.float64); t0 = float(z["t0"])
    except Exception as e:
        return None, "error", str(e)[:120]
    if x.ndim != 2 or x.shape[1] < 0.5 * 86400 * C.FS or not np.any(x):
        return None, "short", ""
    try:
        X, mc = C2.base_clean(x)
        ms = C2.mute_samples(P.z[f"{ns}|mute"].astype(bool), t0, P.d, X.shape[1])
        mc = mc & ~ms[None, :]
        X = np.where(mc, X, 0.0)
        X = C2.notch(X, mc, P.z[f"{ns}|lines"])
        return (t0, X, mc, P.z[f"{ns}|ref"]), "ok", ""
    except Exception as e:
        return None, "error", str(e)[:120]


def add_shifted(SUM, LIVE, v, lv, sh, NS, nw):
    if sh >= 0:
        m_ = min(NS, nw - sh)
        if m_ > 0: SUM[:m_] += v[sh:sh + m_]; LIVE[:m_] += lv[sh:sh + m_]
    else:
        s_ = -sh; m_ = min(NS - s_, nw)
        if m_ > 0: SUM[s_:s_ + m_] += v[:m_]; LIVE[s_:s_ + m_] += lv[:m_]


def site_feats(stS, stL, sts, SITE, noise_S, noise_L):
    """stS/stL (n_st, n) station channel-sum cc and live count at the candidates; noise_* (n_st, NRAND) at random
    origins (data).  Returns dict of per-candidate site features + station cc / z arrays."""
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
    out = dict(n_sites=len(sites), n_sites_live=(SL > 0).sum(0),
               site_cc_max=top(Cs, 0), site_cc_2nd=top(Cs, 1), site_z_max=top(Zs, 0), site_z_2nd=top(Zs, 1),
               site_z_min=np.nanmin(np.where(np.isfinite(Z), Z, np.inf), 0) if len(sites) else np.full(n, np.nan),
               n_z2=(np.nan_to_num(Z, nan=-9) >= 2).sum(0), n_z3=(np.nan_to_num(Z, nan=-9) >= 3).sum(0),
               n_z5=(np.nan_to_num(Z, nan=-9) >= 5).sum(0),
               dom_top=np.nanmax(np.where(np.isfinite(share), share, -np.inf), 0) if len(sites) else np.full(n, np.nan),
               site_noise_mad_med=float(np.nanmedian(mad)) if np.isfinite(mad).any() else np.nan)
    out["site_z_min"] = np.where(np.isfinite(out["site_z_min"]), out["site_z_min"], np.nan)
    out["dom_top"] = np.where(np.isfinite(out["dom_top"]), out["dom_top"], np.nan)
    return out, stc, stz, smed, smad


def day_job(run, date, seg, tpl_tag, prep_dir, CAND, SITE, odir):
    t00 = time.time()
    d = pd.Timestamp(date)
    C2.PREP = prep_dir
    P = C2.DayPrep(d)
    M, TPL = C.load_tpl(tpl_tag)
    fams = sorted(M.fam.unique(), key=lambda s: int(s))
    stas = sorted(M.netsta.unique())
    status = []
    if P.z is None:
        return None, [(run, date, ns, "error", "no prep file for a processed day") for ns in stas]
    PRE = float(M.pre_s.iloc[0]); LEN = float(M.len_s.iloc[0]); L = int(round(LEN * C.FS)); pre_n = int(round(PRE * C.FS))
    DAT = {}
    for ns in stas:
        v, st, msg = load_station(P, d, ns)
        status.append((run, date, ns, st, msg))
        if v is not None:
            DAT[ns] = v
    if not DAT:
        return None, status
    tref = d.timestamp(); NS = int(round(86400 * C.FS)); gap_n = int(round(C.GAP_S * C.FS))
    N = max(v[1].shape[1] for v in DAT.values())
    nfft = int(2 ** np.ceil(np.log2(N + L)))
    CH = {ns: [MF.Chan(X[c], mc[c], float(ref[c]), L, nfft) for c in range(X.shape[0])] for ns, (t0, X, mc, ref) in DAT.items()}
    SHIFT = C.site_shifts(date, [SITE[s] for s in C.STAS], NS)
    R = C.rand_origins(date, NS)
    cand = CAND.get(date, pd.DataFrame(columns=["fam", "t_epoch", "net_mad", "ccsum", "meancc", "nlive", "nchan", "nsta"]))
    cand_f = {f: g for f, g in cand.groupby("fam")}
    rows = []; noise_rows = []; chk = []
    for f in fams:
        sub = M[(M.fam == f) & M.netsta.isin(DAT)]
        g = cand_f.get(f)
        i_d = np.round((g.t_epoch.to_numpy() - tref) * C.FS).astype(np.int64) if g is not None else np.zeros(0, np.int64)
        Q = np.concatenate([i_d, R]); Qg = XP.asarray(Q)
        SUM = XP.zeros(NS); LIVE = XP.zeros(NS); SUMc = XP.zeros(NS); LIVEc = XP.zeros(NS)
        stQS, stQL, stCS, stCL, stsh = {}, {}, {}, {}, {}
        nch = 0; nst = 0
        for r in sub.itertuples():
            T = TPL.get(r.key)
            if T is None: continue
            t0 = DAT[r.netsta][0]; chs = CH[r.netsta]
            if T.shape[0] != len(chs): continue
            sh = int(round((float(r.tt_s) + float(r.best_tau) - PRE - (t0 - tref)) * C.FS))
            shift = SHIFT[SITE[r.netsta]]
            used = 0
            qs = XP.zeros(len(Q)); ql = XP.zeros(len(Q)); cs = XP.zeros(NS); cl = XP.zeros(NS)
            for c, ch in enumerate(chs):
                tp = T[c]
                if not np.all(np.isfinite(tp)) or not np.any(tp): continue
                tg = XP.asarray(tp)
                st_ = float(np.sqrt(max(L * np.sum(tp * tp) - np.sum(tp) ** 2, 0.0)))
                if st_ <= 0: continue
                dot = XP.fft.irfft(ch.F * XP.conj(XP.fft.rfft(tg, nfft)), nfft)[:ch.nw]
                cc = (L * dot - float(np.sum(tp)) * ch.S1) / (st_ * ch.sd)
                ok = ch.ok & XP.isfinite(cc)
                bad = ok & (XP.abs(cc) > 1 + 1e-6); ok = ok & ~bad
                v = XP.where(ok, cc, 0.0); lv = ok.astype(XP.float64)
                add_shifted(SUM, LIVE, v, lv, sh, NS, ch.nw)
                jq = Qg + sh; okq = (jq >= 0) & (jq < ch.nw); jqc = XP.clip(jq, 0, ch.nw - 1)
                qs += XP.where(okq, v[jqc], 0.0); ql += XP.where(okq, lv[jqc], 0.0)
                k = shift % ch.nw
                vr = XP.roll(v, -k); lvr = XP.roll(lv, -k)
                add_shifted(SUMc, LIVEc, vr, lvr, sh, NS, ch.nw)
                add_shifted(cs, cl, vr, lvr, sh, NS, ch.nw)
                nch += 1; used += 1
                del dot, cc, ok, bad, v, lv, vr, lvr
            if used:
                nst += 1
                stQS[r.netsta] = MF.to_np(qs); stQL[r.netsta] = MF.to_np(ql); stCS[r.netsta] = cs; stCL[r.netsta] = cl
                stsh[r.netsta] = sh
        if nch == 0:
            continue
        SUM = MF.to_np(SUM); LIVE = MF.to_np(LIVE); SUMc = MF.to_np(SUMc); LIVEc = MF.to_np(LIVEc)
        ok = LIVE > 0; okc = LIVEc > 0
        if not ok.any(): continue
        med = np.median(SUM[ok]); mad = np.median(np.abs(SUM[ok] - med)) + 1e-30
        nm = np.where(ok, (SUM - med) / mad, -np.inf)
        medc = np.median(SUMc[okc]); madc = np.median(np.abs(SUMc[okc] - medc)) + 1e-30
        nmc = np.where(okc, (SUMc - medc) / madc, -np.inf)
        keep_c = MF.decl_bin(nmc, C.TMIN, gap_n)
        sts = sorted(stQS)
        nd = len(i_d)
        A = np.stack([stQS[s] for s in sts]); B = np.stack([stQL[s] for s in sts])
        nS, nL = A[:, nd:], B[:, nd:]
        # noise table
        for s_ in sts:
            vv = np.where(nL[sts.index(s_)] > 0, nS[sts.index(s_)] / np.maximum(nL[sts.index(s_)], 1), np.nan)
            m_ = np.nanmedian(vv); noise_rows.append(dict(run=run, date=date, fam=f, netsta=s_, sh=stsh[s_],
                                                          noise_med=m_, noise_mad=1.4826 * np.nanmedian(np.abs(vv - m_)),
                                                          n_valid=int(np.isfinite(vv).sum())))
        base = dict(run=run, date=date, fam=f, seg=seg, nchan=nch, nsta=nst, fam_mad=float(mad), fam_mad_c=float(madc))
        for series, idx, src in (("data", i_d, g), ("c1", keep_c, None)):
            if not len(idx):
                continue
            if series == "data":
                sS, sL = A[:, :nd], B[:, :nd]
                # the stored candidate must be reproduced
                ii = np.clip(idx, 0, NS - 1)
                nm_re = np.where((idx >= 0) & (idx < NS), nm[ii], np.nan)
                cs_re = np.where((idx >= 0) & (idx < NS), SUM[ii], np.nan)
                e1 = np.abs(nm_re - src.net_mad.to_numpy()); e2 = np.abs(cs_re - src.ccsum.to_numpy())
                chk.append(dict(fam=f, n=nd, max_err_mad=float(np.nanmax(e1)) if nd else 0.0,
                                max_err_ccsum=float(np.nanmax(e2)) if nd else 0.0, n_nan=int(np.isnan(nm_re).sum())))
                dfm = dict(t_epoch=src.t_epoch.to_numpy(), net_mad=src.net_mad.to_numpy(), ccsum=src.ccsum.to_numpy(),
                           meancc=src.meancc.to_numpy(), nlive=src.nlive.to_numpy(), net_mad_re=nm_re, ccsum_re=cs_re)
            else:
                ig = XP.asarray(idx)
                sS = np.stack([MF.to_np(stCS[s][ig]) for s in sts]); sL = np.stack([MF.to_np(stCL[s][ig]) for s in sts])
                dfm = dict(t_epoch=tref + idx / C.FS, net_mad=nmc[idx], ccsum=SUMc[idx], meancc=SUMc[idx] / np.maximum(LIVEc[idx], 1),
                           nlive=LIVEc[idx].astype(int), net_mad_re=np.nan, ccsum_re=np.nan)
            sf, stc, stz, _, _ = site_feats(sS, sL, sts, SITE, nS, nL)
            R_ = pd.DataFrame(dict(series=series, i=idx, **dfm, **base, **sf))
            for k_, s_ in enumerate(sts):
                R_[f"cc_{s_}"] = stc[k_].astype(np.float32); R_[f"z_{s_}"] = stz[k_].astype(np.float32)
            rows.append(R_)
        del stCS, stCL
        if GPU: MF.cp.get_default_memory_pool().free_all_blocks()
    del CH
    if GPU: MF.cp.get_default_memory_pool().free_all_blocks()
    F = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    if len(F):
        F["k"] = F.groupby("series").cumcount()
    meta = dict(run=run, date=date, seg=seg, tpl=tpl_tag, prep=prep_dir, NS=NS, L=L, pre_n=pre_n,
                stations={ns: dict(t0=float(DAT[ns][0]), N=int(DAT[ns][1].shape[1]), site=SITE[ns]) for ns in DAT},
                shifts={str(k): v for k, v in SHIFT.items()}, check=chk,
                max_err_mad=max([c["max_err_mad"] for c in chk if c["n"]] or [0.0]),
                max_err_ccsum=max([c["max_err_ccsum"] for c in chk if c["n"]] or [0.0]),
                n_data=int((F.series == "data").sum()) if len(F) else 0, n_c1=int((F.series == "c1").sum()) if len(F) else 0,
                n_cand_stored=int(len(cand)), secs=time.time() - t00)
    os.makedirs(odir, exist_ok=True)
    F.to_parquet(f"{odir}/{date}.parquet.tmp", index=False); os.replace(f"{odir}/{date}.parquet.tmp", f"{odir}/{date}.parquet")
    pd.DataFrame(noise_rows).to_parquet(f"{odir}/{date}_noise.parquet", index=False)
    json.dump(meta, open(f"{odir}/{date}.json.tmp", "w")); os.replace(f"{odir}/{date}.json.tmp", f"{odir}/{date}.json")
    return meta, status


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--days", default=""); ap.add_argument("--runs", default="v3b,fwd")
    a = ap.parse_args(); warnings.filterwarnings("ignore")
    D = pd.read_csv(f"{C.OUT}/day_table.csv", dtype={"date": str})
    D = D[D.run.isin(a.runs.split(","))]
    if a.days: D = D[D.date.isin(a.days.split(","))]
    D = D.reset_index(drop=True); D = D[D.index % a.nshard == a.shard]
    todo = [r for r in D.itertuples() if not os.path.exists(f"{C.OUT}/feat/{r.run}/{r.date}.json")]
    C.log(f"[mf shard {a.shard}/{a.nshard}] {len(todo)} of {len(D)} run-days to do, GPU={GPU}")
    SITE = C.sites()
    CANDS = {}
    status = []; t0 = time.time()
    sp = f"{C.OUT}/status/mf_{a.shard}.csv"
    if os.path.exists(sp):
        status = [tuple(x) for x in pd.read_csv(sp).fillna("").itertuples(index=False)]
    for n, r in enumerate(todo):
        key = (r.run, r.seg)
        if key not in CANDS:
            Cd = C.candidates(r.run, None if r.run == "v3b" else r.seg)
            CANDS[key] = {dd: g for dd, g in Cd.groupby("date")}
        try:
            meta, st = day_job(r.run, r.date, r.seg, r.tpl, r.prep, CANDS[key], SITE, f"{C.OUT}/feat/{r.run}")
            status += st
            if meta is not None and (meta["max_err_mad"] > 1e-4 or meta["max_err_ccsum"] > 1e-6):
                C.log(f"  WARNING {r.run} {r.date}: stored candidate not reproduced, max err mad {meta['max_err_mad']:.3g} ccsum {meta['max_err_ccsum']:.3g}")
            if n % 10 == 0 or n < 3:
                C.log(f"  [{a.shard}] {r.run} {r.date} {r.seg}: data {meta and meta['n_data']} / stored {meta and meta['n_cand_stored']}, "
                      f"c1 {meta and meta['n_c1']}, err {meta and meta['max_err_mad']:.2e}, {meta and meta['secs']:.0f}s; "
                      f"{n + 1}/{len(todo)} {time.time() - t0:.0f}s")
        except Exception as e:
            C.log(f"  ERROR {r.run} {r.date}: {e}\n{traceback.format_exc()[-800:]}")
            status += [(r.run, r.date, ns, "error", f"day job: {str(e)[:100]}") for ns in C.STAS]
            if GPU: MF.cp.get_default_memory_pool().free_all_blocks()
        if n % 20 == 0:
            C.write_status("mf", str(a.shard), status)
    C.write_status("mf", str(a.shard), status)
    C.log(f"[mf shard {a.shard}] done {len(todo)} run-days in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
