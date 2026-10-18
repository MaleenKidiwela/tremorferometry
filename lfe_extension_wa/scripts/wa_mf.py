#!/usr/bin/env python
"""Stage 4: permanent-only network matched filter on the FROZEN chain (lfe_extension c_run_month.py
measurement core), with the site-rolled scrambled control computed on the same pass.

Per day:
  data       cache20hz/<NET.STA>/<YYYY>/<DOY>.npz, midnight-aligned (so every station shares t0)
  templates  one set per (family, station EPOCH NODE): a station-day is matched with the template of
             ITS node only; no node template -> the station does not contribute (never silently blind)
  cc         full-day normalised correlation per live channel by rFFT (float64 prefix-sum energy),
             station series = mean over the channels live in both template and data
  guard      c_run_month ROBUST ENERGY FLOOR (per channel vs the station-day median, cc set to 0 at
             masked windows, dilated by the template length) + windows touching the cache validity
             mask `m`
  network    site merge B: within-site w-weighted mean cc, site weight sqrt(sum w^2), w = n_channels x
             SNR weight; network = sum_site W_site x site mean / sum W_site, on the ORIGIN axis
             (template start at station = t + offset_s)
  observable >= 2 SITES with a template and data that day, else the family-day is not scored
  background median and MAD of that family-day's own network series over valid samples
  detections declustered peaks (max over a template-length window) >= med + 8.0 x MAD (both tiers are
             filters on net_mad later: T_ext 8.25, T_core 8.75)
  control    SITE-ROLLED scrambled control: each SITE's series (all its stations together) is rolled by
             its own shift drawn in [300, 3000] s for that day, which destroys between-site coherence
             and keeps each site's own statistics (e4_template_test.py E6 shift_site); same weights,
             same valid masks (rolled with the series), same MAD and declustering.

  python scripts/wa_mf.py --tag val --windows V2008,V2010 --tpl build
  -> results/mf_<tag>/det_<date>.parquet, famday_<date>.parquet, run_<date>.json
"""
from __future__ import annotations
import os, sys, json, time, argparse
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

LOG = f"{W.STUDY}/logs/mf.log"


def load_tpl(tags, all_templates=False):
    """clear templates only (default, the trial rule); all_templates=True (Merlin T1, 2026-09-30): EVERY
    >= 30-event (family, node) template, weighted by its SNR weight w, clear flag kept as a column."""
    T, Q = {}, []
    for tag in tags:
        z = np.load(f"{W.RES}/templates/tpl_{tag}.npz", allow_pickle=False)
        q = pd.read_csv(f"{W.RES}/templates/tpl_qc_{tag}.csv", dtype={"live": str})
        q["live"] = q.live.astype(str).str.zfill(3)
        if all_templates:
            q = q[(q.n_events >= W.MIN_EVENTS) & q.offset_s.notna() & q.w.notna()].copy()
        else:
            q = q[q.clear.fillna(False).astype(bool)].copy()
        q["tag"] = tag
        keys = list(z["keys"])
        tp = z["tpl"]
        kk = {k: i for i, k in enumerate(keys)}
        for r in q.itertuples():
            k = f"{r.fid}|{r.netsta}|{r.node}"
            if (int(r.fid), r.netsta, r.node) in T:      # step 4b: the FIRST tag listed wins (earlier templates kept)
                continue
            T[(int(r.fid), r.netsta, r.node)] = dict(tpl=tp[kk[k]], off=float(r.offset_s), w=float(r.w),
                                                    live=np.array([c == "1" for c in r.live]), tag=tag,
                                                    clear=bool(r.clear))
        Q.append(q)
    return T, pd.concat(Q, ignore_index=True)


def main():
    import torch
    import torch.nn.functional as Fnn
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", required=True)
    ap.add_argument("--windows", required=True)
    ap.add_argument("--tpl", required=True, help="comma list of template tags")
    ap.add_argument("--exclude-nodes", default="", help="csv with `node`: never used (e.g. trial pass 1)")
    ap.add_argument("--fam-chunk", type=int, default=24)
    ap.add_argument("--energy-floor", type=float, default=1e-3)
    ap.add_argument("--seed", type=int, default=20260929)
    ap.add_argument("--days", default="", help="comma list of dates (smoke)")
    ap.add_argument("--no-control", action="store_true")
    ap.add_argument("--n-control", type=int, default=1, help="2 = add a second, independent site-rolled "
                    "control ('control2', own shifts) so the separation statistic can be applied to a control")
    ap.add_argument("--all-templates", action="store_true", help="Merlin T1: every >= 30-event template (w-weighted)")
    ap.add_argument("--drop-dead", type=float, default=0.0, help="> 0: a station-day whose guarded (invalid) "
                    "fraction exceeds this is treated as ABSENT (no data) instead of blanking every family it "
                    "carries (the all-stations-valid network rule); 0 = off (trial behaviour)")
    a = ap.parse_args()
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    od = f"{W.RES}/mf_{a.tag}"
    os.makedirs(od, exist_ok=True)
    S = W.load_stations()
    SITE = dict(zip(S.netsta, S.site))
    TPL, Q = load_tpl(a.tpl.split(","), a.all_templates)
    excl = set(pd.read_csv(a.exclude_nodes).node) if a.exclude_nodes else set()
    fams = W.load_families()
    days = [(w, d, k) for w in a.windows.split(",") for d, k in W.window_days(w)]
    if a.days:
        keep = set(a.days.split(","))
        days = [x for x in days if x[1] in keep]
    n = W.NT
    NV = W.NDAY - n + 1
    NF = int(2 ** np.ceil(np.log2(W.NDAY + n)))
    def guard(X, live, m):
        """frozen c_run_month robust energy floor + cache validity mask, dilated by the template length."""
        xd = X.double()
        EN = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64), xd * xd], 1), 1)
        ew = (EN[:, n:] - EN[:, :-n])[:, :NV]                 # (3, NV) window energy per channel
        del EN, xd
        lv = torch.from_numpy(live).to(dev)
        ewl = ew[lv]
        med_ch = ewl.median(dim=1).values
        e3 = ewl.sum(0)
        med3 = e3.median()
        ref_ch = torch.maximum(med_ch, 1e-2 * med3)
        bad = ((ewl < a.energy_floor * ref_ch[:, None]).any(0) | (e3 < a.energy_floor * med3)
               | (ewl.min(0).values <= 0))
        if float(med3) <= 0:
            bad = torch.ones_like(bad)
        mk = torch.from_numpy(~m).to(dev)
        cmk = torch.cat([torch.zeros(1, device=dev), torch.cumsum(mk.float(), 0)])
        bad |= (cmk[n:n + NV] - cmk[:NV]) > 0
        if bool(bad.any()):
            bad = Fnn.max_pool1d(bad.float()[None, None], kernel_size=2 * n + 1, stride=1, padding=n)[0, 0] > 0
        efl = 1e-10 * ewl.median()
        E = ew.clamp_min(efl).sqrt().float()                  # (3, NV)
        return bad, E

    W.log(LOG, f"[mf:{a.tag}] {len(days)} days, {len(TPL)} templates ({Q.fid.nunique()} families), "
          f"excluded nodes {sorted(excl)}, device {dev}")
    for (win, ds, kind) in days:
        outp = f"{od}/det_{ds}.parquet"
        if os.path.exists(outp):
            continue
        t_day = time.time()
        D0 = pd.Timestamp(ds).timestamp()
        rng = np.random.default_rng(a.seed + int(pd.Timestamp(ds).strftime("%Y%m%d")))
        # ---------------- load the station-days
        D = {}
        for ns in S.netsta:
            p = W.cache_path(ns, ds)
            if not os.path.exists(p):
                continue
            z = np.load(p)
            node = str(z["node"])
            if node in excl:
                continue
            D[ns] = dict(x=z["x"], m=z["m"], live=z["live"], node=node)
        dead = {}
        loaded = list(D)                     # scrambled-control shifts are drawn over the LOADED sites (unchanged)
        if a.drop_dead > 0:
            for ns in list(D):
                if not D[ns]["live"].any():
                    dead[ns] = 1.0; del D[ns]; continue
                bad_, E_ = guard(torch.from_numpy(D[ns]["x"]).to(dev), D[ns]["live"], D[ns]["m"])
                fr_ = float(bad_.float().mean())
                del bad_, E_
                if fr_ > a.drop_dead:
                    dead[ns] = fr_; del D[ns]
        # families x stations present with a template of this node
        pres = {}
        for (f, ns, node), t in TPL.items():
            if ns in D and D[ns]["node"] == node:
                ch = t["live"] & D[ns]["live"]
                if ch.any():
                    pres.setdefault(f, []).append((ns, ch))
        fam_rows, obs = [], {}
        for f in sorted(fams.fid):
            ss = pres.get(f, [])
            sites = sorted({SITE[ns] for ns, _ in ss})
            fam_rows.append(dict(fid=f, date=ds, window=win, kind=kind, n_sta=len(ss), n_sites=len(sites),
                                 sites=",".join(sites), stations=",".join(ns for ns, _ in ss),
                                 observable=len(sites) >= 2))
            if len(sites) >= 2:
                obs[f] = ss
        FL = sorted(obs)
        fi = {f: i for i, f in enumerate(FL)}
        if not FL:
            pd.DataFrame(fam_rows).to_parquet(f"{od}/famday_{ds}.parquet", index=False)
            pd.DataFrame(columns=["fid", "t"]).to_parquet(outp, index=False)
            continue
        # coefficients (site merge B)
        coef = {}
        for f, ss in obs.items():
            wm = {ns: float(ch.sum()) * TPL[(f, ns, D[ns]["node"])]["w"] for ns, ch in ss}
            groups = {}
            for ns, _ in ss:
                groups.setdefault(SITE[ns], []).append(ns)
            Ws = {q: float(np.sqrt(sum(wm[s] ** 2 for s in mem))) for q, mem in groups.items()}
            Wt = sum(Ws.values())
            for q, mem in groups.items():
                wsum = sum(wm[s] for s in mem)
                for s in mem:
                    coef[(f, s)] = Ws[q] * wm[s] / wsum / Wt
        sites_today = sorted({SITE[ns] for ns in loaded})
        shift = {q: int(rng.uniform(*W.ROLL_S) * W.FS) for q in sites_today}
        C2 = (not a.no_control) and a.n_control >= 2
        if C2:
            rng2 = np.random.default_rng(a.seed + 7 + int(pd.Timestamp(ds).strftime("%Y%m%d")))
            shift2 = {q: int(rng2.uniform(*W.ROLL_S) * W.FS) for q in sites_today}
        NET = torch.zeros((len(FL), NV), device=dev)
        VM = torch.ones((len(FL), NV), dtype=torch.bool, device=dev)
        if not a.no_control:
            CTL = torch.zeros((len(FL), NV), device=dev)
            CM = torch.ones((len(FL), NV), dtype=torch.bool, device=dev)
        if C2:
            CTL2 = torch.zeros((len(FL), NV), device=dev)
            CM2 = torch.ones((len(FL), NV), dtype=torch.bool, device=dev)
        n_badw = {}
        by_sta = {}
        for f, ss in obs.items():
            for ns, ch in ss:
                by_sta.setdefault(ns, []).append((f, ch))
        for ns, lst in by_sta.items():
            X = torch.from_numpy(D[ns]["x"]).to(dev)
            DH = torch.fft.rfft(X, n=NF)
            bad, E = guard(X, D[ns]["live"], D[ns]["m"])
            del X; X = None
            n_badw[ns] = int(bad.sum())
            valid = ~bad
            # group families by channel pattern
            pats = {}
            for f, ch in lst:
                pats.setdefault(tuple(ch.tolist()), []).append(f)
            for pat, fl in pats.items():
                chi = [c for c in range(3) if pat[c]]
                for i0 in range(0, len(fl), a.fam_chunk):
                    ff = fl[i0:i0 + a.fam_chunk]
                    Wt_ = torch.from_numpy(np.stack([TPL[(f, ns, D[ns]["node"])]["tpl"][chi] for f in ff])).to(dev)
                    WH = torch.fft.rfft(Wt_, n=NF)
                    num = torch.fft.irfft(WH.conj() * DH[chi][None], n=NF)[:, :, :NV]
                    cc = (num / E[chi][None]).mean(1)              # (nf, NV)
                    del WH, num, Wt_
                    cc[:, bad] = 0.0
                    for j, f in enumerate(ff):
                        c = coef[(f, ns)]
                        d = int(round(TPL[(f, ns, D[ns]["node"])]["off"] * W.FS))
                        a0, b0 = max(0, -d), min(NV, NV - d)
                        tmp = torch.zeros(NV, device=dev)
                        tmp[a0:b0] = cc[j, a0 + d:b0 + d]
                        tv = torch.zeros(NV, dtype=torch.bool, device=dev)
                        tv[a0:b0] = valid[a0 + d:b0 + d]
                        k = fi[f]
                        NET[k] += c * tmp
                        VM[k] &= tv
                        if not a.no_control:
                            sh = shift[SITE[ns]] % NV
                            CTL[k] += c * torch.roll(tmp, sh)
                            CM[k] &= torch.roll(tv, sh)
                        if C2:
                            sh2 = shift2[SITE[ns]] % NV
                            CTL2[k] += c * torch.roll(tmp, sh2)
                            CM2[k] &= torch.roll(tv, sh2)
                        del tmp, tv
                    del cc
            del DH, E, bad, valid
        # ---------------- detections
        rows = []
        famstat = {}
        series = [("data", NET, VM)] + ([] if a.no_control else [("control", CTL, CM)]) + \
            ([("control2", CTL2, CM2)] if C2 else [])
        for lab, SER, MSK in series:
            for f in FL:
                k = fi[f]
                x = SER[k]; vm = MSK[k]
                nvld = int(vm.sum())
                if nvld < 1000:
                    famstat[(f, lab)] = (np.nan, np.nan, nvld, 0)
                    continue
                am = x[vm]
                med = am.median(); mad = (am - med).abs().median().clamp_min(1e-6)
                xx = torch.where(vm, x, torch.full_like(x, -1e9))
                mx = Fnn.max_pool1d(xx[None, None], kernel_size=2 * (n // 2) + 1, stride=1, padding=n // 2)[0, 0]
                idx = torch.nonzero((xx >= mx - 1e-12) & vm & (xx >= med + W.T_REC * mad)).squeeze(1)
                famstat[(f, lab)] = (float(med), float(mad), nvld, int(idx.numel()))
                if idx.numel():
                    v = x[idx]
                    nm = ((v - med) / mad).cpu().numpy()
                    ii = idx.cpu().numpy()
                    for j in range(len(ii)):
                        rows.append((f, D0 + ii[j] / W.FS, float(v[j]), float(nm[j]), lab))
        DET = pd.DataFrame(rows, columns=["fid", "t", "net_cc", "net_mad", "series"])
        DET["date"] = ds; DET["window"] = win; DET["kind"] = kind
        FD = pd.DataFrame(fam_rows)
        labs = ["data", "control"] + (["control2"] if C2 else [])
        for lab in labs:
            FD[f"med_{lab}"] = [famstat.get((f, lab), (np.nan,) * 4)[0] for f in FD.fid]
            FD[f"mad_{lab}"] = [famstat.get((f, lab), (np.nan,) * 4)[1] for f in FD.fid]
            FD[f"nvalid_{lab}"] = [famstat.get((f, lab), (np.nan,) * 4)[2] for f in FD.fid]
        for lab in labs:
            for T, nm in [(W.T_REC, "rec"), (W.T_EXT, "ext"), (W.T_CORE, "core")]:
                c = DET[(DET.series == lab) & (DET.net_mad >= T)].groupby("fid").size()
                FD[f"n_{nm}_{lab}"] = FD.fid.map(c).fillna(0).astype(int)
        FD.to_parquet(f"{od}/famday_{ds}.parquet", index=False)
        DET.to_parquet(outp + ".tmp", index=False)
        os.replace(outp + ".tmp", outp)
        json.dump(dict(date=ds, window=win, n_sta=len(D), nodes={k: v["node"] for k, v in D.items()},
                       n_fam_obs=len(FL), n_det=int((DET.series == "data").sum()),
                       n_ctl=int((DET.series == "control").sum()), n_badw=n_badw,
                       secs=time.time() - t_day, dead=dead, shift_site=shift, shift_site2=shift2 if C2 else {}), open(f"{od}/run_{ds}.json", "w"), indent=1)
        del NET, VM
        if not a.no_control:
            del CTL, CM
        if C2:
            del CTL2, CM2
        if dev == "cuda":
            torch.cuda.empty_cache()
        W.log(LOG, f"[mf:{a.tag}] {ds} {kind}: {len(D)} sta, {len(FL)} fam observable, "
              f"det {int((DET.series=='data').sum())} ctl {int((DET.series=='control').sum())} "
              f"(>=8.75: {int(((DET.series=='data')&(DET.net_mad>=W.T_CORE)).sum())} / "
              f"{int(((DET.series=='control')&(DET.net_mad>=W.T_CORE)).sum())}) {time.time()-t_day:.0f}s")


if __name__ == "__main__":
    main()
