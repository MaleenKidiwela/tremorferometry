#!/usr/bin/env python
"""STEP 2 (Merlin, 2026-09-30) part 3 — WHICH MEMBER FIRED: attribution of group detections to member families.

stage cc  (GPU) — for every multi-member group, every record-tier (>= 8.0 MAD) candidate of the group matched
  filter (results/mf_step2grp; series data and control = the first site-rolled scrambled control), the normalised
  cc of each member's template at the stations of the NEAREST sites only:
    member template  12 s cut from the member's OWN build-half long stack (results/templates/stacks_build.npz,
                     >= 30 events at that station node) at the group template's window mapped into the member frame
                     (start = group offset + tau_m,s), each live channel L2-normalised
    sites            sites (3-km merge B; B013+B943 is one site, PGC+B009-B011 one, LRIV+B005-B007 one) ordered by
                     distance from the group centroid; the first K_MAX = 3 with a station live that day (dead-station
                     rule of the MF) where >= 2 members have a template; every such station of a chosen site is used
    eligible member  a template at EVERY chosen station of that group-day (others: no attribution that day)
    cc               wa_mf chain (1-8 Hz cache, robust energy floor + mask guard, per-channel normalised cc, mean of
                     channels live in both); read at the candidate's aligned position (scrambled control: the site's
                     rolled position, as the MF built it) and at lags -10..+10 samples (+-0.5 s)
    per (candidate, member):  s0_k{1,2,3} = sum of the zero-lag cc over the stations of the first k chosen sites,
                     smax_k3 = sum of the max over lags, mlag_k3 = mean |best lag| (sub-sample, parabolic) and
                     lag_sd_k3 = SD of the best lags over stations (the moveout / lag-difference signature)
stage eval — P(member | group event) = softmax(beta * s0_k3) over the eligible members, beta fitted by maximum
  likelihood on the BUILD half (V2008 + V2010, half 0) candidates matched (+-2 s) to exactly one eligible member's
  array detection; evaluated on the V2010 / V2008 TEST half and A2011: top-1 accuracy vs chance (1 / n eligible) and vs
  "always the most active member", mean P(true), log-loss skill vs uniform, reliability of P.  Variants (columns):
  k = 1, 2 sites, smax, a 2-feature conditional logit (s0_k3, mlag_k3), prior-only and cc + log prior.  Combined with
  the group-level P(real) (STEP 1 lfdr refitted on group candidates, results/step2_groups/step1/candidates.parquet):
  per member family and day  E = sum over the group's candidates of P(real) * P(member).

  python scripts/step2_attrib.py --stage cc [--days 2010-08-20]
  python scripts/step2_attrib.py --stage eval
  -> results/step2_groups/attrib/{cc_<date>.parquet, chosen_<date>.parquet}, results/step2_groups/attrib_*.csv
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_common as W

G = f"{W.RES}/step2_groups"
SET = f"{G}/setup"
AD = os.environ.get("STEP2_ATTRIB_DIR", f"{G}/attrib")   # smoke: results/step2_groups/smoke/attrib
LOG = f"{W.STUDY}/logs/step2.log"
MF = os.environ.get("STEP2_MF", "step2grp")
TAG = "step2grp"
K_MAX = 3
LAG = 10
NL = int(W.LONG_LEN * W.FS)
IPRE = int(W.LONG_PRE * W.FS)
WINS = ["V2008", "V2010", "A2011"]
OUTD = os.environ.get("STEP2_OUT", G)                    # smoke: results/step2_groups/smoke


def load_static():
    GR = pd.read_csv(f"{SET}/families.csv")
    M = pd.read_csv(f"{G}/members.csv")
    ST = pd.read_parquet(f"{G}/align_station.parquet")
    tau = {(int(f), ns, nd): float(t) for f, ns, nd, t in zip(ST.fid, ST.netsta, ST.node, ST.tau_s)}
    dl = dict(zip(M.fid, M.delta_s))
    Qg = pd.read_csv(f"{W.RES}/templates/tpl_qc_{TAG}.csv")
    Qg = Qg[(Qg.n_events >= W.MIN_EVENTS) & Qg.offset_s.notna()]
    offg = {(int(g), ns, nd): float(o) for g, ns, nd, o in zip(Qg.fid, Qg.netsta, Qg.node, Qg.offset_s)}
    Qb = pd.read_csv(f"{W.RES}/templates/tpl_qc_build.csv", dtype={"live": str})
    Qb["live"] = Qb.live.astype(str).str.zfill(3)
    qb = {(int(f), ns, nd): (int(n), l) for f, ns, nd, n, l in zip(Qb.fid, Qb.netsta, Qb.node, Qb.n_events, Qb.live)}
    Z = np.load(f"{W.RES}/templates/stacks_build.npz", allow_pickle=False)
    kidx = {}
    for i, k in enumerate(Z["keys"]):
        f, ns, nd = str(k).split("|")
        kidx[(int(f), ns, nd)] = i
    STK = Z["stacks"]
    S = pd.read_csv(f"{SET}/stations.csv")
    site_xy = S.groupby("site")[["lat", "lon"]].mean()
    SITE = dict(zip(S.netsta, S.site))
    # attribution templates
    ATT = {}
    for r in GR[GR.n_members > 1].itertuples():
        mem = [int(x) for x in r.members.split(";")]
        for (g, ns, nd), og in offg.items():
            if g != r.fid:
                continue
            for m in mem:
                n_l = qb.get((m, ns, nd))
                if n_l is None or n_l[0] < W.MIN_EVENTS or (m, ns, nd) not in kidx:
                    continue
                t_ = tau.get((m, ns, nd), dl[m])
                j0 = IPRE + int(round((og + t_) * W.FS))
                if j0 < 0 or j0 + W.NT > NL:
                    continue
                st = STK[kidx[(m, ns, nd)]].astype(np.float64)
                tp = st[:, j0:j0 + W.NT].copy()
                live = np.array([c == "1" for c in n_l[1]])
                nrm = np.sqrt((tp ** 2).sum(1))
                live &= nrm > 0
                if not live.any():
                    continue
                tp[~live] = 0.0
                tp[live] /= nrm[live, None]
                ATT[(r.fid, m, ns, nd)] = dict(tpl=tp.astype(np.float32), live=live, d=int(round(og * W.FS)))
    return GR, M, ATT, site_xy, SITE, S


def stage_cc(days_arg):
    import torch
    import torch.nn.functional as Fnn
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    t0 = time.time()
    os.makedirs(AD, exist_ok=True)
    GR, M, ATT, site_xy, SITE, S = load_static()
    sites_of = S.groupby("site").netsta.apply(list).to_dict()
    multi = GR[GR.n_members > 1]
    site_order = {}
    for r in multi.itertuples():
        d = W.km(site_xy.lat.to_numpy(), site_xy.lon.to_numpy(), r.lat_use, r.lon_use)
        site_order[r.fid] = [(site_xy.index[i], float(d[i])) for i in np.argsort(d)]
    att_by = {}
    for (g, m, ns, nd), v in ATT.items():
        att_by.setdefault((g, ns, nd), {})[m] = v
    W.log(LOG, f"[attrib:cc] {len(multi)} multi-member groups, {len(ATT)} member attribution templates, device {dev}")
    n = W.NT
    NV = W.NDAY - n + 1
    NF = int(2 ** np.ceil(np.log2(W.NDAY + n)))

    def guard(X, live, m, floor=1e-3):       # == wa_mf.guard
        xd = X.double()
        EN = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64), xd * xd], 1), 1)
        ew = (EN[:, n:] - EN[:, :-n])[:, :NV]
        del EN, xd
        lv = torch.from_numpy(live).to(dev)
        ewl = ew[lv]
        med_ch = ewl.median(dim=1).values
        e3 = ewl.sum(0); med3 = e3.median()
        ref_ch = torch.maximum(med_ch, 1e-2 * med3)
        bad = ((ewl < floor * ref_ch[:, None]).any(0) | (e3 < floor * med3) | (ewl.min(0).values <= 0))
        if float(med3) <= 0:
            bad = torch.ones_like(bad)
        mk = torch.from_numpy(~m).to(dev)
        cmk = torch.cat([torch.zeros(1, device=dev), torch.cumsum(mk.float(), 0)])
        bad |= (cmk[n:n + NV] - cmk[:NV]) > 0
        if bool(bad.any()):
            bad = Fnn.max_pool1d(bad.float()[None, None], kernel_size=2 * n + 1, stride=1, padding=n)[0, 0] > 0
        efl = 1e-10 * ewl.median()
        E = ew.clamp_min(efl).sqrt().float()
        return bad, E

    days = [d for w in WINS for d, _ in W.window_days(w)]
    if days_arg:
        days = days_arg
    lags = torch.arange(-LAG, LAG + 1, device=dev)
    for ds in days:
        outp = f"{AD}/cc_{ds}.parquet"
        if os.path.exists(outp) or not os.path.exists(f"{W.RES}/mf_{MF}/det_{ds}.parquet"):
            continue
        td = time.time()
        D0 = pd.Timestamp(ds).timestamp()
        J = json.load(open(f"{W.RES}/mf_{MF}/run_{ds}.json"))
        nodes = J["nodes"]; sh1 = J["shift_site"]
        DE = pd.read_parquet(f"{W.RES}/mf_{MF}/det_{ds}.parquet")
        DE = DE[DE.series.isin(["data", "control"]) & DE.fid.isin(set(multi.fid))].copy()
        DE["i"] = np.round((DE.t.to_numpy() - D0) * W.FS).astype(np.int64)
        # chosen stations per group-day
        chosen, need = {}, {}
        crow = []
        for r in multi.itertuples():
            if not (DE.fid == r.fid).any():
                continue
            picks = []
            for site, dkm in site_order[r.fid]:
                sts = [ns for ns in sites_of[site] if ns in nodes and len(att_by.get((r.fid, ns, nodes[ns]), {})) >= 2]
                if sts:
                    picks.append((site, dkm, sts))
                if len(picks) == K_MAX:
                    break
            if not picks:
                continue
            allst = [ns for _, _, sts in picks for ns in sts]
            elig = set.intersection(*[set(att_by[(r.fid, ns, nodes[ns])]) for ns in allst])
            if len(elig) < 2:
                crow.append(dict(gid=r.fid, date=ds, n_elig=len(elig), sites="|".join(p[0] for p in picks), stations=",".join(allst)))
                continue
            chosen[r.fid] = (picks, sorted(elig))
            for k, (site, dkm, sts) in enumerate(picks):
                for ns in sts:
                    need.setdefault(ns, []).append((r.fid, k))
            crow.append(dict(gid=r.fid, date=ds, n_elig=len(elig), elig=";".join(map(str, sorted(elig))),
                             sites="|".join(p[0] for p in picks), site_km="|".join(f"{p[1]:.1f}" for p in picks),
                             stations=",".join(allst)))
        # per (gid, station): cc arrays (members x candidates x lags)
        RES = {}
        Cg = {g: x for g, x in DE.groupby("fid")}
        for ns, lst in need.items():
            p = W.cache_path(ns, ds)
            z = np.load(p)
            node = str(z["node"])
            if node != nodes[ns]:
                continue
            X = torch.from_numpy(z["x"]).to(dev)
            live_d = z["live"]
            bad, E = guard(X, live_d, z["m"])
            DH = torch.fft.rfft(X, n=NF)
            del X
            for g, k in lst:
                elig = chosen[g][1]
                A = [ATT[(g, m, ns, node)] for m in elig]
                d = A[0]["d"]
                cand = Cg[g]
                site = SITE[ns]
                shv = np.where(cand.series.to_numpy() == "control", int(sh1.get(site, 0)) % NV, 0)
                pos = ((cand.i.to_numpy() - shv) % NV) + d
                pt = torch.from_numpy(pos).to(dev)
                idx = pt[:, None] + lags[None]                          # (nc, 21)
                okp = (idx >= 0) & (idx < NV)
                idc = idx.clamp(0, NV - 1)
                out = torch.full((len(elig), len(cand), 2 * LAG + 1), float("nan"), device=dev)
                for j, a in enumerate(A):
                    ch = a["live"] & live_d
                    chi = [c for c in range(3) if ch[c]]
                    if not chi:
                        continue
                    WH = torch.fft.rfft(torch.from_numpy(a["tpl"][chi]).to(dev), n=NF)
                    num = torch.fft.irfft(WH.conj() * DH[chi], n=NF)[:, :NV]
                    cc = (num / E[chi]).mean(0)
                    cc[bad] = float("nan")
                    v = cc[idc]
                    v[~okp] = float("nan")
                    out[j] = v
                    del WH, num, cc
                RES[(g, ns)] = (k, out.cpu().numpy())
            del DH, E, bad
        # aggregate
        rows = []
        for g, (picks, elig) in chosen.items():
            cand = Cg[g]
            nc = len(cand)
            arrs = [(k, RES[(g, ns)][1]) for k, (site, dkm, sts) in enumerate(picks) for ns in sts if (g, ns) in RES]
            if not arrs:
                continue
            ne = len(elig)
            cc0 = np.stack([a[:, :, LAG] for _, a in arrs], 2)                 # (ne, nc, nsta)
            ks = np.array([k for k, _ in arrs])
            full = np.stack([a for _, a in arrs], 2)                             # (ne, nc, nsta, 21)
            with np.errstate(invalid="ignore"):
                fz = np.where(np.isfinite(full), full, -9.0)
                kb = fz.argmax(3)
                cmax = np.take_along_axis(fz, kb[..., None], 3)[..., 0]
                cmax = np.where(cmax <= -8.9, np.nan, cmax)
                ym = np.take_along_axis(fz, np.clip(kb - 1, 0, 2 * LAG)[..., None], 3)[..., 0]
                yp = np.take_along_axis(fz, np.clip(kb + 1, 0, 2 * LAG)[..., None], 3)[..., 0]
                den = ym - 2 * cmax + yp
                frac = np.where((kb > 0) & (kb < 2 * LAG) & (den < 0), 0.5 * (ym - yp) / den, 0.0)
                lag = (kb - LAG + frac) / W.FS
                lag = np.where(np.isfinite(cmax), lag, np.nan)
            rec = dict(gid=np.repeat(g, ne * nc), fid=np.repeat(elig, nc), t=np.tile(cand.t.to_numpy(), ne),
                       series=np.tile(cand.series.to_numpy(), ne), net_mad=np.tile(cand.net_mad.to_numpy(), ne),
                       n_elig=ne)
            for kk in range(K_MAX):
                sel = ks <= kk
                rec[f"s0_k{kk+1}"] = np.nansum(cc0[:, :, sel], 2).ravel() if sel.any() else np.nan
                rec[f"n_sta_k{kk+1}"] = np.isfinite(cc0[:, :, sel]).sum(2).ravel() if sel.any() else 0
            rec["smax_k3"] = np.nansum(cmax, 2).ravel()
            with np.errstate(invalid="ignore"):
                rec["mlag_k3"] = np.nanmean(np.abs(lag), 2).ravel()
                rec["lag_sd_k3"] = np.nanstd(lag, 2).ravel()
                rec["lag_k1"] = np.nanmean(lag[:, :, ks == 0], 2).ravel() if (ks == 0).any() else np.nan
            rows.append(pd.DataFrame(rec))
        R = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
        R["date"] = ds
        pd.DataFrame(crow).to_parquet(f"{AD}/chosen_{ds}.parquet", index=False)
        R.to_parquet(outp + ".tmp", index=False)
        os.replace(outp + ".tmp", outp)
        if dev == "cuda":
            torch.cuda.empty_cache()
        W.log(LOG, f"[attrib:cc] {ds}: {len(chosen)} groups attributed, {len(R)} (candidate, member) rows, {time.time()-td:.0f}s")
    W.log(LOG, f"[attrib:cc] done {time.time()-t0:.0f}s")


# ------------------------------------------------------------------ eval
def near_idx(x, ref, tol=W.MATCH_S):
    """index of the nearest ref within tol (-1 if none); ref sorted"""
    x = np.asarray(x, float)
    out = np.full(len(x), -1)
    if len(ref) == 0 or len(x) == 0:
        return out
    j = np.clip(np.searchsorted(ref, x), 0, len(ref) - 1)
    jm = np.clip(j - 1, 0, len(ref) - 1)
    pick = np.where(np.abs(ref[jm] - x) < np.abs(ref[j] - x), jm, j)
    ok = np.abs(ref[pick] - x) <= tol
    out[ok] = pick[ok]
    return out


def softmax_by(logit, key):
    """softmax of logit within groups of key (candidate id); returns P"""
    df = pd.DataFrame(dict(k=key, z=logit))
    mx = df.groupby("k").z.transform("max")
    e = np.exp(df.z - mx)
    return (e / e.groupby(df.k).transform("sum")).to_numpy()


def nll(beta, X, key, y):
    z = X @ np.atleast_1d(beta)
    P = softmax_by(z, key)
    return -np.log(np.clip(P[y], 1e-12, 1)).sum()


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d; h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def stage_eval():
    from scipy.optimize import minimize, minimize_scalar
    from scipy.stats import spearmanr
    t0 = time.time()
    GR = pd.read_csv(f"{SET}/families.csv")
    M = pd.read_csv(f"{G}/members.csv")
    A = pd.concat([pd.read_parquet(p) for p in sorted(glob.glob(f"{AD}/cc_*.parquet"))], ignore_index=True)
    A = A[np.isfinite(A.s0_k3)]
    wk = {d: w for w in WINS for d, _ in W.window_days(w)}
    A["window"] = A.date.map(wk)
    A["half"] = W.half_of(A.t.to_numpy())
    A["cid"] = A.groupby(["gid", "series", "t"]).ngroup()
    nb = dict(zip(M.fid, M.n_build))
    A["n_build"] = A.fid.map(nb).clip(lower=1)
    A["logprior"] = np.log(A.n_build) - A.groupby("cid").n_build.transform(lambda v: np.log(v.sum()))
    # truth (group frame) + scan
    TG = pd.read_parquet(f"{SET}/truth.parquet")
    scan = pd.read_csv(f"{SET}/scan_days.csv")
    sd = {s: set(g.date) for s, g in scan.groupby("det_source")}
    vw = dict(zip(GR.fid, GR.val_windows))
    A["truth_ok"] = [(d in sd.get(f"G{g}", set())) and (w in vw[g].split(",")) for g, d, w in zip(A.gid, A.date, A.window)]
    C = A.drop_duplicates("cid")[["cid", "gid", "t", "series", "date", "window", "half", "truth_ok", "net_mad", "n_elig"]].copy()
    C["tmem"] = ""
    Tg = {g: x.sort_values("t") for g, x in TG.groupby("fid")}
    for g, cg in C[C.truth_ok & (C.series == "data")].groupby("gid"):
        tg = Tg.get(g)
        if tg is None:
            continue
        j = near_idx(cg.t.to_numpy(), tg.t.to_numpy())
        mem = tg.members.to_numpy()
        C.loc[cg.index, "tmem"] = [mem[x] if x >= 0 else "" for x in j]
    A = A.merge(C[["cid", "tmem"]], on="cid", how="left")
    A["is_true"] = [str(f) in s.split(";") if s else False for f, s in zip(A.fid, A.tmem)]
    ntrue = A.groupby("cid").is_true.transform("sum")
    A["matched"] = A.tmem.fillna("") != ""
    A["unique"] = A.matched & (ntrue == 1) & (A.tmem.str.count(";") == 0)
    A["multi_true"] = A.matched & (ntrue >= 1) & ~A.unique
    A["eval_set"] = np.select([A.window.isin(["V2008", "V2010"]) & (A.half == 0), A.window.isin(["V2008", "V2010"]) & (A.half == 1),
                               A.window == "A2011"], ["fit", A.window + "_test", "A2011"], "guard")
    # ---- models
    MODELS = {"cc_k3": ["s0_k3"], "cc_k2": ["s0_k2"], "cc_k1": ["s0_k1"], "ccmax_k3": ["smax_k3"],
              "cc_lag_k3": ["s0_k3", "mlag_k3"], "cc_k3_prior": ["s0_k3", "logprior"]}
    fit = A[(A.eval_set == "fit") & (A.series == "data") & A.unique]
    fitc = fit[fit.cid.isin(fit[fit.is_true].cid)]
    params = {}
    for name, cols in MODELS.items():
        X = np.nan_to_num(fitc[cols].to_numpy(float))
        key = fitc.cid.to_numpy(); y = fitc.is_true.to_numpy()
        if len(cols) == 1:
            r = minimize_scalar(lambda b: nll(b, X, key, y), bounds=(0, 500), method="bounded")
            params[name] = [float(r.x)]
        else:
            r = minimize(lambda b: nll(b, X, key, y), x0=np.array([10.0, 0.0]), method="Nelder-Mead")
            params[name] = [float(v) for v in r.x]
    params["prior_only"] = [1.0]
    MODELS["prior_only"] = ["logprior"]
    for name, cols in MODELS.items():
        z = np.nan_to_num(A[cols].to_numpy(float)) @ np.array(params[name])
        A[f"P_{name}"] = softmax_by(z, A.cid.to_numpy())
    A["P_member"] = A.P_cc_k3
    A["P_unif"] = 1.0 / A.n_elig
    json.dump(dict(params=params, n_fit_candidates=int(fitc.cid.nunique()), models=MODELS), open(f"{OUTD}/attrib_model.json", "w"), indent=1)
    W.log(LOG, f"[attrib:eval] fitted {params} on {fitc.cid.nunique()} build-half candidates ({time.time()-t0:.0f}s)")
    # ---- accuracy tables
    reg = dict(zip(GR.fid, GR.region)); nmem = dict(zip(GR.fid, GR.n_members))
    A["region"] = A.gid.map(reg)
    A["size_class"] = A.gid.map(nmem).map(lambda k: "2" if k == 2 else ("3-6" if k <= 6 else "23"))
    rows = []

    def summarize(sub, label):
        out = []
        for name in list(MODELS) + ["unif"]:
            pc = f"P_{name}"
            d = sub.copy()
            d["rk"] = d.groupby("cid")[pc].rank(ascending=False, method="first")
            top = d[d.rk == 1]
            tru = d[d.is_true]
            ncand = top.cid.nunique()
            k = int(top.is_true.sum())
            lo, hi = wilson(k, ncand)
            chance = float((1 / top.n_elig).mean()) if ncand else np.nan
            ll = -np.log(np.clip(tru[pc], 1e-12, 1)).sum(); llu = -np.log(1 / tru.n_elig).sum()
            # "most active member" baseline
            d["rkp"] = d.groupby("cid").n_build.rank(ascending=False, method="first")
            kp = int(d[(d.rkp == 1)].is_true.sum())
            out.append(dict(label=label, model=name, n_cand=ncand, top1=k / ncand if ncand else np.nan, top1_lo=lo, top1_hi=hi,
                            chance=chance, top1_over_chance=(k / ncand) / chance if ncand else np.nan,
                            most_active_top1=kp / ncand if ncand else np.nan,
                            mean_P_true=float(tru[pc].mean()) if len(tru) else np.nan,
                            mean_P_unif=float((1 / tru.n_elig).mean()) if len(tru) else np.nan,
                            logloss_skill_vs_unif=1 - ll / llu if llu > 0 else np.nan,
                            mean_n_elig=float(top.n_elig.mean()) if ncand else np.nan))
        return out

    D = A[(A.series == "data") & A.unique & A.cid.isin(A[A.is_true].cid)]
    for es in ["fit", "V2010_test", "V2008_test", "A2011"]:
        de = D[D.eval_set == es]
        rows += summarize(de, f"{es}|ALL|ALL|all")
        for rg, x in de.groupby("region"):
            rows += summarize(x, f"{es}|{rg}|ALL|all")
        for sc, x in de.groupby("size_class"):
            rows += summarize(x, f"{es}|ALL|{sc}|all")
        rows += summarize(de[de.net_mad >= W.T_CORE], f"{es}|ALL|ALL|core")
    ACC = pd.DataFrame(rows)
    ACC[["eval_set", "region", "size_class", "tier"]] = ACC.label.str.split("|", expand=True)
    ACC.to_csv(f"{OUTD}/attrib_accuracy.csv", index=False)
    # multi-true candidates (co-detected events): top-1 in the true set vs chance |set| / n
    rows = []
    MT = A[(A.series == "data") & A.multi_true]
    for es, x in MT.groupby("eval_set"):
        x = x.copy()
        x["rk"] = x.groupby("cid").P_member.rank(ascending=False, method="first")
        top = x[x.rk == 1]
        nt = x.groupby("cid").is_true.sum()
        rows.append(dict(eval_set=es, n_cand=top.cid.nunique(), top1_in_set=float(top.is_true.mean()),
                         chance=float((nt / x.groupby("cid").n_elig.first()).mean())))
    pd.DataFrame(rows).to_csv(f"{OUTD}/attrib_multi_true.csv", index=False)
    # per group (test sets)
    rows = []
    Dt = D[D.eval_set.isin(["V2010_test", "V2008_test", "A2011"])].copy()
    Dt["rk"] = Dt.groupby("cid").P_member.rank(ascending=False, method="first")
    for g, x in Dt.groupby("gid"):
        top = x[x.rk == 1]
        rows.append(dict(gid=g, region=reg[g], n_members=nmem[g], n_cand=top.cid.nunique(), top1=float(top.is_true.mean()),
                         chance=float((1 / top.n_elig).mean()), mean_P_true=float(x[x.is_true].P_member.mean()),
                         true_member_counts=";".join(f"{k}:{v}" for k, v in x[x.is_true].fid.value_counts().items())))
    pd.DataFrame(rows).sort_values("n_cand", ascending=False).to_csv(f"{OUTD}/attrib_by_group.csv", index=False)
    # calibration (pairs)
    rows = []
    bins = np.linspace(0, 1, 11)
    for es in ["fit", "V2010_test", "V2008_test", "A2011"]:
        x = D[D.eval_set == es]
        b = np.clip(np.digitize(x.P_member, bins[1:-1]), 0, 9)
        for k, gb in x.groupby(b):
            rows.append(dict(eval_set=es, pbin=f"{bins[k]:.1f}-{bins[k+1]:.1f}", n=len(gb), mean_P=gb.P_member.mean(),
                             frac_true=gb.is_true.mean()))
    CAL = pd.DataFrame(rows)
    CAL.to_csv(f"{OUTD}/attrib_calibration.csv", index=False)
    # scrambled control vs data: how peaked is P
    rows = []
    for es in ["V2010_test", "V2008_test", "A2011"]:
        for s in ["data", "control"]:
            x = A[(A.eval_set == es) & (A.series == s)]
            if not len(x):
                continue
            pm = x.groupby("cid").agg(pmax=("P_member", "max"), n=("n_elig", "first"))
            ent = x.assign(h=-x.P_member * np.log(np.clip(x.P_member, 1e-12, 1))).groupby("cid").h.sum()
            rows.append(dict(eval_set=es, series=s, n_cand=len(pm), mean_Pmax=float(pm.pmax.mean()),
                             mean_Pmax_unif=float((1 / pm.n).mean()),
                             mean_norm_entropy=float((ent / np.log(pm.n)).mean())))
    pd.DataFrame(rows).to_csv(f"{OUTD}/attrib_peakedness.csv", index=False)
    # ---- combine with P(real) (group lfdr refit)
    pcan = os.environ.get("STEP2_PREAL", f"{G}/step1/candidates.parquet")
    if os.path.exists(pcan):
        P = pd.read_parquet(pcan, columns=["fid", "t", "series", "P", "date", "window", "half"])
        P = P[P.series == "data"].rename(columns={"fid": "gid"})
        P["scored"] = (P.half == 1) | (P.window == "A2011")
        att = A[A.series == "data"][["gid", "t", "fid", "P_member", "P_unif"]]
        X = P.merge(att, on=["gid", "t"], how="left")
        memlist = {int(g): [int(v) for v in s.split(";")] for g, s in zip(GR.fid, GR.members)}
        # singletons: P_member = 1; multi-member without attribution that day: uniform over ALL members (flag)
        miss = X.fid.isna()
        fill = []
        for r in X[miss][["gid", "t", "P", "date", "window", "half", "scored"]].itertuples(index=False):
            mm = memlist[int(r.gid)]
            for m in mm:
                fill.append((r.gid, r.t, "data", r.P, r.date, r.window, r.half, r.scored, m, 1.0 / len(mm), 1.0 / len(mm),
                             "singleton" if len(mm) == 1 else "uniform_fallback"))
        Fm = pd.DataFrame(fill, columns=["gid", "t", "series", "P", "date", "window", "half", "scored", "fid", "P_member", "P_unif", "attr_source"])
        X = X[~miss].copy(); X["attr_source"] = "softmax"
        X = pd.concat([X, Fm], ignore_index=True)
        X["fid"] = X.fid.astype(int)
        X["E"] = X.P * X.P_member; X["E_unif"] = X.P * X.P_unif
        X.to_parquet(f"{OUTD}/attrib_expected_candidates.parquet", index=False)
        E = X[X.scored].groupby(["fid", "date"]).agg(E=("E", "sum"), E_unif=("E_unif", "sum"),
                                                     n_softmax=("attr_source", lambda v: int((v == "softmax").sum()))).reset_index()
        # truth per member family-day (member frame, scored half), T1 step-1 sum P for comparison
        T1 = pd.read_parquet(f"{W.RES}/merlin_tests/t1_setup/truth.parquet")
        T1 = T1[(T1.half == 1) | (T1.window == "A2011")]
        own = T1.groupby(["fid", "date"]).size().rename("truth_own").reset_index()
        R1 = pd.read_parquet(f"{W.RES}/step1_prob/rate_family_day.parquet")
        R1 = R1[R1.truth_ok & (R1.kind == "episode")][["fid", "date", "window", "region", "sumP_scored_half"]]
        Y = R1.merge(E, on=["fid", "date"], how="left").merge(own, on=["fid", "date"], how="left")
        for c in ["E", "E_unif", "truth_own", "n_softmax"]:
            Y[c] = Y[c].fillna(0)
        g_of = dict(zip(M.fid, M.gid))
        Y["gid"] = Y.fid.map(g_of); Y["n_members"] = Y.gid.map(nmem)
        Y.to_parquet(f"{OUTD}/attrib_expected_family_day.parquet", index=False)
        rows = []
        for (w, rg), x in [(("ALL", "ALL"), Y)] + [((w, "ALL"), x) for w, x in Y.groupby("window")] + \
                [(("ALL", rg), x) for rg, x in Y.groupby("region")]:
            for lab, xx in [("all families", x), ("members of multi-member groups", x[x.n_members > 1])]:
                if len(xx) < 10:
                    continue
                rows.append(dict(window=w, region=rg, families=lab, n_famdays=len(xx), truth=int(xx.truth_own.sum()),
                                 sum_E=float(xx.E.sum()), sum_E_unif=float(xx.E_unif.sum()), sum_P_T1=float(xx.sumP_scored_half.sum()),
                                 spearman_E=spearmanr(xx.E, xx.truth_own).statistic,
                                 spearman_E_unif=spearmanr(xx.E_unif, xx.truth_own).statistic,
                                 spearman_P_T1=spearmanr(xx.sumP_scored_half, xx.truth_own).statistic))
        pd.DataFrame(rows).to_csv(f"{OUTD}/attrib_expected_count.csv", index=False)
    A.drop(columns=["tmem"]).to_parquet(f"{OUTD}/attrib_pairs.parquet", index=False)
    W.log(LOG, f"[attrib:eval] done {time.time()-t0:.0f}s")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["cc", "eval"], required=True)
    ap.add_argument("--days", default="")
    a = ap.parse_args()
    if a.stage == "cc":
        stage_cc(a.days.split(",") if a.days else None)
    else:
        stage_eval()


if __name__ == "__main__":
    main()
