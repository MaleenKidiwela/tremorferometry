#!/usr/bin/env python
"""WA pre-production step 4 (2026-10-02): cleaning the scrambled control of real activity, and a time-reversed-template
control for comparison.

Scrambled control (wa_mf.py): each SITE's cc series is rolled by its own shift sh (samples, stored in run_<date>.json
shift_site / shift_site2), CTL[i] = sum_site c * cc_site[i - sh_site]  (torch.roll(x, sh)[i] = x[i - sh]).  So a control
candidate at origin index i saw site q's data at index i - sh_q: real origin time t_real = t - sh_q / 20 Hz (mod the
day).  The task text says "t + shift"; the code's convention is used and the opposite sign is run as a check (it must
give chance-level contamination).

  sitez    per candidate (data and every control series) the per-site cc at the candidate position (controls: the
           site's rolled position, exactly as the MF built it) and its z against 1000 random origins of the same
           site-day (step1_features.py method): dominant site = argmax z, its z, the 2nd site, n_z3 (sites with z >= 3)
  contam   control candidate CONTAMINATED if a DATA candidate of any family (net_mad >= 8, i.e. every stored one) lies
           within +-2 s of it AT THE DOMINANT SITE: |(t_real + off(f, q)) - (t2 + off(f2, q))| <= 2 s, off = the family's
           template offset at that site (w-weighted over the site's stations; only families with a template at q).
           Variants (columns): origin-time match (any family, no site offsets); coherent data only (data candidate with
           >= 2 sites at z >= 3); other family only (a data candidate of a DIFFERENT family: the same family's data
           candidate at t_real is by construction the same single-site burst when that site dominates the network);
           chance level (control arrival moved by +-97 s); opposite shift sign.
  select   60 trial-year days (2010-2014, >= 17 core stations cached, not 2011-04): 20 PNSN-quiet (0 tremor windows in
           every WA region table that day and the days either side) / 20 moderate (1-80 windows, NW < 100) / 20 ETS
           (NW_WA >= 100 windows), fixed seed
  mktpl    time-reversed T1 build templates (each channel reversed in time; same offsets, weights, channels)
  report   tables -> results/wa_pre/ctl/

  python scripts/pre_ctlclean.py sitez --mf results/wa_gonogo/mf_core21 --setup results/wa_gonogo/setup --tpl build --out results/wa_pre/ctl/sitez_core21 --workers 12
"""
from __future__ import annotations
import os, sys, glob, json, time, argparse, warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pre_common as P
W = P.W

OD = f"{P.PRE}/ctl"
NRAND = 1000
_G = {}


def _init(setup, tags, cache_dirs):
    warnings.filterwarnings("ignore")
    os.environ["WA_SETUP"] = setup
    W.SETUP = setup
    P.patch_wa(cache_dirs)
    import wa_mf
    T, Q = wa_mf.load_tpl(tags, all_templates=True)
    S = pd.read_csv(f"{setup}/stations.csv")
    _G.update(T=T, SITE=dict(zip(S.netsta, S.site)))


def guard_np(x, live, m, n=W.NT, floor=1e-3):
    import step1_features as S1
    return S1.guard_np(x, live, m, n, floor)


def sitez_day(args):
    mfdir, mf2, ds, out = args
    t0 = time.time()
    T, SITE = _G["T"], _G["SITE"]
    n = W.NT; NV = W.NDAY - n + 1
    D0 = pd.Timestamp(ds).timestamp()
    J = json.load(open(f"{mfdir}/run_{ds}.json"))
    SH = {"data": {}, "control": J["shift_site"], "control2": J.get("shift_site2", {})}
    C = pd.read_parquet(f"{mfdir}/det_{ds}.parquet")
    parts = [C]
    if mf2 and os.path.exists(f"{mf2}/run_{ds}.json"):
        J2 = json.load(open(f"{mf2}/run_{ds}.json"))
        SH["c3"] = J2["shift_site"]; SH["c4"] = J2.get("shift_site2", {})
        C2 = pd.read_parquet(f"{mf2}/det_{ds}.parquet")
        C2 = C2[C2.series != "data"].copy(); C2["series"] = C2.series.map({"control": "c3", "control2": "c4"})
        parts.append(C2)
    C = pd.concat(parts, ignore_index=True)
    if not len(C):
        pd.DataFrame().to_parquet(out); return ds, 0, 0.0
    C["i"] = np.round((C.t.to_numpy() - D0) * W.FS).astype(np.int64)
    FD = pd.read_parquet(f"{mfdir}/famday_{ds}.parquet")
    FD = FD[FD.observable]
    nodes = J["nodes"]
    fam_sta = {int(f): s.split(",") for f, s in zip(FD.fid, FD.stations)}
    rng = np.random.default_rng(1000 + int(pd.Timestamp(ds).strftime("%Y%m%d")))
    R = rng.integers(0, NV, NRAND)
    Cg = {f: g for f, g in C.groupby("fid")}
    CC, RC, NCH, OFF = {}, {}, {}, {}
    sta_fams = {}
    for f, ss in fam_sta.items():
        for ns in ss:
            sta_fams.setdefault(ns, []).append(f)
    for ns, fl in sta_fams.items():
        p = W.cache_path(ns, ds)
        if not os.path.exists(p):
            continue
        z = np.load(p)
        node = str(z["node"])
        if node != nodes.get(ns):
            continue
        x = z["x"]; m = z["m"]; live = z["live"]
        bad, ew = guard_np(x, live, m)
        efl = 1e-10 * np.median(ew[live])
        site = SITE[ns]
        for f in fl:
            key = (f, ns, node)
            if key not in T or f not in Cg:
                continue
            tp = T[key]
            chi = np.nonzero(tp["live"] & live)[0]
            if not len(chi):
                continue
            d = int(round(tp["off"] * W.FS))
            g = Cg[f]
            sh = np.array([int(SH[s].get(site, 0)) % NV if s != "data" else 0 for s in g.series], np.int64)
            j = ((g.i.to_numpy() - sh) % NV) + d
            allj = np.concatenate([j, R + d])
            okj = (allj >= 0) & (allj < NV)
            jj = np.clip(allj, 0, NV - 1)
            idx = jj[:, None] + np.arange(n)[None]
            num = np.einsum("cnk,ck->cn", x[chi[:, None, None], idx[None]], tp["tpl"][chi])
            E = np.sqrt(np.maximum(ew[chi[:, None], jj[None]], efl)).astype(np.float32)
            cc = np.where(okj & ~bad[jj], (num / E).mean(0), np.nan)
            CC.setdefault(f, {})[ns] = cc[:len(j)]
            RC.setdefault(f, {})[ns] = cc[len(j):]
            NCH.setdefault(f, {})[ns] = len(chi)
            OFF.setdefault(f, {})[ns] = tp["off"]
        del x, ew, bad, z
    out_rows = []
    for f, g in Cg.items():
        if f not in CC:
            continue
        sts = [ns for ns in fam_sta.get(f, []) if ns in CC[f]]
        if not sts:
            continue
        wm = {ns: float(NCH[f][ns]) * T[(f, ns, nodes[ns])]["w"] for ns in sts}
        groups = {}
        for ns in sts:
            groups.setdefault(SITE[ns], []).append(ns)
        sites = sorted(groups)
        SC = np.zeros((len(g), len(sites))); SR = np.zeros((NRAND, len(sites)))
        for k, q in enumerate(sites):
            mem = groups[q]; ws = np.array([wm[s] for s in mem])
            SC[:, k] = (np.stack([CC[f][s] for s in mem], 1) * ws).sum(1) / ws.sum()
            SR[:, k] = (np.stack([RC[f][s] for s in mem], 1) * ws).sum(1) / ws.sum()
        med = np.nanmedian(SR, 0); mad = 1.4826 * np.nanmedian(np.abs(SR - med), 0)
        mad = np.where(mad > 0, mad, np.nan)
        Z = np.nan_to_num((SC - med) / mad, nan=-99.0)
        o = np.argsort(-Z, 1)
        z1 = Z[np.arange(len(g)), o[:, 0]]
        z2 = Z[np.arange(len(g)), o[:, 1]] if len(sites) > 1 else np.full(len(g), np.nan)
        sarr = np.array(sites)
        out_rows.append(pd.DataFrame(dict(
            fid=f, t=g.t.to_numpy(), series=g.series.to_numpy(), net_mad=g.net_mad.to_numpy(), n_sites=len(sites),
            dom_site=sarr[o[:, 0]], dom_z=z1, site2=sarr[o[:, 1]] if len(sites) > 1 else "", z2=z2,
            n_z3=(Z >= 3).sum(1),
            sh_dom=[int(SH[s].get(q, 0)) % NV if s != "data" else 0 for s, q in zip(g.series, sarr[o[:, 0]])])))
    F = pd.concat(out_rows, ignore_index=True) if out_rows else pd.DataFrame()
    F["date"] = ds
    F.to_parquet(out + ".tmp", index=False); os.replace(out + ".tmp", out)
    return ds, len(F), time.time() - t0


def sitez(a):
    os.makedirs(a.out, exist_ok=True)
    days = sorted(os.path.basename(p)[4:14] for p in glob.glob(f"{a.mf}/det_*.parquet"))
    if a.days:
        keep = set(open(a.days).read().split()) if os.path.exists(a.days) else set(a.days.split(","))
        days = [d for d in days if d in keep]
    jobs = [(a.mf, a.mf2, d, f"{a.out}/sitez_{d}.parquet") for d in days if not os.path.exists(f"{a.out}/sitez_{d}.parquet")]
    P.log(f"[sitez] {a.mf}: {len(jobs)} days, {a.workers} workers")
    cdirs = [x for x in a.cache_dirs.split(",") if x]
    with ProcessPoolExecutor(a.workers, initializer=_init, initargs=(os.path.abspath(a.setup), a.tpl.split(","), cdirs)) as ex:
        for fu in as_completed([ex.submit(sitez_day, j) for j in jobs]):
            ds, nr, secs = fu.result()
            P.log(f"[sitez] {ds}: {nr} candidates {secs:.0f}s")


# --------------------------------------------------------------------------------------------- offsets at a site
def site_offsets(setup, tags, mfdir):
    """per (date, fid, site) w-weighted template offset over the site's stations live that day."""
    import wa_mf
    T, Q = wa_mf.load_tpl(tags, all_templates=True)
    S = pd.read_csv(f"{setup}/stations.csv"); SITE = dict(zip(S.netsta, S.site))
    tab = {}
    for (f, ns, node), t in T.items():
        if ns in SITE:
            tab.setdefault((node, f), []).append((ns, t["off"], t["w"]))
    return tab, SITE


def near_any(x, ref, tol=W.MATCH_S):
    if len(ref) == 0 or len(x) == 0:
        return np.zeros(len(x), bool)
    j = np.clip(np.searchsorted(ref, x), 1, max(len(ref) - 1, 1))
    if len(ref) == 1:
        return np.abs(ref[0] - x) <= tol
    return np.minimum(np.abs(ref[j] - x), np.abs(ref[j - 1] - x)) <= tol


def near_other(x, fx, ref, rf, tol=W.MATCH_S):
    """bool: a ref element within tol of x whose family differs from fx (ref sorted, rf aligned with ref)."""
    out = np.zeros(len(x), bool)
    if len(ref) == 0:
        return out
    lo = np.searchsorted(ref, x - tol, "left"); hi = np.searchsorted(ref, x + tol, "right")
    for k in np.flatnonzero(hi > lo):
        out[k] = bool((rf[lo[k]:hi[k]] != fx[k]).any())
    return out


def contam(a):
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    tags = a.tpl.split(",")
    import wa_mf
    T, Q = wa_mf.load_tpl(tags, all_templates=True)
    S = pd.read_csv(f"{a.setup}/stations.csv"); SITE = dict(zip(S.netsta, S.site))
    NV = W.NDAY - W.NT + 1
    rows = []
    for p in sorted(glob.glob(f"{a.sitez}/sitez_*.parquet")):
        ds = os.path.basename(p)[6:16]
        Z = pd.read_parquet(p)
        if not len(Z):
            continue
        J = json.load(open(f"{a.mf}/run_{ds}.json")); nodes = J["nodes"]
        D0 = pd.Timestamp(ds).timestamp()
        # family offset at each site (stations live that day with a template of their node)
        off = {}
        for (f, ns, node), t in T.items():
            if nodes.get(ns) == node and ns in SITE:
                off.setdefault((f, SITE[ns]), []).append((t["off"], t["w"]))
        offq = {k: float(np.average([o for o, _ in v], weights=[w for _, w in v])) for k, v in off.items()}
        Dd = Z[Z.series == "data"]
        # arrival times of data candidates per site
        arr, arr_coh, arr_f = {}, {}, {}
        for q in set(SITE.values()):
            fo = {f: offq[(f, q)] for f in Dd.fid.unique() if (f, q) in offq}
            if not fo:
                continue
            m = Dd.fid.isin(fo)
            at = Dd.t[m].to_numpy() + Dd.fid[m].map(fo).to_numpy()
            oo = np.argsort(at)
            arr[q] = at[oo]; arr_f[q] = Dd.fid[m].to_numpy()[oo]
            mc = m & (Dd.n_z3 >= 2)
            arr_coh[q] = np.sort(Dd.t[mc].to_numpy() + Dd.fid[mc].map(fo).to_numpy())
        tor = np.sort(Dd.t.to_numpy()); tor_coh = np.sort(Dd.t[Dd.n_z3 >= 2].to_numpy())
        Cc = Z[Z.series != "data"].copy()
        if not len(Cc):
            continue
        i = np.round((Cc.t.to_numpy() - D0) * W.FS).astype(np.int64)
        sh = Cc.sh_dom.to_numpy(np.int64)
        t_real = D0 + ((i - sh) % NV) / W.FS
        t_flip = D0 + ((i + sh) % NV) / W.FS               # opposite sign (check)
        o_c = np.array([offq.get((f, q), np.nan) for f, q in zip(Cc.fid, Cc.dom_site)])
        res = {k: np.zeros(len(Cc), bool) for k in ["contam", "contam_coh", "contam_other", "chance_p97", "chance_m97",
                                                    "chance_coh_p97", "chance_other_p97", "flip_sign"]}
        fc = Cc.fid.to_numpy()
        qs = Cc.dom_site.to_numpy()
        for q in np.unique(qs):
            mq = qs == q
            ref = arr.get(q, np.zeros(0)); refc = arr_coh.get(q, np.zeros(0))
            x = t_real[mq] + o_c[mq]
            res["contam"][mq] = near_any(x, ref)
            res["contam_coh"][mq] = near_any(x, refc)
            res["contam_other"][mq] = near_other(x, fc[mq], ref, arr_f.get(q, np.zeros(0)))
            res["chance_coh_p97"][mq] = near_any(x + 97.0, refc)
            res["chance_other_p97"][mq] = near_other(x + 97.0, fc[mq], ref, arr_f.get(q, np.zeros(0)))
            res["chance_p97"][mq] = near_any(x + 97.0, ref)
            res["chance_m97"][mq] = near_any(x - 97.0, ref)
            res["flip_sign"][mq] = near_any(t_flip[mq] + o_c[mq], ref)
        for k, v in res.items():
            Cc[k] = v
        Cc["contam_origin"] = near_any(t_real, tor)
        Cc["contam_origin_coh"] = near_any(t_real, tor_coh)
        Cc["off_known"] = np.isfinite(o_c)
        Cc["t_real"] = t_real
        rows.append(Cc.drop(columns=[c for c in ["site2"] if c in Cc]))
    X = pd.concat(rows, ignore_index=True)
    X.to_parquet(a.out, index=False)
    P.log(f"[contam] {a.sitez}: {len(X)} control candidates -> {a.out}")


# --------------------------------------------------------------------------------------------- day selection / templates
def select(a):
    os.makedirs(OD, exist_ok=True)
    out = f"{OD}/days60.csv"
    if os.path.exists(out):
        print(open(out).read()[:500]); return
    import datetime
    from collections import Counter
    cnt = Counter()
    for ns in P.core_list():
        for y in range(2010, 2015):
            d = f"{P.WA_CACHE}/{ns}/{y}"
            if os.path.isdir(d):
                for f in os.listdir(d):
                    if f.endswith(".npz") and f[:3].isdigit():
                        cnt[(datetime.date(y, 1, 1) + datetime.timedelta(int(f[:3]) - 1)).isoformat()] += 1
    t = pd.DataFrame({"date": sorted(cnt)}); t["n_core_cached"] = t.date.map(cnt)
    td = P.tremor_days()
    for r in td:
        t[r] = t.date.map(td[r]).fillna(0).astype(int)
    t["mx"] = t[["NW_WA", "SW_WA", "OLY2010"]].max(axis=1)
    nb = t.date.map(lambda d: max(max(td[r].get((pd.Timestamp(d) + pd.Timedelta(days=k)).strftime("%Y-%m-%d"), 0) for r in td) for k in (-1, 1)))
    t["mx_nb"] = nb
    t = t[(t.n_core_cached >= 17) & (t.date.str[:7] != "2011-04")]
    cls = np.where((t.mx == 0) & (t.mx_nb == 0), "quiet", np.where(t.NW_WA >= 100, "ets",
                   np.where((t.mx >= 1) & (t.mx <= 80) & (t.NW_WA < 100), "moderate", "")))
    t["cls"] = cls
    rng = np.random.default_rng(20261002)
    pick = []
    for c in ["quiet", "moderate", "ets"]:
        g = t[t.cls == c]
        pick.append(g.iloc[np.sort(rng.choice(len(g), min(20, len(g)), replace=False))])
    Pk = pd.concat(pick).sort_values("date")
    Pk.to_csv(out, index=False)
    t.to_csv(f"{OD}/days_candidates.csv", index=False)
    print(Pk.groupby("cls").agg(n=("date", "size"), first=("date", "min"), last=("date", "max"), med_nw=("NW_WA", "median")).to_string())
    print(Pk.date.str[:4].value_counts().sort_index().to_string())


def mktpl(a):
    d = f"{OD}/res/templates"
    os.makedirs(d, exist_ok=True)
    src = f"{P.STUDY}/results/templates"
    for f in ["tpl_build.npz", "tpl_qc_build.csv"]:
        if not os.path.exists(f"{d}/{f}"):
            os.symlink(f"{src}/{f}", f"{d}/{f}")
    if not os.path.exists(f"{d}/tpl_buildrev.npz"):
        z = np.load(f"{src}/tpl_build.npz")
        np.savez(f"{d}/tpl_buildrev.npz", keys=z["keys"], tpl=np.ascontiguousarray(z["tpl"][:, :, ::-1]))
        import shutil
        shutil.copyfile(f"{src}/tpl_qc_build.csv", f"{d}/tpl_qc_buildrev.csv")
    z = np.load(f"{d}/tpl_buildrev.npz"); z0 = np.load(f"{src}/tpl_build.npz")
    assert np.array_equal(z["tpl"][:, :, 0], z0["tpl"][:, :, -1]) and (z["keys"] == z0["keys"]).all()
    assert open(f"{src}/tpl_qc_build.csv", "rb").read() == open(f"{d}/tpl_qc_buildrev.csv", "rb").read(), "qc copy differs"
    print("reversed templates ok", z["tpl"].shape)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("what", choices=["sitez", "contam", "select", "mktpl"])
    ap.add_argument("--mf", default=""); ap.add_argument("--mf2", default="")
    ap.add_argument("--setup", default=f"{P.STUDY}/results/wa_gonogo/setup")
    ap.add_argument("--tpl", default="build")
    ap.add_argument("--out", default=""); ap.add_argument("--sitez", default="")
    ap.add_argument("--days", default=""); ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--cache-dirs", default="")
    a = ap.parse_args()
    if a.what in ("sitez", "contam"):
        os.environ["WA_SETUP"] = os.path.abspath(a.setup); W.SETUP = os.path.abspath(a.setup); P.patch_wa()
    dict(sitez=sitez, contam=contam, select=select, mktpl=mktpl)[a.what](a)
