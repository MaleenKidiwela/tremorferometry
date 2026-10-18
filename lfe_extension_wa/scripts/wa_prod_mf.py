#!/usr/bin/env python
"""WA PRODUCTION matched filter (GPU) -- W5 truth runs, W6 production pass, rich pass.

Detection path = the FROZEN wa_mf.py measurement core, line for line (20 Hz 1-8 Hz cache, per-(family, node) templates,
rFFT normalised cc averaged over the channels live in both, robust energy floor + cache-mask guard dilated by the template
length, drop-dead rule (> 5 % of the day guarded -> station absent), site merge B with w = n_channels x SNR weight, origin
axis, >= 2 sites, per family-day network median / MAD, declustered peaks >= 8.0 MAD).  On the same stations, sites and
templates the 'data' rows reproduce wa_mf.py (smoke-tested before launch).  What is new never changes a data detection:

  networks   several network sums from ONE set of station cc series (core / E1 / E2 / E3 / all / rich, wa_prod_common);
             each network's sites = 3-km single linkage over ITS member stations present that day
  series     data; c1 (and c2 on truth days) = scrambled control: each SITE's series rolled by its own shift of the day
             (wa_prod_common.site_shift, independent of which other sites are present); rev = TIME-REVERSED templates
             (every channel reversed, same offsets / weights / channels / guard / MAD / declustering)
  templates  W4a library (setup/tpl_library.csv: frozen T1 'build' first, D1 rebuild 'v10b', the 2012-14 re-templates)
             + W2-approved carry-over; the rich pass adds the W4b templates (OUT/templates/tpl_w4b.npz)
  features   per candidate, gathered in the same pass: station cc at the candidate (direct dot product at the aligned
             position; scrambled control: the site's rolled position, exactly as the control was built), each station's
             background at 1000 random valid origins of the day (from the full cc series), site cc / site z (within-site
             w-weighted mean), site_cc_max/2nd, site_z_max/2nd/min, n_z2/3/5, frac_z3, dom_top, site_noise_mad_med,
             z_<STA> for the E3 headline stations, n_sta / n_sites / nch_sum / wsum / n_sta_day, net_re (check), the
             family-day HF station and position (raw HF ratio: wa_prod_hf.py), template in-sample flag
Status per station-day -> OUT/status/mf_<run>_s<shard>.csv (ok / dead / excluded / notemplate / nodata / missing / error);
guarded > 3 % by `wa_prod_common.py guard`.  cuFFT plan cache limited and cleared every day; any GPU error is retried once
and otherwise written as 'error' rows (never a silent loss: the DUC_FIX lesson).

  python scripts/wa_prod_mf.py --run truth --days-file OUT/setup/days_truth.json --nets E3,all --series data,c1,c2,rev
  python scripts/wa_prod_mf.py --run prod --days-file OUT/setup/days_prod.json --nets core,all --series data,c1,rev \
         --shard 0 --nshard 3 --wait-cache
"""
from __future__ import annotations
import os, sys, json, glob, time, argparse, warnings, traceback
os.environ.setdefault("PYTORCH_CUDA_ALLOC_CONF", "expandable_segments:True")
import numpy as np
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import wa_prod_common as P
W = P.W

n = W.NT
NV = W.NDAY - n + 1
NF = int(2 ** np.ceil(np.log2(W.NDAY + n)))


# ------------------------------------------------------------------------------------------------ templates
def load_library(rich=False, tags=None, lib_csv=None):
    """{(fid, netsta, tkey): dict(tpl, off, w, live, tag, clear)} ; tkey = node (own templates)."""
    T = {}
    L = pd.read_csv(lib_csv or f"{P.OUT}/setup/tpl_library.csv", dtype={"live": str})
    L["live"] = L.live.astype(str).str.zfill(3)
    if tags:
        L = L[L.tag.isin(tags)]
    for tag, g in L.groupby("tag", sort=False):
        d = f"{P.OUT}/templates" if os.path.exists(f"{P.OUT}/templates/tpl_{tag}.npz") else f"{P.RES}/templates"
        z = np.load(f"{d}/tpl_{tag}.npz", allow_pickle=False)
        kk = {k: i for i, k in enumerate(z["keys"])}
        tp = z["tpl"]
        for r in g.itertuples():
            key = (int(r.fid), r.netsta, r.node)
            if key in T:
                continue
            T[key] = dict(tpl=tp[kk[f"{r.fid}|{r.netsta}|{r.node}"]].astype(np.float32), off=float(r.offset_s), w=float(r.w),
                          live=np.array([c == "1" for c in r.live]), tag=tag, clear=bool(r.clear) if r.clear == r.clear else False)
    # order = the library csv order (first tag wins), as wa_mf.load_tpl
    order = {(int(r.fid), r.netsta, r.node): i for i, r in enumerate(L.itertuples())}
    if rich and os.path.exists(f"{P.OUT}/templates/tpl_qc_w4b.csv"):
        Q = pd.read_csv(f"{P.OUT}/templates/tpl_qc_w4b.csv", dtype={"live": str})
        Q = Q[Q.use.astype(bool)]
        Q["live"] = Q.live.astype(str).str.zfill(3)
        z = np.load(f"{P.OUT}/templates/tpl_w4b.npz", allow_pickle=False)
        kk = {k: i for i, k in enumerate(z["keys"])}
        for r in Q.itertuples():
            key = (int(r.fid), r.netsta, r.tkey)
            if key in T:
                continue
            T[key] = dict(tpl=z["tpl"][kk[f"{r.fid}|{r.netsta}|{r.tkey}"]].astype(np.float32), off=float(r.offset_s),
                          w=float(r.w), live=np.array([c == "1" for c in r.live]), tag="w4b", clear=bool(r.clear))
            order[key] = len(order)
    keys = sorted(T, key=lambda k: order.get(k, 0))
    return {k: T[k] for k in keys}


def node_maps(rich=False):
    """node -> tkey for the pass: own node, W2-approved carry-over (W4a); + W4b nodes in the rich pass."""
    M = {}
    p = f"{P.OUT}/setup/tpl_nodemap.csv"
    if os.path.exists(p):
        C = pd.read_csv(p)
        for r in C.itertuples():
            if r.kind == "carry" and isinstance(r.src_node, str) and r.src_node:
                M[r.node] = r.src_node
    return M


# ------------------------------------------------------------------------------------------------ helpers
def single_linkage(stas, co):
    lab = {s: s for s in stas}
    for i, a in enumerate(stas):
        for b in stas[i + 1:]:
            if W.km(co[a][0], co[a][1], co[b][0], co[b][1]) < W.SITE_KM:
                la, lb = lab[a], lab[b]
                if la != lb:
                    for k in lab:
                        if lab[k] == lb:
                            lab[k] = la
    groups = {}
    for s in stas:
        groups.setdefault(lab[s], []).append(s)
    out = {}
    for mem in groups.values():
        name = "+".join(sorted(x.split(".")[1] for x in mem))
        for s in mem:
            out[s] = name
    return out


def fs_of_node(node):
    try:
        return float(node.split("@")[0].split(".")[-1])
    except Exception:
        return np.nan


def robust(v):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if len(v) < 20:
        return np.nan, np.nan
    m = float(np.median(v)); s = 1.4826 * float(np.median(np.abs(v - m)))
    return m, (s if s > 0 else np.nan)


_ND = {}


def cache_nodata(month):
    """station-days whose raw file the cache stage read but found no usable data in (empty / no band): legitimate absences."""
    if month not in _ND:
        s = set()
        for m in glob.glob(f"{P.OUT}/cache/manifest_p*_{month}.csv") + [f"{P.RES}/cache_manifest.csv"]:
            if os.path.exists(m):
                M = pd.read_csv(m, usecols=lambda c: c in ("netsta", "date", "status"), dtype=str)
                M = M[~M.status.isin(["ok", "cached"]) & ~M.status.str.startswith("unreadable") & ~M.status.str.startswith("crash")]
                s |= set(zip(M.netsta, M.date))
        _ND[month] = s
    return _ND[month]


def cache_ready(month, wait=True, max_wait=36 * 3600):
    """W6 gate: the month's caches are built (phase 1; phase 2 = dl2 roots for months >= 2012-01)."""
    t0 = time.time()
    while True:
        a = os.path.exists(f"{P.OUT}/cache/done_p1_{month}")
        b = month < "2012-01" or os.path.exists(f"{P.OUT}/cache/done_p2_{month}")
        if a and b:
            return True
        if not wait or time.time() - t0 > max_wait:
            return False
        time.sleep(60)


# ------------------------------------------------------------------------------------------------ main
def main():
    import torch
    import torch.nn.functional as Fnn
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True)
    ap.add_argument("--days-file", default=""); ap.add_argument("--days", default="")
    ap.add_argument("--shard", type=int, default=0); ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--nets", default="core,all")
    ap.add_argument("--series", default="data,c1,rev")
    ap.add_argument("--rich", action="store_true", help="add the W4b templates (rich pass)")
    ap.add_argument("--rich-days", default="", help="json {date: [netsta,...]}: only days with a W4b node live")
    ap.add_argument("--wait-cache", action="store_true")
    ap.add_argument("--fam-chunk", type=int, default=24)
    ap.add_argument("--energy-floor", type=float, default=1e-3)
    ap.add_argument("--drop-dead", type=float, default=0.05)
    ap.add_argument("--stations-csv", default="", help="smoke reproduction: network = these stations, sites as given")
    ap.add_argument("--tags", default="", help="smoke reproduction: template tags (e.g. build)")
    ap.add_argument("--no-features", action="store_true")
    ap.add_argument("--no-excl", action="store_true")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    torch.backends.cuda.cufft_plan_cache.max_size = 32
    days = json.load(open(a.days_file)) if a.days_file else [d for d in a.days.split(",") if d]
    days = [d for j, d in enumerate(sorted(days)) if j % a.nshard == a.shard]
    rich_days = json.load(open(a.rich_days)) if a.rich_days else None
    if rich_days is not None:
        days = [d for d in days if d in rich_days]
    OD = f"{P.OUT}/mf/{a.run}"
    P.ensure(OD, f"{P.OUT}/status")
    NETS = a.nets.split(","); SERIES = a.series.split(",")
    CTLS = [s for s in SERIES if s in ("c1", "c2")]
    REV = "rev" in SERIES
    ST = P.stations()
    co = {r.netsta: (float(r.lat), float(r.lon)) for r in ST.itertuples()}
    if a.stations_csv:
        SO = pd.read_csv(a.stations_csv)
        fixed_site = dict(zip(SO.netsta, SO.site.str.replace(".", "")))
        pool = list(SO.netsta)
        for r in SO.itertuples():
            co[r.netsta] = (float(r.lat), float(r.lon))
    else:
        fixed_site = None
        pool = ST.netsta.tolist()
    CS = P.core_sets() if not a.stations_csv else {}
    E3 = CS.get("E3", [])
    TPL = load_library(rich=a.rich, tags=a.tags.split(",") if a.tags else None)
    NMAP = node_maps(a.rich)
    FAMS = P.families() if os.path.exists(f"{P.OUT}/setup/families.csv") else pd.read_csv(f"{P.RES}/merlin_tests/t1_setup/families.csv")
    FIDS = sorted(FAMS.fid.astype(int))
    X_ = pd.read_csv(f"{P.OUT}/setup/exclusions.csv", dtype=str) if (os.path.exists(f"{P.OUT}/setup/exclusions.csv") and not a.no_excl) else pd.DataFrame(columns=["netsta", "date", "reason"])
    EXCL = {(r.netsta, r.date): r.reason for r in X_.itertuples()}
    RAWI = None
    for tag in ("dl2", "nodl2"):
        pth = f"{P.OUT}/setup/raw_inventory_{tag}.parquet"
        if os.path.exists(pth):
            R_ = pd.read_parquet(pth)
            RAWI = set(zip(R_.netsta, R_.date)); break
    # in-sample template days (frozen re-templates built on 2012-14 production days; W4b build days)
    INS = {}
    for tg, win in P.INSAMPLE_TAGS.items():
        INS[tg] = {d for d, _ in W.window_days(win)}
    W4B_DAYS = {}
    if a.rich and os.path.exists(f"{P.OUT}/templates/w4b_build_days.csv"):
        B = pd.read_csv(f"{P.OUT}/templates/w4b_build_days.csv", dtype=str)
        for r in B.itertuples():
            W4B_DAYS.setdefault((r.netsta, r.tkey), set()).add(r.date)
    stas_with_tpl = {k[1] for k in TPL}
    tkeys_by_sta = {}
    for (f, s, tk) in TPL:
        tkeys_by_sta.setdefault(s, set()).add(tk)
    P.log(f"[mf {a.run} {a.shard}/{a.nshard}] {len(days)} days; nets {NETS}; series {SERIES}; templates {len(TPL)} on "
          f"{len(stas_with_tpl)} stations; rich {a.rich}; device {dev}")

    def guard(X, live, m):
        xd = X.double()
        EN = torch.cumsum(torch.cat([torch.zeros((3, 1), device=dev, dtype=torch.float64), xd * xd], 1), 1)
        ew = (EN[:, n:] - EN[:, :-n])[:, :NV]
        del EN, xd
        lv = torch.from_numpy(live).to(dev)
        ewl = ew[lv]
        med_ch = ewl.median(dim=1).values
        e3 = ewl.sum(0)
        med3 = e3.median()
        ref_ch = torch.maximum(med_ch, 1e-2 * med3)
        bad = ((ewl < a.energy_floor * ref_ch[:, None]).any(0) | (e3 < a.energy_floor * med3) | (ewl.min(0).values <= 0))
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

    from concurrent.futures import ThreadPoolExecutor
    EX = ThreadPoolExecutor(4)

    def load_day(ds):
        out = {}
        for ns in pool:
            if ns not in stas_with_tpl:
                continue
            p = P.find_cache(ns, ds)
            if p is None:
                out[ns] = ("nocache", None); continue
            try:
                z = np.load(p)
                x = z["x"]; m = z["m"]; live = z["live"]; node = str(z["node"])
                if x.shape != (3, W.NDAY):
                    out[ns] = ("error", f"cache shape {x.shape}"); continue
                out[ns] = ("ok", dict(x=x, m=m, live=live, node=node))
            except Exception as e:
                out[ns] = ("error", f"cache unreadable {str(e)[:60]}")
        return out

    sp_ = f"{P.OUT}/status/mf_{a.run}_s{a.shard}.csv"
    status = [tuple(x) for x in pd.read_csv(sp_, dtype=str).fillna("").itertuples(index=False)] if os.path.exists(sp_) else []
    todo = [d for d in days if not os.path.exists(f"{OD}/run_{d}.json")]
    fut = {}
    t_all = time.time()
    for nday, ds in enumerate(todo):
        if a.wait_cache and not cache_ready(ds[:7]):
            status.append((ds, "*", "error", "cache month never completed")); continue
        for dd in todo[nday:nday + 2]:
            if dd not in fut:
                fut[dd] = EX.submit(load_day, dd)
        L = fut.pop(ds).result()
        for attempt in (0, 1):
            try:
                st_rows, nrow = run_day(a, ds, L, TPL, NMAP, NETS, SERIES, CTLS, REV, co, fixed_site, CS, E3, FIDS, EXCL,
                                        RAWI, INS, W4B_DAYS, tkeys_by_sta, guard, torch, Fnn, dev, OD, rich_days)
                status += st_rows
                break
            except Exception as e:
                msg = f"{type(e).__name__}: {str(e)[:120]}"
                P.log(f"  [{a.run} {a.shard}] ERROR {ds} attempt {attempt}: {msg}\n{traceback.format_exc()[-1200:]}")
                if dev == "cuda":
                    torch.cuda.empty_cache(); torch.backends.cuda.cufft_plan_cache.clear()
                if attempt == 1:
                    status += [(ds, ns, "error", f"day job: {msg}") for ns in L if L[ns][0] == "ok"]
        if dev == "cuda":
            torch.cuda.empty_cache(); torch.backends.cuda.cufft_plan_cache.clear()
        if nday % 10 == 0:
            P.write_status(f"mf_{a.run}", f"s{a.shard}", status)
            P.log(f"  [{a.run} {a.shard}] {ds} done {nday + 1}/{len(todo)} {time.time() - t_all:.0f}s")
    P.write_status(f"mf_{a.run}", f"s{a.shard}", status)
    P.log(f"[mf {a.run} {a.shard}] done {len(todo)} days {time.time() - t_all:.0f}s")


def run_day(a, ds, L, TPL, NMAP, NETS, SERIES, CTLS, REV, co, fixed_site, CS, E3, FIDS, EXCL, RAWI, INS, W4B_DAYS,
            tkeys_by_sta, guard, torch, Fnn, dev, OD, rich_days):
    t_day = time.time()
    D0 = pd.Timestamp(ds).timestamp()
    era = P.era_of(ds)
    st_rows = []
    D = {}
    early = []
    for ns, (stt, v) in L.items():
        if stt == "nocache":
            if (ns, ds) in cache_nodata(ds[:7]):
                st_rows.append((ds, ns, "nodata", "raw without usable data (cache stage)")); continue
            if RAWI is not None and (ns, ds) in RAWI:
                st_rows.append((ds, ns, "missing", "raw on disk, no cache")); continue
            st_rows.append((ds, ns, "nodata", "")); continue
        if stt == "error":
            st_rows.append((ds, ns, "error", v)); continue
        if (ns, ds) in EXCL:
            st_rows.append((ds, ns, "excluded", EXCL[(ns, ds)])); continue
        if ns.split(".")[0] in ("XU", "XG") and ds < "2012-01-01":
            st_rows.append((ds, ns, "excluded", "truth array station before 2012")); continue
        node = v["node"]
        tk = node if node in tkeys_by_sta.get(ns, ()) else NMAP.get(node)
        if tk is None or tk not in tkeys_by_sta.get(ns, ()):
            st_rows.append((ds, ns, "notemplate", node)); continue
        if not v["live"].any():
            st_rows.append((ds, ns, "dead", "no live channel")); continue
        D[ns] = dict(x=v["x"], m=v["m"], live=v["live"], node=node, tkey=tk)
    # drop-dead rule (frozen T1 rule: > 5 % of the day guarded -> absent)
    GB = {}
    for ns in list(D):
        X = torch.from_numpy(D[ns]["x"]).to(dev)
        bad, E = guard(X, D[ns]["live"], D[ns]["m"])
        fr = float(bad.float().mean())
        if a.drop_dead > 0 and fr > a.drop_dead:
            st_rows.append((ds, ns, "dead", f"guarded {fr:.3f}")); del D[ns], X, bad, E; continue
        GB[ns] = (X, bad, E)
    for ns in D:
        st_rows.append((ds, ns, "ok", D[ns]["node"] + ("" if D[ns]["tkey"] == D[ns]["node"] else f" <- {D[ns]['tkey']}")))
    # ---------------- network membership
    members = {}
    for net in NETS:
        if net == "core":
            mem = [s for s in D if s in set(CS.get(era, []))] if not fixed_site else list(D)
        elif net in ("E1", "E2", "E3"):
            mem = [s for s in D if s in set(CS.get(net, []))]
        else:                                   # all / rich
            mem = list(D)
        members[net] = mem
    SITE = {}
    for net, mem in members.items():
        if fixed_site:
            SITE[net] = {s: fixed_site[s] for s in mem}
        else:
            SITE[net] = single_linkage(sorted(mem), co)
    # families x stations
    pres = {}                                      # wa_mf order: iterate the template library in its (tag, row) order
    for (f, ns, tk), t in TPL.items():
        if ns in D and D[ns]["tkey"] == tk:
            ch = t["live"] & D[ns]["live"]
            if ch.any():
                pres.setdefault(f, []).append((ns, ch))
    obs, fam_rows = {}, []
    for net in NETS:
        ms = set(members[net])
        for f in FIDS:
            ss = [(s, ch) for s, ch in pres.get(f, []) if s in ms]
            sites = sorted({SITE[net][s] for s, _ in ss})
            fam_rows.append(dict(net=net, fid=f, date=ds, n_sta=len(ss), n_sites=len(sites), sites=",".join(sites),
                                 stations=",".join(s for s, _ in ss), observable=len(sites) >= 2))
            if len(sites) >= 2:
                obs[(net, f)] = ss
    keys = sorted(obs)
    ki = {k: i for i, k in enumerate(keys)}
    meta = dict(date=ds, run=a.run, era=era, nets=NETS, series=SERIES, n_loaded=len(D),
                nodes={s: D[s]["node"] for s in D}, tkeys={s: D[s]["tkey"] for s in D},
                members={k: v for k, v in members.items()}, n_obs={net: sum(1 for k in keys if k[0] == net) for net in NETS})
    if not keys:
        pd.DataFrame(fam_rows).to_parquet(f"{OD}/famday_{ds}.parquet", index=False)
        pd.DataFrame(columns=["fid", "t"]).to_parquet(f"{OD}/det_{ds}.parquet", index=False)
        json.dump(meta, open(f"{OD}/run_{ds}.json", "w"), default=str)
        return st_rows, 0
    # coefficients (site merge B) per network
    coef, wmap, wsite = {}, {}, {}
    for (net, f), ss in obs.items():
        wm = {s: float(ch.sum()) * TPL[(f, s, D[s]["tkey"])]["w"] for s, ch in ss}
        groups = {}
        for s, _ in ss:
            groups.setdefault(SITE[net][s], []).append(s)
        Ws = {q: float(np.sqrt(sum(wm[s] ** 2 for s in mem))) for q, mem in groups.items()}
        Wt = sum(Ws.values())
        for q, mem in groups.items():
            wsum = sum(wm[s] for s in mem)
            for s in mem:
                coef[(net, f, s)] = Ws[q] * wm[s] / wsum / Wt
        wmap[(net, f)] = wm; wsite[(net, f)] = (groups, Ws, Wt)
    shifts = {}
    for net in NETS:
        for k_, c in enumerate(CTLS, 1):
            for q in set(SITE[net].values()):
                shifts[(net, c, q)] = P.site_shift(ds, q, NV, k_)
    nk = len(keys)
    ACC = {ser: torch.zeros((nk, NV), device=dev) for ser in ["data"] + CTLS + (["rev"] if REV else [])}
    VM = torch.ones((nk, NV), dtype=torch.bool, device=dev)
    CM = {c: torch.ones((nk, NV), dtype=torch.bool, device=dev) for c in CTLS}
    rng = np.random.default_rng(1000 + int(ds.replace("-", "")))
    RI = rng.integers(0, NV, P.NRAND)
    RIt = torch.from_numpy(RI).to(dev)
    BG = {}                                         # (f, s, kind) -> station cc at the random origins (np, NaN = invalid)
    by_sta = {}
    for (net, f), ss in obs.items():
        for s, ch in ss:
            by_sta.setdefault(s, {})[f] = ch
    for ns, fch in by_sta.items():
        X, bad, E = GB[ns]
        DH = torch.fft.rfft(X, n=NF)
        valid = ~bad
        pats = {}
        for f, ch in fch.items():
            pats.setdefault(tuple(ch.tolist()), []).append(f)
        tk = D[ns]["tkey"]
        for pat, fl in pats.items():
            chi = [c for c in range(3) if pat[c]]
            for i0 in range(0, len(fl), a.fam_chunk):
                ff = fl[i0:i0 + a.fam_chunk]
                Tn = np.stack([TPL[(f, ns, tk)]["tpl"][chi] for f in ff])
                Wt_ = torch.from_numpy(Tn).to(dev)
                WH = torch.fft.rfft(Wt_, n=NF)
                num = torch.fft.irfft(WH.conj() * DH[chi][None], n=NF)[:, :, :NV]
                cc = (num / E[chi][None]).mean(1)
                del WH, num, Wt_
                cc[:, bad] = 0.0
                if REV:
                    Wr = torch.from_numpy(np.ascontiguousarray(Tn[:, :, ::-1])).to(dev)
                    WHr = torch.fft.rfft(Wr, n=NF)
                    numr = torch.fft.irfft(WHr.conj() * DH[chi][None], n=NF)[:, :, :NV]
                    ccr = (numr / E[chi][None]).mean(1)
                    del WHr, numr, Wr
                    ccr[:, bad] = 0.0
                for j, f in enumerate(ff):
                    d = int(round(TPL[(f, ns, tk)]["off"] * W.FS))
                    a0, b0 = max(0, -d), min(NV, NV - d)
                    tmp = torch.zeros(NV, device=dev)
                    tmp[a0:b0] = cc[j, a0 + d:b0 + d]
                    tv = torch.zeros(NV, dtype=torch.bool, device=dev)
                    tv[a0:b0] = valid[a0 + d:b0 + d]
                    jr = RIt + d
                    okr = (jr >= 0) & (jr < NV)
                    jrc = jr.clamp(0, NV - 1)
                    v = torch.where(okr & valid[jrc], cc[j, jrc], torch.full_like(cc[j, jrc], float("nan")))
                    BG[(f, ns, "n")] = v.cpu().numpy()
                    if REV:
                        tmpr = torch.zeros(NV, device=dev)
                        tmpr[a0:b0] = ccr[j, a0 + d:b0 + d]
                        vr = torch.where(okr & valid[jrc], ccr[j, jrc], torch.full_like(ccr[j, jrc], float("nan")))
                        BG[(f, ns, "r")] = vr.cpu().numpy()
                    for net in NETS:
                        if (net, f) not in obs or (net, f, ns) not in coef:
                            continue
                        k = ki[(net, f)]
                        c = coef[(net, f, ns)]
                        ACC["data"][k] += c * tmp
                        VM[k] &= tv
                        if REV:
                            ACC["rev"][k] += c * tmpr
                        for cser in CTLS:
                            sh = shifts[(net, cser, SITE[net][ns])]
                            ACC[cser][k] += c * torch.roll(tmp, sh)
                            CM[cser][k] &= torch.roll(tv, sh)
                    del tmp, tv
                    if REV:
                        del tmpr
                del cc
                if REV:
                    del ccr
        del DH, valid
    # ---------------- detections (frozen rule)
    rows = []
    famstat = {}
    for ser in ACC:
        MSK = VM if ser in ("data", "rev") else CM[ser]
        SER = ACC[ser]
        for (net, f) in keys:
            k = ki[(net, f)]
            x = SER[k]; vm = MSK[k]
            nvld = int(vm.sum())
            if nvld < 1000:
                famstat[(net, f, ser)] = (np.nan, np.nan, nvld, 0); continue
            am = x[vm]
            med = am.median(); mad = (am - med).abs().median().clamp_min(1e-6)
            xx = torch.where(vm, x, torch.full_like(x, -1e9))
            mx = Fnn.max_pool1d(xx[None, None], kernel_size=2 * (n // 2) + 1, stride=1, padding=n // 2)[0, 0]
            idx = torch.nonzero((xx >= mx - 1e-12) & vm & (xx >= med + W.T_REC * mad)).squeeze(1)
            famstat[(net, f, ser)] = (float(med), float(mad), nvld, int(idx.numel()))
            if idx.numel():
                v = x[idx]
                nm = ((v - med) / mad).cpu().numpy()
                ii = idx.cpu().numpy()
                vv = v.cpu().numpy()
                for j in range(len(ii)):
                    rows.append((net, ser, f, int(ii[j]), D0 + ii[j] / W.FS, float(vv[j]), float(nm[j])))
    del ACC, VM, CM
    if dev == "cuda":
        torch.cuda.empty_cache()
    DET = pd.DataFrame(rows, columns=["net", "series", "fid", "i", "t", "net_cc", "net_mad"])
    DET["date"] = ds
    # ---------------- features (pass 2: direct dot products at the candidate positions)
    if len(DET) and not a.no_features:
        DET = features(DET, D, GB, TPL, obs, coef, wmap, wsite, SITE, shifts, BG, CTLS, E3, ds, torch, dev, famstat)
    # in-sample flag: a contributing station's template was built on this day
    if len(DET):
        ins = {}
        for (net, f), ss in obs.items():
            flag = False
            for s, _ in ss:
                t = TPL[(f, s, D[s]["tkey"])]
                if t["tag"] in INS and ds in INS[t["tag"]]:
                    flag = True; break
                if t["tag"] == "w4b" and ds in W4B_DAYS.get((s, D[s]["tkey"]), ()):
                    flag = True; break
            ins[(net, f)] = flag
        DET["tpl_insample"] = [ins.get((nt, f), False) for nt, f in zip(DET.net, DET.fid)]
    DET["run"] = a.run
    DET["k"] = np.arange(len(DET))
    FD = pd.DataFrame(fam_rows)
    for ser in ["data"] + CTLS + (["rev"] if REV else []):
        FD[f"med_{ser}"] = [famstat.get((nt, f, ser), (np.nan,) * 4)[0] for nt, f in zip(FD.net, FD.fid)]
        FD[f"mad_{ser}"] = [famstat.get((nt, f, ser), (np.nan,) * 4)[1] for nt, f in zip(FD.net, FD.fid)]
        FD[f"nvalid_{ser}"] = [famstat.get((nt, f, ser), (np.nan,) * 4)[2] for nt, f in zip(FD.net, FD.fid)]
        FD[f"nrec_{ser}"] = [famstat.get((nt, f, ser), (np.nan,) * 4)[3] for nt, f in zip(FD.net, FD.fid)]
    FD.to_parquet(f"{OD}/famday_{ds}.parquet", index=False)
    DET.to_parquet(f"{OD}/det_{ds}.parquet.tmp", index=False)
    os.replace(f"{OD}/det_{ds}.parquet.tmp", f"{OD}/det_{ds}.parquet")
    meta.update(n_det={f"{nt}|{s}": int(((DET.net == nt) & (DET.series == s)).sum()) for nt in NETS for s in SERIES},
                shifts={f"{k[0]}|{k[1]}|{k[2]}": v for k, v in shifts.items()}, secs=round(time.time() - t_day, 1),
                sites={net: SITE[net] for net in NETS})
    if "net_re" in DET and len(DET):
        dd = DET[DET.series == "data"]
        meta["max_abs_net_re_minus_net_cc"] = float(np.nanmax(np.abs(dd.net_re - dd.net_cc))) if len(dd) else np.nan
    json.dump(meta, open(f"{OD}/run_{ds}.json.tmp", "w"), default=str)
    os.replace(f"{OD}/run_{ds}.json.tmp", f"{OD}/run_{ds}.json")
    for ns in list(GB):
        del GB[ns]
    P.log(f"  [{a.run}] {ds}: {len(D)} sta; obs {meta['n_obs']}; det " +
          " ".join(f"{k}:{v}" for k, v in meta["n_det"].items() if v) + f"; {meta['secs']}s")
    return st_rows, len(DET)


def features(DET, D, GB, TPL, obs, coef, wmap, wsite, SITE, shifts, BG, CTLS, E3, ds, torch, dev, famstat):
    """station cc at each candidate's aligned position -> site / network features (CPU aggregation)."""
    ar = torch.arange(n, device=dev)
    # requests per station: (row index in DET, f, kind, position j)
    DET = DET.reset_index(drop=True)
    req = {}
    for (net, ser, f), g in DET.groupby(["net", "series", "fid"]).groups.items():
        ss = obs.get((net, f))
        if ss is None:
            continue
        ii = DET.loc[g, "i"].to_numpy(np.int64)
        kind = "r" if ser == "rev" else "n"
        for s, ch in ss:
            d = int(round(TPL[(f, s, D[s]["tkey"])]["off"] * W.FS))
            if ser in ("c1", "c2"):
                sh = shifts[(net, ser, SITE[net][s])]
                j = ((ii - sh) % NV) + d
            else:
                j = ii + d
            req.setdefault(s, []).append((np.asarray(g), f, kind, j, ch))
    CCV = {}                                         # (row, s) -> cc ; stored as dict of arrays per (s)
    for s, lst in req.items():
        X, bad, E = GB[s]
        tk = D[s]["tkey"]
        rows_all = np.concatenate([r for r, _, _, _, _ in lst])
        J = np.concatenate([j for _, _, _, j, _ in lst])
        tix = []
        Tl, CHl = [], []
        tid = {}
        for r, f, kind, j, ch in lst:
            key = (f, kind)
            if key not in tid:
                tp = TPL[(f, s, tk)]["tpl"]
                tid[key] = len(Tl)
                Tl.append(tp[:, ::-1].copy() if kind == "r" else tp)
                CHl.append(ch.astype(np.float32))
            tix.append(np.full(len(r), tid[key], np.int64))
        TI = np.concatenate(tix)
        TT = torch.from_numpy(np.stack(Tl)).to(dev)
        CH = torch.from_numpy(np.stack(CHl)).to(dev)
        out = np.full(len(J), np.nan, np.float32)
        okj = (J >= 0) & (J < NV)
        for b0 in range(0, len(J), 40000):
            sl = slice(b0, b0 + 40000)
            jj = torch.from_numpy(np.clip(J[sl], 0, NV - 1)).to(dev)
            ti = torch.from_numpy(TI[sl]).to(dev)
            win = X[:, (jj[:, None] + ar[None]).reshape(-1)].reshape(3, len(jj), n).permute(1, 0, 2)
            num = (win * TT[ti]).sum(-1)                                   # (R, 3)
            den = E[:, jj].T                                               # (R, 3)
            chm = CH[ti]
            v = ((num / den) * chm).sum(1) / chm.sum(1)
            v = torch.where(bad[jj], torch.full_like(v, float("nan")), v)
            out[sl] = v.cpu().numpy()
            del win, num, den, chm, v, jj, ti
        out[~okj] = np.nan
        CCV[s] = (rows_all, out)
    # assemble station cc per row (dense per station)
    nr = len(DET)
    STAV = {}
    for s, (rr, vv) in CCV.items():
        a_ = np.full(nr, np.nan, np.float64); a_[rr] = vv; STAV[s] = a_
    cols = {c: np.full(nr, np.nan) for c in ["net_re", "site_cc_max", "site_cc_2nd", "site_z_max", "site_z_2nd", "site_z_min",
                                             "n_z2", "n_z3", "n_z5", "frac_z3", "dom_top", "site_noise_mad_med", "wsum",
                                             "nch_sum", "n_sta", "n_sites", "n_sta_day", "hf_pos", "hf_fs"]}
    zc = {s: np.full(nr, np.nan, np.float32) for s in E3}
    hf_sta = np.full(nr, "", dtype=object)
    for (net, ser, f), g in DET.groupby(["net", "series", "fid"]).groups.items():
        ss = obs.get((net, f))
        if ss is None:
            continue
        g = np.asarray(g)
        groups, Ws, Wt = wsite[(net, f)]
        wm = wmap[(net, f)]
        sites = sorted(groups)
        kind = "r" if ser == "rev" else "n"
        SC = np.full((len(g), len(sites)), np.nan); SR = np.full((P.NRAND, len(sites)), np.nan)
        for q_, q in enumerate(sites):
            mem = groups[q]; ws = np.array([wm[s] for s in mem])
            C_ = np.stack([STAV[s][g] if s in STAV else np.full(len(g), np.nan) for s in mem], 1)
            Rr = np.stack([BG[(f, s, kind)] for s in mem], 1)
            SC[:, q_] = (C_ * ws).sum(1) / ws.sum()
            SR[:, q_] = (Rr * ws).sum(1) / ws.sum()
            for s_i, s in enumerate(mem):
                if s in zc:
                    m_, sd_ = robust(BG[(f, s, kind)])
                    zc[s][g] = (C_[:, s_i] - m_) / sd_
        net_re = (SC * np.array([Ws[q] for q in sites])).sum(1) / Wt
        med = np.nanmedian(SR, 0); mad = 1.4826 * np.nanmedian(np.abs(SR - med), 0)
        mad = np.where(mad > 0, mad, np.nan)
        Z = (SC - med) / mad
        contrib = SC * np.array([Ws[q] for q in sites]) / Wt
        Zs = np.sort(np.nan_to_num(Z, nan=-99), 1)[:, ::-1]
        Cs = np.sort(np.nan_to_num(SC, nan=-99), 1)[:, ::-1]
        with np.errstate(invalid="ignore", divide="ignore"):
            dom = np.nanmax(contrib, 1) / np.where(net_re > 0, net_re, np.nan)
            zmin = np.nanmin(np.where(np.isfinite(Z), Z, np.inf), 1)
        zmin[~np.isfinite(zmin)] = np.nan
        ns_ = len(sites)
        cols["net_re"][g] = net_re
        cols["site_cc_max"][g] = np.where(Cs[:, 0] > -98, Cs[:, 0], np.nan)
        cols["site_cc_2nd"][g] = np.where(Cs[:, 1] > -98, Cs[:, 1], np.nan) if ns_ > 1 else np.nan
        cols["site_z_max"][g] = np.where(Zs[:, 0] > -98, Zs[:, 0], np.nan)
        cols["site_z_2nd"][g] = np.where(Zs[:, 1] > -98, Zs[:, 1], np.nan) if ns_ > 1 else np.nan
        cols["site_z_min"][g] = zmin
        cols["n_z2"][g] = (np.nan_to_num(Z, nan=-9) >= 2).sum(1)
        cols["n_z3"][g] = (np.nan_to_num(Z, nan=-9) >= 3).sum(1)
        cols["n_z5"][g] = (np.nan_to_num(Z, nan=-9) >= 5).sum(1)
        cols["frac_z3"][g] = cols["n_z3"][g] / max(ns_, 1)
        cols["dom_top"][g] = dom
        cols["site_noise_mad_med"][g] = float(np.nanmedian(mad)) if np.isfinite(mad).any() else np.nan
        cols["wsum"][g] = float(sum(wm.values()))
        cols["nch_sum"][g] = float(sum(int(ch.sum()) for _, ch in ss))
        cols["n_sta"][g] = len(ss); cols["n_sites"][g] = ns_
        cols["n_sta_day"][g] = len(set(SITE[net]))
        # HF station of the family-day (wa_hf rule: clear first, then w; native >= 40 Hz); position incl. control roll
        cand = sorted([s for s, _ in ss if fs_of_node(D[s]["node"]) >= 40],
                      key=lambda s: (not TPL[(f, s, D[s]["tkey"])]["clear"], -TPL[(f, s, D[s]["tkey"])]["w"]))
        if cand:
            s = cand[0]
            off = TPL[(f, s, D[s]["tkey"])]["off"]
            ii = DET.loc[g, "i"].to_numpy(np.int64)
            if ser in ("c1", "c2"):
                ii = (ii - shifts[(net, ser, SITE[net][s])]) % NV
            hf_sta[g] = s
            cols["hf_pos"][g] = ii / W.FS + off
            cols["hf_fs"][g] = fs_of_node(D[s]["node"])
    for c, v in cols.items():
        DET[c] = v
    DET["hf_sta"] = hf_sta
    for s in E3:
        DET[f"z_{s.split('.')[1]}"] = zc[s]
    return DET


if __name__ == "__main__":
    main()
