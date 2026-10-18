#!/usr/bin/env python
"""Locate INDIVIDUAL detections of the 26 candidate-distinct Oregon families from cross-station
S times, as Savard & Bostock (2015, BSSA 105:2128) locate their detections.

S&B: S arrival times per station from MCCC (VanDecar & Crosson 1990) alignment of the horizontal
waveforms, with a per-station timing uncertainty sigma; >= 4 stations; HYPOINVERSE epicentre with
the depth FIXED (32 km on Vancouver Island); keep solutions with nominal erh < 3 km.

Per detection (det46 time t = the family template's window start), the CODE path is reused unchanged
from xs_detect.py (a literal port of lfedetect.m / cleanwf.m / mccc2.m / ccsearch.m / gettwin.m):
  cleanwf_batch (per-window detrend/demean/unit max + spectral + envelope admission) ->
  stage1_batch (envelope cc, gettwin lag limits) -> lfedetect_one (channels with any envelope cc >
  0.43, envelope MCCC, +-4 s mask around the envelope-stack maximum) -> finish_detect (ccsearch:
  Cmin decreasing by 0.01, waveform lags +-1 s, max sigma < 0.3 s; t_S = tabs - tdel - tdel0;
  Hyp2000 weight = fix(sigma/0.1)).
Differences from main_scc (logged in NOTES.md):
  * the 18 s window is placed PER STATION, centred on the family template's expected S time there
    (a detection is already known, so no scanning).  Expected S = the 1-D (USGS_OR, straight ray)
    S time from a family centre, with the origin anchored to the template stack's clear S picks
    (hypo26/picks.csv, rS >= 1.5; if a family has none, the seed stations' stack maxima).
    Pass 1 centre = the ALP0/ALS0/BWO0 centroid; pass 2 centre = the family's pass-1 median
    epicentre (USGS_OR, erh < 3 km, >= 5 events), windows re-cut.  gettwin limits then apply to the
    RESIDUAL lags after this pre-alignment (they are generous: >= 0.15 s same station, 1.5 D/3.6).
  * no hour normalisation (every later step is per-window amplitude-normalised).
  * min_sta = 3 in the code path so a 3-station tally can be reported; only >= 4 are located.
Location: Hyp2000 (hyp1.40) S-only with Savard's writeHYPO bogus P (3 s before the earliest S,
weight code 3, and WET 1 .5 .2 0 so code 3 has zero weight, as his comment says); depth fixed with
ZTR <z> T at the McCrory (2012) interface depth under the family centre; models USGS_OR (Vp/Vs 1.71,
headline) and GSC (1.73, comparison); ERR default 0.15 s.  All solutions saved; erh is a column.

  python scripts/xs_event_locate.py --out results/chain_paper/event_locs [--workers 12]
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import xs_common as C                                                        # noqa: E402
import xs_detect as XD                                                       # noqa: E402

R = f"{C.RES}/chain_paper"
SEED = ["XA.ALP0", "XA.ALS0", "XA.BWO0"]
ZONLY = {"UW.BRO", "UW.HSO", "UW.MPO", "UW.RNO", "UW.SSO"}
STAS = [s for s in C.STATIONS.netsta if os.path.isdir(f"{C.CACHE}/{s}") and s not in ZONLY]
USGS = pd.read_csv(f"{C.PROJ}/data/velocity_models/usgs_cascadia_v1.7/profiles/oregon_TB2015.csv")
NW = XD.NW
RMIN_TPL = 1.5
LONG_PRE = 12.0                     # stacks_loc.npz / picks.csv sample 0 = t - 12 s (xs_grow)


def log(*a):
    print(f"[{time.strftime('%H:%M:%S')}]", *a, flush=True)


# ------------------------------------------------------------------ geometry ---------------
_G = None


def iface_depth(lat, lon):
    global _G
    if _G is None:
        _G = pd.read_csv(f"{C.PROJ}/data/slab_models/mccrory2012/Modeled_slab_XYZ_Grid_GCS_WGS84.txt")
    d = np.hypot((_G.Longitude - lon) * np.cos(np.radians(lat)), _G.Latitude - lat)
    return float(-_G.Depth_km.iloc[int(np.argmin(d))])


def mean_slowness_s(z):
    """depth-averaged S slowness of USGS_OR from 0 to z km (straight-ray travel time)."""
    zt = USGS.ztop_km.values
    zb = np.minimum(USGS.zbot_km.values, z)
    th = np.clip(zb - zt, 0, None)
    return float((th / USGS.vs.values).sum() / max(th.sum(), 1e-9))


def ts_pred(lat, lon, z, stas):
    from pyproj import Proj
    p = Proj(proj="utm", zone=10)
    x0, y0 = p(lon, lat)
    u = mean_slowness_s(z)
    out = {}
    for s in stas:
        x, y = C.utm_km(s)
        D = np.hypot(x - x0 / 1e3, y - y0 / 1e3)
        out[s] = float(np.hypot(D, z) * u)
    return out


def km_xy(lat, lon):
    from pyproj import Proj
    x, y = Proj(proj="utm", zone=10)(np.asarray(lon), np.asarray(lat))
    return np.asarray(x) / 1e3, np.asarray(y) / 1e3


# ------------------------------------------------------------------ window offsets ---------
FULL = {r.netsta.split(".")[1]: r.netsta for r in C.STATIONS.itertuples()}


def clear_picks(picks, g):
    pg = picks[(picks.grp == g) & picks.tS.notna()]
    pg = pg.assign(ns=pg.sta.map(FULL))
    return pg[(pg.rS >= RMIN_TPL) & pg.ns.isin(STAS)], pg


_IFG = None


def iface_grid():
    global _IFG
    if _IFG is None:
        g = pd.read_csv(f"{C.PROJ}/data/slab_models/mccrory2012/Modeled_slab_XYZ_Grid_GCS_WGS84.txt")
        _IFG = g[g.Latitude.between(43.4, 45.4) & g.Longitude.between(-124.3, -122.3)].reset_index(drop=True)
    return _IFG


def fit_centre(clear):
    """grid search over interface nodes: min L1 misfit of clear S picks with free origin."""
    if len(clear) < 3:
        return None
    G = iface_grid()
    tau = (clear.tS.values - LONG_PRE)
    mis = []
    for r in G.itertuples():
        T = ts_pred(r.Latitude, r.Longitude, -r.Depth_km, list(clear.ns))
        res = tau - np.array([T[s] for s in clear.ns])
        mis.append(np.abs(res - np.median(res)).sum())
    k = int(np.argmin(mis))
    return float(G.Latitude[k]), float(G.Longitude[k])


def family_offsets(centres, picks):
    """{grp: {sta: window start offset (s) after t}}: window centred on the template's clear S pick
    (rS >= 1.5) where there is one, else on the 1-D prediction from the centre, origin anchored to the
    clear picks (or the seed stations' stack maxima when a family has no clear pick)."""
    offs, info = {}, []
    for g, (la, lo) in centres.items():
        z = iface_depth(la, lo)
        T = ts_pred(la, lo, z, STAS)
        clear, pg = clear_picks(picks, g)
        src = "clear"
        anc_set = clear
        if not len(clear):
            anc_set = pg[pg.ns.isin(SEED)]
            src = "seed_max"
        anc = float(np.median([(r.tS - LONG_PRE) - T[r.ns] for r in anc_set.itertuples()]))
        o = {s: anc + T[s] - XD.TLEN / 2 for s in STAS}
        for r in clear.itertuples():
            o[r.ns] = (r.tS - LONG_PRE) - XD.TLEN / 2
        offs[g] = o
        info.append(dict(grp=g, c_lat=la, c_lon=lo, z_iface=z, anchor_src=src, n_clear=len(clear),
                         clear_stas=",".join(clear.sta), anchor_s=anc))
    return offs, pd.DataFrame(info)


# ------------------------------------------------------------------ S times per detection --
PAD = 10.0
SOS = {}


def site_count(stas, r_km=3.0):
    """number of distinct sites (stations < 3 km apart = one site; cf. gettwin clusterflag)."""
    xy = [C.utm_km(s) for s in stas]
    lab = list(range(len(stas)))
    for i in range(len(stas)):
        for j in range(i):
            if np.hypot(xy[i][0] - xy[j][0], xy[i][1] - xy[j][1]) < r_km:
                lab[i] = lab[j]
    return len(set(lab))


def measure_day(args):
    date, dd, offs, band, shift = args
    from scipy.signal import butter, sosfiltfilt
    sos = butter(4, band, "band", fs=C.FS, output="sos") if band else None
    D = {s: C.load_day(s, date) for s in STAS}
    ns = len(STAS)
    twin = XD.gettwin(STAS)
    npad = int(PAD * C.FS)
    rows, trows = [], []
    for r in dd.itertuples():
        g = int(r.grp)
        tdet = r.t + shift
        N = np.zeros((1, ns, NW))
        E = np.zeros((1, ns, NW))
        valid = np.zeros((1, ns), bool)
        wst = np.zeros(ns)
        for i, s in enumerate(STAS):
            if D[s] is None:
                continue
            t0, x, m = D[s]
            ws = tdet + offs[g][s]
            wst[i] = ws
            a = int(round((ws - t0) * C.FS))
            if a - npad < 0 or a + NW + npad > x.shape[1]:
                continue
            if not m[a:a + NW].all():
                continue
            seg = x[1:3, a - npad:a + NW + npad].astype(np.float64)
            if sos is not None:
                seg = sosfiltfilt(sos, seg, axis=1)
            N[0, i] = seg[0, npad:npad + NW]
            E[0, i] = seg[1, npad:npad + NW]
            valid[0, i] = True
        keep, Nc, Ec = XD.cleanwf_batch(N, E, valid)
        rec = dict(det_id=int(r.det_id), grp=g, t=float(r.t), date=date, n_valid=int(valid.sum()),
                   n_clean=int(keep.sum()), n_env_sta=0, n_sta=0, n_site=0, reason="")
        ks = np.flatnonzero(keep[0])
        if ks.size < 3:
            rec["reason"] = "clean<3"
            rows.append(rec)
            continue
        sub = np.concatenate([ks, ks + ns])
        tw = twin[np.ix_(sub, sub)]
        Rr, pairs, _ = XD.stage1_batch(Nc[:, ks], Ec[:, ks], np.ones((1, ks.size), bool), tw)
        o = XD.lfedetect_one(Nc[0, ks], Ec[0, ks], None, tw, XD.CENV, 3, Rr[0].astype(float), pairs)
        if o is None:
            rec["reason"] = "env<3"
            rows.append(rec)
            continue
        rec["n_env_sta"] = int(np.unique(np.concatenate([o["ian"], o["iae"]])).size)
        why = {}
        fd = XD.finish_detect(o, ks.size, 3, XD.SIGR, False, why)
        if fd is None:
            rec["reason"] = why.get("reason", "cs")
            rows.append(rec)
            continue
        st_used = [STAS[ks[si]] for (si, *_) in fd["sta"]]
        rec.update(n_sta=len(fd["sta"]), n_site=site_count(st_used), n_ch=fd["n_ch"],
                   sigr_max=fd["sigr_max"], rmean=fd["rmean_mean"], n_seed=sum(s in SEED for s in st_used),
                   reason="ok")
        rows.append(rec)
        for (si, comp, ts, wt, sg, rm) in fd["sta"]:
            gi = ks[si]
            trows.append(dict(det_id=int(r.det_id), grp=g, sta=STAS[gi], tS=wst[gi] + ts,
                              tS_rel=wst[gi] + ts - tdet, comp=comp, wt=wt, sigr=sg, rmean=rm))
    return rows, trows


def measure(det, offs, workers, band, shift):
    jobs = [(d, dd.iloc[k:k + 250], offs, band, shift) for d, dd in det.groupby("date")
            for k in range(0, len(dd), 250)]
    R1, R2 = [], []
    with ProcessPoolExecutor(workers) as ex:
        for k, (a, b) in enumerate(ex.map(measure_day, jobs)):
            R1 += a
            R2 += b
            if k % 20 == 0:
                log(f"  measured {k + 1}/{len(jobs)} jobs")
    return pd.DataFrame(R1), pd.DataFrame(R2)


# ------------------------------------------------------------------ Hyp2000 ---------------
def _dm(v):
    d = int(abs(v))
    return d, (abs(v) - d) * 60.0


def write_run(wd, TT, stas, model, vpvs, z):
    import ab_hypo as H
    zt, vp = H.MODELS[model]
    os.makedirs(wd, exist_ok=True)
    with open(f"{wd}/sta.sta", "w") as f:
        for s in stas:
            la, lo = C.COORD[s]
            ld, lm = _dm(la)
            od, om = _dm(lo)
            net, st = s.split(".")
            f.write(f"{st:<5s} {net:<2s}  BHZ  {ld:2d} {lm:7.4f} {od:3d} {om:7.4f}W{0:4d}\n")
    with open(f"{wd}/p.crh", "w") as f:
        f.write(f"{model} 1-D P MODEL\n")
        for zz, v in zip(zt, vp):
            f.write(f"{v:5.2f}{zz:5.2f}\n")
    n = 0
    with open(f"{wd}/pha.arc", "w") as f:
        for did, gg in TT.groupby("det_id"):
            tmin = gg.tS.min()
            base_t = pd.Timestamp(np.floor((tmin - 12.0) / 60.0) * 60.0, unit="s")
            b = base_t.timestamp()
            base = base_t.strftime("%Y%m%d%H%M")
            imin = int(gg.tS.values.argmin())
            lines = []
            for k, r in enumerate(gg.itertuples()):
                net, st = r.sta.split(".")
                sS = f"{int(round((r.tS - b) * 100)):5d}"
                if k == imin:                                  # writeHYPO.m bogus P, weight 3
                    sP, pp, wP = f"{int(round((r.tS - 3.0 - b) * 100)):5d}", "IP", "3"
                else:
                    sP, pp, wP = "    0", "  ", "4"
                lines.append(f"{st:<5s}{net:<2s}  BHZ {pp} {wP}{base}{sP}{'':7s}{sS}ES {int(r.wt)}")
            f.write(f"{base}{int(round((tmin - 10.0 - b) * 100)):4d}{'':150s}\n")
            f.write("\n".join(lines) + "\n")
            f.write(f"{'':62s}{int(did):10d}\n")
            n += 1
    with open(f"{wd}/run.hyp", "w") as f:
        f.write(f"""200 T 1900 0
LET 5 2 3 0 0
H71 1 1 3
STA 'sta.sta'
CRH 1 'p.crh'
POS {vpvs:.2f}
ZTR {z:.2f} T
WET 1. .5 .2 0.
MIN 4
ERF F
TOP F
LST 0 0 0
KPR 0
REP F F
COP 3
CAR 1
PHS 'pha.arc'
PRT 'out.prt'
SUM 'out.sum'
ARC 'out.arc'
LOC
STO
""")
    return n


def _loc_one(args):
    import ab_hypo as H
    import ab_hypo_cvm  # noqa: F401  registers USGS_OR
    wd, TT, model, vpvs, z = args
    shutil.rmtree(wd, ignore_errors=True)
    write_run(wd, TT, STAS, model, vpvs, z)
    h, _ = H.run_hyp(wd)
    return h.rename(columns={"grp": "det_id"})


def locate(TT, fams, info, out, tag, workers):
    good = TT.groupby("det_id").sta.transform("size") >= 4
    TT = TT[good]
    jobs = []
    for model, vpvs in (("USGS_OR", 1.71), ("GSC", 1.73)):
        for g in fams:
            tt = TT[TT.grp == g]
            if not len(tt):
                continue
            z = float(info.set_index("grp").at[g, "z_iface"])
            jobs.append((f"{out}/wd_{tag}/{model}/{g}", tt, model, vpvs, z))
    res = []
    with ProcessPoolExecutor(workers) as ex:
        for j, h in zip(jobs, ex.map(_loc_one, jobs)):
            if len(h):
                res.append(h.assign(model=j[2], vpvs=j[3], z_fix=j[4]))
    Hh = pd.concat(res, ignore_index=True)
    gmap = TT.drop_duplicates("det_id").set_index("det_id").grp
    Hh["grp"] = Hh.det_id.map(gmap)
    return Hh


# ------------------------------------------------------------------ summary ---------------
def robust_sigma(v):
    v = np.asarray(v, float)
    return 1.4826 * np.median(np.abs(v - np.median(v)))


def shared_flags(det, win=5.0):
    """True where a detection is within win s of a detection of ANOTHER family."""
    d = det.sort_values("t")
    t = d.t.values
    g = d.grp.values
    sh = np.zeros(len(d), bool)
    j0 = 0
    for i in range(len(d)):
        while t[j0] < t[i] - win:
            j0 += 1
        j = j0
        while j < len(d) and t[j] <= t[i] + win:
            if g[j] != g[i]:
                sh[i] = True
                break
            j += 1
    return pd.Series(sh, index=d.det_id.values)


def summarise(Hh, det, meas, model="USGS_OR", erh_max=3.0):
    h = Hh[(Hh.model == model) & (Hh.erh < erh_max)].copy()
    x, y = km_xy(h.lat.values, h.lon.values)
    h["x"], h["y"] = x, y
    fam_rows = []
    for g, dg in det.groupby("grp"):
        m = meas[meas.grp == g]
        hg = h[h.grp == g]
        allg = Hh[(Hh.model == model) & (Hh.grp == g)]
        r = dict(grp=g, n_det=len(dg), n_sta3=int((m.n_sta >= 3).sum()), n_sta4=int((m.n_sta >= 4).sum()),
                 n_hyp=len(allg), n_erh3=len(hg))
        if len(hg) >= 3:
            mx, my = np.median(hg.x), np.median(hg.y)
            sx, sy = robust_sigma(hg.x), robust_sigma(hg.y)
            rr = np.hypot(hg.x - mx, hg.y - my)
            la, lo = np.median(hg.lat), np.median(hg.lon)
            hu = hg[~hg.shared]
            r.update(med_lat=la, med_lon=lo, med_x=mx, med_y=my, sig_x=sx, sig_y=sy,
                     r50_km=float(np.median(rr)), se_x=1.2533 * sx / np.sqrt(len(hg)),
                     se_y=1.2533 * sy / np.sqrt(len(hg)), med_erh=float(hg.erh.median()),
                     med_rms=float(hg.rms.median()), frac_shared=float(hg.shared.mean()),
                     n_unshared=len(hu),
                     med_x_unsh=float(np.median(hu.x)) if len(hu) >= 3 else np.nan,
                     med_y_unsh=float(np.median(hu.y)) if len(hu) >= 3 else np.nan)
        # activity timing (all 46-day detections of the family)
        tt = pd.to_datetime(dg.t, unit="s")
        r.update(t_p05=tt.quantile(0.05), t_p50=tt.quantile(0.5),
                 peak_day=tt.dt.floor("D").value_counts().idxmax())
        fam_rows.append(r)
    F = pd.DataFrame(fam_rows).sort_values("n_erh3", ascending=False)
    for c in ("med_lat", "med_lon", "med_x", "med_y", "se_x", "se_y", "r50_km", "med_x_unsh", "med_y_unsh"):
        if c not in F:
            F[c] = np.nan
    L = F.dropna(subset=["med_x"]).reset_index(drop=True)
    pr = []
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            a, b = L.iloc[i], L.iloc[j]
            d = float(np.hypot(a.med_x - b.med_x, a.med_y - b.med_y))
            se = float(np.hypot(np.hypot(a.se_x, a.se_y), np.hypot(b.se_x, b.se_y)))
            sp = float(np.hypot(a.r50_km, b.r50_km))
            du = float(np.hypot(a.med_x_unsh - b.med_x_unsh, a.med_y_unsh - b.med_y_unsh))
            pr.append(dict(g1=int(a.grp), g2=int(b.grp), d_km=d, z_sep=d / se if se > 0 else np.inf,
                           d_over_spread=d / sp if sp > 0 else np.inf, d_unshared_km=du))
    P = pd.DataFrame(pr)
    return h, F, P


def components(nodes, edges):
    par = {n: n for n in nodes}

    def f(a):
        while par[a] != a:
            par[a] = par[par[a]]
            a = par[a]
        return a
    for a, b in edges:
        par[f(a)] = f(b)
    return len({f(n) for n in nodes})


def plot_map(h, F, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    G = pd.read_csv(f"{C.PROJ}/data/slab_models/mccrory2012/Modeled_slab_XYZ_Grid_GCS_WGS84.txt")
    G = G[(G.Latitude.between(43.2, 45.6)) & (G.Longitude.between(-124.6, -122.2))]
    fig, ax = plt.subplots(figsize=(9, 10))
    tc = ax.tricontour(G.Longitude, G.Latitude, -G.Depth_km, levels=np.arange(10, 70, 5),
                       colors="0.6", linewidths=0.7)
    ax.clabel(tc, fmt="%d km", fontsize=7)
    L = F.dropna(subset=["med_lat"]).sort_values("grp")
    cmap = plt.get_cmap("tab20")
    rest = h[~h.grp.isin(L.grp)]
    ax.scatter(rest.lon, rest.lat, s=14, color="0.5", alpha=0.7, lw=0, label="families with < 3 located")
    for k, r in enumerate(L.itertuples()):
        c = cmap(k % 20)
        mk = "o" if k < 20 else "s"
        hg = h[h.grp == r.grp]
        ax.scatter(hg.lon, hg.lat, s=16, color=c, alpha=0.8, marker=mk, lw=0)
    for k, r in enumerate(L.itertuples()):
        c = cmap(k % 20)
        ax.scatter(r.med_lon, r.med_lat, s=90, color=c, edgecolor="k", lw=1.0, zorder=5,
                   marker="o" if k < 20 else "s")
        ax.annotate(f"{r.grp} ({r.n_erh3})", (r.med_lon, r.med_lat), fontsize=7, zorder=6,
                    xytext=(4, 4), textcoords="offset points")
    for s in C.STATIONS.itertuples():
        ax.plot(s.lon, s.lat, "k^", ms=6, mfc="w" if s.netsta not in SEED else "r")
        ax.annotate(s.sta, (s.lon, s.lat), fontsize=6, xytext=(3, -8), textcoords="offset points")
    ax.legend(loc="lower left", fontsize=7)
    ax.set_xlim(-124.3, -122.3)
    ax.set_ylim(43.4, 45.4)
    ax.set_aspect(1 / np.cos(np.radians(44.4)))
    ax.set_xlabel("longitude")
    ax.set_ylabel("latitude")
    ax.set_title("Oregon 26 families: individual detections located from cross-station S times\n"
                 "USGS_OR Vp/Vs 1.71, depth fixed to McCrory 2012 interface, erh < 3 km; "
                 "large = family median (label: grp (n))", fontsize=9)
    fig.tight_layout()
    fig.savefig(f"{out}/map.png", dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------ main ------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{R}/event_locs")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--band", default="1,8", help="per-window band-pass (Hz); '' = cache band 1-10 Hz")
    ap.add_argument("--shift", type=float, default=0.0, help="NULL control: move every detection by this many s")
    ap.add_argument("--max-det", type=int, default=0, help="smoke test: detections per family")
    a = ap.parse_args()
    band = tuple(float(v) for v in a.band.split(",")) if a.band else None
    os.makedirs(a.out, exist_ok=True)
    sys.path.insert(0, f"{C.PROJ}/lfe_array_discovery/scripts")
    fam = pd.read_csv(f"{R}/scan26/families.csv")
    fams = sorted(fam[fam.kept].grp.astype(int))
    det = pd.read_parquet(f"{R}/scan26/det46_families.parquet")
    det = det[det.grp.isin(fams)].sort_values("t").reset_index(drop=True)
    det["det_id"] = np.arange(len(det))
    if a.max_det:
        det = det.sort_values("cc", ascending=False).groupby("grp").head(a.max_det).sort_values("t")
    if a.shift:                     # NULL: drop shifted times within 20 s of any real detection
        tr = np.sort(det.t.values)
        ts = det.t.values + a.shift
        k = np.searchsorted(tr, ts)
        dmin = np.minimum(np.abs(ts - tr[np.clip(k - 1, 0, len(tr) - 1)]), np.abs(ts - tr[np.clip(k, 0, len(tr) - 1)]))
        det = det[dmin > 20.0]
        det = det[(det.t + a.shift) < pd.to_datetime(det.date).astype("datetime64[ns]").astype("int64") / 1e9 + 86400 - 60]
    sh = shared_flags(det)
    det["shared"] = det.det_id.map(sh)
    picks = pd.read_csv(f"{R}/hypo26/picks.csv")
    la0 = np.mean([C.COORD[s][0] for s in SEED])
    lo0 = np.mean([C.COORD[s][1] for s in SEED])
    centres = {}
    for g in fams:
        c = fit_centre(clear_picks(picks, g)[0])
        centres[g] = c if c is not None else (la0, lo0)
    log(f"band {band}, shift {a.shift} s")
    for p in (1, 2):
        log(f"pass {p}: {len(det)} detections, {len(fams)} families, {len(STAS)} stations")
        offs, info = family_offsets(centres, picks)
        info.to_csv(f"{a.out}/centres_pass{p}.csv", index=False)
        meas, TT = measure(det, offs, a.workers, band, a.shift)
        meas.to_parquet(f"{a.out}/meas_pass{p}.parquet", index=False)
        TT.to_parquet(f"{a.out}/stimes_pass{p}.parquet", index=False)
        log(f"  n_sta>=3: {(meas.n_sta >= 3).sum()}, >=4: {(meas.n_sta >= 4).sum()}, n_site>=4: "
            f"{(meas.n_site >= 4).sum()}; reasons {meas.reason.value_counts().to_dict()}")
        if not (meas.n_sta >= 4).any():
            log("  nothing to locate")
            break
        Hh = locate(TT, fams, info, a.out, f"pass{p}", a.workers)
        mi = meas.set_index("det_id")
        Hh["shared"] = Hh.det_id.map(det.set_index("det_id").shared)
        for c in ("n_sta", "n_site", "n_seed", "sigr_max", "rmean"):
            Hh[c] = Hh.det_id.map(mi[c])
        Hh.to_csv(f"{a.out}/hyp_pass{p}.csv", index=False)
        log(f"  located {Hh.groupby('model').size().to_dict()}; erh<3 "
            f"{Hh[Hh.erh < 3].groupby('model').size().to_dict()}")
        h, F, P = summarise(Hh, det, meas)
        F.to_csv(f"{a.out}/family_summary_pass{p}.csv", index=False)
        new = {}
        for r in F.itertuples():
            ok = r.n_erh3 >= 5 and np.isfinite(getattr(r, "med_lat", np.nan))
            new[r.grp] = (float(r.med_lat), float(r.med_lon)) if ok else centres[r.grp]
        centres = new
    else:
        F.to_csv(f"{a.out}/family_summary.csv", index=False)
        P.to_csv(f"{a.out}/family_pairs.csv", index=False)
        hg, Fg, _ = summarise(Hh, det, meas, model="GSC")
        if "med_lat" in Fg:
            Fg[["grp", "n_erh3", "med_lat", "med_lon", "r50_km"]].rename(
                columns=lambda c: c if c == "grp" else f"gsc_{c}").merge(F, on="grp").to_csv(
                f"{a.out}/family_summary_with_gsc.csv", index=False)
        if "med_lat" in F:
            plot_map(h, F, a.out)
            nodes = list(F.dropna(subset=["med_x"]).grp)
            if not len(P):
                P = pd.DataFrame(columns=["g1", "g2", "z_sep", "d_over_spread"])
            for lab, e in (("z_sep<3 (medians not significantly apart)", P[P.z_sep < 3]),
                           ("d<spread (clouds overlap)", P[P.d_over_spread < 1])):
                log(f"  distinct groups by {lab}: {components(nodes, zip(e.g1, e.g2))} of {len(nodes)}")
    open(f"{a.out}/EVLOC_DONE", "w").write("ok\n")


if __name__ == "__main__":
    main()
