#!/usr/bin/env python
# FROZEN COPY for avenue A4 (2026-09-25) of the file named without the a4f_ prefix (md5 at copy time in results/a4_seed/A4_LOG.md); do not edit.
"""STAGE 6 -- Royer & Bostock (2014) sec.5 / sec.6 location: Hilbert-rotated templates, zero-phase P
and S picks, Hyp2000 (Klein 2002), then hypoDD (Waldhauser 2001) with Table S1 EXACTLY.

  sec.6 (i)  "each 3-component LFE template seismogram is normalized to unit, maximum absolute
             amplitude and phase rotated by 90 degrees to aid in identification and picking of
             zero-phase, primary arrivals"
  sec.6 (ii) "locations are determined from travel-time picks corresponding to those P- and
             S-arrivals judged to display unambiguous, zero-phase waveforms"
  sec.5      Hyp2000 locations, then hypoDD "to a combination of ordinary phase picks from our LFE
             catalogue and differential travel times from phase correlation of P- and S-waves",
             Table S1: IDAT 3; 5 sets x NITER 8; WTCCP/WTCCS .01,.01,1,1,1; WRCC -9,-9,-9,5,5;
             WDCC -9,-9,2,2,.5; WTCTP/WTCTS 1,1,.01,.01,.01; WRCT -9,-9,6,6,6; WDCT -9,8,8,8,8;
             DAMP 50,50,40,40,40.

Programs (built from the public sources, tools/README): /home/jovyan/tools/bin/hyp1.40
(Hypoinverse-2000 v1.40, USGS) and /home/jovyan/tools/bin/hypoDD + ph2dt (HYPODD 1.3, Waldhauser).

Not given in the paper, resolved from the authors' group (headline in brackets):
  * velocity model -- Royer's thesis: "the 1-D, fixed VP/VS model used for LFE location" (unnamed);
    Bostock's group used the GSC routine-location 1-D model -> [GSC onshore] (see MODELS), with
    Crosson (1976, PNSN) as a column;
  * Vp/Vs [1.73] (the GSC uniform value) / 1.78 / 1.85;
  * the "unambiguous" test for a zero-phase pick: the rotated trace's extremum must exceed every
    other extremum in the search window (outside +-0.5 s) by the ratio R: 1.2 / [1.5] / 2.0;
  * P/S search windows: S = largest horizontal extremum inside the stack, P = largest vertical
    extremum 1.5-14 s before S (S-P for 25-60 km deep sources within ~100 km);
  * ph2dt parameters: the program's documented defaults; hypoDD ISTART 2 (from the Hyp2000
    locations), ISOLV 2 (LSQR, implied by DAMP), OBSCC 0, OBSCT 8, DIST 200 km;
  * cc differential times: +-1.0 s (P) / +-1.5 s (S) windows about the picks, lag search +-0.5 s,
    weight cc^2, kept for cc >= 0.7 (the hypoDD manual's example practice).
Hyp2000's own distance/residual weighting: program defaults.

Station bootstrap (per family, Hyp2000): stations resampled with replacement --boot times; the
68th-percentile epicentral and depth deviations are reported (as for the calibrated locator).

  python scripts/ab_hypo.py --sub NW --det .../grow/det_it3.parquet --out .../hypo
-> <out>/stacks_loc.npz, picks.csv, station_snr.csv, hyp_<vpvs>_R<r>.csv, dd_<vpvs>_R<r>.csv,
   loc_summary.csv
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy.signal import hilbert

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import a4f_common as A                                                       # noqa: E402
import a4f_grow as G                                                         # noqa: E402

BIN = "/home/jovyan/tools/bin"
# 1-D models.  HEADLINE = the GSC routine-location ONSHORE model, the model Bostock's group used
# (Bostock et al. 2015 sec.3: "predicted times from the 1-D model used for routine earthquake
# location by the Geological Survey of Canada"; Savard 2018 UBC PhD, "the 1D GSC velocity model
# used for routine locations" with HYPOINVERSE; Royer 2014 PhD: "the 1-D, fixed VP/VS model used
# for LFE location").  Layers digitised EXACTLY from the vector graphics of Savard (2018) Fig.4.7a
# ("GSC onshore model"; thesis PDF page 70, circle.library.ubc.ca bitstream 333260):
#   0-1 km 5.01, 1-6 km 6.01, 6-30 km 6.71, 30-46 km 7.11, >46 km 7.76 km/s.
# GSC scales Vp by a uniform Vp/Vs of 1.73.  Crosson (1976, PNSN) is kept as a column.
MODELS = {
    "GSC": ([0.0, 1.0, 6.0, 30.0, 46.0], [5.01, 6.01, 6.71, 7.11, 7.76]),
    "CRO": ([0.0, 4.0, 9.0, 16.0, 20.0, 25.0, 32.0], [5.40, 6.38, 6.59, 6.73, 6.86, 7.00, 7.90]),
}
ZTOP, VP = MODELS["GSC"]
VPVS = [1.73, 1.78, 1.85]
RMIN = [1.2, 1.5, 2.0]
HEAD = ("GSC", 1.73, 1.5)
TABLE_S1 = [  # NITER WTCCP WTCCS WRCC WDCC WTCTP WTCTS WRCT WDCT DAMP
    (8, 0.01, 0.01, -9, -9, 1, 1, -9, -9, 50),
    (8, 0.01, 0.01, -9, -9, 1, 1, -9, 8, 50),
    (8, 1, 1, -9, 2, 0.01, 0.01, 6, 8, 40),
    (8, 1, 1, 5, 2, 0.01, 0.01, 6, 8, 40),
    (8, 1, 1, 5, 0.5, 0.01, 0.01, 6, 8, 40)]
T0 = pd.Timestamp("2000-01-01 00:00:00")
PRE = G.LONG_PRE                  # stack sample 0 = detection time - PRE


def rotate(st):
    """sec.6(i): 90-degree phase rotation (Hilbert) then unit max across 3 components."""
    r = np.imag(hilbert(st, axis=-1))
    m = np.abs(r).max()
    return r / m if m > 0 else r


def pick_one(x, t0, t1, excl=0.5):
    """Largest |extremum| of x in samples [t0, t1); ratio to the largest outside +-excl s."""
    t0, t1 = max(int(t0), 0), min(int(t1), x.size)
    if t1 - t0 < 4:
        return None
    w = np.abs(x[t0:t1])
    k = int(np.argmax(w))
    e = int(round(excl * A.FS))
    rest = np.r_[w[:max(k - e, 0)], w[k + e + 1:]]
    ratio = float(w[k] / rest.max()) if rest.size and rest.max() > 0 else np.inf
    return t0 + k, ratio, float(np.sign(x[t0 + k]))


def make_picks(ST):
    rows = []
    for g, T in ST.items():
        for s, st in T.items():
            if not np.isfinite(st).all() or not (st ** 2).sum() > 0:
                continue
            r = rotate(st)
            n = r.shape[1]
            has_h = (r[1] ** 2).sum() > 0 and (r[2] ** 2).sum() > 0
            env = np.sqrt((np.abs(hilbert(st, axis=-1)) ** 2).sum(0))
            snr = float(env.max() / np.median(env)) if np.median(env) > 0 else np.nan
            row = dict(grp=int(g), sta=s, snr=snr, three_c=bool(has_h))
            if has_h:
                h = np.where(np.abs(r[1]) >= np.abs(r[2]), r[1], r[2])
                sp = pick_one(h, 0.1 * n, 0.9 * n)
                if sp:
                    kS, rS, pS = sp
                    row.update(tS=kS / A.FS, rS=rS, polS=pS)
                    pp = pick_one(r[0], kS - 14 * A.FS, kS - 1.5 * A.FS)
                    if pp:
                        row.update(tP=pp[0] / A.FS, rP=pp[1], polP=pp[2])
            else:                                   # vertical-only: P only, whole-stack search
                pp = pick_one(r[0], 0.1 * n, 0.9 * n)
                if pp:
                    row.update(tP=pp[0] / A.FS, rP=pp[1], polP=pp[2])
            rows.append(row)
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------ Hyp2000 ------------
def _dm(v):
    d = int(abs(v))
    return d, (abs(v) - d) * 60.0


def write_hyp_inputs(wd, P, stas, vpvs, rmin, boot_idx=None, model="GSC"):
    zt, vp = MODELS[model]
    os.makedirs(wd, exist_ok=True)
    with open(f"{wd}/sta.sta", "w") as f:
        for s in stas:
            la, lo, el = A.COORD[s]
            ld, lm = _dm(la)
            od, om = _dm(lo)
            f.write(f"{s:<5s} {A.NET[s]:<2s}  BHZ  {ld:2d} {lm:7.4f} {od:3d} {om:7.4f}W"
                    f"{int(round(el)):4d}\n")
    with open(f"{wd}/p.crh", "w") as f:
        f.write(f"{model} 1-D P MODEL\n")
        for z, v in zip(zt, vp):
            f.write(f"{v:5.2f}{z:5.2f}\n")
    ids = []
    with open(f"{wd}/pha.arc", "w") as f:
        for g, gg in P.groupby("grp"):
            if boot_idx is not None:
                gg = gg.iloc[boot_idx[g]] if g in boot_idx else gg.iloc[:0]
            lines = []
            ot = T0 + pd.Timedelta(minutes=2 * int(g))
            base = ot.strftime("%Y%m%d%H%M")
            for r in gg.itertuples():
                okP = np.isfinite(getattr(r, "tP", np.nan)) and getattr(r, "rP", 0) >= rmin
                okS = np.isfinite(getattr(r, "tS", np.nan)) and getattr(r, "rS", 0) >= rmin
                if not (okP or okS):
                    continue
                sP = f"{int(round((10 + r.tP) * 100)):5d}" if okP else "    0"
                sS = f"{int(round((10 + r.tS) * 100)):5d}" if okS else "    0"
                wP = "0" if okP else "4"
                wS = "0" if okS else "4"
                ln = (f"{r.sta:<5s}{A.NET[r.sta]:<2s}  BHZ "
                      f"{'IP' if okP else '  '} {wP}{base}{sP}"
                      f"{'':7s}{sS}{'ES' if okS else '  '} {wS}")
                lines.append(ln)
            if len(lines) < 2:
                continue
            ids.append(int(g))
            f.write(f"{base}{1000:4d}{'':150s}\n")
            f.write("\n".join(lines) + "\n")
            f.write(f"{'':62s}{int(g):10d}\n")
    with open(f"{wd}/run.hyp", "w") as f:
        f.write(f"""200 T 1900 0
LET 5 2 3 0 0
H71 1 1 3
STA 'sta.sta'
CRH 1 'p.crh'
POS {vpvs:.2f}
ZTR 35
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
    return ids


def run_hyp(wd):
    r = subprocess.run([f"{BIN}/hyp1.40"], input="@run.hyp\n", cwd=wd, text=True,
                       capture_output=True, timeout=1800)
    rows = []
    p = f"{wd}/out.sum"
    if not os.path.exists(p):
        return pd.DataFrame(), r.stdout[-2000:]
    for ln in open(p):
        if len(ln) < 94:
            continue
        try:
            la = float(ln[16:18]) + float(ln[19:23]) / 6000.0
            lo = -(float(ln[23:26]) + float(ln[27:31]) / 6000.0)
            dp = float(ln[31:36]) / 100.0
            num = int(ln[39:42] or 0)
            gap = int(ln[42:45] or 0)
            rms = float(ln[48:52] or 0) / 100.0
            erh = float(ln[85:89] or 0) / 100.0
            erz = float(ln[89:93] or 0) / 100.0
            yr, mo, dy, hh, mm = ln[0:4], ln[4:6], ln[6:8], ln[8:10], ln[10:12]
            sec = float(ln[12:16]) / 100.0
            ot = pd.Timestamp(f"{yr}-{mo}-{dy} {hh}:{mm}:00") + pd.Timedelta(seconds=sec)
            grp = int(ln[136:146])                  # event ID = template id (terminator)
            rows.append(dict(grp=grp, lat=la, lon=lo, dep=dp, num=num, gap=gap, rms=rms,
                             erh=erh, erz=erz, ot_s=(ot - T0).total_seconds()))
        except ValueError:
            continue
    return pd.DataFrame(rows), r.stdout[-2000:]


# ------------------------------------------------------------------------ hypoDD -------------
def cc_dt(P, ROT, pairs, maxlag=0.5):
    """Differential times from phase correlation of the ROTATED P (Z) and S (N+E) waveforms, for
    the event pairs ph2dt linked (dt.ct) at every common station.  Returns rows g1 g2 sta ph cc
    dt_pick lag (lag = shift of event 2 that best aligns it on event 1)."""
    out = []
    Pk = P.set_index(["grp", "sta"])
    L = int(round(maxlag * A.FS))
    have = set(Pk.index)
    for g1, g2 in pairs:
        for s in P.sta[P.grp == g1]:
            if (g2, s) not in have:
                continue
            for ph, half, comp in (("P", 1.0, [0]), ("S", 1.5, [1, 2])):
                t1 = Pk.loc[(g1, s)].get("t" + ph, np.nan)
                t2 = Pk.loc[(g2, s)].get("t" + ph, np.nan)
                if not (np.isfinite(t1) and np.isfinite(t2)):
                    continue
                x1, x2 = ROT[(g1, s)][comp], ROT[(g2, s)][comp]
                k1, k2 = int(round(t1 * A.FS)), int(round(t2 * A.FS))
                h = int(round(half * A.FS))
                if (k1 - h < 0 or k1 + h >= x1.shape[1] or k2 - h - L < 0
                        or k2 + h + L >= x2.shape[1]):
                    continue
                a = x1[:, k1 - h:k1 + h + 1]
                b = np.stack([x2[:, k2 - h + l:k2 + h + 1 + l] for l in range(-L, L + 1)])
                den = np.sqrt((a ** 2).sum() * (b ** 2).sum((1, 2)))
                c = np.where(den > 0, (b * a[None]).sum((1, 2)) / np.where(den > 0, den, 1), 0)
                i = int(np.argmax(c))
                out.append(dict(g1=g1, g2=g2, sta=s, ph=ph, cc=float(c[i]),
                                dt_pick=(k1 - k2) / A.FS, lag=(i - L) / A.FS))
    return pd.DataFrame(out)


def run_dd(wd, P, H, ROT, stas, vpvs, rmin, ccmin=0.7, model="GSC"):
    zt, vp = MODELS[model]
    """ph2dt (catalogue) + dt.cc (phase correlation) + hypoDD with Table S1."""
    os.makedirs(wd, exist_ok=True)
    H = H.set_index("grp")
    with open(f"{wd}/station.dat", "w") as f:
        for s in stas:
            f.write(f"{s:<7s} {A.COORD[s][0]:10.5f} {A.COORD[s][1]:11.5f}\n")
    with open(f"{wd}/phase.dat", "w") as f:
        for g, gg in P.groupby("grp"):
            if g not in H.index:
                continue
            h = H.loc[g]
            ot = T0 + pd.Timedelta(seconds=float(h.ot_s))
            f.write(f"# {ot.year} {ot.month} {ot.day} {ot.hour} {ot.minute} "
                    f"{ot.second + ot.microsecond / 1e6:.2f} {h.lat:.5f} {h.lon:.5f} {h.dep:.3f} "
                    f"0.0 {h.erh:.2f} {h.erz:.2f} {h.rms:.2f} {int(g)}\n")
            tref = (T0 + pd.Timedelta(minutes=2 * int(g)) - ot).total_seconds() + 10.0
            for r in gg.itertuples():
                if np.isfinite(getattr(r, "tP", np.nan)) and r.rP >= rmin:
                    f.write(f"{r.sta:<7s} {tref + r.tP:8.3f} 1.000 P\n")
                if np.isfinite(getattr(r, "tS", np.nan)) and r.rS >= rmin:
                    f.write(f"{r.sta:<7s} {tref + r.tS:8.3f} 1.000 S\n")
    with open(f"{wd}/ph2dt.inp", "w") as f:
        f.write("station.dat\nphase.dat\n"
                "* MINWGHT MAXDIST MAXSEP MAXNGH MINLNK MINOBS MAXOBS (documented defaults)\n"
                "0 200 10 10 8 8 20\n")
    subprocess.run([f"{BIN}/ph2dt", "ph2dt.inp"], cwd=wd, capture_output=True, text=True,
                   timeout=1800)
    # cross-correlation differential times for the pairs ph2dt linked
    pairs = []
    if os.path.exists(f"{wd}/dt.ct"):
        for ln in open(f"{wd}/dt.ct"):
            if ln.startswith("#"):
                x = ln.split()
                pairs.append((int(x[1]), int(x[2])))
    CC = cc_dt(P, ROT, pairs)
    CC.to_csv(f"{wd}/cc_dt.csv", index=False)
    ok = P.set_index(["grp", "sta"])
    with open(f"{wd}/dt.cc", "w") as f:
        if len(CC):
            C = CC[(CC.cc >= ccmin) & CC.g1.isin(H.index) & CC.g2.isin(H.index)]
            for (g1, g2), cc in C.groupby(["g1", "g2"]):
                ot1, ot2 = float(H.loc[g1].ot_s), float(H.loc[g2].ot_s)
                off1 = 120.0 * g1 + 10.0 - ot1
                off2 = 120.0 * g2 + 10.0 - ot2
                f.write(f"# {int(g1)} {int(g2)} 0.0\n")
                for r in cc.itertuples():
                    rr = ok.loc[(g1, r.sta)]
                    rq = ok.loc[(g2, r.sta)]
                    if r.ph == "P" and not (rr.get("rP", 0) >= rmin and rq.get("rP", 0) >= rmin):
                        continue
                    if r.ph == "S" and not (rr.get("rS", 0) >= rmin and rq.get("rS", 0) >= rmin):
                        continue
                    # tt1 - tt2 = (pick1 + off1) - (pick2 + lag-corrected pick2 + off2)
                    dt = (r.dt_pick - r.lag) + off1 - off2
                    f.write(f"{r.sta:<7s} {dt:9.4f} {r.cc ** 2:6.3f} {r.ph}\n")
    rows = "\n".join(" ".join(str(v) for v in t) for t in TABLE_S1)
    with open(f"{wd}/hypoDD.inp", "w") as f:
        f.write(f"""dt.cc
dt.ct
event.sel
station.dat
hypoDD.loc
hypoDD.reloc
hypoDD.sta
hypoDD.res
hypoDD.src
* IDAT IPHA DIST
3 3 200
* OBSCC OBSCT
0 8
* ISTART ISOLV NSET
2 2 {len(TABLE_S1)}
* NITER WTCCP WTCCS WRCC WDCC WTCTP WTCTS WRCT WDCT DAMP  (Royer & Bostock 2014 Table S1)
{rows}
* NLAY RATIO
{len(zt)} {vpvs:.2f}
{' '.join(f'{z:.1f}' for z in zt)}
{' '.join(f'{v:.2f}' for v in vp)}
* CID
0
* ID
""")
    if not os.path.exists(f"{wd}/event.sel") and os.path.exists(f"{wd}/event.dat"):
        shutil.copy(f"{wd}/event.dat", f"{wd}/event.sel")
    r = subprocess.run([f"{BIN}/hypoDD", "hypoDD.inp"], cwd=wd, capture_output=True, text=True,
                       timeout=7200)
    out = []
    p = f"{wd}/hypoDD.reloc"
    if os.path.exists(p) and os.path.getsize(p) > 0:
        D = np.loadtxt(p, ndmin=2)
        for d in D:
            out.append(dict(grp=int(d[0]), dd_lat=d[1], dd_lon=d[2], dd_dep=d[3],
                             dd_ex_m=d[7], dd_ey_m=d[8], dd_ez_m=d[9]))
    return pd.DataFrame(out), (r.stdout + r.stderr)[-3000:]


def bootstrap(wd, P, stas, vpvs, rmin, H, nboot, seed=0, model="GSC"):
    rng = np.random.default_rng(seed)
    res = {int(g): [] for g in H.grp}
    for b in range(nboot):
        idx = {}
        for g, gg in P.groupby("grp"):
            n = len(gg)
            idx[g] = rng.integers(0, n, n) if n else np.zeros(0, int)
        bw = f"{wd}/boot"
        write_hyp_inputs(bw, P, stas, vpvs, rmin, boot_idx=idx, model=model)
        B, _ = run_hyp(bw)
        for r in B.itertuples():
            if r.grp in res:
                res[r.grp].append((r.lat, r.lon, r.dep))
    rows = []
    Hi = H.set_index("grp")
    for g, v in res.items():
        if len(v) >= 10:
            v = np.array(v)
            h = Hi.loc[g]
            rows.append(dict(grp=g, boot_km=float(np.quantile(
                A.epi_km(v[:, 0], v[:, 1], h.lat, h.lon), 0.68)),
                boot_dz=float(np.quantile(np.abs(v[:, 2] - h.dep), 0.68)), boot_n=len(v)))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", required=True, choices=["NW", "SW", "NVI"])
    ap.add_argument("--det", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--pool", default="", help="stations stacked/picked; default A.POOL[sub]; "
                    "'ext:<csv>:<R>' = extension table within R km")
    ap.add_argument("--boot", type=int, default=60)
    ap.add_argument("--vpvs", default=",".join(map(str, VPVS)))
    ap.add_argument("--models", default="GSC,CRO")
    ap.add_argument("--rmin", default=",".join(map(str, RMIN)))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    stas = A.pool_from_arg(a.pool, a.sub)
    det = pd.read_parquet(a.det)
    A.log(f"[hypo {a.sub}] {det.grp.nunique()} templates, {len(det)} detections, "
          f"{len(stas)} stations")
    sp = f"{a.out}/stacks_loc.npz"
    if os.path.exists(sp):
        ST = G.load_stacks(sp)
    else:
        ST = G.restack(A.DayCache(stas), stas, det)
        np.savez_compressed(sp, **{f"{g}_{s}": ST[g][s].astype(np.float32)
                                   for g in ST for s in ST[g]})
    P = make_picks(ST)
    P.to_csv(f"{a.out}/picks.csv", index=False)
    # which stations carry each family: per-station SNR (permanent = non-XU) -- a COLUMN table
    Q = P[["grp", "sta", "snr", "three_c"]].copy()
    Q["net"] = [A.NET[s] for s in Q.sta]
    Q["permanent"] = ~Q.net.isin(["XU", "XG", "ZK", "Z3", "7D", "ZU", "YW"])
    Q.sort_values(["grp", "snr"], ascending=[True, False]).to_csv(
        f"{a.out}/station_snr.csv", index=False)
    A.log(f"[hypo] {len(P)} station stacks; P picks {P.get('tP', pd.Series()).notna().sum()}, "
          f"S picks {P.get('tS', pd.Series()).notna().sum()}")
    ROT = {(g, s): rotate(ST[g][s]) for g in ST for s in ST[g]}
    summ = []
    combos = [(m, float(v), float(r)) for m in a.models.split(",")
              for v in a.vpvs.split(",") for r in a.rmin.split(",")]
    for model, vpvs, rmin in combos:
        if True:
            tag = f"{model}_{vpvs:.2f}_R{rmin:.1f}"
            wd = f"{a.out}/run_{tag}"
            ids = write_hyp_inputs(wd, P, stas, vpvs, rmin, model=model)
            H, log = run_hyp(wd)
            if not len(H):
                A.log(f"[hypo {tag}] Hyp2000 produced nothing; tail:\n{log}")
                summ.append(dict(model=model, vpvs=vpvs, rmin=rmin, n_input=len(ids), n_hyp=0))
                continue
            is_head = (model == HEAD[0] and abs(vpvs - HEAD[1]) < 1e-6
                       and abs(rmin - HEAD[2]) < 1e-6)
            if a.boot and is_head:
                Bt = bootstrap(wd, P, stas, vpvs, rmin, H, a.boot, model=model)
                H = H.merge(Bt, on="grp", how="left")
            D, dlog = run_dd(f"{wd}/dd", P, H, ROT, stas, vpvs, rmin, model=model)
            open(f"{wd}/dd/hypoDD_stdout.txt", "w").write(dlog)
            M = H.merge(D, on="grp", how="left") if len(D) else H.assign(dd_dep=np.nan)
            M["dd_kept"] = M.dd_dep.notna()
            M["model"], M["vpvs"], M["rmin"], M["headline"] = model, vpvs, rmin, is_head
            M.to_csv(f"{a.out}/loc_{tag}.csv", index=False)
            if is_head:
                M.to_csv(f"{a.out}/loc_headline.csv", index=False)
            lo, hi = A.SUBARRAY[a.sub]["depth_pub"]
            summ.append(dict(
                model=model, vpvs=vpvs, rmin=rmin, headline=is_head, n_input=len(ids), n_hyp=len(H),
                hyp_dep_med=float(H.dep.median()), hyp_dep_p05=float(H.dep.quantile(.05)),
                hyp_dep_p95=float(H.dep.quantile(.95)),
                hyp_frac_in_band=float(H.dep.between(lo, hi).mean()),
                hyp_erh_mean=float(H.erh.mean()), hyp_erz_mean=float(H.erz.mean()),
                n_dd=int(M.dd_kept.sum()), dd_keep_frac=float(M.dd_kept.mean()),
                dd_dep_med=float(M.dd_dep.median()) if M.dd_kept.any() else np.nan,
                dd_frac_in_band=float(M.dd_dep[M.dd_kept].between(lo, hi).mean())
                if M.dd_kept.any() else np.nan,
                boot_km_med=float(M.boot_km.median()) if "boot_km" in M else np.nan,
                published_band=f"{lo}-{hi}",
                published_sig_h=A.SUBARRAY[a.sub]["sig_h_pub"],
                published_sig_z=A.SUBARRAY[a.sub]["sig_z_pub"],
                published_dd_keep=A.SUBARRAY[a.sub]["dd_keep_pub"]))
            A.log(f"[hypo {tag}] Hyp2000 {len(H)}/{len(ids)} located, depth median "
                  f"{H.dep.median():.1f} km; hypoDD kept {int(M.dd_kept.sum())}")
    S = pd.DataFrame(summ)
    S.to_csv(f"{a.out}/loc_summary.csv", index=False)
    A.log("\n" + S.to_string(index=False))


if __name__ == "__main__":
    main()
