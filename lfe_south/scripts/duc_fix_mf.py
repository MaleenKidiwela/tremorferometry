#!/usr/bin/env python
"""DUC_FIX matched filter -- Ducellier's detector (find_all_LFEs_parallel.py analyze_data, verbatim
where it is stated) on the FIXED data path.

HIS RULES (his code):
  * ONE-HOUR blocks from each hour: data [Tstart, Tstart + 3600 + duration]; the correlation is the
    per-CHANNEL normalised cc of each of his templates; a channel enters the hour only if it has a
    COMPLETE record of it (his `len(subdata[0].data) == ndata`; here: every sample of the hour span
    recorded -- the cache's pre-fill valid mask -- and the channel not dead that day);
  * meancc = mean over those nchannel channels; MAD = median(|meancc - mean(meancc)|) of the HOUR;
    detections where meancc >= 8 x MAD, declustered by clean_LFEs (runs within 1 s -> max sample);
    time = Tstart + lag*dt (the template's first sample);
  * horizontals rotated to the template orientation when the data azimuth differs (his rotate_data);
  * catalogue (his analysis scripts): cc * nchannel >= threshold_perm[family]  -> column pass_his.
FIXES (duc_fix_clean = the sh_v2 fixes): flat-run / gap guard, burst muting, line notching, dead-channel
and dead-station floors, EXACT block-local window sums (no hour/day-long cumsum), relative energy guard
(window energy > 1e-3 x the channel's reference, scaled to the window length), and any |cc| > 1 left
after that is counted and zeroed.  A muted / flat / sub-guard window contributes cc = 0 but the channel
still counts in nchannel (his denominator), so muting can only LOWER a detection statistic.

Templates: his (STA, CHA) for epochs whose status is 'his'; for 'new' epochs, the re-templated
(fam, epoch, cha) from results/duc_fix/templates_new.{csv,npz} when kept, else the station is absent.

  python scripts/duc_fix_mf.py --t0 2008-04-01 --t1 2008-05-11 --out results/duc_fix/repro_mf.parquet
"""
import os, sys, time, argparse, warnings, pickle
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from duc_fix_common import *
import duc_fix_clean as CL
try:
    import cupy as cp; XP = cp; GPU = True
    cp.get_default_memory_pool().set_limit(size=12 * 1024 ** 3)      # share the GPU
except Exception:
    XP = np; GPU = False
def to_np(a): return cp.asnumpy(a) if GPU else a

THR = 8.0; FREQ0 = 1.0
NH = int(round(3600 * FS)); NDAY = 24 * NH
TAIL = int(round(60 * FS))          # next-day samples needed by hour 23 (max duration 45 s)
NREF = int(round(20.0 * FS))


def template_orient():
    O = {}
    base = f"{DUC}/src/templates_v1"
    for fam in os.listdir(base):
        for fn in os.listdir(f"{base}/{fam}"):
            if fn.endswith(".pkl") and fn[-5] in "EN12":
                k = fn[:-4]
                if k in O: continue
                tr, o = pickle.load(open(f"{base}/{fam}/{fn}", "rb")); O[tuple(k.split("_", 1))] = float(o["azimuth"])
    return O


def new_templates():
    p = f"{FIX}/templates_new.csv"
    if not os.path.exists(p): return {}
    M = pd.read_csv(p, dtype={"fam": str}); M = M[M.keep == 1]
    Z = np.load(f"{FIX}/templates_new.npz")
    return {(r.fam, r.epoch, r.cha): Z[r.key].astype(np.float64) for r in M.itertuples() if r.key in Z.files}


def decluster(idx, cc, maxdiff):
    out = []; run = [idx[0]]
    for i in idx[1:]:
        if i - run[-1] <= maxdiff: run.append(i)
        else: out.append(run); run = [i]
    out.append(run)
    return np.array([r[int(np.argmax(cc[r]))] for r in out])


class DayData:
    """cleaned day d (+ TAIL samples of d+1) for the stations requested; rotation applied.
    `cache` keeps DayPrep objects and loaded station-days across consecutive days."""
    def __init__(self, d, want, look, orient, cache):
        self.d = pd.Timestamp(d); self.ch = {}
        d1 = self.d + pd.Timedelta(days=1)
        for ns in want:
            sta = ns.split(".")[1]
            ep = look(sta, self.d)
            if ep is None: continue
            v = self._get(cache, ns, self.d)
            if v is None: continue
            t0, X, mc, ref, chs, rec = v
            if abs(t0 - self.d.timestamp()) > 1e-3: continue
            nxt = None
            ep1 = look(sta, d1)
            if ep1 is not None and ep1.epoch == ep.epoch:
                nxt = self._get(cache, ns, d1)
                if nxt is not None and abs(nxt[0] - d1.timestamp()) > 1e-3: nxt = None
            X, mc, rec, ref = self._rot(ep, chs, X, mc, rec, ref, orient, sta)
            if nxt is not None:
                Xn, mcn, recn, _ = self._rot(ep, nxt[4], nxt[1], nxt[2], nxt[5], nxt[3], orient, sta)
            N = X.shape[1]
            for i, c in enumerate(chs):
                x = np.zeros(NDAY + TAIL); m = np.zeros(NDAY + TAIL, bool); r = np.zeros(NDAY + TAIL, bool)
                n = min(N, NDAY); x[:n] = X[i, :n]; m[:n] = mc[i, :n]; r[:n] = rec[i, :n]
                if nxt is not None and c in nxt[4]:
                    j = nxt[4].index(c)
                    x[NDAY:] = Xn[j, :TAIL]; m[NDAY:] = mcn[j, :TAIL]; r[NDAY:] = recn[j, :TAIL]
                self.ch[(sta, c)] = (x, m, r, float(ref[i]), ep)

    @staticmethod
    def _get(cache, ns, d):
        k = ("L", ns, d)
        if k not in cache:
            P = cache.get(("P", d))
            if P is None:
                P = CL.DayPrep(d); cache[("P", d)] = P
            cache[k] = P.load(ns)
        return cache[k]

    @staticmethod
    def _rot(ep, chs, X, mc, rec, ref, orient, sta):
        """his rotate_data: horizontals to the template orientation when the data azimuth differs."""
        if ep.status != "his" or len(chs) != 3 or not isinstance(ep.az, str) or "/" not in ep.az:
            return X, mc, rec, ref
        cE, cN = chs[0], chs[1]
        tE, tN = orient.get((sta, cE)), orient.get((sta, cN))
        dE, dN = [float(v) for v in ep.az.split("/")]
        if tE is None or tN is None or (abs(dE - tE) < 0.5 and abs(dN - tN) < 0.5):
            return X, mc, rec, ref
        r = np.radians
        E, Nn = X[0].copy(), X[1].copy()
        X = X.copy(); mc = mc.copy(); rec = rec.copy(); ref = np.array(ref, float).copy()
        X[0] = np.cos(r(dE - tE)) * E + np.cos(r(dN - tE)) * Nn
        X[1] = np.cos(r(dE - tN)) * E + np.cos(r(dN - tN)) * Nn
        both = mc[0] & mc[1]; mc[0] = both; mc[1] = both
        rb = rec[0] & rec[1]; rec[0] = rb; rec[1] = rb
        rr = 0.0 if min(ref[0], ref[1]) <= 0 else 0.5 * (ref[0] + ref[1]); ref[0] = rr; ref[1] = rr
        X = np.where(mc, X, 0.0)
        return X, mc, rec, ref


def run(days, out, fams_only=None, diag_out=None, exclude=()):
    warnings.filterwarnings("ignore")
    fp = families(); thr = thresholds(); HT = his_templates(); NT = new_templates()
    orient = template_orient(); look = epoch_lookup()
    H = his_rows(); NETOF = dict(zip(H.sta, H.net))
    famst = {r.family: r.stations.split(",") for r in fp.itertuples()}
    dur = {r.family: float(r.duration) for r in fp.itertuples()}
    fams = [f for f in fp.family if (fams_only is None or f in fams_only)]
    want = sorted({f"{NETOF[s]}.{s}" for f in fams for s in famst[f] if s in NETOF and s not in set(exclude)})
    maxdiff = int(np.floor(1.0 / ((1 / FS) * FREQ0)))
    rows = []; diag = []; t00 = time.time(); cache = {}
    for d in days:
        d = pd.Timestamp(d)
        for k in [k for k in cache if k[-1] < d]: cache.pop(k)
        DD = DayData(d, want, look, orient, cache)
        if not DD.ch: continue
        nfft = int(2 ** np.ceil(np.log2(NDAY + TAIL + 1000)))
        G = {}        # (sta, cha) -> dict(F, x, m, rec_ok per hour, ref)
        def gch(k):
            if k in G: return G[k]
            x, m, r, ref, ep = DD.ch[k]
            xg = XP.asarray(x)
            # hour participation: every sample of [h*NH, h*NH + NH + L) recorded (L added per family)
            cr = np.concatenate([[0], np.cumsum(~r)])
            G[k] = dict(F=XP.fft.rfft(xg, nfft), x=xg, x2=xg * xg, m=XP.asarray(m), cr=cr, ref=ref, ep=ep, S={})
            return G[k]
        lastL = None
        for f in sorted(fams, key=lambda f_: dur[f_]):
            L = int(round(dur[f] * FS)) + 1
            if L != lastL:
                for g_ in G.values(): g_["S"].clear()
                lastL = L
            # templates available today for this family
            TT = []
            for sta in famst[f]:
                for (s2, c), (x, m, r, ref, ep) in DD.ch.items():
                    if s2 != sta: continue
                    if ep.status == "his": tp = HT.get(f, {}).get((sta, c))
                    else: tp = NT.get((f, ep.epoch, c))
                    if tp is None or tp.size != L or not np.all(np.isfinite(tp)) or not np.any(tp): continue
                    TT.append(((sta, c), tp))
            if not TT: continue
            nw = NDAY + 1            # cc samples 0 .. 24*NH (the last one = hour 23's i = NH)
            SUM = XP.zeros((24, NH + 1)); LIVE = XP.zeros((24, NH + 1)); NCH = np.zeros(24)
            ngt1 = 0
            hidx = XP.asarray((np.arange(24)[:, None] * NH + np.arange(NH + 1)[None, :]))
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
            if not NCH.any(): continue
            SUMn = to_np(SUM); LIVEn = to_np(LIVE)
            thrp = float(thr.get(f, 0.0))
            ndet = 0; npass = 0; mads = []
            for h in range(24):
                if NCH[h] == 0: continue
                mcc = SUMn[h] / NCH[h]
                mad = float(np.median(np.abs(mcc - mcc.mean()))); mads.append(mad)
                idx = np.where(mcc >= THR * mad)[0]
                if idx.size == 0 or mad <= 0: continue
                keep = decluster(idx, mcc, maxdiff)
                t0h = d.timestamp() + h * 3600.0
                for i in keep:
                    c_ = float(mcc[i]); nc = int(NCH[h])
                    p = int(thrp > 0 and c_ * nc >= thrp); npass += p; ndet += 1
                    rows.append((f, t0h + i / FS, c_, nc, int(LIVEn[h, i]), mad, c_ * nc, thrp, p))
            diag.append(dict(fam=f, day=str(d.date()), n_tpl=len(TT), nch_mean=float(NCH[NCH > 0].mean()),
                             hours=int((NCH > 0).sum()), mad_med=float(np.median(mads)) if mads else np.nan,
                             n_ccgt1=ngt1, n_det=ndet, n_pass=npass,
                             new_tpl=int(sum(1 for k, _ in TT if DD.ch[k][4].status == "new"))))
        del G
        if GPU: cp.get_default_memory_pool().free_all_blocks()
    D = pd.DataFrame(rows, columns=["fam", "t_epoch", "cc", "nchannel", "nlive", "mad", "ccsum", "thr_perm", "pass_his"])
    D["t"] = pd.to_datetime(D.t_epoch, unit="s").astype("datetime64[ns]")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)
    D.to_parquet(out)
    dg = pd.DataFrame(diag); dg.to_csv(diag_out or out.replace(".parquet", "_diag.csv"), index=False)
    log(f"[duc-mf-fix] {len(days)} days: {len(D)} detections (pass_his {int(D.pass_his.sum()) if len(D) else 0}), "
        f"family-days {len(dg)}, |cc|>1 zeroed {int(dg.n_ccgt1.sum()) if len(dg) else 0}, "
        f"{time.time() - t00:.0f}s -> {out}")
    return D


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--t0"); ap.add_argument("--t1")
    ap.add_argument("--month", default=""); ap.add_argument("--out", required=True)
    ap.add_argument("--fams", default=""); ap.add_argument("--exclude", default="")
    a = ap.parse_args()
    if a.month:
        m0 = pd.Timestamp(a.month + "-01"); days = pd.date_range(m0, m0 + pd.offsets.MonthBegin(1), freq="D", inclusive="left")
    else:
        days = pd.date_range(a.t0, a.t1, freq="D", inclusive="left")
    run(days, a.out, set(a.fams.split(",")) if a.fams else None, exclude=tuple(x for x in a.exclude.split(",") if x))


if __name__ == "__main__":
    main()
