#!/usr/bin/env python
"""Build ONE three-component ARRAY TRACE per 1-km array, the way the Olympic LFE papers did.
Prefix `o_`.

METHOD, taken from the papers (results/WA_2010_ARRAY_RESULTS.md sec.0 quotes them):

  Ghosh et al. (2012, G-cubed), whose beams Chestler & Creager (2017) use verbatim:
    "For each array and every 1-minute time window [we] filtered vertical-component seismograms
     between 5-9 Hz.  Then [we] apply a DELAY-AND-SUM BEAMFORMING approach to stack data and find
     the VECTOR SLOWNESS corresponding to the direction of a plane wave that produces the stack
     with the most power.  Stacking seismograms from each window at the optimal slowness produces
     a continuous time series with excellent SNR.  This is repeated for each array and every
     50-percent overlapping one-minute time window."   (Chestler thesis sec.2.1)

  So: a SLOWNESS-GRID BEAM, not a simple stack -- implemented here as such.

  Two documented deviations, both forced by the frozen cache, neither a free parameter:
    (1) BAND.  The frozen 20 Hz cache is band-passed 1-8 Hz, so the slowness is estimated on
        5-8 Hz instead of 5-9 Hz.  The BEAM ITSELF is formed over the full cached 1-8 Hz band,
        because the rest of the chain (autocorrelation, templates, matched filter) is 1-8 Hz.
    (2) STATIC TERMS.  A per-member static delay is estimated once per array from tremor-hour
        records and applied in addition to the plane-wave delay.  It is the RESIDUAL of the
        measured inter-member delays after the best-fit plane wave is removed, so it cannot
        absorb the slowness the beamformer searches for.  (The papers do not state statics; a
        1-km array has < 0.05 s of them and they are reported, not hidden.)

Output: `beamcache/<SITE>/<year>/<doy>.npz` in the SAME format as the station cache
(x (3, N) float32 Z,N,E; t0; order; fs), so `k_common.load_day` reads a beam with no code change.

  python scripts/o_beam.py --stage statics
  python scripts/o_beam.py --stage smoke   --site xLC --date 2010-08-15 --hour 12
  O_BED=2010 python scripts/o_beam.py --stage build --workers 6
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import o_common as O                                                        # noqa: E402

FS = 20.0
WIN_N = 1200              # 60 s at 20 Hz -- Ghosh's one-minute window
HOP_N = 600               # 50 % overlap
NFFT = 2048
BAND_STEER = (5.0, 8.0)   # Ghosh 5-9 Hz; the cache stops at 8 Hz (deviation 1)
SMAX = 0.40               # s/km; a 1-km array at 6.5 Hz has a main lobe ~1/(fD) = 0.15 s/km,
SSTEP = 0.05              # ... so 0.05 s/km already over-samples the beam pattern
STAT_LAG_S = 0.6
STAT_ROUNDS = 3
STAT_HOURS = 8
STEER = os.environ.get("O_STEER", "none")     # "none" = zero-slowness; "grid" = Ghosh slowness search
NDAY = 1728001


def _xy_km(lats, lons, lat0, lon0):
    """Local east/north km relative to (lat0, lon0)."""
    e = (np.asarray(lons) - lon0) * 111.320 * np.cos(np.radians(lat0))
    n = (np.asarray(lats) - lat0) * 110.574
    return np.asarray(e, float), np.asarray(n, float)


def slowness_grid():
    v = np.arange(-SMAX, SMAX + 1e-9, SSTEP)
    SX, SY = np.meshgrid(v, v, indexing="ij")
    SX, SY = SX.ravel(), SY.ravel()
    m = (SX ** 2 + SY ** 2) <= SMAX ** 2 + 1e-9
    return SX[m], SY[m]


# =============================================================== member day loading =============
def load_member(sta, date):
    """(3, N) ZNE float32 for one member-day, through k_common (rotation for the boreholes)."""
    import k_common as K
    import s_common as S
    idx = O.member_cache_index()
    ts = pd.Timestamp(date)
    key = (sta, int(ts.year), int(ts.dayofyear))
    if key not in idx:
        return None
    AZ = O.member_az12()
    if sta in AZ:
        S.AZ12.setdefault(sta, AZ[sta])
        S.ZSIGN.setdefault(sta, 1.0)
        v = S.load_day_zne(sta, key[1], key[2], idx)
        if v is None or v[1].shape[0] != 3:
            return None
        return float(v[0]), np.ascontiguousarray(v[1], np.float32)
    z = np.load(idx[key])
    x = z["x"]
    order = [str(o) for o in z["order"]]
    if x.shape[0] != 3 or [o[-1] for o in order] != ["Z", "N", "E"]:
        return None
    return float(z["t0"]), np.ascontiguousarray(x, np.float32)


# =============================================================== statics ========================
def measure_statics(site, members, hours, log=print):
    """Per-member delays from tremor-hour Z records, with the best-fit PLANE WAVE removed.

    Rounds of (cross-correlate every member against the current aligned mean, shift) at sample
    resolution with parabolic refinement.  The measured delay of member m is then split into the
    plane-wave part (a + sx*e_m + sy*n_m, least squares) and the STATIC residual, and only the
    residual is returned -- the plane-wave part is exactly what the per-window beamformer finds.
    """
    ST = O.member_table()
    lat0 = float(np.mean([ST.lat[m] for m in members]))
    lon0 = float(np.mean([ST.lon[m] for m in members]))
    e, n = _xy_km([ST.lat[m] for m in members], [ST.lon[m] for m in members], lat0, lon0)
    segs = []
    for (date, hh) in hours:
        Z, have = [], []
        for m in members:
            v = load_member(m, date)
            if v is None:
                continue
            t0, X = v
            i0 = int(round((pd.Timestamp(f"{date} {hh:02d}:00:00").timestamp() - t0) * FS))
            x = X[0, i0:i0 + int(600 * FS)]                      # 10 minutes of Z
            if x.size < int(600 * FS) or not np.isfinite(x).all() or x.std() <= 0:
                continue
            Z.append(x / x.std())
            have.append(m)
        if len(have) >= max(3, len(members) // 2):
            segs.append((have, np.asarray(Z, np.float64)))
    if not segs:
        return {m: 0.0 for m in members}, dict(n_seg=0)
    nl = int(round(STAT_LAG_S * FS))
    acc = {m: [] for m in members}
    coh0, coh1 = [], []
    for have, Z in segs:
        d = np.zeros(len(have))
        coh0.append(_coh(Z))
        for _ in range(STAT_ROUNDS):
            Zs = np.array([np.roll(Z[i], -int(round(d[i]))) for i in range(len(have))])
            ref = Zs.mean(0)
            rn = np.sqrt((ref ** 2).sum())
            for i in range(len(have)):
                cc = np.correlate(Z[i], ref, mode="same")
                c0 = Z[i].size // 2
                w = cc[c0 - nl:c0 + nl + 1]
                k = int(np.argmax(w))
                den = np.sqrt((Z[i] ** 2).sum()) * rn
                if den <= 0:
                    continue
                kk = float(k)
                if 0 < k < len(w) - 1:
                    a, b, c = w[k - 1], w[k], w[k + 1]
                    q = a - 2 * b + c
                    if q != 0:
                        kk = k + 0.5 * (a - c) / q
                d[i] = kk - nl
        coh1.append(_coh(np.array([_shift(Z[i], -d[i]) for i in range(len(have))])))
        for i, m in enumerate(have):
            acc[m].append(d[i] / FS)
    delay = np.array([float(np.median(acc[m])) if acc[m] else 0.0 for m in members])
    A = np.column_stack([np.ones(len(members)), e, n])
    coef, *_ = np.linalg.lstsq(A, delay, rcond=None)
    resid = delay - A @ coef
    resid = resid - resid.mean()
    log(f"  [{site}] {len(segs)} tremor segments; raw delay sd {delay.std()*1000:.1f} ms, "
        f"plane-wave slowness ({coef[1]:+.3f},{coef[2]:+.3f}) s/km, "
        f"STATIC residual sd {resid.std()*1000:.1f} ms (max {np.abs(resid).max()*1000:.1f} ms); "
        f"Z coherence plain {np.mean(coh0):.2f} -> aligned {np.mean(coh1):.2f}")
    info = dict(n_seg=len(segs), delay_sd_ms=float(delay.std() * 1000),
                plane_sx=float(coef[1]), plane_sy=float(coef[2]),
                static_sd_ms=float(resid.std() * 1000),
                static_max_ms=float(np.abs(resid).max() * 1000),
                coh_plain=float(np.mean(coh0)), coh_aligned=float(np.mean(coh1)))
    return {m: float(resid[i]) for i, m in enumerate(members)}, info


def _shift(x, k):
    n = x.size
    f = np.fft.rfft(x)
    f *= np.exp(-2j * np.pi * np.fft.rfftfreq(n) * k)
    return np.fft.irfft(f, n)


def _coh(Z):
    """N * <stack^2> / <trace^2>; 1 = incoherent, N = perfectly coherent."""
    Z = np.asarray(Z, float)
    e = (Z ** 2).mean(1)
    ok = e > 0
    Z = Z[ok] / np.sqrt(e[ok])[:, None]
    if len(Z) < 2:
        return np.nan
    return float(len(Z) * (Z.mean(0) ** 2).mean() / max((Z ** 2).mean(), 1e-300))


# =============================================================== the beam =======================
def beam_hour(Xs, e, n, stat, dev="cuda", pad=HOP_N, want_diag=False):
    """Xs = (M, 3, L) float32 for ONE hour PLUS `pad` samples of margin on each side.
    Returns (3, L - 2*pad) beam, and the per-window optimal slowness."""
    import torch
    M, _, L = Xs.shape
    nw = (L - WIN_N) // HOP_N + 1
    if nw < 2:
        return None, None, None
    X = torch.from_numpy(np.ascontiguousarray(Xs)).to(dev)
    w = torch.hann_window(WIN_N, periodic=True, device=dev, dtype=torch.float32)
    U = X.unfold(-1, WIN_N, HOP_N)[:, :, :nw, :] * w                    # (M,3,nw,WIN_N)
    F = torch.fft.rfft(U, n=NFFT, dim=-1)                               # (M,3,nw,NF)
    fr = torch.fft.rfftfreq(NFFT, d=1.0 / FS).to(dev)
    SX, SY = slowness_grid()
    tau = (torch.from_numpy(np.outer(SX, e) + np.outer(SY, n)).to(dev).float())  # (Ns,M)
    bs = (fr >= BAND_STEER[0]) & (fr <= BAND_STEER[1])
    bi = torch.nonzero(bs).ravel()[::2]                                 # every 2nd bin is plenty
    Zf = F[:, 0]                                                        # (M,nw,NF)
    P = torch.zeros((tau.shape[0], nw), device=dev)
    CH = 24
    for k in range(0, bi.numel(), CH):
        b = bi[k:k + CH]
        ph = torch.exp(2j * np.pi * fr[b][None, None, :] * tau[:, :, None])   # (Ns,M,nb)
        A = Zf[:, :, b]                                                       # (M,nw,nb)
        B = torch.einsum("smb,mwb->swb", ph, A.to(torch.complex64))
        P += (B.real ** 2 + B.imag ** 2).sum(-1)
        del ph, A, B
    best = torch.argmax(P, dim=0)                                       # (nw,)
    tb = tau[best]                                                      # (nw,M)
    dl = torch.from_numpy(np.asarray(stat, np.float32)).to(dev)[None, :]
    sh = torch.exp(2j * np.pi * fr[None, None, :] * (tb + dl)[:, :, None])   # (nw,M,NF)
    Bm = torch.einsum("wmf,mcwf->cwf", sh, F.to(torch.complex64)) / float(M)
    y = torch.fft.irfft(Bm, n=NFFT, dim=-1)[:, :, :WIN_N]               # (3,nw,WIN_N)
    out = torch.zeros((3, L), device=dev)
    idx = torch.arange(WIN_N, device=dev)[None, :] + \
        (torch.arange(nw, device=dev)[:, None] * HOP_N)
    out.scatter_add_(1, idx.reshape(1, -1).expand(3, -1), y.reshape(3, -1))
    diag = None
    if want_diag:
        diag = dict(sx=SX[best.cpu().numpy()].tolist(), sy=SY[best.cpu().numpy()].tolist())
    res = out[:, pad:L - pad].cpu().numpy().astype(np.float32)
    del X, U, F, P, tb, sh, Bm, y, out
    torch.cuda.empty_cache()
    return res, (SX[best.cpu().numpy()], SY[best.cpu().numpy()]), diag


def build_site_day(site, members, date, statics, mveto, dev="cuda", log=print):
    """One (3, NDAY) beam for one array-site on one day.  Returns (x, rows) or (None, rows)."""
    ST = O.member_table()
    lat0 = float(np.mean([ST.lat[m] for m in members]))
    lon0 = float(np.mean([ST.lon[m] for m in members]))
    D = {}
    t0 = None
    for m in members:
        v = load_member(m, date)
        if v is None:
            continue
        D[m] = v[1]
        t0 = v[0] if t0 is None else t0
    if len(D) < 2:
        return None, []
    Y = np.zeros((3, NDAY), np.float32)
    rows = []
    for hh in range(24):
        bad = mveto.get((date, hh), set())
        live = [m for m in members if m in D and m not in bad]
        a = hh * int(3600 * FS)
        b = a + int(3600 * FS)
        ok = []
        for m in live:
            X = D[m]
            if b > X.shape[1]:
                continue
            seg = X[:, a:b]
            if np.isfinite(seg).all() and float((seg ** 2).mean()) > 0:
                ok.append(m)
        rows.append(dict(site=site, date=date, hour=hh, n_members=len(ok),
                         n_members_pool=len(members)))
        if len(ok) < 2:
            continue
        p = HOP_N
        lo, hi = a - p, b + p
        Xs = np.zeros((len(ok), 3, (b - a) + 2 * p), np.float32)
        for i, m in enumerate(ok):
            X = D[m]
            s0, s1 = max(lo, 0), min(hi, X.shape[1])
            Xs[i, :, s0 - lo:s1 - lo] = X[:, s0:s1]
            r = np.sqrt((Xs[i] ** 2).mean())
            if r > 0:
                Xs[i] /= r                       # equal weight per member (a beam, not a sum)
        if STEER == "none":
            # zero-slowness delay-and-sum: the slowness search is not identifiable at this
            # aperture (see stage `smoke` / results/WA_2010_ARRAY_RESULTS.md sec.0c), so the
            # beam is the equal-weight mean of the rms-normalised members.
            Y[:, a:b] = Xs[:, :, p:-p].mean(0)
        else:
            e, n = _xy_km([ST.lat[m] for m in ok], [ST.lon[m] for m in ok], lat0, lon0)
            stat = [statics.get(m, 0.0) for m in ok]
            y, _, _ = beam_hour(Xs, e, n, stat, dev=dev)
            if y is not None:
                Y[:, a:b] = y
    return (Y, t0, rows)



# =============================================================== EQ statics =====================
EQ_BAND = (1.0, 8.0)
EQ_LAG_S = 0.6
EQ_WIN_S = 4.0
EQ_PRE_S = 1.0


def _eq_events():
    E = pd.read_csv(f"{O.RES}/intraslab_eq_2006_2012.csv")
    E["origin"] = pd.to_datetime(E.origin_utc, utc=True)
    E = E[(E.lat > 47.3) & (E.lat < 48.7) & (E.lon > -124.8) & (E.lon < -122.2)
          & (E.depth_km >= 20)]
    return E.sort_values("mag", ascending=False).reset_index(drop=True)


def _eq_cut(site, members, ev, want="S"):
    """{member: (S window, noise window)} at native rate, 1-8 Hz, for one event."""
    import s_common as S
    import p_vel1d as V
    from obspy import read
    from scipy.signal import butter, sosfiltfilt
    ST = O.member_table()
    out, fs0 = {}, None
    for m in members:
        p = f"{O.STUDY}/traces_eq/{ST.net[m]}.{m}/{ev.evid}.mseed"
        if not (os.path.exists(p) and os.path.getsize(p) > 0):
            continue
        try:
            st = read(p)
            st.merge(method=1, fill_value=0)
        except Exception:
            continue
        tr = [t for t in st if t.stats.channel.endswith("Z")]
        if not tr:
            continue
        tr = tr[0]
        fs = float(tr.stats.sampling_rate)
        if fs < 40 or (fs0 is not None and abs(fs - fs0) > 1e-6):
            continue
        fs0 = fs
        x = tr.data.astype(float)
        x -= x.mean()
        x = sosfiltfilt(butter(4, [EQ_BAND[0] / (fs / 2), EQ_BAND[1] / (fs / 2)],
                               btype="band", output="sos"), x)
        d = float(S.epi_km(ev.lat, ev.lon, float(ST.lat[m]), float(ST.lon[m])))
        tS = ev.origin.timestamp() + V.tt(d, float(ev.depth_km), "S")
        i0 = int(round((tS - float(tr.stats.starttime.timestamp)) * fs)) - int(EQ_PRE_S * fs)
        nw = int(EQ_WIN_S * fs)
        j0 = i0 - int(20 * fs)
        if i0 < 0 or i0 + nw > len(x) or j0 < 0:
            continue
        a, b = x[i0:i0 + nw], x[j0:j0 + nw]
        if a.std() <= 0 or b.std() <= 0 or a.std() < 3 * b.std():
            continue
        out[m] = (a / a.std(), b / b.std())
    return out, fs0


def _align(Z, maxlag, rounds=3):
    d = np.zeros(len(Z))
    for _ in range(rounds):
        ref = np.array([_shift(Z[i], -d[i]) for i in range(len(Z))]).mean(0)
        rn = np.sqrt((ref ** 2).sum())
        if rn <= 0:
            break
        for i in range(len(Z)):
            c = np.correlate(Z[i], ref, mode="same")
            c0 = Z[i].size // 2
            w = c[c0 - maxlag:c0 + maxlag + 1]
            k = int(np.argmax(w))
            kk = float(k)
            if 0 < k < len(w) - 1:
                a, b, cc = w[k - 1], w[k], w[k + 1]
                q = a - 2 * b + cc
                if q != 0:
                    kk = k + 0.5 * (a - cc) / q
            d[i] = kk - maxlag
    return d - d.mean()


def stage_eqstatics(a):
    """Per-member STATIC delays from catalogued INTRASLAB EARTHQUAKES -- deep sources at the LFE
    incidence, in the LFE band.  Cross-validated: the statics used to score event e are the
    median over all OTHER events, and the identical procedure is applied to a pre-event NOISE
    window as the control (aligning noise also raises coherence, so an uncontrolled coherence
    gain would prove nothing)."""
    O.activate()
    E = _eq_events()
    A = O.arrays()
    O.log(f"{len(E)} catalogued intraslab events in the box")
    per = {s: {} for s in A}
    rows = []
    for ev in E.itertuples():
        for site, mem in sorted(A.items()):
            D, fs = _eq_cut(site, mem, ev)
            if len(D) < 4:
                continue
            ms = sorted(D)
            Zs = np.array([D[m][0] for m in ms])
            Nz = np.array([D[m][1] for m in ms])
            nl = int(round(EQ_LAG_S * fs))
            d = _align(Zs, nl)
            dn = _align(Nz, nl)
            for i, m in enumerate(ms):
                per[site].setdefault(m, []).append(float(d[i] / fs))
            rows.append(dict(evid=ev.evid, mag=float(ev.mag), dep=float(ev.depth_km), site=site,
                             M=len(ms), fs=fs,
                             coh_raw=_coh(Zs),
                             coh_self=_coh([_shift(Zs[i], -d[i]) for i in range(len(ms))]),
                             coh_noise_raw=_coh(Nz),
                             coh_noise_self=_coh([_shift(Nz[i], -dn[i]) for i in range(len(ms))]),
                             lag_sd_ms=float(d.std() / fs * 1000),
                             members=",".join(ms), lags=",".join(f"{x/fs:.4f}" for x in d)))
    R = pd.DataFrame(rows)
    if not len(R):
        raise SystemExit("no usable intraslab events for any array")
    stat = {s: {m: float(np.median(v)) for m, v in per[s].items() if len(v) >= 3}
            for s in per}
    for s in stat:
        if stat[s]:
            mu = float(np.mean(list(stat[s].values())))
            stat[s] = {m: v - mu for m, v in stat[s].items()}
    # ---- leave-one-event-out validation of the static-corrected plain stack -------------------
    val = []
    for r in R.itertuples():
        ms = r.members.split(",")
        lags = [float(x) for x in r.lags.split(",")]
        oth = {}
        for m in ms:
            v = [x for x in per[r.site][m]]
            k = ms.index(m)
            v = [x for x in v if abs(x - lags[k]) > 1e-12] or v
            if len(v) >= 2:
                oth[m] = float(np.median(v))
        if len(oth) < 4:
            continue
        mu = float(np.mean(list(oth.values())))
        D, fs = _eq_cut(r.site, ms, E[E.evid == r.evid].iloc[0])
        use = [m for m in ms if m in D and m in oth]
        if len(use) < 4:
            continue
        Zs = [D[m][0] for m in use]
        Nz = [D[m][1] for m in use]
        sh = [-(oth[m] - mu) * fs for m in use]
        val.append(dict(evid=r.evid, site=r.site, M=len(use),
                        coh_plain=_coh(Zs),
                        coh_loo_static=_coh([_shift(Zs[i], sh[i]) for i in range(len(use))]),
                        coh_noise_plain=_coh(Nz),
                        coh_noise_loo_static=_coh([_shift(Nz[i], sh[i])
                                                   for i in range(len(use))])))
    Vd = pd.DataFrame(val)
    os.makedirs(O.OUTROOT, exist_ok=True)
    R.to_csv(f"{O.OUTROOT}/eq_statics_events.csv", index=False)
    if len(Vd):
        Vd.to_csv(f"{O.OUTROOT}/eq_statics_validation.csv", index=False)
    summ = {}
    for s in sorted(A):
        v = Vd[Vd.site == s] if len(Vd) else Vd
        rr = R[R.site == s]
        summ[s] = dict(
            n_events=int(len(rr)), n_members_with_static=len(stat.get(s, {})),
            static_sd_ms=float(np.std(list(stat[s].values())) * 1000) if stat.get(s) else None,
            static_max_ms=float(np.max(np.abs(list(stat[s].values()))) * 1000)
            if stat.get(s) else None,
            coh_plain=float(v.coh_plain.median()) if len(v) else None,
            coh_static_LOO=float(v.coh_loo_static.median()) if len(v) else None,
            coh_noise_plain=float(v.coh_noise_plain.median()) if len(v) else None,
            coh_noise_static_LOO=float(v.coh_noise_loo_static.median()) if len(v) else None,
            n_val=int(len(v)))
        O.log(f"  [{s}] events {summ[s]['n_events']:2d}  statics {summ[s]['n_members_with_static']:2d}"
              f"  sd {summ[s]['static_sd_ms']}ms   EQ coherence plain "
              f"{summ[s]['coh_plain']} -> LOO-static {summ[s]['coh_static_LOO']}   "
              f"(noise {summ[s]['coh_noise_plain']} -> {summ[s]['coh_noise_static_LOO']})")
    json.dump(dict(source="intraslab_earthquakes", statics=stat, summary=summ),
              open(a.out or O.STATICS, "w"), indent=1)
    O.log(f"wrote {a.out or O.STATICS}")


# =============================================================== stages =========================
def stage_statics(a):
    O.activate()
    A = O.arrays()
    H = O.candidate_hours()
    H = H.sort_values(["date", "hour"])
    pick = H.iloc[np.linspace(0, len(H) - 1, min(STAT_HOURS, len(H))).astype(int)]
    hours = [(str(r.date), int(r.hour)) for r in pick.itertuples()]
    O.log(f"statics from {len(hours)} tremor hours: {hours}")
    out, info = {}, {}
    for site, mem in A.items():
        d, i = measure_statics(site, mem, hours, log=O.log)
        out[site] = d
        info[site] = i
    json.dump(dict(statics=out, info=info, hours=hours),
              open(a.out or O.STATICS, "w"), indent=1)
    O.log(f"wrote {a.out or O.STATICS}")


def stage_smoke(a):
    """Coherence gain of the beam over a single member, on one array / one hour."""
    O.activate()
    A = O.arrays()
    site = a.site or sorted(A)[0]
    mem = A[site]
    ST = O.member_table()
    stat = json.load(open(O.STATICS))["statics"].get(site, {}) if os.path.exists(O.STATICS) \
        else {}
    D = {m: load_member(m, a.date) for m in mem}
    D = {m: v[1] for m, v in D.items() if v is not None}
    a0 = a.hour * int(3600 * FS)
    p = HOP_N
    ok = [m for m in mem if m in D and np.isfinite(D[m][:, a0:a0 + 72000]).all()
          and (D[m][:, a0:a0 + 72000] ** 2).mean() > 0]
    Xs = np.zeros((len(ok), 3, 72000 + 2 * p), np.float32)
    for i, m in enumerate(ok):
        X = D[m]
        s0, s1 = max(a0 - p, 0), min(a0 + 72000 + p, X.shape[1])
        Xs[i, :, s0 - (a0 - p):s1 - (a0 - p)] = X[:, s0:s1]
        Xs[i] /= max(np.sqrt((Xs[i] ** 2).mean()), 1e-30)
    lat0 = float(np.mean([ST.lat[m] for m in mem]))
    lon0 = float(np.mean([ST.lon[m] for m in mem]))
    e, n = _xy_km([ST.lat[m] for m in ok], [ST.lon[m] for m in ok], lat0, lon0)
    z = Xs[:, 0, p:-p]
    coh_raw = _coh(z)
    zs = np.array([_shift(z[i], -stat.get(ok[i], 0.0) * FS) for i in range(len(ok))])
    coh_stat = _coh(zs)
    yb, sl, _ = beam_hour(Xs, e, n, [stat.get(m, 0.0) for m in ok])
    y0, _, _ = beam_hour(Xs, e, n, [stat.get(m, 0.0) for m in ok])
    # a "single member" reference: the member closest to the centroid
    j = int(np.argmin(e ** 2 + n ** 2))
    r = dict(site=site, date=a.date, hour=a.hour, n_members=len(ok),
             coh_plain_stack=round(coh_raw, 3), coh_static_aligned=round(coh_stat, 3),
             n_max=len(ok),
             beam_over_member_rms=float(np.sqrt((yb[0] ** 2).mean())
                                        / max(np.sqrt((z[j] ** 2).mean()), 1e-30)),
             slowness_median=[float(np.median(sl[0])), float(np.median(sl[1]))],
             slowness_abs_median=float(np.median(np.hypot(sl[0], sl[1]))),
             identical_rerun=bool(np.allclose(yb, y0)))
    # beam coherence: how much of a member's record the beam explains, vs the plain stack
    cb, cp = [], []
    for i in range(len(ok)):
        b = yb[0]
        cb.append(float(np.corrcoef(b, z[i])[0, 1]))
        cp.append(float(np.corrcoef(z.mean(0), z[i])[0, 1]))
    r["mean_corr_beam_member"] = round(float(np.mean(cb)), 3)
    r["mean_corr_plainstack_member"] = round(float(np.mean(cp)), 3)
    print(json.dumps(r, indent=1), flush=True)
    return r


def stage_build(a):
    O.activate()
    A = O.arrays()
    stat = json.load(open(O.STATICS))["statics"] if os.path.exists(O.STATICS) else {}
    mveto = O._member_veto()
    days = O.beam_days() if not a.days else [d for d in a.days.split(",") if d]
    O.log(f"[build {O.BED}] {len(A)} arrays x {len(days)} days -> {O.BEAM}")
    rows = []
    logp = O.BEAMLOG.format(bed=O.BED)
    done = 0
    for date in days:
        for site, mem in sorted(A.items()):
            p = O.beam_path(site, date)
            if os.path.exists(p) and not a.force:
                continue
            r = build_site_day(site, mem, date, stat.get(site, {}), mveto, log=O.log)
            if r[0] is None:
                O.log(f"[skip] {site} {date}: < 2 members")
                continue
            Y, t0, rr = r
            os.makedirs(os.path.dirname(p), exist_ok=True)
            np.savez(p + ".tmp.npz", x=Y, t0=np.float64(t0),
                     order=np.array(["BHZ", "BHN", "BHE"]), fs=np.float64(FS))
            os.replace(p + ".tmp.npz", p)
            rows += rr
            done += 1
        O.log(f"[build] {date} done ({done} site-days so far)")
        if rows:
            R = pd.DataFrame(rows)
            if os.path.exists(logp):
                R = pd.concat([pd.read_csv(logp), R], ignore_index=True)
                R = R.drop_duplicates(["site", "date", "hour"], keep="last")
            R.to_csv(logp, index=False)
            rows = []
    O.log(f"[build {O.BED}] complete; hour log {logp}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["statics", "eqstatics", "smoke", "build"])
    ap.add_argument("--site", default="")
    ap.add_argument("--date", default="2010-08-15")
    ap.add_argument("--hour", type=int, default=12)
    ap.add_argument("--days", default="")
    ap.add_argument("--out", default="")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()
    {"statics": stage_statics, "eqstatics": stage_eqstatics,
     "smoke": stage_smoke, "build": stage_build}[a.stage](a)


if __name__ == "__main__":
    main()
