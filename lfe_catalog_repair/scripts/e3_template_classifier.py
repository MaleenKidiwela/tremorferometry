#!/usr/bin/env python
"""E3 step 4 — Bostock family templates as a CLASSIFIER on the associated events.

The question: do the PyOcto-associated events look like KNOWN LFE families, more often than the
shifted-pick chance events of step 1 do?  This is a classifier on ~100-300 candidates/day, not a
detector on 86,400 s/day, so D3's false-alarm-floor arithmetic does not apply.

Network cross-correlation.  For a candidate origin time T and a family f, the template at station
s is the family's own SeisBench stack, the segment from S-1 s to S+7 s (160 samples at 20 Hz) --
the stacks are P-aligned at sample 800, so the segment's offset from the family origin is
ttP(f,s) + (S_sample - 800)/20, i.e. the family's own stack S-P plus the Crosson-1D P moveout to
that station.  The 3 components are concatenated into one 480-sample vector, L2-normalised, and
dotted with the identically normalised data window; the network cc is the mean over the stations
that have data.  A single common lag (the origin-time shift) is scanned over +-3 s at 0.05 s, so
the moveout is held fixed and only the origin time floats -- which is what makes this a test of
network coherence and not of 6 independent single-station correlations.

Event sets
  assoc  vi2010 cascadia thr 0.3, N >= 4
  null   the step-1 shifted-pick chance events, same threshold and N cut
  bostock  Bostock's own catalogued events in the window, at the equivalent origin time from the
           PER-FAMILY offsets of step 5, each scored against ITS OWN family template -- this is
           the "how many MADs does Bostock's own event register with his own template on these six
           stations" measurement that decides whether D3's 3.5 MAD was a template-quality artifact
  random   uniform random times, the per-template per-day background whose median and MAD
           normalise everything else

Outputs results/e3_classifier_events.parquet, results/e3_classifier_summary.csv
"""
from __future__ import annotations

import argparse
import os
import sys

import h5py
import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
TRACES = "lfe_catalog_repair/traces"
STACKDIR = f"{RES}/external/other/seisbench_lfestacks_cascadia"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import d3_vel1d as v1d  # noqa: E402
import e4_chanorder as chanorder  # noqa: E402

LEGACY_CHAN = bool(int(os.environ.get("E3_LEGACY_CHAN_ORDER", "0")))  # 1 = reproduce the pre-fix bug
_ORDER_SEEN = {}

OVERLAP = ["PGC", "LZB", "VGZ", "YOUB", "SNB", "PFB"]
NETSTA = {"PGC": "CN.PGC", "LZB": "CN.LZB", "VGZ": "CN.VGZ", "YOUB": "CN.YOUB",
          "SNB": "CN.SNB", "PFB": "CN.PFB"}
FS = 20.0
TL = 160            # template samples = 8 s at 20 Hz (S-1 s .. S+7 s)
PRE = 1.0           # template starts PRE seconds before S
LAGS = np.arange(-3.0, 3.0001, 0.05)


def bandpass(x, fs, f1=1.0, f2=8.0):
    from scipy.signal import butter, sosfiltfilt
    sos = butter(4, [f1 / (fs / 2), f2 / (fs / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def load_templates():
    """(fam, sta) -> (3, TL) template + the family hypocentre."""
    m = pd.read_csv(f"{STACKDIR}/metadata.csv")
    m = m[m.station_code.isin(OVERLAP) & m.trace_S_arrival_sample.notna()]
    T, meta = {}, {}
    with h5py.File(f"{STACKDIR}/waveforms.hdf5", "r") as f:
        for r in m.itertuples():
            b, idx = r.trace_name.split("$")
            parts = idx.split(",")
            w = np.asarray(f["data"][b][int(parts[0]), :3, :int(parts[2].lstrip(":"))], float)
            w = bandpass(w, 20.0)
            s0 = int(round(r.trace_S_arrival_sample - PRE * FS))
            if s0 < 0 or s0 + TL > w.shape[1]:
                continue
            seg = w[:, s0:s0 + TL]
            v = seg.ravel()
            nrm = np.linalg.norm(v)
            if nrm <= 0:
                continue
            T[(float(r.source_family_id), r.station_code)] = (v / nrm).astype(np.float32)
            # offset of the template start from the family ORIGIN time
            meta.setdefault(float(r.source_family_id),
                            dict(lat=r.source_latitude_deg, lon=r.source_longitude_deg,
                                 dep=r.source_depth_km))
            meta[float(r.source_family_id)][f"sp_{r.station_code}"] = \
                (r.trace_S_arrival_sample - 800.0) / 20.0
    return T, meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--min-stations", type=int, default=4)
    ap.add_argument("--threshold", default="0.3")
    ap.add_argument("--null-seed", type=int, default=1)
    ap.add_argument("--max-events", type=int, default=4000)
    ap.add_argument("--n-random", type=int, default=4000)
    ap.add_argument("--batch", type=int, default=48)
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    import torch
    dev = args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu"
    print(f"[env] {dev} {torch.cuda.get_device_name(0) if dev == 'cuda' else ''}", flush=True)

    Tpl, meta = load_templates()
    fams = sorted(meta)
    print(f"[tpl] {len(Tpl)} (family, station) templates, {len(fams)} families, "
          f"{len({k[1] for k in Tpl})} stations", flush=True)

    stas = pd.read_csv(f"{RES}/e2_stations.csv").set_index("sta")
    # per (family, station): template start offset from the family origin
    off = np.full((len(fams), len(OVERLAP)), np.nan)
    have = np.zeros((len(fams), len(OVERLAP)), bool)
    for i, fid in enumerate(fams):
        mm = meta[fid]
        for j, s in enumerate(OVERLAP):
            if (fid, s) not in Tpl or f"sp_{s}" not in mm:
                continue
            dx = np.hypot((float(stas.loc[s, "lon"]) - mm["lon"]) * 111.195
                          * np.cos(np.radians(48.6)),
                          (float(stas.loc[s, "lat"]) - mm["lat"]) * 111.195)
            ttP = v1d.tt(dx, mm["dep"], "P")
            off[i, j] = ttP + mm[f"sp_{s}"] - PRE
            have[i, j] = True
    print(f"[tpl] template-start offsets from origin: "
          f"{np.nanmin(off):.1f}..{np.nanmax(off):.1f} s; "
          f"stations per family median {have.sum(1).mean():.1f}", flush=True)

    # stacked template tensor (nfam, nsta, 480), zeros where missing
    TT = np.zeros((len(fams), len(OVERLAP), 3 * TL), np.float32)
    for i, fid in enumerate(fams):
        for j, s in enumerate(OVERLAP):
            if have[i, j]:
                TT[i, j] = Tpl[(fid, s)]
    TTt = torch.from_numpy(TT).to(dev)
    HAVE = torch.from_numpy(have).to(dev)

    # ---------------- event sets ----------------
    thr = args.threshold
    def load_ev(path):
        e = pd.read_csv(path)
        return e[e["n_stations"] >= args.min_stations].time.to_numpy()
    sets = {}
    sets["assoc"] = load_ev(f"{RES}/e2_events_vi2010_cascadia_t{thr}.csv")
    sets["null"] = load_ev(f"{RES}/e3_nullevents_vi2010_cascadia_t{thr}_s{args.null_seed}.csv")
    bo = pd.read_parquet(f"{RES}/e2_gt_bostock.parquet")
    fo = pd.read_csv(f"{RES}/e3_bostock_family_offsets.csv")
    acc = fo[fo.accept.astype(bool)]
    omap = dict(zip(acc.bkey.astype(int), acc.offset_s.astype(float)))
    bo = bo[bo.bkey.astype(int).isin(omap)].copy()
    bo["t"] = (bo.time.astype("datetime64[ns]").astype("int64") / 1e9
               + bo.bkey.astype(int).map(omap).to_numpy())
    t0w, t1w = sets["assoc"].min(), sets["assoc"].max()
    bo = bo[(bo.t >= t0w) & (bo.t <= t1w)]
    rng = np.random.default_rng(11)
    if len(bo) > args.max_events:
        bo = bo.iloc[np.sort(rng.choice(len(bo), args.max_events, replace=False))]
    sets["bostock"] = bo.t.to_numpy()
    sets["random"] = np.sort(rng.uniform(t0w, t1w, args.n_random))
    for k in sets:
        if len(sets[k]) > args.max_events:
            sets[k] = np.sort(rng.choice(sets[k], args.max_events, replace=False))
    print("[sets] " + ", ".join(f"{k}={len(v)}" for k, v in sets.items()), flush=True)
    # the family each Bostock event belongs to, for the own-template score
    bo_fam = np.array([fams.index(f) if f in fams else -1
                       for f in bo.fam_id.to_numpy()], int) if len(bo) else np.zeros(0, int)

    # ---------------- load the 6 stations' 20 Hz data for the window ----------------
    days = pd.date_range(pd.Timestamp(t0w, unit="s").normalize(),
                         pd.Timestamp(t1w, unit="s").normalize(), freq="D")
    from obspy import read
    from scipy.signal import resample_poly
    DATA = {}
    for s in OVERLAP:
        for d in days:
            p = f"{TRACES}/{NETSTA[s]}/{d.year}/{d.dayofyear:03d}.mseed"
            if not os.path.exists(p):
                continue
            st = read(p); st.merge(method=1, fill_value=0)
            fs0 = float(st[0].stats.sampling_rate)
            npts = min(tr.stats.npts for tr in st)
            ch = sorted({tr.stats.channel for tr in st})
            # FIX 2026-09-11 (E4 step 0): old line was [Z]+sorted(non-Z) = Z,E,N for CN stations,
            # but the stack templates are Z,N,E (component_order "Z12", 1=N -- verified by
            # e4_polarization_check.py, 7.9 deg vs 37.3 deg).  Every E3 step-4 cc was 90 deg rotated.
            order = chanorder.resolve_order(ch, station=f"{s}/{d.date()}", legacy=LEGACY_CHAN)
            if s not in _ORDER_SEEN:
                print(f"  [chanorder] {s}: {order}", flush=True); _ORDER_SEEN[s] = order
            X = np.vstack([st.select(channel=c)[0].data[:npts].astype(float) for c in order])
            X = bandpass(X, fs0)
            q = int(round(fs0 / FS))
            X = resample_poly(X, 1, q, axis=-1) if q > 1 else X
            DATA[(s, d)] = (st[0].stats.starttime.timestamp, X.astype(np.float32))
        print(f"  [data] {s}: {sum(1 for k in DATA if k[0]==s)} days", flush=True)

    # ---------------- the network cc, as a grouped conv1d ----------------
    # For each (event, family, station) the data segment starts at the template's nominal position
    # minus 3 s and is 280 samples long; a grouped convolution with the 129 unit-norm templates as
    # filters gives the numerator at all 121 integer lags at once, and the same convolution of the
    # SQUARED data with a ones-kernel gives the window norms.  (The first version looped over
    # families in Python -- 12 M iterations -- and could not finish.)
    NLAG = int(round(2 * 3.0 * FS)) + 1          # 121 lags = +-3.0 s at 0.05 s
    WIN = TL + NLAG - 1
    LAG0 = -3.0
    nf = len(fams)
    Wt = TTt.permute(1, 0, 2).reshape(len(OVERLAP), nf, 3, TL).contiguous()   # (nsta, nfam, 3, TL)
    Wones = torch.ones_like(Wt)
    import torch.nn.functional as Fnn
    rows = []
    # BUGFIX 2026-09-11 (first run): the MAD scale was taken from the RANDOM-time distribution of
    # the BEST-of-129 cc, which is a maximum over 129 families and therefore far above any single
    # template's background.  Scoring Bostock's own-template cc against it gave a nonsensical
    # -5.7 MAD.  The full (n_events, n_fam) matrix is now kept so every cc can be expressed in the
    # MAD of ITS OWN template's random-time distribution, which is what the acceptance asks for.
    PERFAM, TIMES = {}, {}
    for setname, times in sets.items():
        times = np.sort(times)
        own = None
        if setname == "bostock":
            o = np.argsort(sets["bostock"])
            own = bo_fam[o]
        num = np.zeros((len(times), nf), np.float32)
        cnt = np.zeros((len(times), nf), np.float32)
        best_lag = np.zeros((len(times), nf), np.float32)
        acc = torch.zeros((0,))
        dayidx = pd.Series(pd.to_datetime(times, unit="s")).dt.normalize()
        for d, gidx in pd.Series(np.arange(len(times))).groupby(dayidx.values):
            ii = gidx.to_numpy()
            tb = times[ii]
            A = torch.zeros((len(ii), nf, NLAG), device=dev)
            C = torch.zeros((len(ii), nf), device=dev)
            for j, stn in enumerate(OVERLAP):
                key = (stn, pd.Timestamp(d))
                if key not in DATA:
                    continue
                tstart, X = DATA[key]
                i0 = np.rint((tb[:, None] + off[None, :, j] + LAG0 - tstart) * FS).astype(np.int64)
                ok = have[None, :, j] & (i0 >= 0) & (i0 + WIN <= X.shape[1])
                i0c = np.clip(i0, 0, max(X.shape[1] - WIN, 0))
                idx = i0c[:, :, None] + np.arange(WIN)[None, None, :]
                seg = X[:, idx]                                   # (3, nb, nfam, WIN)
                seg = np.transpose(seg, (1, 2, 0, 3))             # (nb, nfam, 3, WIN)
                D = torch.from_numpy(np.ascontiguousarray(seg)).to(dev)
                nb_ = D.shape[0]
                Dr = D.reshape(nb_, nf * 3, WIN)
                n = Fnn.conv1d(Dr, Wt[j].reshape(nf, 3, TL), groups=nf)        # (nb, nf, NLAG)
                e = Fnn.conv1d(Dr * Dr, Wones[j].reshape(nf, 3, TL), groups=nf)
                cc = n / e.clamp_min(1e-20).sqrt()
                M = torch.from_numpy(ok).to(dev).unsqueeze(-1)
                A += torch.where(M, cc, torch.zeros_like(cc))
                C += M.squeeze(-1).float()
            net = A / C.clamp_min(1).unsqueeze(-1)
            net = torch.where((C > 0).unsqueeze(-1), net, torch.full_like(net, -9.0))
            v, l = net.max(dim=2)
            num[ii] = v.cpu().numpy()
            cnt[ii] = C.cpu().numpy()
            best_lag[ii] = (l.cpu().numpy() * (1.0 / FS) + LAG0)
        per_fam = np.where(cnt > 0, num, np.nan)
        with np.errstate(all="ignore"):
            bi = np.nanargmax(np.nan_to_num(per_fam, nan=-9.0), axis=1)
        best_cc = per_fam[np.arange(len(times)), bi]
        srt = np.sort(np.nan_to_num(per_fam, nan=-9.0), axis=1)[:, ::-1]
        own_cc = np.full(len(times), np.nan)
        if own is not None:
            m = own >= 0
            own_cc[m] = per_fam[np.where(m)[0], own[m]]
        PERFAM[setname] = per_fam.astype(np.float32)     # (n_events, n_fam), for the PER-TEMPLATE
        TIMES[setname] = times                            # background (see the note below)
        rows.append(pd.DataFrame(dict(
            set=setname, time=times, best_fam=[fams[k] for k in bi], best_cc=best_cc,
            best_lag=best_lag[np.arange(len(times)), bi], own_cc=own_cc,
            cc2=srt[:, 1] if srt.shape[1] > 1 else np.nan,
            cc8=srt[:, 7] if srt.shape[1] > 7 else np.nan,
            n_sta=np.nanmax(cnt, axis=1))))
        print(f"  [{setname}] {len(times)} events scored", flush=True)
    out = pd.concat(rows, ignore_index=True)
    out = out[np.isfinite(out.best_cc) & (out.best_cc > -5)]
    out.to_parquet(f"{RES}/e3_classifier_events.parquet", index=False)
    np.savez_compressed(f"{RES}/e3_classifier_perfam.npz",
                        fams=np.asarray(fams, float),
                        **{f"cc_{k}": v for k, v in PERFAM.items()},
                        **{f"t_{k}": TIMES[k] for k in PERFAM})

    # ---------------- acceptance ----------------
    rnd = out[out.set == "random"].best_cc.to_numpy()
    p99 = float(np.percentile(out[out.set == "null"].best_cc.to_numpy(), 99))
    med = float(np.median(rnd)); mad = float(np.median(np.abs(rnd - med)))
    # PER-TEMPLATE background from the random times: one median and MAD per family
    Rm = PERFAM["random"]
    fmed = np.nanmedian(Rm, axis=0)
    fmad = np.nanmedian(np.abs(Rm - fmed[None, :]), axis=0)
    print("[bg] per-template random-time cc: median of the family medians "
          f"{np.nanmedian(fmed):.4f}, median MAD {np.nanmedian(fmad):.4f}")
    def in_own_mad(M, fam_index):
        z = (M[np.arange(len(fam_index)), fam_index] - fmed[fam_index]) / \
            np.clip(fmad[fam_index], 1e-9, None)
        return z
    mrows = []
    for k, M in PERFAM.items():
        Z = (M - fmed[None, :]) / np.clip(fmad[None, :], 1e-9, None)
        bz = np.nanmax(Z, axis=1)
        mrows.append(dict(set=k, n=len(M), best_mad_median=float(np.nanmedian(bz)),
                          best_mad_p90=float(np.nanpercentile(bz, 90)),
                          frac_best_mad_ge6=float(np.nanmean(bz >= 6)),
                          frac_best_mad_ge8=float(np.nanmean(bz >= 8))))
    if len(bo):
        m = bo_fam >= 0
        ob = PERFAM["bostock"]
        ordr = np.argsort(sets["bostock"])
        fi = bo_fam[ordr]
        ok = fi >= 0
        z = (ob[np.where(ok)[0], fi[ok]] - fmed[fi[ok]]) / np.clip(fmad[fi[ok]], 1e-9, None)
        mrows.append(dict(set="bostock_own_template", n=int(ok.sum()),
                          best_mad_median=float(np.nanmedian(z)),
                          best_mad_p90=float(np.nanpercentile(z, 90)),
                          frac_best_mad_ge6=float(np.nanmean(z >= 6)),
                          frac_best_mad_ge8=float(np.nanmean(z >= 8))))
    MAD = pd.DataFrame(mrows)
    MAD.to_csv(f"{RES}/e3_classifier_mad.csv", index=False)
    print("\n[per-TEMPLATE MAD, the scale D3's '3.5 MAD' was quoted in]")
    print(MAD.round(3).to_string(index=False))
    srows = []
    for k, g in out.groupby("set"):
        v = g.best_cc.to_numpy()
        srows.append(dict(set=k, n=len(g), cc_median=float(np.median(v)),
                          cc_p90=float(np.percentile(v, 90)),
                          cc_p99=float(np.percentile(v, 99)),
                          frac_above_null_p99=float((v > p99).mean()),
                          mad_median=float((np.median(v) - med) / max(mad, 1e-12)),
                          mad_p90=float((np.percentile(v, 90) - med) / max(mad, 1e-12))))
    ob = out[(out.set == "bostock") & out.own_cc.notna()]
    if len(ob):
        v = ob.own_cc.to_numpy()
        srows.append(dict(set="bostock_own_template", n=len(ob),
                          cc_median=float(np.median(v)), cc_p90=float(np.percentile(v, 90)),
                          cc_p99=float(np.percentile(v, 99)),
                          frac_above_null_p99=float((v > p99).mean()),
                          mad_median=float((np.median(v) - med) / max(mad, 1e-12)),
                          mad_p90=float((np.percentile(v, 90) - med) / max(mad, 1e-12))))
    S = pd.DataFrame(srows)
    S["null_p99"] = p99
    S["random_median"] = med
    S["random_mad"] = mad
    S.to_csv(f"{RES}/e3_classifier_summary.csv", index=False)
    print("\n[classifier] random-time background: median cc {:.4f}, MAD {:.4f}; "
          "null-event best-cc 99th pct {:.4f}".format(med, mad, p99))
    print(S.round(4).to_string(index=False))
    print(f"\n[write] {RES}/e3_classifier_events.parquet, e3_classifier_summary.csv")


if __name__ == "__main__":
    main()
