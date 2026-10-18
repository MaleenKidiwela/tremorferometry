#!/usr/bin/env python
"""D3 stage 2a — build 3-COMPONENT family templates.

The project's family templates (`data/<sta>_disc.npz`) are Z-ONLY: each is the mean of the family's
member cuts, every cut demeaned + unit-normalised, window [template_pre, template_post) = [-1, +1) s
around the member (candidate) time, 40 Hz, 2-8 Hz bandpass (see scripts/discover_gpu.py `_cut_day`
and `main`). 3-C matched filtering needs the SAME stack on the two horizontals.

This script reproduces that construction on all three components using the member times in
`<sta>_disc.members.parquet`. Member times are spread over the station's whole record and are NOT
inside the D3 pilot download, so short 3-C cutouts (+-PAD s) are fetched from FDSN with
`get_waveforms_bulk` and cached under <outroot>/cutouts/<NET>.<STA>/.

Verification built in: the rebuilt Z template is cross-correlated with the stored Z template; the
distribution of that cc is printed and written to the report CSV (it should be ~1.0).

Usage:
  python lfe_catalog_repair/scripts/d3_build_3c_templates.py --net UW --sta WISH \
     --npz data/wish_disc.npz --members data/wish_disc.members.parquet \
     --summary data/wish_disc.summary.csv --chan 'BH?' \
     --out lfe_catalog_repair/results/d3_templates3c_WISH.npz
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd

PAD = 40.0          # seconds of pad each side of a member time (filter transient + pre-window RMS)
FS = 40.0
FMIN, FMAX = 2.0, 8.0
T_PRE, T_POST = -1.0, 1.0
M = int(round((T_POST - T_PRE) * FS))    # 80


def _client(net):
    from obspy.clients.fdsn import Client
    provs = (["NCEDC", "https://service.ncedc.org"] if net in ("BK", "NC")
             else ["IRIS", "https://service.iris.edu", "EARTHSCOPE"])
    for prov in provs:
        for _ in range(3):
            try:
                return Client(prov, timeout=120)
            except Exception:
                time.sleep(3)
    raise RuntimeError("no FDSN service reachable")


def fetch_cutouts(net, sta, chan, times, cachedir, chunk=150):
    """Fetch +-PAD s 3-C cutouts for `times` (list of pd.Timestamp), cached per chunk as mseed."""
    from obspy import UTCDateTime, read
    cachedir = Path(cachedir)
    cachedir.mkdir(parents=True, exist_ok=True)
    cl = None
    streams = []
    for i0 in range(0, len(times), chunk):
        part = times[i0:i0 + chunk]
        f = cachedir / f"chunk_{i0:06d}.mseed"
        nd = cachedir / f"chunk_{i0:06d}.nd"
        if f.exists():
            try:
                streams.append(read(str(f)))
                continue
            except Exception:
                f.unlink(missing_ok=True)
        if nd.exists():
            continue
        bulk = [(net, sta, "*", chan,
                 UTCDateTime(pd.Timestamp(t).to_pydatetime()) - PAD,
                 UTCDateTime(pd.Timestamp(t).to_pydatetime()) + PAD) for t in part]
        st = None
        for _ in range(3):
            try:
                if cl is None:
                    cl = _client(net)
                st = cl.get_waveforms_bulk(bulk)
                break
            except Exception as e:
                if "No data" in str(e) or "204" in str(e):
                    st = None
                    break
                time.sleep(5)
                cl = None
        if st is None or len(st) == 0:
            nd.touch()
            continue
        st.write(str(f) + ".tmp", format="MSEED")
        os.replace(str(f) + ".tmp", str(f))
        streams.append(st)
        print(f"  cutouts {i0 + len(part)}/{len(times)} ({len(st)} traces)", flush=True)
    return streams


def _prep(tr):
    """resample_poly -> FS, demean, bandpass, in place on a copy. Returns (data, t0_epoch) or None."""
    from scipy.signal import resample_poly
    tr = tr.copy()
    nat = float(tr.stats.sampling_rate)
    if abs(nat - FS) > 1e-6:
        fr = Fraction(int(round(FS)), int(round(nat))).limit_denominator(1000)
        try:
            tr.data = resample_poly(tr.data.astype(np.float64), fr.numerator, fr.denominator).astype(np.float32)
            tr.stats.sampling_rate = FS
        except Exception:
            return None
    tr.detrend("demean")
    try:
        tr.filter("bandpass", freqmin=FMIN, freqmax=FMAX, corners=4, zerophase=True)
    except Exception:
        return None
    return tr.data.astype(np.float32), float(tr.stats.starttime.timestamp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--net", required=True)
    ap.add_argument("--sta", required=True)
    ap.add_argument("--npz", required=True, help="stored Z-only template npz (read-only)")
    ap.add_argument("--members", required=True)
    ap.add_argument("--summary", required=True)
    ap.add_argument("--chan", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cutoutroot", default="lfe_catalog_repair/traces/cutouts")
    ap.add_argument("--tmin", default=None, help="drop member times before this ISO date (sensor-era guard)")
    ap.add_argument("--tmax", default=None, help="drop member times at/after this ISO date")
    ap.add_argument("--tag", default="", help="suffix for the cutout cache dir (per sensor era)")
    a = ap.parse_args()

    dz = np.load(a.npz, allow_pickle=True)
    summ = pd.read_csv(a.summary)
    mem = pd.read_parquet(a.members)
    mem["time"] = pd.to_datetime(mem["time"])
    fams = [f for f in summ.family_id.tolist() if f in dz.files]
    mem = mem[mem.family_id.isin(fams)].drop_duplicates(["family_id", "time"]).reset_index(drop=True)
    n0 = len(mem)
    if a.tmin:
        mem = mem[mem.time >= pd.Timestamp(a.tmin)]
    if a.tmax:
        mem = mem[mem.time < pd.Timestamp(a.tmax)]
    mem = mem.reset_index(drop=True)
    if len(mem) != n0:
        print(f"[3c {a.sta}] era filter {a.tmin}..{a.tmax}: {n0} -> {len(mem)} member times, "
              f"{mem.family_id.nunique()} families keep >=1 member", flush=True)
    print(f"[3c {a.sta}] {len(fams)} families, {len(mem)} distinct member times", flush=True)

    times = sorted(mem.time.unique())
    cachedir = Path(a.cutoutroot) / f"{a.net}.{a.sta}{a.tag}"
    streams = fetch_cutouts(a.net, a.sta, a.chan, list(times), cachedir)
    print(f"[3c {a.sta}] {len(streams)} cutout chunks cached in {cachedir}", flush=True)

    # index preprocessed traces by (component, start, end)
    from obspy import Stream
    allst = Stream()
    for s in streams:
        allst += s
    # NOTE: deliberately NOT calling Stream.merge() -- with 10k-30k independent 80 s cutouts obspy's
    # merge is quadratic and effectively never returns (killed after 20 min on WISH's 8,835 traces).
    # Each cutout is already contiguous; `cut()` below picks the right trace by binary search and
    # tries its neighbours, so a cutout that FDSN split in two simply loses that one member.
    comps = sorted({tr.stats.channel[-1] for tr in allst})
    print(f"[3c {a.sta}] components present: {comps}", flush=True)
    prepped = {c: [] for c in comps}          # (t0, t1, data, t0_epoch)
    for tr in allst:
        c = tr.stats.channel[-1]
        r = _prep(tr)
        if r is None:
            continue
        data, t0 = r
        prepped[c].append((t0, t0 + data.size / FS, data, t0))
    for c in comps:
        prepped[c].sort(key=lambda x: x[0])
    starts = {c: np.array([p[0] for p in prepped[c]]) for c in comps}

    def cut(c, tepoch):
        arr = starts.get(c)
        if arr is None or arr.size == 0:
            return None
        k = int(np.searchsorted(arr, tepoch, side="right")) - 1
        for kk in (k, k - 1, k + 1, k + 2):
            if kk < 0 or kk >= len(prepped[c]):
                continue
            t0, t1, data, _ = prepped[c][kk]
            i_peak = int(round((tepoch - t0) * FS))
            i0 = i_peak + int(round(T_PRE * FS))
            i1 = i0 + M
            if i0 < 0 or i1 > data.size:
                continue
            seg = data[i0:i1].astype(np.float64)
            seg = seg - seg.mean()
            nrm = float(np.linalg.norm(seg))
            if nrm == 0 or not np.isfinite(nrm):
                continue
            return (seg / nrm).astype(np.float32)
        return None

    out = {}
    rows = []
    mem["ep"] = mem.time.astype("datetime64[ns]").astype("int64") / 1e9
    for fam, g in mem.groupby("family_id", sort=False):
        stacks = {}
        for c in comps:
            cs = [cut(c, e) for e in g.ep.values]
            cs = [x for x in cs if x is not None]
            if len(cs) == 0:
                continue
            T = np.mean(np.asarray(cs), axis=0)
            n = float(np.linalg.norm(T))
            if n == 0 or not np.isfinite(n):
                continue
            stacks[c] = ((T / n).astype(np.float32), len(cs))
        if not stacks:
            rows.append(dict(family_id=fam, n_members=len(g), cc_z_vs_stored=np.nan,
                             **{f"n_{c}": 0 for c in comps}))
            continue
        Zs = np.asarray(dz[fam], np.float64)
        cczs = np.nan
        if "Z" in stacks:
            Tz = stacks["Z"][0].astype(np.float64)
            a1 = Zs - Zs.mean()
            a2 = Tz - Tz.mean()
            d = np.linalg.norm(a1) * np.linalg.norm(a2)
            cczs = float(a1 @ a2 / d) if d > 0 else np.nan
        for c, (T, n) in stacks.items():
            out[f"{fam}|{c}"] = T
        rows.append(dict(family_id=fam, n_members=len(g), cc_z_vs_stored=cczs,
                         **{f"n_{c}": stacks.get(c, (None, 0))[1] for c in comps}))

    rep = pd.DataFrame(rows)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    np.savez(a.out, **out)
    rep.to_csv(str(Path(a.out).with_suffix("")) + "_report.csv", index=False)
    ok = rep.cc_z_vs_stored.dropna()
    print(f"[3c {a.sta}] wrote {len(out)} component-templates for "
          f"{rep[[c for c in rep.columns if c.startswith('n_')]].gt(0).any(axis=1).sum()} families -> {a.out}")
    if len(ok):
        print(f"  Z-rebuild vs stored Z: cc median {ok.median():.4f}, "
              f"p05 {ok.quantile(0.05):.4f}, frac>0.99 {(ok > 0.99).mean():.3f}, n={len(ok)}")
    for c in comps:
        k = f"n_{c}"
        if k in rep:
            print(f"  component {c}: templates for {(rep[k] > 0).sum()} families, "
                  f"median members {rep.loc[rep[k] > 0, k].median()}")


if __name__ == "__main__":
    sys.exit(main())
