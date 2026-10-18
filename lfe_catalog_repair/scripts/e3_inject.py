#!/usr/bin/env python
"""E3 step 3 — injection recall-vs-SNR and the false-alarm rate on REAL noise.

Bostock family stacks (SeisBench `LFEStacksCascadiaBostock2015`, 20 Hz, 3-C, restituted velocity,
P at sample 800) recorded AT the host station are band-passed 1-8 Hz, scaled to a target
peak-amplitude / host-RMS ratio, resampled to the host sampling rate and added into REAL
quiet-day 3-C day files at known random times.  The picker (both weight sets) is then run on the
injected day AND on the same day un-injected, so the baseline -- chance picks plus any residual
real LFE -- is measured, not assumed, and subtracted.

Host days
  CN.PGC   the 6 quietest 2010-05..07 days (0 Bostock + 0-4 Lin events, 0 PNSN tremor within
           60 km): genuinely quiet, the clean measurement.
  CN.LZB, CN.YOUB  the quietest days available on disk, which all lie INSIDE the 2010-08..09 ETS
           window (only the 34-day vi2010 window was downloaded), so their false-alarm rate F is
           an UPPER bound and their baseline subtraction is doing real work.  Stated, not hidden.

Two SNR measures are recorded for every injection, because the project's existing recall-vs-
magnitude curve is calibrated in the second one:
  `snr_peak` = peak |amplitude| of the injected 1-8 Hz signal / RMS of the 1-8 Hz host in a +-15 s
               window -- the knob the `--levels` are set in;
  `snr_rms`  = RMS of the injected signal over the 6 s from S-1 s / the same host RMS -- the
               definition behind D3's "median Lin LFE at PGC has s/n = 0.078" (which came from
               rms(signal window)/rms(null window) = sqrt(1 + N (s/n)^2) at L = 6 s).

Outputs results/e3_inject_picks.parquet (per injection per model per threshold),
        results/e3_inject_recall.csv, results/e3_inject_falsealarm.csv.
"""
from __future__ import annotations

import argparse
import os
import sys
import time

import h5py
import numpy as np
import pandas as pd

RES = "lfe_catalog_repair/results"
TRACES = "lfe_catalog_repair/traces"
STACKDIR = f"{RES}/external/other/seisbench_lfestacks_cascadia"
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import e2_pick as EP  # noqa: E402
import e4_chanorder as chanorder  # noqa: E402

LEGACY_CHAN = bool(int(os.environ.get("E3_LEGACY_CHAN_ORDER", "0")))  # 1 = reproduce the pre-fix bug
_ORDER_SEEN = {}

LEVELS = (0.1, 0.2, 0.3, 0.5, 1.0, 2.0)
SWEEP = EP.SWEEP
HOSTS = {
    "PGC":  ("CN.PGC",  ["2010-05-20", "2010-05-02", "2010-06-03", "2010-05-30",
                         "2010-06-13", "2010-05-01"]),
    "LZB":  ("CN.LZB",  ["2010-09-12", "2010-09-09", "2010-09-11", "2010-08-12",
                         "2010-08-11", "2010-08-13"]),
    "YOUB": ("CN.YOUB", ["2010-09-12", "2010-09-09", "2010-09-11", "2010-08-12",
                         "2010-08-11", "2010-08-13"]),
}


def load_stacks(station):
    """(n, 3, nt) 20 Hz stacks recorded at `station`, with their S sample."""
    m = pd.read_csv(f"{STACKDIR}/metadata.csv")
    m = m[(m.station_code == station) & (m.trace_S_arrival_sample.notna())]
    out = []
    with h5py.File(f"{STACKDIR}/waveforms.hdf5", "r") as f:
        for r in m.itertuples():
            b, idx = r.trace_name.split("$")
            parts = idx.split(",")
            i = int(parts[0])
            n = int(parts[2].lstrip(":"))
            w = np.asarray(f["data"][b][i, :3, :n], dtype=float)
            out.append((w, float(r.trace_S_arrival_sample), float(r.source_family_id)))
    return out


def bandpass(x, fs, f1=1.0, f2=8.0, corners=4):
    from scipy.signal import butter, sosfiltfilt
    sos = butter(corners, [f1 / (fs / 2), f2 / (fs / 2)], btype="band", output="sos")
    return sosfiltfilt(sos, x, axis=-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stations", default="PGC,LZB,YOUB")
    ap.add_argument("--models", default="cascadia,jcms")
    ap.add_argument("--per-day-per-level", type=int, default=50)
    ap.add_argument("--device", default="cuda")
    ap.add_argument("--torch-threads", type=int, default=4)
    ap.add_argument("--seed", type=int, default=7)
    args = ap.parse_args()

    import torch
    torch.set_num_threads(args.torch_threads)
    from obspy import read
    from scipy.signal import resample_poly
    dev = args.device if (args.device != "cuda" or torch.cuda.is_available()) else "cpu"
    print(f"[env] device={dev} "
          f"{torch.cuda.get_device_name(0) if dev == 'cuda' else ''}", flush=True)
    models = {mn: EP.load_model(mn, dev) for mn in args.models.split(",")}

    rng = np.random.default_rng(args.seed)
    inj_rows, pick_rows, fa_rows = [], [], []
    for stn in args.stations.split(","):
        netsta, daylist = HOSTS[stn]
        stacks = load_stacks(stn)
        if not stacks:
            print(f"[warn] no stacks recorded at {stn}; skipped", flush=True)
            continue
        print(f"[{stn}] {len(stacks)} family stacks, {len(daylist)} host days", flush=True)
        for dstr in daylist:
            d = pd.Timestamp(dstr)
            p = f"{TRACES}/{netsta}/{d.year}/{d.dayofyear:03d}.mseed"
            if not os.path.exists(p):
                print(f"[warn] missing {p}", flush=True)
                continue
            st = read(p)
            st.merge(method=1, fill_value=0)
            fs = float(st[0].stats.sampling_rate)
            t0 = st[0].stats.starttime.timestamp
            npts = min(tr.stats.npts for tr in st)
            # channel order Z, N/1, E/2 to match the stack's Z12
            # FIX 2026-09-11 (E4 step 0): the old three lines below sorted the horizontals
            # alphabetically (BHE before BHN) -> the data were Z,E,N while the stack is Z,N,E
            # (1=N, verified by e4_polarization_check.py: 7.9 deg vs 37.3 deg).  The injected LFE
            # was therefore fed to the picker with its horizontals swapped.
            #   zc = [c for c in chans if c.endswith("Z")]
            #   hc = [c for c in chans if not c.endswith("Z")]
            #   order = zc + hc
            chans = sorted({tr.stats.channel for tr in st})
            order = chanorder.resolve_order(chans, station=f"{stn}/{d.date()}", legacy=LEGACY_CHAN)
            if stn not in _ORDER_SEEN:
                print(f"  [chanorder] {stn}: {order}", flush=True); _ORDER_SEEN[stn] = order
            X = np.vstack([st.select(channel=c)[0].data[:npts].astype(float) for c in order])
            Xf = bandpass(X, fs)

            # ---- injection plan: all levels interleaved on one day
            n_inj = args.per_day_per_level * len(LEVELS)
            slots = np.arange(n_inj) * (npts / fs / n_inj) + 300.0
            slots = slots + rng.uniform(-60, 60, n_inj)
            lev = np.tile(np.array(LEVELS), args.per_day_per_level)
            rng.shuffle(lev)
            Xi = X.copy()
            plan = []
            for ts, L in zip(slots, lev):
                w, s_samp, fid = stacks[rng.integers(len(stacks))]
                wf = bandpass(w, 20.0)
                up = int(round(fs / 20.0))
                wr = resample_poly(wf, up, 1, axis=-1) if up > 1 else wf
                s_r = s_samp * up
                # host RMS in +-15 s of the injected S, over the 3 components
                i_s = int(round(ts * fs))
                a, b = max(i_s - int(15 * fs), 0), min(i_s + int(15 * fs), npts)
                rms_host = float(np.sqrt((Xf[:, a:b] ** 2).mean()))
                if rms_host <= 0 or not np.isfinite(rms_host):
                    continue
                # peak of the stack in the S window (S-2 .. S+8 s)
                sa, sb = int(s_r - 2 * fs), int(s_r + 8 * fs)
                sa, sb = max(sa, 0), min(sb, wr.shape[1])
                pk = float(np.abs(wr[:, sa:sb]).max())
                if pk <= 0:
                    continue
                g = L * rms_host / pk
                wr = wr * g
                i0 = i_s - int(round(s_r))
                j0, j1 = max(i0, 0), min(i0 + wr.shape[1], npts)
                if j1 - j0 < wr.shape[1] * 0.9:
                    continue
                Xi[:, j0:j1] += wr[:, j0 - i0:j1 - i0]
                # the s/n the project's curve is calibrated in: RMS over 6 s from S-1 s
                ra, rb = int(s_r - 1 * fs), int(s_r + 5 * fs)
                snr_rms = float(np.sqrt((wr[:, max(ra, 0):min(rb, wr.shape[1])] ** 2).mean())
                                / rms_host)
                plan.append(dict(station=stn, day=dstr, t_S=t0 + ts, level=float(L),
                                 fam_id=fid, snr_peak=float(L), snr_rms=snr_rms,
                                 rms_host=rms_host))
            plan = pd.DataFrame(plan)
            inj_rows.append(plan)
            print(f"  [{stn} {dstr}] fs={fs:.0f} Hz, {len(plan)} injections", flush=True)

            # ---- pick the un-injected and the injected day with every weight set
            for tagname, DATA in (("base", X), ("inj", Xi)):
                stw = st.copy()
                for c, row in zip(order, DATA):
                    tr = stw.select(channel=c)[0]
                    tr.data = np.ascontiguousarray(row[:npts], dtype=np.float32)
                    tr.stats.npts = npts
                for mn, mdl in models.items():
                    t1 = time.time()
                    pk = EP.picks_for_day(mdl, stw, 1024, 600)
                    df = pd.DataFrame(pk)
                    if len(df) == 0:
                        df = pd.DataFrame(columns=["phase", "time", "prob", "thr_first"])
                    df = df[df.phase == "S"]
                    tt = (pd.to_datetime(df.time).astype("datetime64[ns]").astype("int64")
                          .to_numpy() / 1e9) if len(df) else np.zeros(0)
                    # BUGFIX 2026-09-11: a pick exists at threshold t iff thr_first <= t <= prob
                    # (e2_pick docstring).  The first version filtered on thr_first alone, which
                    # made every threshold identical to the 0.1 pick set.
                    pick_rows.append(dict(station=stn, day=dstr, model=mn, run=tagname,
                                          n_picks=len(df),
                                          times=tt, thr=df.thr_first.to_numpy() if len(df)
                                          else np.zeros(0),
                                          prob=df.prob.to_numpy() if len(df) else np.zeros(0),
                                          seconds=round(time.time() - t1, 1)))
                    for thr in SWEEP:
                        nS = int(((df.thr_first <= thr) & (df.prob >= thr)).sum()) if len(df) else 0
                        fa_rows.append(dict(station=stn, day=dstr, model=mn, run=tagname,
                                            threshold=thr, n_S=nS, day_seconds=npts / fs))
    inj = pd.concat(inj_rows, ignore_index=True)
    inj.to_parquet(f"{RES}/e3_inject_plan.parquet", index=False)
    pd.DataFrame(fa_rows).to_csv(f"{RES}/e3_inject_falsealarm.csv", index=False)

    # ---- recall: matched against the injected S, minus the SAME query on the un-injected day
    P = {(r["station"], r["day"], r["model"], r["run"]): (r["times"], r["thr"], r["prob"])
         for r in pick_rows}
    rows = []
    for (stn, day, mn), _ in pd.DataFrame(pick_rows).groupby(["station", "day", "model"]):
        pass
    keys = sorted({(r["station"], r["day"], r["model"]) for r in pick_rows})
    det = []
    for stn, day, mn in keys:
        g = inj[(inj.station == stn) & (inj.day == day)]
        if len(g) == 0:
            continue
        for run in ("inj", "base"):
            tt, th, pr = P[(stn, day, mn, run)]
            o = np.argsort(tt); tt, th, pr = tt[o], th[o], pr[o]
            for thr in SWEEP:
                s = np.sort(tt[(th <= thr) & (pr >= thr)])
                if len(s) == 0:
                    for tol in (1.0, 3.0):
                        det.append(dict(station=stn, day=day, model=mn, run=run, threshold=thr,
                                        tol=tol, hit=np.zeros(len(g), bool),
                                        dt=np.full(len(g), np.nan), level=g.level.to_numpy(),
                                        snr_rms=g.snr_rms.to_numpy()))
                    continue
                q = g.t_S.to_numpy()
                i = np.searchsorted(s, q)
                dbest = np.full(len(q), np.inf)
                for off in (-1, 0):
                    j = np.clip(i + off, 0, len(s) - 1)
                    dd = s[j] - q
                    m = np.abs(dd) < np.abs(dbest)
                    dbest = np.where(m, dd, dbest)
                for tol in (1.0, 3.0):
                    det.append(dict(station=stn, day=day, model=mn, run=run, threshold=thr,
                                    tol=tol, hit=np.abs(dbest) <= tol, dt=dbest,
                                    level=g.level.to_numpy(), snr_rms=g.snr_rms.to_numpy()))
    D = pd.concat([pd.DataFrame(dict(station=d["station"], day=d["day"], model=d["model"],
                                     run=d["run"], threshold=d["threshold"], tol=d["tol"],
                                     hit=d["hit"], dt=d["dt"], level=d["level"],
                                     snr_rms=d["snr_rms"])) for d in det], ignore_index=True)
    D.to_parquet(f"{RES}/e3_inject_detections.parquet", index=False)
    piv = D.groupby(["station", "model", "threshold", "tol", "level", "run"]).agg(
        n=("hit", "size"), frac=("hit", "mean"),
        snr_rms=("snr_rms", "median")).reset_index()
    w = piv.pivot_table(index=["station", "model", "threshold", "tol", "level", "snr_rms"],
                        columns="run", values=["n", "frac"]).reset_index()
    w.columns = ["_".join([c for c in col if c]).strip("_") for col in w.columns]
    w["recall_raw"] = w["frac_inj"]
    w["baseline"] = w["frac_base"]
    w["recall"] = (w.frac_inj - w.frac_base) / (1.0 - w.frac_base).clip(lower=1e-9)
    bias = D[(D.run == "inj") & (D.tol == 1.0)].groupby(
        ["station", "model", "threshold", "level"]).apply(
        lambda g: pd.Series(dict(bias_s=np.nanmedian(g.dt[g.hit]),
                                 scatter_s=(np.nanpercentile(g.dt[g.hit], 84)
                                            - np.nanpercentile(g.dt[g.hit], 16)) / 2
                                 if g.hit.sum() > 10 else np.nan)),
        include_groups=False).reset_index()
    w = w.merge(bias, on=["station", "model", "threshold", "level"], how="left")
    w.to_csv(f"{RES}/e3_inject_recall.csv", index=False)
    print("\n[recall, tol 1 s, thr 0.3]")
    print(w[(w.tol == 1.0) & (w.threshold == 0.3)][
        ["station", "model", "level", "snr_rms", "recall_raw", "baseline", "recall",
         "bias_s", "scatter_s"]].round(4).to_string(index=False))
    print(f"\n[write] {RES}/e3_inject_recall.csv, e3_inject_falsealarm.csv, "
          f"e3_inject_detections.parquet, e3_inject_plan.parquet")


if __name__ == "__main__":
    main()
