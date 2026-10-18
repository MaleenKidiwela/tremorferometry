#!/usr/bin/env python
"""STEP 1c' STAGE 1 -- D8: is the 2005 bed limited by the TIME BASE, the CHAIN, or the PHYSICS?

Merlin's prescription step 1 (results/merlin_review_step1_2026-09-20.md), over the 50 busiest KEPT
family-hours:

  (a) ENVELOPE STACK per site: stack every catalogued event of the family in that hour at
      t_origin + fam_tt(family, station), 3-C, per-trace rms normalised; report the envelope peak's
      OFFSET from the prediction and its height in background MAD.
  (b) PER-EVENT S-ENERGY EXCESS per site: the 3-s window at the predicted S over the MEDIAN 3-s
      energy of its own +- 60 s neighbourhood, with a +30 s shifted CONTROL.
  (c) CHAIN COMPARISON on >= 5 of those hours: the frozen 20 Hz 1-8 Hz cache against the raw
      100 Hz mseed on disk band-passed 2-8 Hz, compared by linear-stack coherence
      (N x <stack^2> / <trace^2>; 1.0 = incoherent).  NOTHING is downloaded.

VERDICT (pre-registered):
  ACCEPT            median envelope-peak offset within +- 3 s at >= 4 sites AND control excess in
                    [0.9, 1.1] AND chain variants within x1.2 in stack coherence.
  ABANDON_CHAIN     chain variants differ by > x1.5  -> "re-cache at 2-8 Hz first".
  ABANDON_TIMEBASE  envelope peaks absent or displaced > 5 s -> "time base wrong: fix fam_tt /
                    offsets for 2005".
The hidden families are not used anywhere in this stage.
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
import w_common as W                                                        # noqa: E402
import z_common as Z                                                        # noqa: E402

FS = W.FS
CUT_PRE_S = 65.0            # envelope-stack cut starts 65 s before the predicted S ...
CUT_LEN_S = 130.0           # ... and is 130 s long, so the +- 60 s energy neighbourhood and the
                            # +30 s control both fit inside it and the MAD background is the
                            # hour's own tremor
PEAK_SEARCH_S = 10.0        # envelope peak searched within +- 10 s of the prediction
EX_WIN_S = 3.0              # S-energy window
EX_NEIGH_S = 60.0           # +- 60 s neighbourhood
EX_CONTROL_S = 30.0         # +30 s control
CHAIN_WIN_S = 10.0          # chain comparison: the engine's own 10-s window at the predicted S
ACC_PEAK_S = 3.0
ACC_PEAK_SITES = 4
ABANDON_PEAK_S = 5.0
ACC_CONTROL = (0.90, 1.10)
ACC_CHAIN = 1.20
ABANDON_CHAIN = 1.50


def envelope_and_excess(cut, det, mem, order, FT, hours, out):
    rows_env, rows_ex = [], []
    for (bk, d, h, nev) in hours:
        ev = hours_events[(bk, d, h)]
        for s in det:
            tt = FT.get((bk, s))
            if tt is None:
                continue
            segs, ex, exc = [], [], []
            for t0 in ev:
                tS = t0 + tt
                seg = cut.cut(s, tS - CUT_PRE_S, CUT_LEN_S)
                if seg is None:
                    continue
                segs.append(seg)
                # per-event S-energy excess, from the SAME cut
                e = (seg ** 2).sum(0)
                nw = int(round(EX_WIN_S * FS))
                cs = np.concatenate([[0.0], np.cumsum(e)])
                win = cs[nw:] - cs[:-nw]                       # energy of every 3-s window
                iS = int(round(CUT_PRE_S * FS))
                iC = int(round((CUT_PRE_S + EX_CONTROL_S) * FS))
                lo = max(iS - int(EX_NEIGH_S * FS), 0)
                hi = min(iS + int(EX_NEIGH_S * FS), len(win))
                med = float(np.median(win[lo:hi])) if hi > lo else np.nan
                if np.isfinite(med) and med > 0:
                    if iS < len(win):
                        ex.append(float(win[iS] / med))
                    if iC < len(win):
                        exc.append(float(win[iC] / med))
            if len(segs) < 5:
                continue
            Sg = Z.norm_rows(np.stack(segs))
            if len(Sg) < 5:
                continue
            st = Sg.mean(0)
            env = Z.envelope_1d(st)
            nom = int(round(CUT_PRE_S * FS))
            a = max(nom - int(PEAK_SEARCH_S * FS), 0)
            b = min(nom + int(PEAK_SEARCH_S * FS), len(env))
            i = int(a + np.argmax(env[a:b]))
            rows_env.append(dict(bkey=bk, date=d, hour=h, sta=s, site=W.sites()[1].get(s, s),
                                 n_stacked=int(len(Sg)), peak_offset_s=float((i - nom) / FS),
                                 peak_mad=Z.mad_height(env, i),
                                 coherence=Z.stack_coherence(np.stack(segs))))
            rows_ex.append(dict(bkey=bk, date=d, hour=h, sta=s, site=W.sites()[1].get(s, s),
                                n_events=len(ex),
                                excess_med=float(np.median(ex)) if ex else np.nan,
                                excess_p90=float(np.percentile(ex, 90)) if ex else np.nan,
                                control_med=float(np.median(exc)) if exc else np.nan))
        W.log(f"[d8] family {bk} {d} {h:02d}h ({nev} events) done")
    E = pd.DataFrame(rows_env)
    X = pd.DataFrame(rows_ex)
    E.to_csv(out + "_envelope.csv", index=False)
    X.to_csv(out + "_excess.csv", index=False)
    return E, X


def chain_compare(cut, det, FT, hours, out, n_hours=5, max_events=60):
    """Frozen 20 Hz 1-8 Hz cache vs raw 100 Hz 2-8 Hz (and 1-8 Hz), same windows, same statistic."""
    import obspy
    from scipy.signal import butter, sosfiltfilt
    rows = []
    done = 0
    for (bk, d, h, nev) in hours:
        ev = hours_events[(bk, d, h)]
        got = False
        for s in det:
            p = Z.raw_path(s, d)
            tt = FT.get((bk, s))
            if p is None or tt is None:
                continue
            try:
                stz = obspy.read(p)
            except Exception as e:                                          # noqa: BLE001
                W.log(f"[d8/chain] {s} {d}: unreadable ({e})")
                continue
            stz.merge(method=1, fill_value=0)
            chan = {tr.stats.channel[-1]: tr for tr in stz}
            if not all(c in chan for c in "ZNE"):
                continue
            fs_n = float(chan["Z"].stats.sampling_rate)
            t_file = float(chan["Z"].stats.starttime.timestamp)
            x = np.stack([chan[c].data.astype(np.float64) for c in "ZNE"], 0)
            x = x - x.mean(1, keepdims=True)

            def band(lo, hi):
                sos = butter(4, [lo / (fs_n / 2), hi / (fs_n / 2)], btype="band", output="sos")
                return np.ascontiguousarray(sosfiltfilt(sos, x, axis=-1))

            # the variant TAG must not carry the native rate: the bed mixes 100 Hz (C8/PO) and
            # 40 Hz (CN) raw channels, and a rate-tagged column makes the per-station comparison
            # pivot empty.  The rate is recorded in its own column instead.
            variants = {"native_2-8Hz": (band(2.0, 8.0), fs_n, t_file),
                        "native_1-8Hz": (band(1.0, 8.0), fs_n, t_file)}
            # frozen cache windows
            segs = []
            for t0 in ev[:max_events]:
                seg = cut.cut(s, t0 + tt - 1.0, CHAIN_WIN_S)
                if seg is not None:
                    segs.append(seg)
            if len(segs) < 5:
                continue
            rows.append(dict(bkey=bk, date=d, hour=h, sta=s, variant="frozen_20Hz_1-8Hz",
                             fs=FS, n=len(segs), coherence=Z.stack_coherence(np.stack(segs))))
            for tag, (arr, fsx, tf) in variants.items():
                nw = int(round(CHAIN_WIN_S * fsx))
                sg = []
                for t0 in ev[:max_events]:
                    i0 = int(round((t0 + tt - 1.0 - tf) * fsx))
                    if i0 < 0 or i0 + nw > arr.shape[1]:
                        continue
                    sg.append(arr[:, i0:i0 + nw])
                if len(sg) >= 5:
                    rows.append(dict(bkey=bk, date=d, hour=h, sta=s, variant=tag, fs=fsx,
                                     n=len(sg), coherence=Z.stack_coherence(np.stack(sg))))
            got = True
            W.log(f"[d8/chain] {bk} {d} {h:02d}h {s}: "
                  f"{ {r['variant']: round(r['coherence'], 3) for r in rows[-3:]} }")
        if got:
            done += 1
        if done >= n_hours:
            break
    C = pd.DataFrame(rows)
    C.to_csv(out + "_chain.csv", index=False)
    return C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=f"{Z.OUT}/d8")
    ap.add_argument("--n-hours", type=int, default=50)
    ap.add_argument("--chain-hours", type=int, default=5)
    ap.add_argument("--min-events", type=int, default=10)
    a = ap.parse_args()
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    mem, of, order, det = W.sites()
    ST = W.station_table()
    FT = Z.fam_tt_map(det, ST)
    days = {d.strftime("%Y-%m-%d") for d in W.bed_days()}
    B = Z.bostock_year(kept_only=True, window_days=days)
    hours = Z.busiest_family_hours(B, a.n_hours, a.min_events)
    global hours_events
    hours_events = {(bk, d, h): np.sort(B[(B.bkey == bk) & (B.date == d) & (B.hh == h)]
                                        .t_origin.to_numpy(float)) for (bk, d, h, _) in hours}
    W.log(f"[d8] {len(hours)} busiest kept family-hours, "
          f"{sum(len(v) for v in hours_events.values())} events, {len(det)} stations")
    cut = Z.Cutter(maxn=40)
    E, X = envelope_and_excess(cut, det, mem, order, FT, hours, a.out)
    C = chain_compare(cut, det, FT, hours, a.out, n_hours=a.chain_hours)

    # ---------------------------------------------------------------- verdict
    site_off = (E.groupby("site").peak_offset_s.median() if len(E) else pd.Series(dtype=float))
    site_mad = (E.groupby("site").peak_mad.median() if len(E) else pd.Series(dtype=float))
    n_ok = int((site_off.abs() <= ACC_PEAK_S).sum())
    worst_site = float(site_off.abs().max()) if len(site_off) else np.nan
    site_ex = (X.groupby("site").excess_med.median() if len(X) else pd.Series(dtype=float))
    site_p90 = (X.groupby("site").excess_p90.median() if len(X) else pd.Series(dtype=float))
    ctrl = float(X.control_med.median()) if len(X) else np.nan
    ratio = np.nan
    if len(C):
        piv = C.pivot_table(index=["bkey", "date", "hour", "sta"], columns="variant",
                            values="coherence")
        piv = piv.dropna()
        if len(piv):
            r = piv.max(axis=1) / piv.min(axis=1).clip(lower=1e-9)
            ratio = float(r.median())
            ratio_max = float(r.max())
        else:
            ratio_max = np.nan
    else:
        ratio_max = np.nan

    verdict = "ACCEPT"
    reason = ""
    if np.isfinite(ratio) and ratio > ABANDON_CHAIN:
        verdict, reason = "ABANDON_CHAIN", ("chain variants differ by x%.2f (> %.2f): re-cache at "
                                            "2-8 Hz first" % (ratio, ABANDON_CHAIN))
    elif (not len(E)) or (np.isfinite(worst_site) and n_ok < ACC_PEAK_SITES
                          and worst_site > ABANDON_PEAK_S):
        verdict, reason = "ABANDON_TIMEBASE", (
            "envelope peaks absent or displaced > %.1f s: time base wrong, fix s_common.fam_tt / "
            "offsets for 2005 before scoring anything" % ABANDON_PEAK_S)
    elif n_ok < ACC_PEAK_SITES:
        verdict, reason = "BELOW_ACCEPT", (f"only {n_ok} sites with |median peak offset| <= "
                                           f"{ACC_PEAK_S} s (need {ACC_PEAK_SITES})")
    elif not (ACC_CONTROL[0] <= ctrl <= ACC_CONTROL[1]):
        verdict, reason = "BELOW_ACCEPT", f"+30 s control excess {ctrl:.3f} outside {ACC_CONTROL}"
    elif np.isfinite(ratio) and ratio > ACC_CHAIN:
        verdict, reason = "BELOW_ACCEPT", f"chain variants differ by x{ratio:.2f} (> {ACC_CHAIN})"

    out = dict(
        n_family_hours=len(hours), n_events=int(sum(len(v) for v in hours_events.values())),
        n_sites_with_peak=int(len(site_off)),
        median_peak_offset_s_by_site={k: float(v) for k, v in site_off.items()},
        median_peak_mad_by_site={k: float(v) for k, v in site_mad.items()},
        n_sites_within_3s=n_ok, worst_site_abs_offset_s=worst_site,
        median_excess_by_site={k: float(v) for k, v in site_ex.items()},
        p90_excess_by_site={k: float(v) for k, v in site_p90.items()},
        control_excess_median=ctrl,
        chain_variants=sorted(C.variant.unique().tolist()) if len(C) else [],
        chain_ratio_median=ratio, chain_ratio_max=ratio_max,
        n_chain_rows=int(len(C)),
        accept=dict(peak_within_s=ACC_PEAK_S, peak_sites_ge=ACC_PEAK_SITES,
                    control_range=list(ACC_CONTROL), chain_ratio_le=ACC_CHAIN,
                    abandon_chain_gt=ABANDON_CHAIN, abandon_peak_gt_s=ABANDON_PEAK_S),
        verdict=verdict, reason=reason)
    json.dump(out, open(a.out + ".json", "w"), indent=1, default=str)
    print(json.dumps(out, indent=1, default=str))


hours_events = {}

if __name__ == "__main__":
    main()
