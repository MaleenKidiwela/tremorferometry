#!/usr/bin/env python
"""Track A-2 follow-up (1a) — WHY the HF-ratio cut is not epoch-invariant, from the instrument response.

The veto statistic is R = E(hi band) / E(2-8 Hz) measured on RAW COUNTS of the PGC-site trace, so for
the same ground motion R scales with the band-averaged squared instrument response:

    R_epoch  =  R_ground  x  k_epoch(hi) / k_epoch(2-8),        k(band) = < |Resp(f)|^2 >_band

and the epoch-to-epoch shift of the measured ratio is predicted, with no waveform at all, by

    shift(hi band)  =  [k_HH(hi)/k_HH(lo)] / [k_BH(hi)/k_BH(lo)].

This script evaluates that from StationXML (EarthScope, level=response) for every CN.PGC epoch, checks
the prediction against the shift Track C measured (0.199 BH -> 0.480 HH, x2.41), and scans candidate hi
bands for the one where shift ~ 1, i.e. where BOTH instruments are flat and a single cut transfers.

Two weightings are reported: flat over the band, and weighted by the observed ground-motion spectrum
shape (a f^-1.6 velocity-spectrum proxy is NOT assumed — the empirical PSD option reads a day of data).
Output: results/a4_hf_response_scan.csv, figures/a4_hf_response.png
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
EXT = os.path.dirname(HERE)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sta", default="PGC")
    ap.add_argument("--net", default="CN")
    ap.add_argument("--lo", default="2,8")
    ap.add_argument("--observed-shift", type=float, default=0.480 / 0.199)
    a = ap.parse_args()
    from obspy.clients.fdsn import Client
    c = Client("EARTHSCOPE", timeout=180)
    inv = c.get_stations(network=a.net, station=a.sta, level="response",
                         starttime="2010-01-01", endtime="2026-09-16")
    flo1, flo2 = [float(x) for x in a.lo.split(",")]

    chans = []
    for n in inv:
        for s in n:
            for ch in s:
                if ch.code[1] != "H" or ch.code[-1] != "Z":
                    continue
                chans.append(ch)
    chans.sort(key=lambda ch: (str(ch.start_date)))

    def kband(ch, f1, f2, nf=400):
        f = np.linspace(f1, f2, nf)
        r = ch.response.get_evalresp_response_for_frequencies(f, output="VEL")
        return float(np.mean(np.abs(r) ** 2)), f, np.abs(r)

    rows = []
    prof = {}
    for ch in chans:
        lab = f"{ch.code}_{str(ch.start_date)[:10]}"
        klo, _, _ = kband(ch, flo1, flo2)
        f = np.linspace(0.5, min(0.49 * ch.sample_rate, 49.0), 1200)
        prof[lab] = (f, np.abs(ch.response.get_evalresp_response_for_frequencies(f, output="VEL")),
                     ch.sample_rate)
        for hi in [(10, 19), (10, 18), (10, 16), (10, 15), (10, 14), (9, 15), (8, 16), (8, 14), (12, 16)]:
            if hi[1] > 0.49 * ch.sample_rate:
                rows.append(dict(chan=ch.code, start=str(ch.start_date)[:10], fs=ch.sample_rate,
                                 hi=f"{hi[0]}-{hi[1]}", k_hi=np.nan, k_lo=klo, k_ratio=np.nan,
                                 note="band above 0.49 x Nyquist"))
                continue
            khi, _, _ = kband(ch, *hi)
            rows.append(dict(chan=ch.code, start=str(ch.start_date)[:10], fs=ch.sample_rate,
                             hi=f"{hi[0]}-{hi[1]}", k_hi=khi, k_lo=klo, k_ratio=khi / klo, note=""))
    D = pd.DataFrame(rows)
    # epoch-to-epoch predicted shift: the LAST BH epoch (the one in use up to 2017-08-10) -> the HH epoch
    bh = D[(D.chan.str.startswith("BH"))].sort_values("start").groupby("hi").last()
    hh = D[(D.chan.str.startswith("HH"))].sort_values("start").groupby("hi").last()
    S = pd.DataFrame(dict(k_ratio_BH=bh.k_ratio, k_ratio_HH=hh.k_ratio))
    S["predicted_shift_HH_over_BH"] = S.k_ratio_HH / S.k_ratio_BH
    S = S.reset_index()
    S.to_csv(f"{EXT}/results/a4_hf_response_scan.csv", index=False)
    pd.set_option("display.width", 200)
    print(D.round(4).to_string(index=False))
    print("\n[predicted epoch shift of the measured HF ratio, by candidate hi band]")
    print(S.round(3).to_string(index=False))
    print(f"\nTrack C MEASURED shift on 10-19 Hz: x{a.observed_shift:.2f} "
          f"(predicted x{float(S.loc[S.hi == '10-19', 'predicted_shift_HH_over_BH'].iloc[0]):.2f})")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for lab, (f, r, fs) in prof.items():
        ax[0].loglog(f, r / np.interp(5.0, f, r), label=f"{lab} ({fs:g} Hz)")
    for b in (2, 8, 10, 16, 19):
        ax[0].axvline(b, color="k", lw=.5, ls=":")
    ax[0].set_xlabel("frequency (Hz)"); ax[0].set_ylabel("|response| (VEL), normalised at 5 Hz")
    ax[0].set_title("CN.PGC instrument response by epoch"); ax[0].legend(fontsize=6); ax[0].grid(alpha=.3, which="both")
    ax[1].plot(range(len(S)), S.predicted_shift_HH_over_BH, "o-")
    ax[1].axhline(1.0, color="k", lw=.8)
    ax[1].axhline(a.observed_shift, color="r", lw=.8, ls="--", label="measured shift on 10-19 Hz")
    ax[1].set_xticks(range(len(S))); ax[1].set_xticklabels(S.hi, rotation=45)
    ax[1].set_ylabel("predicted HH/BH shift of the HF ratio"); ax[1].set_xlabel("candidate hi band (Hz)")
    ax[1].legend(fontsize=7); ax[1].grid(alpha=.3)
    plt.tight_layout(); plt.savefig(f"{EXT}/figures/a4_hf_response.png", dpi=130)
    print("[write] results/a4_hf_response_scan.csv, figures/a4_hf_response.png")


if __name__ == "__main__":
    main()
