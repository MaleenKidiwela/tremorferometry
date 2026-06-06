"""Overlay all station-metadata changes (FDSN channel epochs) on a dv/v series,
so dv/v features can be checked against instrument/response changes.

Major changes (sensor type or sample-rate change) are drawn red+solid+labeled;
minor changes (same sensor, new response epoch) are grey+dashed. Reusable across
stations.
"""
from __future__ import annotations

import argparse

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from obspy import UTCDateTime
from obspy.clients.fdsn import Client


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--dvv-csv", required=True)
    p.add_argument("--network", default="UW")
    p.add_argument("--station", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--smooth-days", type=int, default=60)
    p.add_argument("--eq-date", default="2001-02-28", help="reference EQ (Nisqually); '' to skip")
    return p.parse_args()


def metadata_changes(network, station):
    """Return sorted list of (date, label, is_major) for every vertical-channel epoch
    boundary, flagging sensor/sample-rate changes as major."""
    c = Client("EARTHSCOPE")
    inv = c.get_stations(network=network, station=station, level="channel",
                         channel="*HZ", starttime="1985-01-01", endtime="2027-01-01")
    eps = []
    for net in inv:
        for sta in net:
            for ch in sta:
                sens = (ch.sensor.description if ch.sensor else "").split(",")[0]
                eps.append((ch.start_date, ch.code, round(ch.sample_rate), sens))
    eps.sort(key=lambda r: r[0])
    out, prev = [], None
    for s, code, sr, sens in eps:
        if prev is None:
            prev = (code, sr, sens)
            continue
        major = (sens != prev[2]) or (sr != prev[1]) or (code[0] != prev[0][0])
        label = f"{code} {sr}Hz {sens}" if major else f"{code} epoch"
        out.append((pd.Timestamp(s.datetime), label, major))
        prev = (code, sr, sens)
    # also emit the very first epoch as the install
    if eps:
        s0, code0, sr0, sens0 = eps[0]
        out.insert(0, (pd.Timestamp(s0.datetime), f"{code0} {sr0}Hz {sens0} (install)", True))
    return out


def main():
    args = parse_args()
    d = pd.read_csv(args.dvv_csv)
    d["date"] = pd.to_datetime(d["date"])
    med = d.groupby("date")["dvv"].median() * 100
    med = med.rolling(args.smooth_days, center=True, min_periods=args.smooth_days // 3).median()

    changes = metadata_changes(args.network, args.station)
    t0, t1 = med.index.min(), med.index.max()
    changes = [(t, l, m) for (t, l, m) in changes if t0 <= t <= t1]

    fig, ax = plt.subplots(figsize=(13, 5))
    # raw daily median (light) + smoothed (bold)
    raw = d.groupby("date")["dvv"].median() * 100
    ax.plot(raw.index, raw.values, color="0.8", lw=0.4, zorder=1)
    ax.plot(med.index, med.values, color="k", lw=1.6, zorder=4,
            label=f"cross-patch {args.smooth_days}-d median")
    ax.axhline(0, color="0.6", lw=0.6)

    for t, lab, major in changes:
        ax.axvline(t, color=("red" if major else "0.55"),
                   lw=(1.8 if major else 1.0), ls=("-" if major else ":"), zorder=3)
        ax.text(t, ax.get_ylim()[1], " " + lab, rotation=90, va="top", ha="left",
                fontsize=6.5, color=("red" if major else "0.4"), zorder=5)

    if args.eq_date:
        ax.axvline(pd.Timestamp(args.eq_date), color="purple", lw=1.3, ls="--", zorder=3)
        ax.text(pd.Timestamp(args.eq_date), ax.get_ylim()[0], " Nisqually M6.8",
                rotation=90, va="bottom", ha="right", fontsize=7, color="purple")

    ax.set_ylabel("dv/v (%)")
    ax.set_xlabel("date")
    ax.xaxis.set_major_locator(mdates.YearLocator(2))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    ax.set_title(f"{args.station} coda dv/v with station-metadata changes "
                 "(red = sensor/rate change; grey = response epoch)", fontsize=10)
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"wrote {args.out}")
    print("metadata changes overlaid:")
    for t, lab, m in changes:
        print(f"  {t.date()}  {'MAJOR' if m else 'minor'}  {lab}")


if __name__ == "__main__":
    main()
