#!/usr/bin/env python
"""E4 step 1 (Merlin-amended) — turn the template test into a calibrated SENSITIVITY MODEL.

Inputs  results/e4_tt_events_<tag>.parquet, e4_tt_station_<tag>.csv, e4_tt_netfam_<tag>.csv
Outputs results/e4_step1_registration.csv   per (station, L, magnitude bin): measured excess,
                                            predicted excess, ratio, distance
        results/e4_step1_magslope.csv        MAD-vs-magnitude slope with family fixed effects
        results/e4_step1_channels.csv        single-channel vs 3-C vs network gain
        results/e4_step1_network.csv         network by magnitude bin and template length
        results/e4_step1_model.md            the numbers in prose

PREDICTION (Merlin vet F1).  Per (family, station) the stack metadata gives the implied
single-event amplitude SNR  r1 = 10^(trace_S_snr_db/20) / sqrt(trace_stack_count).  A matched
filter against a template of quality q returns, at an event of amplitude ratio r, a correlation
excess of about q*r, and the cc noise scale is sigma = 1/sqrt(2 B L) / sqrt(n_chan) (B = 7 Hz).
Expressed in the units this test uses (the MAD of the template's own continuous cc series),
    predicted excess [MAD] = q * r1 * A(M) / MAD_scale_measured
with q = 0.89 for the SeisBench stacks and A(M) the amplitude of an event of magnitude M relative
to the family's median magnitude.  A(M) is MEASURED here (the MAD-vs-magnitude slope), not assumed.
"""
from __future__ import annotations
import numpy as np, pandas as pd, sys, os

RES = "/home/jovyan/tremorferometry/lfe_catalog_repair/results"
SD = f"{RES}/external/other/seisbench_lfestacks_cascadia"
RING = ["PGC", "LZB", "VGZ", "SNB", "YOUB", "PFB"]
MB = [(-9, 1.5, "M<1.5"), (1.5, 1.8, "M1.5-1.8"), (1.8, 2.1, "M1.8-2.1"), (2.1, 9, "M>=2.1")]
B_HZ = 7.0
Q_STACK = 0.89
R_KM_DEG = 111.195


def load(tag):
    e = pd.read_parquet(f"{RES}/e4_tt_events_{tag}.parquet")
    st = pd.read_csv(f"{RES}/e4_tt_station_{tag}.csv")
    nf = pd.read_csv(f"{RES}/e4_tt_netfam_{tag}.csv")
    return e, st, nf


def r1_table():
    m = pd.read_csv(f"{SD}/metadata.csv")
    m = m[m.station_code.isin(RING) & m.trace_S_arrival_sample.notna()].copy()
    m["r1"] = 10 ** (m.trace_S_snr_db / 20.0) / np.sqrt(m.trace_stack_count)
    stas = pd.read_csv(f"{RES}/e2_stations.csv").set_index("sta")
    cs = np.cos(np.radians(48.6))
    m["epi_km"] = np.hypot(
        (m.station_code.map(stas.lon) - m.source_longitude_deg) * R_KM_DEG * cs,
        (m.station_code.map(stas.lat) - m.source_latitude_deg) * R_KM_DEG)
    m["hypo_km"] = np.hypot(m.epi_km, m.source_depth_km)
    return m[["source_family_id", "station_code", "r1", "trace_S_snr_db", "trace_stack_count",
              "epi_km", "hypo_km"]].rename(columns={"source_family_id": "fam",
                                                    "station_code": "sta"})


def main(tag="ring_p2"):
    e, st, nf = load(tag)
    R1 = r1_table()
    # per (L, station) measured MAD scale and random-time background
    bg = st.groupby(["L", "sta"]).agg(mad_scale=("mad_scale", "median"),
                                      rand=("rand3_median", "median"),
                                      randZ=("randZ_median", "median"),
                                      fa6=("fa6_per_day", "median"),
                                      fa8=("fa8_per_day", "median")).reset_index()
    nbg = nf.groupby("L").agg(net_rand=("net_rand_med", "median"),
                              net_fa6=("net_fa6", "median"), net_fa8=("net_fa8", "median"),
                              shift_fa6=("shift_fa6", "median"), shift_fa8=("shift_fa8", "median"),
                              n_ch=("n_ch", "median")).reset_index()

    # ---------------- long form: one row per (event, L, station)
    rows = []
    for s in RING:
        c = ["doy", "L", "t_origin", "fam", "mag", f"sta_{s}", f"staZ_{s}", f"lag_{s}"]
        g = e[c].rename(columns={f"sta_{s}": "mad3", f"staZ_{s}": "madZ", f"lag_{s}": "lag"})
        g = g[g.mad3.notna()]
        g["sta"] = s
        rows.append(g)
    LN = pd.concat(rows, ignore_index=True)
    LN = LN.merge(bg[["L", "sta", "rand", "randZ", "mad_scale"]], on=["L", "sta"], how="left")
    LN = LN.merge(R1, on=["fam", "sta"], how="left")
    LN["exc3"] = LN.mad3 - LN["rand"]

    # ---------------- MAD-vs-magnitude slope with FAMILY fixed effects
    ms = []
    for (L, s), g in LN.groupby(["L", "sta"]):
        g = g.dropna(subset=["mag", "exc3"])
        if len(g) < 300:
            continue
        fm = g.groupby("fam")[["mag", "exc3"]].transform("mean")
        dx = (g.mag - fm.mag).to_numpy(); dy = (g.exc3 - fm.exc3).to_numpy()
        k = np.isfinite(dx) & np.isfinite(dy)
        sl = np.polyfit(dx[k], dy[k], 1)[0] if k.sum() > 50 else np.nan
        # between-family slope (confounded with distance/site, reported for contrast)
        fg = g.groupby("fam").agg(m=("mag", "mean"), x=("exc3", "median"))
        slb = np.polyfit(fg.m, fg.x, 1)[0] if len(fg) > 5 else np.nan
        # implied A(M): what the excess WOULD be if amplitude scaled as 10^(1.5 dM)
        base = g.exc3.median()
        ms.append(dict(L=L, sta=s, n=len(g), within_family_slope=sl, between_family_slope=slb,
                       median_excess=base,
                       predicted_slope_if_moment=base * np.log(10) * 1.5))
    MS = pd.DataFrame(ms)
    MS.to_csv(f"{RES}/e4_step1_magslope.csv", index=False)

    # ---------------- registration table
    reg = []
    for (L, s), g in LN.groupby(["L", "sta"]):
        mad_scale = float(bg[(bg.L == L) & (bg.sta == s)].mad_scale.iloc[0])
        sigma_pred = 1.0 / np.sqrt(2 * B_HZ * L) / np.sqrt(3)
        mad_pred = 0.6745 * sigma_pred
        r1m = float(g.r1.median())
        d = float(g.epi_km.median())
        for lo, hi, lab in MB + [(-9, 9, "all")]:
            h = g[(g.mag > lo) & (g.mag <= hi)]
            if len(h) < 50:
                continue
            meas = float(h.exc3.median())
            pred = Q_STACK * r1m / mad_scale
            reg.append(dict(L=L, sta=s, magbin=lab, n=len(h),
                            epi_km=d, r1_median=r1m,
                            mad_scale_meas=mad_scale, mad_scale_pred=mad_pred,
                            excess_meas=meas, excess_pred=pred,
                            ratio=meas / pred if pred > 0 else np.nan,
                            excZ=float((h.madZ - h["randZ"]).median()),
                            lag_median=float(h.lag.median())))
    REG = pd.DataFrame(reg)
    REG.to_csv(f"{RES}/e4_step1_registration.csv", index=False)

    # ---------------- channel scaling: 1 channel (Z) vs 3-C vs network
    ch = []
    for L in sorted(e.L.unique()):
        nb = float(nbg[nbg.L == L].net_rand.iloc[0])
        g = e[e.L == L]
        per = {}
        for s in RING:
            gg = LN[(LN.L == L) & (LN.sta == s)]
            if len(gg) < 200:
                continue
            # Z-only background: the random-time median is not stored per component, so use the
            # station's 3-C random median as the reference for both (they share the same day and
            # the same scan; the Z MAD scale is its own).  Reported as a RATIO, which is what the
            # sqrt(channels) check needs.
            # FIX: the single-channel excess must be referenced to the Z series' OWN random-time
            # background (the max-over-scan background of a 1-channel cc series is not that of the
            # 3-C mean), which the step-1 run now records as randZ_median.
            per[s] = dict(e1=float(gg.madZ.median() - gg["randZ"].median()),
                          e3=float(gg.exc3.median()))
        if not per:
            continue
        e1 = np.median([v["e1"] for v in per.values()])
        e3 = np.median([v["e3"] for v in per.values()])
        enet = float(g.net_mad.median() - nb)
        # ideal network from the per-station values: sqrt(sum of squares) weighted by channels
        ideal = float(np.sqrt(np.sum([v["e3"] ** 2 for v in per.values()])))
        # the network actually averages ~2.8 stations per family, not 6 -> also compute the
        # per-family ideal from the stations that family actually has
        ch.append(dict(L=L, e_1chan=e1, e_3chan=e3, gain_3over1=e3 / e1 if e1 else np.nan,
                       e_network=enet, ideal_6sta=ideal,
                       ratio_net_ideal6=enet / ideal if ideal else np.nan,
                       net_rand_bg=nb))
    CH = pd.DataFrame(ch)
    CH.to_csv(f"{RES}/e4_step1_channels.csv", index=False)

    # ---------------- network by magnitude and template length
    net = []
    for L in sorted(e.L.unique()):
        nb = float(nbg[nbg.L == L].net_rand.iloc[0])
        row0 = nbg[nbg.L == L].iloc[0]
        g = e[e.L == L]
        for lo, hi, lab in MB + [(-9, 9, "all")]:
            h = g[(g.mag > lo) & (g.mag <= hi)]
            if len(h) < 50:
                continue
            net.append(dict(L=L, magbin=lab, n=len(h), n_ch_median=float(row0.n_ch),
                            net_mad_median=float(h.net_mad.median()),
                            net_mad_p90=float(h.net_mad.quantile(.9)),
                            net_rand_bg=nb, excess=float(h.net_mad.median() - nb),
                            frac_ge6=float((h.net_mad >= 6).mean()),
                            frac_ge8=float((h.net_mad >= 8).mean()),
                            shift_frac_ge6=float((h.net_mad_shift >= 6).mean()),
                            shift_frac_ge8=float((h.net_mad_shift >= 8).mean()),
                            fa6_per_tpl_day=float(row0.net_fa6),
                            fa8_per_tpl_day=float(row0.net_fa8),
                            shift_fa8_per_tpl_day=float(row0.shift_fa8)))
    NT = pd.DataFrame(net)
    NT.to_csv(f"{RES}/e4_step1_network.csv", index=False)

    pd.set_option("display.width", 200)
    print("\n=== REGISTRATION (excess in MAD over the random-time background) ===")
    print(REG[REG.magbin == "M1.8-2.1"].round(3).to_string(index=False))
    print("\n=== MAD-vs-magnitude slope (within family) ===")
    print(MS.round(3).to_string(index=False))
    print("\n=== CHANNEL SCALING ===")
    print(CH.round(3).to_string(index=False))
    print("\n=== NETWORK ===")
    print(NT.round(4).to_string(index=False))
    print(f"\n[write] {RES}/e4_step1_{{registration,magslope,channels,network}}.csv")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "ring_p2")
