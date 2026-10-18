#!/usr/bin/env python
"""s11_ball2_meta -- S11 ballistic round 2, §1 metadata check (round2/PLAN.md).
(1) EarthScope MUSTANG daily timing metrics for the 7 core stations -> tables/meta/mustang_<NET.STA>_<metric>.csv
(2) channel / sampling rate per data row from the cut-store index -> tables/meta/row_channel.parquet

    python s11_ball2_meta.py fetch
    python s11_ball2_meta.py channels
"""
from __future__ import annotations
import os, sys, io, time, glob, urllib.request
import numpy as np
import pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import s10_common as S
import s10_multi_common as M

OUT = f"{S.ROOT}/svi_science/results/s11_depth/ballistic/round2/tables/meta"
STA = ["CN.PGC", "PB.B009", "PB.B010", "PB.B011", "PB.B926", "CN.VGZ", "CN.SNB"]
METRICS = ["timing_quality", "clock_locked", "suspect_time_tag", "timing_correction", "num_gaps", "max_gap"]
URL = "https://service.earthscope.org/mustang/measurements/1/query?metric={m}&net={n}&sta={s}&timewindow={a},{b}&format=text"


def fetch():
    os.makedirs(OUT, exist_ok=True)
    for ns in STA:
        n, s = ns.split(".")
        for m in METRICS:
            p = f"{OUT}/mustang_{ns}_{m}.csv"
            if os.path.exists(p):
                continue
            parts = []
            for y0 in range(2003, 2027, 4):
                u = URL.format(m=m, n=n, s=s, a=f"{y0}-01-01T00:00:00", b=f"{min(y0 + 4, 2027)}-01-01T00:00:00")
                for _ in range(3):
                    try:
                        txt = urllib.request.urlopen(u, timeout=120).read().decode()
                        break
                    except Exception as e:
                        txt = ""; time.sleep(5)
                lines = txt.splitlines()
                if len(lines) > 2:
                    parts.append(pd.read_csv(io.StringIO("\n".join(lines[1:]))))
            D = pd.concat(parts) if parts else pd.DataFrame(columns=["value", "target", "start", "end", "lddate"])
            D.to_csv(p, index=False)
            print(ns, m, len(D), flush=True)


def channels():
    os.makedirs(OUT, exist_ok=True)
    X = pd.read_parquet(f"{M.TB}/rows.parquet", columns=["row_id", "det_id", "series", "fam"])
    X = X[(X.series == "data") & (X.det_id >= 0)]
    RS = pd.read_parquet(f"{M.TB}/row_station.parquet", columns=["row_id", "netsta", "date"])
    RS = RS[RS.netsta.isin(STA)].merge(X, on="row_id")
    RS["month"] = RS.date.str[:7]
    out = []
    for mo, g in RS.groupby("month"):
        p = f"{S.CUT}/index_{mo}.parquet"
        if not os.path.exists(p):
            continue
        import pyarrow.parquet as pq
        I = pq.read_table(p, columns=["det_id", "netsta", "cha", "fs", "t0"],
                          filters=[("det_id", "in", list(map(int, g.det_id.unique())))]).to_pandas()
        I = I[I.cha.str[-1] == "Z"]
        out.append(g.merge(I, on=["det_id", "netsta"], how="inner")[["row_id", "netsta", "date", "fam", "cha", "fs", "t0"]])
    R = pd.concat(out)
    R["t0_frac_ms"] = ((R.t0 * R.fs) % 1.0) / R.fs * 1e3
    R.to_parquet(f"{OUT}/row_channel.parquet", index=False)
    print(R.groupby(["netsta", "cha", "fs"]).size())


def compare():
    """per (station, episode): MUSTANG daily metrics (Z channel, mean over the episode's days) and channel band fractions,
    joined to the round-2 EP timing terms (quartet instrument terms and far-station terms)."""
    T = f"{S.ROOT}/svi_science/results/s11_depth/ballistic/round2/tables"
    U = pd.read_csv(f"{M.TB}/episode_units.csv")
    K = pd.read_csv(f"{T}/clock_terms.csv"); K = K[K.bin == "EP"]
    F = pd.read_csv(f"{T}/far_terms.csv"); F = F[F.bin == "EP"]
    terms = pd.concat([K.rename(columns={"lab": "episode"})[["station", "episode", "K", "sigK", "n_fam", "glitch"]]
                       .rename(columns={"K": "term", "sigK": "sig", "glitch": "flag"}),
                       F.rename(columns={"lab": "episode", "far": "station", "F": "term", "sigF": "sig"})
                       [["station", "episode", "term", "sig", "n_fam", "flag"]]])
    terms = terms.groupby(["station", "episode"]).agg(term=("term", "mean"), sig=("sig", "mean"), n_fam=("n_fam", "sum"),
                                                     flag=("flag", "max")).reset_index().merge(U[["episode", "start", "end"]], on="episode")
    met = {}
    for ns in STA:
        for m in METRICS:
            p = f"{OUT}/mustang_{ns}_{m}.csv"
            if not os.path.exists(p):
                continue
            D = pd.read_csv(p)
            if not len(D):
                continue
            D = D[D.target.astype(str).str.contains(r"\.[BHE]HZ\.|\.EHZ\.", regex=True)]
            D["day"] = pd.to_datetime(D.start.str[:10], format="%Y/%m/%d")
            met[(ns, m)] = D.groupby("day").value.mean()
    RC = pd.read_parquet(f"{OUT}/row_channel.parquet")
    RC["day"] = pd.to_datetime(RC.date)
    rows = []
    for r in terms.itertuples():
        a, b = pd.Timestamp(r.start), pd.Timestamp(r.end)
        o = dict(station=r.station, episode=r.episode, term_ms=r.term, sig_ms=r.sig, n_fam=r.n_fam, flag=r.flag)
        for m in METRICS:
            s_ = met.get((r.station, m))
            o[m] = float(s_[(s_.index >= a) & (s_.index <= b)].mean()) if s_ is not None else np.nan
        c = RC[(RC.netsta == r.station) & (RC.day >= a) & (RC.day <= b)]
        o["frac_BH"] = float((c.cha.str[:2] == "BH").mean()) if len(c) else np.nan
        o["fs_dev_ppm_max"] = float((np.abs(c.fs - c.fs.round()) / c.fs.round()).max() * 1e6) if len(c) else np.nan
        rows.append(o)
    R = pd.DataFrame(rows); R.to_csv(f"{OUT}/meta_vs_terms.csv", index=False)
    with pd.option_context("display.width", 250, "display.max_rows", 200):
        print(R[R.flag == True].round(3).to_string())
        print(R.groupby(["station", "flag"])[METRICS + ["frac_BH"]].mean().round(3).to_string())


if __name__ == "__main__":
    {"fetch": fetch, "channels": channels, "compare": compare}[sys.argv[1]]()
