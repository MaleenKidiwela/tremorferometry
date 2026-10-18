#!/usr/bin/env python
"""S11 tomo/eq (PLAN.md §6.3): parse cached ComCat phase-data QuakeML into data/eq_picks/picks.parquet.

  python -I s11_eq_parse.py [--raw DIR] [--out FILE]
One row per (event, pick associated to the preferred origin). Stdlib XML only; files are treated as data.
Also writes data/eq_picks/download_log.md (counts per source / class / region).
"""
import argparse
import glob
import os
import site
import sys
import xml.etree.ElementTree as ET

# run under python -I (no script dir / cwd on sys.path); the user site holds pyarrow, so add it back explicitly
if site.getusersitepackages() not in sys.path:
    sys.path.append(site.getusersitepackages())

import numpy as np
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
SEL = f"{ROOT}/svi_science/results/s11_depth/tomo/eq/selection.csv"
NS = "{http://quakeml.org/xmlns/bed/1.2}"


def tx(el, path, cast=str):
    if el is None:
        return None
    e = el.find(path)
    if e is None or e.text is None:
        return None
    try:
        return cast(e.text.strip())
    except ValueError:
        return None


def parse_file(path, cc):
    root = ET.parse(path).getroot()
    ev = root.find(f".//{NS}event")
    if ev is None:
        return []
    pref = tx(ev, f"{NS}preferredOriginID")
    origins = ev.findall(f"{NS}origin")
    org = next((o for o in origins if o.get("publicID") == pref), None)
    if org is None:          # choose the origin with most arrivals
        org = max(origins, key=lambda o: len(o.findall(f"{NS}arrival")), default=None)
    if org is None:
        return []
    ot = pd.Timestamp(tx(org, f"{NS}time/{NS}value"))
    ot = ot.tz_convert(None) if ot.tzinfo else ot
    olat = tx(org, f"{NS}latitude/{NS}value", float)
    olon = tx(org, f"{NS}longitude/{NS}value", float)
    odep = tx(org, f"{NS}depth/{NS}value", float)
    picks = {}
    for p in ev.findall(f"{NS}pick"):
        w = p.find(f"{NS}waveformID")
        t = pd.Timestamp(tx(p, f"{NS}time/{NS}value"))
        t = t.tz_convert(None) if t.tzinfo else t
        picks[p.get("publicID")] = dict(
            time=t, unc=tx(p, f"{NS}time/{NS}uncertainty", float),
            net=w.get("networkCode") if w is not None else None, sta=w.get("stationCode") if w is not None else None,
            cha=w.get("channelCode") if w is not None else None, loc=w.get("locationCode") if w is not None else None,
            mode=tx(p, f"{NS}evaluationMode"), hint=tx(p, f"{NS}phaseHint"),
            quality=tx(p, f"{NS}onset"))
    rows = []
    for a in org.findall(f"{NS}arrival"):
        pk = picks.get(tx(a, f"{NS}pickID"))
        if pk is None:
            continue
        dist = tx(a, f"{NS}distance", float)
        rows.append(dict(cc_id=cc, origin_time=ot, o_lat=olat, o_lon=olon,
                         o_depth_km=odep / 1000 if odep is not None else np.nan,
                         net=pk["net"], sta=pk["sta"], cha=pk["cha"], loc=pk["loc"],
                         phase=tx(a, f"{NS}phase"), time=pk["time"], unc=pk["unc"], mode=pk["mode"],
                         onset=pk["quality"], time_weight=tx(a, f"{NS}timeWeight", float),
                         time_residual=tx(a, f"{NS}timeResidual", float),
                         dist_deg=dist, dist_km=dist * 111.19 if dist is not None else np.nan,
                         azimuth=tx(a, f"{NS}azimuth", float)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", default=f"{ROOT}/data/eq_picks/raw/comcat")
    ap.add_argument("--out", default=f"{ROOT}/data/eq_picks/picks.parquet")
    a = ap.parse_args()
    files = sorted(glob.glob(f"{a.raw}/*.phase.xml"))
    rows, bad = [], []
    for f in files:
        cc = os.path.basename(f)[: -len(".phase.xml")]
        try:
            rows += parse_file(f, cc)
        except Exception as e:
            bad.append((cc, str(e)[:80]))
    P = pd.DataFrame(rows)
    S = pd.read_csv(SEL)
    S = S[S.source == "comcat"].drop_duplicates("cc_id")
    if len(P):
        P["origin_time"] = P.origin_time.astype("datetime64[ns]")
        P["time"] = P.time.astype("datetime64[ns]")
        P["t_rel"] = (P.time - P.origin_time).dt.total_seconds()
        P["phase1"] = P.phase.str.upper().str[0]
        P["stakey"] = P.net.fillna("") + "." + P.sta.fillna("")
        P = P.merge(S[["cc_id", "id", "tier", "kind", "reg", "zone", "along", "lat", "lon", "depth", "mag", "t_PT"]]
                    .rename(columns=dict(id="event", lat="cat_lat", lon="cat_lon", depth="cat_depth_km")),
                    on="cc_id", how="left")
        P.to_parquet(a.out, index=False)
    # download log
    E = os.path.dirname(a.out)
    st = {k: set(os.path.basename(x).split(".")[0] for x in glob.glob(f"{a.raw}/*.{k}"))
          for k in ("done", "nophase", "fail")}
    S["status"] = np.select([S.cc_id.isin(st["done"]), S.cc_id.isin(st["nophase"]), S.cc_id.isin(st["fail"])],
                            ["picks", "no phase-data", "failed"], "not yet")
    Sall = pd.read_csv(SEL)
    lines = ["# Earthquake-pick download log (S11 tomo/eq, PLAN.md §6)", "",
             f"Written by s11_eq_parse.py at {pd.Timestamp.now(tz='America/Los_Angeles'):%Y-%m-%d %H:%M} PT.", "",
             "## Per source", "",
             f"- USGS ComCat phase-data: {len(S)} events requested; " +
             ", ".join(f"{k} {v}" for k, v in S.status.value_counts().items()),
             f"- NRCan FDSN (includearrivals=true): returns no arrivals (verified 2026-10-13); "
             f"{int((Sall.source == 'none_public').sum())} NRCan-only selected events have no public pick source.",
             "- PNSN: the ComCat `uw` phase-data product is the PNSN pick set (no separate fetch).", "",
             f"- Parsed rows: {len(P)} (P {int((P.phase1 == 'P').sum()) if len(P) else 0}, "
             f"S {int((P.phase1 == 'S').sum()) if len(P) else 0}), events with rows "
             f"{P.cc_id.nunique() if len(P) else 0}, stations {P.stakey.nunique() if len(P) else 0}; "
             f"unparseable files {len(bad)}", "", "## Status by region / class", "",
             "```", S.groupby(["reg", "kind", "status"]).size().unstack(fill_value=0).to_string(), "```"]
    if len(P):
        ps = P.groupby(["cc_id", "stakey"]).phase1.agg(lambda v: {"P", "S"} <= set(v))
        n_ps = ps.groupby(level=0).sum()
        k = P.drop_duplicates("cc_id").set_index("cc_id").reindex(n_ps.index)
        lines += ["", "## Stations with both P and S per event (median / events with >= 4)", "", "```",
                  pd.DataFrame(dict(reg=k.reg, kind=k.kind, n=n_ps)).groupby(["reg", "kind"]).n
                  .agg(events="size", median="median", ge4=lambda v: int((v >= 4).sum())).to_string(), "```"]
    open(f"{E}/download_log.md", "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))
    if bad:
        print("unparseable:", bad[:10])


if __name__ == "__main__":
    main()
