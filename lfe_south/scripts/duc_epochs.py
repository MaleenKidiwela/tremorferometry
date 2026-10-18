#!/usr/bin/env python
"""DUC_FIX stage 0 -- station epochs 2008-01 .. 2026-09 for HIS permanent stations (FDSN StationXML,
resp/<NET>.<STA>.xml, fetched at response level 2000-01-01..2026-09-30).

An EPOCH is a run of days with one constant KEY = (location, channel codes, sample rate, sensor
model, horizontal azimuths).  Pure gain/response updates are not epoch changes (response removal
handles them).  At each day the channel group used is:
  1. HIS row (stations_permanent.txt: loc + channels, on..off; his off = 2025-01-01 is his
     placeholder for 'still running' and is treated as open) when all his channels are live;
  2. otherwise the replacement at the same station: same component set (1-C vertical / 3-C),
     band preference his band > EH > HH > SH > BH, location his > '' > 00 > 01 > 02 > 03.
status 'his'  = a KEY he used on this station inside his catalogue period (2004-2011) through one
                of his rows -> HIS template kept (sensor swaps inside his period included: he ran
                his templates through them, response removal only).
status 'new'  = any other KEY -> re-templated by duc_retemplate.py from our own high-confidence
                detections made WITHOUT this station, or dropped for that epoch (decided there).
Writes results/duc_fix/station_epochs.csv.
"""
import os, sys, warnings
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from duc_fix_common import *

BANDS = ["EH", "HH", "SH", "BH"]
LOCS = ["", "00", "01", "02", "03"]
KEY0 = pd.Timestamp("2004-01-01")


def chan_epochs(net, sta):
    from obspy import read_inventory
    inv = read_inventory(f"{RESP}/{net}.{sta}.xml")
    rows = []
    for n in inv:
        for s in n:
            if s.code != sta: continue
            for c in s:
                if c.code[:2] not in BANDS: continue
                sens = (c.sensor.description or "") if c.sensor else ""
                rows.append(dict(loc=c.location_code, cha=c.code, start=pd.Timestamp(str(c.start_date)[:10]),
                                 end=pd.Timestamp(str(c.end_date)[:10]) if c.end_date and c.end_date.year < 2900
                                 else pd.Timestamp("2100-01-01"),
                                 sr=float(c.sample_rate), az=float(c.azimuth or 0.0),
                                 sensor=sens.split(",")[0].strip()))
    return pd.DataFrame(rows)


def pick_group(live, his_loc, his_chans):
    """live: channel-epoch rows alive on the day.  Returns (loc, [chans], source) or None."""
    have = {(r.loc, r.cha) for r in live.itertuples()}
    if all((his_loc, c) in have for c in his_chans):
        return his_loc, list(his_chans), "his_row"
    band0 = his_chans[0][:2]
    bands = [band0] + [b for b in BANDS if b != band0]
    locs = [his_loc] + [l for l in LOCS if l != his_loc] + sorted({r.loc for r in live.itertuples()} - set(LOCS))
    threec = len(his_chans) == 3
    for b in bands:
        for l in locs:
            if not threec:
                if (l, b + "Z") in have: return l, [b + "Z"], "replacement"
            else:
                for hz in (("E", "N"), ("1", "2")):
                    cs = [b + hz[0], b + hz[1], b + "Z"]
                    if all((l, c) in have for c in cs): return l, cs, "replacement"
    return None


def main():
    warnings.filterwarnings("ignore")
    H = his_rows()
    out = []
    for (net, sta), hr in H.groupby(["net", "sta"], sort=False):
        CE = chan_epochs(net, sta)
        hr = hr.copy()
        hr["on_t"] = pd.to_datetime(hr.on)
        hr["off_t"] = pd.to_datetime(hr.off).where(hr.off != "2025-01-01", pd.Timestamp("2100-01-01"))
        bps = set([KEY0, pd.Timestamp(T_END)])
        for r in CE.itertuples(): bps |= {r.start, r.end}
        for r in hr.itertuples(): bps |= {r.on_t, r.off_t}
        bps = sorted(b for b in bps if KEY0 <= b <= pd.Timestamp(T_END))
        segs = []
        for a, b in zip(bps[:-1], bps[1:]):
            if b <= a: continue
            live = CE[(CE.start <= a) & (CE.end > a)]
            act = hr[(hr.on_t <= a) & (hr.off_t > a)]
            if len(act): row = act.iloc[-1]
            else:  # his row for the station's most recent listing (chans / loc template)
                prev = hr[hr.on_t <= a]; row = (prev.iloc[-1] if len(prev) else hr.iloc[0])
            hl, hc = row["loc"], row["chans"].split(",")
            g = pick_group(live, hl, hc) if len(live) else None
            if g is None:
                segs.append(dict(a=a, b=b, key=None)); continue
            loc, chans, src = g
            if src == "his_row" and not len(act): src = "replacement"
            L = live[(live["loc"] == loc) & live.cha.isin(chans)].drop_duplicates("cha").set_index("cha")
            srv = float(L.loc[chans[-1], "sr"]); sens = str(L.loc[chans[-1], "sensor"])
            az = tuple(round(float(L.loc[c, "az"]), 1) for c in chans if c[-1] != "Z")
            key = (loc, ",".join(chans), srv, sens, az)
            segs.append(dict(a=a, b=b, key=key, src=src))
        his_keys = {s["key"] for s in segs if s["key"] is not None and s.get("src") == "his_row"
                    and s["a"] < pd.Timestamp("2012-01-01")}
        # merge contiguous equal keys
        M = []
        for s in segs:
            if M and M[-1]["key"] == s["key"] and M[-1]["b"] == s["a"]:
                M[-1]["b"] = s["b"]
            else:
                M.append(dict(s))
        prevkey = None; eid = 0
        for s in M:
            if s["b"] <= pd.Timestamp(T_START) or s["key"] is None:
                if s["key"] is not None: prevkey = s["key"]
                continue
            loc, chans, srv, sens, az = s["key"]
            a = max(s["a"], pd.Timestamp(T_START))
            flags = []
            if prevkey is not None and prevkey != s["key"]:
                for i, nm in enumerate(["loc", "channels", "sample_rate", "sensor", "azimuth"]):
                    if prevkey[i] != s["key"][i]: flags.append(f"{nm}:{prevkey[i]}->{s['key'][i]}")
            status = "his" if s["key"] in his_keys else "new"
            eid += 1
            out.append(dict(net=net, sta=sta, epoch=f"{sta}#{eid}", start=a.date(), end=s["b"].date(),
                            loc=loc, chans=chans, sr=srv, sensor=sens, az="/".join(map(str, az)),
                            source=s["src"], status=status, n_days=(s["b"] - a).days,
                            change=";".join(flags)))
            prevkey = s["key"]
    E = pd.DataFrame(out)
    # sensor swaps inside a 'his' run are invisible in the key only if the model name is equal; list
    # every his-status epoch that follows a sensor change so the report can say it was kept, as he did
    os.makedirs(FIX, exist_ok=True)
    E.to_csv(f"{FIX}/station_epochs.csv", index=False)
    fam = families(); used = set(",".join(fam.stations).split(","))
    E["in_families"] = E.sta.isin(used)
    print(E.groupby("status").agg(epochs=("epoch", "size"), stations=("sta", "nunique"),
                                   days=("n_days", "sum")).to_string())
    pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500)
    print(E[["epoch", "start", "end", "loc", "chans", "sr", "sensor", "status", "change"]].to_string(index=False))


if __name__ == "__main__":
    main()
