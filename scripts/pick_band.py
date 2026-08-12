#!/usr/bin/env python
"""Pick ONE vertical channel band for a station's fleet run: the Z band spanning the most of the tremor era
(2010-2026), preferring native-40 BHZ (no decimation). Multi-band stations get their record-spanning band for
era-1; a second band (if it covers a distinct later era) is a per-era follow-up. Prints e.g. 'BH?'.
Usage: python scripts/pick_band.py <NET> <STA>"""
import sys
from obspy.clients.fdsn import Client

net, sta = sys.argv[1], sys.argv[2]
inv = None
provs = (["NCEDC", "https://service.ncedc.org"] if net in ("BK", "NC")   # N. California networks live at NCEDC
         else ["IRIS", "https://service.iris.edu", "EARTHSCOPE"])
for prov in provs:
    try:
        inv = Client(prov, timeout=30).get_stations(network=net, station=sta, level="channel"); break
    except Exception:
        pass
from collections import defaultdict
# aggregate span PER BAND across all metadata epochs (a band's HHZ can be split into several epochs; scoring
# each epoch separately let a short single-epoch BHZ beat a long multi-epoch HHZ -> picked the short era).
band_span, band_n40 = defaultdict(float), {}
if inv:
    for ch in inv[0][0].channels:
        # require [BHES]H? = high-gain SEISMOMETER (2nd char H). Excludes N=accelerometer (HNZ/BNZ), L=low-gain.
        if not ch.code.endswith("Z") or ch.code[0] not in "BHES" or ch.code[1] != "H":
            continue
        s = ch.start_date.year
        e = ch.end_date.year if ch.end_date else 2026
        b = ch.code[:2]
        band_span[b] += max(0, min(e, 2026) - max(s, 2010))
        band_n40[b] = abs(ch.sample_rate - 40) < 1


def _score(b):
    return band_span[b] + (2 if band_n40.get(b) else 0) + (100 if b[0] in "BH" else 0)   # prefer broadband; native-40 tiebreak


print((max(band_span, key=_score) + "?") if band_span else "BH?,HH?,EH?,SH?")
