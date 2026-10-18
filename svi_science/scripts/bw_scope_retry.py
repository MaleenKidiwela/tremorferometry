#!/usr/bin/env python
"""bw_scope_retry -- SVI backward scope, step 2b: serial re-try of every NRCan sample that did not return data.

NRCan/CNDC FDSN (www.earthquakescanada.nrcan.gc.ca) answers the SAME request alternately 404 "Not found" and
200 (seen 2026-10-07: PGC 2008-03-15 BHZ 404 then 200 x5; BH? 404 x6), so one miss is not "no data".
For every NRCan row of samples_raw.csv with http != 200, re-request up to 5 times (3 parallel lanes, each serial),
1.5 s apart, cycling the channel pattern over the exact Z channel(s) active that day in the NRCan metadata
(e.g. BHZ, HHZ) and the band wildcard.  Z-only answers count as 'delivered' (availability, not 3-C QC).
Also re-tries EarthScope 204s for CN once (to confirm the 204 is stable).

Out: svi_science/results/backward_scope/samples_retry.csv
"""
import io, time
import numpy as np, pandas as pd, requests

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/backward_scope"
BASE = {"NRCAN": "https://www.earthquakescanada.nrcan.gc.ca/fdsnws", "EARTHSCOPE": "https://service.earthscope.org/fdsnws"}


def one(centre, net, sta, day, cha):
    url = (f"{BASE[centre]}/dataselect/1/query?net={net}&sta={sta}&cha={cha}"
           f"&starttime={day}T10:00:00&endtime={day}T10:10:00")
    try:
        r = requests.get(url, timeout=90)
    except Exception as e:
        return "exc", 0, ""
    if r.status_code != 200 or not r.content:
        return r.status_code, 0, ""
    try:
        from obspy import read
        st = read(io.BytesIO(r.content))
        return 200, len(r.content), ";".join(sorted({f"{t.stats.channel}@{t.stats.sampling_rate:g}" for t in st}))
    except Exception:
        return "parse", len(r.content), ""


def main():
    S = pd.read_csv(f"{OUT}/samples_raw.csv")
    C = pd.read_csv(f"{OUT}/inventory_channels.csv")
    C["loc"] = C["loc"].fillna("")
    todo = S[((S.centre == "NRCAN") & (S.http.astype(str) != "200") & (S.http.astype(str) != "no_epoch")) |
             ((S.centre == "EARTHSCOPE") & (S.net == "CN") & (S.http.astype(str) == "204"))]
    print(len(todo), "to retry", flush=True)
    from concurrent.futures import ThreadPoolExecutor
    rows = []

    def work(r):
        c = C[(C.centre == r.centre) & (C.net == r.net) & (C.sta == r.sta) & (C.start.str[:10] <= r.date) &
              (C.end.str[:10] >= r.date) & C.cha.str.endswith(("Z", "3"))]
        pats = sorted(set(c.cha)) or ["HHZ", "BHZ", "EHZ"]
        pats = pats + [p[:2] + "?" for p in pats]
        ntry = 5 if r.centre == "NRCAN" else 1
        got = None; log = []
        for k in range(ntry):
            p = pats[k % len(pats)]
            code, nb, ch = one(r.centre, r.net, r.sta, r.date, p)
            log.append(f"{p}:{code}")
            if code == 200:
                got = (nb, ch); break
            time.sleep(1.5)
        return dict(centre=r.centre, net=r.net, sta=r.sta, date=r.date, first_http=r.http,
                    retry_ok=got is not None, bytes=got[0] if got else 0, chans=got[1] if got else "",
                    attempts=len(log), log=" ".join(log))

    with ThreadPoolExecutor(3) as ex:          # 3 lanes; each lane is strictly serial with 1.5-s pauses
        for i, rec in enumerate(ex.map(work, list(todo.itertuples()))):
            rows.append(rec)
            if i % 100 == 0:
                print(i, time.strftime("%H:%M:%S"), flush=True)
                pd.DataFrame(rows).to_csv(f"{OUT}/samples_retry.partial.csv", index=False)
    pd.DataFrame(rows).to_csv(f"{OUT}/samples_retry.csv", index=False)
    print("done", flush=True)


if __name__ == "__main__":
    main()
