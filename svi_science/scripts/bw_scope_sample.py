#!/usr/bin/env python
"""bw_scope_sample -- SVI backward-extension scope, step 2: small FDSN dataselect samples.

For every permanent-network station in inventory_stations.csv (temporary arrays XG/XU/YW/4P and strong-motion-only
stations are listed but not sampled) and every year 2003..2013: four 10-minute requests
(YYYY-02-15, -05-15, -08-15, -11-15, 10:00-10:10 UTC; seismic velocity bands HH?,BH?,EH?,SH?) at EVERY centre whose
metadata lists the station (NRCan/CNDC and/or EarthScope).  Nothing is written to disk except the result table:
waveforms are parsed in memory (obspy) and dropped.

Out: svi_science/results/backward_scope/samples_raw.csv (one row per request)
"""
import io, sys, time, threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np, pandas as pd, requests

ROOT = "/home/jovyan/tremorferometry"
OUT = f"{ROOT}/svi_science/results/backward_scope"
CENTRES = {"NRCAN": "https://www.earthquakescanada.nrcan.gc.ca/fdsnws",
           "EARTHSCOPE": "https://service.earthscope.org/fdsnws"}
SEM = {"NRCAN": threading.Semaphore(4), "EARTHSCOPE": threading.Semaphore(8)}
YEARS = range(2003, 2014)
MD = ["02-15", "05-15", "08-15", "11-15"]
TEMP = {"XG", "XU", "YW", "4P", "X6"}


def fetch(centre, net, sta, day):
    t0 = f"{day}T10:00:00"; t1 = f"{day}T10:10:00"
    url = (f"{CENTRES[centre]}/dataselect/1/query?net={net}&sta={sta}&cha=HH?,BH?,EH?,SH?"
           f"&starttime={t0}&endtime={t1}")
    rec = dict(centre=centre, net=net, sta=sta, date=day, http="", bytes=0, chans="", fs="",
               cov_min=0.0, flat_frac_max="", rms_z="", note="")
    r = None
    with SEM[centre]:
        for k in range(4):
            try:
                r = requests.get(url, timeout=120)
                if r.status_code in (200, 204):
                    break
                if r.status_code == 404:          # FDSN 'no data' (nodata=404 default at some centres)
                    # NRCan has returned transient 404s for data that exist; retry once more after a pause
                    if k >= 1:
                        break
            except Exception as e:
                rec["note"] = f"exc:{type(e).__name__}"
            time.sleep(4 * (k + 1))
    if r is None:
        rec["http"] = "fail"; return rec
    rec["http"] = r.status_code; rec["bytes"] = len(r.content)
    if r.status_code != 200 or not r.content:
        return rec
    try:
        from obspy import read
        st = read(io.BytesIO(r.content))
        by = {}
        for tr in st:
            by.setdefault(tr.stats.channel, []).append(tr)
        chans, fss, covs, flats = [], [], [], []
        rmsz = ""
        for ch, trs in sorted(by.items()):
            fs = float(trs[0].stats.sampling_rate)
            n = sum(t.stats.npts for t in trs)
            cov = min(n / (600.0 * fs), 1.0)
            x = np.concatenate([t.data.astype(float) for t in trs])
            flat = float(np.mean(np.diff(x) == 0)) if x.size > 1 else 1.0
            chans.append(ch); fss.append(f"{fs:g}"); covs.append(cov); flats.append(flat)
            if ch.endswith("Z") and x.size > 10:
                rmsz = f"{np.std(x - x.mean()):.4g}"
        rec.update(chans=";".join(chans), fs=";".join(sorted(set(fss))), cov_min=round(min(covs), 4),
                   flat_frac_max=round(max(flats), 4), rms_z=rmsz)
    except Exception as e:
        rec["note"] = f"parse:{type(e).__name__}"
    return rec


def main():
    S = pd.read_csv(f"{OUT}/inventory_stations.csv")
    S = S[~S.net.isin(TEMP)]
    S = S[S.bands.str.contains("HH|BH|EH|SH")]
    jobs = []
    for _, s in S.iterrows():
        first = str(s["first"])[:10]; last = str(s["last"])[:10]
        for c in s.centres.split(";"):
            for y in YEARS:
                for md in MD:
                    d = f"{y}-{md}"
                    if d < first[:8] + "01" and d < first:      # outside metadata window -> still record
                        jobs.append((c, s.net, s.sta, d, "no_epoch")); continue
                    if d > last:
                        jobs.append((c, s.net, s.sta, d, "no_epoch")); continue
                    jobs.append((c, s.net, s.sta, d, ""))
    print(len(jobs), "jobs", sum(j[4] == "" for j in jobs), "requests", flush=True)
    rows = [dict(centre=c, net=n, sta=st, date=d, http="no_epoch") for c, n, st, d, f in jobs if f]
    todo = [j for j in jobs if not j[4]]
    with ThreadPoolExecutor(12) as ex:
        futs = [ex.submit(fetch, c, n, st, d) for c, n, st, d, _ in todo]
        for i, f in enumerate(as_completed(futs)):
            rows.append(f.result())
            if i % 200 == 0:
                print(i, time.strftime("%H:%M:%S"), flush=True)
                pd.DataFrame(rows).to_csv(f"{OUT}/samples_raw.partial.csv", index=False)
    D = pd.DataFrame(rows).sort_values(["net", "sta", "centre", "date"])
    D.to_csv(f"{OUT}/samples_raw.csv", index=False)
    print("done", len(D), flush=True)


if __name__ == "__main__":
    main()
