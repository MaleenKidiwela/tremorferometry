#!/usr/bin/env python
"""g2_local_dl (2026-10-11) -- Study G2 step 0: download the local earthquake catalogues (small queries).

ComCat (USGS FDSN, includes PNSN 'uw') CSV by year, all event types, M >= 1.0, box 45.9-50.5 N, 126.8-120.4 W.
NRCan FDSN event (QuakeML, for event type) by month, same box, M >= 1.0; flaky -> up to 8 retries per request.
Raw files -> results/g2_local/dl/{comcat,nrcan}/ ; parsed NRCan -> dl/nrcan_events.csv.  No computation.

    python -I svi_science/scripts/g2_local_dl.py
"""
import os, re, time, subprocess
import pandas as pd

ROOT = "/home/jovyan/tremorferometry"
DL = f"{ROOT}/svi_science/results/g2_local/dl"
BOX = dict(minlatitude=45.9, maxlatitude=50.5, minlongitude=-126.8, maxlongitude=-120.4, minmagnitude=1.0)


def get(url, out, tries=8):
    for k in range(tries):
        r = subprocess.run(["curl", "-s", "-m", "180", "-o", out, "-w", "%{http_code}", url], capture_output=True, text=True)
        code = r.stdout.strip()
        if code in ("200", "204") and os.path.exists(out):
            return code
        time.sleep(3 + 4 * k)
    return code


def q(d):
    return "&".join(f"{k}={v}" for k, v in d.items())


def comcat():
    os.makedirs(f"{DL}/comcat", exist_ok=True)
    for y in range(2003, 2027):
        out = f"{DL}/comcat/comcat_{y}.csv"
        if os.path.exists(out) and os.path.getsize(out) > 100:
            continue
        url = ("https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&orderby=time-asc&"
               + q(dict(BOX, starttime=f"{y}-01-01", endtime=f"{y + 1}-01-01" if y < 2026 else "2026-10-11")))
        print(y, get(url, out), flush=True)


EV = re.compile(r'<event publicID="smi:local/([^"]+)">(.*?)</event>', re.S)


def _f(pat, s):
    m = re.search(pat, s, re.S)
    return m.group(1) if m else None


def nrcan():
    os.makedirs(f"{DL}/nrcan", exist_ok=True)
    rows, bad = [], []
    for y in range(2003, 2027):
        for mth in range(1, 13):
            if y == 2026 and mth > 10:
                break
            s = pd.Timestamp(year=y, month=mth, day=1); e = s + pd.offsets.MonthBegin(1)
            out = f"{DL}/nrcan/nrcan_{y}_{mth:02d}.xml"
            if not (os.path.exists(out) and os.path.getsize(out) > 0):
                url = ("https://www.earthquakescanada.nrcan.gc.ca/fdsnws/event/1/query?format=xml&"
                       + q(dict(BOX, starttime=s.strftime("%Y-%m-%d"), endtime=e.strftime("%Y-%m-%d"))))
                code = get(url, out)
                if code not in ("200", "204"):
                    bad.append((y, mth, code)); print("FAIL", y, mth, code, flush=True)
                    continue
            txt = open(out, encoding="utf-8", errors="replace").read() if os.path.exists(out) else ""
            for eid, body in EV.findall(txt):
                rows.append(dict(id=eid, type=_f(r"^<type>([^<]*)</type>", body),
                                 time=_f(r"<time><value>([^<]+)</value>", body),
                                 lat=_f(r"<latitude><value>([^<]+)</value>", body),
                                 lon=_f(r"<longitude><value>([^<]+)</value>", body),
                                 depth=_f(r"<depth><value>([^<]+)</value>", body),
                                 depthType=_f(r"<depthType>([^<]+)</depthType>", body),
                                 mag=_f(r"<magnitude[^>]*><mag><value>([^<]+)</value>", body),
                                 magType=_f(r"</mag><type>([^<]+)</type>", body),
                                 place=_f(r"<description><text>([^<]*)</text>", body)))
        print("nrcan", y, len(rows), flush=True)
    D = pd.DataFrame(rows)
    D.to_csv(f"{DL}/nrcan_events.csv", index=False)
    pd.DataFrame(bad, columns=["year", "month", "code"]).to_csv(f"{DL}/nrcan_failed_months.csv", index=False)
    print("NRCan rows", len(D), "failed months", len(bad))


if __name__ == "__main__":
    comcat()
    nrcan()
