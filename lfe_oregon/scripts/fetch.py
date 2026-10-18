"""Fetch raw waveforms for the Thomas & Bostock (2015) central-Oregon reproduction.
Stations: notes/thomas2015_candidate_stations.csv (30 HH/BH/EH stations, 2009-08-03..09-17).
Output raw/<NET.STA>/<YYYY-MM-DD>.mseed (one file per station-day, all 3 components of the best band);
manifest results/fetch_manifest.csv. Resumable: skips existing files; retries twice.
"""
import csv, os, sys, time
from concurrent.futures import ThreadPoolExecutor
from obspy import UTCDateTime as T
from obspy.clients.fdsn import Client

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAS = list(csv.DictReader(open(os.path.join(ROOT, "..", "notes", "thomas2015_candidate_stations.csv"))))
DAYS = [T("2009-08-03") + 86400 * i for i in range(46)]
PREF = ["HH", "BH", "EH"]

def job(args):
    r, d = args
    band = next(b for b in PREF if b in r["bands"].split("/"))
    out = os.path.join(ROOT, "raw", f'{r["net"]}.{r["sta"]}', d.strftime("%Y-%m-%d") + ".mseed")
    if os.path.exists(out) or os.path.exists(out + ".nd"):
        return r["net"], r["sta"], d.date, band, "present", 0
    os.makedirs(os.path.dirname(out), exist_ok=True)
    for k in range(3):
        try:
            st = Client("EARTHSCOPE", timeout=120).get_waveforms(r["net"], r["sta"], "*", band + "?", d, d + 86400)
            st.write(out, format="MSEED")
            return r["net"], r["sta"], d.date, band, "data", len(st)
        except Exception as e:
            if "No data" in str(e) or "204" in str(e):
                open(out + ".nd", "w").close()
                return r["net"], r["sta"], d.date, band, "nodata", 0
            time.sleep(5 * (k + 1))
    return r["net"], r["sta"], d.date, band, "error", 0

if __name__ == "__main__":
    jobs = [(r, d) for r in STAS for d in DAYS]
    man = os.path.join(ROOT, "results", "fetch_manifest.csv")
    with open(man, "w", newline="") as f, ThreadPoolExecutor(6) as ex:
        w = csv.writer(f); w.writerow(["net", "sta", "date", "band", "status", "n_traces"])
        for i, row in enumerate(ex.map(job, jobs), 1):
            w.writerow(row); f.flush()
            if i % 100 == 0:
                print(f"{time.strftime('%H:%M:%S')} {i}/{len(jobs)}", flush=True)
    open(os.path.join(ROOT, "results", "FETCH_DONE"), "w").write(time.strftime("%FT%TZ\n"))
    print("done", flush=True)
