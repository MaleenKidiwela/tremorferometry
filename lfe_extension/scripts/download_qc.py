#!/usr/bin/env python
"""Track B step 4 - per-month download QC + the station x month coverage figure.

  python scripts/download_qc.py --month 2014-01      one month  -> results/download_qc_2014-01.csv
  python scripts/download_qc.py --all                every month present in the manifest
  python scripts/download_qc.py --fig-only           just redraw figures/download_coverage.png

Per station and month: band/channels/fs, days with data, no-data days, error/missing days,
median & min coverage fraction, days with coverage < 0.9, fragments per channel, MB.
Days served from ../lfe_catalog_repair/traces (status present_elsewhere) are read from
their recorded path so the month is scored on ALL the data Track C will see.
"""
import os
import csv
import glob
import argparse
import statistics as stats
from collections import defaultdict
from datetime import date, timedelta

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results")
FIGS = os.path.join(HERE, "figures")
MANIFEST = os.path.join(RESULTS, "download_manifest.csv")
QUEUE = os.path.join(RESULTS, "download_queue.csv")


def read_manifest():
    rows = {}
    if not os.path.exists(MANIFEST):
        return rows
    with open(MANIFEST) as f:
        for r in csv.DictReader(f):
            if not r.get("date"):
                continue
            rows[(r["sta"], r["date"])] = r            # last row for a key wins
    return rows


def queue_index():
    q = defaultdict(list)
    with open(QUEUE) as f:
        for r in csv.DictReader(f):
            if r["status"] == "no_epoch":
                continue
            q[(r["sta"], r["date"][:7])].append(r)
    return q


def enrich_present_elsewhere(rec):
    """present_elsewhere rows carry no coverage - read the foreign file's headers once."""
    p = rec.get("note", "")
    if not p or not os.path.exists(p):
        return rec
    try:
        import sys
        sys.path.insert(0, os.path.join(HERE, "scripts"))
        from download_day import verify
        n, fs, cov, nf, _ = verify(p, rec.get("band") or p[-20:-15])
        rec = dict(rec)
        rec.update(n_channels=n, fs=f"{fs:g}", coverage_frac=f"{cov:.4f}", n_frag=nf,
                   bytes=os.path.getsize(p))
    except Exception:
        pass
    return rec


def month_days(month):
    y, m = int(month[:4]), int(month[5:7])
    d = date(y, m, 1)
    out = []
    while d.month == m:
        out.append(d.isoformat())
        d += timedelta(days=1)
    return out


def qc_month(month, man, qidx, verbose=True):
    out = []
    stas = sorted({k[0] for k in qidx if k[1] == month})
    for sta in stas:
        qrows = qidx[(sta, month)]
        exp = len(qrows)
        net = qrows[0]["net"]
        bands = sorted({r["band"] for r in qrows if r["band"]})
        recs = []
        n_pe = 0
        for r in qrows:
            rec = man.get((sta, r["date"]))
            if rec is None:
                recs.append(dict(status="missing"))
                continue
            if rec["status"] == "present_elsewhere":
                rec = enrich_present_elsewhere(rec)
                n_pe += 1
            recs.append(rec)
        cov = [float(r["coverage_frac"]) for r in recs
               if r.get("status") in ("data", "present_elsewhere") and r.get("coverage_frac")]
        frag = [int(r["n_frag"]) for r in recs if r.get("n_frag")]
        nch = sorted({int(r["n_channels"]) for r in recs if r.get("n_channels")})
        fss = sorted({r["fs"] for r in recs if r.get("fs")})
        mb = sum(int(r["bytes"]) for r in recs if r.get("bytes")) / 1e6
        out.append(dict(
            month=month, net=net, sta=sta, band=",".join(bands), fs=",".join(fss),
            n_channels=",".join(str(x) for x in nch),
            days_expected=exp,
            days_data=sum(1 for r in recs if r.get("status") in ("data", "present_elsewhere")),
            days_present_elsewhere=n_pe,
            days_nodata=sum(1 for r in recs if r.get("status") == "nodata"),
            days_error=sum(1 for r in recs if r.get("status") == "error"),
            days_missing=sum(1 for r in recs if r.get("status") == "missing"),
            cov_median=round(stats.median(cov), 4) if cov else "",
            cov_min=round(min(cov), 4) if cov else "",
            days_cov_lt_0p9=sum(1 for c in cov if c < 0.9),
            frag_median=round(stats.median(frag), 1) if frag else "",
            frag_max=max(frag) if frag else "",
            mb=round(mb, 1)))
    p = os.path.join(RESULTS, f"download_qc_{month}.csv")
    if out:
        tmp = p + ".tmp"
        with open(tmp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
            w.writeheader()
            w.writerows(out)
        os.replace(tmp, p)
        if verbose:
            print(f"{p}: {len(out)} stations, "
                  f"{sum(r['days_data'] for r in out)} station-days with data, "
                  f"{round(sum(r['mb'] for r in out)/1000,2)} GB")
    return out


def figure(man=None):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    rows = []
    for p in sorted(glob.glob(os.path.join(RESULTS, "download_qc_*.csv"))):
        with open(p) as f:
            rows += list(csv.DictReader(f))
    if not rows:
        print("no QC csvs yet - no figure")
        return
    stas = sorted({f"{r['net']}.{r['sta']}" for r in rows})
    months = sorted({r["month"] for r in rows})
    M = np.full((len(stas), len(months)), np.nan)
    for r in rows:
        i, j = stas.index(f"{r['net']}.{r['sta']}"), months.index(r["month"])
        exp = max(1, int(r["days_expected"]))
        cm = float(r["cov_median"]) if r["cov_median"] else 0.0
        M[i, j] = cm * int(r["days_data"]) / exp        # fraction of the month actually usable
    fig, ax = plt.subplots(figsize=(max(8, 0.22 * len(months) + 3), 0.32 * len(stas) + 2))
    im = ax.imshow(M, aspect="auto", vmin=0, vmax=1, cmap="viridis",
                   interpolation="nearest")
    ax.set_yticks(range(len(stas)))
    ax.set_yticklabels(stas, fontsize=7)
    step = max(1, len(months) // 40)
    ax.set_xticks(range(0, len(months), step))
    ax.set_xticklabels(months[::step], rotation=90, fontsize=6)
    ax.set_title("lfe_extension Track B - usable fraction of each station-month "
                 "(median coverage x days with data / days in epoch)", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.7, label="usable fraction")
    fig.tight_layout()
    os.makedirs(FIGS, exist_ok=True)
    p = os.path.join(FIGS, "download_coverage.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print("figure ->", p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--month")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--fig-only", action="store_true")
    a = ap.parse_args()
    if a.fig_only:
        figure()
        return
    man, qidx = read_manifest(), queue_index()
    months = [a.month] if a.month else sorted({k[1] for k in qidx}) if a.all else []
    if not months:
        ap.error("give --month YYYY-MM, --all or --fig-only")
    for m in months:
        qc_month(m, man, qidx)
    figure()


if __name__ == "__main__":
    main()
