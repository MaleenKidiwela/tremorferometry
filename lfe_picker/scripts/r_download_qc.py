#!/usr/bin/env python
"""P3D step 4 - per-station QC of a priority class + the station x month coverage figure.

  python scripts/r_download_qc.py --priority P0     -> results/p3d/download_qc_P0.csv
  python scripts/r_download_qc.py --all             -> one csv per priority class
  python scripts/r_download_qc.py --fig-only        -> just redraw figures/r_download_coverage.png

Adapted from scripts/o_download_qc.py, but scored PER STATION over a priority class (the P3D
question is "did this station actually deliver data in the SSE windows?", not "is this month
complete"), and the verdict column is explicit: metadata "open" is not delivery.

  delivered = YES   >=50 % of the in-epoch days of the class came back with data
              PARTIAL  1..50 %
              NO    0 days with data although the station has an epoch over the class
              N/A   the station has no epoch inside the class at all
"""
import os
import csv
import argparse
import statistics as stats
from collections import defaultdict

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(HERE, "results", "p3d")
FIGS = os.path.join(HERE, "figures")
MANIFEST = os.path.join(RESULTS, "download_manifest.csv")
QUEUE = os.path.join(RESULTS, "download_queue.csv")


def read_manifest():
    rows = {}
    with open(MANIFEST) as f:
        for r in csv.DictReader(f):
            if r.get("date"):
                rows[(r["sta"], r["date"])] = r        # last row for a key wins
    return rows


def read_queue():
    q = []
    with open(QUEUE) as f:
        q = list(csv.DictReader(f))
    return q


def qc_priority(prio, man, q, verbose=True):
    out = []
    sel = [r for r in q if r["priority"] == prio]
    stas = sorted({(r["net"], r["sta"]) for r in sel})
    for net, sta in stas:
        rows = [r for r in sel if r["sta"] == sta]
        inep = [r for r in rows if r["status"] == "pending"]
        noep = len(rows) - len(inep)
        recs = []
        for r in inep:
            rec = man.get((sta, r["date"]))
            recs.append(rec if rec is not None else dict(status="missing"))
        cov = [float(r["coverage_frac"]) for r in recs
               if r.get("status") == "data" and r.get("coverage_frac")]
        frag = [int(r["n_frag"]) for r in recs if r.get("status") == "data" and r.get("n_frag")]
        nch = sorted({int(r["n_channels"]) for r in recs
                      if r.get("status") == "data" and r.get("n_channels")})
        fss = sorted({r["fs"] for r in recs if r.get("status") == "data" and r.get("fs")})
        mb = sum(int(r["bytes"]) for r in recs if r.get("bytes"))
        dd = sum(1 for r in recs if r.get("status") == "data")
        dates = sorted(r["date"] for r in recs if r.get("status") == "data")
        exp = len(inep)
        if exp == 0:
            verdict = "N/A"
        elif dd == 0:
            verdict = "NO"
        elif dd >= 0.5 * exp:
            verdict = "YES"
        else:
            verdict = "PARTIAL"
        out.append(dict(
            priority=prio, net=net, sta=sta,
            band=",".join(sorted({r["band"] for r in inep if r["band"]})),
            fs=",".join(fss), n_channels=",".join(str(x) for x in nch),
            days_in_epoch=exp, days_no_epoch=noep, days_data=dd,
            days_nodata=sum(1 for r in recs if r.get("status") == "nodata"),
            days_error=sum(1 for r in recs if r.get("status") == "error"),
            days_missing=sum(1 for r in recs if r.get("status") == "missing"),
            frac_delivered=round(dd / exp, 4) if exp else "",
            delivered=verdict,
            first_data=dates[0] if dates else "", last_data=dates[-1] if dates else "",
            cov_median=round(stats.median(cov), 4) if cov else "",
            cov_min=round(min(cov), 4) if cov else "",
            days_cov_lt_0p9=sum(1 for c in cov if c < 0.9),
            frag_median=round(stats.median(frag), 1) if frag else "",
            frag_max=max(frag) if frag else "",
            gb=round(mb / 1e9, 3),
            mb_per_day=round(mb / 1e6 / dd, 1) if dd else ""))
    p = os.path.join(RESULTS, f"download_qc_{prio}.csv")
    if out:
        tmp = p + ".tmp"
        with open(tmp, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
            w.writeheader()
            w.writerows(out)
        os.replace(tmp, p)
        if verbose:
            print(f"{p}")
            hdr = ("net.sta", "band", "in_ep", "data", "nodata", "err", "miss",
                   "frac", "verdict", "covmed", "fragmed", "GB")
            print("  {:<10}{:<6}{:>7}{:>7}{:>8}{:>6}{:>6}{:>7}  {:<8}{:>8}{:>8}{:>8}"
                  .format(*hdr))
            for r in out:
                print("  {:<10}{:<6}{:>7}{:>7}{:>8}{:>6}{:>6}{:>7}  {:<8}{:>8}{:>8}{:>8}".format(
                    f"{r['net']}.{r['sta']}", r["band"], r["days_in_epoch"], r["days_data"],
                    r["days_nodata"], r["days_error"], r["days_missing"],
                    str(r["frac_delivered"]), r["delivered"], str(r["cov_median"]),
                    str(r["frag_median"]), r["gb"]))
            print(f"  TOTAL {sum(r['days_data'] for r in out)} station-days with data, "
                  f"{round(sum(r['gb'] for r in out), 2)} GB")
    return out


def figure(q, prios=None, tag=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    man = read_manifest()
    # expected in-epoch days per station-month from the queue
    exp = defaultdict(int)
    got = defaultdict(list)
    prio_of = {}
    for r in q:
        key = (f"{r['net']}.{r['sta']}", r["date"][:7])
        prio_of.setdefault(r["date"][:7], set()).add(r["priority"])
        if r["status"] != "pending":
            continue
        if prios and r["priority"] not in prios:
            continue
        exp[key] += 1
        rec = man.get((r["sta"], r["date"]))
        if rec is not None and rec.get("status") == "data" and rec.get("coverage_frac"):
            got[key].append(float(rec["coverage_frac"]))
    stas = sorted({k[0] for k in exp})
    months = sorted({k[1] for k in exp})
    M = np.full((len(stas), len(months)), np.nan)
    for (s, m), e in exp.items():
        v = got.get((s, m), [])
        M[stas.index(s), months.index(m)] = (sum(v) / e) if e else np.nan
    fig, ax = plt.subplots(figsize=(max(9, 0.16 * len(months) + 3), 0.34 * len(stas) + 2.4))
    im = ax.imshow(M, aspect="auto", vmin=0, vmax=1, cmap="viridis", interpolation="nearest")
    ax.set_yticks(range(len(stas)))
    ax.set_yticklabels(stas, fontsize=8)
    step = max(1, len(months) // 46)
    ax.set_xticks(range(0, len(months), step))
    ax.set_xticklabels(months[::step], rotation=90, fontsize=6)
    # mark the months that contain P0 (SSE) days
    for j, m in enumerate(months):
        if "P0" in prio_of.get(m, set()):
            ax.plot([j], [-0.75], marker="v", color="crimson", ms=3, clip_on=False)
    lab = ("+".join(sorted(prios)) if prios else "all priorities")
    ax.set_title(f"P3D download [{lab}] - delivered fraction of each station-month\n"
                 f"(sum of daily coverage / in-epoch {lab} days queued; "
                 "red = month contains SSE days); white = nothing queued / no epoch", fontsize=9)
    fig.colorbar(im, ax=ax, shrink=0.7, label="delivered fraction")
    fig.tight_layout()
    os.makedirs(FIGS, exist_ok=True)
    p = os.path.join(FIGS, f"r_download_coverage{tag}.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    print("figure ->", p)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--priority", default="P0")
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--fig-only", action="store_true")
    ap.add_argument("--no-fig", action="store_true")
    a = ap.parse_args()
    q = read_queue()
    if a.fig_only:
        figure(q, {a.priority} if a.priority else None,
               f"_{a.priority}" if a.priority else "")
        figure(q, None, "_all")
        return
    man = read_manifest()
    for prio in (sorted({r["priority"] for r in q}) if a.all else [a.priority]):
        qc_priority(prio, man, q)
    if not a.no_fig:
        figure(q, {a.priority}, "")          # the headline figure: the QC'd class only
        figure(q, None, "_all")              # every priority, for the running fill


if __name__ == "__main__":
    main()
